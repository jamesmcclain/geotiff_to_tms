/*
 * Copyright (c) 2017, James McClain
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 * 3. All advertising materials mentioning features or use of this
 *    software must display the following acknowledgement: This product
 *    includes software developed by Dr. James W. McClain.
 * 4. Neither the names of the authors nor the names of the
 *    contributors may be used to endorse or promote products derived
 *    from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHORS ''AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <arpa/inet.h>
#include <float.h>
#include <math.h>

#include "ansi.h"
#include "constants.h"
#include "landsat.h"
#include "load.h"
#include "pngwrite.h"

#include "gdal.h"
#include "cpl_conv.h"
#include "ogr_srs_api.h"
#include "proj_api.h"


const char * webmercator = "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext +no_defs";
projPJ webmercator_pj = NULL;
landsat_scene scene = {.filename = "/tmp/LC08_L1TP_139045_20170304_20170316_01_T1_B%d.TIF"};
/* landsat_scene scene = {.filename = "/vsicurl/https://landsat-pds.s3.amazonaws.com/c1/L8/139/045/LC08_L1TP_139045_20170304_20170316_01_T1/LC08_L1TP_139045_20170304_20170316_01_T1_B%d.TIF"}; */

int fetch(double xmin, double xmax, double ymin, double ymax, int verbose);
uint8_t sigmoidal(uint16_t v);
void world_to_image(double * xy);
void zxy_approx(int fd, int z, int _x, int _y, int verbose);
void zxy_exact(int fd, int z, int _x, int _y, int verbose);


void load(int verbose)
{
  char * wkt = NULL, * dstProj4 = NULL;
  OGRSpatialReferenceH srs = NULL;
  char r_filename[(1<<8)];
  char g_filename[(1<<8)];
  char b_filename[(1<<8)];

  sprintf(r_filename, scene.filename, 4);
  sprintf(g_filename, scene.filename, 3);
  sprintf(b_filename, scene.filename, 2);

  GDALAllRegister();

  /* Datasets and bands */
  scene.r_dataset = GDALOpen(r_filename, GA_ReadOnly);
  if(scene.r_dataset == NULL) {
    fprintf(stderr, ANSI_COLOR_RED "GDALOpen problem (red band)" ANSI_COLOR_RESET "\n");
    exit(-1);
  }
  scene.g_dataset = GDALOpen(g_filename, GA_ReadOnly);
  if(scene.g_dataset == NULL) {
    fprintf(stderr, ANSI_COLOR_RED "GDALOpen problem (green)" ANSI_COLOR_RESET "\n");
    exit(-1);
  }
  scene.b_dataset = GDALOpen(b_filename, GA_ReadOnly);
  if(scene.b_dataset == NULL) {
    fprintf(stderr, ANSI_COLOR_RED "GDALOpen problem (blue)" ANSI_COLOR_RESET "\n");
    exit(-1);
  }
  scene.r_band = GDALGetRasterBand(scene.r_dataset, 1);
  scene.g_band = GDALGetRasterBand(scene.g_dataset, 1);
  scene.b_band = GDALGetRasterBand(scene.b_dataset, 1);
  scene.width  = GDALGetRasterXSize(scene.r_dataset);
  scene.height = GDALGetRasterYSize(scene.r_dataset);

  /* SRS.  This is from the red band, but is assumed to be valid for
     all of the bands. */
  srs = OSRNewSpatialReference(NULL);
  wkt = calloc(STRING_BUFFER_SIZE, sizeof(char)); // No memory leak, this is freed from within `OSRImportFromWkt`!
  strncpy(wkt, GDALGetProjectionRef(scene.r_dataset), STRING_BUFFER_SIZE);
  if (verbose)
    fprintf(stderr, ANSI_COLOR_GREEN "WKT: " ANSI_COLOR_CYAN "%s" ANSI_COLOR_RESET "\n", wkt);
  OSRImportFromWkt(srs, &wkt);
  OSRExportToProj4(srs, &dstProj4);
  if (verbose)
    fprintf(stderr, ANSI_COLOR_GREEN "Proj4: " ANSI_COLOR_CYAN "%s" ANSI_COLOR_RESET "\n", dstProj4);

  webmercator_pj = pj_init_plus(webmercator);
  scene.destination_pj = pj_init_plus(dstProj4);

  /* Transform.  This is from the red band, but is assumed to be valid
     for all of the bands. */
  GDALGetGeoTransform(scene.r_dataset, scene.transform);

  CPLFree(dstProj4);
  OSRRelease(srs);
}

int fetch(double xmin, double xmax, double ymin, double ymax, int verbose)
{
  int startx = (int)xmin, starty = (int)ymin;
  int rsizex = (int)(xmax-xmin), rsizey = (int)(ymax-ymin);
  int wsizex = TEXTURE_BUFFER_SIZE, wsizey = TEXTURE_BUFFER_SIZE;
  int deltax = 0, deltay = 0;

  memset(scene.tile, 0, sizeof(scene.tile));
  memset(scene.r_texture, 0, sizeof(scene.r_texture));
  memset(scene.g_texture, 0, sizeof(scene.g_texture));
  memset(scene.b_texture, 0, sizeof(scene.b_texture));

  if (startx < 0) {
    rsizex += startx;
    startx = 0;
    deltax = TEXTURE_BUFFER_SIZE-(int)(wsizex * (rsizex/(xmax-xmin)));
  }
  if (starty < 0) {
    rsizey += starty;
    starty = 0;
    deltay = TEXTURE_BUFFER_SIZE-(int)(wsizey * (rsizey/(ymax-ymin)));
  }
  if (startx + rsizex > scene.width) {
    rsizex = scene.width - startx;
    rsizex = rsizex > 0? rsizex : 0;
  }
  if (starty + rsizey > scene.height) {
    rsizey = scene.height - starty;
    rsizey = rsizey > 0? rsizey : 0;
  }

  wsizex = (int)(wsizex * (rsizex/(xmax-xmin)));
  wsizey = (int)(wsizey * (rsizey/(ymax-ymin)));

  if (!rsizex || !rsizey || !wsizex || !wsizey)
    return 0;

  if (verbose) {
    fprintf(stderr,
            ANSI_COLOR_GREEN "X "
            ANSI_COLOR_CYAN "start: %d, rsize: %d, wsize: %d, delta: %d"
            ANSI_COLOR_RESET "\n",
            startx, rsizex, wsizex, deltax);
    fprintf(stderr,
            ANSI_COLOR_GREEN "Y "
            ANSI_COLOR_CYAN "start: %d, rsize: %d, wsize: %d, delta: %d"
            ANSI_COLOR_RESET "\n",
            starty, rsizey, wsizey, deltay);
  }

  if (GDALRasterIO(scene.r_band, GF_Read,
                   startx, starty, rsizex, rsizey,
                   scene.r_texture + deltax + (TEXTURE_BUFFER_SIZE*deltay),
                   wsizex, wsizey,
                   GDT_UInt16, 0, TEXTURE_BUFFER_SIZE*sizeof(uint16_t))) {
    fprintf(stderr, ANSI_COLOR_RED "GDALRasterIO problem (red band)" ANSI_COLOR_RESET "\n");
    exit(-1);
  }
  if (GDALRasterIO(scene.g_band, GF_Read,
                   startx, starty, rsizex, rsizey,
                   scene.g_texture + deltax + (TEXTURE_BUFFER_SIZE*deltay),
                   wsizex, wsizey,
                   GDT_UInt16, 0, TEXTURE_BUFFER_SIZE*sizeof(uint16_t))) {
    fprintf(stderr, ANSI_COLOR_RED "GDALRasterIO problem (green band)" ANSI_COLOR_RESET "\n");
    exit(-1);
  }
  if (GDALRasterIO(scene.b_band, GF_Read,
                   startx, starty, rsizex, rsizey,
                   scene.b_texture + deltax + (TEXTURE_BUFFER_SIZE*deltay),
                   wsizex, wsizey,
                   GDT_UInt16, 0, TEXTURE_BUFFER_SIZE*sizeof(uint16_t))) {
    fprintf(stderr, ANSI_COLOR_RED "GDALRasterIO problem (blue band)" ANSI_COLOR_RESET "\n");
    exit(-1);
  }

  return 1;
}

void zxy(int fd, int z, int _x, int _y, int verbose)
{
  if (z < 7)
    zxy_exact(fd, z, _x, _y, verbose);
  else
    zxy_approx(fd, z, _x, _y, verbose);
}

void zxy_exact(int fd, int z, int _x, int _y, int verbose)
{
  double xmin = DBL_MAX, ymin = DBL_MAX;
  double xmax = DBL_MIN, ymax = DBL_MIN;

  if (verbose)
    fprintf(stderr, ANSI_COLOR_YELLOW "start: z=%d x=%d, y=%d pid=%d" ANSI_COLOR_RESET "\n", z, _x, _y, getpid());

  /*
    TMS to Pseudo Web Mercator
    Source: https://wiki.openstreetmap.org/wiki/Slippy_map_tilenames
  */
  for (int j = 0; j < TILE_SIZE; ++j) {
    for (int i = 0; i < TILE_SIZE; ++i) {
      int index = (i + j*TILE_SIZE)<<1;
      scene.patch[index+0] = _x + (i/((double)TILE_SIZE));           // tile space
      scene.patch[index+0] /= pow(2.0, z);                           // 0-1 scaled, translated Web Mercator
      scene.patch[index+0] = (2*scene.patch[index+0] - 1) * M_PI * RADIUS; // Web Mercator
      scene.patch[index+1] = _y + (j/((double)TILE_SIZE));
      scene.patch[index+1] /= pow(2.0, z);
      scene.patch[index+1] = (1 - 2*scene.patch[index+1]) * M_PI * RADIUS;
    }
  }

  /* Web Mercator to world coordinates */
  pj_transform(webmercator_pj, scene.destination_pj,
               TILE_SIZE*TILE_SIZE, 2,
               scene.patch, scene.patch+1, NULL);

  /* world coordinates to image coordinates */
  for (int j = 0; j < TILE_SIZE; ++j) {
    for (int i = 0; i < TILE_SIZE; ++i) {
      int index = (i + j*TILE_SIZE)<<1;
      world_to_image(scene.patch+index);
      xmin = fmin(scene.patch[index+0], xmin);
      xmax = fmax(scene.patch[index+0], xmax);
      ymin = fmin(scene.patch[index+1], ymin);
      ymax = fmax(scene.patch[index+1], ymax);
    }
  }

  xmin = floor(xmin);
  xmax = ceil(xmax);
  ymin = floor(ymin);
  ymax = ceil(ymax);

  /* Read clipped textures from image */
  if (!fetch(xmin, xmax, ymin, ymax, verbose)) goto done;

  /* Sample from the texture */
  for (int j = 0; j < TILE_SIZE; ++j) {
    for (int i = 0; i < TILE_SIZE; ++i) {
      double _u, _v;
      uint8_t byte = 0;
      int index = (i + j*TILE_SIZE)<<1;
      _u = scene.patch[index+0];
      _v = scene.patch[index+1];
      if (!isnan(_u) && !isnan(_v)) {
        int u = (int)(((_u - xmin)/(xmax-xmin)) * TEXTURE_BUFFER_SIZE);
        int v = (int)(((_v - ymin)/(ymax-ymin)) * TEXTURE_BUFFER_SIZE);
        byte |= scene.tile[4*i + 4*j*TILE_SIZE + 0] = sigmoidal(scene.r_texture[u + v*TEXTURE_BUFFER_SIZE]);
        byte |= scene.tile[4*i + 4*j*TILE_SIZE + 1] = sigmoidal(scene.g_texture[u + v*TEXTURE_BUFFER_SIZE]);
        byte |= scene.tile[4*i + 4*j*TILE_SIZE + 2] = sigmoidal(scene.b_texture[u + v*TEXTURE_BUFFER_SIZE]);
        scene.tile[4*i + 4*j*TILE_SIZE + 3] = (byte ? -1 : 0);
      }
    }
  }

 done:
  write_png(fd, scene.tile, TILE_SIZE, TILE_SIZE, 0);

  if (verbose)
    fprintf(stderr, ANSI_COLOR_YELLOW "finish: z=%d x=%d, y=%d" ANSI_COLOR_RESET "\n", z, _x, _y);
}

void zxy_approx(int fd, int z, int _x, int _y, int verbose)
{
  double xmin = DBL_MAX, ymin = DBL_MAX;
  double xmax = DBL_MIN, ymax = DBL_MIN;

  if (verbose)
    fprintf(stderr, ANSI_COLOR_YELLOW "start: z=%d x=%d, y=%d pid=%d" ANSI_COLOR_RESET "\n", z, _x, _y, getpid());

  /*
    TMS to Pseudo Web Mercator
    Source: https://wiki.openstreetmap.org/wiki/Slippy_map_tilenames
  */
  for (int i = 0; i < (TILE_SIZE<<1); i+=2) {
    // top, bottom longitudes
    scene.top[i+0] = _x + (i/(TILE_SIZE*2.0));     // tile space
    scene.top[i+0] /= pow(2.0, z);                 // 0-1 scaled, translated Web Mercator
    scene.top[i+0] = (2*scene.top[i+0] - 1) * M_PI;      // Web Mercator in radians
    scene.top[i+0] = scene.bot[i+0] = scene.top[i+0] * RADIUS; // Web Mercator in radians*radius

    // top, bottom latitudes
    scene.top[i+1] = (1 - 2*((_y+0) / pow(2.0, z))) * M_PI * RADIUS;
    scene.bot[i+1] = (1 - 2*((_y+1) / pow(2.0, z))) * M_PI * RADIUS;

    // left, right longitudes
    scene.left[i+0]  = (2*((_x+0) / pow(2.0, z)) - 1) * M_PI * RADIUS;
    scene.right[i+0] = (2*((_x+1) / pow(2.0, z)) - 1) * M_PI * RADIUS;

    // left, right latitudes
    scene.left[i+1] = _y + (i/(TILE_SIZE*2.0));
    scene.left[i+1] /= pow(2.0, z);
    scene.left[i+1] = scene.right[i+1] = (1 - 2*scene.left[i+1]) * M_PI * RADIUS;
  }

  /* Web Mercator to world coordinates */
  pj_transform(webmercator_pj, scene.destination_pj, TILE_SIZE, 2, scene.top,   scene.top+1, NULL);
  pj_transform(webmercator_pj, scene.destination_pj, TILE_SIZE, 2, scene.bot,   scene.bot+1, NULL);
  pj_transform(webmercator_pj, scene.destination_pj, TILE_SIZE, 2, scene.left,  scene.left+1, NULL);
  pj_transform(webmercator_pj, scene.destination_pj, TILE_SIZE, 2, scene.right, scene.right+1, NULL);

  /* world coordinates to image coordinates */
  for (int i = 0; i < (TILE_SIZE<<1); i+=2) {
    world_to_image(scene.top+i);
    world_to_image(scene.bot+i);
    world_to_image(scene.left+i);
    world_to_image(scene.right+i);
    xmin = fmin(scene.right[i], fmin(scene.left[i], fmin(scene.bot[i], fmin(scene.top[i], xmin))));
    xmax = fmax(scene.right[i], fmax(scene.left[i], fmax(scene.bot[i], fmax(scene.top[i], xmax))));
    ymin = fmin(scene.right[i+1], fmin(scene.left[i+1], fmin(scene.bot[i+1], fmin(scene.top[i+1], ymin))));
    ymax = fmax(scene.right[i+1], fmax(scene.left[i+1], fmax(scene.bot[i+1], fmax(scene.top[i+1], ymax))));
  }
  xmin = floor(xmin);
  xmax = ceil(xmax);
  ymin = floor(ymin);
  ymax = ceil(ymax);

  /* Read clipped textures from image */
  if (!fetch(xmin, xmax, ymin, ymax, verbose)) goto done;

  /* Sample from the texture */
  for (int j = 0; j < TILE_SIZE; ++j) {
    for (int i = 0; i < TILE_SIZE; ++i) {
      double _u, _v;
      uint8_t byte = 0;
      _u = scene.top[2*i]*((double)j/TILE_SIZE) + scene.bot[2*i]*(1-((double)j/TILE_SIZE));
      _v = scene.left[2*j+1]*((double)i/TILE_SIZE) + scene.right[2*j+1]*(1-((double)i/TILE_SIZE));
      if (!isnan(_u) && !isnan(_v)) {
        int u = (int)(((_u - xmin)/(xmax-xmin)) * TEXTURE_BUFFER_SIZE);
        int v = (int)(((_v - ymin)/(ymax-ymin)) * TEXTURE_BUFFER_SIZE);
        byte |= scene.tile[4*i + 4*j*TILE_SIZE + 0] = sigmoidal(scene.r_texture[u + v*TEXTURE_BUFFER_SIZE]);
        byte |= scene.tile[4*i + 4*j*TILE_SIZE + 1] = sigmoidal(scene.g_texture[u + v*TEXTURE_BUFFER_SIZE]);
        byte |= scene.tile[4*i + 4*j*TILE_SIZE + 2] = sigmoidal(scene.b_texture[u + v*TEXTURE_BUFFER_SIZE]);
        scene.tile[4*i + 4*j*TILE_SIZE + 3] = (byte ? -1 : 0);
      }
    }
  }

 done:
  write_png(fd, scene.tile, TILE_SIZE, TILE_SIZE, 0);

  if (verbose)
    fprintf(stderr, ANSI_COLOR_YELLOW "finish: z=%d x=%d, y=%d" ANSI_COLOR_RESET "\n", z, _x, _y);
}

uint8_t sigmoidal(uint16_t _u) {
  if (!_u) return 0;

  double u = ((double)_u) / 23130.235294118;
  double beta = 10, alpha = 0.50;
  double numer = 1/(1+exp(beta*(alpha-u))) - 1/(1+exp(beta));
  double denom = 1/(1+exp(beta*(alpha-1))) - 1/(1+exp(beta*alpha));
  double gu = fmax(0.0, fmin(1.0, numer / denom));
  return ((1<<8)-1)*gu;
}

void world_to_image(double * xy)
{
  double world_x = xy[0], world_y = xy[1];

  // Source: http://www.gdal.org/classGDALDataset.html#a5101119705f5fa2bc1344ab26f66fd1d
  xy[0] = (-world_x*scene.transform[5] + world_y*scene.transform[2] + scene.transform[0]*scene.transform[5] - scene.transform[2]*scene.transform[3])/(scene.transform[2]*scene.transform[4] - scene.transform[1]*scene.transform[5]);
  xy[1] = (world_x*scene.transform[4] - world_y*scene.transform[1] + scene.transform[0]*scene.transform[4] + scene.transform[1]*scene.transform[3])/(scene.transform[2]*scene.transform[4] - scene.transform[1]*scene.transform[5]);
}
