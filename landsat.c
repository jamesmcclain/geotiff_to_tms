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
landsat_scene * scene = NULL;

void load_scene(landsat_scene * s, int verbose);
int fetch(double xmin, double xmax, double ymin, double ymax, landsat_scene * s, int verbose);
void zxy_approx(int fd, int z, int _x, int _y, landsat_scene * s, int verbose);
void zxy_exact(int fd, int z, int _x, int _y, landsat_scene * s, int verbose);
uint8_t sigmoidal(uint16_t _u);
void world_to_image(double * xy, landsat_scene * s);


void load(int verbose)
{
  GDALAllRegister();

  scene = calloc(sizeof(landsat_scene), 1);
  (scene + 0)->filename = "/tmp/LC08_L1TP_139045_20170304_20170316_01_T1_B%d.TIF";
  /* (scene + 0)->filename = "/vsicurl/https://landsat-pds.s3.amazonaws.com/c1/L8/139/045/LC08_L1TP_139045_20170304_20170316_01_T1/LC08_L1TP_139045_20170304_20170316_01_T1_B%d.TIF"; */

  webmercator_pj = pj_init_plus(webmercator);
  load_scene(scene + 0, verbose);
}

void load_scene(landsat_scene * s, int verbose)
{
  char * wkt = NULL, * dstProj4 = NULL;
  OGRSpatialReferenceH srs = NULL;
  char r_filename[(1<<8)];
  char g_filename[(1<<8)];
  char b_filename[(1<<8)];

  sprintf(r_filename, s->filename, 4);
  sprintf(g_filename, s->filename, 3);
  sprintf(b_filename, s->filename, 2);

  /* Datasets and bands */
  s->r_dataset = GDALOpen(r_filename, GA_ReadOnly);
  if(s->r_dataset == NULL) {
    fprintf(stderr, ANSI_COLOR_RED "GDALOpen problem (red band)" ANSI_COLOR_RESET "\n");
    exit(-1);
  }
  s->g_dataset = GDALOpen(g_filename, GA_ReadOnly);
  if(s->g_dataset == NULL) {
    fprintf(stderr, ANSI_COLOR_RED "GDALOpen problem (green)" ANSI_COLOR_RESET "\n");
    exit(-1);
  }
  s->b_dataset = GDALOpen(b_filename, GA_ReadOnly);
  if(s->b_dataset == NULL) {
    fprintf(stderr, ANSI_COLOR_RED "GDALOpen problem (blue)" ANSI_COLOR_RESET "\n");
    exit(-1);
  }
  s->r_band = GDALGetRasterBand(s->r_dataset, 1);
  s->g_band = GDALGetRasterBand(s->g_dataset, 1);
  s->b_band = GDALGetRasterBand(s->b_dataset, 1);
  s->width  = GDALGetRasterXSize(s->r_dataset);
  s->height = GDALGetRasterYSize(s->r_dataset);

  /* SRS.  This is from the red band, but is assumed to be valid for
     all of the bands. */
  srs = OSRNewSpatialReference(NULL);
  wkt = calloc(STRING_BUFFER_SIZE, sizeof(char)); // No memory leak, this is freed from within `OSRImportFromWkt`!
  strncpy(wkt, GDALGetProjectionRef(s->r_dataset), STRING_BUFFER_SIZE);
  if (verbose)
    fprintf(stderr, ANSI_COLOR_GREEN "WKT: " ANSI_COLOR_CYAN "%s" ANSI_COLOR_RESET "\n", wkt);
  OSRImportFromWkt(srs, &wkt);
  OSRExportToProj4(srs, &dstProj4);
  if (verbose)
    fprintf(stderr, ANSI_COLOR_GREEN "Proj4: " ANSI_COLOR_CYAN "%s" ANSI_COLOR_RESET "\n", dstProj4);

  s->destination_pj = pj_init_plus(dstProj4);

  /* Transform.  This is from the red band, but is assumed to be valid
     for all of the bands. */
  GDALGetGeoTransform(s->r_dataset, s->transform);

  CPLFree(dstProj4);
  OSRRelease(srs);
}

int fetch(double xmin, double xmax, double ymin, double ymax, landsat_scene * s, int verbose)
{
  int startx = (int)xmin, starty = (int)ymin;
  int rsizex = (int)(xmax-xmin), rsizey = (int)(ymax-ymin);
  int wsizex = TEXTURE_BUFFER_SIZE, wsizey = TEXTURE_BUFFER_SIZE;
  int deltax = 0, deltay = 0;

  memset(s->tile, 0, sizeof(s->tile));
  memset(s->r_texture, 0, sizeof(s->r_texture));
  memset(s->g_texture, 0, sizeof(s->g_texture));
  memset(s->b_texture, 0, sizeof(s->b_texture));

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
  if (startx + rsizex > s->width) {
    rsizex = s->width - startx;
    rsizex = rsizex > 0? rsizex : 0;
  }
  if (starty + rsizey > s->height) {
    rsizey = s->height - starty;
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

  if (GDALRasterIO(s->r_band, GF_Read,
                   startx, starty, rsizex, rsizey,
                   s->r_texture + deltax + (TEXTURE_BUFFER_SIZE*deltay),
                   wsizex, wsizey,
                   GDT_UInt16, 0, TEXTURE_BUFFER_SIZE*sizeof(uint16_t))) {
    fprintf(stderr, ANSI_COLOR_RED "GDALRasterIO problem (red band)" ANSI_COLOR_RESET "\n");
    exit(-1);
  }
  if (GDALRasterIO(s->g_band, GF_Read,
                   startx, starty, rsizex, rsizey,
                   s->g_texture + deltax + (TEXTURE_BUFFER_SIZE*deltay),
                   wsizex, wsizey,
                   GDT_UInt16, 0, TEXTURE_BUFFER_SIZE*sizeof(uint16_t))) {
    fprintf(stderr, ANSI_COLOR_RED "GDALRasterIO problem (green band)" ANSI_COLOR_RESET "\n");
    exit(-1);
  }
  if (GDALRasterIO(s->b_band, GF_Read,
                   startx, starty, rsizex, rsizey,
                   s->b_texture + deltax + (TEXTURE_BUFFER_SIZE*deltay),
                   wsizex, wsizey,
                   GDT_UInt16, 0, TEXTURE_BUFFER_SIZE*sizeof(uint16_t))) {
    fprintf(stderr, ANSI_COLOR_RED "GDALRasterIO problem (blue band)" ANSI_COLOR_RESET "\n");
    exit(-1);
  }

  return 1;
}

void zxy(int fd, int z, int _x, int _y, int verbose)
{
  landsat_scene * s = scene;

  if (z < 7)
    zxy_exact(fd, z, _x, _y, s, verbose);
  else
    zxy_approx(fd, z, _x, _y, s, verbose);
}

void zxy_exact(int fd, int z, int _x, int _y, landsat_scene * s, int verbose)
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
      s->patch[index+0] = _x + (i/((double)TILE_SIZE));              // tile space
      s->patch[index+0] /= pow(2.0, z);                              // 0-1 scaled, translated Web Mercator
      s->patch[index+0] = (2*s->patch[index+0] - 1) * M_PI * RADIUS; // Web Mercator
      s->patch[index+1] = _y + (j/((double)TILE_SIZE));
      s->patch[index+1] /= pow(2.0, z);
      s->patch[index+1] = (1 - 2*s->patch[index+1]) * M_PI * RADIUS;
    }
  }

  /* Web Mercator to world coordinates */
  pj_transform(webmercator_pj, s->destination_pj,
               TILE_SIZE*TILE_SIZE, 2,
               s->patch, s->patch+1, NULL);

  /* world coordinates to image coordinates */
  for (int j = 0; j < TILE_SIZE; ++j) {
    for (int i = 0; i < TILE_SIZE; ++i) {
      int index = (i + j*TILE_SIZE)<<1;
      world_to_image(s->patch+index, s);
      xmin = fmin(s->patch[index+0], xmin);
      xmax = fmax(s->patch[index+0], xmax);
      ymin = fmin(s->patch[index+1], ymin);
      ymax = fmax(s->patch[index+1], ymax);
    }
  }

  xmin = floor(xmin);
  xmax = ceil(xmax);
  ymin = floor(ymin);
  ymax = ceil(ymax);

  /* Read clipped textures from image */
  if (!fetch(xmin, xmax, ymin, ymax, s, verbose)) goto done;

  /* Sample from the textures */
  for (int j = 0; j < TILE_SIZE; ++j) {
    for (int i = 0; i < TILE_SIZE; ++i) {
      double _u, _v;
      uint8_t byte = 0;
      int index = (i + j*TILE_SIZE)<<1;
      _u = s->patch[index+0];
      _v = s->patch[index+1];
      if (!isnan(_u) && !isnan(_v)) {
        int u = (int)(((_u - xmin)/(xmax-xmin)) * TEXTURE_BUFFER_SIZE);
        int v = (int)(((_v - ymin)/(ymax-ymin)) * TEXTURE_BUFFER_SIZE);
        byte |= s->tile[4*i + 4*j*TILE_SIZE + 0] = sigmoidal(s->r_texture[u + v*TEXTURE_BUFFER_SIZE]);
        byte |= s->tile[4*i + 4*j*TILE_SIZE + 1] = sigmoidal(s->g_texture[u + v*TEXTURE_BUFFER_SIZE]);
        byte |= s->tile[4*i + 4*j*TILE_SIZE + 2] = sigmoidal(s->b_texture[u + v*TEXTURE_BUFFER_SIZE]);
        s->tile[4*i + 4*j*TILE_SIZE + 3] = (byte ? -1 : 0);
      }
    }
  }

 done:
  write_png(fd, s->tile, TILE_SIZE, TILE_SIZE, 0);

  if (verbose)
    fprintf(stderr, ANSI_COLOR_YELLOW "finish: z=%d x=%d, y=%d" ANSI_COLOR_RESET "\n", z, _x, _y);
}

void zxy_approx(int fd, int z, int _x, int _y, landsat_scene * s, int verbose)
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
    s->top[i+0] = _x + (i/(TILE_SIZE*2.0));           // tile space
    s->top[i+0] /= pow(2.0, z);                       // 0-1 scaled, translated Web Mercator
    s->top[i+0] = (2*s->top[i+0] - 1) * M_PI;         // Web Mercator in radians
    s->top[i+0] = s->bot[i+0] = s->top[i+0] * RADIUS; // Web Mercator in radians*radius

    // top, bottom latitudes
    s->top[i+1] = (1 - 2*((_y+0) / pow(2.0, z))) * M_PI * RADIUS;
    s->bot[i+1] = (1 - 2*((_y+1) / pow(2.0, z))) * M_PI * RADIUS;

    // left, right longitudes
    s->left[i+0]  = (2*((_x+0) / pow(2.0, z)) - 1) * M_PI * RADIUS;
    s->right[i+0] = (2*((_x+1) / pow(2.0, z)) - 1) * M_PI * RADIUS;

    // left, right latitudes
    s->left[i+1] = _y + (i/(TILE_SIZE*2.0));
    s->left[i+1] /= pow(2.0, z);
    s->left[i+1] = s->right[i+1] = (1 - 2*s->left[i+1]) * M_PI * RADIUS;
  }

  /* Web Mercator to world coordinates */
  pj_transform(webmercator_pj, s->destination_pj, TILE_SIZE, 2, s->top,   s->top+1, NULL);
  pj_transform(webmercator_pj, s->destination_pj, TILE_SIZE, 2, s->bot,   s->bot+1, NULL);
  pj_transform(webmercator_pj, s->destination_pj, TILE_SIZE, 2, s->left,  s->left+1, NULL);
  pj_transform(webmercator_pj, s->destination_pj, TILE_SIZE, 2, s->right, s->right+1, NULL);

  /* world coordinates to image coordinates */
  for (int i = 0; i < (TILE_SIZE<<1); i+=2) {
    world_to_image(s->top+i, s);
    world_to_image(s->bot+i, s);
    world_to_image(s->left+i, s);
    world_to_image(s->right+i, s);
    xmin = fmin(s->right[i], fmin(s->left[i], fmin(s->bot[i], fmin(s->top[i], xmin))));
    xmax = fmax(s->right[i], fmax(s->left[i], fmax(s->bot[i], fmax(s->top[i], xmax))));
    ymin = fmin(s->right[i+1], fmin(s->left[i+1], fmin(s->bot[i+1], fmin(s->top[i+1], ymin))));
    ymax = fmax(s->right[i+1], fmax(s->left[i+1], fmax(s->bot[i+1], fmax(s->top[i+1], ymax))));
  }
  xmin = floor(xmin);
  xmax = ceil(xmax);
  ymin = floor(ymin);
  ymax = ceil(ymax);

  /* Read clipped textures from image */
  if (!fetch(xmin, xmax, ymin, ymax, s, verbose)) goto done;

  /* Sample from the textures */
  for (int j = 0; j < TILE_SIZE; ++j) {
    for (int i = 0; i < TILE_SIZE; ++i) {
      double _u, _v;
      uint8_t byte = 0;
      _u = s->top[2*i]*((double)j/TILE_SIZE) + s->bot[2*i]*(1-((double)j/TILE_SIZE));
      _v = s->left[2*j+1]*((double)i/TILE_SIZE) + s->right[2*j+1]*(1-((double)i/TILE_SIZE));
      if (!isnan(_u) && !isnan(_v)) {
        int u = (int)(((_u - xmin)/(xmax-xmin)) * TEXTURE_BUFFER_SIZE);
        int v = (int)(((_v - ymin)/(ymax-ymin)) * TEXTURE_BUFFER_SIZE);
        byte |= s->tile[4*i + 4*j*TILE_SIZE + 0] = sigmoidal(s->r_texture[u + v*TEXTURE_BUFFER_SIZE]);
        byte |= s->tile[4*i + 4*j*TILE_SIZE + 1] = sigmoidal(s->g_texture[u + v*TEXTURE_BUFFER_SIZE]);
        byte |= s->tile[4*i + 4*j*TILE_SIZE + 2] = sigmoidal(s->b_texture[u + v*TEXTURE_BUFFER_SIZE]);
        s->tile[4*i + 4*j*TILE_SIZE + 3] = (byte ? -1 : 0);
      }
    }
  }

 done:
  write_png(fd, s->tile, TILE_SIZE, TILE_SIZE, 0);

  if (verbose)
    fprintf(stderr, ANSI_COLOR_YELLOW "finish: z=%d x=%d, y=%d" ANSI_COLOR_RESET "\n", z, _x, _y);
}

uint8_t sigmoidal(uint16_t _u)
{
  if (!_u) return 0;

  double u = ((double)_u) / 23130.235294118;
  double beta = 10, alpha = 0.50;
  double numer = 1/(1+exp(beta*(alpha-u))) - 1/(1+exp(beta));
  double denom = 1/(1+exp(beta*(alpha-1))) - 1/(1+exp(beta*alpha));
  double gu = fmax(0.0, fmin(1.0, numer / denom));
  return ((1<<8)-1)*gu;
}

void world_to_image(double * xy, landsat_scene * s)
{
  double world_x = xy[0], world_y = xy[1];

  // Source: http://www.gdal.org/classGDALDataset.html#a5101119705f5fa2bc1344ab26f66fd1d
  xy[0] = (-world_x*s->transform[5] + world_y*s->transform[2] + s->transform[0]*s->transform[5] - s->transform[2]*s->transform[3])/(s->transform[2]*s->transform[4] - s->transform[1]*s->transform[5]);
  xy[1] = (world_x*s->transform[4] - world_y*s->transform[1] + s->transform[0]*s->transform[4] + s->transform[1]*s->transform[3])/(s->transform[2]*s->transform[4] - s->transform[1]*s->transform[5]);
}
