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
int scene_count = 3;
uint8_t tile[TILE_SIZE * TILE_SIZE * 4]; // RGBA ergo 4

void load_scene(landsat_scene * s, int verbose);
int fetch(landsat_scene * s, int verbose);
void zxy_approx(int z, int _x, int _y, landsat_scene * s, int verbose);
void zxy_exact(int z, int _x, int _y, landsat_scene * s, int verbose);
uint8_t sigmoidal(uint16_t _u);
void world_to_image(double * xy, landsat_scene * s);

#define LARGER(a,b) (a = a > (b) ? a : (b))

void load(int verbose)
{
  GDALAllRegister();

  scene = calloc(sizeof(landsat_scene), scene_count);
  /* (scene + 0)->filename = "/tmp/LC08_L1TP_139043_20170304_20170316_01_T1_B%d.TIF"; */
  /* (scene + 1)->filename = "/tmp/LC08_L1TP_139044_20170304_20170316_01_T1_B%d.TIF"; */
  /* (scene + 2)->filename = "/tmp/LC08_L1TP_139045_20170304_20170316_01_T1_B%d.TIF"; */
  /* (scene + 0)->filename = "/vsicurl/https://s3-us-west-2.amazonaws.com/landsat-pds/c1/L8/139/044/LC08_L1TP_139044_20170304_20170316_01_T1/LC08_L1TP_139044_20170304_20170316_01_T1_B%d.TIF"; */
  /* (scene + 1)->filename = "/vsicurl/https://s3-us-west-2.amazonaws.com/landsat-pds/c1/L8/139/045/LC08_L1TP_139045_20170304_20170316_01_T1/LC08_L1TP_139045_20170304_20170316_01_T1_B%d.TIF"; */
  /* (scene + 2)->filename = "/vsicurl/https://s3-us-west-2.amazonaws.com/landsat-pds/c1/L8/139/043/LC08_L1TP_139043_20170304_20170316_01_T1/LC08_L1TP_139043_20170304_20170316_01_T1_B%d.TIF"; */
  /* //  docker run --rm -it --name my-apache-app -p 8080:80 -v /tmp/:/usr/local/apache2/htdocs/:ro httpd:2.4 */
  /* (scene + 0)->filename = "/vsicurl/http://localhost:8080/LC08_L1TP_139043_20170304_20170316_01_T1_B%d.TIF"; */
  /* (scene + 1)->filename = "/vsicurl/http://localhost:8080/LC08_L1TP_139044_20170304_20170316_01_T1_B%d.TIF"; */
  /* (scene + 2)->filename = "/vsicurl/http://localhost:8080/LC08_L1TP_139045_20170304_20170316_01_T1_B%d.TIF"; */
  (scene + 0)->filename = "/home/ec2-user/mnt/c1/L8/139/044/LC08_L1TP_139044_20170304_20170316_01_T1/LC08_L1TP_139044_20170304_20170316_01_T1_B%d.TIF";
  (scene + 1)->filename = "/home/ec2-user/mnt/c1/L8/139/045/LC08_L1TP_139045_20170304_20170316_01_T1/LC08_L1TP_139045_20170304_20170316_01_T1_B%d.TIF";
  (scene + 2)->filename = "/home/ec2-user/mnt/c1/L8/139/043/LC08_L1TP_139043_20170304_20170316_01_T1/LC08_L1TP_139043_20170304_20170316_01_T1_B%d.TIF";

  webmercator_pj = pj_init_plus(webmercator);
  for (int i = 0; i < scene_count; ++i)
    load_scene(scene + i, verbose);
}

void load_scene(landsat_scene * s, int verbose)
{
  char * wkt = NULL, * dstProj4 = NULL;
  OGRSpatialReferenceH srs = NULL;
  char r_filename[(1<<8)];
  char g_filename[(1<<8)];
  char b_filename[(1<<8)];

  if (verbose)
    fprintf(stderr, ANSI_COLOR_YELLOW "%s" ANSI_COLOR_RESET "\n", s->filename);

  sprintf(r_filename, s->filename, 4);
  sprintf(g_filename, s->filename, 3);
  sprintf(b_filename, s->filename, 2);

  /* Datasets and bands */
  s->r_dataset = GDALOpen(r_filename, GA_ReadOnly);
  if(s->r_dataset == NULL) {
    fprintf(stderr, ANSI_COLOR_RED "GDALOpen problem (red band band)" ANSI_COLOR_RESET "\n");
    exit(-1);
  }
  s->g_dataset = GDALOpen(g_filename, GA_ReadOnly);
  if(s->g_dataset == NULL) {
    fprintf(stderr, ANSI_COLOR_RED "GDALOpen problem (green band)" ANSI_COLOR_RESET "\n");
    exit(-1);
  }
  s->b_dataset = GDALOpen(b_filename, GA_ReadOnly);
  if(s->b_dataset == NULL) {
    fprintf(stderr, ANSI_COLOR_RED "GDALOpen problem (blue band)" ANSI_COLOR_RESET "\n");
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
    fprintf(stderr, ANSI_COLOR_GREEN "WKT=%s" ANSI_COLOR_RESET "\n", wkt);
  OSRImportFromWkt(srs, &wkt);
  OSRExportToProj4(srs, &dstProj4);
  if (verbose)
    fprintf(stderr, ANSI_COLOR_GREEN "Proj4=%s" ANSI_COLOR_RESET "\n", dstProj4);

  s->destination_pj = pj_init_plus(dstProj4);

  /* Transform.  This is from the red band, but is assumed to be valid
     for all of the bands. */
  GDALGetGeoTransform(s->r_dataset, s->transform);

  CPLFree(dstProj4);
  OSRRelease(srs);

  if (verbose)
    fprintf(stderr, ANSI_COLOR_BLUE "pid=%d" ANSI_COLOR_RESET "\n", getpid());
}

int fetch(landsat_scene * s, int verbose)
{
  double xmin = s->xmin, xmax = s->xmax, ymin = s->ymin, ymax = s->ymax;
  int startx = (int)xmin, starty = (int)ymin;
  unsigned int rsizex = (int)(xmax-xmin), rsizey = (int)(ymax-ymin);
  unsigned int wsizex = TEXTURE_BUFFER_SIZE, wsizey = TEXTURE_BUFFER_SIZE;
  unsigned int deltax = 0, deltay = 0;

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
            ANSI_COLOR_CYAN "startx=%d rsize=%d wsize=%d delta=%d pid=%d"
            ANSI_COLOR_RESET "\n",
            startx, rsizex, wsizex, deltax, getpid());
    fprintf(stderr,
            ANSI_COLOR_CYAN "starty=%d rsize=%d wsize=%d delta=%d pid=%d"
            ANSI_COLOR_RESET "\n",
            starty, rsizey, wsizey, deltay, getpid());
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
  int exact = (z < 7) ? 1 : 0;

  #pragma omp parallel for schedule(dynamic, 1)
  for (int i = -1; i < scene_count; ++i) {
    if (i == -1)
      memset(tile, 0, sizeof(tile));
    else if (exact)
      zxy_exact(z, _x, _y, scene + i, verbose);
    else if (!exact)
      zxy_approx(z, _x, _y, scene + i, verbose);
  }

  // Sample from textures to produce tile
  for (int i = 0; i < scene_count; ++i)
    {

      landsat_scene * s = scene + i;

      if (!s->dirty) continue;

      double * patch = s->coordinates.patch;
      double * top   = s->coordinates.periphery.top;
      double * bot   = s->coordinates.periphery.bot;
      double * left  = s->coordinates.periphery.left;
      double * right = s->coordinates.periphery.right;
      double xmin = s->xmin, xmax = s->xmax, ymin = s->ymin, ymax = s->ymax;

      for (int j = 0; j < TILE_SIZE; ++j) {
        for (int i = 0; i < TILE_SIZE; ++i) {
          double _u, _v;
          uint8_t red, byte = 0;
          int index = (i + j*TILE_SIZE)<<1;

          if (z < 7) {
            _u = patch[index+0];
            _v = patch[index+1];
          }
          else {
            _u = top[2*i]*((double)j/TILE_SIZE) + bot[2*i]*(1-((double)j/TILE_SIZE));
            _v = left[2*j+1]*((double)i/TILE_SIZE) + right[2*j+1]*(1-((double)i/TILE_SIZE));
          }

          if (!isnan(_u) && !isnan(_v)) {
            int u = (int)(((_u - xmin)/(xmax-xmin)) * TEXTURE_BUFFER_SIZE);
            int v = (int)(((_v - ymin)/(ymax-ymin)) * TEXTURE_BUFFER_SIZE);
            int index = 4*i + 4*j*TILE_SIZE;

            byte |= red = sigmoidal(s->r_texture[u + v*TEXTURE_BUFFER_SIZE]);
            if (tile[index + 3] == 0 || tile[index + 0] < red) { // write into empty pixels
              tile[index + 0] = red;
              byte |= tile[index + 1] = sigmoidal(s->g_texture[u + v*TEXTURE_BUFFER_SIZE]);
              byte |= tile[index + 2] = sigmoidal(s->b_texture[u + v*TEXTURE_BUFFER_SIZE]);
              tile[index + 3] = (byte ? -1 : 0);
            }
          }
        }
      }
    }

  write_png(fd, tile, TILE_SIZE, TILE_SIZE, 0);
}

void zxy_exact(int z, int _x, int _y, landsat_scene * s, int verbose)
{
  double * patch = s->coordinates.patch;
  double xmin = DBL_MAX, ymin = DBL_MAX;
  double xmax = DBL_MIN, ymax = DBL_MIN;

  if (verbose)
    fprintf(stderr,
            ANSI_COLOR_YELLOW "z=%d x=%d y=%d pid=%d" ANSI_COLOR_RESET "\n",
            z, _x, _y, getpid());

  /*
    TMS to Pseudo Web Mercator
    Source: https://wiki.openstreetmap.org/wiki/Slippy_map_tilenames
  */
  for (int j = 0; j < TILE_SIZE; ++j) {
    for (int i = 0; i < TILE_SIZE; ++i) {
      int index = (i + j*TILE_SIZE)<<1;
      patch[index+0] = _x + (i/((double)TILE_SIZE));           // tile space
      patch[index+0] /= pow(2.0, z);                           // 0-1 scaled, translated Web Mercator
      patch[index+0] = (2*patch[index+0] - 1) * M_PI * RADIUS; // Web Mercator
      patch[index+1] = _y + (j/((double)TILE_SIZE));
      patch[index+1] /= pow(2.0, z);
      patch[index+1] = (1 - 2*patch[index+1]) * M_PI * RADIUS;
    }
  }

  /* Web Mercator to world coordinates */
  pj_transform(webmercator_pj, s->destination_pj,
               TILE_SIZE*TILE_SIZE, 2,
               patch, patch+1, NULL);

  /* world coordinates to image coordinates */
  for (int j = 0; j < TILE_SIZE; ++j) {
    for (int i = 0; i < TILE_SIZE; ++i) {
      int index = (i + j*TILE_SIZE)<<1;
      world_to_image(patch+index, s);
      xmin = fmin(patch[index+0], xmin);
      xmax = fmax(patch[index+0], xmax);
      ymin = fmin(patch[index+1], ymin);
      ymax = fmax(patch[index+1], ymax);
    }
  }

  s->xmin = xmin = floor(xmin);
  s->xmax = xmax = ceil(xmax);
  s->ymin = ymin = floor(ymin);
  s->ymax = ymax = ceil(ymax);

  s->dirty = fetch(s, verbose);
}

void zxy_approx(int z, int _x, int _y, landsat_scene * s, int verbose)
{
  double * top   = s->coordinates.periphery.top;
  double * bot   = s->coordinates.periphery.bot;
  double * left  = s->coordinates.periphery.left;
  double * right = s->coordinates.periphery.right;
  double xmin = DBL_MAX, ymin = DBL_MAX;
  double xmax = DBL_MIN, ymax = DBL_MIN;

  if (verbose)
    fprintf(stderr,
            ANSI_COLOR_YELLOW "z=%d x=%d y=%d pid=%d" ANSI_COLOR_RESET "\n",
            z, _x, _y, getpid());

  /*
    TMS to Pseudo Web Mercator
    Source: https://wiki.openstreetmap.org/wiki/Slippy_map_tilenames
  */
  for (int i = 0; i < (TILE_SIZE<<1); i+=2) {
    // top, bottom longitudes
    top[i+0] = _x + (i/(TILE_SIZE*2.0));           // tile space
    top[i+0] /= pow(2.0, z);                       // 0-1 scaled, translated Web Mercator
    top[i+0] = (2*top[i+0] - 1) * M_PI;         // Web Mercator in radians
    top[i+0] = bot[i+0] = top[i+0] * RADIUS; // Web Mercator in radians*radius

    // top, bottom latitudes
    top[i+1] = (1 - 2*((_y+0) / pow(2.0, z))) * M_PI * RADIUS;
    bot[i+1] = (1 - 2*((_y+1) / pow(2.0, z))) * M_PI * RADIUS;

    // left, right longitudes
    left[i+0]  = (2*((_x+0) / pow(2.0, z)) - 1) * M_PI * RADIUS;
    right[i+0] = (2*((_x+1) / pow(2.0, z)) - 1) * M_PI * RADIUS;

    // left, right latitudes
    left[i+1] = _y + (i/(TILE_SIZE*2.0));
    left[i+1] /= pow(2.0, z);
    left[i+1] = right[i+1] = (1 - 2*left[i+1]) * M_PI * RADIUS;
  }

  /* Web Mercator to world coordinates */
  pj_transform(webmercator_pj, s->destination_pj, TILE_SIZE, 2, top,   top+1, NULL);
  pj_transform(webmercator_pj, s->destination_pj, TILE_SIZE, 2, bot,   bot+1, NULL);
  pj_transform(webmercator_pj, s->destination_pj, TILE_SIZE, 2, left,  left+1, NULL);
  pj_transform(webmercator_pj, s->destination_pj, TILE_SIZE, 2, right, right+1, NULL);

  /* world coordinates to image coordinates */
  for (int i = 0; i < (TILE_SIZE<<1); i+=2) {
    world_to_image(top+i, s);
    world_to_image(bot+i, s);
    world_to_image(left+i, s);
    world_to_image(right+i, s);
    xmin = fmin(right[i], fmin(left[i], fmin(bot[i], fmin(top[i], xmin))));
    xmax = fmax(right[i], fmax(left[i], fmax(bot[i], fmax(top[i], xmax))));
    ymin = fmin(right[i+1], fmin(left[i+1], fmin(bot[i+1], fmin(top[i+1], ymin))));
    ymax = fmax(right[i+1], fmax(left[i+1], fmax(bot[i+1], fmax(top[i+1], ymax))));
  }
  s->xmin = xmin = floor(xmin);
  s->xmax = xmax = ceil(xmax);
  s->ymin = ymin = floor(ymin);
  s->ymax = ymax = ceil(ymax);

  s->dirty = fetch(s, verbose);
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
