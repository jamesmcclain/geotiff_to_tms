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
#include "gdal.h"
#include "cpl_conv.h"
#include "ogr_srs_api.h"
#include "proj_api.h"

#include "load.h"
#include "ansi.h"
#include "pngwrite.h"


GDALDatasetH r_dataset = NULL;
GDALDatasetH g_dataset = NULL;
GDALDatasetH b_dataset = NULL;
GDALRasterBandH r_band = NULL;
GDALRasterBandH g_band = NULL;
GDALRasterBandH b_band = NULL;
const char * latlng = "+proj=longlat +datum=WGS84 +no_defs ";
const char * filename = "/tmp/LC08_L1TP_139045_20170304_20170316_01_T1_B%d.TIF";
/* const char * filename = "/vsicurl/https://landsat-pds.s3.amazonaws.com/c1/L8/139/045/LC08_L1TP_139045_20170304_20170316_01_T1/LC08_L1TP_139045_20170304_20170316_01_T1_B%d.TIF"; */
projPJ src = NULL, dst = NULL;
double t[6];
uint32_t width, height;

void world_to_image(double * xy);

#define STRING_BUFFER_SIZE (1<<10)
#define TILE_SIZE (1<<8)
#define TEXTURE_BUFFER_SIZE ((int)(ceil(sqrt(2)*TILE_SIZE)))


void load()
{
  char * wkt = NULL, * dstProj4 = NULL;
  OGRSpatialReferenceH srs = NULL;
  char r_filename[(1<<8)];
  char g_filename[(1<<8)];
  char b_filename[(1<<8)];

  sprintf(r_filename, filename, 4);
  sprintf(g_filename, filename, 3);
  sprintf(b_filename, filename, 2);

  GDALAllRegister();

  /* Dataset */
  r_dataset = GDALOpen(r_filename, GA_ReadOnly);
  g_dataset = GDALOpen(g_filename, GA_ReadOnly);
  b_dataset = GDALOpen(b_filename, GA_ReadOnly);
  if(r_dataset == NULL) {
    fprintf(stderr, ANSI_COLOR_RED "GDALOpen problem" ANSI_COLOR_RESET "\n");
    exit(-1);
  }
  r_band = GDALGetRasterBand(r_dataset, 1);
  g_band = GDALGetRasterBand(g_dataset, 1);
  b_band = GDALGetRasterBand(b_dataset, 1);
  width = GDALGetRasterXSize(r_dataset);
  height = GDALGetRasterYSize(r_dataset);

  /* SRS */
  srs = OSRNewSpatialReference(NULL);
  wkt = calloc(STRING_BUFFER_SIZE, sizeof(char));
  strncpy(wkt, GDALGetProjectionRef(r_dataset), STRING_BUFFER_SIZE);
  /* fprintf(stderr, ANSI_COLOR_GREEN "WKT: " ANSI_COLOR_CYAN "%s" ANSI_COLOR_RESET "\n", wkt); */
  OSRImportFromWkt(srs, &wkt);
  OSRExportToProj4(srs, &dstProj4);
  /* fprintf(stderr, ANSI_COLOR_GREEN "Proj4: " ANSI_COLOR_CYAN "%s" ANSI_COLOR_RESET "\n", dstProj4); */

  src = pj_init_plus(latlng);
  dst = pj_init_plus(dstProj4);

  /* Transform */
  GDALGetGeoTransform(r_dataset, t);

  CPLFree(dstProj4);
  OSRRelease(srs);
}

void zxy(int fd, int z, int _x, int _y)
{
  double * top;
  double * bot;
  double * left;
  double * right;
  double n = pow(2.0, z);
  double xmin = DBL_MAX, ymin = DBL_MAX;
  double xmax = DBL_MIN, ymax = DBL_MIN;
  uint16_t * r_texture = NULL;
  uint16_t * g_texture = NULL;
  uint16_t * b_texture = NULL;
  uint16_t * tile = NULL;

  top   = calloc((TILE_SIZE<<1), sizeof(double));
  bot   = calloc((TILE_SIZE<<1), sizeof(double));
  left  = calloc((TILE_SIZE<<1), sizeof(double));
  right = calloc((TILE_SIZE<<1), sizeof(double));

  /*
    TMS to longitude, latitude pairs
    Source: https://wiki.openstreetmap.org/wiki/Slippy_map_tilenames
  */
  for (int i = 0; i < (TILE_SIZE<<1); i+=2) {
    top[i]   = _x + (i/(TILE_SIZE*2.0));
    top[i]   = bot[i] = ((2*M_PI*top[i])/n)-M_PI;
    top[i+1] = atan(sinh(M_PI*(1-(2*_y)/n)));
    bot[i+1] = atan(sinh(M_PI*(1-(2*(_y+1))/n)));
    left[i+1] = _y + (i/(TILE_SIZE*2.0));
    left[i+1] = right[i+1] = atan(sinh(M_PI*(1-(2*left[i+1])/n)));
    left[i]   = ((2*M_PI*_x)/n)-M_PI;
    right[i]  = ((2*M_PI*(_x+1))/n)-M_PI;
  }

  /* longitude, latitude pairs to world coordinates */
  pj_transform(src, dst, TILE_SIZE, 2, top, top+1, NULL);
  pj_transform(src, dst, TILE_SIZE, 2, bot, bot+1, NULL);
  pj_transform(src, dst, TILE_SIZE, 2, left, left+1, NULL);
  pj_transform(src, dst, TILE_SIZE, 2, right, right+1, NULL);

  /* world coordinates to image coordinates */
  for (int i = 0; i < (TILE_SIZE<<1); i+=2) {
    world_to_image(top+i);
    world_to_image(bot+i);
    world_to_image(left+i);
    world_to_image(right+i);
    xmin = fmin(right[i], fmin(left[i], fmin(bot[i], fmin(top[i], xmin))));
    xmax = fmax(right[i], fmax(left[i], fmax(bot[i], fmax(top[i], xmax))));
    ymin = fmin(right[i+1], fmin(left[i+1], fmin(bot[i+1], fmin(top[i+1], ymin))));
    ymax = fmax(right[i+1], fmax(left[i+1], fmax(bot[i+1], fmax(top[i+1], ymax))));
  }
  xmin = floor(xmin);
  xmax = ceil(xmax);
  ymin = floor(ymin);
  ymax = ceil(ymax);

  /* Read clipped textures from image */
  {
    int startx = (int)xmin, starty = (int)ymin;
    int rsizex = (int)(xmax-xmin), rsizey = (int)(ymax-ymin);
    int wsizex = TEXTURE_BUFFER_SIZE, wsizey = TEXTURE_BUFFER_SIZE;
    int deltax = 0, deltay = 0;

    r_texture = calloc(TEXTURE_BUFFER_SIZE * TEXTURE_BUFFER_SIZE, sizeof(*r_texture));
    g_texture = calloc(TEXTURE_BUFFER_SIZE * TEXTURE_BUFFER_SIZE, sizeof(*g_texture));
    b_texture = calloc(TEXTURE_BUFFER_SIZE * TEXTURE_BUFFER_SIZE, sizeof(*b_texture));

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
    if (startx + rsizex >= width) {
      rsizex = width - startx - 1;
      rsizex = rsizex > 0? rsizex : 0;
    }
    if (starty + rsizey >= height) {
      rsizey = height - starty - 1;
      rsizey = rsizey > 0? rsizey : 0;
    }
    wsizex = (int)(wsizex * (rsizex/(xmax-xmin)));
    wsizey = (int)(wsizey * (rsizey/(ymax-ymin)));

    if (GDALRasterIO(r_band, GF_Read,
                     startx, starty, rsizex, rsizey,
                     r_texture + deltax + (TEXTURE_BUFFER_SIZE*deltay),
                     wsizex, wsizey,
                     GDT_UInt16, 0, TEXTURE_BUFFER_SIZE*sizeof(uint16_t))) {
      fprintf(stderr, ANSI_COLOR_RED "GDALRasterIO problem (red band)" ANSI_COLOR_RESET "\n");
      exit(-1);
    }
    if (GDALRasterIO(g_band, GF_Read,
                     startx, starty, rsizex, rsizey,
                     g_texture + deltax + (TEXTURE_BUFFER_SIZE*deltay),
                     wsizex, wsizey,
                     GDT_UInt16, 0, TEXTURE_BUFFER_SIZE*sizeof(uint16_t))) {
      fprintf(stderr, ANSI_COLOR_RED "GDALRasterIO problem (green band)" ANSI_COLOR_RESET "\n");
      exit(-1);
    }
    if (GDALRasterIO(b_band, GF_Read,
                     startx, starty, rsizex, rsizey,
                     b_texture + deltax + (TEXTURE_BUFFER_SIZE*deltay),
                     wsizex, wsizey,
                     GDT_UInt16, 0, TEXTURE_BUFFER_SIZE*sizeof(uint16_t))) {
      fprintf(stderr, ANSI_COLOR_RED "GDALRasterIO problem (blue band)" ANSI_COLOR_RESET "\n");
      exit(-1);
    }
  }

  /* Sample from the texture */
  tile = calloc(TILE_SIZE*TILE_SIZE*4, sizeof(*tile));
  for (int j = 0; j < TILE_SIZE; ++j) {
    for (int i = 0; i < TILE_SIZE; ++i) {
      int u, v;
      double _u, _v;
      uint16_t dword = 0;
      _u = top[2*i]*((double)j/TILE_SIZE) + bot[2*i]*(1-((double)j/TILE_SIZE));
      _u = (_u - xmin)/(xmax-xmin);
      _v = left[2*j+1]*((double)i/TILE_SIZE) + right[2*j+1]*(1-((double)i/TILE_SIZE));
      _v = (_v - ymin)/(ymax-ymin);
      u = (int)(_u * TEXTURE_BUFFER_SIZE);
      v = (int)(_v * TEXTURE_BUFFER_SIZE);
      dword |= tile[4*i + 4*j*TILE_SIZE + 0] = htons(r_texture[u + v*TEXTURE_BUFFER_SIZE]);
      dword |= tile[4*i + 4*j*TILE_SIZE + 1] = htons(g_texture[u + v*TEXTURE_BUFFER_SIZE]);
      dword |= tile[4*i + 4*j*TILE_SIZE + 2] = htons(b_texture[u + v*TEXTURE_BUFFER_SIZE]);
      tile[4*i + 4*j*TILE_SIZE + 3] = (dword ? -1 : 0);
    }
  }

  write_png(fd, tile, TILE_SIZE, TILE_SIZE);

  /* Cleanup */
  free(tile);
  free(b_texture);
  free(g_texture);
  free(r_texture);
  free(right);
  free(left);
  free(bot);
  free(top);
}

void world_to_image(double * xy)
{
  double world_x = xy[0], world_y = xy[1];

  xy[0] = (-world_x*t[5] + world_y*t[2] + t[0]*t[5] - t[2]*t[3])/(t[2]*t[4] - t[1]*t[5]);
  xy[1] = (world_x*t[4] - world_y*t[1] + t[0]*t[4] + t[1]*t[3])/(t[2]*t[4] - t[1]*t[5]);
}
