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
#include <float.h>
#include <math.h>
#include "ansi.h"
#include "gdal.h"
#include "cpl_conv.h"
#include "ogr_srs_api.h"
#include "proj_api.h"


GDALDatasetH dataset = NULL;
GDALRasterBandH band = NULL;
const char * latlng = "+proj=longlat +datum=WGS84 +no_defs ";
const char * filename = "/tmp/LC08_L1TP_139045_20170304_20170316_01_T1_B1.TIF";
/* const char * filename = "/vsicurl/https://landsat-pds.s3.amazonaws.com/c1/L8/139/045/LC08_L1TP_139045_20170304_20170316_01_T1/LC08_L1TP_139045_20170304_20170316_01_T1_B1.TIF"; */
projPJ src, dst = NULL;
double t[6];
uint32_t width, height;

#define STRING_BUFFER_SIZE (1<<10)
#define TILE_SIZE (1<<8)
#define BUFFER_SIZE (ceil(sqrt(2)*TILE_SIZE))


void world_to_image(double * xy)
{
  double world_x = xy[0], world_y = xy[1];

  xy[0] = (-world_x*t[5] + world_y*t[2] + t[0]*t[5] - t[2]*t[3])/(t[2]*t[4] - t[1]*t[5]);
  xy[1] = (world_x*t[4] - world_y*t[1] + t[0]*t[4] + t[1]*t[3])/(t[2]*t[4] - t[1]*t[5]);
}

void init()
{
  GDALAllRegister();
}

void load()
{
  char * wkt = NULL, * dstProj4 = NULL;
  OGRSpatialReferenceH srs = NULL;

  /* Dataset */
  dataset = GDALOpen(filename, GA_ReadOnly);
  if(dataset == NULL) {
    fprintf(stderr, ANSI_COLOR_RED "Exiting" ANSI_COLOR_RESET);
    exit(-1);
  }
  band = GDALGetRasterBand(dataset, 1);
  width = GDALGetRasterXSize(dataset);
  height = GDALGetRasterYSize(dataset);

  /* SRS */
  srs = OSRNewSpatialReference(NULL);
  wkt = calloc(STRING_BUFFER_SIZE, sizeof(char));
  strncpy(wkt, GDALGetProjectionRef(dataset), STRING_BUFFER_SIZE);
  fprintf(stderr, ANSI_COLOR_GREEN "WKT: " ANSI_COLOR_CYAN "%s\n" ANSI_COLOR_RESET, wkt);
  OSRImportFromWkt(srs, &wkt);
  OSRExportToProj4(srs, &dstProj4);
  fprintf(stderr, ANSI_COLOR_GREEN "Proj4: " ANSI_COLOR_CYAN "%s\n" ANSI_COLOR_RESET, dstProj4);

  src = pj_init_plus(latlng);
  dst = pj_init_plus(dstProj4);

  /* Transform */
  GDALGetGeoTransform(dataset, t);
  for (int i = 0; i < 6; ++i) {
    fprintf(stderr, "%lf ", t[i]);
  }
  fprintf(stderr, "\n");

  CPLFree(dstProj4);
  OSRRelease(srs);
}

void zxy(int z, int _x, int _y)
{
  double * top;
  double * bot;
  double * left;
  double * right;
  double n = pow(2.0, z);
  double xmin = DBL_MAX, ymin = DBL_MAX;
  double xmax = DBL_MIN, ymax = DBL_MIN;
  uint16_t * texture = NULL;

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
  ymax = floor(ymax);

  texture = calloc(BUFFER_SIZE * BUFFER_SIZE, sizeof(uint16_t));
  GDALRasterIO(band, GF_Read,
               (int)xmin, (int)ymin, (int)(xmax-xmin), (int)(ymax-ymin),
               texture, BUFFER_SIZE, BUFFER_SIZE,
               GDT_UInt16, 0, 0);

  fprintf(stderr, "%lf %lf\n", xmin, xmax);
  fprintf(stderr, "%lf %lf\n", ymin, ymax);
}
