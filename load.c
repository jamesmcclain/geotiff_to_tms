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
#include <arpa/inet.h>
#include <png.h>
#include "ansi.h"
#include "gdal.h"
#include "cpl_conv.h"
#include "ogr_srs_api.h"
#include "proj_api.h"


GDALDatasetH r_dataset = NULL;
GDALDatasetH g_dataset = NULL;
GDALDatasetH b_dataset = NULL;
GDALRasterBandH r_band = NULL;
GDALRasterBandH g_band = NULL;
GDALRasterBandH b_band = NULL;
const char * latlng = "+proj=longlat +datum=WGS84 +no_defs ";
const char * filename = "/tmp/LC08_L1TP_139045_20170304_20170316_01_T1_B%d.TIF";
/* const char * filename = "/vsicurl/https://landsat-pds.s3.amazonaws.com/c1/L8/139/045/LC08_L1TP_139045_20170304_20170316_01_T1/LC08_L1TP_139045_20170304_20170316_01_T1_B%d.TIF"; */
projPJ src, dst = NULL;
double t[6];
uint32_t width, height;

#define STRING_BUFFER_SIZE (1<<10)
#define TILE_SIZE (1<<8)
#define BUFFER_SIZE (ceil(sqrt(2)*TILE_SIZE))


void write_png(char *file_name,
               const uint16_t * r_texture,
               const uint16_t * g_texture,
               const uint16_t * b_texture,
               int width, int height)
{
  FILE *fp;
  png_structp png_ptr;
  png_infop info_ptr;
  png_color_8 sig_bit;
  uint16_t * tile = NULL;

  tile = calloc(width*height*3, sizeof(*tile));
  for (int j = 0; j < height; ++j) {
    for (int i = 0; i < width; ++i) {
      tile[3*i + 3*j*width + 0] = htons(r_texture[i + j*width]);
      tile[3*i + 3*j*width + 1] = htons(g_texture[i + j*width]);
      tile[3*i + 3*j*width + 2] = htons(b_texture[i + j*width]);
    }
  }

  fp = fopen(file_name, "wb");

  png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  info_ptr = png_create_info_struct(png_ptr);
  if (setjmp(png_jmpbuf(png_ptr))) {
    fprintf(stderr, ANSI_COLOR_RED "libpng issue\n" ANSI_COLOR_RESET);
    exit(-1);
  }

  png_init_io(png_ptr, fp);
  png_set_IHDR(png_ptr, info_ptr,
               width, height, 8*sizeof(*tile),
               PNG_COLOR_TYPE_RGB,
               PNG_INTERLACE_NONE,
               PNG_COMPRESSION_TYPE_BASE,
               PNG_FILTER_TYPE_BASE);

  sig_bit.red   = 8*sizeof(*r_texture);
  sig_bit.green = 8*sizeof(*g_texture);
  sig_bit.blue  = 8*sizeof(*b_texture);
  /* sig_bit.gray = 8*sizeof(*texture); */
  /* sig_bit.alpha = 8; */
  png_set_sBIT(png_ptr, info_ptr, &sig_bit);
  png_write_info(png_ptr, info_ptr);
  /* png_set_invert_mono(png_ptr); */
  /* png_set_shift(png_ptr, &sig_bit); */
  /* png_set_packing(png_ptr); */
  /* png_set_swap_alpha(png_ptr); */
  /* png_set_filler(png_ptr, 0, PNG_FILLER_BEFORE); */
  /* png_set_bgr(png_ptr); */
  /* png_set_swap(png_ptr); */
  /* png_set_packswap(png_ptr); */

  png_bytep row_pointers[height];
  for (png_uint_32 i = 0; i < height; i++)
    row_pointers[i] = (png_bytep)(tile + i*width*3);
  png_write_image(png_ptr, row_pointers);
  png_write_end(png_ptr, info_ptr);
  png_destroy_write_struct(&png_ptr, &info_ptr);
  fclose(fp);
  free(tile);
}

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
  char r_filename[(1<<8)];
  char g_filename[(1<<8)];
  char b_filename[(1<<8)];

  sprintf(r_filename, filename, 1);
  sprintf(g_filename, filename, 2);
  sprintf(b_filename, filename, 3);

  /* Dataset */
  r_dataset = GDALOpen(r_filename, GA_ReadOnly);
  g_dataset = GDALOpen(g_filename, GA_ReadOnly);
  b_dataset = GDALOpen(b_filename, GA_ReadOnly);
  if(r_dataset == NULL) {
    fprintf(stderr, ANSI_COLOR_RED "GDALOpen issue\n" ANSI_COLOR_RESET);
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
  fprintf(stderr, ANSI_COLOR_GREEN "WKT: " ANSI_COLOR_CYAN "%s\n" ANSI_COLOR_RESET, wkt);
  OSRImportFromWkt(srs, &wkt);
  OSRExportToProj4(srs, &dstProj4);
  fprintf(stderr, ANSI_COLOR_GREEN "Proj4: " ANSI_COLOR_CYAN "%s\n" ANSI_COLOR_RESET, dstProj4);

  src = pj_init_plus(latlng);
  dst = pj_init_plus(dstProj4);

  /* Transform */
  GDALGetGeoTransform(r_dataset, t);
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
  uint16_t * r_texture = NULL;
  uint16_t * g_texture = NULL;
  uint16_t * b_texture = NULL;

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

  r_texture = calloc(BUFFER_SIZE * BUFFER_SIZE, sizeof(*r_texture));
  g_texture = calloc(BUFFER_SIZE * BUFFER_SIZE, sizeof(*g_texture));
  b_texture = calloc(BUFFER_SIZE * BUFFER_SIZE, sizeof(*b_texture));
  if (GDALRasterIO(r_band, GF_Read,
                   (int)xmin, (int)ymin, (int)(xmax-xmin), (int)(ymax-ymin),
                   r_texture, BUFFER_SIZE, BUFFER_SIZE,
                   GDT_UInt16, 0, 0)) {
    fprintf(stderr, ANSI_COLOR_RED "%lf %lf %lf %lf" ANSI_COLOR_RESET "\n", xmin, ymin, xmax, ymax);
    fprintf(stderr, ANSI_COLOR_RED "GDALRasterIO issue (red band)\n" ANSI_COLOR_RESET);
    exit(-1);
  }
  if (GDALRasterIO(g_band, GF_Read,
                   (int)xmin, (int)ymin, (int)(xmax-xmin), (int)(ymax-ymin),
                   g_texture, BUFFER_SIZE, BUFFER_SIZE,
                   GDT_UInt16, 0, 0)) {
    fprintf(stderr, ANSI_COLOR_RED "GDALRasterIO issue (green band)\n" ANSI_COLOR_RESET);
    exit(-1);
  }
  if (GDALRasterIO(b_band, GF_Read,
                   (int)xmin, (int)ymin, (int)(xmax-xmin), (int)(ymax-ymin),
                   b_texture, BUFFER_SIZE, BUFFER_SIZE,
                   GDT_UInt16, 0, 0)) {
    fprintf(stderr, ANSI_COLOR_RED "GDALRasterIO issue (blue band)\n" ANSI_COLOR_RESET);
    exit(-1);
  }
  write_png("/tmp/tile.png", r_texture, g_texture, b_texture, BUFFER_SIZE, BUFFER_SIZE);

  free(b_texture);
  free(g_texture);
  free(r_texture);
  free(right);
  free(left);
  free(bot);
  free(top);
}

/* void whole() */
/* { */
/*   uint16_t * texture = NULL; */

/*   texture = calloc(width * height, sizeof(*texture)); */
/*   if (GDALRasterIO(r_band, GF_Read, */
/*                    0, 0, width, height, */
/*                    texture, width, height, */
/*                    GDT_UInt16, 0, 0)) { */
/*     fprintf(stderr, ANSI_COLOR_RED "GDALRasterIO issue" ANSI_COLOR_RESET); */
/*     exit(-1); */
/*   } */
/*   write_png("/tmp/whole.png", texture, width, height); */
/*   free(texture); */
/* } */
