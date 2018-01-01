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

#ifndef __LANDSAT_H__
#define __LANDSAT_H__

#include "gdal.h"
#include "cpl_conv.h"
#include "ogr_srs_api.h"
#include "proj_api.h"

#include "constants.h"
#include "lesser_landsat_scene.h"


struct periphery_struct {
  double   top[TILE_SIZE<<1]; // two coordinates ergo shift
  double   bot[TILE_SIZE<<1];
  double  left[TILE_SIZE<<1];
  double right[TILE_SIZE<<1];
};

union coordinates_struct {
  struct periphery_struct periphery;
  double patch [TILE_SIZE2<<1]; // two coordinates ergo shift
};

typedef struct greater_landsat_scene_struct {
  struct lesser_landsat_scene_struct lesser;

  // Datasets
  GDALDatasetH r_dataset;
  GDALDatasetH g_dataset;
  GDALDatasetH b_dataset;

  // Bands
  GDALRasterBandH r_band;
  GDALRasterBandH g_band;
  GDALRasterBandH b_band;

  // Textures
  uint16_t r_texture[TILE_SIZE2];
  uint16_t g_texture[TILE_SIZE2];
  uint16_t b_texture[TILE_SIZE2];

  // Coordinates
  union coordinates_struct coordinates;

  // bounding box for the current tile in source image coordinates
  double xmin, xmax, ymin, ymax;

  // Information about the source image.  The bounding box given on
  // the second and third lines is the intersection of the source
  // image with the tile.
  uint32_t src_width, src_height;
  uint32_t src_window_xmin, src_window_ymin;
  uint32_t src_window_width, src_window_height;

  // The starting point of the texture within the tile, as well as its
  // height and width.  Here, the red, green, and blue textures are
  // referred to in the singular.
  uint32_t tile_window_xmin, tile_window_ymin;
  uint32_t tile_window_width, tile_window_height;

  int dirty;

} landsat_scene;

#endif
