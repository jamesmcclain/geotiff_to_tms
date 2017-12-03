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

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <algorithm>
#include <limits>
#include <memory>

#include <arpa/inet.h>

#include <boost/geometry/algorithms/intersection.hpp>
#include <boost/geometry/algorithms/intersects.hpp>
#include <boost/math/constants/constants.hpp>

#include "ansi.h"
#include "greater_landsat_scene.h"
#include "lesser_landsat_scene.h"
#include "load.h"
#include "pngwrite.h"
#include "projection.h"

#include "rtree.hpp"
#include "textures.hpp"

const char * indexfile = nullptr;
const char * prefix = nullptr;
const double pi = boost::math::constants::pi<double>();
const rtree_t * rtree_ptr = nullptr;
projPJ webmercator = nullptr;
bi::managed_mapped_file * file = nullptr;

#define FUDGE (0.75)
#define ENV_INDEX "TMS_INDEX"
#define ENV_PREFIX "TMS_PREFIX"

uint8_t tile[TILE_SIZE2*4]; // RGBA ergo 4

void zxy_read(int z, int x, int y, const value_t & pair, texture_data & data);
void fetch(const value_t & pair, const box_t & tile_bb, texture_data & data);
void zxy_commit(const std::vector<texture_data> & data);
uint8_t sigmoidal(uint16_t _u);

// http://localhost:8001/{z}/{x}/{y}.png

// Global
void preload(int verbose, void * extra)
{
  const char * variable = nullptr;

  GDALAllRegister();

  // index file
  variable = std::getenv(ENV_INDEX);
  if (variable)
    indexfile = variable;
  else
    indexfile = DEFAULT_INDEXFILE;

  // prefix
  variable = std::getenv(ENV_PREFIX);
  if (variable)
    prefix = variable;
  else
    prefix = DEFAULT_READ_PREFIX;

  fprintf(stderr, ANSI_COLOR_BLUE "index file\t = " ANSI_COLOR_GREEN "%s" ANSI_COLOR_RESET "\n", indexfile);
  fprintf(stderr, ANSI_COLOR_BLUE "prefix\t\t = " ANSI_COLOR_GREEN "%s" ANSI_COLOR_RESET "\n", prefix);

  webmercator = pj_init_plus(WEBMERCATOR);

  file = new bi::managed_mapped_file(bi::open_only, indexfile);
  rtree_ptr = file->find_or_construct<rtree_t>("rtree")(params_t(), indexable_t(), equal_to_t(), allocator_t(file->get_segment_manager()));
}

// Local
void load(int verbose, void * extra)
{
}

void zxy(int fd, int z, int x, int y, int verbose, void * extra)
{
  std::vector<value_t> scene_list;
  std::vector<texture_data> texture_list;
  double z2 = pow(2.0, z);
  double xmin = (2*((x+0) / z2) - 1) * pi * RADIUS;
  double xmax = (2*((x+1) / z2) - 1) * pi * RADIUS;
  double ymin = (1 - 2*((y+0) / z2)) * pi * RADIUS;
  double ymax = (1 - 2*((y+1) / z2)) * pi * RADIUS;
  point_t smaller = point_t(std::min(xmin, xmax), std::min(ymin,ymax));
  point_t larger = point_t(std::max(xmin, xmax), std::max(ymin,ymax));
  box_t box = box_t(smaller, larger);

  rtree_ptr->query(bgi::intersects(box), std::back_inserter(scene_list));
  texture_list.resize(scene_list.size());

  for (int i = -1; i < (int)scene_list.size(); ++i) {
    if (i == -1)
      memset(tile, 0, sizeof(tile));
    else
      zxy_read(z, x, y, scene_list[i], texture_list[i]);
  }

  zxy_commit(texture_list);
  png_write(fd, tile, TILE_SIZE, TILE_SIZE, 0);
}

void fetch(const value_t & scene, const box_t & tile_bounding_box, texture_data & data)
{
  GDALDatasetH handles[3];
  GDALRasterBandH bands[3];
  char pattern[STRING_LEN];
  char filename[STRING_LEN];

  box_t image_bounding_box = box_t(point_t(0, 0), point_t(scene.second.width-1, scene.second.height-1));

  // If no intersection, short circuit
  if (!bg::intersects(tile_bounding_box, image_bounding_box))
    return;

  // Compute the texutre bounding box by intersecting the tile
  // bounding box with the image (scene) bounding box.
  bg::intersection(tile_bounding_box, image_bounding_box, data.bounding_box);

  // Compute dimensions and scales
  double texture_bounding_box_width = XMAX(data.bounding_box) - XMIN(data.bounding_box);
  double texture_bounding_box_height = YMAX(data.bounding_box) - YMIN(data.bounding_box);
  double tile_bounding_box_width = XMAX(tile_bounding_box) - XMIN(tile_bounding_box);
  double tile_bounding_box_height = YMAX(tile_bounding_box) - YMIN(tile_bounding_box);
  double w = TILE_SIZE * (texture_bounding_box_width  / tile_bounding_box_width);
  double h = TILE_SIZE * (texture_bounding_box_height / tile_bounding_box_height);

  data.texture_width  = std::max(SMALL_TILE_SIZE, static_cast<int>(std::round(w)));
  data.texture_height = std::max(SMALL_TILE_SIZE, static_cast<int>(std::round(h)));

  if (((XMIN(data.bounding_box) == XMIN(image_bounding_box) &&
        XMAX(data.bounding_box) == XMAX(image_bounding_box)) ||
       (YMIN(data.bounding_box) == YMIN(image_bounding_box) &&
        YMAX(data.bounding_box) == YMAX(image_bounding_box))) &&
      data.texture_width == data.texture_height &&
      data.texture_width == SMALL_TILE_SIZE) { // If previews are usable ...
    data.xscale = (double)data.texture_width  / scene.second.width;
    data.yscale = (double)data.texture_height / scene.second.height;
    data.bounding_box = image_bounding_box; // adjust the texture bounding box

    for (int i = 0; i < 3; ++i)
      data.textures[i] = static_cast<const uint16_t *>(scene.second.rgb[i]);
  }
  else { // If previews are not usable ...
    data.xscale = (double)(data.texture_width - FUDGE)/texture_bounding_box_width;
    data.yscale = (double)(data.texture_height - FUDGE)/texture_bounding_box_height;

    // Open handles and bands
    sprintf(pattern, "%s%s", prefix, scene.second.filename);
    for (int i = 4; i > 1; --i) {
      sprintf(filename, pattern, i);
      if ((handles[4-i] = GDALOpen(filename, GA_ReadOnly)) == NULL) {
        fprintf(stderr, ANSI_COLOR_RED "Could not open %s" ANSI_COLOR_RESET "\n", filename);
        exit(-1);
      }
      bands[4-i] = GDALGetRasterBand(handles[4-i], 1);
    }

    // Fetch textures
    // Reference: http://www.gdal.org/classGDALRasterBand.html#a30786c81246455321e96d73047b8edf1
    #pragma omp parallel for schedule(static) num_threads(3)
    for (int i = 0; i < 3; ++i) {
      data.textures[i] = std::shared_ptr<uint16_t>(new uint16_t[data.texture_width * data.texture_height],
                                                   [](uint16_t * p){ delete[] p;});
      if (GDALRasterIO(bands[i],
                       GF_Read,
                       static_cast<int>(std::floor(XMIN(data.bounding_box))),
                       static_cast<int>(std::floor(YMIN(data.bounding_box))),
                       static_cast<int>(std::ceil(XMAX(data.bounding_box)-XMIN(data.bounding_box))),
                       static_cast<int>(std::ceil(YMAX(data.bounding_box)-YMIN(data.bounding_box))),
                       std::get<std::shared_ptr<uint16_t>>(data.textures[i]).get(),
                       data.texture_width, data.texture_height,
                       GDT_UInt16, 0, 0)) {
        fprintf(stderr, ANSI_COLOR_RED "Could not fetch texture" ANSI_COLOR_RESET "\n");
        exit(-1);
      }
    }

    // Close handles
    for (int i = 0; i < 3; ++i)
      GDALClose(handles[i]);
  }

  return;
}

void zxy_read(int z, int x, int y, const value_t & scene, texture_data & data)
{
  double xmin, xmax, ymin, ymax;

  data.xs.resize(TILE_SIZE * TILE_SIZE);
  data.ys.resize(TILE_SIZE * TILE_SIZE);

  /*
    TMS to Web Mercator
    Source: https://wiki.openstreetmap.org/wiki/Slippy_map_tilenames
    Source: https://software.intel.com/en-us/node/524530
  */
  double z2 = pow(2.0, z);

  #pragma omp simd collapse(2)
  for (int j = 0; j < TILE_SIZE; ++j) {
    for ( int i = 0; i < TILE_SIZE; ++i) {
      int index = (i + j*TILE_SIZE);
      double u, v;

      u = x + (i/((double)(TILE_SIZE+1))); // tile space
      u /= z2;                             // 0-1 scaled, translated Web Mercator
      u = (2*u - 1) * pi * RADIUS;         // Web Mercator
      data.xs[index] = u;

      v = y + (j/((double)(TILE_SIZE+1))); // tile space
      v /= z2;                             // 0-1 scaled, translated Web Mercator
      v = (1 - 2*v) * pi * RADIUS;         // Web Mercator
      data.ys[index] = v;
    }
  }

  /* Web Mercator to world coordinates */
  projPJ projection = pj_init_plus(scene.second.proj4);
  pj_transform(webmercator, projection, TILE_SIZE*TILE_SIZE, 1, &(data.xs[0]), &(data.ys[0]), NULL);
  pj_free(projection);

  /* World coordinates to image coordinates */
  for (int j = 0; j < TILE_SIZE; ++j) {
    for (int i = 0; i < TILE_SIZE; ++i) {
      int index = (i + j*TILE_SIZE);

      world_to_image(&data.xs[index], &data.ys[index], scene.second.transform);
    }
  }
  xmin = *std::min_element(std::begin(data.xs), std::end(data.xs));
  xmax = *std::max_element(std::begin(data.xs), std::end(data.xs));
  ymin = *std::min_element(std::begin(data.ys), std::end(data.ys));
  ymax = *std::max_element(std::begin(data.ys), std::end(data.ys));

  /* Bounding box of the tile in image coordinates */
  box_t tile_bounding_box = box_t(point_t(std::floor(xmin), std::floor(ymin)),
                                  point_t(std::ceil(xmax), std::ceil(ymax)));

  fetch(scene, tile_bounding_box, data);
}

void zxy_commit(const std::vector<texture_data> & texture_list)
{
  for (unsigned int k = 0; k < texture_list.size(); ++k) { // For each scene
    const auto texture = texture_list[k];
    const uint16_t * rgb[3] = {nullptr, nullptr, nullptr};

    if (std::holds_alternative<const uint16_t *>(texture.textures[0])) {
      for (int i = 0; i < 3; ++i)
        rgb[i] = std::get<const uint16_t *>(texture.textures[i]);
    }
    else if (std::holds_alternative<std::shared_ptr<uint16_t>>(texture.textures[0])) {
      for (int i = 0; i < 3; ++i)
        rgb[i] = std::get<std::shared_ptr<uint16_t>>(texture.textures[i]).get();
    }

    #pragma omp simd collapse(2)
    for (int j = 0; j < TILE_SIZE; ++j) { // tile coordinate
      for (int i = 0; i < TILE_SIZE; ++i) { // tile coordinate
        int tile_index = (i + j*TILE_SIZE)*4;
        double x = texture.xs[tile_index/4], y = texture.ys[tile_index/4]; // scene image coordinates
        int u = static_cast<int>(std::round(texture.xscale*(x-XMIN(texture.bounding_box)))); // texture coordinate
        int v = static_cast<int>(std::round(texture.yscale*(y-YMIN(texture.bounding_box)))); // texture coordinate

        if (0 <= u && u < static_cast<int>(texture.texture_width) &&
            0 <= v && v < static_cast<int>(texture.texture_height)) {
          int texture_index = u + v*(texture.texture_width);
          uint8_t red, byte = 0;

          byte |= red = sigmoidal(rgb[0][texture_index]);
          if (tile[tile_index + 3] == 0 /* alpha channel */ ||
              tile[tile_index + 0] < red /* red channel */) { // write into empty pixels
            tile[tile_index + 0] = red;
            byte |= tile[tile_index + 1] = sigmoidal(rgb[1][texture_index]);
            byte |= tile[tile_index + 2] = sigmoidal(rgb[2][texture_index]);
            tile[tile_index + 3] = (byte ? -1 : 0);
          }
        }
      }
    }
  }
}

uint8_t sigmoidal(uint16_t _u)
{
  if (!_u) return 0;

  double u = ((double)_u) / 23130.235294118;
  double beta = 10, alpha = 0.50;
  double numer = 1/(1+exp(beta*(alpha-u))) - 1/(1+exp(beta));
  double denom = 1/(1+exp(beta*(alpha-1))) - 1/(1+exp(beta*alpha));
  double gu = std::max(0.0, std::min(1.0, numer / denom));
  return ((1<<8)-1)*gu;
}
