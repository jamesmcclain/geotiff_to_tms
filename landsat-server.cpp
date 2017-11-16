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
#include <cmath>
#include <limits>

#include <arpa/inet.h>

#include <boost/compute/detail/lru_cache.hpp>
#include <boost/geometry/algorithms/intersection.hpp>
#include <boost/geometry/algorithms/intersects.hpp>

#include "ansi.h"
#include "greater_landsat_scene.h"
#include "lesser_landsat_scene.h"
#include "load.h"
#include "pngwrite.h"
#include "projection.h"

#include "landsat_scene_handles.hpp"
#include "rtree.hpp"

namespace bcd = boost::compute::detail;

typedef bcd::lru_cache<const char *, landsat_scene_handles> cache;

const char * indexfile = nullptr;
const char * prefix = nullptr;
rtree_t * rtree_ptr = nullptr;
projPJ webmercator = nullptr;
bi::managed_mapped_file * file = nullptr;
cache * lru = nullptr;

#define MAX_LEN (1<<10)
#define XMIN(b) ((b).min_corner().get<0>())
#define YMIN(b) ((b).min_corner().get<1>())
#define XMAX(b) ((b).max_corner().get<0>())
#define YMAX(b) ((b).max_corner().get<1>())

// #define OVERZOOM_FACTOR (2)
// #define OVERZOOM (overzoom ? OVERZOOM_FACTOR : 1)

uint8_t tile[TILE_SIZE*4]; // RGBA ergo 4

int fetch(const value_t & pair, const box_t & tile_bb, uint16_t * textures[3], ibox_t & texture_box);
uint8_t sigmoidal(uint16_t _u);
// void load_scene(landsat_scene * s, int verbose);
// void zxy_approx_commit();
// void zxy_approx_read(int z, int _x, int _y, landsat_scene * s, int verbose);
// void zxy_exact_commit(int overzoom);
void zxy_exact_read(int z, int x, int y, const value_t & pair);


void preload(int verbose, void * extra)
{
  GDALAllRegister();

  indexfile = DEFAULT_INDEXFILE;
  prefix = DEFAULT_PREFIX;
  webmercator = pj_init_plus(WEBMERCATOR);

  file = new bi::managed_mapped_file(bi::open_only, indexfile);
  rtree_ptr = file->find_or_construct<rtree_t>("rtree")(params_t(), indexable_t(), equal_to_t(), allocator_t(file->get_segment_manager()));
  lru = new cache(1<<16);
}

void load(int verbose, void * extra)
{
}

void zxy(int fd, int z, int x, int y, int verbose, void * extra)
{
  std::vector<value_t> results;

  {
    double z2 = pow(2.0, z);
    double xmin = (2*((x+0) / z2) - 1) * M_PI * RADIUS;
    double xmax = (2*((x+1) / z2) - 1) * M_PI * RADIUS;
    double ymin = (1 - 2*((y+0) / z2)) * M_PI * RADIUS;
    double ymax = (1 - 2*((y+1) / z2)) * M_PI * RADIUS;
    point_t smaller = point_t(std::min(xmin, xmax), std::min(ymin,ymax));
    point_t larger = point_t(std::max(xmin, xmax), std::max(ymin,ymax));
    box_t box = box_t(smaller, larger);

    rtree_ptr->query(bgi::intersects(box), std::back_inserter(results));
  }

  int exact = (z < 7) ? 1 : 0;
  // int overzoom = (z < 3) ? 1 : 0;

  // #pragma omp parallel for schedule(dynamic, 1)
  for (int i = -1; i < (int)results.size(); ++i) {
    if (i == -1)
      memset(tile, 0, sizeof(tile));
    else /*if (exact)*/
      zxy_exact_read(z, x, y, results[i]);
    // else if (!exact)
    //   zxy_approx_read(z, x, y, results[i], verbose);
  }

  // // Sample from textures to produce tile
  // if (exact)
  //   zxy_exact_commit(overzoom);
  // else
  //   zxy_approx_commit();

  // write_png(fd, tile, TILE_SIZE, TILE_SIZE, 0);
}

int fetch(const value_t & pair,
          const box_t & tile_bb,
          uint16_t * textures[3],
          ibox_t & texture_box)
{
  GDALDatasetH handles[3];
  GDALRasterBandH bands[3];
  char pattern[MAX_LEN];
  char filename[MAX_LEN];

  box_t image_bb = box_t(point_t(0, 0), point_t(pair.second.width-1, pair.second.height-1));

  // If no intersection, short circuit
  if (!bg::intersects(tile_bb, image_bb))
    return 0;

  // Compute image box
  box_t image_box;
  bg::intersection(tile_bb, image_bb, image_box);

  // Compute texture box
  double xscale = TILE_SIZE / (XMAX(tile_bb) - XMIN(tile_bb));
  double yscale = TILE_SIZE / (YMAX(tile_bb) - YMIN(tile_bb));
  texture_box.min_corner().set<0>(static_cast<int>(floor(xscale*(XMIN(image_box)-XMIN(tile_bb)))));
  texture_box.min_corner().set<1>(static_cast<int>(floor(yscale*(YMIN(image_box)-YMIN(tile_bb)))));
  texture_box.max_corner().set<0>(static_cast<int>(ceil(xscale*(XMAX(image_box)-XMIN(tile_bb)))));
  texture_box.max_corner().set<1>(static_cast<int>(ceil(yscale*(YMAX(image_box)-YMIN(tile_bb)))));

  // Open reg, green, blue datasets
  sprintf(pattern, "%s%s", DEFAULT_PREFIX, pair.second.filename);
  for (int i = 4; i > 1; --i) {
    sprintf(filename, pattern, i);
    if ((handles[4-i] = GDALOpen(filename, GA_ReadOnly)) == NULL) exit(-1); // XXX
    bands[4-i] = GDALGetRasterBand(handles[4-i], 1);
  }

  // Fetch textures
  // Reference: http://www.gdal.org/classGDALRasterBand.html#a30786c81246455321e96d73047b8edf1
  for (int i = 0; i < 3; ++i) {
    if (GDALRasterIO(bands[i],
                     GF_Read,
                     static_cast<int>(floor(XMIN(image_box))),
                     static_cast<int>(floor(YMIN(image_box))),
                     static_cast<int>(ceil(XMAX(image_box)-XMIN(image_box))),
                     static_cast<int>(ceil(YMAX(image_box)-YMIN(image_box))),
                     textures[i],
                     XMAX(texture_box)-XMIN(texture_box), YMAX(texture_box)-YMIN(texture_box),
                     GDT_UInt16, 0, 0)) exit(-1);
  }

  return 1;
}

void zxy_exact_read(int z, int x, int y, const value_t & pair)
{
  double xmin = std::numeric_limits<double>::max();
  double ymin = std::numeric_limits<double>::max();
  double xmax = std::numeric_limits<double>::min();
  double ymax = std::numeric_limits<double>::min();
  auto xs = std::vector<double>(TILE_SIZE * TILE_SIZE);
  auto ys = std::vector<double>(TILE_SIZE * TILE_SIZE);

  /*
    TMS to Pseudo Web Mercator
    Source: https://wiki.openstreetmap.org/wiki/Slippy_map_tilenames
    Source: https://software.intel.com/en-us/node/524530
  */
  {
    double z2 = pow(2.0, z);

    // #pragma omp simd collapse(2)
    for (int j = 0; j < TILE_SIZE; ++j) {
      for ( int i = 0; i < TILE_SIZE; ++i) {
        int index = (i + j*TILE_SIZE);
        double u, v;

        u = x + (i/((double)(TILE_SIZE+1))); // tile space
        u /= z2;                             // 0-1 scaled, translated Web Mercator
        u = (2*u - 1) * M_PI * RADIUS;       // Web Mercator
        xs[index] = u;

        v = y + (j/((double)(TILE_SIZE+1))); // tile space
        v /= z2;                             // 0-1 scaled, translated Web Mercator
        v = (1 - 2*v) * M_PI * RADIUS;       // Web Mercator
        ys[index] = v;
      }
    }
  }

  /* Web Mercator to world coordinates */
  projPJ projection = pj_init_plus(pair.second.proj4);
  pj_transform(webmercator, projection, TILE_SIZE*TILE_SIZE, 1, &xs[0], &ys[0], NULL);
  pj_free(projection);

  /* World coordinates to image coordinates */
  for (int j = 0; j < TILE_SIZE; ++j) {
    for (int i = 0; i < TILE_SIZE; ++i) {
      int index = (i + j*TILE_SIZE);
      double uv[2] = {xs[index], ys[index]};

      world_to_image(uv, pair.second.transform);
      xmin = fmin(uv[0], xmin);
      xmax = fmax(uv[0], xmax);
      ymin = fmin(uv[1], ymin);
      ymax = fmax(uv[1], ymax);
    }
  }

  /* Bounding box of the tile in image coordinates */
  box_t bounding_box = box_t(point_t(round(xmin), round(ymin)),
                             point_t(round(xmax), round(ymax)));

  /* Textures and texture bounding box (the latter in tile coordinates) */
  ibox_t texture_box;
  uint16_t * textures[3];
  for (int i = 0; i < 3; ++i)
    textures[i] = static_cast<uint16_t *>(calloc(TILE_SIZE * TILE_SIZE, sizeof(uint16_t)));

  fetch(pair, bounding_box, textures, texture_box);
  fprintf(stderr, "XXX\n");
  // s->dirty = fetch(s, overzoom, verbose);
}

// void zxy_approx_commit()
// {
//   for (int i = 0; i < scene_count; ++i) {

//     landsat_scene * s = scene + i;

//     if (!s->dirty) continue;

//     for (unsigned int j = 0; (j < s->tile_window_height); ++j) {
//       for (unsigned int i = 0; (i < s->tile_window_width); ++i) {
//         uint8_t red, byte = 0;
//         int tile_index = ((i + s->tile_window_xmin) + (j + s->tile_window_ymin)*TILE_SIZE)<<2;
//         int texture_index = (i + j*s->tile_window_width);

//         byte |= red = sigmoidal(s->r_texture[texture_index]);
//         if (tile[tile_index + 3] == 0 || tile[tile_index + 0] < red) { // write into empty pixels
//           tile[tile_index + 0] = red;
//           byte |= tile[tile_index + 1] = sigmoidal(s->g_texture[texture_index]);
//           byte |= tile[tile_index + 2] = sigmoidal(s->b_texture[texture_index]);
//           tile[tile_index + 3] = (byte ? -1 : 0);
//         }
//       }
//     }
//   }
// }

// void zxy_approx_read(int z, int _x, int _y, landsat_scene * s, int verbose)
// {
//   double * top   = s->coordinates.periphery.top;
//   double * bot   = s->coordinates.periphery.bot;
//   double * left  = s->coordinates.periphery.left;
//   double * right = s->coordinates.periphery.right;
//   double xmin = DBL_MAX, ymin = DBL_MAX;
//   double xmax = DBL_MIN, ymax = DBL_MIN;

//   if (verbose) {
//     fprintf(stderr,
//             ANSI_COLOR_YELLOW "z=%d x=%d y=%d pid=%d" ANSI_COLOR_RESET "\n",
//             z, _x, _y, getpid());
//   }

//   /*
//     TMS to Pseudo Web Mercator
//     Source: https://wiki.openstreetmap.org/wiki/Slippy_map_tilenames
//   */
//   for (int i = 0; i < (TILE_SIZE<<1); i+=2) {
//     // top, bottom x-values
//     top[i+0] = _x + (i/(TILE_SIZE*2.0));     // tile space
//     top[i+0] /= pow(2.0, z);                 // 0-1 scaled, translated Web Mercator
//     top[i+0] = (2*top[i+0] - 1) * M_PI;      // Web Mercator in radians
//     top[i+0] = bot[i+0] = top[i+0] * RADIUS; // Web Mercator in radians*radius

//     // top, bottom y-values
//     top[i+1] = (1 - 2*((_y+0) / pow(2.0, z))) * M_PI * RADIUS;
//     bot[i+1] = (1 - 2*((_y+1) / pow(2.0, z))) * M_PI * RADIUS;

//     // left, right x-values
//     left[i+0]  = (2*((_x+0) / pow(2.0, z)) - 1) * M_PI * RADIUS;
//     right[i+0] = (2*((_x+1) / pow(2.0, z)) - 1) * M_PI * RADIUS;

//     // left, right y-values
//     left[i+1] = _y + (i/(TILE_SIZE*2.0));
//     left[i+1] /= pow(2.0, z);
//     left[i+1] = right[i+1] = (1 - 2*left[i+1]) * M_PI * RADIUS;
//   }

//   /* Web Mercator to world coordinates */
//   pj_transform(webmercator_pj, s->lesser.projection,
//                TILE_SIZE<<2, 2, top, top+1, NULL); //top, bot, left, right ergo shift

//   /* World coordinates to image coordinates */
//   for (int i = 0; i < (TILE_SIZE<<1); i+=2) {
//     world_to_image(top+i, s->lesser.transform);
//     world_to_image(bot+i, s->lesser.transform);
//     world_to_image(left+i, s->lesser.transform);
//     world_to_image(right+i, s->lesser.transform);
//     xmin = fmin(right[i], fmin(left[i], fmin(bot[i], fmin(top[i], xmin))));
//     xmax = fmax(right[i], fmax(left[i], fmax(bot[i], fmax(top[i], xmax))));
//     ymin = fmin(right[i+1], fmin(left[i+1], fmin(bot[i+1], fmin(top[i+1], ymin))));
//     ymax = fmax(right[i+1], fmax(left[i+1], fmax(bot[i+1], fmax(top[i+1], ymax))));
//   }

//   /* Bounding box of the tile */
//   s->xmin = round(xmin);
//   s->xmax = round(xmax);
//   s->ymin = round(ymin);
//   s->ymax = round(ymax);

//   s->dirty = fetch(s, 0, verbose);
// }

// void zxy_exact_commit(int overzoom)
// {
//   for (int i = 0; i < scene_count; ++i) {
//     landsat_scene * s = scene + i;

//     if (!s->dirty) continue;

//     const double * patch = s->coordinates.patch;

//     for (unsigned int j = 0; j < TILE_SIZE; ++j) {
//       for (unsigned int i = 0; i < TILE_SIZE; ++i) {
//         int patch_index = (i + j*TILE_SIZE)<<1;
//         int tile_index = patch_index<<1;
//         double u = patch[patch_index + 0];
//         double v = patch[patch_index + 1];

//         u = round((((u - s->xmin)/(s->xmax - s->xmin)) * (OVERZOOM*TILE_SIZE - 1))) - OVERZOOM*s->tile_window_xmin;
//         v = round((((v - s->ymin)/(s->ymax - s->ymin)) * (OVERZOOM*TILE_SIZE - 1))) - OVERZOOM*s->tile_window_ymin;

//         if ((0 <= u && u < OVERZOOM*s->tile_window_width) &&
//             (0 <= v && v < OVERZOOM*s->tile_window_height)) {
//           uint8_t red, byte = 0;
//           int texture_index = (u + v*OVERZOOM*s->tile_window_width);

//           byte |= red = sigmoidal(s->r_texture[texture_index]);
//           if (tile[tile_index + 3] == 0 || tile[tile_index + 0] < red) { // write into empty pixels
//             tile[tile_index + 0] = red;
//             byte |= tile[tile_index + 1] = sigmoidal(s->g_texture[texture_index]);
//             byte |= tile[tile_index + 2] = sigmoidal(s->b_texture[texture_index]);
//             tile[tile_index + 3] = (byte ? -1 : 0);
//           }
//         }
//       }
//     }
//   }
// }

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
