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
#include <cfloat>
#include <cmath>

#include <arpa/inet.h>

#include <boost/compute/detail/lru_cache.hpp>

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
projPJ webmercator_pj = nullptr;
bi::managed_mapped_file * file = nullptr;
cache * lru = nullptr;

// #define OVERZOOM_FACTOR (2)
// #define OVERZOOM (overzoom ? OVERZOOM_FACTOR : 1)

// const char * webmercator = "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext +no_defs";
// projPJ webmercator_pj = NULL;
// landsat_scene * scene = NULL;
// int scene_count = -1;
// uint8_t tile[TILE_SIZE2<<2]; // RGBA ergo 4

// int fetch(landsat_scene * s, int overzoom, int verbose);
// uint8_t sigmoidal(uint16_t _u);
// void load_scene(landsat_scene * s, int verbose);
// void zxy_approx_commit();
// void zxy_approx_read(int z, int _x, int _y, landsat_scene * s, int verbose);
// void zxy_exact_commit(int overzoom);
// void zxy_exact_read(int z, int _x, int _y, int overzoom, landsat_scene * s, int verbose);


void preload(int verbose, void * extra)
{
  GDALAllRegister();

  indexfile = DEFAULT_INDEXFILE;
  prefix = DEFAULT_PREFIX;
  webmercator_pj = pj_init_plus(WEBMERCATOR);

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

  fprintf(stderr, "%ld %ld\n", rtree_ptr->size(), results.size());

  // int exact = (z < 7) ? 1 : 0;
  // int overzoom = (z < 3) ? 1 : 0;

  // #pragma omp parallel for schedule(static, 1)
  // for (int i = -1; i < scene_count; ++i) {
  //   if (i == -1)
  //     memset(tile, 0, sizeof(tile));
  //   else if (exact)
  //     zxy_exact_read(z, _x, _y, overzoom, scene + i, verbose);
  //   else if (!exact)
  //     zxy_approx_read(z, _x, _y, scene + i, verbose);
  // }

  // // Sample from textures to produce tile
  // if (exact)
  //   zxy_exact_commit(overzoom);
  // else
  //   zxy_approx_commit();

  // write_png(fd, tile, TILE_SIZE, TILE_SIZE, 0);
}

// int fetch(landsat_scene * s, int overzoom, int verbose)
// {
//   int src_window_xmax, src_window_ymax;

//   // If the tile is disjoint from the source image, exit.
//   if ((s->xmin >= s->src_width-1) ||
//       (s->ymin >= s->src_height-1) ||
//       (s->xmax <= 0) ||
//       (s->ymax <= 0)) {
//     return 0;
//   }

//   // Clip window to source data
//   if (s->xmin < 0) s->src_window_xmin = 0; else s->src_window_xmin = s->xmin;
//   if (s->ymin < 0) s->src_window_ymin = 0; else s->src_window_ymin = s->ymin;
//   if (s->xmax > s->src_width-1) src_window_xmax = s->src_width-1; else src_window_xmax = s->xmax;
//   if (s->ymax > s->src_height-1) src_window_ymax = s->src_height-1; else src_window_ymax = s->ymax;
//   s->src_window_width = src_window_xmax - s->src_window_xmin;
//   s->src_window_height = src_window_ymax - s->src_window_ymin;
//   s->tile_window_xmin = floor(TILE_SIZE*((s->src_window_xmin - s->xmin)/(s->xmax - s->xmin)));
//   s->tile_window_ymin = floor(TILE_SIZE*((s->src_window_ymin - s->ymin)/(s->ymax - s->ymin)));
//   s->tile_window_width = ceil(TILE_SIZE*((s->src_window_width)/(s->xmax - s->xmin)));
//   s->tile_window_height = ceil(TILE_SIZE*((s->src_window_height)/(s->ymax - s->ymin)));

//   if (verbose) {
//     fprintf(stderr,
//             ANSI_COLOR_CYAN "width=%d height=%d "
//             ANSI_COLOR_MAGENTA "tx=%d ty=%d twidth=%d theight=%d "
//             ANSI_COLOR_BLUE "pid=%d"
//             ANSI_COLOR_RESET "\n",
//             s->src_window_width, s->src_window_height,
//             s->tile_window_xmin, s->tile_window_ymin,
//             s->tile_window_width, s->tile_window_height,
//             getpid());
//   }

//   // Reference: http://www.gdal.org/classGDALRasterBand.html#a30786c81246455321e96d73047b8edf1
//   if (GDALRasterIO(s->r_band, GF_Read,
//                    s->src_window_xmin, s->src_window_ymin, s->src_window_width, s->src_window_height,
//                    s->r_texture, OVERZOOM * s->tile_window_width, OVERZOOM * s->tile_window_height,
//                    GDT_UInt16, 0, 0)) {
//     fprintf(stderr, ANSI_COLOR_RED "GDALRasterIO problem (red band)" ANSI_COLOR_RESET "\n");
//     exit(-1);
//   }
//   if (GDALRasterIO(s->g_band, GF_Read,
//                    s->src_window_xmin, s->src_window_ymin, s->src_window_width, s->src_window_height,
//                    s->g_texture, OVERZOOM * s->tile_window_width, OVERZOOM * s->tile_window_height,
//                    GDT_UInt16, 0, 0)) {
//     fprintf(stderr, ANSI_COLOR_RED "GDALRasterIO problem (green band)" ANSI_COLOR_RESET "\n");
//     exit(-1);
//   }
//   if (GDALRasterIO(s->b_band, GF_Read,
//                    s->src_window_xmin, s->src_window_ymin, s->src_window_width, s->src_window_height,
//                    s->b_texture, OVERZOOM * s->tile_window_width, OVERZOOM * s->tile_window_height,
//                    GDT_UInt16, 0, 0)) {
//     fprintf(stderr, ANSI_COLOR_RED "GDALRasterIO problem (blue band)" ANSI_COLOR_RESET "\n");
//     exit(-1);
//   }

//   return 1;
// }

// void zxy_exact_read(int z, int _x, int _y, int overzoom, landsat_scene * s, int verbose)
// {
//   double * patch = s->coordinates.patch;
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
//   for (int j = 0; j < TILE_SIZE; ++j) {
//     for (int i = 0; i < TILE_SIZE; ++i) {
//       int index = (i + j*TILE_SIZE)<<1;
//       patch[index+0] = _x + (i/((double)(TILE_SIZE+1)));       // tile space
//       patch[index+0] /= pow(2.0, z);                           // 0-1 scaled, translated Web Mercator
//       patch[index+0] = (2*patch[index+0] - 1) * M_PI * RADIUS; // Web Mercator
//       patch[index+1] = _y + (j/((double)(TILE_SIZE+1)));
//       patch[index+1] /= pow(2.0, z);
//       patch[index+1] = (1 - 2*patch[index+1]) * M_PI * RADIUS;
//     }
//   }

//   /* Web Mercator to world coordinates */
//   pj_transform(webmercator_pj, s->lesser.projection,
//                TILE_SIZE*TILE_SIZE, 2,
//                patch, patch+1, NULL);

//   /* World coordinates to image coordinates */
//   for (int j = 0; j < TILE_SIZE; ++j) {
//     for (int i = 0; i < TILE_SIZE; ++i) {
//       int index = (i + j*TILE_SIZE)<<1;
//       world_to_image(patch+index, s->lesser.transform);
//       xmin = fmin(patch[index+0], xmin);
//       xmax = fmax(patch[index+0], xmax);
//       ymin = fmin(patch[index+1], ymin);
//       ymax = fmax(patch[index+1], ymax);
//     }
//   }

//   /* Bounding box of the tile */
//   s->xmin = round(xmin);
//   s->xmax = round(xmax);
//   s->ymin = round(ymin);
//   s->ymax = round(ymax);

//   s->dirty = fetch(s, overzoom, verbose);
// }

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
