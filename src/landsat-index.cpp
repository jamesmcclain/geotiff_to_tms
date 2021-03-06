/*
 * Copyright (c) 2017-2018, James McClain
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
#include <cstring>
#include <algorithm>
#include <limits>

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "gdal.h"
#include "cpl_conv.h"
#include "ogr_srs_api.h"
#include "proj_api.h"

#include "ansi.h"
#include "lesser_landsat_scene.h"
#include "projection.h"

#include "rtree.hpp"

#define INCREMENTS (1<<8)
#define POSTFIX "index.html"

const char * webmercator = WEBMERCATOR;
projPJ webmercator_pj = NULL;


/* Reference: https://stackoverflow.com/questions/4986673/c11-rvalues-and-move-semantics-confusion-return-statement */
box_t bounding_box(lesser_landsat_scene_struct scene)
{
  double xmin, xmax, ymin, ymax;
  std::vector<double> t_x = std::vector<double>(INCREMENTS);
  std::vector<double> t_y = std::vector<double>(INCREMENTS);
  std::vector<double> b_x = std::vector<double>(INCREMENTS);
  std::vector<double> b_y = std::vector<double>(INCREMENTS);
  std::vector<double> l_x = std::vector<double>(INCREMENTS);
  std::vector<double> l_y = std::vector<double>(INCREMENTS);
  std::vector<double> r_x = std::vector<double>(INCREMENTS);
  std::vector<double> r_y = std::vector<double>(INCREMENTS);
  projPJ projection_pj = NULL;

  if (scene.width == BAD || scene.height == BAD)
    return box_t(point_t(std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest()),
                 point_t(std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest()));

  // Get world coordinates around the periphery
  #pragma omp simd
  for (int i = 0; i < INCREMENTS; ++i) {
    t_x[i] = b_x[i] = 0.5 + static_cast<double>(i*scene.width)/INCREMENTS;  // x values across the top and bottom
    t_y[i] = l_x[i] = 0.5;                                                  // y values across the top and x values on the left
    b_y[i] = 0.5 + static_cast<double>(scene.height);                       // y values across the bottom
    r_x[i] = 0.5 + static_cast<double>(scene.width);                        // x values on the right
    l_y[i] = r_y[i] = 0.5 + static_cast<double>(i*scene.height)/INCREMENTS; // y values on the left and right

    image_to_world(&t_x[i], &t_y[i], scene.transform);
    image_to_world(&b_x[i], &b_y[i], scene.transform);
    image_to_world(&l_x[i], &l_y[i], scene.transform);
    image_to_world(&r_x[i], &r_y[i], scene.transform);
  }

  // Get WebMercator coordinates
  projection_pj = pj_init_plus(scene.proj4);
  pj_transform(projection_pj, webmercator_pj, INCREMENTS, 1, &t_x[0], &t_y[0], NULL);
  pj_transform(projection_pj, webmercator_pj, INCREMENTS, 1, &b_x[0], &b_y[0], NULL);
  pj_transform(projection_pj, webmercator_pj, INCREMENTS, 1, &l_x[0], &l_y[0], NULL);
  pj_transform(projection_pj, webmercator_pj, INCREMENTS, 1, &r_x[0], &r_y[0], NULL);
  pj_free(projection_pj);

  // Get WebMercator bounding box
  xmin = std::min(*std::min_element(std::begin(r_x), std::end(r_x)),
                  std::min(*std::min_element(std::begin(l_x), std::end(l_x)),
                           std::min(*std::min_element(std::begin(b_x), std::end(b_x)),
                                    *std::min_element(std::begin(t_x), std::end(t_x)))));
  xmax = std::max(*std::max_element(std::begin(r_x), std::end(r_x)),
                  std::max(*std::max_element(std::begin(l_x), std::end(l_x)),
                           std::max(*std::max_element(std::begin(b_x), std::end(b_x)),
                                    *std::max_element(std::begin(t_x), std::end(t_x)))));
  ymin = std::min(*std::min_element(std::begin(r_y), std::end(r_y)),
                  std::min(*std::min_element(std::begin(l_y), std::end(l_y)),
                           std::min(*std::min_element(std::begin(b_y), std::end(b_y)),
                                    *std::min_element(std::begin(t_y), std::end(t_y)))));
  ymax = std::max(*std::max_element(std::begin(r_y), std::end(r_y)),
                  std::max(*std::max_element(std::begin(l_y), std::end(l_y)),
                           std::max(*std::max_element(std::begin(b_y), std::end(b_y)),
                                    *std::max_element(std::begin(t_y), std::end(t_y)))));

  return box_t(point_t(xmin, ymin), point_t(xmax, ymax));
}

void metadata(const std::string & prefix, struct lesser_landsat_scene_struct & s, int verbose)
{
  char * wkt = NULL, * proj4 = NULL;
  OGRSpatialReferenceH srs = NULL;
  GDALDatasetH handles[3];
  char pattern[STRING_LEN];
  char filename[STRING_LEN];

  s.width = s.height = 0;

  // Open the bands
  sprintf(pattern, "%s%s", prefix.c_str(), s.filename);
  for (int i = 0; i < 3; ++ i) {
    sprintf(filename, pattern, 4-i);
    if (verbose) fprintf(stderr, ANSI_COLOR_YELLOW "%s" ANSI_COLOR_RESET "\n", filename);
    for (int j = 0; ((handles[i] = GDALOpen(filename, GA_ReadOnly)) == NULL); ++j) {
      if (j >= RETRIES) {
        fprintf(stderr, ANSI_COLOR_RED "Failed: %s:%d (handle)" ANSI_COLOR_RESET "\n", s.filename, i);
        s.width = s.height = BAD;
        break;
      }
      else {
        fprintf(stderr, ANSI_COLOR_RED "Retrying: %s:%d (handle)" ANSI_COLOR_RESET "\n", s.filename, i);
        sleep(3);
      }
    }
  }

  if (s.width == BAD || s.height == BAD) goto no_handles;

  // Get projection
  srs = OSRNewSpatialReference(NULL);
  // This is never freed, but freeing it after OSRImportFromWkt is not
  // valid either.  The problem seems to originate from within GDAL.
  wkt = static_cast<char *>(CPLMalloc(STRING_LEN * sizeof(char)));
  strncpy(wkt, GDALGetProjectionRef(handles[0]), STRING_LEN);
  OSRImportFromWkt(srs, &wkt);
  OSRExportToProj4(srs, &proj4);
  strncpy(s.proj4, proj4, 1<<8);

  // Get transform
  GDALGetGeoTransform(handles[0], s.transform);

  // Get dimensions
  s.width  = GDALGetRasterXSize(handles[0]);
  s.height = GDALGetRasterYSize(handles[0]);

  // Collect previews
  for (int i = 0; i < 3; ++i) {
    GDALRasterBandH band = GDALGetRasterBand(handles[i], 1);
    int w = GDALGetRasterXSize(handles[i]);
    int h = GDALGetRasterYSize(handles[i]);

    for (int j = 0; GDALRasterIO(band,
                                 GF_Read,
                                 0, 0, w, h,
                                 &(s.rgb[i]),
                                 SMALL_TILE_SIZE, SMALL_TILE_SIZE,
                                 GDT_UInt16, 0, 0); ++j) {
      if (j >= RETRIES) {
        fprintf(stderr, ANSI_COLOR_RED "Failed: %s:%d (preview)" ANSI_COLOR_RESET "\n", s.filename, i);
        s.width = s.height = BAD;
        break;
      }
      else {
        fprintf(stderr, ANSI_COLOR_RED "Retrying: %s:%d (preview)" ANSI_COLOR_RESET "\n", s.filename, i);
        sleep(3);
      }
    }
  }

  if (s.width == BAD || s.height == BAD) goto no_previews;

  // Calculate minima, maxima
  s.max[0] = s.max[1] = s.max[2] = 0;
  for (int j = 0; j < SMALL_TILE_SIZE; ++j) {
    for (int i = 0; i < SMALL_TILE_SIZE; ++i) {
      uint16_t r, g, b;
      r = s.rgb[0][i + j*SMALL_TILE_SIZE];
      g = s.rgb[1][i + j*SMALL_TILE_SIZE];
      b = s.rgb[2][i + j*SMALL_TILE_SIZE];
      if (r | g | b) {
        s.max[0] = std::max(s.max[0], r);
        s.max[1] = std::max(s.max[1], g);
        s.max[2] = std::max(s.max[2], b);
      }
    }
  }

  // Cleanup
 no_previews:
  CPLFree(proj4);
  OSRRelease(srs);
  for (int i = 0; i < 3; ++i)
    GDALClose(handles[i]);

 no_handles:
  return;
}

int main(int argc, const char ** argv)
{
  std::vector<value_t> scene_metadata_list;
  std::vector<const char *> scene_filename_list;
  char buffer[STRING_LEN];
  char product_id[STRING_LEN];
  char infix[STRING_LEN];
  std::string list_prefix = std::string(DEFAULT_LIST_PREFIX);
  std::string read_prefix = std::string(DEFAULT_READ_PREFIX);
  std::string stem = std::string(DEFAULT_STEM);
  int order_of_magnitude = 20;

  // Arguments from command line
  if (argc > 1) stem = std::string(argv[1]);
  if (argc > 2) sscanf(argv[2], "%d", &order_of_magnitude);
  if (argc > 3) read_prefix = std::string(argv[3]);
  if (argc > 4) list_prefix = std::string(argv[4]);
  const std::string indexfile = stem + std::string(INDEX_EXTENSION);
  const std::string bulkfile = stem + std::string(BULK_EXTENSION);

  fprintf(stderr, ANSI_COLOR_BLUE "index file \t\t =" ANSI_COLOR_GREEN " %s" ANSI_COLOR_RESET "\n", indexfile.c_str());
  fprintf(stderr, ANSI_COLOR_BLUE "bulk file \t\t =" ANSI_COLOR_GREEN " %s" ANSI_COLOR_RESET "\n", bulkfile.c_str());
  fprintf(stderr, ANSI_COLOR_BLUE "order of magnitude \t =" ANSI_COLOR_GREEN " %d" ANSI_COLOR_RESET "\n", order_of_magnitude);
  fprintf(stderr, ANSI_COLOR_BLUE "read_prefix \t\t =" ANSI_COLOR_GREEN " %s" ANSI_COLOR_RESET "\n", read_prefix.c_str());
  fprintf(stderr, ANSI_COLOR_BLUE "list_prefix \t\t =" ANSI_COLOR_GREEN " %s" ANSI_COLOR_RESET "\n", list_prefix.c_str());

  // Initialize
  webmercator_pj = pj_init_plus(webmercator);
  GDALAllRegister();

  // Read the scene list
  while (fgets(buffer, STRING_LEN, stdin) != NULL) {
    char * filename = new char[1<<8];

    sscanf(buffer, "%[^,]", product_id);
    sscanf(strstr(buffer, list_prefix.c_str()) + strlen(list_prefix.c_str()), "%s", infix);
    *(strstr(infix, POSTFIX)) = '\0';
    scene_metadata_list.push_back(std::make_pair(box_t(point_t(0, 0), point_t(1, 1)), -1));
    sprintf(filename, "%s%s_B%%d.TIF", infix, product_id);
    scene_filename_list.push_back(filename);
  }
  fprintf(stderr, ANSI_COLOR_BLUE "scenes \t\t\t =" ANSI_COLOR_GREEN " %ld" ANSI_COLOR_RESET "\n", scene_metadata_list.size());

  // Persistence
  int fd = open(bulkfile.c_str(), O_CREAT|O_WRONLY, S_IRUSR|S_IWUSR|S_IRGRP|S_IWGRP|S_IROTH|S_IWOTH);
  if (fd == -1) {
    fprintf(stderr, ANSI_COLOR_RED "open(2) problem" ANSI_COLOR_RESET "\n");
    exit(-1);
  }
  bi::managed_mapped_file file(bi::create_only, indexfile.c_str(), 1<<order_of_magnitude);

  // Get data/metadata
  #pragma omp parallel for
  for (unsigned int i = 0; i < scene_metadata_list.size(); ++i) {
    lesser_landsat_scene_struct scene;

    // Read and write scene information
    memcpy(&scene.filename, scene_filename_list[i], sizeof(scene.filename));
    delete scene_filename_list[i];
    metadata(read_prefix, scene, true);
    if (lseek(fd, sizeof(scene)*i, SEEK_SET) == -1) {
      fprintf(stderr, ANSI_COLOR_RED "lseek64(2) problem" ANSI_COLOR_RESET "\n");
      exit(-1);
    }
    if (write(fd, &scene, sizeof(scene)) == -1) {
      fprintf(stderr, ANSI_COLOR_RED "write(2) problem" ANSI_COLOR_RESET "\n");
      exit(-1);
    }

    // Record scene bouinding box information
    scene_metadata_list[i].first = bounding_box(scene);
    scene_metadata_list[i].second = static_cast<uint64_t>(i);
  }

  close(fd);

  // Build R-Tree
  // Resource: http://www.boost.org/doc/libs/1_65_1/libs/geometry/doc/html/geometry/spatial_indexes/rtree_examples/index_stored_in_mapped_file_using_boost_interprocess.html
  allocator_t alloc(file.get_segment_manager());
  rtree_t * rtree_ptr = file.find_or_construct<rtree_t>("rtree")(params_t(), indexable_t(), equal_to_t(), alloc);
  rtree_ptr->insert(scene_metadata_list);
  bi::managed_mapped_file::shrink_to_fit(indexfile.c_str());

  return 0;
}
