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
#include <cstring>
#include <cmath>

#include <boost/interprocess/managed_mapped_file.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

#include "gdal.h"
#include "cpl_conv.h"
#include "ogr_srs_api.h"
#include "proj_api.h"

#include "ansi.h"
#include "landsat_scene.h"
#include "projection.h"

#define MAX_LEN (1<<10)
#define INCREMENTS (1<<8)
#define PREFIX "https://s3-us-west-2.amazonaws.com/landsat-pds/"
#define POSTFIX "index.html"
#define DEFAULT_PREFIX "/vsicurl/https://s3-us-west-2.amazonaws.com/landsat-pds/"
#define DEFAULT_INDEXFILE "/tmp/index.data"

const char * webmercator = "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext +no_defs";
projPJ webmercator_pj = NULL;

namespace bi  = boost::interprocess;
namespace bg  = boost::geometry;
namespace bgm = boost::geometry::model;
namespace bgi = boost::geometry::index;

// Resource: http://www.boost.org/doc/libs/1_63_0/libs/geometry/doc/html/geometry/spatial_indexes/rtree_examples/index_stored_in_mapped_file_using_boost_interprocess.html
typedef bgm::point<double, 2, bg::cs::cartesian> point_t;
typedef bgm::box<point_t> box_t;
typedef std::pair<box_t, landsat_scene> value_t;
typedef bgi::linear<32,8> params_t;
typedef bgi::indexable<value_t> indexable_t;
typedef bgi::equal_to<value_t> equal_to_t;
typedef bi::allocator<value_t, bi::managed_mapped_file::segment_manager> allocator_t;
typedef bgi::rtree<value_t, params_t, indexable_t, equal_to_t, allocator_t> rtree_t;


void bounding_box(std::pair<box_t, landsat_scene> & pair)
{
  double xmin = DBL_MAX;
  double ymin = DBL_MAX;
  double xmax = DBL_MIN;
  double ymax = DBL_MIN;
  double * t = static_cast<double *>(calloc(INCREMENTS, 2*sizeof(double)));
  double * b = static_cast<double *>(calloc(INCREMENTS, 2*sizeof(double)));
  double * l = static_cast<double *>(calloc(INCREMENTS, 2*sizeof(double)));
  double * r = static_cast<double *>(calloc(INCREMENTS, 2*sizeof(double)));

  // Get world coordinates around the periphery
  for (int i = 0; i < INCREMENTS; ++i) {
    t[2*i + 0] = b[2*i + 0] = 0.5 + static_cast<double>(i*pair.second.width)/INCREMENTS;  // x values across the top and bottom
    t[2*i + 1] = l[2*i + 0] = 0.0;                                                        // y values across the top and x values on the left// XXX
    b[2*i + 1] = pair.second.height;                                                      // y values across the bottom
    r[2*i + 0] = pair.second.width;                                                       // x values on the right
    l[2*i + 1] = r[2*i + 1] = 0.5 + static_cast<double>(i*pair.second.height)/INCREMENTS; // y values on the left and right

    image_to_world(t + 2*i, pair.second.transform);
    image_to_world(b + 2*i, pair.second.transform);
    image_to_world(l + 2*i, pair.second.transform);
    image_to_world(r + 2*i, pair.second.transform);
  }

  // Get WebMercator coordinates
  pj_transform(webmercator_pj, pair.second.projection, INCREMENTS, 2, t+0, t+1, NULL);
  pj_transform(webmercator_pj, pair.second.projection, INCREMENTS, 2, b+0, b+1, NULL);
  pj_transform(webmercator_pj, pair.second.projection, INCREMENTS, 2, l+0, l+1, NULL);
  pj_transform(webmercator_pj, pair.second.projection, INCREMENTS, 2, r+0, r+1, NULL);

  // Get WebMercator bounding box
  for (int i = 0; i < INCREMENTS; ++i) {
    xmin = fmin(r[2*i + 0], fmin(l[2*i + 0], fmin(b[2*i + 0], fmin(t[2*i + 0], xmin))));
    xmax = fmax(r[2*i + 0], fmax(l[2*i + 0], fmax(b[2*i + 0], fmax(t[2*i + 0], xmax))));
    ymin = fmin(r[2*i + 1], fmin(l[2*i + 1], fmin(b[2*i + 1], fmin(t[2*i + 1], ymin))));
    ymax = fmax(r[2*i + 1], fmax(l[2*i + 1], fmax(b[2*i + 1], fmax(t[2*i + 1], ymax))));
  }
  pair.first = box_t(point_t(xmin, ymin), point_t(xmax, ymax));

  free(r);
  free(l);
  free(b);
  free(t);
}

void metadata(const char * prefix, struct landsat_scene & s, int verbose)
{
  char * wkt = NULL, * proj4 = NULL;
  OGRSpatialReferenceH srs = NULL;
  GDALDatasetH dataset;
  char pattern[MAX_LEN];
  char filename[MAX_LEN];

  // Open the red band of the scene
  sprintf(pattern, "%s%s", prefix, s.filename);
  sprintf(filename, pattern, 4); // 4 stands for red
  if ((dataset = GDALOpen(filename, GA_ReadOnly)) == NULL) exit(-1); // XXX
  if (verbose) fprintf(stderr, ANSI_COLOR_YELLOW "%s" ANSI_COLOR_RESET "\n", filename);

  // Get projection
  srs = OSRNewSpatialReference(NULL);
  // This is never freed, but freeing wkt is not valid either.  The
  // problem seems to originate from within GDAL.
  wkt = static_cast<char *>(CPLMalloc(MAX_LEN * sizeof(char)));
  strncpy(wkt, GDALGetProjectionRef(dataset), MAX_LEN);
  OSRImportFromWkt(srs, &wkt);
  OSRExportToProj4(srs, &proj4);
  s.projection = pj_init_plus(proj4);

  // Get transform
  GDALGetGeoTransform(dataset, s.transform);

  // Get dimensions
  s.width  = GDALGetRasterXSize(dataset);
  s.height = GDALGetRasterYSize(dataset);

  // Cleanup
  CPLFree(proj4);
  OSRRelease(srs);
  GDALClose(dataset);
}

int main(int argc, const char ** argv)
{
  std::vector<value_t> scene_list;
  char buffer[MAX_LEN];
  char product_id[MAX_LEN];
  char infix[MAX_LEN];
  const char * prefix = DEFAULT_PREFIX;
  const char * indexfile = DEFAULT_INDEXFILE;
  int order_of_magnitude = 28;

  // Arguments from command line
  if (argc > 1) indexfile = argv[1];
  if (argc > 2) sscanf(argv[2], "%d", &order_of_magnitude);
  if (argc > 3) prefix = argv[3];
  fprintf(stderr, ANSI_COLOR_BLUE "index file \t\t =" ANSI_COLOR_GREEN " %s" ANSI_COLOR_RESET "\n", indexfile);
  fprintf(stderr, ANSI_COLOR_BLUE "order of magnitude \t =" ANSI_COLOR_GREEN " %d" ANSI_COLOR_RESET "\n", order_of_magnitude);
  fprintf(stderr, ANSI_COLOR_BLUE "prefix \t\t\t =" ANSI_COLOR_GREEN " %s" ANSI_COLOR_RESET "\n", prefix);

  // Initialize
  webmercator_pj = pj_init_plus(webmercator);
  GDALAllRegister();

  // Read the scene list
  while (fgets(buffer, MAX_LEN, stdin) != NULL) {
    sscanf(buffer, "%[^,]", product_id);
    sscanf(strstr(buffer, PREFIX), PREFIX "%s", infix);
    *(strstr(infix, POSTFIX)) = '\0';
    scene_list.push_back(std::make_pair(box_t(point_t(0, 0), point_t(1, 1)), landsat_scene()));
    sprintf(scene_list.back().second.filename, "%s%s_B%%d.TIF", infix, product_id);
  }
  fprintf(stderr, ANSI_COLOR_BLUE "scenes \t\t\t =" ANSI_COLOR_GREEN " %ld" ANSI_COLOR_RESET "\n", scene_list.size());

  // Get metadata
  #pragma omp parallel for
  for (unsigned int i = 0; i < scene_list.size(); ++i) {
    metadata(prefix, scene_list[i].second, true);
    bounding_box(scene_list[i]);
  }

  // Build R-Tree
  bi::managed_mapped_file file(bi::create_only, indexfile, 1<<order_of_magnitude);
  allocator_t alloc(file.get_segment_manager());
  rtree_t * rtree_ptr = file.find_or_construct<rtree_t>("rtree")(params_t(), indexable_t(), equal_to_t(), alloc);
  rtree_ptr->insert(scene_list);

  return 0;
}
