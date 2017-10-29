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

#include <boost/interprocess/managed_mapped_file.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>

#include "landsat_scene.h"

#define MAX_LEN (1<<10)
#define PREFIX "https://s3-us-west-2.amazonaws.com/landsat-pds/"
#define POSTFIX "index.html"
#define DEFAULT "/vsicurl/https://s3-us-west-2.amazonaws.com/landsat-pds/"


int main(int argc, const char ** argv)
{
  std::vector<landsat_scene> scene_list;
  char buffer[MAX_LEN];
  char product_id[MAX_LEN];
  char infix[MAX_LEN];
  const char * prefix = DEFAULT;

  // Prefix from command line
  if (argc > 1) prefix = argv[1];

  while (fgets(buffer, MAX_LEN, stdin) != NULL) {
    struct landsat_scene scene;

    sscanf(buffer, "%[^,]", product_id);
    sscanf(strstr(buffer, PREFIX), PREFIX "%s", infix);
    *(strstr(infix, POSTFIX)) = '\0';

    sprintf(scene.filename, "%s%s%s_B%%d.TIF", prefix, infix, product_id);
    scene_list.push_back(scene);
  }

  fprintf(stdout, "%ld\n", scene_list.size());
  for (auto i = scene_list.begin(); i != scene_list.end(); ++i) {
    fprintf(stdout, "%s\n", i->filename);
  }

  return 0;
}
