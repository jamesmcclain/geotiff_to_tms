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
#include "ansi.h"
#include "gdal.h"
#include "cpl_conv.h"
#include "ogr_srs_api.h"

#define BUFFERSIZE (1<<10)

void init()
{
  GDALAllRegister();
}

void load()
{
  const char * filename = "/tmp/LC08_L1TP_139045_20170304_20170316_01_T1_B1.TIF";
  char * wkt;
  char * proj4;
  GDALDatasetH * dataset = NULL;
  double transform[6];
  OGRSpatialReferenceH srs = NULL;

  dataset = GDALOpen(filename, GA_ReadOnly);
  if(dataset == NULL) {
    fprintf(stderr, ANSI_COLOR_RED "Exiting" ANSI_COLOR_RESET);
    exit(-1);
  }

  srs = OSRNewSpatialReference(NULL);
  wkt = calloc(BUFFERSIZE, sizeof(char));
  strncpy(wkt, GDALGetProjectionRef(dataset), BUFFERSIZE);
  OSRImportFromWkt(srs, &wkt);
  OSRExportToProj4(srs, &proj4);
  fprintf(stdout, "%s\n", proj4);
  CPLFree(proj4);

  GDALGetGeoTransform(dataset, transform);
  for (int i = 0; i < 6; ++i) {
    fprintf(stdout, "%lf ", transform[i]);
  }
  fprintf(stdout, "\n");

}
