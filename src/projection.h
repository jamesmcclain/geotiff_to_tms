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

#ifndef __PROJECTION_H__
#define __PROJECTION_H__


#ifdef __cplusplus
extern "C" {
#endif

  void inline image_to_world(double * __restrict__ x, double * __restrict__ y, const double * __restrict__ transform)
  {
    // Source: http://www.gdal.org/classGDALDataset.html#a5101119705f5fa2bc1344ab26f66fd1d
    double image_x = *x, image_y = *y;
    *x = transform[0] + image_x*transform[1] + image_y*transform[2];
    *y = transform[3] + image_x*transform[4] + image_y*transform[5];
  }

  void inline world_to_image(double * __restrict__ x, double * __restrict__ y, const double * __restrict__ transform)
  {
    // Source: http://www.gdal.org/classGDALDataset.html#a5101119705f5fa2bc1344ab26f66fd1d
    double world_x = *x, world_y = *y;
    *x = (-world_x*transform[5] + world_y*transform[2] + transform[0]*transform[5] - transform[2]*transform[3])/(transform[2]*transform[4] - transform[1]*transform[5]);
    *y = ( world_x*transform[4] - world_y*transform[1] + transform[0]*transform[4] + transform[1]*transform[3])/(transform[2]*transform[4] - transform[1]*transform[5]);
  }

#ifdef __cplusplus
}
#endif

#endif
