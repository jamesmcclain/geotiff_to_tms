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

#ifndef __LANDSAT_SCENE_HANDLES_HPP__
#define __LANDSAT_SCENE_HANDLES_HPP__

#include "lesser_landsat_scene.h"
#include "gdal.h"
#include "cpl_conv.h"
#include "ogr_srs_api.h"
#include "proj_api.h"


// Resource: https://en.wikipedia.org/wiki/Rule_of_three_(C++_programming)
class landsat_scene_handles {

 public:
  landsat_scene_handles() :
    r(nullptr), g(nullptr), b(nullptr), p(nullptr) {}

  landsat_scene_handles(const landsat_scene_handles & other) :
    r(other.r), g(other.g), b(other.b), p(other.p) {}

  landsat_scene_handles(landsat_scene_handles && other) noexcept :
    r(other.r), g(other.g), b(other.b), p(other.p)
  {
    other.r = other.g = other.b = other.p = nullptr;
  }

  landsat_scene_handles & operator=(const landsat_scene_handles & other)
  {
    landsat_scene_handles tmp(other);
    *this = std::move(tmp);
    return *this;
  }

  landsat_scene_handles & operator=(landsat_scene_handles && other) noexcept
  {
    pj_free(p);
    GDALClose(b);
    GDALClose(g);
    GDALClose(r);

    r = other.r;
    g = other.g;
    b = other.b;
    p = other.p;

    other.r = other.g = other.b = other.p = nullptr;

    return *this;
  }

  landsat_scene_handles(GDALDatasetH red, GDALDatasetH green, GDALDatasetH blue, projPJ proj) :
    r(red), g(green), b(blue), p(proj) {}

  ~landsat_scene_handles() noexcept
  {
    pj_free(p);
    GDALClose(b);
    GDALClose(g);
    GDALClose(r);
  }

  bool operator!=(const landsat_scene_handles & other) const
  {
    return ((r != other.r) || (g != other.g) || (b != other.b) || (p != other.p));
  }

private:
  GDALDatasetH r;
  GDALDatasetH g;
  GDALDatasetH b;
  projPJ p;

};

#endif
