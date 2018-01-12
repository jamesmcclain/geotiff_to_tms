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

#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__

#define BAD (0xBAADF00D)
#define BULK_EXTENSION ".bulk"
#define DEFAULT_LIST_PREFIX "https://s3-us-west-2.amazonaws.com/landsat-pds/"
#define DEFAULT_READ_PREFIX "/vsicurl/https://s3-us-west-2.amazonaws.com/landsat-pds/"
#define DEFAULT_STEM "/tmp/landsat"
#define INDEX_EXTENSION ".index"
#define RADIUS (6378137.0)
#define RETRIES (3)
#define SMALL_TILE_SIZE2 (SMALL_TILE_SIZE * SMALL_TILE_SIZE)
#define SMALL_TILE_SIZE (int)(1.414213562373095048801688724209698078569671875376948073176*(1<<SMALL_TILE_ZOOM))
#define SMALL_TILE_ZOOM (8)
#define STRING_BUFFER_SIZE (1<<10)
#define STRING_LEN (1<<10)
#define TILE_SIZE (1<<8)
#define TILE_SIZE2 (TILE_SIZE * TILE_SIZE)
#define WEBMERCATOR "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext +no_defs"

#endif
