#!/usr/bin/env python

"""
Copyright (c) 2017, James McClain
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the
   distribution.
3. All advertising materials mentioning features or use of this
   software must display the following acknowledgement: This product
   includes software developed by Dr. James W. McClain.
4. Neither the names of the authors nor the names of the
   contributors may be used to endorse or promote products derived
   from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE AUTHORS ''AS IS'' AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

import sys
import csv
import re
import gzip


# bulk_pluck.py metadata.csv scene_list.gz
# metadata.csv comes from https://landsat.usgs.gov/landsat-bulk-metadata-service
# scene_list.gz comes from http://landsat-pds.s3.amazonaws.com/c1/L8/scene_list.gz

scenes = {}

def open2(filename):
    gzipped = re.compile("\.gz$")
    if (gzipped.search(filename)):
        return gzip.open(filename)
    else:
        return open(filename)

with open2(sys.argv[1]) as f:
    reader = csv.DictReader(f, delimiter=',')

    for row in reader:
        product_id = row['LANDSAT_PRODUCT_ID']
        cloud_cover = float(row['CLOUD_COVER_LAND'])
        is_day = (row['dayOrNight'] == 'DAY')
        if (0.0 <= cloud_cover and cloud_cover <= 10.0 and is_day):
            pair = (row['row'], row['path'])
            if pair in scenes:
                elevation_1 = float(row['sunElevation'])
                elevation_2 = float(scenes[pair]['sunElevation'])
                if (elevation_1 > elevation_2):
                    scenes[pair] = row
            else:
                scenes[pair] = row

desired = set([row['LANDSAT_PRODUCT_ID'] for row in scenes.values()])

with open2(sys.argv[2]) as f:
    reader = csv.DictReader(f, delimiter=',')
    writer = csv.DictWriter(sys.stdout, fieldnames=reader.fieldnames)

    for csv_row in reader:
        if csv_row['productId'] in desired:
            writer.writerow(csv_row)
