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

scene_file = sys.argv[2]
desired = set([])
seen = set([])

with open(sys.argv[1]) as f:
    reader = csv.DictReader(f, delimiter=',')
    for row in reader:
        product_id = row['LANDSAT_PRODUCT_ID']
        desired.add(product_id)

with gzip.open(sys.argv[2]) as f:
    product_id_expr = re.compile('^([^,]*)')
    path_col_expr = re.compile('c1/L8/(\d+)/(\d+)/')
    for line in f:
        line = line[0:len(line)-1]
        match = re.search(product_id_expr, line)
        if (match != None) and (match.group(1) in desired):
            match = re.search(path_col_expr, line)
            path = match.group(1)
            row = match.group(2)
            pair = (path, row)
            if pair not in seen:
                seen.add(pair)
                print(line)
