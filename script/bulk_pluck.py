#!/usr/bin/env python

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
