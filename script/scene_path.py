#!/usr/bin/env python

import sys
import re


path = sys.argv[1]
expr = re.compile("c1/L8/{}/(\d+)/".format(path))
seen = set([])

for line in sys.stdin:
    line = line[0:len(line)-1]
    match = re.search(expr, line)

    if match != None:
        row = match.group(1)
        if row not in seen:
            print(line)
            seen.add(row)
