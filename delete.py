#!/usr/bin/env python
# -*- coding:utf-8 -*-
from __future__ import division  # very important !!!!
import os
import gzip
import sys
import re

out = open('chr1A.1.raw', 'w')
with open('/home/wangzh/share/annotation/IWGSC_annotation_blastp', 'r') as mydata:
    for line in mydata.readlines()[1:]:
        line = line.strip().split('\t')
        line1 = line[0].strip('"')
        out.write(line1 + "\n")
out.close()


