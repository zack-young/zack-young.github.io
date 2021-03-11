#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@author:Levy
@file:finalgeno.py
@time:2018/7/2417:35
"""

from __future__ import division # very important !!!!
from collections import defaultdict
import os
import gzip
import sys
import re
import math

infile = open("/data/user/yangzz/mapping/S15/ReseqPipeNew/09.filterVCF/GT_AD/s11_AD/chr1A.gt", "r")
PART = infile.split("_")
part = PART[0]
art = part.lstrip('s')
print(art)
infile.close
