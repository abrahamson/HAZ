#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import datetime
import glob
import os

import checksum


PATH = '../PEER Verification Tests'
PATTERN = os.path.join(PATH, '*', '*')
fp = open('ref_checksums.txt', 'wt')
fp.write('# Created: {}\n'.format(datetime.datetime.now()))

for dirpath in glob.glob(PATTERN):
    if not os.path.isdir(dirpath):
        continue
    digest = checksum.hash_directory(dirpath)
    fp.write('\t'.join([
        os.path.relpath(dirpath, PATH),
        digest
    ]) + '\n')

fp.close()
