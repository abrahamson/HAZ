#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import glob
import os

import checksum

PATH_SRC = '../PEER Verification Tests/'
PATH_TEST = '../tests'

REF_CHECKSUMS = checksum.load_ref_checksums()


def test_sets():
    pattern = os.path.join(PATH_SRC, '*', '*')
    for dirpath in glob.glob(pattern):
        if not os.path.isdir(dirpath):
            continue


        try:
            os.maked


