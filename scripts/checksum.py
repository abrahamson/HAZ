#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import hashlib
import os


def hash_directory(path):
    """Create a SHA-1 hash including all files in a directory.
    """
    m = hashlib.sha1()

    for root, _, fnames in os.walk(path):
        for fname in fnames:
            fpath = os.path.join(root, fname)
            size = os.path.getsize(fpath)
            with open(fpath, 'rb') as fp:
                while fp.tell() != size:
                    # Read 256k at a time
                    m.update(fp.read(0x40000))
    return m.hexdigest()


def load_ref_checksums():
    """Load reference checksums.
    """
    with open('ref_checksums.txt') as fp:
        next(fp)
        checksums = dict(line.strip().split('\t') for line in fp)

    return checksums
