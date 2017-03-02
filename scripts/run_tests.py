#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import datetime
import functools
import glob
import multiprocessing
import pathlib
import os
import re
import shutil
import subprocess
import sys
import time

from typing import List

import numpy as np

import checksum
import io_tools


def check_dictionary(actual, expected, rtol: float, atol: float):
    """Check each value in a dictionary.

    Parameters
    ----------
    actual: dict
        Dictionary of actual values.
    expected: dict
        Dictionary of expected (reference) values to be tested against.
    rtol: float
        Relative tolerance of the float comparisons, see :func:`numpy.allclose`.
    atol: float
        Absolute tolerance of the float comparisons, see :func:`numpy.allclose`.
    Returns
    -------
    dict
        Only contains keys of the actual dictionary that include mismatches
        with the expected values.
    """
    dist_keys = ['dist_avg']
    dist_atol = 1e-03

    result = dict()
    for key, value in actual.items():
        try:
            # Check distance metrics with a large *atol* value.
            r = check_value(value, expected[key], rtol,
                            (dist_atol if key in dist_keys else atol))
            if r is not None:
                result[key] = r
        except KeyError:
            result[key] = 'Reference not found!'

    return result or None


def check_value(actual, expected, rtol, atol):
    """Test a value with with an expected value using approximate float tests.

    Parameters
    ----------
    actual: str, int, float, dict, list, tuple
        Value to be tested.
    expected: : str, int, float, dict, list, tuple
        Value to test against (reference). Type should match the acutal type.
    rtol: float
        Relative tolerance of the float comparisons, see :func:`numpy.allclose`.
    atol: float
        Absolute tolerance of the float comparisons, see :func:`numpy.allclose`.
    Returns
    -------
    None or tuple(actual, expected):
        If value agress then *None* is returned. If values disagree, then a
        tuple of the value is returned.
    """
    assert type(actual) == type(expected)

    if isinstance(actual, (str, int)):
        if actual != expected:
            return actual, expected
    elif isinstance(actual, float):
        if not np.isclose(actual, expected, rtol=rtol, atol=atol):
            return actual, expected
    elif isinstance(actual, dict):
        return check_dictionary(actual, expected, rtol, atol)
    elif isinstance(actual, (list, tuple)):
        if isinstance(actual[0], float):
            # Mask NaN values
            actual = np.asarray(actual)
            expected = np.asarray(expected)
            if not np.allclose(actual, expected, rtol=rtol, atol=atol,
                               equal_nan=True):
                return actual.tolist(), expected.tolist()
        else:
            results = []
            for ae in zip(actual, expected):
                r = check_value(ae[0], ae[1], rtol, atol)
                if r is not None:
                    results.append(r)
            return results or None
    else:
        raise NotImplementedError


def print_errors(name, errors):
    if isinstance(errors, dict):
        for k, v in errors.items():
            print_errors(name + '[%s]' % k, v)
    elif isinstance(errors, list):
        for i, e in enumerate(errors):
            print_errors(name + '[%d]' % i, e)
    elif isinstance(errors, tuple) and len(errors) == 2:
        indent = '    '
        print(indent + name)
        if isinstance(errors[0], list):
            if isinstance(errors[0][0], float):
                # List of floats
                print(indent + '    {:12s}{:12s}{:12s}'.format(
                    'Actual', 'Expected', 'Rel.Diff (%)'
                ))
                actual = np.array(errors[0])
                expected = np.array(errors[1])
                rel_diff = 100. * (actual - expected) / expected
                rows = np.c_[actual, expected, rel_diff].tolist()
                for i, row in enumerate(rows):
                    print(indent + '{:4d}{:12g}{:12g}{:12g}'.format(i, *row))
            else:
                print(indent +
                      '    {:.12s}{:.12s}'.format('Actual', 'Expected'))
                for i, ae in enumerate(zip(*errors)):
                    print(indent +
                          '{:4d}{:.12g}{:.12g}'.format(i, *ae))
        else:
            print(indent + '  Actual: {}'.format(errors[0]))
            print(indent + 'Expected: {}'.format(errors[1]))
    else:
        raise NotImplementedError


def iter_cases(path_src: str,
               patterns: List[str],
               all_cases: bool=False
               ) -> pathlib.PurePath:
    """Iterate over test cases.
    """
    # Excluded filenames
    excluded = [
        # FIXME: Currently calculation of fractiles is not supported
        'Run_Fractiles.txt',
    ]
    # These cases take a couple hours to run
    long_cases = [
        'Set1/S1Test05',
        'Set1/S1Test06',
        'Set1/S1Test07',
        'Set2/S2Test2a',
        'Set2/S2Test2b',
        'Set2/S2Test2c',
        'Set2/S2Test2d',
    ]

    pattern_run = r'^Run_\S+\.txt$'

    for root, dirnames, fnames in os.walk(path_src):
        for fname in fnames:
            if fname in excluded:
                continue
            if not re.match(pattern_run, fname):
                continue

            if patterns and not any(re.search(p, root) for p in patterns):
                continue

            path = pathlib.Path(root).parent
            # Check cases with long runtimes
            if not all_cases and any(path.match(lc) for lc in long_cases):
                print('Skipping:', path)
                continue

            yield path


def run_haz(path: pathlib.PurePath, haz_bin: str):
    # Use the first Run_ filename
    fpath = next(path.glob('Run_*.txt'))
    print('Running HAZ on:', fpath)
    haz_bin = os.path.abspath(haz_bin)
    with open(os.devnull, 'w') as fp:
        p = subprocess.Popen([haz_bin],
                             stdout=fp,
                             stdin=subprocess.PIPE,
                             cwd=str(fpath.parent))
        # Process one file
        p.communicate(bytes('0\n0\n%s\n' % fpath.name, 'ascii'))
        p.wait()


def test_path(path_ref: pathlib.PurePath,
             force: bool=True,
             haz_bin: str='HAZ',
             root_ref: str='',
             root_test: str='',
             rtol: float=1E-3) -> bool:

    print(path_ref)
    path_test = pathlib.Path(root_test, path_ref.relative_to(root_ref))
    if not path_test.exists() or force:
        try:
            shutil.rmtree(path_test)
            # Wait for my slow computer :-/
            time.sleep(1)
        except FileNotFoundError:
            pass

        # Copy files over
        shutil.copytree(path_ref.joinpath('Input'), path_test)
        # Run HAZ and track the duration
        start = datetime.datetime.now()
        run_haz(path_test, haz_bin)
        time_diff = datetime.datetime.now() - start

        print(
            'Calculation time: {} {}'.format(
                path_ref.relative_to(root_ref), time_diff)
        )

    ok = True
    for fpath_test in path_test.iterdir():
        ext = fpath_test.suffix
        fpath_ref = path_ref.joinpath('Output', fpath_test.name)

        if not fpath_ref.exists():
            continue

        if ext == '.out3':
            expected = io_tools.read_out3(str(fpath_ref))
            actual = io_tools.read_out3(str(fpath_test))
        elif ext == '.out4':
            expected = io_tools.read_out4(str(fpath_ref))
            actual = io_tools.read_out4(str(fpath_test))
        else:
            continue

        # Check for errors
        errors = check_value(actual, expected, rtol, atol=1e-08)
        ok &= (not errors)
        if errors:
            print('Errors in:', name)
            print_errors('%s: ' % fname, errors)
    return ok

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Perform and report tests for HAZ PSHA code.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('-a', '--all_cases', action='store_true',
                        help='Perform all test cases, ' +
                             'which might take many hours.')
    parser.add_argument('-b', '--haz_bin', type=str,
                        default='../build/HAZ.exe',
                        help='Name of HAZ binary. ' +
                             'Path is required if not in PATH variable.')
    parser.add_argument('-c', '--cores', type=int, default=1,
                        help='Number of cores to use.')
    parser.add_argument('-f', '--force', action='store_true',
                        help='Force HAZ to rerun; otherwise it only runs if '
                             'test case directory is empty.')
    parser.add_argument('-r', '--rtol', type=float, default=2E-3,
                        help='Relative tolerance used for float comparisons.')
    parser.add_argument('-s', '--root_ref', type=str,
                        default='../PEER_Verification_Tests/',
                        help='Root path of test cases')
    parser.add_argument('-t', '--root_test', type=str,
                        default='../tests',
                        help='Root path used in testing; created if needed.')
    parser.add_argument('patterns', default=None, nargs='*',
                        help='Only process test cases matching the pattern')
    args = parser.parse_args()

    if args.cores == 1:
        # Single thread
        ok = True
        for name in iter_cases(args.root_ref, args.patterns, args.all_cases):
            ok &= test_path(name, force=args.force, root_ref=args.root_ref,
                           root_test=args.root_test, haz_bin=args.haz_bin,
                           rtol=args.rtol)
    else:
        # Multi-threaded
        processes = min(max(1, args.cores), multiprocessing.cpu_count())

        with multiprocessing.Pool(processes) as pool:
            results = pool.map_async(
                functools.partial(test_path, force=args.force,
                                  root_ref=args.root_ref,
                                  root_test=args.root_test,
                                  haz_bin=args.haz_bin, rtol=args.rtol),
                iter_cases(args.root_ref, args.patterns, args.all_cases))
            ok = all(results.get())

    sys.exit(0 if ok else 1)
