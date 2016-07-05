#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import functools
import glob
import multiprocessing
import os
import shutil
import subprocess
import sys
import time

import numpy as np

import checksum
import io_tools

REF_CHECKSUMS = checksum.load_ref_checksums()


def check_dictionary(actual, expected, rtol):
    """Check each value in a dictionary.

    Parameters
    ----------
    actual: dict
        Dictionary of actual values.
    expected: dict
        Dictionary of expected (reference) values to be tested against.
    rtol: float
        Relative tolerance of the float comparisons, see :func:`numpy.allclose`.
    Returns
    -------
    dict
        Only contains keys of the actual dictionary that include mismatches
        with the expected values.
    """
    result = dict()
    for key, value in actual.items():
        try:
            r = check_value(value, expected[key], rtol)
            if r is not None:
                result[key] = r
        except KeyError:
            result[key] = 'Reference not found!'

    return result or None


def check_value(actual, expected, rtol):
    """Test a value with with an expected value using approximate float tests.

    Parameters
    ----------
    actual: str, int, float, dict, list, tuple
        Value to be tested.
    expected: : str, int, float, dict, list, tuple
        Value to test against (reference). Type should match the acutal type.
    rtol: float
        Relative tolerance of the float comparisons, see :func:`numpy.allclose`.
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
        if not np.isclose(actual, expected, rtol=rtol):
            return actual, expected
    elif isinstance(actual, dict):
        return check_dictionary(actual, expected, rtol)
    elif isinstance(actual, (list, tuple)):
        if isinstance(actual[0], float):
            # Mask NaN values
            actual = np.asarray(actual)
            expected = np.asarray(expected)
            mask = np.isfinite(actual)
            if (np.any(mask) and
                    not np.allclose(actual[mask], expected[mask], rtol=rtol)):
                return actual.tolist(), expected.tolist()
        else:
            results = []
            for ae in zip(actual, expected):
                r = check_value(ae[0], ae[1], rtol)
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


def iter_cases(path):
    pattern = os.path.join(path, '*', '*')
    for dirpath in sorted(glob.glob(pattern)):
        if not os.path.isdir(dirpath):
            continue
        name = os.path.relpath(dirpath, path).replace(os.sep, '/')
        # Skip long cases if not required
        yield name


def run_haz(path, haz_bin):
    fname = glob.glob(os.path.join(path, 'Run_*'))[0]
    print('Running HAZ on:', fname)
    dirname, basename = os.path.split(fname)
    haz_bin = os.path.abspath(haz_bin)
    with open(os.devnull, 'w') as fp:
        p = subprocess.Popen([haz_bin],
                             stdout=fp,
                             stdin=subprocess.PIPE,
                             cwd=dirname)
        # Process one file
        p.communicate(bytes('0\n0\n%s\n' % basename, 'ascii'))
        p.wait()


def test_set(name, all_cases, force=True, haz_bin='HAZ',
             root_src='', root_test='', rtol=1E-3):
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

    # fixme check this
    # if not REF_CHECKSUMS[name] == checksum.hash_directory(dirpath):
    #     logging.critical('Reference case has changed!')
    #     assert False
    ok = True
    path_test = os.path.join(root_test, name)
    if not os.path.exists(path_test) or force:
        if name in long_cases and not all_cases:
            return True

        try:
            shutil.rmtree(path_test)
            # Wait for my slow computer :-/
            time.sleep(1)
        except FileNotFoundError:
            pass

        shutil.copytree(os.path.join(root_src, name, 'Input'),
                        path_test)
        run_haz(path_test, haz_bin)

    for fname in glob.glob(os.path.join(path_test, '*')):
        ext = os.path.splitext(fname)[1]
        fname_expected = os.path.join(root_src, name, 'Output',
                                      os.path.basename(fname))
        if ext == '.out3':
            expected = io_tools.read_out3(fname_expected)
            actual = io_tools.read_out3(fname)
        elif ext == '.out4':
            expected = io_tools.read_out4(fname_expected)
            actual = io_tools.read_out4(fname)
        else:
            continue

        # Check for errors
        errors = check_value(actual, expected, rtol)
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
    parser.add_argument('-r', '--rtol', type=float, default=1E-2,
                        help='Relative tolerance used for float comparisons.')
    parser.add_argument('-s', '--root_src', type=str,
                        default='../PEER Verification Tests/',
                        help='Root path of test cases')
    parser.add_argument('-t', '--root_test', type=str,
                        default='../tests',
                        help='Root path used in testing; created if needed.')
    args = parser.parse_args()

    processes = min(max(1, args.cores), multiprocessing.cpu_count())
    with multiprocessing.Pool(processes) as pool:
        results = pool.map_async(
            functools.partial(test_set,
                              force=args.force, all_cases=args.all_cases,
                              root_src=args.root_src, root_test=args.root_test,
                              haz_bin=args.haz_bin, rtol=args.rtol),
            iter_cases(args.root_src)
        )
        ok = all(results.get())
    # Single thread
    # ok = True
    # for name in iter_cases(args.root_src):
    #     ok &= test_set(name, force=args.force, all_cases=args.all_cases,
    #              root_src=args.root_src, root_test=args.root_test,
    #              haz_bin=args.haz_bin, rtol=args.rtol)
    sys.exit(0 if ok else 1)
