# Python Scripts Supporting HAZ

This folder contains scripts for reading HAZ files, and running tests on the
`HAZ`.

## Requirements

The scripts are written for Python 3.5, and currently does not support
previous versions. Additionally, `numpy` is required for the numeric tests.
The Python installer is available from
[here](https://www.python.org/downloads/). After Python is installed, `numpy`
may be installed using `pip`:

    pip install numpy

Alternatively, Python and `numpy` can be installed via `conda`. After
installing [miniconda3](http://conda.pydata.org/miniconda.html), the `numpy`
library can be installed with:

    conda install numpy

## Overview

The main program is `run_tests.py`. The commands line options can be printed
with:

    $> python run_tests.py --help
    usage: run_tests.py [-h] [-a] [-b HAZ_BIN] [-c CORES] [-f] [-r RTOL]
                        [-s ROOT_SRC] [-t ROOT_TEST]

    Perform and report tests for HAZ PSHA code.

    optional arguments:
      -h, --help            show this help message and exit
      -a, --all_cases       Perform all test cases, which might take many hours.
                            (default: False)
      -b HAZ_BIN, --haz_bin HAZ_BIN
                            Name of HAZ binary. Path is required if not in PATH
                            variable. (default: ../build/HAZ.exe)
      -c CORES, --cores CORES
                            Number of cores to use. (default: 1)
      -f, --force           Force HAZ to rerun; otherwise it only runs if test
                            case directory is empty. (default: False)
      -r RTOL, --rtol RTOL  Relative tolerance used for float comparisons.
                            (default: 0.002)
      -s ROOT_SRC, --root_src ROOT_SRC
                            Root path of test cases (default:
                            ../PEER_Verification_Tests/)
      -t ROOT_TEST, --root_test ROOT_TEST
                            Root path used in testing; created if needed.
                            (default: ../tests)

Depending the on configuration of the system calling the Python 3.5
interpreter may require changing how `run_tests.py` is invoked.

The standard tests can then be run with:

    $> python run_tests.py --cores 2 --haz_bin ../build/HAZ --rtol 0.005
    Running HAZ on: ../tests\Set1/S1Test05_fast\Run_S1Test5_fast.txt
    Running HAZ on: ../tests\Set1/S1Test06_fast\Run_S1Test6_fast.txt
    Running HAZ on: ../tests\Set1/S1Test07_fast\Run_S1Test7_fast.txt
    Running HAZ on: ../tests\Set1/S1Test08a\Run_S1Test8a.txt
    Running HAZ on: ../tests\Set1/S1Test08b\Run_S1Test8b.txt
    Running HAZ on: ../tests\Set1/S1Test08c\Run_S1Test8c.txt
    Running HAZ on: ../tests\Set2/S2Test2a_fast\Run_S2Test2a_fast.txt
    Running HAZ on: ../tests\Set2/S2Test2b_fast\Run_S2Test2b_fast.txt
    Running HAZ on: ../tests\Set2/S2Test3b\Run_S2Test3b.txt
    Running HAZ on: ../tests\Set2/S2Test3c\Run_S2Test3c.txt
    Running HAZ on: ../tests\Set2/S2Test3d\Run_S2Test3d.txt
    Running HAZ on: ../tests\Set2/S2Test4a\Run_S2Test4a.txt
    Running HAZ on: ../tests\Set1/S1Test01\Run_S1Test1.txt
    Running HAZ on: ../tests\Set1/S1Test02\Run_S1Test2.txt
    Running HAZ on: ../tests\Set1/S1Test04\Run_S1Test4.txt
    Running HAZ on: ../tests\Set1/S1Test10\Run_S1Test10.txt
    Running HAZ on: ../tests\Set1/S1Test11\Run_S1Test11.txt
    Running HAZ on: ../tests\Set2/S2Test1\Run_S2Test1.txt
    Running HAZ on: ../tests\Set2/S2Test2c_fast\Run_S2Test2c_fast.txt
    Running HAZ on: ../tests\Set2/S2Test2d_fast\Run_S2Test2d_fast.txt
    Running HAZ on: ../tests\Set2/S2Test3a\Run_S2Test3a.txt
    Running HAZ on: ../tests\Set2/S2Test4b\Run_S2Test4b.txt
    Running HAZ on: ../tests\Set2/S2Test5a\Run_S2Test5a.txt
    Running HAZ on: ../tests\Set2/S2Test5b\Run_S2Test5b.txt

By default, only the tests with fast calculation time are performed. All
tests, can be run with the `-a` or `--all_tests` flag, such as:

    $> python run_tests.py -a -c 2 -b ../build/HAZ -r 0.005
