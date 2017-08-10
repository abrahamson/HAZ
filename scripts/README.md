# Python Scripts Supporting HAZ

This folder contains scripts for reading HAZ files, and running tests on the
`HAZ`.

## Requirements

The scripts are written for Python 3.5, and currently do not support
earlier versions. Additionally, the `numpy` library is required for the
numeric tests.  The Python installer is available from
[here](https://www.python.org/downloads/). After Python is installed, the
`numpy` library may be installed using `pip`:

    pip install numpy

Alternatively, Python and the `numpy` library can be installed via the `conda`
package manager. After installing
[miniconda3](http://conda.pydata.org/miniconda.html), the `numpy` library can
be installed with:

    conda install numpy

## Overview

The main program is `run_tests.py`. The commands line options can be printed
with:

    $> python run_tests.py --help
    usage: run_tests.py [-h] [-a] [-b HAZ_BIN] [-c CORES] [-f] [-r RTOL]
                        [-s ROOT_SRC] [-t ROOT_TEST]

    Perform and report tests for HAZ PSHA code.

    positional arguments:
      patterns              Only process test cases matching the pattern 
                            (default: None)

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
      -s ROOT_REF, --root_ref ROOT_REF
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

## Test Duration

During the testing process, the time for running `HAZ` on a given input file is
reported. If the test output is redirected to a file, the test duration can be
extracted as follows:

    $> python run_tests.py -a -c 2 -b ../build/HAZ -r 0.005 > test.log
    $> cat haz_tests.log | sed -En 's/.*: (\S+) (\S+)$/\1\t\2/p' | sort

    Set1/S1Test01           0:00:00.938700
    Set1/S1Test02           0:00:33.722160
    Set1/S1Test04           0:00:45.615901
    Set1/S1Test05           1:15:58.920194
    Set1/S1Test06           1:18:25.974351
    Set1/S1Test07           1:18:23.266588
    Set1/S1Test08a          0:00:33.947122
    Set1/S1Test08b          0:00:35.287663
    Set1/S1Test08c          0:00:34.911442
    Set1/S1Test10           0:00:05.554910
    Set1/S1Test11           0:00:24.390950
    Set2/S2Test1            0:02:35.734937
    Set2/S2Test2a           1:13:19.895491
    Set2/S2Test2b           1:21:30.252216
    Set2/S2Test2c           1:19:56.661891
    Set2/S2Test2d           1:13:25.263971
    Set2/S2Test3a           0:00:07.209199
    Set2/S2Test3b           0:00:07.187167
    Set2/S2Test3c           0:00:10.300279
    Set2/S2Test3d           0:00:07.637015
    Set2/S2Test4a           0:00:22.463504
    Set2/S2Test4b           0:00:29.631422
    Set2/S2Test5a           0:00:05.331427
    Set2/S2Test5b           0:00:04.926723
    Set3/S3Test1a           0:00:01.603328
    Set3/S3Test1b           0:00:01.516753
    Set3/S3Test2/Fractiles  0:07:55.357953
    Set3/S3Test2/Hazard     0:08:17.658489

The calculation time ranges from less than a second to almost to 90 minutes for
this test computer.
