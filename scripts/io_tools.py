import re
from typing import Callable, Dict, List, Tuple, TypeVar

import numpy as np

T = TypeVar("T", dict, int, float, str)


def fixed_split(
    line: str, pairs: List[Tuple[int, Callable]], keep_tail: bool = False
) -> List[T]:
    """Read a fixed width line of text.

    Parameters
    ----------
    line : str
        text line to parse
    pairs : list[(int, func)]
        list of tuples specifying the width and the string parser
    keep_tail : bool
        remove entries without values

    Returns
    -------
    values : list
        values read
    """
    line = line.rstrip()

    i = 0
    values = []
    for w, p in pairs:
        try:
            v = p(line[i : i + w])
        except ValueError:
            v = None

        values.append(v)
        i += w

    if not keep_tail:
        while values and values[-1] is None:
            values.pop()

    return tuple(values)


def parse_blocks(lines: List[str], pattern: str) -> List[str]:
    """Separate lines into a series of blocks using an identifying pattern.

    Parameters
    ----------
    lines: list[str]
        List of strings
    pattern: str
        Regex pattern used to identify the blocks.
    Returns
    -------
    blocks : list[str]
        Parsed blocks
    """
    starts = [i for i, l in enumerate(lines) if re.search(pattern, l)]
    starts.append(len(lines))
    blocks = [lines[s : starts[i + 1]] for i, s in enumerate(starts[:-1])]
    return blocks


def parse_site_line(line: str) -> Dict[str, T]:
    """Parse the site line in out3 and out4 files.

    Parameters
    ----------
    line : str
        Line of text
    Returns
    -------
    site : dict
        Site dictionary containing keys:
            id: int
                identifier
            lon: float
                longitude (E)
            lat: float
                latitude (N)
    """
    parts = line.split()
    s = {
        "id": int(parts[1]),
        "lon": float(parts[3]),
        "lat": float(parts[4]),
    }
    return s


def pop_lines(lines: List[str], count: int) -> List[str]:
    """Pop a given number of lines from the start of a list.

    Parameters
    ----------
    lines: list[str]
        List of strings that is modified inplace
    count: int
        Number of items to remove

    Returns
    -------
    popped : list[str]
        Lines removed
    """
    popped = []
    for _ in range(count):
        popped.append(lines.pop(0))
    return popped


def pop_until(lines: List[str], pattern: str, include_match: bool = True) -> List[str]:
    """Pop lines from the start of a list until a regex pattern is found.

    Parameters
    ----------
    lines: list[str]
        List of strings that is modified inplace
    pattern: str
        Regex patten to match against, see `func`:re.match:.
    include_match: bool
        If the match should be removed.

    Returns
    -------
    popped : list[str]
        Lines removed
    """
    popped = []
    while not re.match(pattern, lines[0]):
        popped.append(lines.pop(0))

    if lines and include_match:
        popped.append(lines.pop(0))

    return popped


def read_out3(fname: str) -> Dict[str, T]:
    """Read a HAZ out3 formatted text file.

    Parameters
    ----------
    fname: str
        File name to open

    Returns
    -------
    site : dict
        Dictionary containing site information
    """
    assert fname.endswith(".out3")

    with open(fname) as fp:
        # Skip header line
        next(fp)
        lines = list(fp)

    parts = lines.pop(0).split()
    n_faults = int(parts[0])
    n_cases = int(parts[1])

    pop_until(lines, "^ SITE", include_match=False)

    site = parse_site_line(lines.pop(0))
    site["models"] = list()
    n_atten = int(lines.pop(0).split()[0])
    for i in range(n_atten):
        parts = lines.pop(0).split()
        m = dict(id=int(parts[1]), weight=float(parts[2]))
        n_ampl = int(lines.pop(0).split()[0])

        # Read the amplitudes
        m["ampl"] = fixed_split(lines.pop(0)[61:], n_ampl * [(12, float)])
        m["faults"] = []
        for _ in range(n_faults):
            values = fixed_split(
                lines.pop(0),
                [(40, str)]
                + 2 * [(6, float)]
                + [(9, float)]
                + n_ampl * [(12, float)],
            )
            f = {
                "name": values[0].strip(),
                "wt_segment": values[1],
                "wt_aleatory": values[2],
                "dist_min": values[3],
                "rate": values[4:],
            }
            m["faults"].append(f)

        names = ["rate", "prob", "mag_avg", "dist_avg", "epsilon_avg"]
        for name, line in zip(names, lines):
            m[name] = fixed_split(line[61:], n_ampl * [(12, float)])

        site["models"].append(m)

    return site


def test_read_out3():
    """Test reading of a HAZ out3 formatted text file."""
    site = read_out3(
        "../PEER_Verification_Tests/" + "Set2/S2Test1/Output/Set2Test1_Site1.out3"
    )

    np.testing.assert_allclose(site["lat"], 0.000)
    np.testing.assert_allclose(site["lon"], -65.000)

    model = site["models"][0]

    # Test the faults
    expected_faults = [
        (
            "01_AreaSource",
            1.000,
            1.000,
            5.0,
            [
                0.3944e-01,
                0.2283e-01,
                0.3926e-02,
                0.1337e-02,
                0.6207e-03,
                0.3293e-03,
                0.1888e-03,
                0.1142e-03,
                0.7184e-04,
                0.4663e-04,
                0.3105e-04,
                0.2114e-04,
                0.1466e-04,
                0.1034e-04,
                0.5371e-05,
                0.2928e-05,
                0.1662e-05,
                0.9771e-06,
            ],
        ),
        (
            "02_FaultB",
            1.000,
            1.000,
            50.0,
            [
                0.1221e-01,
                0.1018e-01,
                0.2854e-02,
                0.5060e-03,
                0.8081e-04,
                0.1412e-04,
                0.2798e-05,
                0.6268e-06,
                0.1569e-06,
                0.4335e-07,
                0.1308e-07,
                0.4272e-08,
                0.1498e-08,
                0.5611e-09,
                0.9376e-10,
                0.1926e-10,
                0.4726e-11,
                0.1350e-11,
            ],
        ),
        (
            "03_FaultC",
            1.000,
            1.000,
            25.0,
            [
                0.5912e-02,
                0.5858e-02,
                0.3982e-02,
                0.2153e-02,
                0.1052e-02,
                0.4863e-03,
                0.2224e-03,
                0.1030e-03,
                0.4876e-04,
                0.2371e-04,
                0.1185e-04,
                0.6082e-05,
                0.3204e-05,
                0.1730e-05,
                0.5397e-06,
                0.1827e-06,
                0.6642e-07,
                0.2574e-07,
            ],
        ),
    ]
    keys = ["name", "wt_segment", "wt_aleatory", "dist_min", "rate"]
    for i, ef in enumerate(expected_faults):
        for key, expected in zip(keys, ef):
            actual = model["faults"][i][key]
            if isinstance(expected, (str, int)):
                np.testing.assert_equal(actual, expected, err_msg="Key: %s" % key)
            else:
                np.testing.assert_allclose(actual, expected, err_msg="Key: %s" % key)

    # Test the total probabilities
    test_values = [
        (
            "ampl",
            [
                0.001000,
                0.010000,
                0.050000,
                0.100000,
                0.150000,
                0.200000,
                0.250000,
                0.300000,
                0.350000,
                0.400000,
                0.450000,
                0.500000,
                0.550000,
                0.600000,
                0.700000,
                0.800000,
                0.900000,
                1.000000,
            ],
        ),
        (
            "rate",
            [
                0.5755e-01,
                0.3887e-01,
                0.1076e-01,
                0.3996e-02,
                0.1754e-02,
                0.8297e-03,
                0.4140e-03,
                0.2178e-03,
                0.1207e-03,
                0.7038e-04,
                0.4291e-04,
                0.2723e-04,
                0.1787e-04,
                0.1207e-04,
                0.5911e-05,
                0.3111e-05,
                0.1729e-05,
                0.1003e-05,
            ],
        ),
        (
            "prob",
            [
                0.5593e-01,
                0.3812e-01,
                0.1071e-01,
                0.3988e-02,
                0.1752e-02,
                0.8293e-03,
                0.4139e-03,
                0.2178e-03,
                0.1207e-03,
                0.7038e-04,
                0.4291e-04,
                0.2723e-04,
                0.1787e-04,
                0.1207e-04,
                0.5911e-05,
                0.3111e-05,
                0.1729e-05,
                0.1003e-05,
            ],
        ),
        (
            "mag_avg",
            [
                0.553e01,
                0.565e01,
                0.602e01,
                0.611e01,
                0.613e01,
                0.611e01,
                0.608e01,
                0.603e01,
                0.598e01,
                0.593e01,
                0.588e01,
                0.584e01,
                0.581e01,
                0.578e01,
                0.573e01,
                0.570e01,
                0.568e01,
                0.566e01,
            ],
        ),
        (
            "dist_avg",
            [
                0.599e02,
                0.511e02,
                0.340e02,
                0.271e02,
                0.237e02,
                0.218e02,
                0.203e02,
                0.190e02,
                0.177e02,
                0.165e02,
                0.154e02,
                0.144e02,
                0.135e02,
                0.128e02,
                0.116e02,
                0.108e02,
                0.103e02,
                0.984e01,
            ],
        ),
        (
            "epsilon_avg",
            [
                -0.487e01,
                -0.171e01,
                -0.521e00,
                0.128e00,
                0.588e00,
                0.972e00,
                0.128e01,
                0.152e01,
                0.171e01,
                0.185e01,
                0.197e01,
                0.207e01,
                0.216e01,
                0.224e01,
                0.240e01,
                0.254e01,
                0.269e01,
                0.282e01,
            ],
        ),
    ]
    for key, expected in test_values:
        np.testing.assert_allclose(model[key], expected, err_msg="Key: %s" % key)


def read_out4(fname: str) -> Dict[str, T]:
    """Read a HAZ out4 formatted text file.

    Parameters
    ----------
    fname: str
        File name to open

    Returns
    -------
    site : dict
        Dictionary containing site information
    """
    assert fname.endswith(".out4")

    with open(fname) as fp:
        # Skip header line
        next(fp)
        lines = list(fp)

    pop_until(lines, "^ SITE", include_match=False)
    site = parse_site_line(lines.pop(0))
    n_prob = int(lines.pop(0).split()[0])
    n_ampl = int(lines.pop(0))

    site["amplitudes"] = fixed_split(lines.pop(0)[44:], n_ampl * [(12, float)])
    pop_lines(lines, 3)

    keys = ["mag_avg", "dist_avg", "epsilon_avg", "xcost_avg"]
    for key in keys:
        site[key] = fixed_split(lines.pop(0)[44:], n_ampl * [(12, float)])

    # Move two sections down
    pop_until(lines, "^------------------")
    pop_until(lines, "^------------------")
    pop_lines(lines, 4)

    site["bins"] = []
    for l in pop_until(lines, "^------------------")[:-2]:
        values = fixed_split(l[2:], 6 * [(7, float)] + n_ampl * [(12, float)])
        b = {
            "epsilon_range": (values[0], values[1]),
            "mag_range": (values[2], values[3]),
            "dist_range": (values[4], values[5]),
            "probabilities": values[6:],
        }
        site["bins"].append(b)

    return site


def test_read_out4():
    """Test reading of a HAZ out4 formatted text file."""
    site = read_out4(
        "../PEER_Verification_Tests/" + "Set2/S2Test1/Output/Set2Test1_Site1.out4"
    )

    np.testing.assert_allclose(site["lat"], 0.000)
    np.testing.assert_allclose(site["lon"], -65.000)

    test_averages = [
        (
            "amplitudes",
            [
                0.001,
                0.010,
                0.050,
                0.100,
                0.150,
                0.200,
                0.250,
                0.300,
                0.350,
                0.400,
                0.450,
                0.500,
                0.550,
                0.600,
                0.700,
                0.800,
                0.900,
                1.000,
            ],
        ),
        (
            "mag_avg",
            [
                0.553e01,
                0.565e01,
                0.602e01,
                0.611e01,
                0.613e01,
                0.611e01,
                0.608e01,
                0.603e01,
                0.598e01,
                0.593e01,
                0.588e01,
                0.584e01,
                0.581e01,
                0.578e01,
                0.573e01,
                0.570e01,
                0.568e01,
                0.566e01,
            ],
        ),
        (
            "dist_avg",
            [
                0.599e02,
                0.511e02,
                0.340e02,
                0.271e02,
                0.237e02,
                0.218e02,
                0.203e02,
                0.190e02,
                0.177e02,
                0.165e02,
                0.154e02,
                0.144e02,
                0.135e02,
                0.128e02,
                0.116e02,
                0.108e02,
                0.103e02,
                0.984e01,
            ],
        ),
        (
            "epsilon_avg",
            [
                -0.487e01,
                -0.171e01,
                -0.521e00,
                0.128e00,
                0.588e00,
                0.972e00,
                0.128e01,
                0.152e01,
                0.171e01,
                0.185e01,
                0.197e01,
                0.207e01,
                0.216e01,
                0.224e01,
                0.240e01,
                0.254e01,
                0.269e01,
                0.282e01,
            ],
        ),
        (
            "xcost_avg",
            [
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
            ],
        ),
    ]
    for key, actual in test_averages:
        if key in ["xcost_avg"]:
            # Don't test nan values
            continue
        np.testing.assert_allclose(actual, site[key])

    # Test the bins
    test_bins = [
        (
            (-12.8, 10.0),
            (5.0, 5.1),
            (0.0, 20.0),
            [
                0.459e-02,
                0.678e-02,
                0.181e-01,
                0.251e-01,
                0.297e-01,
                0.346e-01,
                0.401e-01,
                0.463e-01,
                0.525e-01,
                0.586e-01,
                0.642e-01,
                0.693e-01,
                0.738e-01,
                0.777e-01,
                0.843e-01,
                0.897e-01,
                0.943e-01,
                0.986e-01,
            ],
        ),
        (
            (-12.8, 10.0),
            (6.9, 7.0),
            (100.0, 200.0),
            [
                0.000e00,
                0.000e00,
                0.000e00,
                0.000e00,
                0.000e00,
                0.000e00,
                0.000e00,
                0.000e00,
                0.000e00,
                0.000e00,
                0.000e00,
                0.000e00,
                0.000e00,
                0.000e00,
                0.000e00,
                0.000e00,
                0.000e00,
                0.000e00,
            ],
        ),
    ]
    keys = ["epsilon_range", "mag_range", "dist_range", "probabilities"]
    for idx, test_bin in zip([0, -1], test_bins):
        for key, expected in zip(keys, test_bin):
            actual = site["bins"][idx][key]
            np.testing.assert_allclose(actual, expected)


def read_fractiles(fname):
    with open(fname) as fp:
        lines = list(fp)

    period = float(lines.pop(0)[:12])

    intensity = np.array(fixed_split(lines.pop(0), 18 * [(12, float)]))

    block = pop_until(lines, r"^\s+$", include_match=False)

    values = [fixed_split(b, [(7, float)] + 17 * [(12, float)]) for b in block]

    fractiles = np.array([v[0] for v in values])
    hazard = np.array([v[1:] for v in values])

    mean = np.array(fixed_split(lines.pop(0)[7:], 17 * [(12, float)]))
    pop_until(lines, r"^-+$", include_match=False)
    pop_lines(lines, 4)
    hazard = float(lines.pop(0)[35:])
    pop_lines(lines, 2)

    summary = np.array(
        [fixed_split(l, [(11, float)] + 5 * [(12, float)]) for l in lines],
        dtype=np.dtype(
            [
                ("period", "<f8"),
                ("5th", "<f8"),
                ("10th", "<f8"),
                ("50th", "<f8"),
                ("90th", "<f8"),
                ("95th", "<f8"),
            ]
        ),
    )

    # FIXME


def test_read_fractiles():
    fname = (
        "../PEER_Verification_Tests/"
        "Set3/S3Test2/Fractiles/Output/Fract_Set3Test2_Site1.out"
    )
    # FIXME
