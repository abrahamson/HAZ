import re
import numpy as np


def fixed_split(line, pairs, keep_tail=False):
    """Read a fixed width line of text.

    Args:
        line (str): text line to parse
        pairs (list): list of tuples specifying the width and the string parser
        keep_tail (bool): remove entries without values

    Returns:
        list: values read
    """
    line = line.rstrip()

    i = 0
    values = []
    for w, p in pairs:
        try:
            v = p(line[i: i + w])
        except ValueError:
            v = None

        values.append(v)
        i += w

    if not keep_tail:
        while values[-1] is None:
            values.pop()

    return values


def parse_blocks(lines, pattern):
    """Separate lines into a series of blocks using an identifying pattern.

    Parameters
    ----------
    lines: list[str]
        List of strings
    pattern: str
        Regex pattern used to identify the blocks.
    Returns
    -------
    list[str]
    """
    starts = [i for i, l in enumerate(lines)
              if re.search(pattern, l)]
    starts.append(len(lines))
    blocks = [lines[s: starts[i + 1]] for i, s in enumerate(starts[:-1])]
    return blocks


def parse_site_line(line):
    """Parse the site line in out3 and out4 files.

    Parameters
    ----------
    line: str
        Line of text
    Returns
    -------
    dict
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
        'id': int(parts[1]),
        'lon': np.float(parts[3]),
        'lat': np.float(parts[4]),
    }
    return s


def pop_lines(lines, count):
    """Pop a given number of lines from the start of a list.

    Parameters
    ----------
    lines: list[str]
        List of strings that is modified inplace
    count: int
        Number of items to remove

    Returns
    -------
    list[str]
        Lines removed
    """
    popped = []
    for _ in range(count):
        popped.append(lines.pop(0))
    return popped


def pop_until(lines, pattern, include_match=True):
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
    list[str]
        Lines removed
    """
    popped = []
    while not re.match(pattern, lines[0]):
        popped.append(lines.pop(0))

    if lines and include_match:
        popped.append(lines.pop(0))

    return popped


def read_out3(fname):
    """Read a HAZ out3 formatted text file.

    Parameters
    ----------
    fname: str
        File name to open

    Returns
    -------
    dict
    """
    assert fname.endswith('.out3')

    with open(fname) as fp:
        lines = list(fp)

    parts = lines.pop(0).split()
    n_faults = int(parts[0])
    n_cases = int(parts[1])

    pop_until(lines, '^ SITE', include_match=False)

    site = parse_site_line(lines.pop(0))
    site['models'] = list()
    n_atten = int(lines.pop(0).split()[0])
    for i in range(n_atten):
        parts = lines.pop(0).split()
        m = dict(
            id=int(parts[1]),
            weight=float(parts[2])
        )
        n_ampl = int(lines.pop(0).split()[0])

        # Read the amplitudes
        m['ampl'] = fixed_split(lines.pop(0)[61:],
                                n_ampl * [(12, np.float)])
        m['faults'] = []
        for _ in range(n_faults):
            values = fixed_split(
                lines.pop(0), [(40, str)] + 2 * [(6, np.float)] +
                              [(9, np.float)] + n_ampl * [(12, np.float)])
            f = {
                'name': values[0].strip(),
                'wt_segment': values[1],
                'wt_aleatory': values[2],
                'dist_min': values[3],
                'rate': values[4:]
            }
            m['faults'].append(f)

        names = ['rate', 'prob', 'mag_avg', 'dist_avg', 'epsilon_avg']
        for name, line in zip(names, lines):
            m[name] = fixed_split(line[61:], n_ampl * [(12, np.float)])

        site['models'].append(m)

    return site


def test_read_out3():
    """Test reading of a HAZ out3 formatted text file."""
    site = read_out3('../PEER Verification Tests/' +
                     'Set2/S2Test1/Output/Set2Test1_Site1.out3')

    np.testing.assert_allclose(site['lat'], 0.000)
    np.testing.assert_allclose(site['lon'], -65.000)

    model = site['models'][0]

    # Test the faults
    expected_faults = [
        ('01_AreaSource', 1.000, 1.000, 0.7, [
            0.3944E-01, 0.2283E-01, 0.3926E-02, 0.1337E-02, 0.6207E-03,
            0.3293E-03, 0.1888E-03, 0.1142E-03, 0.7184E-04, 0.4663E-04,
            0.3105E-04, 0.2114E-04, 0.1466E-04, 0.1034E-04, 0.5371E-05,
            0.2928E-05, 0.1662E-05, 0.9771E-06]),
        ('01_FaultB', 1.000, 1.000, 50.0, [
            0.1221E-01, 0.1018E-01, 0.2854E-02, 0.5060E-03, 0.8081E-04,
            0.1412E-04, 0.2798E-05, 0.6268E-06, 0.1569E-06, 0.4335E-07,
            0.1308E-07, 0.4272E-08, 0.1498E-08, 0.5611E-09, 0.9376E-10,
            0.1926E-10, 0.4726E-11, 0.1350E-11]),
        ('01_FaultC', 1.000, 1.000, 25.0, [
            0.5912E-02, 0.5858E-02, 0.3982E-02, 0.2153E-02, 0.1052E-02,
            0.4863E-03, 0.2224E-03, 0.1030E-03, 0.4876E-04, 0.2371E-04,
            0.1185E-04, 0.6082E-05, 0.3204E-05, 0.1730E-05, 0.5397E-06,
            0.1827E-06, 0.6642E-07, 0.2574E-07]),
    ]
    keys = ['name', 'wt_segment', 'wt_aleatory', 'dist_min', 'rate']
    for i, ef in enumerate(expected_faults):
        for key, expected in zip(keys, ef):
            actual = model['faults'][i][key]
            if isinstance(expected, (str, int)):
                np.testing.assert_equal(actual, expected)
            else:
                np.testing.assert_allclose(actual, expected)

    # Test the total probabilities
    test_values = [
        ('ampl', [0.001000, 0.010000, 0.050000, 0.100000, 0.150000, 0.200000,
                  0.250000, 0.300000, 0.350000, 0.400000, 0.450000, 0.500000,
                  0.550000, 0.600000, 0.700000, 0.800000, 0.900000, 1.000000]
         ),
        ('rate', [0.5755E-01, 0.3887E-01, 0.1076E-01, 0.3996E-02, 0.1754E-02,
                  0.8297E-03, 0.4140E-03, 0.2178E-03, 0.1207E-03, 0.7038E-04,
                  0.4291E-04, 0.2723E-04, 0.1787E-04, 0.1207E-04, 0.5911E-05,
                  0.3111E-05, 0.1729E-05, 0.1003E-05]),
        ('prob', [0.5593E-01, 0.3812E-01, 0.1071E-01, 0.3988E-02, 0.1752E-02,
                  0.8293E-03, 0.4139E-03, 0.2178E-03, 0.1207E-03, 0.7038E-04,
                  0.4291E-04, 0.2723E-04, 0.1787E-04, 0.1207E-04, 0.5911E-05,
                  0.3111E-05, 0.1729E-05, 0.1003E-05]),
        ('mag_avg',
         [0.553E+01, 0.565E+01, 0.602E+01, 0.611E+01, 0.613E+01, 0.611E+01,
          0.608E+01, 0.603E+01, 0.598E+01, 0.593E+01, 0.588E+01, 0.584E+01,
          0.581E+01, 0.578E+01, 0.573E+01, 0.570E+01, 0.568E+01, 0.566E+01]),
        ('dist_avg',
         [0.599E+02, 0.511E+02, 0.340E+02, 0.271E+02, 0.237E+02, 0.218E+02,
          0.203E+02, 0.190E+02, 0.177E+02, 0.165E+02, 0.154E+02, 0.144E+02,
          0.135E+02, 0.128E+02, 0.116E+02, 0.108E+02, 0.103E+02, 0.984E+01]),
        ('epsilon_avg',
         [-0.487E+01, -0.171E+01, -0.521E+00, 0.128E+00, 0.588E+00, 0.972E+00,
          0.128E+01, 0.152E+01, 0.171E+01, 0.185E+01, 0.197E+01, 0.207E+01,
          0.216E+01, 0.224E+01, 0.240E+01, 0.254E+01, 0.269E+01, 0.282E+01]),
    ]
    for key, expected in test_values:
        np.testing.assert_allclose(model[key], expected)


def read_out4(fname):
    """Read a HAZ out4 formatted text file.

    Parameters
    ----------
    fname: str
        File name to open

    Returns
    -------
    dict
    """
    assert fname.endswith('.out4')

    with open(fname) as fp:
        lines = list(fp)

    pop_until(lines, '^ SITE', include_match=False)
    site = parse_site_line(lines.pop(0))
    n_prob = int(lines.pop(0).split()[0])
    n_ampl = int(lines.pop(0))

    site['amplitudes'] = fixed_split(lines.pop(0)[44:],
                                     n_ampl * [(12, np.float)])
    pop_lines(lines, 3)

    keys = ['mag_avg', 'dist_avg', 'epsilon_avg', 'xcost_avg']
    for key in keys:
        site[key] = fixed_split(lines.pop(0)[44:], n_ampl * [(12, np.float)])

    # Move two sections down
    pop_until(lines, '^------------------')
    pop_until(lines, '^------------------')
    pop_lines(lines, 4)

    site['bins'] = []
    for l in pop_until(lines, '^------------------')[:-2]:
        values = fixed_split(
            l[2:], 6 * [(7, np.float)] + n_ampl * [(12, np.float)])
        b = {
            'epsilon_range': (values[0], values[1]),
            'mag_range': (values[2], values[3]),
            'dist_range': (values[4], values[5]),
            'probabilities': values[6:]
        }
        site['bins'].append(b)

    return site


def test_read_out4():
    """Test reading of a HAZ out4 formatted text file."""
    site = read_out4('../PEER Verification Tests/' +
                     'Set2/S2Test1/Output/Set2Test1_Site1.out4')

    np.testing.assert_allclose(site['lat'], 0.000)
    np.testing.assert_allclose(site['lon'], -65.000)

    test_averages = [
        ('amplitudes',
         [0.001, 0.010, 0.050, 0.100, 0.150, 0.200, 0.250, 0.300, 0.350, 0.400,
          0.450, 0.500, 0.550, 0.600, 0.700, 0.800, 0.900, 1.000]),
        ('mag_avg',
         [0.553E+01, 0.565E+01, 0.602E+01, 0.611E+01, 0.613E+01, 0.611E+01,
          0.608E+01, 0.603E+01, 0.598E+01, 0.593E+01, 0.588E+01, 0.584E+01,
          0.581E+01, 0.578E+01, 0.573E+01, 0.570E+01, 0.568E+01, 0.566E+01]),
        ('dist_avg',
         [0.599E+02, 0.511E+02, 0.340E+02, 0.271E+02, 0.237E+02, 0.218E+02,
          0.203E+02, 0.190E+02, 0.177E+02, 0.165E+02, 0.154E+02, 0.144E+02,
          0.135E+02, 0.128E+02, 0.116E+02, 0.108E+02, 0.103E+02, 0.984E+01]),
        ('epsilon_avg',
         [-0.487E+01, -0.171E+01, -0.521E+00, 0.128E+00, 0.588E+00, 0.972E+00,
          0.128E+01, 0.152E+01, 0.171E+01, 0.185E+01, 0.197E+01, 0.207E+01,
          0.216E+01, 0.224E+01, 0.240E+01, 0.254E+01, 0.269E+01, 0.282E+01]),
        ('xcost_avg',
         [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
          np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
          np.nan, np.nan]),
    ]
    for key, actual in test_averages:
        if key in ['xcost_avg']:
            # Don't test nan values
            continue
        np.testing.assert_allclose(actual, site[key])

    # Test the bins
    test_bins = [
        ((-12.8, 10.0), (5.0, 5.1), (0.0, 20.0),
         [0.459E-02, 0.678E-02, 0.181E-01, 0.251E-01, 0.297E-01, 0.346E-01,
          0.401E-01, 0.463E-01, 0.525E-01, 0.586E-01, 0.642E-01, 0.693E-01,
          0.738E-01, 0.777E-01, 0.843E-01, 0.897E-01, 0.943E-01, 0.986E-01]),
        ((-12.8, 10.0), (6.9, 7.0), (100.0, 200.0),
         [0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
          0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
          0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00])
    ]
    keys = ['epsilon_range', 'mag_range', 'dist_range', 'probabilities']
    for idx, test_bin in zip([0, -1], test_bins):
        for key, expected in zip(keys, test_bin):
            actual = site['bins'][idx][key]
            np.testing.assert_allclose(actual, expected)