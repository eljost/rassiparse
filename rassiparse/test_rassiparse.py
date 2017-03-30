import pytest

from rassiparse.rassiparse import parse_rassi, parse_subspaces


def test_bodypi_rassi():
    # 4 JOB.Mix
    # C2v symmetry
    # RAS(24,2,2;9,6,9)

    fn = "../logs/bodypi_rassi.out"
    irreps = 4
    with open(fn) as handle:
        text = handle.read()

    sf_states, trans_dict = parse_rassi(text)

    assert(len(sf_states) == 8)
    assert(trans_dict[(1,3)] == pytest.approx(1.0902241))
    assert(trans_dict[(1,2)] == pytest.approx(0.39636729E-03))

    irrep_inds = parse_subspaces(text, irreps)
    assert(irrep_inds == [[3], [0, 4, 7], [1, 5, 8], [2, 6, 9]])
