from rassiparse.rassiparse import conf_diff

def test_conf_diff():
    c1 = "22222 00000"
    c2 = "2222u d0000"
    assert(set(conf_diff(c1, c2)) == set([(4, 5)]))

    c2 = "22220 du000"
    assert(set(conf_diff(c1, c2)) == set([(4, 5), (4, 6)]))

    # Ambiguous transitions
    c2 = "22200 dudu0"
    assert(set(conf_diff(c1, c2)) == set([]))

    # Spinflip
    c2 = "2222u u0000"
    assert(set(conf_diff(c1, c2)) == set([(4, 5)]))
