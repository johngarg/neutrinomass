#!/usr/bin/env python3

from neutrinomass.database.database import *


def test_conjugate_term():
    test_terms = [
        ["L.conj", "F,10,3,1/6,1", "F,10,3,7/6,1"],
        ["F,11,1,1/2,0", "F,11,2,0,0", "F,20,0,1/3,2", "F,20,1,5/6,2"],
        ["Q", "S,02,1,7/6,-2", "S,02,2,5/3,-2"],
        ["S,10,0,2/3,1", "S,11,0,1,0", "S,11,1,1/2,0"],
    ]
    conj_terms = [
        ["L", "F,01,3,-1/6,-1", "F,01,3,-7/6,-1"],
        ["F,11,1,-1/2,0", "F,11,2,0,0", "F,02,0,-1/3,-2", "F,02,1,-5/6,-2"],
        ["Q.conj", "S,20,1,-7/6,2", "S,20,2,-5/3,2"],
        ["S,01,0,-2/3,-1", "S,11,0,-1,0", "S,11,1,-1/2,0"],
    ]
    proc_terms = [conjugate_term(t) for t in test_terms]

    for i, sorted_conj in enumerate(proc_terms):
        assert list(sorted_conj) == sorted(conj_terms[i])
