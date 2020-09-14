#!/usr/bin/env python3

from neutrinomass.database.database import *

# db_path = "/Users/johngargalionis/Desktop/operators/"
# mvdb = ModelDatabase(db_path)


# def test_query():
#     no_scalars = lambda m: not m.contains_field("S,*")

#     func = lambda m: no_scalars(m) and ModelDatabase.no_seesaws(m)
#     for k, v in mvdb.query(func).items():
#         assert not v


# def test_process():
#     mvdb.process(filter_seesaws=True)
#     assert len(mvdb.data["8p"]) == 8


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


def test_query():
    import pickle
    import pandas as pd

    data = pickle.load(
        open("/Users/johngargalionis/Desktop/raw_democratic_data.p", "rb")
    )
    proc_data = {
        k: [LazyCompletion(head=x, tail=y) for x, y in v] for k, v in data.items()
    }
    mvdb = ModelDatabase(path="", data=proc_data)
    return mvdb
