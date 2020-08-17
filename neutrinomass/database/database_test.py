#!/usr/bin/env python3

from neutrinomass.database import *


def test_query():
    db_path = "/Users/johngargalionis/Desktop/operators/"
    mvdb = ModelDatabase(db_path)

    no_seesaws = (
        lambda m: not m.contains_field("F,00,2,0,0")
        and not m.contains_field("F,00,0,0,0")
        and not m.contains_field("S,00,2,1,0")
    )

    no_scalars = lambda m: not m.contains_field("S,*")

    for k, v in mvdb.query(lambda m: no_scalars(m) and no_seesaws(m)).items():
        assert not v
