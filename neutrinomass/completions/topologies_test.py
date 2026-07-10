#!/usr/bin/env python3

from pathlib import Path

from neutrinomass.completions.topologies import (
    canonical_topology_name,
    get_topology_data,
    paper_topology_name,
)


def test_published_5s2f_topology_names():
    assert paper_topology_name("5s2f_18") == "5s2f_10"
    assert paper_topology_name("5s2f_19") == "5s2f_11"
    assert paper_topology_name("5s2f_20") == "5s2f_13"
    assert paper_topology_name("5s2f_21") == "5s2f_14"

    assert canonical_topology_name("5s2f_18") == "5s2f_3"
    assert canonical_topology_name("5s2f_19") == "5s2f_4"
    assert canonical_topology_name("5s2f_20") == "5s2f_5"
    assert canonical_topology_name("5s2f_21") == "5s2f_6"

    for number in range(1, 25):
        canonical = f"5s2f_{number}"
        assert canonical_topology_name(paper_topology_name(canonical)) == canonical


def test_5s2f_topology_data_uses_published_order_and_matching_files():
    topology_data = get_topology_data(5, 2)

    assert [data["topology"] for data in topology_data] == [
        f"5s2f_{number}" for number in range(1, 25)
    ]
    assert [data["canonical_topology"] for data in topology_data] == [
        "5s2f_1",
        "5s2f_10",
        "5s2f_11",
        "5s2f_12",
        "5s2f_13",
        "5s2f_14",
        "5s2f_15",
        "5s2f_16",
        "5s2f_17",
        "5s2f_18",
        "5s2f_19",
        "5s2f_2",
        "5s2f_20",
        "5s2f_21",
        "5s2f_22",
        "5s2f_23",
        "5s2f_24",
        "5s2f_3",
        "5s2f_4",
        "5s2f_5",
        "5s2f_6",
        "5s2f_7",
        "5s2f_8",
        "5s2f_9",
    ]

    for data in topology_data:
        assert Path(data["partition_file"]).stem == data["canonical_topology"]
        assert Path(data["diagram_file"]).stem == data["canonical_topology"]
        assert Path(data["graph_file"]).stem == data["canonical_topology"]
