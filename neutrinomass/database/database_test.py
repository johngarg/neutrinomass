#!/usr/bin/env python3

from neutrinomass.database.database import *
import pytest


class DummyLazyCompletion:
    def __init__(self, value):
        self.value = value
        self.head = {"terms": []}
        self.quantum_numbers = set()

    def force(self):
        return self.value


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


def test_force_preserves_every_operator():
    database = ModelDatabase(
        path=None,
        data={
            "operator-a": [DummyLazyCompletion("a1"), DummyLazyCompletion("a2")],
            "operator-b": [DummyLazyCompletion("b1")],
        },
    )

    database.force()

    assert database.data == {
        "operator-a": ["a1", "a2"],
        "operator-b": ["b1"],
    }
    assert database.is_forced


def test_path_database_starts_unordered(tmp_path):
    database = ModelDatabase(str(tmp_path))

    assert database.is_ordered is False

    database.filter()

    assert database.is_ordered is True


def test_legacy_executable_formats_require_explicit_trust(tmp_path):
    legacy = tmp_path / "op_test.dat"
    legacy.write_text("not_python()\n", encoding="utf-8")

    with pytest.raises(ValueError, match="execute Python"):
        read_completions(str(legacy))

    completion = LazyCompletion(
        head={"operator_name": "test", "quantum_numbers": []},
        tail="not_python()",
    )
    with pytest.raises(ValueError, match="execute Python"):
        completion.force()


def test_legacy_completions_can_be_streamed_after_explicit_trust(tmp_path):
    legacy = tmp_path / "op_test.dat"
    legacy.write_text(
        "LazyCompletion(head={'operator_name': 'test', "
        "'quantum_numbers': []}, tail='1')\n"
        "LazyCompletion(head={'operator_name': 'other', "
        "'quantum_numbers': []}, tail='2')\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="execute Python"):
        list(iter_completions(str(legacy)))

    streamed = list(iter_completions(str(legacy), trusted=True))

    assert [item.operator_name for item in streamed] == ["test", "other"]
    assert [item.force() for item in streamed] == [1, 2]
