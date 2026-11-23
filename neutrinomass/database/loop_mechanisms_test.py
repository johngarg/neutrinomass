from types import SimpleNamespace
from sympy import Rational

from neutrinomass.completions.core import ComplexScalar,MajoranaFermion
from neutrinomass.database.loop_mechanisms import candidate_loop_witness,concrete_vertices
from neutrinomass.database.serialization import loads_completion


def scotogenic(baryon=0):
    return SimpleNamespace(exotics={
        ComplexScalar("eta","i700",charges={"y":Rational(1,2),"3b":baryon}),
        MajoranaFermion("N","u700",charges={"y":0,"3b":0})})


def test_one_scalar_copy_can_occupy_both_triangle_edges():
    witness=candidate_loop_witness(scotogenic())
    assert witness["verified"]
    assert len(witness["required_species"]) == 2
    assert witness["numerator_kinds"].count("mass") == 1
    record=loads_completion(witness["record"])
    assert record.canonical_topology == "T3"
    assert all(term.safe_simplify() != 0 for term in record.terms)
    assert all(sum(field.charges.get("3b",0) for field in term.indexed_fields)==0
               for term in record.terms)
    assert set(witness["species_assignment"]) == {"eta","N"}


def test_loop_filter_uses_candidate_baryon_assignments():
    assert candidate_loop_witness(scotogenic(baryon=3)) is None


def test_antisymmetric_scalar_vertex_needs_two_admitted_copies():
    eta=ComplexScalar("eta","i701",charges={"y":Rational(-1,2),"3b":0})
    other=ComplexScalar("other","i702",charges=eta.charges)
    rho=ComplexScalar("rho","",charges={"y":1,"3b":0})
    assert not concrete_vertices((eta,eta,rho))
    assert concrete_vertices((eta,other,rho))
