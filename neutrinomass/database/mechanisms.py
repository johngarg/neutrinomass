"""Candidate-specific higher-scale tree searches with concrete UV witnesses.

The search domain is the published scalar/fermion operator inventory. Absence
of an additional tree in this domain is not a proof of low-scale genuineness.
"""

from contextlib import closing

from neutrinomass.completions.species import candidate_completion_stream
from neutrinomass.completions.amplitudes import audit_amplitude_symmetrisation
from neutrinomass.database.global_equivalence import species_descriptors
from neutrinomass.database.serialization import dumps_completion,loads_completion
from neutrinomass.database.rebuild import operator_registry,operator_scale_tev,validate_completion


def candidate_tree_witness(candidate,higher_operators,*,topology=None,progress=None):
    """Exhaust every admitted species assignment, stopping at a proved witness.

    Copies are assigned before each vertex's symmetry check. Every propagator
    comes from candidate.exotics with its original baryon number. Couplings
    are generic; the candidate's listed interaction subset is not a symmetry.
    """
    operators=operator_registry()
    available=species_descriptors(candidate)
    for name in higher_operators:
        operator=operators[name]
        generated=0
        with closing(candidate_completion_stream(operator,candidate.exotics,topology=topology)) as stream:
            for completion in stream:
                generated+=1
                decoded=loads_completion(dumps_completion(completion))
                validate_completion(decoded,check_vanishing=False)
                if any(term.safe_simplify() == 0 for term in decoded.terms):
                    continue
                audit=audit_amplitude_symmetrisation(decoded)
                if audit.status == "unsupported":
                    raise ValueError(f"UNSUPPORTED_CANDIDATE_MECHANISM: {audit.reason}")
                if audit.is_zero:
                    continue
                required=species_descriptors(decoded)
                if any(available.get(label) != descriptor for label,descriptor in required.items()):
                    raise ValueError("CANDIDATE_MECHANISM_CHANGED_SPECIES")
                return dict(code="HIGHER_SCALE_CANDIDATE_TREE",verified=True,
                    operator=name,scale_tev=operator_scale_tev(operator),
                    required_species=tuple(sorted(required.values(),key=repr)),
                    species_assignment=dict(sorted(required.items())),
                    record=dumps_completion(decoded),amplitude_witness=audit.witness)
        if progress is not None:
            progress(dict(operator=name,complete=True,generated=generated))
    return None
