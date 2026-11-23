"""Exact interaction-multiset equivalence, independent of EFT provenance.

Operator names, topology, momentum expansions and EFT projections are not
part of this relation. Heavy particle nature, copies, charges, tensor
contractions and interaction multiplicities are. Hashes never decide equality.
"""

from collections import Counter, OrderedDict, defaultdict
from itertools import permutations, product

from neutrinomass.completions.completions import exotic_species_kind
from neutrinomass.completions.equivalence import (
    equivalent_lagrangians, field_label_parts,
)
from neutrinomass.database.serialization import dumps_completion, loads_completion


def canonical_species_descriptor(field):
    """Include zero-hypercharge conjugation and the full Dirac orbit."""
    states = [field, field.conj]
    if exotic_species_kind(field) == "dirac_fermion":
        states += [field.dirac_partner(), field.dirac_partner().conj]
    descriptions = [
        (tuple(f.sm_irrep), tuple(sorted((k, str(v)) for k, v in f.charges.items())))
        for f in states
    ]
    return exotic_species_kind(field), min(descriptions, key=repr)


def species_descriptors(completion):
    species = {}
    for field in completion.exotics:
        label = field_label_parts(field.label)[0]
        descriptor = canonical_species_descriptor(field)
        if label in species and species[label] != descriptor:
            raise ValueError("Inconsistent species identity in global comparison")
        species[label] = descriptor
    return species


def global_lagrangian_bucket(completion):
    """A necessary invariant only; contraction isomorphism remains required."""
    return (
        len(completion.terms),
        tuple(sorted(Counter(species_descriptors(completion).values()).items(), key=repr)),
    )


def equivalent_global_lagrangians(left, right):
    if global_lagrangian_bucket(left) != global_lagrangian_bucket(right):
        return False
    a, b = species_descriptors(left), species_descriptors(right)
    groups_a, groups_b = defaultdict(list), defaultdict(list)
    for label, descriptor in a.items():
        groups_a[descriptor].append(label)
    for label, descriptor in b.items():
        groups_b[descriptor].append(label)
    descriptors = sorted(groups_a, key=repr)
    sources = [sorted(groups_a[key]) for key in descriptors]
    targets = [permutations(sorted(groups_b[key])) for key in descriptors]
    for assignment in product(*targets):
        mapping = {x: y for xs, ys in zip(sources, assignment) for x, y in zip(xs, ys)}
        labels = sorted(mapping)
        orientations = [
            [(mapping[label], conjugate, partner)
             for conjugate in (False, True)
             for partner in ((False, True) if a[label][0] == "dirac_fermion" else (False,))]
            for label in labels
        ]
        for choices in product(*orientations):
            if equivalent_lagrangians(left.terms, right.terms, dict(zip(labels, choices))):
                return True
    return False


class DecodedCache:
    """Bound payload bytes and object count; cached objects are read-only."""

    def __init__(self, max_bytes=8 * 1024**2, max_records=32):
        self.max_bytes, self.max_records = max_bytes, max_records
        self.bytes = 0
        self.entries = OrderedDict()

    def get(self, payload):
        if payload in self.entries:
            self.entries.move_to_end(payload)
            return self.entries[payload]
        completion = loads_completion(payload)
        size = len(payload.encode("utf-8"))
        if size <= self.max_bytes and self.max_records:
            self.entries[payload] = completion
            self.bytes += size
            while self.bytes > self.max_bytes or len(self.entries) > self.max_records:
                old, _ = self.entries.popitem(last=False)
                self.bytes -= len(old.encode("utf-8"))
        return completion


class GlobalLagrangianIndex:
    """Disk-backed classes with exact comparison and unioned filter status."""

    def __init__(self, connection):
        self.connection = connection
        self.cache = DecodedCache()
        connection.executescript("""
            CREATE TABLE IF NOT EXISTS global_lagrangians (
                id INTEGER PRIMARY KEY, bucket TEXT NOT NULL,
                payload TEXT NOT NULL, filtered INTEGER NOT NULL);
            CREATE INDEX IF NOT EXISTS global_lagrangian_buckets
                ON global_lagrangians(bucket);
        """)

    def add(self, completion, *, filtered=False):
        bucket = repr(global_lagrangian_bucket(completion))
        for identifier, payload, was_filtered in self.connection.execute(
            "SELECT id, payload, filtered FROM global_lagrangians WHERE bucket=? ORDER BY id",
            (bucket,),
        ):
            if equivalent_global_lagrangians(completion, self.cache.get(payload)):
                if filtered and not was_filtered:
                    self.connection.execute(
                        "UPDATE global_lagrangians SET filtered=1 WHERE id=?", (identifier,))
                return identifier, False
        cursor = self.connection.execute(
            "INSERT INTO global_lagrangians(bucket,payload,filtered) VALUES (?,?,?)",
            (bucket, dumps_completion(completion), int(filtered)),
        )
        return cursor.lastrowid, True

    def counts(self):
        return self.connection.execute(
            "SELECT COUNT(*), COALESCE(SUM(filtered),0) FROM global_lagrangians"
        ).fetchone()
