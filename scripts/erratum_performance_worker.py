"""Use a certified global-matching backend with an unchanged pilot producer.

The original source still generates, validates and packages completions.
Only global interaction comparison uses the separately recorded backend.
Partial reduction snapshots keep their cursor and representatives; the original
manifest and immutable database are archived before an atomic contract change.
"""
import importlib.util
import json
import os
from pathlib import Path
import sys

from neutrinomass.database.checkpoints import (
    CheckpointStore, atomic_json, digest_file, identity_hash,
)
from neutrinomass.database.resumable_reduction import input_contract


def verified(reference):
    path = Path(reference["path"])
    if digest_file(path) != reference["sha256"]:
        raise ValueError("PERFORMANCE_REFERENCE_CHANGED")
    return path


def verify_backend(spec):
    from neutrinomass.completions import equivalence
    from neutrinomass.database import global_equivalence, resumable_reduction
    for name, module in (("equivalence", equivalence), ("global_equivalence", global_equivalence),
                         ("resumable_reduction", resumable_reduction)):
        if digest_file(module.__file__) != spec["baseline_files"][name]:
            raise ValueError("PERFORMANCE_BASELINE_CHANGED")
    backend = verified(spec["module"])
    certificate = json.loads(verified(spec["validation"]).read_text())
    if not (certificate.get("passed") is True and certificate["backend_sha256"] == digest_file(backend)
            and certificate["baseline_files"] == spec["baseline_files"]):
        raise ValueError("PERFORMANCE_BACKEND_NOT_VALIDATED")
    for reference in certificate["evidence"]:
        verified(reference)
    if not (certificate.get("full_local_suites") == [0, 1]
            and certificate.get("full_spartan_suites") == [0, 1]
            and certificate.get("full_ordered_71p_replays") == [0, 1]):
        raise ValueError("PERFORMANCE_VALIDATION_INCOMPLETE")
    return backend


def adopt_partial(directory, old_identity, new_identity, archive_root, spec):
    """Archive and rebind only a verified, compatible global reduction."""
    directory = Path(directory)
    manifest = directory/"checkpoint.json"
    if not manifest.exists():
        return
    state = json.loads(manifest.read_text())
    if state["contract"] == new_identity:
        with CheckpointStore(directory, new_identity):
            pass
        return
    if old_identity.get("relation") != "global-lagrangian":
        raise ValueError("PERFORMANCE_TRANSITION_NOT_GLOBAL")
    expected = dict(old_identity, performance_backend=spec)
    if expected != new_identity:
        raise ValueError("PERFORMANCE_TRANSITION_CHANGED_SCIENCE")
    with CheckpointStore(directory, old_identity) as original:
        if original.state["complete"]:
            raise ValueError("PERFORMANCE_TRANSITION_COMPLETED_STAGE")
        reference = dict(path=str(manifest.resolve()), sha256=digest_file(manifest))
        archive = Path(archive_root)/reference["sha256"]
        archive.mkdir(parents=True, exist_ok=True)
        for artifact in original.state["artifacts"].values():
            source = original.verify(artifact)
            target = archive/source.name
            if not target.exists():
                os.link(source, target)
            if digest_file(target) != artifact["sha256"]:
                raise ValueError("PERFORMANCE_ARCHIVE_CORRUPT")
        archived = archive/"checkpoint.json"
        if archived.exists() and digest_file(archived) != reference["sha256"]:
            raise ValueError("PERFORMANCE_ARCHIVE_MANIFEST_CHANGED")
        atomic_json(archived, original.state)
        if digest_file(archived) != reference["sha256"]:
            raise ValueError("PERFORMANCE_ARCHIVE_MANIFEST_CHANGED")
        updated = dict(original.state, contract=new_identity, contract_sha256=identity_hash(new_identity),
            performance_transition=dict(original_checkpoint=dict(path=str(archived),
                sha256=digest_file(archived)), original_producer=old_identity.get("source_commit"),
                validation=spec["validation"], preserved_cursor=True))
        atomic_json(manifest, updated)


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def install(config):
    from neutrinomass.database import global_equivalence, pipeline, resumable_reduction
    spec = config["performance_backend"]
    backend = load_module("certified_interaction_backend", verify_backend(spec))
    global_equivalence.equivalent_lagrangians = backend.equivalent_lagrangians
    original_reduce = resumable_reduction.reduce_resumable

    def reduce(paths, directory, contract, **kwargs):
        if kwargs.get("relation", "physical") != "global-lagrangian":
            return original_reduce(paths, directory, contract, **kwargs)
        paths = tuple(paths)
        manifest = Path(directory)/"checkpoint.json"
        old_identity = dict(contract, algorithm="exact-reduction-v1", relation="global-lagrangian",
                            inputs=input_contract(paths), ranking=None, maximum_rank=None)
        if kwargs.get("representative_rank") is not None or kwargs.get("maximum_rank") is not None:
            raise ValueError("PERFORMANCE_UNVALIDATED_RANKING")
        # An already completed original stage remains immutable with its actual
        # producer. Only unfinished/new stages acquire the backend contract.
        if manifest.exists():
            state = json.loads(manifest.read_text())
            if state["complete"] and state["contract"] == old_identity:
                return original_reduce(paths, directory, contract, **kwargs)
        new_contract = dict(contract, performance_backend=spec)
        new_identity = dict(old_identity, performance_backend=spec)
        adopt_partial(directory, old_identity, new_identity, config["performance_archive"], spec)
        return original_reduce(paths, directory, new_contract, **kwargs)

    pipeline.reduce_resumable = reduce
    resumable_reduction.reduce_resumable = reduce
    original_json = pipeline.atomic_json

    def publish(path, data):
        if Path(path).name in (".report.json", "report.json") and "stage_manifest" in data:
            data = dict(data, stage_manifest=dict(data["stage_manifest"],
                        performance_backends=dict(global_lagrangians=spec)))
        return original_json(path, data)

    pipeline.atomic_json = publish


def main():
    config = json.loads(Path(sys.argv[2]).read_text())
    install(config)
    worker = load_module("original_erratum_worker", verified(config["base_worker"]))
    original_provenance = worker.runtime_provenance

    def provenance():
        return dict(original_provenance(), performance_backend=config["performance_backend"],
                    performance_driver=dict(path=str(Path(__file__).resolve()), sha256=digest_file(__file__)))

    worker.runtime_provenance = provenance
    return worker.main()


if __name__ == "__main__":
    raise SystemExit(main())
