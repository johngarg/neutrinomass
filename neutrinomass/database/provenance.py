"""Source, import and pinned-environment provenance for stage contracts."""

import importlib.metadata
import os
from pathlib import Path
import platform
import subprocess
import sys
import hashlib

from neutrinomass.database.checkpoints import digest_file


def runtime_provenance(*,allow_dirty=False):
    import neutrinomass
    root=Path(__file__).resolve().parents[2]
    commit=subprocess.check_output(["git","rev-parse","HEAD"],cwd=root,text=True).strip()
    dirty=bool(subprocess.check_output(["git","status","--porcelain"],cwd=root,text=True).strip())
    if dirty and not allow_dirty:
        raise ValueError("A reproducible worker requires a clean detached checkout")
    dirty_digest=None
    if dirty:
        digest=hashlib.sha256(subprocess.check_output(["git","diff","HEAD","--binary"],cwd=root))
        untracked=subprocess.check_output(["git","ls-files","--others","--exclude-standard","-z"],cwd=root).split(b"\0")
        for name in sorted(filter(None,untracked)):
            digest.update(name)
            digest.update(digest_file(root/os.fsdecode(name)).encode())
        dirty_digest=digest.hexdigest()
    imported=Path(neutrinomass.__file__).resolve()
    if imported.parent != root/"neutrinomass":
        raise ValueError("Imported neutrinomass is outside the recorded source checkout")
    lock=root/"requirements-test.lock"
    pins={}
    for line in lock.read_text().splitlines():
        line=line.strip()
        if not line or line.startswith("#"):
            continue
        name,version=line.split("==")
        installed=importlib.metadata.version(name)
        if installed != version:
            raise ValueError(f"Dependency pin mismatch: {name}: {installed} != {version}")
        pins[name]=installed
    return dict(source_commit=commit,dirty_worktree=dirty,working_tree_digest=dirty_digest,lockfile_sha256=digest_file(lock),
        dependencies=pins,python_version=platform.python_version(),interpreter=sys.executable,
        imported_package=str(imported),source_directory=str(root),
        modules=os.environ.get("LOADEDMODULES",""),hash_seed=os.environ.get("PYTHONHASHSEED",""),
        platform=platform.platform(),slurm_job_id=os.environ.get("SLURM_JOB_ID"))


def scientific_contract(provenance):
    return {key:provenance[key] for key in (
        "source_commit","dirty_worktree","working_tree_digest","lockfile_sha256","dependencies","python_version","hash_seed")}
