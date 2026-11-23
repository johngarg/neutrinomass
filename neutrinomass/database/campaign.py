"""Conservative allocated-core accounting and bounded submission reservations."""

from contextlib import contextmanager
import fcntl
import json
from pathlib import Path
import time

from neutrinomass.database.checkpoints import atomic_json


class CampaignBudgetError(RuntimeError):
    pass


class CampaignLedger:
    def __init__(self,path):
        self.path=Path(path)

    @contextmanager
    def transaction(self):
        self.path.parent.mkdir(parents=True,exist_ok=True)
        with self.path.with_suffix(".lock").open("a+") as lock:
            fcntl.flock(lock,fcntl.LOCK_EX)
            state=json.loads(self.path.read_text()) if self.path.exists() else dict(
                schema_version=1,limit_core_seconds=3000*3600,
                pilot_limit_core_seconds=200*3600,max_workers=4,
                jobs={},production_gate=None)
            yield state
            atomic_json(self.path,state)

    @staticmethod
    def charged(state,phase=None):
        return sum(job.get("allocated_core_seconds",job["reserved_core_seconds"])
                   for job in state["jobs"].values()
                   if phase is None or job["phase"] == phase)

    def reserve(self,work_id,*,phase="validation",cpus=1,wall_seconds=7200,retry_of=None):
        if phase not in {"validation","production","finalisation"}:
            raise ValueError("Unknown campaign phase")
        if cpus != 1 or not 0 < wall_seconds <= 7200:
            raise ValueError("Workers require one CPU and at most two hours")
        with self.transaction() as state:
            if state.get("accounting_anomaly"):
                raise CampaignBudgetError("Unresolved Slurm allocation accounting anomaly")
            if work_id in state["jobs"]:
                raise ValueError("Stable submission identifier already reserved")
            if phase != "validation" and state["production_gate"] is None:
                raise CampaignBudgetError("71p production gate has not passed")
            live=[j for j in state["jobs"].values() if "allocated_core_seconds" not in j]
            if len(live) >= state["max_workers"]:
                raise CampaignBudgetError("Four-worker concurrency limit")
            amount=cpus*wall_seconds
            if self.charged(state)+amount > state["limit_core_seconds"]:
                raise CampaignBudgetError("Campaign allocation would exceed 3000 core-hours")
            if phase == "validation" and self.charged(state,"validation")+amount > state["pilot_limit_core_seconds"]:
                raise CampaignBudgetError("Initial validation allocation would exceed 200 core-hours")
            retries=0
            if retry_of is not None:
                old=state["jobs"][retry_of]
                if old.get("failure_class") != "infrastructure" or "allocated_core_seconds" not in old:
                    raise ValueError("Only diagnosed terminal infrastructure failures permit automatic retries")
                retries=old["retries"]+1
                if retries > 2:
                    raise ValueError("Two automatic infrastructure retries exhausted")
            state["jobs"][work_id]=dict(phase=phase,cpus=cpus,wall_seconds=wall_seconds,
                reserved_core_seconds=amount,status="reserved",created_unix=time.time(),
                retries=retries,retry_of=retry_of)
            return dict(state["jobs"][work_id])

    def submitted(self,work_id,job_id):
        with self.transaction() as state:
            job=state["jobs"][work_id]
            if job.get("job_id") not in {None,str(job_id)}:
                raise ValueError("Submission already bound to another Slurm job")
            job.update(job_id=str(job_id),status="submitted")

    def finished(self,work_id,*,elapsed_seconds,allocated_cpus,status,failure_class=None):
        if elapsed_seconds < 0 or allocated_cpus < 1:
            raise ValueError("Invalid Slurm allocation accounting")
        if failure_class not in {None,"infrastructure","scientific","corruption","no-progress"}:
            raise ValueError("Unknown failure diagnosis")
        with self.transaction() as state:
            job=state["jobs"][work_id]
            actual=elapsed_seconds*allocated_cpus
            if "allocated_core_seconds" in job and job["allocated_core_seconds"] != actual:
                raise ValueError("Previously recorded accounting changed")
            job.update(status=status,allocated_core_seconds=actual,failure_class=failure_class,
                       finished_unix=time.time())
            if actual > job["reserved_core_seconds"]:
                state["accounting_anomaly"]=dict(work_id=work_id,actual=actual)

    def summary(self):
        with self.transaction() as state:
            charged=self.charged(state)
            return dict(state,charged_core_hours=charged/3600,
                        remaining_core_hours=(state["limit_core_seconds"]-charged)/3600,
                        validation_charged_core_hours=self.charged(state,"validation")/3600)
