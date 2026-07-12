#!/usr/bin/env python3

"""Fetch and safely extract the published Priority-4 legacy inputs."""

import argparse
import hashlib
from pathlib import Path
import shutil
import tempfile
from time import sleep
from urllib.request import Request, urlopen
import zipfile


ARCHIVE_URL = (
    "https://zenodo.org/api/records/4054618/files/"
    "raw_completions.zip/content"
)
ARCHIVE_SIZE = 198_689_883
ARCHIVE_MD5 = "f7a199f7718607e3740137e85c1d488b"
EXPECTED_FILES = 243
CHUNK_SIZE = 1024 * 1024


def md5_file(path):
    digest = hashlib.md5(usedforsecurity=False)
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(CHUNK_SIZE), b""):
            digest.update(block)
    return digest.hexdigest()


def archive_is_valid(path):
    path = Path(path)
    return (
        path.is_file()
        and path.stat().st_size == ARCHIVE_SIZE
        and md5_file(path) == ARCHIVE_MD5
    )


def download_archive(destination):
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if archive_is_valid(destination):
        print(f"Verified cached Zenodo archive: {destination}", flush=True)
        return destination

    temporary = destination.with_name(f".{destination.name}.part")
    request = Request(ARCHIVE_URL, headers={"User-Agent": "neutrinomass/1.0"})
    for attempt in range(1, 4):
        digest = hashlib.md5(usedforsecurity=False)
        size = 0
        try:
            print(f"Downloading {ARCHIVE_URL} (attempt {attempt}/3)", flush=True)
            with urlopen(request, timeout=120) as response, temporary.open(
                "wb"
            ) as out:
                while True:
                    block = response.read(CHUNK_SIZE)
                    if not block:
                        break
                    out.write(block)
                    digest.update(block)
                    size += len(block)
                    if size % (64 * CHUNK_SIZE) == 0:
                        print(f"Downloaded {size // CHUNK_SIZE} MiB", flush=True)
            if size != ARCHIVE_SIZE or digest.hexdigest() != ARCHIVE_MD5:
                raise ValueError(
                    "Downloaded Zenodo archive failed its published "
                    "size/checksum"
                )
            temporary.replace(destination)
            return destination
        except (OSError, ValueError):
            temporary.unlink(missing_ok=True)
            if attempt == 3:
                raise
            print("Download failed; retrying", flush=True)
            sleep(2**attempt)


def raw_members(archive):
    members = [
        item
        for item in archive.infolist()
        if not item.is_dir()
        and Path(item.filename).name.startswith("op_")
        and Path(item.filename).name.endswith(".dat")
    ]
    names = [Path(item.filename).name for item in members]
    if len(names) != EXPECTED_FILES or len(set(names)) != EXPECTED_FILES:
        raise ValueError(
            f"Expected {EXPECTED_FILES} unique op_*.dat files, found "
            f"{len(names)} files and {len(set(names))} unique names"
        )
    return members


def extract_archive(archive_path, destination):
    archive_path = Path(archive_path)
    destination = Path(destination)
    existing = sorted(destination.glob("op_*.dat")) if destination.exists() else []
    if len(existing) == EXPECTED_FILES:
        print(f"Found {EXPECTED_FILES} legacy inputs in {destination}", flush=True)
        return destination
    if existing:
        raise ValueError(
            f"Refusing to overwrite partial legacy directory {destination}: "
            f"found {len(existing)} op_*.dat files"
        )

    destination.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive_path) as archive:
        members = raw_members(archive)
        uncompressed_size = sum(item.file_size for item in members)
        free = shutil.disk_usage(destination.parent).free
        if free < uncompressed_size + 1024**3:
            raise ValueError(
                "Insufficient free space for the legacy archive: need at least "
                f"{(uncompressed_size + 1024**3) / 1024**3:.1f} GiB"
            )
        with tempfile.TemporaryDirectory(
            prefix=".raw-completions-extract-", dir=destination.parent
        ) as temporary:
            staged = Path(temporary) / destination.name
            staged.mkdir()
            for index, member in enumerate(members, 1):
                target = staged / Path(member.filename).name
                with archive.open(member) as source, target.open("wb") as output:
                    shutil.copyfileobj(source, output, length=CHUNK_SIZE)
                if index % 25 == 0 or index == len(members):
                    print(
                        f"Extracted {index}/{len(members)} legacy inputs",
                        flush=True,
                    )
            if destination.exists():
                destination.rmdir()
            staged.replace(destination)
    return destination


def ensure_legacy_archive(destination):
    destination = Path(destination).resolve()
    existing = sorted(destination.glob("op_*.dat")) if destination.exists() else []
    if len(existing) == EXPECTED_FILES:
        print(f"Found {EXPECTED_FILES} legacy inputs in {destination}", flush=True)
        return destination
    if existing:
        raise ValueError(
            f"Refusing to overwrite partial legacy directory {destination}: "
            f"found {len(existing)} op_*.dat files"
        )
    if destination.exists() and any(destination.iterdir()):
        raise ValueError(
            f"Refusing to overwrite non-empty legacy directory {destination}"
        )
    archive = destination.parent / "raw_completions.zip"
    download_archive(archive)
    return extract_archive(archive, destination)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    ensure_legacy_archive(args.destination)


if __name__ == "__main__":
    main()
