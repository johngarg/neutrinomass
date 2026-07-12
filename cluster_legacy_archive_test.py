import hashlib
import zipfile

import pytest

from cluster import fetch_legacy_archive


def test_extract_archive_flattens_only_complete_raw_inventory(
    tmp_path, monkeypatch
):
    archive_path = tmp_path / "raw.zip"
    with zipfile.ZipFile(archive_path, "w") as archive:
        archive.writestr("raw_completions/op_1.dat", b"one")
        archive.writestr("nested/op_2.dat", b"two")
        archive.writestr("README.txt", b"ignored")
    monkeypatch.setattr(fetch_legacy_archive, "EXPECTED_FILES", 2)
    destination = tmp_path / "raw_completions"

    fetch_legacy_archive.extract_archive(archive_path, destination)

    assert sorted(path.name for path in destination.iterdir()) == [
        "op_1.dat",
        "op_2.dat",
    ]
    assert (destination / "op_1.dat").read_bytes() == b"one"


def test_download_archive_rejects_wrong_checksum(tmp_path, monkeypatch):
    payload = b"not-the-published-archive"
    source = tmp_path / "source.zip"
    source.write_bytes(payload)
    monkeypatch.setattr(fetch_legacy_archive, "ARCHIVE_URL", source.as_uri())
    monkeypatch.setattr(fetch_legacy_archive, "ARCHIVE_SIZE", len(payload))
    monkeypatch.setattr(fetch_legacy_archive, "ARCHIVE_MD5", "wrong")
    monkeypatch.setattr(fetch_legacy_archive, "sleep", lambda _: None)

    with pytest.raises(ValueError, match="published size/checksum"):
        fetch_legacy_archive.download_archive(tmp_path / "download.zip")

    assert not (tmp_path / "download.zip").exists()
    assert not (tmp_path / ".download.zip.part").exists()


def test_archive_checksum_constants_match_known_payload(tmp_path, monkeypatch):
    payload = b"known"
    archive = tmp_path / "raw_completions.zip"
    archive.write_bytes(payload)
    monkeypatch.setattr(fetch_legacy_archive, "ARCHIVE_SIZE", len(payload))
    monkeypatch.setattr(
        fetch_legacy_archive,
        "ARCHIVE_MD5",
        hashlib.md5(payload, usedforsecurity=False).hexdigest(),
    )

    assert fetch_legacy_archive.archive_is_valid(archive)
