import io
import tarfile

import pytest

from bcflow.download import safe_extract


def make_archive(path, member):
    with tarfile.open(path, "w") as tar:
        info = tarfile.TarInfo(member)
        info.size = 2
        tar.addfile(info, io.BytesIO(b"ok"))


def test_safe_archive_extracts(tmp_path):
    path = tmp_path / "safe.tar"
    make_archive(path, "nested/data.csv")
    safe_extract(path, tmp_path / "dest")
    assert (tmp_path / "dest/nested/data.csv").read_text() == "ok"


@pytest.mark.parametrize("member", ["../escape.txt", "/absolute.txt"])
def test_archive_traversal_is_rejected_before_extraction(tmp_path, member):
    path = tmp_path / "bad.tar"
    make_archive(path, member)
    with pytest.raises(ValueError, match="Unsafe archive"):
        safe_extract(path, tmp_path / "dest")
    assert not (tmp_path / "escape.txt").exists()


def test_archive_links_are_rejected(tmp_path):
    path = tmp_path / "links.tar"
    with tarfile.open(path, "w") as tar:
        info = tarfile.TarInfo("link")
        info.type = tarfile.SYMTYPE
        info.linkname = "../outside"
        tar.addfile(info)
    with pytest.raises(ValueError, match="Unsafe archive"):
        safe_extract(path, tmp_path / "dest")

