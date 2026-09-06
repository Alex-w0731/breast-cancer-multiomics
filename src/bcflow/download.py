"""Bounded, checksummed downloads and safe extraction for public research data."""

import hashlib
import tarfile
from pathlib import Path

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from .core import mkdir, sha256, write_json


def session():
    s = requests.Session()
    s.headers["User-Agent"] = "BCFlow/0.1 public-research-reproducibility"
    s.mount("https://", HTTPAdapter(max_retries=Retry(total=4, backoff_factor=1,
             status_forcelist=[429, 500, 502, 503, 504], allowed_methods=["GET", "HEAD"])))
    return s


def download(url, target, expected_md5=None, max_bytes=2_000_000_000):
    if not url.startswith("https://"):
        raise ValueError("HTTPS is required")
    target = Path(target)
    mkdir(target.parent)
    if target.exists():
        h = hashlib.md5(usedforsecurity=False)
        with target.open("rb") as f:
            for b in iter(lambda: f.read(1024 * 1024), b""):
                h.update(b)
        if expected_md5 and h.hexdigest() != expected_md5:
            raise ValueError(f"Existing {target} failed expected checksum; quarantine it first")
        if not expected_md5:
            raise FileExistsError(f"Refusing an unverified cached file: {target}")
        return target
    part = target.with_suffix(target.suffix + ".part")
    md5 = hashlib.md5(usedforsecurity=False)
    n = 0
    with session().get(url, stream=True, timeout=(15, 180)) as response:
        response.raise_for_status()
        with part.open("wb") as f:
            for chunk in response.iter_content(1024 * 1024):
                n += len(chunk)
                if n > max_bytes:
                    raise ValueError("Download exceeds the declared size bound")
                md5.update(chunk)
                f.write(chunk)
    if expected_md5 and md5.hexdigest() != expected_md5:
        raise ValueError("Checksum mismatch; .part retained for inspection")
    part.replace(target)
    write_json(target.with_suffix(target.suffix + ".provenance.json"), {
        "url": url, "bytes": n, "md5": md5.hexdigest(), "sha256": sha256(target),
        "expected_md5_verified": expected_md5 is not None})
    return target


def safe_extract(archive, destination):
    """Reject links, special files and all traversal before extracting any member."""
    destination = mkdir(destination).resolve()
    with tarfile.open(archive) as tar:
        members = tar.getmembers()
        for member in members:
            resolved = (destination / member.name).resolve()
            if not resolved.is_relative_to(destination) or not (member.isfile() or member.isdir()):
                raise ValueError(f"Unsafe archive member: {member.name}")
        tar.extractall(destination, members=members, filter="data")

