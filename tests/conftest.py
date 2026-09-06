import os
from pathlib import Path

cache = Path("work/test_matplotlib")
cache.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(cache.resolve()))
