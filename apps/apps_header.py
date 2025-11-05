import sys

from pathlib import Path

HERE = Path(__file__).parent.resolve()
SRC = (HERE / ".." / "src").resolve()
sys.path.insert(0, str(SRC))