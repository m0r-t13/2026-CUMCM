"""ASCII entry point for the Windows launcher."""
from pathlib import Path
import runpy

if __name__ == "__main__":
    runpy.run_path(
        str(Path(__file__).with_name("\u68c0\u6d4b\u70b9\u7ed8\u56fe.py")),
        run_name="__main__",
    )
