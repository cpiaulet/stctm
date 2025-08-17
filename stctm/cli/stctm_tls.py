# stctm/cli/stctm_tls.py
from __future__ import annotations
import argparse
import sys
from pathlib import Path

def main(argv: list[str] | None = None) -> int:
    """
    Console entry point for:  stctm_TLS your_ini_file.ini
    """
    argv = sys.argv[1:] if argv is None else argv

    p = argparse.ArgumentParser(
        prog="stctm_TLS",
        description="Run TLS stellar contamination retrievals using an INI config."
    )
    p.add_argument("ini_path", help="Path to your INI file")
    args = p.parse_args(argv)

    ini = Path(args.ini_path)
    if not ini.is_file():
        p.error(f"INI file not found: {ini}")


    from stctm.runners.tls_runner import main as runner_main
    # The original runner expects argv with the ini path at index 1
    return int(runner_main(["stctm_TLS", str(ini)]) or 0)
