# stctm/cli/stctm_tls.py
import argparse
import sys
from pathlib import Path

def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv

    p = argparse.ArgumentParser(
        prog="stctm_TLS",
        description="Run TLS stellar contamination retrievals using an INI config."
    )
    p.add_argument("ini_path", help="Path to your INI file")

    # IMPORTANT: don't error on extra -key=value overrides
    args, unknown = p.parse_known_args(argv)

    ini = Path(args.ini_path)
    if not ini.is_file():
        p.error(f"INI file not found: {ini}")

    # Hand argv (including overrides) through to the runner
    from stctm.runners.tls_runner import main as runner_main
    return int(runner_main(["stctm_TLS", str(ini), *unknown]) or 0)
