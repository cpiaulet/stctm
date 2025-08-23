# stctm/cli/stctm_exotune.py
import argparse
import sys
from pathlib import Path

def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    p = argparse.ArgumentParser(
        prog="stctm_exotune",
        description="Run ExoTune retrievals using an INI config."
    )
    p.add_argument("ini_path", help="Path to your INI file")
    # to allow -key=value overrides to pass through
    args, unknown = p.parse_known_args(argv)

    ini = Path(args.ini_path)
    if not ini.is_file():
        p.error(f"INI file not found: {ini}")

    from stctm.runners.exotune_runner import main as runner_main
    # forward any overrides so the inner parser in xtu can read them
    return int(runner_main(["stctm_exotune", str(ini), *unknown]) or 0)
