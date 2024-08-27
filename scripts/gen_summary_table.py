"""
Script to generate the summary table for OTTER searches

This needs to be run everytime this repo is downloaded right now because
the summary table stores absolute paths in it which is user dependent
"""

import os
import argparse
from otter import Otter


def main():
    p = argparse.ArgumentParser()
    p.add_argument(
        "--otterroot", help="Root directory of the otter installation", required=False
    )
    p.add_argument(
        "--datadir", help="Directory where the otter json are", required=False
    )
    args = p.parse_args()

    if args.otterroot is not None and args.datadir is None:
        otterpath = os.path.join(args.otterroot, "otterdb", ".otter")
    elif args.otterroot is None and args.datadir is not None:
        otterpath = args.datadir
    else:
        raise ValueError("Incorrect arguments provided!")

    # find the otter directory and then generate the summary table
    db = Otter(datadir=otterpath)
    db.generate_summary_table(save=True)


if __name__ == "__main__":
    main()
