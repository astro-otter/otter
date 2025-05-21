import os
import glob
import json
import otter
import argparse


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--otterdir", help="Directory where the otter json files will go")
    args = p.parse_args()

    for f in glob.glob(os.path.join(args.otterdir, "*.json")):
        with open(f, "r") as f2:
            j = json.load(f2)

        t = otter.Transient(j)
        t = t._derive_classification_flags(t)

        with open(f, "w") as f2:
            json.dump(dict(t), f2, indent=4)


if __name__ == "__main__":
    main()
