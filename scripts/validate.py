"""
Validate all of the OTTER JSON files
"""

import os
import glob
import json
import otter
import argparse
import logging
from pydantic import ValidationError

logger = logging.getLogger(__name__)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--otterdir", help="Directory where the otter json files will go")
    args = p.parse_args()

    validation_errors = ""
    for f in glob.glob(os.path.join(args.otterdir, "*.json")):
        with open(f, "r") as f2:
            j = json.load(f2)

        try:
            otter.schema.OtterSchema.model_validate(j)
        except ValidationError as e:
            validation_errors += f"""
{f}
-------------------------
{e}\n\n
"""
    logger.error(validation_errors)


if __name__ == "__main__":
    main()
