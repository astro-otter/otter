"""
Drop some unconfirmed events that were historically proposed to be TDEs
but have no published classification
"""

import os
import otter


def main():
    OBJECTS_TO_DROP = [  # noqa
        "2020zeb",
        "ASASSN-20jr",
        "ASASSN-20il",
        "AT2021msu",
        "ASASSN-20lk",
    ]

    USER = "root"  # noqa
    ROOT_PASSWORD = os.environ.get("ARANGO_ROOT_PASSWORD")  # noqa
    URL = os.environ.get("ARANGO_URL")  # noqa

    if ROOT_PASSWORD is None:
        raise ValueError("Missing root password!")

    if URL is None:
        raise ValueError("Missing ARANGO_URL")

    db = otter.Otter(url=URL, username=USER, password=ROOT_PASSWORD)

    res = db.query(names=OBJECTS_TO_DROP)
    coll = db["transients"]
    for t in res:
        print(t)
        doc = coll.fetchDocument(t["_key"])
        doc.delete()


if __name__ == "__main__":
    main()
