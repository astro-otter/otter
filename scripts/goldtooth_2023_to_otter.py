"""
Convert Table 1 from Goldtooth et al. (2023) to the OTTER format
and merge with the rest of the OTTER data (including adding Host info)
"""

import argparse
import re
import json
import os

import pandas as pd
import numpy as np

from astropy.coordinates import SkyCoord

from otter import Otter, Transient
from otter.util import bibcode_to_hrn


def read_goldtooth_catalog(filename):
    """
    Function to read the goldtooth TSV file for their
    catalog and output a pandas dataframe
    """

    # read in the TSV
    df = pd.read_csv(filename, sep="\t")
    del df["Unnamed: 7"]
    df = df.drop(0)

    # clean up the data
    df = df.replace({"cdots": np.nan})
    df["TDE Name"] = df["TDE Name"].str.replace(r"\^a", "")
    df["TDE Name"] = df["TDE Name"].str.replace(r"\^b", "")
    df["TDE Name"] = df["TDE Name"].str.replace(r"\^c", "")
    df["TDE Name"] = df["TDE Name"].str.strip()

    return df


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--otterdir", help="Directory where the otter json files will go")
    p.add_argument("--indir", help="Directory where dirty files are")
    args = p.parse_args()

    goldtooth_bibcode = "2023PASP..135c4101G"

    # get the reference mapping
    df = pd.read_csv(os.path.join(args.indir, "goldtooth_2023_citations_to_search.csv"))
    goldtooth_ref_map = {
        f'{row["name"]} {row.year}': row.bibcode for _, row in df.iterrows()
    }

    # connect to otter
    db = Otter(args.otterdir)

    # parse and clean the Table 1 dataset
    to_search = {"name": [], "year": [], "tde": []}

    data = read_goldtooth_catalog(os.path.join(args.indir, "goldtooth_2023_data.tsv"))

    for i, row in data.iterrows():
        tde = row["TDE Name"]
        host = row["Host Name"]

        try:
            disc_date = "-".join(row["Discovery Date"].split("(")[0].split("/"))
            disc_date_format = "iso"
        except Exception:
            disc_date = None

        # handle references
        bibcodes = [goldtooth_bibcode]
        refs = row.References
        vals = re.findall(r"(\w+(?: et al\.)? \(\d{4}[a-z]?(?:, \d{4}[a-z]?)*\))", refs)
        for val in vals:
            author, years = val.split("(")
            author = author.strip("et al.").strip()
            years = years.strip(")").strip()
            if "," in years:
                years = years.split(",")
                to_search["name"] += [author] * len(years)
                to_search["year"] += years
                to_search["tde"] += [tde] * len(years)

                for y in years:
                    k = f"{author} {y.strip()}"
                    bibcodes.append(goldtooth_ref_map[k])

            else:
                to_search["name"].append(author)
                to_search["year"].append(years)
                to_search["tde"].append(tde)
                k = f"{author} {years}"
                bibcodes.append(goldtooth_ref_map[k])

        # host info
        host_ra = row["Host R.A."]
        host_dec = row["Host Decl."]

        if not pd.isna(host_ra):
            host_coord = SkyCoord(host_ra, host_dec, unit=("hour", "deg"))
        else:
            host_coord = None

        if pd.isna(host):
            host = None

        # check otter for this object
        t = db.query(names=tde)
        if len(t) == 0:
            # this one didn't have a name matching the name given in Goldtooth et al.
            # let's try searching with the host coord
            t = db.cone_search(host_coord, radius=10)

        if len(t) > 0:
            t = t[0]
        else:
            t = None  # this is a new object for OTTER!

        # reformat the data into a Transient object
        tnew = Transient(
            {
                "schema_version": {
                    "value": 0,
                    "comments": "From Goldtooth et al. (2023)",
                },
                "name": {
                    "default_name": tde,
                    "alias": [{"value": tde, "reference": goldtooth_bibcode}],
                },
                "classification": [
                    dict(object_class="TDE", confidence=1, reference=bibcodes)
                ],
                "distance": [
                    dict(
                        value=row["z"],
                        reference=bibcodes,
                        computed=False,
                        default=True,
                        distance_type="redshift",
                    )
                ],
            }
        )

        if t is None:
            # just use the host coordinates, this should be good for TDEs
            tnew["coordinate"] = [
                dict(
                    ra=host_coord.ra.value,
                    dec=host_coord.dec.value,
                    ra_units="deg",
                    dec_units="deg",
                    coordinate_type="equitorial",
                    reference=bibcodes,
                    default=True,
                )
            ]
        else:
            # add coordinates from the existing objects
            tnew["coordinate"] = [t["coordinate"][0]]

        if disc_date is not None:
            tnew["date_reference"] = [
                dict(
                    value=disc_date,
                    date_format=disc_date_format,
                    date_type="discovery",
                    reference=bibcodes,
                    default=True,
                )
            ]

        if host_coord is not None:
            tnew["host"] = [
                dict(
                    host_name=host,
                    host_ra=host_coord.ra.value,
                    host_dec=host_coord.dec.value,
                    host_ra_units="deg",
                    host_dec_units="deg",
                    reference=bibcodes,
                )
            ]

        tnew["reference_alias"] = [
            {"name": b, "human_readable_name": bibcode_to_hrn(b)} for b in bibcodes
        ]

        if t is not None:
            # then we need to merge with existing data
            t_to_save = t + tnew

            # then overwrite existing file by finding the original file path
            tab = db.generate_summary_table()
            outfile = tab.json_path[tab["name"] == t.default_name].values[0]

        else:
            t_to_save = tnew

            # this doesn't exist yet!
            outfile = os.path.join(
                args.otterdir, t_to_save.default_name + ".json"
            ).replace(" ", "-")

        # print(json.dumps(dict(t_to_save), indent=4))
        # print()

        # print(tde, outfile)
        # print()

        with open(outfile, "w") as f:
            json.dump(dict(t_to_save), f, indent=4)


if __name__ == "__main__":
    main()
