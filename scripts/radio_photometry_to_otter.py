"""
Convert the Radio Photometry that Noah, Collin, and Kate gathered to the Otter format
"""

# import
import os
import argparse
import otter
import pandas as pd
import numpy as np


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--otterdir", help="Directory where the otter json files will go")
    p.add_argument("--indir", help="Directory where dirty files are")
    args = p.parse_args()

    db = otter.Otter(args.otterdir)

    # read in the metadata and photometry files
    meta = pd.read_csv(os.path.join(args.indir, "tde_radio_data_meta.csv"))
    phot = pd.read_csv(os.path.join(args.indir, "tde_radio_data_photometry.csv"))

    # drop duplicated names in meta and keep the first
    meta = meta.drop_duplicates(subset="name", keep="first")

    # merge the meta and phot data
    data = pd.merge(phot, meta, on="name", how="inner")

    # perform some data checks
    assert (
        len(data[pd.isna(data.ra)].name.unique()) == 0
    ), "Missing some RA and Decs, please check the input files!"
    for name in meta.name:
        assert len(data[data.name == name]) == len(
            phot[phot.name == name]
        ), f"failed on {name}"

    # actually do the data conversion to OTTER
    all_jsons = []

    for name, tde in data.groupby("name"):
        ########### FOR NOW JUST CAUSE IT'S BROKEN ###################
        if name in {"AT2020opy", "AT2019azh", "AT2020vwl", "CNSS J0019+00"}:
            continue

        json = {}
        tde = tde.reset_index()

        # name first
        json["name"] = dict(
            default_name=name,
            alias=[dict(value=name, reference=[tde.coord_bibcode[0]])],
        )

        # coordinates
        json["coordinate"] = [
            dict(
                ra=tde.ra[0],
                dec=tde.dec[0],
                ra_units=tde.ra_unit[0],
                dec_units=tde.dec_unit[0],
                reference=[tde.coord_bibcode[0]],
                coordinate_type="equitorial",
            )
        ]

        # redshift
        if not np.any(pd.isna(tde["redshift"])):
            json["distance"] = [
                dict(
                    value=tde.redshift[0],
                    reference=[tde.redshift_ref[0]],
                    compute=False,
                    distance_type="redshift",
                )
            ]

        json["classification"] = [
            dict(
                object_class="TDE",
                confidence=1,  # we know this is at least an tde
                reference=[tde.bibcode[0]],
            )
        ]

        # discovery date
        # print(tde)
        if not np.any(pd.isna(tde.discovery_date)):
            json["date_reference"] = [
                dict(
                    value=tde.discovery_date.tolist()[0].strip(),
                    date_format=tde.discovery_date_format.tolist()[0].lower(),
                    reference=tde.discovery_date_ref.tolist()[0],
                    computed=False,
                    date_type="discovery",
                )
            ]

        # now the radio photometry
        phot_sources = []
        json["photometry"] = []
        for (src, tele), phot in tde.groupby(["bibcode", "telescope"]):
            if src not in phot_sources:
                phot_sources.append(src)

            if len(np.unique(phot.flux_unit)) == 1:
                raw_units = phot.flux_unit.tolist()[0]
            else:
                raw_units = phot.flux_unit.values

            json_phot = dict(
                reference=src,
                raw=phot.flux.tolist(),
                raw_err=phot.flux_err.tolist(),
                raw_units=raw_units,
                date=phot.date.tolist(),
                date_format=phot.date_format.tolist(),
                upperlimit=phot.upperlimit.tolist(),
                filter_key=phot["band_name"].tolist(),
                obs_type="radio",
                telescope=tele,
            )

            corrs = ["corr_k", "corr_s", "corr_host", "corr_av", "corr_hostav"]
            for c in corrs:
                json_phot[c] = False if np.all(pd.isna(phot[c])) else phot[c].tolist()
                if np.any(json_phot[c]):
                    v = c.replace("corr", "val")
                    json_phot[v] = phot[v].tolist()

            json["photometry"].append(json_phot)

        # filter alias
        # radio filters first
        filter_keys1 = ["band_name", "band_eff_freq", "band_eff_freq_unit"]
        filter_map = (
            tde[filter_keys1].drop_duplicates().set_index("band_name")
        )  # .to_dict(orient='index')
        try:
            filter_map_radio = filter_map.to_dict(orient="index")
        except Exception:
            print(filter_map)
            print(name)
            raise Exception
        json["filter_alias"] = [
            dict(
                filter_key=filt,
                freq_eff=val["band_eff_freq"],
                freq_units=val["band_eff_freq_unit"],
            )
            for filt, val in filter_map_radio.items()
        ]

        # reference alias
        json["reference_alias"] = [
            dict(
                name=tde.coord_bibcode[0],
                human_readable_name=otter.util.bibcode_to_hrn(tde.coord_bibcode[0]),
            )
        ] + [dict(name=val, human_readable_name=val) for val in phot_sources]

        curr_refs = {v["name"] for v in json["reference_alias"]}
        if tde.redshift_ref[0] not in curr_refs and not np.any(pd.isna(tde.redshift)):
            json["reference_alias"].append(
                dict(
                    name=tde.redshift_ref[0],
                    human_readable_name=otter.util.bibcode_to_hrn(tde.redshift_ref[0]),
                )
            )

        curr_refs = {v["name"] for v in json["reference_alias"]}
        if tde.discovery_date_ref[0] not in curr_refs and not np.any(
            pd.isna(tde.discovery_date)
        ):
            json["reference_alias"].append(
                dict(
                    name=tde.discovery_date_ref[0],
                    human_readable_name=otter.util.bibcode_to_hrn(
                        tde.discovery_date_ref[0]
                    ),
                )
            )

        # print(jj.dumps(json, indent=4))
        # print('#####################################################################')
        # print()
        # print()

        all_jsons.append(json)

    db.save(all_jsons, testing=False)


if __name__ == "__main__":
    main()
