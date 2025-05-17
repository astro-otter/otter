"""
Convert the Radio Photometry that Noah, Collin, and Kate gathered to the Otter format
"""

# import
import os
import argparse
import otter
import pandas as pd
import numpy as np
from astropy import units as u


def main():
    pp = argparse.ArgumentParser()
    pp.add_argument("--otterdir", help="Directory where the otter json files will go")
    pp.add_argument("--indir", help="Directory where dirty files are")
    pp.add_argument("--debug", action=argparse.BooleanOptionalAction)
    args = pp.parse_args()

    db = otter.Otter(datadir=args.otterdir, gen_summary=True)

    # read in the metadata and photometry files
    meta = pd.read_csv(os.path.join(args.indir, "meta.csv"))
    phot = pd.read_csv(os.path.join(args.indir, "photometry.csv"))

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

        class_flag = 0
        if "classification_flag" in tde:
            class_flag = tde.classification_flag[0]

        json["classification"] = dict(
            value=dict(
                object_class="TDE",
                confidence=class_flag,
                reference=[tde.bibcode[0]],
            )
        )

        # discovery date
        # print(tde)
        if not np.any(pd.isna(tde.discovery_date)):
            json["date_reference"] = [
                dict(
                    value=str(tde.discovery_date.tolist()[0]).strip(),
                    date_format=tde.discovery_date_format.tolist()[0].lower(),
                    reference=tde.discovery_date_ref.tolist(),
                    computed=False,
                    date_type="discovery",
                )
            ]

        # host information
        if not np.any(pd.isna(tde.host_ref)):
            host_info = dict(
                host_name=tde.host_name.tolist()[0].strip(),
                host_ra=tde.host_ra.tolist()[0],
                host_dec=tde.host_dec.tolist()[0],
                host_ra_units=tde.host_ra_unit.tolist()[0],
                host_dec_units=tde.host_dec_unit.tolist()[0],
                reference=[tde.host_ref.tolist()[0]],
            )

            if not pd.isna(tde.host_redshift.tolist()[0]):
                host_info["host_z"] = tde.host_redshift.tolist()[0]

            if "host" in json:
                json["host"].append(host_info)
            else:
                json["host"] = [host_info]

        # now the radio photometry
        tde["obs_type"] = [
            otter.util.freq_to_obstype(vv * u.Unit(uu))
            for vv, uu in zip(
                tde.band_eff_freq.astype(float).values, tde.band_eff_freq_unit.values
            )
        ]

        phot_sources = []
        unique_filter_keys = []
        index_for_match = []
        json["photometry"] = []
        for (src, tele, obstype), p in tde.groupby(
            ["bibcode", "telescope", "obs_type"], dropna=False
        ):
            if src not in phot_sources:
                phot_sources.append(src)

            if len(np.unique(p.flux_unit)) == 1:
                raw_units = p.flux_unit.tolist()[0]
            else:
                raw_units = p.flux_unit.values

            # add a column to phot with the unique filter key
            if obstype == "radio":
                filter_uq_key = (
                    p.band_eff_freq.astype(str) + p.band_eff_freq_unit
                ).tolist()
            elif obstype == "uvoir":
                filter_uq_key = p.band_name.astype(str).tolist()
            else:
                raise ValueError("not prepared for this obstype!")

            unique_filter_keys += filter_uq_key
            index_for_match += p.index.tolist()

            json_phot = dict(
                reference=src,
                raw=p.flux.astype(float).tolist(),
                raw_err=p.flux_err.astype(float).tolist(),
                raw_units=raw_units,
                date=p.date.tolist(),
                date_format=p.date_format.tolist(),
                upperlimit=p.upperlimit.tolist(),
                filter_key=filter_uq_key,
                obs_type=obstype,
            )

            if not pd.isna(tele):
                json_phot["telescope"] = tele

            raw_err_detail = {}
            for key in ["statistical_err", "systematic_err", "iss_err"]:
                if not np.all(pd.isna(p[key])):
                    k = key.split("_")[0]

                    # fill the nan values
                    # this is to match with the official json format
                    # and works with arangodb document structure
                    p[key].fillna(0, inplace=True)

                    raw_err_detail[k] = p[key].tolist()

            if len(raw_err_detail) > 0:
                json_phot["raw_err_detail"] = raw_err_detail

            corrs = ["corr_k", "corr_s", "corr_host", "corr_av", "corr_hostav"]
            for c in corrs:
                json_phot[c] = False if np.all(pd.isna(p[c])) else p[c].tolist()
                if np.any(json_phot[c]):
                    v = c.replace("corr", "val")

                    # fill the nan values
                    # this is to match with the official json format
                    # and works with arangodb document structure
                    p[v].fillna("null", inplace=True)

                    json_phot[v] = p[v].tolist()

            json["photometry"].append(json_phot)

        tde["filter_uq_key"] = pd.Series(unique_filter_keys, index=index_for_match)

        # filter alias
        # radio filters first
        filter_keys1 = ["filter_uq_key", "band_eff_freq", "band_eff_freq_unit"]
        filter_map = (
            tde[filter_keys1].drop_duplicates().set_index("filter_uq_key")
        )  # .to_dict(orient='index')
        try:
            filter_map_radio = filter_map.to_dict(orient="index")
        except Exception:
            print(filter_map)
            print(name)
            raise Exception

        json["filter_alias"] = []
        for filt, val in filter_map_radio.items():
            obs_type = otter.util.freq_to_obstype(
                float(val["band_eff_freq"]) * u.Unit(val["band_eff_freq_unit"])
            )
            if obs_type == "radio":
                filter_name = otter.util.freq_to_band(
                    float(val["band_eff_freq"]) * u.Unit(val["band_eff_freq_unit"])
                )
            else:
                filter_name = filt

            json["filter_alias"].append(
                dict(
                    filter_key=filt,
                    filter_name=filter_name,
                    freq_eff=float(val["band_eff_freq"]),
                    freq_units=val["band_eff_freq_unit"],
                )
            )

        # reference alias
        # gather all the bibcodes
        all_bibcodes = [tde.coord_bibcode[0]] + phot_sources
        if tde.redshift_ref[0] not in all_bibcodes and not np.any(
            pd.isna(tde.redshift)
        ):
            all_bibcodes.append(tde.redshift_ref[0])

        if tde.discovery_date_ref[0] not in all_bibcodes and not np.any(
            pd.isna(tde.discovery_date)
        ):
            all_bibcodes.append(tde.discovery_date_ref[0])

        if tde.host_ref[0] not in all_bibcodes and not np.any(pd.isna(tde.host_ref)):
            all_bibcodes.append(tde.host_ref[0])

        # find the hrn's for all of these bibcodes
        uq_bibcodes, all_hrns = otter.util.bibcode_to_hrn(all_bibcodes)

        # package these into the reference alias
        json["reference_alias"] = [
            dict(name=name, human_readable_name=hrn)
            for name, hrn in zip(uq_bibcodes, all_hrns)
        ]

        # print(jj.dumps(json, indent=4))
        # print('#####################################################################')
        # print()
        # print()

        all_jsons.append(json)

    db.save(all_jsons, testing=args.debug)


if __name__ == "__main__":
    main()
