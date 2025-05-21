"""
Convert the Dataset provided in Masterson, Megan et al. (2024) to the OTTER format
"""

import os
import glob
import otter
import pandas as pd
import numpy as np

from astropy import units as u
from astropy import constants as const

class_conf_map = {
    "WTP 14abnpgk": 3.2,  # silver
    "WTP 14acnjbu": 3.2,  # silver
    "WTP 14adbjsh": 3.3,  # gold
    "WTP 14adbwvs": 1,  # the spectrum in the paper is archival
    "WTP 14adeqka": 3.3,  # gold
    "WTP 15abymdq": 3.3,  # gold
    "WTP 15acbgpn": 1,  # uses archival spectrum
    "WTP 15acbuuv": 3.3,  # gold
    "WTP 16aaqrcr": 3.3,  # gold
    "WTP 16aatsnw": 1,  # archival spec
    "WTP 17aaldjb": 1,  # archival spec
    "WTP 17aalzpx": 3.3,  # gold
    "WTP 17aamoxe": 3.3,  # gold
    "WTP 17aamzew": 1,  # archival spec
    "WTP 17aanbso": 1,  # archival spec
    "WTP 18aajkmk": 3.3,  # gold
    "WTP 18aamced": 1,  # arcihval spec
    "WTP 18aampwj": 3.3,  # gold
}


# we need the minimum, maximum, and effective wavelengths
def xray_to_wave(xraycode):
    min_energy, max_energy = [float(val) for val in xraycode.split("-")]
    wave_eff = (const.h * const.c / (max_energy - min_energy) / u.keV).to(u.nm).value

    # these variable names seem backwards but remember that energy and
    # wavelength are inversely proportional!
    wave_min = (const.h * const.c / max_energy / u.keV).to(u.nm).value
    wave_max = (const.h * const.c / min_energy / u.keV).to(u.nm).value

    return {
        "filter_key": xraycode,
        "filter_name": xraycode,
        "wave_eff": wave_eff,
        "wave_min": wave_min,
        "wave_max": wave_max,
        "wave_units": "nm",
    }


def pkl_to_csv(pklfile, *args, **kwargs):
    df = pd.read_pickle(pklfile, *args, **kwargs)
    outpath = pklfile.replace("pkl", "csv")
    df.to_csv(outpath, index=False)


def clean_meta_data(path):
    meta = pd.read_csv(path, sep="\t")
    del meta["Unnamed: 9"]

    # replace with appropriate bibcodes
    meta["bibcode"] = meta["First Presented"]
    del meta["First Presented"]
    meta["bibcode"] = meta.bibcode.str.replace("This work", "2024ApJ...961..211M")
    meta["bibcode"] = meta.bibcode.str.replace(
        "Panagiotou et al. (2023)", "2023ApJ...948L...5P"
    )
    meta["bibcode"] = meta.bibcode.str.replace(
        "Dai et al. (2020)", "2020ApJ...896L..27D"
    )
    meta["bibcode"] = meta.bibcode.str.replace(
        "Jiang et al. (2021a)", "2021ApJS..252...32J"
    )
    meta["bibcode"] = meta.bibcode.str.replace(
        "Wang et al. (2022b)", "2022ApJS..258...21W"
    )
    meta["bibcode"] = meta.bibcode.str.replace(
        "Arcavi et al. (2018)", "2018ATel11953....1A"
    )
    meta["bibcode"] = meta.bibcode.str.replace(
        "Frederick et al. (2019)", "2019ApJ...883...31F"
    )
    meta["bibcode"] = meta.bibcode.str.replace(
        "Falco et al. (2018)", "2018TNSCR..81....1F"
    )
    meta["bibcode"] = meta.bibcode.str.split("; ")
    meta["bibcode"] = [
        blist + ["2024ApJ...961..211M"] if "2024ApJ...961..211M" not in blist else blist
        for blist in meta.bibcode
    ]

    return meta


def clean_uvoir_data(csvdir):
    # ASASSN Data First
    asassn = pd.read_csv(os.path.join(csvdir, "Masterson24_asassn_opt.csv"))
    asassn_refflux = pd.read_csv(
        os.path.join(csvdir, "Masterson24_asassn_referenceflux.csv")
    )

    # reformat the refflux
    out_refflux = {"name": [], "band": [], "refflux": [], "refflux_unit": []}
    for i, row in asassn_refflux.iterrows():
        out_refflux["name"] += [row["name"], row["name"]]
        out_refflux["band"] += ["g", "V"]
        out_refflux["refflux"] += [row.g_refflux, row.V_refflux]
        out_refflux["refflux_unit"] += ["mJy", "mJy"]

    asassn_refflux_reformated = pd.DataFrame(out_refflux)

    asassn = asassn.merge(asassn_refflux_reformated, on=["name", "band"])

    # now clean up the asassn data a little
    asassn["upperlimit"] = asassn.apply(lambda row: row.flux_err == -99, axis=1)

    # convert upperlimits from 10 sigma to 3 sigma
    asassn["flux_err"] = asassn.apply(
        lambda row: 0 if row.upperlimit else row.flux_err, axis=1
    )
    asassn["flux"] = asassn.apply(
        lambda row: 3 * (row.flux / 10) if row.upperlimit else row.flux, axis=1
    )

    # add some descriptive columns
    asassn["time_format"] = "decimalyear"
    asassn["flux_unit"] = "mJy"

    asassn["telescope"] = "ASASSN"

    # WISE Next
    wise = pd.read_csv(os.path.join(csvdir, "Masterson24_wise_ir.csv"))
    wise_refflux = pd.read_csv(
        os.path.join(csvdir, "Masterson24_wise_referenceflux.csv")
    )

    # reformat the refflux
    out_refflux = {"name": [], "band": [], "refflux": [], "refflux_unit": []}
    for i, row in wise_refflux.iterrows():
        out_refflux["name"] += [row["name"], row["name"]]
        out_refflux["band"] += ["W1", "W2"]
        out_refflux["refflux"] += [row.W1_refflux, row.W2_refflux]
        out_refflux["refflux_unit"] += ["mJy", "mJy"]

    wise_refflux_reformated = pd.DataFrame(out_refflux)

    wise = wise.merge(wise_refflux_reformated, on=["name", "band"])

    # now clean up the wise data a little
    wise["upperlimit"] = wise.apply(lambda row: row.flux_err == -99, axis=1)

    # convert upperlimits from 10 sigma to 3 sigma
    wise["flux_err"] = wise.apply(
        lambda row: 0 if row.upperlimit else row.flux_err, axis=1
    )
    wise["flux"] = wise.apply(
        lambda row: 3 * (row.flux / 10) if row.upperlimit else row.flux, axis=1
    )

    # add some descriptive columns
    wise["time_format"] = "mjd"
    wise["flux_unit"] = "mJy"

    wise["telescope"] = "WISE"

    # combine and return the datasets
    return pd.concat([asassn, wise])


def clean_xray_data(csvdir):
    # eROSITA Data
    ero = pd.read_csv(os.path.join(csvdir, "Masterson24_eRO_Xray.csv"))

    out_ero = {"name": [], "band": [], "date": [], "flux": [], "flux_err": []}
    for _, row in ero.iterrows():
        out_ero["name"] += [row["name"]] * 5
        out_ero["band"] += [row.band] * 5
        out_ero["date"] += [row[f"year_{i}"] for i in range(1, 6)]
        out_ero["flux"] += [row[f"flux_{i}"] for i in range(1, 6)]
        out_ero["flux_err"] += [row[f"flux_err_{i}"] for i in range(1, 6)]

    ero = pd.DataFrame(out_ero)
    ero = ero[~pd.isna(ero.date)]

    ero["upperlimit"] = ero.apply(lambda row: row.flux_err == -1, axis=1)
    ero["flux_err"] = ero.apply(
        lambda row: 0 if row.upperlimit else row.flux_err, axis=1
    )
    ero["flux_unit"] = "erg/s/cm^2"
    ero["date_format"] = "decimalyear"

    ero["band_name"] = ero.band.str.replace(" keV", "")
    del ero["band"]

    ero["telescope"] = "eROSITA"

    # Chandra Data
    chandra = pd.read_csv(os.path.join(csvdir, "Masterson24_chandra_Xray.csv"))

    chandra["flux_err"] = np.mean(chandra[["flux_lowerr", "flux_uperr"]], axis=1)
    chandra["upperlimit"] = chandra.apply(lambda row: row.flux_err == -1, axis=1)
    chandra["flux_err"] = chandra.apply(
        lambda row: 0 if row.upperlimit else row.flux_err, axis=1
    )
    chandra["flux_lowerr"] = chandra.apply(
        lambda row: 0 if row.upperlimit else row.flux_lowerr, axis=1
    )
    chandra["flux_uperr"] = chandra.apply(
        lambda row: 0 if row.upperlimit else row.flux_uperr, axis=1
    )
    chandra["flux_unit"] = "erg/s/cm^2"

    chandra["date"] = chandra.mjd
    chandra["date_format"] = "mjd"
    del chandra["mjd"]

    chandra["band_name"] = chandra.band.str.replace("keV", "")
    del chandra["band"]

    chandra["telescope"] = "Chandra"

    # Swift Data
    swift = pd.read_csv(os.path.join(csvdir, "Masterson24_swift_Xray.csv"))

    swift["flux_err"] = np.mean(swift[["flux_lowerr", "flux_uperr"]], axis=1)
    swift["upperlimit"] = swift.apply(lambda row: row.flux_err == -1, axis=1)
    swift["flux_err"] = swift.apply(
        lambda row: 0 if row.upperlimit else row.flux_err, axis=1
    )
    swift["flux_lowerr"] = swift.apply(
        lambda row: 0 if row.upperlimit else row.flux_lowerr, axis=1
    )
    swift["flux_uperr"] = swift.apply(
        lambda row: 0 if row.upperlimit else row.flux_uperr, axis=1
    )
    swift["flux_unit"] = "erg/s/cm^2"

    swift["date"] = swift.mjd
    swift["date_format"] = "mjd"
    del swift["mjd"]

    swift["band_name"] = swift.band.str.replace("keV", "")
    del swift["band"]

    swift["telescope"] = "Swift"

    # XMM Data
    xmm = pd.read_csv(os.path.join(csvdir, "Masterson24_xmm_Xray.csv"))

    xmm["upperlimit"] = xmm.apply(lambda row: row.flux_err == -1, axis=1)
    xmm["flux_err"] = xmm.apply(
        lambda row: 0 if row.upperlimit else row.flux_err, axis=1
    )
    xmm["flux_unit"] = "erg/s/cm^2"

    xmm["date"] = xmm.time
    xmm["date_format"] = "decimalyear"
    del xmm["time"]

    xmm["band_name"] = xmm.band.str.replace("keV", "")
    del xmm["band"]

    xmm["telescope"] = "XMM-Slew"

    # merge and return the data
    return pd.concat([ero, chandra, swift, xmm])


def main():
    import argparse

    pp = argparse.ArgumentParser()
    pp.add_argument("--otterdir", help="Directory where the otter json files will go")
    pp.add_argument("--indir", help="Directory where dirty files are")
    args = pp.parse_args()

    db = otter.Otter(datadir=args.otterdir, gen_summary=True)

    # define some paths we need
    pklfiles = glob.glob(os.path.join(args.indir, "pkl", "*.pkl"))
    csvdir = os.path.join(args.indir, "csv")
    metapath = os.path.join(args.indir, "meta.txt")

    # make sure all of the pickle files are converted
    for f in pklfiles:
        pkl_to_csv(f)

    # get the three datasets
    meta = clean_meta_data(metapath)
    xray = clean_xray_data(csvdir)
    uvoir = clean_uvoir_data(csvdir)

    # reformat the data
    all_jsons = []
    for idx, row in meta.iterrows():
        otter_json = {}

        # name info
        tns_used = False
        if not pd.isna(row["TNS Name"]):
            otter_json["name"] = dict(
                default_name=row["TNS Name"].replace(" ", ""),
                alias=[
                    dict(value=row["TNS Name"].replace(" ", ""), reference="TNS"),
                    dict(value=row["WTP Name"].replace(" ", ""), reference=row.bibcode),
                ],
            )
            tns_used = True
        else:
            otter_json["name"] = dict(
                default_name=row["WTP Name"].replace(" ", ""),
                alias=[
                    dict(value=row["WTP Name"].replace(" ", ""), reference=row.bibcode)
                ],
            )

        # coordinates
        otter_json["coordinate"] = [
            dict(
                ra=row["R.A."],
                dec=row["Decl."],
                ra_units="deg",
                dec_units="deg",
                reference=row.bibcode,
                coordinate_type="equatorial",
            )
        ]

        # distance
        otter_json["distance"] = [
            dict(value=row.Redshift, reference=row.bibcode, distance_type="redshift"),
            dict(
                value=row["D _L"],
                reference=row.bibcode,
                unit="Mpc",
                distance_type="luminosity",
            ),
        ]

        # date_reference
        otter_json["date_reference"] = [
            dict(
                value=row["t _disrupt"],
                reference=row.bibcode,
                date_format="mjd",
                date_type="explosion",
            )
        ]

        # classification
        # some of these have spectra post-tde, use the map at the top of this script
        # to flag these classifications
        otter_json["classification"] = dict(
            value=[
                dict(
                    object_class="TDE",
                    confidence=class_conf_map[row["WTP Name"]],
                    reference=row.bibcode,
                )
            ]
        )

        # photometry
        # uvoir first
        otter_json["photometry"] = []
        wtp_name = row["WTP Name"].replace(" ", "")
        uvoir_filters_used = []
        for keys, grp in uvoir[uvoir["name"] == wtp_name].groupby("telescope"):
            otter_json["photometry"].append(
                dict(
                    date=grp.time.tolist(),
                    date_format=grp.time_format.tolist(),
                    filter_key=grp.band.tolist(),
                    raw=grp.flux.tolist(),
                    raw_err=grp.flux_err.tolist(),
                    raw_units=grp.flux_unit.tolist(),
                    telescope=keys,
                    upperlimit=grp.upperlimit.tolist(),
                    corr_k=False,
                    corr_s=False,
                    corr_host=True,
                    corr_av=False,
                    corr_hostav=False,
                    val_host=grp.refflux.tolist(),
                    obs_type="uvoir",
                    reference=["2024ApJ...961..211M"],
                )
            )

            for filt in grp.band.unique():
                if filt not in uvoir_filters_used:
                    uvoir_filters_used.append(filt)

        # now x ray photometry
        xray_filters_used = []

        xray_models = {
            "eROSITA": dict(
                model_name="Absorbed Powerlaw",
                param_names=["Gamma", "N_H"],
                param_values=[2.0, 3e20],
                param_units=["None", "cm^-2"],
                min_energy=0.2,
                max_energy=2.3,
                energy_units="keV",
                model_reference="2024ApJ...961..211M",
            ),
            "Chandra": dict(
                model_name="Powerlaw",
                param_names=["Gamma"],
                param_values=[2.0],
                param_units=["None"],
                min_energy=0.5,
                max_energy=7,
                energy_units="keV",
                model_reference="2024ApJ...961..211M",
            ),
            "Swift": dict(
                model_name="Powerlaw",
                param_names=["Gamma"],
                param_values=[2.0],
                param_units=["None"],
                min_energy=0.3,
                max_energy=10,
                energy_units="keV",
                model_reference="2024ApJ...961..211M",
            ),
            "XMM-Slew": dict(
                model_name="Powerlaw",
                param_names=["Gamma"],
                param_values=[2.0],
                param_units=["None"],
                min_energy=0.2,
                max_energy=2,
                energy_units="keV",
                model_reference="2024ApJ...961..211M",
            ),
        }

        for keys, grp in xray[xray["name"] == wtp_name].groupby("telescope"):
            to_app = dict(
                date=grp.date.tolist(),
                date_format=grp.date_format.tolist(),
                filter_key=grp.band_name.tolist(),
                raw=grp.flux.tolist(),
                raw_err=grp.flux_err.tolist(),
                raw_units=grp.flux_unit.tolist(),
                telescope=keys,
                upperlimit=grp.upperlimit.tolist(),
                corr_k=False,
                corr_s=False,
                corr_host=False,
                corr_av=False,
                corr_hostav=False,
                obs_type="xray",
                xray_model=[xray_models[keys]] * len(grp),
                reference=["2024ApJ...961..211M"],
            )

            if not np.all(pd.isna(grp.flux_uperr)) and not np.all(
                pd.isna(grp.flux_lowerr)
            ):
                to_app["raw_err_detail"] = dict(
                    upper=grp.flux_uperr.tolist(), lower=grp.flux_lowerr.tolist()
                )

            otter_json["photometry"].append(to_app)

            for filt in grp.band_name.unique():
                if filt not in xray_filters_used:
                    xray_filters_used.append(filt)

        # deal with the filter alias
        otter_json["filter_alias"] = []
        for filt in uvoir_filters_used:
            otter_json["filter_alias"].append(
                dict(
                    filter_key=filt,
                    filter_name=filt,
                    wave_eff=otter.util.FILTER_MAP_WAVE[filt],
                    wave_units="nm",
                )
            )

        for filt in xray_filters_used:
            otter_json["filter_alias"].append(xray_to_wave(filt))

        # reference alias
        uq_bibcodes, all_hrns = otter.util.bibcode_to_hrn(row.bibcode)

        # package these into the reference alias
        otter_json["reference_alias"] = [
            dict(name=name, human_readable_name=hrn)
            for name, hrn in zip(uq_bibcodes, all_hrns)
        ]

        if tns_used:
            otter_json["reference_alias"].append(
                dict(name="TNS", human_readable_name="TNS")
            )

        # print(json.dumps(otter_json, indent=4))
        # print('####################################################################')
        # print()
        # print()

        all_jsons.append(otter_json)

    db.save(all_jsons, testing=False)


if __name__ == "__main__":
    main()
