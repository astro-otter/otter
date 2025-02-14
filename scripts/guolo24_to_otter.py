"""
Guolo+24 Xray modeling paper to the Otter format
"""

import os
import otter
import re

import pandas as pd
import numpy as np
from astropy.table import Table
from astropy.io import fits
from astropy import units as u
from astropy import constants as const

guolo_bibcode = "2024ApJ...966..160G"

# from the paper text
xraycodes_reduction = {
    "Swift/XRT": "0.3 - 10",
    "XMM-Newton": "0.3 - 10.0",
    "eRO": "0.3 - 2.0",
}

xraycodes = {"XRT": "0.3 - 10", "XMM": "0.2 - 12", "eRO": "0.3 - 2.0"}

# generate a mapping to be used as we convert to OTTER
tele_map = {"XRT": "Swift", "XMM": "XMM-Newton", "eRO": "eROSITA"}

name_map = {
    "at2018zr": "AT 2018zr",
    "at2018hyz": "AT 2018hyz",
    "at2019azh": "AT 2019azh",
    "at2019dsg": "AT 2019dsg",
    "at2019ehz": "AT 2019ehz",
    "at2019qiz": "",  # just don't use this one, use the one from Matt Nicholl's paper
    "at2019teq": "AT 2019teq",
    "at2019vcb": "AT 2019vcb",
    "at2020ddv": "AT 2020ddv",
    "at2020ksf": "AT 2020ksf",
    "at2020ocn": "AT 2020ocn",
    "at2021ehb": "AT 2021ehb",
    "at2021yzv": "AT 2021yzv",
    "asassn-14li": "ASASSN-14li",
    "asassn-15oi": "ASASSN-15oi",
    "at2018fyk": "AT 2018fyk",
}

tele_map_models = {"XRT": "Swift/XRT", "XMM": "XMM-Newton"}


def get_param(pstr):
    pstr = pstr.replace("(fixed)", "").strip()

    regex_pattern = r"[-+]?\d*\.?\d+"
    regex_result = re.findall(regex_pattern, pstr)

    params_and_errs = [float(v) for v in regex_result[:3]]
    if len(regex_result) < 4:
        return params_and_errs

    f = 10 ** float(regex_result[-1])
    return [v * f for v in params_and_errs]


def get_closest_model(date, name, tele, model_params, refs):
    # rint(name, date, tele)

    # read in clean the model file

    # find the closest model to the input data points
    test_set = model_params[
        (model_params.Source == name) * (model_params.Instrument == tele)
    ]
    closest_model_idx = np.argmin(np.abs(test_set.center_mjd.astype(float) - date))

    closest_model = test_set.iloc[closest_model_idx]

    # two possible models, setup to parse both
    min_energy, max_energy = xraycodes_reduction[tele].split(" - ")
    if closest_model["f _sc"] == "cdots":
        T_p, T_p_lo, T_p_up = get_param(closest_model["T _p (K)"])  # noqa: N806
        R_p, R_p_lo, R_p_up = get_param(closest_model["R _p (cm)"])  # noqa: N806

        model_dict = dict(
            model_name="tdediscspec",
            param_names=["T_p", "R_p"],
            param_values=[T_p, R_p],
            param_value_upper_err=[T_p_up, R_p_up],
            param_value_lower_err=[T_p_lo, R_p_lo],
            param_upperlimit=[False, False],
            param_units=["K", "cm"],
            min_energy=float(min_energy),
            max_energy=float(max_energy),
            energy_units="keV",
            model_referece=["2023MNRAS.519.5828M", guolo_bibcode],
        )
    else:
        T_p, T_p_lo, T_p_up = get_param(closest_model["T _p (K)"])  # noqa: N806
        R_p, R_p_lo, R_p_up = get_param(closest_model["R _p (cm)"])  # noqa: N806
        f_sc, f_sc_lo, f_sc_up = get_param(closest_model["f _sc"])
        y_sc, y_sc_lo, y_sc_up = get_param(closest_model["Gamma_sc"])

        model_dict = dict(
            model_name="simPL+tdediscspec",
            param_names=["T_p", "R_p", "f_sc", "Gamma_sc"],
            param_values=[T_p, R_p, f_sc, y_sc],
            param_value_upper_err=[T_p_up, R_p_up, f_sc_up, y_sc_up],
            param_value_lower_err=[T_p_lo, R_p_lo, f_sc_lo, y_sc_lo],
            param_upperlimit=[False, False, False, False],
            param_units=["K", "cm", None, None],
            min_energy=float(min_energy),
            max_energy=float(max_energy),
            energy_units="keV",
            model_referece=["2023MNRAS.519.5828M", guolo_bibcode],
        )

        refs.append("2023MNRAS.519.5828M")

    if name in {"ASASSN-14li", "AT 2019dsg"}:
        model_dict["model_name"] += "+absorption"
        model_dict["param_names"].append("N_H")
        model_dict["param_value_upper_err"].append(0)
        model_dict["param_value_lower_err"].append(0)
        model_dict["param_upperlimit"].append(False)
        model_dict["param_units"].append("cm^-2")

        if name == "ASASSN-14li":
            model_dict["param_values"].append(5e20)
        if name == "AT 2019dsg":
            model_dict["param_values"].append(3e20)

    return model_dict


# we need the minimum, maximum, and effective wavelengths
def xray_to_wave(xraycode):
    min_energy, max_energy = [float(val) for val in xraycode.split(" - ")]
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


filter_mapping = {
    "XRT": xray_to_wave(xraycodes["XRT"]),
    "XMM": xray_to_wave(xraycodes["XMM"]),
    "eRO": xray_to_wave(xraycodes["eRO"]),
}


def main():
    import argparse

    pp = argparse.ArgumentParser()
    pp.add_argument("--otterdir", help="Directory where the otter json files will go")
    pp.add_argument("--indir", help="Directory where dirty files are")
    args = pp.parse_args()

    model_params = pd.read_csv(
        os.path.join(args.indir, "guolo23_model_params.txt"), sep="\t"
    ).ffill()
    model_params["center_mjd"] = model_params.MJD.apply(
        lambda val: float(get_param(val)[0])
    )
    model_params.drop(
        model_params[
            (model_params.Source == "AT 2019teq") * (model_params.MJD == "58915")
        ].index,
        inplace=True,
    )

    phot_file = os.path.join(args.indir, "light_curves.fits")
    meta_file = os.path.join(args.indir, "source_info.fits")

    _, meta_table = fits.open(meta_file)
    meta = Table(meta_table.data).to_pandas()

    # all this data is already in OTTER so pull the coordinates from there
    db = otter.Otter(datadir=args.otterdir)
    meta_coord_data = {"name": [], "ra": [], "dec": [], "ref": []}
    for name in list(meta.name.unique()):
        res = db.query(names=name)
        if len(res) == 0:
            print(name, ", No coordinates found!")

        ra, dec = res[0].get_skycoord().to_string("hmsdms").split(" ")
        ref = res[0]._get_default("coordinate")["reference"][
            0
        ]  # the first one is fine for referencing

        meta_coord_data["name"].append(name)
        meta_coord_data["ra"].append(ra)
        meta_coord_data["dec"].append(dec)
        meta_coord_data["ref"].append(ref)

    meta_coord_data = pd.DataFrame(meta_coord_data)

    meta = meta.merge(meta_coord_data, on="name")

    all_lc_data = []
    for lc, m in zip(fits.open(phot_file)[1:], meta):
        lc_data = Table(lc.data).to_pandas()
        lc_data["name"] = lc.header["EXTNAME"]
        all_lc_data.append(lc_data)

    all_lc_data = pd.concat(all_lc_data)

    data = all_lc_data.merge(
        meta, left_on=all_lc_data.name.str.lower(), right_on=meta.name.str.lower()
    )

    # convert this to OTTER JSON format now
    problem_children = []
    new_data = []
    for tdename, grp in data.groupby("key_0"):
        # We will use a different X-ray dataset from Matt Nicholl for this one
        if tdename == "at2019qiz":
            continue

        otter_json = {}
        refs = [guolo_bibcode]

        # names
        otter_json["name"] = dict(
            default_name=tdename, alias=[dict(value=tdename, reference=[guolo_bibcode])]
        )

        possible_alias = grp.alt_name.values[0]
        if len(possible_alias) > 0:
            otter_json["name"]["alias"].append(
                dict(value=possible_alias, reference=[guolo_bibcode])
            )

        # coordinate
        otter_json["coordinate"] = [
            dict(
                ra=grp.ra.iloc[0],
                dec=grp.dec.iloc[0],
                ra_units="hour",
                dec_units="deg",
                reference=[grp.ref.iloc[0]],
                coordinate_type="equitorial",
            )
        ]
        refs.append(grp.ref.iloc[0])

        # redshift
        z = grp.z.iloc[0]
        otter_json["distance"] = [
            dict(value=float(z), reference=[guolo_bibcode], distance_type="redshift")
        ]

        # classification
        otter_json["classification"] = [
            dict(object_class="TDE", confidence=1, reference=[guolo_bibcode])
        ]

        # date peak
        otter_json["date_reference"] = [
            dict(
                value=float(grp.mjd_peak.iloc[0]),
                reference=[guolo_bibcode],
                date_type="peak",
                date_format="MJD",
            )
        ]

        # photometry
        otter_json["photometry"] = []
        otter_json["filter_alias"] = []

        for inst, subgrp in grp.groupby("instrument"):
            otter_json["photometry"].append(
                dict(
                    raw=subgrp.Fx.tolist(),
                    raw_err=list(
                        np.sqrt(subgrp.Fx_perr**2 + subgrp.Fx_nerr**2).fillna(0)
                    ),
                    raw_units=["erg cm^-2 s^-1"] * len(subgrp.Fx),
                    raw_err_detail={
                        "upper": subgrp.Fx_perr.fillna(0).tolist(),
                        "lower": subgrp.Fx_nerr.fillna(0).tolist(),
                    },
                    upperlimit=[pd.isna(v) for v in subgrp.Fx_perr],
                    date=subgrp.mjd.tolist(),
                    date_format=["mjd"] * len(subgrp),
                    filter_key=subgrp.instrument.replace(xraycodes).tolist(),
                    telescope=tele_map[inst],
                    instrument=inst,
                    reference=[guolo_bibcode]
                    if inst != "eRO"
                    else [guolo_bibcode, "2020ATel14246....1G"],
                    obs_type="xray",
                    xray_model=[],
                )
            )

            if inst == "eRO":
                refs.append("2020ATel14246....1G")

            for time in subgrp.mjd:
                try:
                    otter_json["photometry"][-1]["xray_model"].append(
                        get_closest_model(
                            date=time,
                            name=name_map[tdename],
                            tele=tele_map_models[inst],
                            model_params=model_params,
                            refs=refs,
                        )
                    )
                except (ValueError, KeyError):
                    otter_json["photometry"][-1]["xray_model"].append(None)
                    problem_children.append((tdename, time, inst))

            otter_json["filter_alias"].append(filter_mapping[inst])

        # references
        ref_details = list(zip(*otter.util.bibcode_to_hrn(np.unique(refs))))
        otter_json["reference_alias"] = [
            dict(name=bib, human_readable_name=hrn) for bib, hrn in ref_details
        ]

        # print(json.dumps(otter_json, indent=4))
        # print()
        # print()

        # validate this dataset
        otter.schema.OtterSchema(**otter_json)

        new_data.append(otter.Transient(otter_json))

    missing_models = pd.DataFrame(
        problem_children, columns=["name", "mjd", "telescope"]
    )
    missing_models.to_csv("guolo23_missing_models.csv")

    db.save(new_data, testing=False)


if __name__ == "__main__":
    main()
