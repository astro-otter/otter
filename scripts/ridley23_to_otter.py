import os
import otter
import pandas as pd
import numpy as np

from astropy import units as u


def main():
    import argparse

    pp = argparse.ArgumentParser()
    pp.add_argument("--otterdir", help="Directory where the otter json files will go")
    pp.add_argument("--indir", help="Directory where dirty files are")
    args = pp.parse_args()

    # Metadata first
    ridley23_bibcode = "2024MNRAS.531.1905R"

    meta = dict(
        name=dict(
            default_name="AT2017bcc",
            alias=[dict(value="AT2017bcc", reference=[ridley23_bibcode])],
        ),
        coordinate=[
            dict(
                reference=["2017TNSTR.220....1C"],  # from TNS Discovery Report
                ra=172.970800581,
                dec=29.9958695719,
                ra_units="deg",
                dec_units="deg",
                default=True,
                coordinate_type="equatorial",
            )
        ],
        distance=[
            dict(
                value=0.148,
                reference=["2017TNSCR.221....1B"],  # from TNS Classification report
                distance_type="redshift",
            )
        ],
        classification=dict(
            value=[
                dict(
                    object_class="ANT",
                    confidence=3.3,  # they have a spectrum
                    default=True,
                    reference=ridley23_bibcode,
                )
            ]
        ),
        date_reference=[
            dict(
                value=2457802.9139931,
                date_format="JD",
                date_type="discovery",
                reference=["2017TNSTR.220....1C"],
            )
        ],
        reference_alias=[
            dict(name=ridley23_bibcode, human_readable_name="Ridley et al. (2023)"),
            dict(
                name="2017TNSTR.220....1C", human_readable_name="Chambers et al. (2017)"
            ),
            dict(
                name="2017TNSCR.221....1B", human_readable_name="Barbarino et al (2017)"
            ),
        ],
    )

    # uvoir dataset

    uvoir = pd.read_csv(os.path.join(args.indir, "at2017bcc_uvoir_phot.csv"))

    uvoir_phot = dict(
        reference=[ridley23_bibcode],
        raw=uvoir.mag.tolist(),
        raw_err=uvoir.magerr.tolist(),
        raw_units="mag(AB)",
        corr_k=False,
        corr_s=False,
        corr_av=False,
        corr_host=False,
        corr_hostav=False,
        upperlimit=[False] * len(uvoir),
        date=uvoir.mjd.tolist(),
        date_format="MJD",
        filter_key=uvoir["filter"].tolist(),
        obs_type="uvoir",
    )

    filts, idxs = np.unique(uvoir["filter"], return_index=True)
    filter_alias = [
        dict(filter_key=k, filter_name=k, wave_eff=w, wave_units="AA")
        for k, w in zip(filts, uvoir.wavelength.iloc[idxs])
    ]

    # test the otter dataset format based on our schema
    otter.schema.PhotometrySchema(**uvoir_phot)
    for f in filter_alias:
        otter.schema.FilterSchema(**f)

    # radio dataset
    radio = pd.read_csv(os.path.join(args.indir, "at2017bcc_radio_phot.csv"), sep="\t")

    radio_phot = dict(
        reference=[ridley23_bibcode],
        raw=radio["F(μJy)"].tolist(),
        raw_err=radio["dF(μJy)"].tolist(),
        raw_units="uJy",
        corr_k=False,
        corr_s=False,
        corr_av=False,
        corr_host=False,
        corr_hostav=False,
        upperlimit=[False] * len(radio),
        date=radio.MJD.tolist(),
        date_format="MJD",
        filter_key=(
            radio["Frequency(GHz)"].astype(str) + ["GHz"] * len(radio)
        ).tolist(),
        obs_type="radio",
        telescope="VLA",
    )

    filter_alias += [
        dict(
            filter_key=str(f) + "GHz",
            filter_name=otter.util.freq_to_band(f * u.GHz),
            freq_eff=f,
            freq_units="GHz",
        )
        for f in radio["Frequency(GHz)"].unique()
    ]

    # test the otter dataset format based on our schema
    otter.schema.PhotometrySchema(**radio_phot)
    for f in filter_alias:
        otter.schema.FilterSchema(**f)

    # X-ray dataset
    xray = pd.read_csv(os.path.join(args.indir, "at2017bcc_xray_phot.csv"), sep=",")

    xray_phot = dict(
        reference=[ridley23_bibcode],
        raw=xray["Counts[s^−1]"].tolist(),
        raw_err=((xray.Counts_upper + xray.Counts_lower.abs()) / 2).tolist(),
        raw_units="ct",
        corr_k=False,
        corr_s=False,
        corr_av=False,
        corr_host=False,
        corr_hostav=False,
        upperlimit=[False] * len(xray),
        date=xray.MJD.tolist(),
        date_format="MJD",
        filter_key="0.3 - 10",
        obs_type="xray",
        telescope="swift",
        raw_err_detail=dict(
            upper=xray.Counts_upper.tolist(), lower=xray.Counts_lower.tolist()
        ),
        xray_model=[
            dict(
                model_name="Absorbed Powerlaw",
                param_names=["Gamma", "N_H"],
                param_values=[1.54, 1.9e20],
                param_units=["None", "cm^-2"],
                param_value_err_upper=[0.14, 0],
                param_value_err_lower=[-0.11, 0],
                model_reference=ridley23_bibcode,
                min_energy=0.3,
                max_energy=10,
                energy_units="keV",
            )
            for _ in range(len(xray))
        ],
    )

    filter_alias += [
        dict(
            filter_key="0.3 - 10",
            filter_name="0.3 - 10",
            wave_eff=((10 - 0.3) * u.keV).to(u.nm, equivalencies=u.spectral()).value,
            wave_units="nm",
            wave_max=(0.3 * u.keV).to(u.nm, equivalencies=u.spectral()).value,
            wave_min=(10 * u.keV).to(u.nm, equivalencies=u.spectral()).value,
        )
    ]

    # test the otter dataset format based on our schema
    otter.schema.PhotometrySchema(**xray_phot)
    for f in filter_alias:
        otter.schema.FilterSchema(**f)

    # combine everything and save it
    otter_json = meta | {
        "photometry": [uvoir_phot, radio_phot, xray_phot],
        "filter_alias": filter_alias,
    }

    otter.schema.OtterSchema(**otter_json)  # checks the validity of the json file

    t = otter.Transient(otter_json)
    db = otter.Otter(datadir=args.otterdir, gen_summary=True)
    db.save([t], testing=False)


if __name__ == "__main__":
    main()
