import os
import pandas as pd
import numpy as np
import otter

from astropy import units as u

from dataclasses import dataclass


@dataclass
class Param:
    name: float | str | int
    unit: str
    description: str = None


month_map = dict(
    Jan="01",
    Feb="02",
    Mar="03",
    Apr="04",
    May="05",
    Jun="06",
    Jul="07",
    Aug="08",
    Sep="09",
    Oct="10",
    Nov="11",
    Dec="12",
)

hajela24_bibcode = "2024arXiv240719019H"

radio_srcmap = {
    "this": [hajela24_bibcode],
    "VLASS": [hajela24_bibcode, "2020PASP..132c5001L"],
}

xray_srcmap = {
    "this": [hajela24_bibcode],
    "Gezari17": [hajela24_bibcode, "2017ApJ...851L..47G"],
    "Holoien18": [hajela24_bibcode, "2016MNRAS.463.3813H"],
}

band_map = {"swift": "0.3 - 10", "XMM-Newton": "0.2 - 12"}


gPL = Param(  # noqa: N816
    name="gamma_PL", unit="None", description="Photon Index for Powerlaw model"
)
FPL = Param(  # noqa: N816
    name="F_PL",
    unit="10^-14 erg/cm^2/s",
    description="Unabsorbed flux for Powerlaw model",
)
kT = Param(  # noqa: N816
    name="kT_BB", unit="10^-2 keV", description="k*T for the Blackbody model"
)
FBB = Param(  # noqa: N816
    name="F_BB",
    unit="10^-13 erg/cm^2/s",
    description="Unabsorbed flux for Blackbody model",
)
RBB = Param(  # noqa: N816
    name="R_BB", unit="10^12 cm", description="Blackbody radius"
)
all_ = [gPL, FPL, kT, FBB, RBB]


def compute_flux(row):
    if (
        not pd.isna(row.F_PL)
        and not pd.isna(row.F_BB)
        and not row.F_PL_upperlimit
        and not row.F_BB_upperlimit
    ):
        return (
            row.F_PL * 1e-14 + row.F_BB * 1e-13,
            (
                np.sqrt((row.F_PL_upper * 1e-14) ** 2 + (row.F_BB_upper * 1e-13) ** 2)
                + np.sqrt((row.F_PL_lower * 1e-14) ** 2 + (row.F_BB_lower * 1e-13) ** 2)
            )
            / 2,
            False,
            np.sqrt((row.F_PL_upper * 1e-14) ** 2 + (row.F_BB_upper * 1e-13) ** 2),
            np.sqrt((row.F_PL_lower * 1e-14) ** 2 + (row.F_BB_lower * 1e-13) ** 2),
        )

    elif (
        (pd.isna(row.F_PL) or row.F_PL_upperlimit)
        and not pd.isna(row.F_BB)
        and not row.F_BB_upperlimit
    ):
        return (
            row.F_BB * 1e-13,
            (row.F_BB_upper + abs(row.F_BB_lower)) / 2 * 1e-13,
            False,
            row.F_BB_upper,
            abs(row.F_BB_lower),
        )

    elif (
        (pd.isna(row.F_BB) or row.F_BB_upperlimit)
        and not pd.isna(row.F_PL)
        and not row.F_PL_upperlimit
    ):
        return (
            row.F_PL * 1e-14,
            (row.F_PL_upper + abs(row.F_PL_lower)) / 2 * 1e-14,
            False,
            row.F_PL_upper,
            abs(row.F_PL_lower),
        )

    elif row.F_PL_upperlimit:
        return row.F_PL * 1e-14, 0, True, 0, 0

    else:
        print(row)
        raise ValueError()


def main():
    import argparse

    pp = argparse.ArgumentParser()
    pp.add_argument("--otterdir", help="Directory where the otter json files will go")
    pp.add_argument("--indir", help="Directory where dirty files are")
    args = pp.parse_args()

    meta = dict(
        name=dict(
            default_name="ASASSN-15oi",
            alias=[dict(value="ASASSN-15oi", reference=[hajela24_bibcode])],
        ),
        coordinate=[
            dict(
                reference=["2016MNRAS.463.3813H"],
                ra="20:39:09.12",
                dec="-30:45:20.84",
                ra_units="hourangle",
                dec_units="deg",
                default=True,
                coordinate_type="equatorial",
            )
        ],
        reference_alias=[
            dict(name=hajela24_bibcode, human_readable_name="Hajela et al. (2024)"),
            dict(
                name="2016MNRAS.463.3813H", human_readable_name="Holoien et al. (2016)"
            ),
            dict(name="2020PASP..132c5001L", human_readable_name="Lucy et al. (2020)"),
            dict(
                name="2017ApJ...851L..47G", human_readable_name="Gezari et al. (2017)"
            ),
        ],
    )

    # UV/Opt/IR Photometry
    uvoir = pd.read_csv(
        os.path.join(args.indir, "asassn15oi_uvoir.txt"), sep=" ", index_col=False
    )

    uvoir_phot = dict(
        reference=[hajela24_bibcode],
        raw=uvoir.Magnitude.tolist(),
        raw_err=uvoir.Magnitude_Error.tolist(),
        raw_units="mag(AB)",  # this is based on a table A1 note in the paper
        corr_k=False,
        corr_s=False,
        corr_av=True,
        corr_host=False,
        corr_hostav=False,
        val_av=0.185,
        upperlimit=[False] * len(uvoir),
        date=uvoir.MJD.tolist(),
        date_format="MJD",
        filter_key=uvoir["Filter"].tolist(),
        obs_type="uvoir",
    )

    filts, idxs = np.unique(uvoir["Filter"], return_index=True)
    print(filts)
    filter_alias = [
        dict(
            filter_key=f,
            filter_name=f,
            wave_eff=otter.util.FILTER_MAP_WAVE[f],
            wave_units="nm",
        )
        for f in filts
    ]

    # test the otter dataset format based on our schema
    otter.schema.PhotometrySchema(**uvoir_phot)
    for f in filter_alias:
        otter.schema.FilterSchema(**f)

    # Radio photometry
    radio = pd.read_csv(os.path.join(args.indir, "asassn15oi_radio.txt"), sep=" ")

    radio["iso"] = radio.Date.replace(month_map, regex=True)

    radio_phot = [
        dict(
            reference=radio_srcmap[src],
            raw=grp.Fν_mJy.tolist(),
            raw_err=grp.Fv_err.tolist(),
            raw_units="mJy",
            corr_k=False,
            corr_s=False,
            corr_av=False,
            corr_host=False,
            corr_hostav=False,
            upperlimit=grp.upperlimit.tolist(),
            date=grp.iso.tolist(),
            date_format="iso",
            filter_key=(grp.ν_GHz.astype(str) + ["GHz"] * len(grp)).tolist(),
            obs_type="radio",
            telescope=tele,
        )
        for (tele, src), grp in radio.groupby(["Observatory", "Source"])
    ]

    filter_alias += [
        dict(
            filter_key=str(f) + "GHz",
            filter_name=otter.util.freq_to_band(f * u.GHz),
            freq_eff=f,
            freq_units="GHz",
        )
        for f in radio.ν_GHz.unique()
    ]

    # test the otter dataset format based on our schema
    for r in radio_phot:
        otter.schema.PhotometrySchema(**r)
    for f in filter_alias:
        otter.schema.FilterSchema(**f)

    # Finally X-ray photometry
    xray = pd.read_csv(os.path.join(args.indir, "asassn15oi_xray.txt"), sep=" ")
    for p in all_:
        xray[[p.name, p.name + "_lower", p.name + "_upper"]].astype(float)

    xray["model_name"] = [
        "Absorbed Powerlaw + Blackbody",
        "Absorbed Powerlaw + Blackbody",
        "Absorbed Powerlaw + Blackbody",
        "Absorbed Powerlaw + Blackbody",
        "Absorbed Powerlaw + Blackbody",
        "Absorbed Powerlaw",
        "Absorbed Powerlaw",
        "Powerlaw + Blackbody",
        "Absorbed Powerlaw + Blackbody",
        "Absorbed Powerlaw + Blackbody",
        "Powerlaw + Blackbody",
        "Absorbed Powerlaw + Blackbody",
        "Absorbed Powerlaw + Blackbody",
        "Absorbed Powerlaw + Blackbody",
        "Absorbed Powerlaw",
    ]

    xray["param_names"] = [
        all_,
        all_,
        all_,
        all_,
        [gPL, FPL, kT, FBB],
        [gPL, FPL],
        [gPL, FPL],
        [gPL, kT, FBB, RBB],
        all_,
        all_,
        [gPL, kT, FBB, RBB],
        all_,
        all_,
        [gPL, FPL, kT, FBB],
        [gPL, FPL],
    ]

    (
        xray["flux"],
        xray["flux_err"],
        xray["upperlimit"],
        xray["flux_upper"],
        xray["flux_lower"],
    ) = list(zip(*xray.apply(compute_flux, axis=1).tolist()))

    # N_H = 5.6e20  # cm^-2

    xray_phot = [
        dict(
            reference=xray_srcmap[src],
            raw=grp.flux.tolist(),
            raw_err=grp.flux_err.tolist(),
            raw_units="erg/s/cm^2",
            corr_k=False,
            corr_s=False,
            corr_av=False,
            corr_host=False,
            corr_hostav=False,
            val_host=(grp.F_abs * 1e-14).fillna("null").tolist(),
            upperlimit=grp.upperlimit.tolist(),
            date=(grp.dt + 57248.0).tolist(),
            date_format="MJD",
            filter_key=band_map[tele],
            obs_type="xray",
            telescope=tele,
            raw_err_detail=dict(
                upper=xray.flux_upper.tolist(), lower=xray.flux_lower.tolist()
            ),
            xray_model=[
                dict(
                    model_name=row.model_name,
                    param_names=[p.name for p in row.param_names],
                    param_values=[row[name.name] for name in row.param_names],
                    param_units=[p.unit for p in row.param_names],
                    param_value_err_upper=[
                        row[name.name + "_upper"] for name in row.param_names
                    ],
                    param_value_err_lower=[
                        row[name.name + "_lower"] for name in row.param_names
                    ],
                    param_upperlimit=[
                        row[name.name + "_upperlimit"] for name in row.param_names
                    ],
                    param_descriptions=[p.description for p in row.param_names],
                    model_reference=[hajela24_bibcode],
                    min_energy=float(band_map[tele].split("-")[0].strip()),
                    max_energy=float(band_map[tele].split()[-1].strip()),
                    energy_units="keV",
                )
                for _, row in grp.iterrows()
            ],
        )
        for (src, tele), grp in xray.groupby(["source", "telescope"])
    ]

    filter_alias += [
        dict(
            filter_key="0.3 - 10",
            filter_name="0.3 - 10",
            wave_eff=((10 - 0.3) * u.keV).to(u.nm, equivalencies=u.spectral()).value,
            wave_units="nm",
            wave_max=(0.3 * u.keV).to(u.nm, equivalencies=u.spectral()).value,
            wave_min=(10 * u.keV).to(u.nm, equivalencies=u.spectral()).value,
        ),
        dict(
            filter_key="0.2 - 12",
            filter_name="0.2 - 12",
            wave_eff=((12 - 0.2) * u.keV).to(u.nm, equivalencies=u.spectral()).value,
            wave_units="nm",
            wave_max=(0.2 * u.keV).to(u.nm, equivalencies=u.spectral()).value,
            wave_min=(12 * u.keV).to(u.nm, equivalencies=u.spectral()).value,
        ),
    ]

    # test the otter dataset format based on our schema
    for r in xray_phot:
        otter.schema.PhotometrySchema(**r)
    for f in filter_alias:
        otter.schema.FilterSchema(**f)

    # then merge everything
    otter_json = meta | {
        "photometry": [uvoir_phot] + radio_phot + xray_phot,
        "filter_alias": filter_alias,
    }

    otter.schema.OtterSchema(**otter_json)

    t = otter.Transient(otter_json)

    db = otter.Otter(datadir=args.otterdir)

    db.save([t], testing=False)


if __name__ == "__main__":
    main()
