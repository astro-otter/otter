"""
Converted the "Curated Optical TDE Catalog" from Mummery et al. (2023) and
provided to us by Sjoert van Velzen via https://github.com/sjoertvv/manyTDE
"""

import os
import glob
import json
import argparse
import numpy as np

from otter import Otter, Transient, util


def fix_old_bibcodes(ref):
    # build in some temporary checks until the catalog is fixed
    if ref == "Wevers 2020MNRAS.497L...1W":
        ref = ref.replace("Wevers ", "")

    if ref == "2018A%26A...610A..14K":
        ref = "2018A&A...610A..14K"

    if ref == "2024arXiv240111773C":
        ref = "2024A&A...689A.350C"

    if ref == "SDSS":
        # this is just for now until I figure out what this means in more detail
        # this reference is to the SDSS DR18 paper which briefly describes
        # updates to their velocity dispersion calculations
        ref = "2023ApJS..267...44A"

    if ref == "arXiv:2405.11343":
        ref = "2024ApJ...976...34Y"

    return ref


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--otterdir", help="Directory where the otter json files will go")
    p.add_argument("--indir", help="input directory for the raw data")
    p.add_argument("--debug", action=argparse.BooleanOptionalAction)
    args = p.parse_args()

    # define some variables from the inputs
    datapath = args.indir
    otterpath = args.otterdir

    # read in the OTTER database
    db = Otter(datadir=otterpath, gen_summary=True)

    # perform the conversion
    injsons = glob.glob(os.path.join(datapath, "*.json"))
    catalog_bibcode = "2024MNRAS.527.2452M"
    catalog_hrn = "Mummery et al. (2023)"

    new_transients = []
    for path in injsons:
        # print(path)
        with open(path, "r") as f:
            indata = json.load(f)

        otterjson = {}

        # handle references first
        bibcodes = [fix_old_bibcodes(b) for b in indata["paper_ref"].split(",")]
        hrn, correct_bibcodes = util.bibcode_to_hrn(bibcodes)

        otterjson["reference_alias"] = [
            dict(name=b, human_readable_name=h) for b, h in zip(correct_bibcodes, hrn)
        ]

        # add a citation to Mummery et al. (2023)
        # since that's where this catalog is from
        otterjson["reference_alias"] += [
            dict(name=catalog_bibcode, human_readable_name=catalog_hrn)
        ]

        # now do the names
        otterjson["name"] = dict(
            default_name=indata["name"],
            alias=[
                dict(
                    value=indata["name"],
                    reference=indata["paper_ref"].split(",") + [catalog_bibcode],
                )
            ],
        )

        # and now the classification
        otterjson["classification"] = dict(
            value=[
                dict(
                    object_class="TDE",
                    confidence=3.3,  # these are all spectroscopically confirmed
                    reference=indata["paper_ref"].split(",") + [catalog_bibcode],
                    default=True,
                )
            ]
        )

        # now the coordinate
        otterjson["coordinate"] = [
            dict(
                ra=indata["ra"],
                dec=indata["dec"],
                ra_units="deg",
                dec_units="deg",
                reference=indata["paper_ref"].split(",") + [catalog_bibcode],
                coordinate_type="equatorial",
            )
        ]

        # now redshift
        otterjson["distance"] = [
            dict(
                value=indata["z"],
                reference=indata["paper_ref"].split(",") + [catalog_bibcode],
                distance_type="redshift",
            )
        ]

        # now the velocity dispersion
        if "vel_disp_km/s" in indata["host"]:
            ref = indata["host"]["vel_disp_source"]
            ref = fix_old_bibcodes(ref)

            _, hrn = util.bibcode_to_hrn(ref)

            curr_bibcodes = {x["name"] for x in otterjson["reference_alias"]}
            if ref not in curr_bibcodes:
                otterjson["reference_alias"].append(
                    dict(name=ref, human_readable_name=hrn[0])
                )

            otterjson["distance"].append(
                dict(
                    value=indata["host"]["vel_disp_km/s"],
                    error=indata["host"]["e_vel_disp_km/s"],
                    reference=[ref, catalog_bibcode],
                    unit="km/s",
                    distance_type="dispersion_measure",
                )
            )

        # now peak date
        otterjson["date_reference"] = [
            dict(
                value=indata["peak_mjd"],
                date_format="mjd",
                reference=indata["paper_ref"].split(",") + [catalog_bibcode],
                computed=False,
                date_type="peak",
            )
        ]

        # now the photometry
        lc = indata["lightcurve"]
        mjd, filter, flux, fluxerr = np.array(lc["data"]).T
        # filtermap = {filt:freq for filt,freq in zip(lc['filters'],lc['frequency_Hz'])}

        flux = flux.astype(float)
        fluxerr = fluxerr.astype(float)

        # clean out NaN values from the flux/fluxerr arrays
        notna = np.argwhere(~np.isnan(flux))
        flux = flux[notna]
        fluxerr = fluxerr[notna]
        mjd = mjd[notna]
        filter = filter[notna]

        # package the data into the photometry array
        otterjson["photometry"] = [
            dict(
                raw=list(flux),
                raw_err=list(fluxerr),
                raw_units="Jy",
                date=list(mjd),
                date_format="mjd",
                filter_key=list(filter),
                obs_type="uvoir",
                reference=indata["paper_ref"].split(",") + [catalog_bibcode],
                upperlimit=list(flux < 3 * fluxerr),
                corr_k=False,
                corr_s=False,
                corr_av=False,
                corr_host=False,
                corr_hostav=False,
                # 3.2 factor conversion value from E(B-V) to A_V comes from
                # https://astronomy.swin.edu.au/cosmos/i/interstellar+reddening
                val_av=3.2 * indata["extinction"]["e_bv"],
            )
        ]

        otterjson["filter_alias"] = [
            dict(filter_key=filt, filter_name=filt, freq_eff=freq, freq_units="Hz")
            for filt, freq in zip(lc["filters"], lc["frequency_Hz"])
        ]

        # print(filtermap)

        # print(indata.keys())
        # print(json.dumps(otterjson, indent=4))

        newt = Transient(otterjson)
        new_transients.append(newt)

    print()
    print(f"Found {len(new_transients)} TDEs that are usable!")

    db.save(new_transients, testing=args.debug)


if __name__ == "__main__":
    main()
