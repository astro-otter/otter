"""
Download and convert the TNS data to OTTER format, then merge!
"""

# imports
import os
import zipfile
import requests
import json
import time
from copy import deepcopy
from collections import OrderedDict
import warnings

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord, search_around_sky
from astropy import units as u
from astropy.time import Time

from otter.io.otter import Otter, Transient
from otter.util import filter_to_obstype, FILTER_MAP_WAVE, bibcode_to_hrn


# helper functions
def set_bot_tns_marker():
    tns_marker = (
        'tns_marker{"tns_id": '
        + str(os.environ["TNS_BOT_ID"])
        + ', "type": "bot", "name": "'
        + str(os.environ["TNS_BOT_NAME"])
        + '"}'
    )
    return tns_marker


def format_to_json(source):
    parsed = json.loads(source, object_pairs_hook=OrderedDict)
    result = json.dumps(parsed, indent=4)
    return result


def is_string_json(string):
    try:
        json_object = json.loads(string)
    except Exception:
        return False
    return json_object


def search(
    ra: str,
    dec: str,
    radius: str = "10",
    radius_units: str = "arcsec",
    test: bool = False,
) -> dict:
    # get the url depending on if we are testing or not
    if test:
        url_tns_api = "https://sandbox.wis-tns.org/api/get"
    else:
        url_tns_api = "https://www.wis-tns.org/api/get"

    # setup the search object json format
    search_obj = [("ra", ra), ("dec", dec), ("radius", radius), ("units", radius_units)]

    # setup some important info for the search
    search_url = url_tns_api + "/search"
    tns_marker = set_bot_tns_marker()
    headers = {"user-agent": tns_marker}
    json_file = OrderedDict(search_obj)
    api_key = os.environ["TNS_API_KEY"]
    search_data = {"api_key": api_key, "data": json.dumps(json_file)}

    # actually perform the search
    print(headers, search_data)
    response = requests.post(search_url, headers=headers, data=search_data)
    return response


def get(objname: str, objid: int, phot: str = "1", spec: str = "0", test: bool = False):
    # get the url depending on if we are testing or not
    if test:
        url_tns_api = "https://sandbox.wis-tns.org/api/get"
    else:
        url_tns_api = "https://www.wis-tns.org/api/get"

    get_url = url_tns_api + "/object"
    tns_marker = set_bot_tns_marker()
    headers = {"User-Agent": tns_marker}

    get_obj = [
        ("objname", objname),
        ("objid", objid),
        ("photometry", phot),
        ("spectra", spec),
    ]
    json_file = OrderedDict(get_obj)

    get_data = {"api_key": os.environ["TNS_API_KEY"], "data": json.dumps(json_file)}
    response = requests.post(get_url, headers=headers, data=get_data)
    return response


def download_daily_csv(outfile="tns_public_objects.csv"):
    """
    Download the TNS daily CSV
    """

    # url to the daily staging of all TNS objects
    url = "https://www.wis-tns.org/system/files/tns_public_objects/"
    url += "tns_public_objects.csv.zip"

    # create the download header
    tns_marker = set_bot_tns_marker()
    headers = {"User-Agent": tns_marker}

    # put the api key in the data
    data = {"api_key": os.environ["TNS_API_KEY"]}

    # request the download
    response = requests.post(url, headers=headers, data=data)

    # write output to a zip file
    if response.status_code == 200:
        # this means it succeeded!
        outzip = outfile + ".zip"
        with open(outzip, "wb") as f:
            f.write(response.content)

    else:
        raise ValueError(
            f"Downloading the TNS Daily CSV Failed with Code {response.status_code}"
        )

    # unzip the outfile
    with zipfile.ZipFile(outzip, "r") as z:
        z.extractall(os.getcwd())

    # return the output file for ease later
    return outfile


def search_tns_daily_csv(db: Otter, tnsdaily: str, sep=5) -> pd.DataFrame:
    """
    Compare the TNS Daily CSV with an OTTER database and return a
    merged pandas dataframe of them
    """

    # read in the data
    tns_df = pd.read_csv(tnsdaily, skiprows=1)
    otter_df = db.generate_summary_table()

    # convert to skycoords
    tns_skycoords = SkyCoord(tns_df.ra, tns_df.declination, unit="deg")
    otter_skycoords = SkyCoord(otter_df.ra, otter_df.dec)

    # perform the crossmatching
    tns_idx, otter_idx, _, _ = search_around_sky(
        tns_skycoords, otter_skycoords, sep * u.arcsec
    )

    # now join the two dataframes on these indices
    tns_df_matches = tns_df.iloc[tns_idx].reset_index().add_suffix("_tns")
    otter_df_matches = otter_df.iloc[otter_idx].reset_index().add_suffix("_otter")
    alldata = pd.concat([tns_df_matches, otter_df_matches], axis=1)

    return alldata


# conversion functions
def tns_phot_to_otter_phot(photlist):
    data_otter_format = []
    filters_added_to_alias = []
    filter_alias = []
    badfilternames = {"Other", ""}

    for point in photlist:
        # first get flux and fluxerr
        flux, fluxerr = point["flux"], point["fluxerr"]
        if fluxerr == "":
            fluxerr = 0

        # now check if they are empty, if so we need to instead store
        # the limflux with upperlimit = True
        if isinstance(flux, str) and len(flux) == 0:
            # upperlimit!
            upperlimit = True
            flux, fluxerr = point["limflux"], 0
        else:
            upperlimit = False

        # now handle the flux units
        flux_unit = point["flux_unit"]["name"]
        if flux_unit == "ABMag":
            flux_unit = "mag(AB)"
        if flux_unit == "VegaMag":
            flux_unit = "vega"

        # fliters
        filterused = point["filters"]["name"]
        if filterused in badfilternames:
            continue  # we don't trust these

        # now handle the dates
        jd_date = Time(point["jd"], format="jd")
        date = jd_date.to_value("mjd")

        # we've extracted all the relevant info so we can convert to the OTTER format!
        # package all this info
        otter_format = dict(
            raw=[flux],
            raw_err=[fluxerr],
            upperlimit=upperlimit,
            raw_units=[flux_unit],
            date=str(date),  # always strings since some dates are formatted differently
            date_format="mjd",
            filter_key=filterused,
            reference=["TNS"],  # CORRECT THIS LATER? SEE GITHUB ISSUE
            obs_type=filter_to_obstype(filterused),
            # We assume that a sane observational astronomer would upload raw
            # but host subtracted photometry...
            corr_k=False,
            corr_s=False,
            corr_av=False,
            corr_host=True,
            corr_hostav=False,
        )

        # now some other optional info
        if "telescope" in point:
            otter_format["telescope"] = point["telescope"]["name"]
        if "instrument" in point:
            otter_format["instrument"] = point["instrument"]["name"]

        data_otter_format.append(otter_format)

        # now deal with the filter alias
        if filterused not in filters_added_to_alias:  # no double adding
            fa = dict(
                filter_key=filterused,
                filter_name=filterused,
                wave_eff=FILTER_MAP_WAVE[filterused],
                wave_units="nm",
            )
            filter_alias.append(fa)

        filters_added_to_alias.append(filterused)

    return data_otter_format, filter_alias


def tns_response_to_otter(row, photlist):
    """
    Convert all the TNS data we've pulled to the correct
    format for the OTTER schema

    Args:
        row [pandas.Series]: A row of the merged dataframe
        photlist [list[dict]]: A list of photometry from the TNS get API
    Returns:
        otter.Transient object in the correct format
    """

    out = deepcopy(Transient())

    # package the photometry
    out["photometry"], out["filter_alias"] = tns_phot_to_otter_phot(photlist)

    # get TNS coordinates and cite the discovery paper
    out["coordinate"] = [
        dict(
            ra=row.ra_tns,
            dec=row.declination_tns,
            epoch="J2000",
            frame="ICRS",
            coordinate_type="equitorial",
            ra_units="deg",
            dec_units="deg",
            reference=np.unique(
                [b.strip() for b in row.Discovery_ADS_bibcode_tns.split(",")]
            ).tolist(),
        )
    ]

    # Deal with TNS name
    out["name"] = dict(
        default_name=row.name_tns,
        alias=[
            dict(value=row.name_tns, reference=["TNS"]),
        ],
    )

    if not pd.isna(row.internal_names_tns):  # add more aliases
        alias_list = row.internal_names_tns.replace(",", "").split()
        for val in alias_list:
            out["name/alias"].append(dict(value=val, reference=["TNS"]))

    # add in discovery date
    if not (pd.isna(row.discoverydate_tns) and pd.isna(row.Discovery_ADS_bibcode_tns)):
        out["date_reference"] = [
            dict(
                value=Time(row.discoverydate_tns).mjd,
                date_format="mjd",
                date_type="discovery",
                reference=np.unique(
                    [b.strip() for b in row.Discovery_ADS_bibcode_tns.split(",")]
                ).tolist(),
                computed=False,
            )
        ]

    # add in the TNS redshift and classification
    # but only do this if there is a value for both redshift
    # and class ADS bibcde
    if not pd.isna(row.Class_ADS_bibcodes_tns):
        if not pd.isna(row.redshift_tns):
            out["distance"] = [
                dict(
                    value=row.redshift_tns,
                    computed=False,
                    distance_type="redshift",
                    reference=np.unique(
                        [b.strip() for b in row.Class_ADS_bibcodes_tns.split(",")]
                    ).tolist(),
                )
            ]

        if not pd.isna(row.type_tns):
            out["classification"] = [
                dict(
                    object_class=row.type_tns,
                    confidence=0.5,
                    reference=np.unique(
                        [b.strip() for b in row.Class_ADS_bibcodes_tns.split(",")]
                    ).tolist(),
                )
            ]

    # Deal with references
    out["reference_alias"] = [
        dict(
            name=row.Discovery_ADS_bibcode_tns,
            human_readable_name=bibcode_to_hrn(row.Discovery_ADS_bibcode_tns),
        )
    ]

    if not pd.isna(row.Class_ADS_bibcodes_tns):
        for bib in row.Class_ADS_bibcodes_tns.split(","):
            try:
                toadd = dict(name=bib, human_readable_name=bibcode_to_hrn(bib))
                out["reference_alias"].append(toadd)
            except ValueError:
                warnings.warn(f"Skipping the bibcode {bib}")
    return out


# now we can actually do all the data downloading and conversion
def main():
    import argparse

    p = argparse.ArgumentParser()
    p.add_argument("--otterdir", help="Directory where the otter json files will go")
    p.add_argument(
        "--sep",
        help="sky separation for comparing to otter in arcsec",
        default=1,
        type=float,
    )
    p.add_argument("--debug", action=argparse.BooleanOptionalAction)
    args = p.parse_args()

    # Create a connection with the database
    db = Otter(datadir=args.otterdir, gen_summary=True)

    # then download the TNS daily csv and crossmatch with OTTER
    # the crossmatching does all the photometry searching too
    outfile = download_daily_csv()
    tns_matches = search_tns_daily_csv(db, outfile, sep=args.sep)
    print(f"Found {len(tns_matches)} matches in TNS!")
    print()

    # merge the data!
    all_transients = []
    for idx, data in tns_matches.iterrows():
        print(f"Merging TNS:{data.name_tns} with OTTER:{data.name_otter}")

        # use that name and objid to grab the associated photometry
        response = get(data.name_tns, data.objid_tns)
        remaining = int(response.headers.get("x-rate-limit-remaining"))
        if remaining == 0:
            # this means we have to wait until the reset time
            # because we are out of calls to the TNS API
            time_to_reset = int(
                response.headers.get("x-rate-limit-reset")
            )  # in seconds
            print("###################################################################")
            print(f"Out of API Calls! Waiting {time_to_reset}s before the next one!")
            print("###################################################################")
            time.sleep(time_to_reset)

        outjson = response.json()
        outdata = outjson["data"]["reply"]
        phot = outdata["photometry"]

        # now clean up the photometry
        tns_transient = tns_response_to_otter(data, phot)
        all_transients.append(tns_transient)

    # now merge this data with the rest of otter and write it to a file
    print(all_transients)
    db.save(all_transients, testing=args.debug)


if __name__ == "__main__":
    main()
