"""
Converts the data from tde.space to the otter json format while also incorporating it
with existing files
"""

# imports
import os
import glob
from copy import deepcopy
import numpy as np
import pandas as pd
import json
from otter import Otter, Transient
from otter import util as otter_const
from otter import util as otter_helper
import astropy.units as u
import astropy.constants as const

# some useful mappings
bandwavelengths = {
    "u": 354.0,
    "g": 475.0,
    "r": 622.0,
    "i": 763.0,
    "z": 905.0,
    "u'": 354.0,
    "g'": 475.0,
    "r'": 622.0,
    "i'": 763.0,
    "z'": 905.0,
    "u_SDSS": 354.3,
    "g_SDSS": 477.0,
    "r_SDSS": 623.1,
    "i_SDSS": 762.5,
    "z_SDSS": 913.4,
    "U": 365.0,
    "B": 445.0,
    "V": 551.0,
    "R": 658.0,
    "I": 806.0,
    "Y": 1020.0,
    "J": 1220.0,
    "H": 1630.0,
    "K": 2190.0,
    "M2": 260.0,
    "W1": 224.6,
    "W2": 192.8,
    "w": 622.0,
}

bandreps = {
    "Ks": ["K_s"],
    "M2": ["uvm2", "UVM2", "UVm2", "Um2", "m2", "um2"],
    "W1": ["uvw1", "UVW1", "UVw1", "Uw1", "w1", "uw1"],
    "W2": ["uvw2", "UVW2", "UVw2", "Uw2", "w2", "uw2"],
}


# we need the minimum, maximum, and effective wavelengths
def xray_to_wave(xraycode):
    min_energy, max_energy = [float(val) for val in xraycode.split(" - ")]
    wave_eff = (const.h * const.c / (max_energy - min_energy) / u.keV).to(u.nm).value

    # these variable names seem backwards but remember that energy and
    # wavelength are inversely proportional!
    wave_min = (const.h * const.c / max_energy / u.keV).to(u.nm).value
    wave_max = (const.h * const.c / min_energy / u.keV).to(u.nm).value

    return (wave_min, wave_eff, wave_max)


xraycodes = {
    "0.3 - 10": xray_to_wave("0.3 - 10"),
    "0.5 - 8": xray_to_wave("0.5 - 8"),
    "0.3 - 2.0": xray_to_wave("0.3 - 2.0"),
    "0.2 - 2.0": xray_to_wave("0.2 - 2.0"),
}

# from https://imagine.gsfc.nasa.gov/science/toolbox/spectrum_chart.html
bandoptions = {
    "radio": {"min": 1e-3, "max": 1e-100},
    "x-ray": {"min": 1e-11, "max": 1e-8},
    "optical": {"min": 4e-7, "max": 7e-7},
    "infared": {"min": 7e-7, "max": 1e-3},
    "uv": {"min": 1e-8, "max": 4e-7},
}

# now the functions


def clean_schema(schema):
    """
    Clean out Nones and empty lists from the given subschema
    """
    protected_keys = {
        "reference",
        "corr_k",
        "corr_s",
        "corr_av",
        "corr_host",
        "corr_hostav",
    }
    for key, val in list(schema.items()):
        if val is None or (isinstance(val, (list, dict)) and len(val) == 0):
            if key not in protected_keys:
                del schema[key]
    return schema


def main():
    import argparse

    p = argparse.ArgumentParser()
    p.add_argument("--indir", help="Directory for the basic radio data")
    p.add_argument("--outdir", help="Directory where the otter json files will go")
    args = p.parse_args()

    if not os.path.exists(args.outdir):
        os.mkdir(args.outdir)

    exc = {"ASAS-SN Supernovae", "MAST", "Swift TOO"}

    def mappedsrc(source_map, src):
        out = []
        for key in src:
            if key in source_map:
                toappend = source_map[key]
                out.append(toappend)
        return out

    allschemas = []
    allsrcs = 0
    badsrcs = 0
    allphot = 0
    badphot = 0

    for file in glob.glob(os.path.join(args.indir, "*.json")):
        # theres this one bugged file from tde.space called AT.json
        # that has data from multiple transients in it
        # so we'll just skip that one
        if os.path.basename(file) == "AT.json":
            continue

        print(f"Reformating {file}")
        with open(file, "r") as f:
            j = json.load(f)[0]

        # copy over the references
        schema = Transient(deepcopy(otter_const.schema))
        source_map = {}

        for src in j["sources"]:
            allsrcs += 1  # another source, double are okay I think...
            if "bibcode" not in src and src["name"] in exc:
                bib = src["name"]
            elif "bibcode" not in src:
                print(f"No bibcode for {src}, skipping!")
                badsrcs += 1
                continue
            else:
                bib = src["bibcode"]
            sub = deepcopy(otter_const.subschema["reference_alias"])

            if bib == "Transient Name Server":
                sub["name"] = "TNS"
            else:
                sub["name"] = bib

            if "reference" in src:
                sub["human_readable_name"] = src["reference"]
            else:
                sub["human_readable_name"] = src["name"]
            schema["reference_alias"].append(sub)

            source_map[src["alias"]] = bib

        # copy over names and aliases
        schema["name/default_name"] = j["name"]
        for val in j["alias"]:
            sub = deepcopy(otter_const.subschema["name/alias"])
            sub["value"] = val["value"]
            ref = mappedsrc(source_map, val["source"].split(","))
            if "Transient Name Server" in val["source"].split(","):
                ref.append("TNS")
            sub["reference"] = ref
            schema["name/alias"].append(sub)

        # copy over coordinates
        if "ra" in j and "dec" in j:
            key1, key2 = "ra", "dec"
        elif "hostra" in j and "hostdec" in j:
            # see the Table 1 note for why this is okay to do in this paper
            # https://iopscience.iop.org/article/10.3847/1538-4357/aa633b
            key1, key2 = "hostra", "hostdec"
        else:
            print(f"Skipping {file} because no ra and dec associated with it!")
            continue

        for ra, dec in zip(j[key1], j[key2]):
            sub = deepcopy(otter_const.subschema["coordinate"])
            if ra["source"] != dec["source"]:
                src = list(set(ra["source"].split(",")) or set(ra["source"].split(",")))
            else:
                src = ra["source"].split(",")

            if not any(s in source_map for s in src):
                print(
                    f"Skipping {(ra, dec)} because it does not"
                    + "have a reliable reference!"
                )
                continue

        if ra["u_value"] == "hours" or ":" in ra["value"]:
            ra_u = "hour"
        elif ra["u_value"] == "degrees" or ra["u_value"] == "floatdegrees":
            ra_u = "deg"
        else:
            ra_u = ra["u_value"]

        if dec["u_value"] == "degrees" or dec["u_value"] == "floatdegrees":
            dec_u = "deg"
        else:
            dec_u = dec["u_value"]

        sub["ra"] = ra["value"]
        sub["dec"] = dec["value"]
        sub["ra_units"] = ra_u
        sub["dec_units"] = dec_u
        sub["reference"] = mappedsrc(source_map, src)
        sub["coordinate_type"] = "equitorial"
        schema["coordinate"].append(clean_schema(sub))

        # copy over distance measurements
        # first redshift
        if "redshift" in j or "lumdist" in j or "comovingdist" in j:
            if "redshift" in j:
                for ii, z in enumerate(j["redshift"]):
                    sub = deepcopy(otter_const.subschema["distance"])
                    src = z["source"].split(",")
                    if not any(s in source_map for s in src):
                        print(
                            f"Skipping {z} because it does not have a "
                            + "reliable reference!"
                        )
                        continue

                    sub["value"] = z["value"]
                    sub["reference"] = mappedsrc(source_map, src)
                    sub["computed"] = False
                    sub["distance_type"] = "redshift"
                    if ii == 0:
                        sub["default"] = True
                    schema["distance"].append(clean_schema(sub))
                del j["redshift"]

            if "lumdist" in j:
                for ii, d in enumerate(j["lumdist"]):
                    sub = deepcopy(otter_const.subschema["distance"])
                    src = d["source"].split(",")
                    if not any(s in source_map for s in src):
                        print(
                            f"Skipping {d} because it does not have a "
                            + "reliable reference!"
                        )
                        continue

                    sub["value"] = d["value"]
                    sub["reference"] = mappedsrc(source_map, src)
                    sub["computed"] = False
                    sub["unit"] = d["u_value"]
                    sub["distance_type"] = "luminosity"
                    if ii == 0:
                        sub["default"] = True
                    schema["distance"].append(clean_schema(sub))
                del j["lumdist"]

            if "comovingdist" in j:
                for ii, d in enumerate(j["comovingdist"]):
                    sub = deepcopy(otter_const.subschema["distance"])
                    src = d["source"].split(",")
                    if not any(s in source_map for s in src):
                        print(
                            f"Skipping {d} because it does not have a "
                            + "reliable reference!"
                        )
                        continue

                    sub["value"] = d["value"]
                    sub["reference"] = mappedsrc(source_map, src)
                    sub["computed"] = False
                    sub["unit"] = d["u_value"]
                    sub["distance_type"] = "comoving"
                    if ii == 0:
                        sub["default"] = True
                    schema["distance"].append(clean_schema(sub))
                del j["comovingdist"]

            if len(schema["distance"]) == 0:
                del schema["distance"]

        # dates
        if "discoverdate" in j:
            for ii, d in enumerate(j["discoverdate"]):
                sub = deepcopy(otter_const.subschema["date_reference"])
                src = d["source"].split(",")
                if not any(s in source_map for s in src):
                    print(
                        f"Skipping {d} because it does not have a reliable reference!"
                    )
                    continue

                if "/" in d["value"]:
                    indate = d["value"].replace("/", "-")
                    sub["date_format"] = "iso"
                    day = d["value"].split("/")[-1]
                    if "." in day:
                        # shorten the day so it doesn't have a float
                        indate = indate.replace(day, day.split(".")[0])

                        # then this has a time of day associated with it
                        intime = float(day) % 1  # just get the decimal point
                        hrs = int(np.floor(intime * 24))
                        min_leftover = (intime * 24) % 1
                        minutes = int(np.floor(min_leftover * 60))
                        sec_leftover = (min_leftover * 60) % 1
                        sec = sec_leftover * 60  # astropy allows this to be a float!

                        indate += f" {hrs}:{minutes}:{sec}"

                elif (
                    isinstance(d["value"], str) and len(d["value"]) == 4
                ):  # this is just a year
                    indate = f'{d["value"]}-01-01'
                    sub["date_format"] = "iso"
                else:
                    print(j["discoverdate"])
                    raise ValueError("New Type of date found, please fix!")

                sub["value"] = indate
                sub["reference"] = mappedsrc(source_map, src)

                if "derived" in d:
                    sub["computed"] = d["derived"]
                else:
                    sub["computed"] = False

                sub["computed"] = False
                if ii == 0:
                    sub["default"] = True
                sub["date_type"] = "discovery"
                schema["date_reference"].append(clean_schema(sub))
            del j["discoverdate"]

        # copy classification
        if "claimedtype" in j:
            # print(json.dumps(j, indent=4))
            for ii, d in enumerate(j["claimedtype"]):
                src = d["source"].split(",")
                if not any(s in source_map for s in src):
                    print(
                        f"Skipping {d} because it does not have a reliable reference!"
                    )
                    continue

                # c = d["value"]
                conf = 0  # It is unclear where this classification came from

                sub = deepcopy(otter_const.subschema["classification"])
                sub["object_class"] = "TDE"
                sub["reference"] = mappedsrc(source_map, src)
                sub["confidence"] = conf
                schema["classification"]["value"].append(clean_schema(sub))
            del j["claimedtype"]

        # copy photometry
        if "photometry" in j:
            phot = pd.DataFrame(j["photometry"])
            allphot += len(phot)

            if "telescope" in phot and "u_time" in phot:
                gby = ["source", "telescope", "u_time"]
            elif "u_time" in phot:
                gby = ["source", "u_time"]
            else:
                gby = ["source"]

            for grouped_by, group in phot.groupby(gby):
                if len(grouped_by) == 3:
                    ref, telescope, tfmt = grouped_by
                elif len(grouped_by) == 2:
                    ref, tfmt = grouped_by
                    telescope = None
                else:
                    ref = grouped_by
                    telescope = None
                    raise ValueError(
                        "Time format uncertain"
                    )  # for now, hopefully we never come to this anyways...

                # clean up group
                group = group.dropna(axis=1, how="any")
                # print(group)
                sub = deepcopy(otter_const.subschema["photometry"])
                src = ref.split(",")
                if not any(s in source_map for s in src):
                    print(
                        f"Skipping {d} because it does not have a reliable reference!"
                    )
                    badphot += len(group)
                    continue

                sub["reference"] = mappedsrc(source_map, src)
                sub["telescope"] = telescope

                usedfluxforraw = False
                if "magnitude" in group:
                    sub["raw"] = list(group.magnitude.astype(float))
                    if "e_magnitude" in group:
                        sub["raw_err"] = list(group.e_magnitude.astype(float))
                    if "system" in group:
                        sub["raw_units"] = list(group.system)
                    else:
                        sub["raw_units"] = (
                            "mag(AB)"  # hopefully this assumption is okay!
                        )
                    if "band" in group:
                        sub["filter_key"] = list(group.band)
                    else:
                        print(
                            f"Skipping {group} because no filter given so unreliable!"
                        )
                        badphot += len(group)
                        continue  # don't add this one
                elif "countrate" in group:
                    sub["raw"] = list(group.countrate.astype(float))
                    if "e_countrate" in group:
                        sub["raw_err"] = list(group.e_countrate.astype(float))
                    if "u_countrate" in group:
                        raw_units = []
                        for unit in group.u_countrate:
                            if unit == "s^-1":
                                raw_units.append("ct")
                            else:
                                raw_units.append(unit)

                        sub["raw_units"] = raw_units
                    else:
                        print(f"Skipping {group} because no units given so unreliable!")
                        badphot += len(group)
                        continue  # don't add this one

                    if "energy" in group:
                        sub["filter_key"] = [
                            " - ".join(e) if isinstance(e, list) else e
                            for e in group.energy
                        ]
                    else:
                        raise ValueError()
                        print(
                            f"Skipping {group} because no filter given so unreliable!"
                        )
                        continue  # don't add this one

                    # add the xray model, if information on it is given
                    sub["xray_model"] = []

                    has_abs_pl = "photonindex" in group and "nhmw" in group
                    has_pl = "photonindex" in group and "nhmw" not in group
                    for _, row in group.iterrows():
                        if (
                            has_abs_pl
                            and not pd.isna(row.photonindex)
                            and not pd.isna(row.nhmw)
                        ):
                            sub["xray_model"].append(
                                dict(
                                    model_name="Absorbed Powerlaw",
                                    param_names=["photonindex", "nhmw"],
                                    param_values=[row.photonindex, row.nhmw],
                                    param_units=["None", "10^22 cm^-2"],
                                    min_energy=float(row.energy[0]),
                                    max_energy=float(row.energy[1]),
                                    energy_units=row.u_energy,
                                )
                            )
                        elif (
                            has_abs_pl
                            and not pd.isna(row.photonindex)
                            and pd.isna(row.nhmw)
                        ):
                            sub["xray_model"].append(
                                dict(
                                    model_name="Powerlaw",
                                    param_names=["photonindex"],
                                    param_values=[row.photonindex],
                                    param_units=["None"],
                                    min_energy=float(row.energy[0]),
                                    max_energy=float(row.energy[1]),
                                    energy_units=row.u_energy,
                                )
                            )

                        elif has_pl and not pd.isna(row.photonindex):
                            sub["xray_model"].append(
                                dict(
                                    model_name="Powerlaw",
                                    param_names=["photonindex"],
                                    param_values=[row.photonindex],
                                    param_units=["None"],
                                    min_energy=float(row.energy[0]),
                                    max_energy=float(row.energy[1]),
                                    energy_units=row.u_energy,
                                )
                            )

                        else:
                            # the x-ray model is unknown :(
                            sub["xray_model"].append("Unknown X-ray Model!")

                elif (
                    "flux" in group and "u_flux" in group
                ):  # we will just use the flux for the raw data
                    usedfluxforraw = True
                    sub["raw"] = list(group.flux.astype(float))
                    sub["raw_units"] = list(group.u_flux)
                    if "e_flux" in group:
                        sub["raw_err"] = list(group.e_flux.astype(float))
                    if "band" in group:
                        sub["filter_key"] = list(group.band)
                    elif "energy" in group:
                        sub["filter_key"] = [
                            " - ".join(e) if isinstance(e, list) else e
                            for e in group.energy
                        ]
                    else:
                        print(
                            f"Skipping {group} because no filter given so unreliable!"
                        )
                        raise ValueError()
                        continue  # don't add this one

                elif (
                    "fluxdensity" in group and "u_fluxdensity" in group
                ):  # we will just use the flux for the raw data
                    usedfluxforraw = True
                    sub["raw"] = list(group.fluxdensity.astype(float))
                    sub["raw_units"] = list(group.u_fluxdensity)
                    if "e_fluxdensity" in group:
                        sub["raw_err"] = list(group.e_fluxdensity.astype(float))
                    if "frequency" in group and "u_frequency" in group:
                        sub["filter_key"] = [
                            freq + freq_unit
                            for freq, freq_unit in zip(
                                group.frequency, group.u_frequency
                            )
                        ]
                    else:
                        print(f"Skipping {group} because no frequency,  so unreliable!")
                        raise ValueError(
                            "Something didn't work with the tde.space radio data"
                        )
                        continue  # don't add this one

                else:
                    print(
                        "Skipping this photometry point because "
                        + "it is an unknown type! Please fix!"
                    )
                    badphot += len(group)

                if "flux" in group and "u_flux" in group and not usedfluxforraw:
                    sub["value"] = list(group.flux.astype(float))
                    sub["value_units"] = list(group.u_flux)
                    if "e_flux" in group:
                        sub["value_err"] = list(group.e_flux.astype(float))

                if "upperlimit" in group:
                    sub["upperlimit"] = list(group["upperlimit"])
                else:
                    sub["upperlimit"] = [False] * len(
                        group
                    )  # james tended to only store this if it is True

                sub["date"] = [
                    t[0] if isinstance(t, list) else t for t in group["time"]
                ]
                sub["date_format"] = tfmt

                # add the observation type

                if sub["filter_key"] is None:
                    continue
                sub["obs_type"] = []
                for filt in sub["filter_key"]:
                    if filt in bandwavelengths:
                        sub["obs_type"].append("uvoir")
                    elif filt in xraycodes:
                        sub["obs_type"].append("xray")
                    else:
                        sub["obs_type"].append(otter_helper.filter_to_obstype(filt))

                # add info about the corrections applied to the photometry
                # we don't have any right now so for now we will just put None
                sub["corr_k"] = None
                sub["corr_s"] = None
                sub["corr_av"] = None
                sub["corr_host"] = None
                sub["corr_hostav"] = None

                schema["photometry"].append(clean_schema(sub))

            del j["photometry"]

        # create the filter alias
        schema["filter_alias"] = []
        for val in schema["photometry"]:
            filter_keys = np.unique(val["filter_key"])
            for key in filter_keys:
                curr_filts = [filt["filter_key"] for filt in schema["filter_alias"]]
                if key in curr_filts:
                    continue

                sub = deepcopy(otter_const.subschema["filter_alias"])
                sub["filter_key"] = key
                sub["filter_name"] = key
                if key in bandwavelengths:
                    sub["wave_eff"] = bandwavelengths[key]
                    sub["wave_units"] = "nm"

                elif key in otter_const.FILTER_MAP_WAVE:
                    sub["wave_units"] = "nm"
                    sub["wave_eff"] = otter_const.FILTER_MAP_WAVE[key]

                elif key in xraycodes:
                    sub["wave_min"], sub["wave_eff"], sub["wave_max"] = xraycodes[key]
                    sub["wave_units"] = "nm"

                elif "Hz" in key:
                    freq = float("".join(filter(str.isdigit, key)))
                    sub["freq_eff"] = freq
                    sub["filter_name"] = otter_helper.freq_to_band(freq)
                    sub["freq_units"] = key.replace(str(freq), "")

                else:
                    raise ValueError(
                        "Can not add filter {key} because we dont know wave_eff"
                    )

                for key, val in list(sub.items()):
                    if val is None:
                        del sub[key]

                schema["filter_alias"].append(sub)

        # ADD HOST INFO ONCE WE HAVE A BETTER IDEA OF FORMATTING

        # remove everything we've added from j to help keep track
        del j["name"]
        del j["alias"]
        del j["sources"]
        del j[key1]
        del j[key2]

        # print(j.keys())
        # print(json.dumps(j, indent=4))

        # some checks of the outputs
        assert len(schema["coordinate"]) > 0
        assert len(schema["name/alias"]) > 0

        # get rid of dispersion measure cause none of these have it
        if len(schema["distance"]) == 0:
            del schema["distance"]

        del schema["spectra"]

        if len(schema["classification"]) == 0:
            # give it a low classification score
            schema["classification"] = {
                "value": [
                    {
                        "object_class": "TDE",
                        "confidence": 0,  # we don't trust this classification very much
                        "reference": ["2017ApJ...835...64G"],
                        "default": True,
                    }
                ]
            }

        if len(schema["photometry"]) == 0:
            del schema["photometry"]
            del schema["filter_alias"]  # this will be empty too

        if "date_reference" in schema and len(schema["date_reference"]) == 0:
            del schema["date_reference"]

        if "distance" in schema and len(schema["distance"]) == 0:
            del schema["distance"]

        allschemas.append(schema)

        # json_schema = json.dumps(dict(schema), indent=4) # clean_schema(dict(schema))
        print(schema)
        print()

    # print some useful stats
    print()
    print("#######################################################################")
    print(
        f"Skipped {badsrcs/allsrcs * 100 : .2f}% "
        + "of sources because we did not recognize them"
    )
    print(
        f"Skipped {badphot/allphot * 100 : .2f}% "
        + "of photometry points because we did not trust them"
    )
    print("#######################################################################")
    print()

    # add the data to the output directory
    p = args.outdir
    otter = Otter(datadir=p, gen_summary=True)

    otter.save(allschemas, testing=False)


if __name__ == "__main__":
    main()
