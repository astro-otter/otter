"""
Keep the local JSON file with references up to date
"""

import ads
import otter
import json


def query_ads(bibcodes):
    query = f"bibcode:{bibcodes[0]}"
    if len(bibcodes) > 1:
        for b in bibcodes[1:]:
            query += f" OR {b}"

    qobj = ads.SearchQuery(q=query)
    qobj.execute()  # do the query
    adsquery = list(qobj)

    hrns = []
    bibcodes_to_return = []
    for res in adsquery:
        authors = res.author
        year = res.year

        if len(authors) == 0:
            raise ValueError("This ADS bibcode does not exist!")
        elif len(authors) == 1:
            author = authors[0]
        elif len(authors) == 2:
            author = authors[0] + " & " + authors[1]
        else:  # longer than 2
            author = authors[0] + " et al."

        # generate the human readable name
        hrn = author + " (" + year + ")"
        hrns.append(hrn)
        bibcodes_to_return.append(res.bibcode)

    missing = list(set(bibcodes) - set(bibcodes_to_return))

    return bibcodes_to_return, hrns, missing


def main():
    db = otter.Otter()

    all_transients = db.get_meta()

    bibcodes = []
    for t in all_transients:
        for ref in t["reference_alias"]:
            if ref["name"].strip() not in bibcodes:
                bibcodes.append(ref["name"].strip())

    dt = 10
    cursor = 0
    bibcode_map_local = {}
    tobreak = False
    all_missing = []
    while True:
        print(cursor, "--->", cursor + dt)

        if cursor + dt > len(bibcodes):
            bibs, hrns, missing = query_ads(bibcodes[cursor:])
            tobreak = True
        else:
            bibs, hrns, missing = query_ads(bibcodes[cursor : cursor + dt])

        for b, h in zip(bibs, hrns):
            bibcode_map_local[b] = h

        cursor += dt

        if tobreak:
            break

        all_missing += missing

    old_bibcode, successful_bibcodes, successful_hrns = [], [], []
    toskip = {"ASAS-SN Supernovae", "MAST", "TNS"}

    for missing in all_missing:
        if missing in toskip:
            old_bibcode.append(missing)
            successful_bibcodes.append(missing)
            successful_hrns.append(missing)
            continue

        query = f"bibcode:{missing} OR {missing}"
        print(missing)
        qobj = ads.SearchQuery(q=query)
        qobj.execute()  # do the query
        adsquery = list(qobj)
        print(adsquery)
        res = adsquery[0]
        authors = res.author
        year = res.year

        if len(authors) == 0:
            raise ValueError("This ADS bibcode does not exist!")
        elif len(authors) == 1:
            author = authors[0]
        elif len(authors) == 2:
            author = authors[0] + " & " + authors[1]
        else:  # longer than 2
            author = authors[0] + " et al."

        # generate the human readable name
        hrn = author + " (" + year + ")"
        bibcode_map_local[res.bibcode] = hrn
        bibcode_map_local[missing] = hrn

    with open("reference_map_local.json", "w+") as f:
        json.dump(bibcode_map_local, f, indent=4)


if __name__ == "__main__":
    main()
