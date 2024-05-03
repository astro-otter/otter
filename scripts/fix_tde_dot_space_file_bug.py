'''
tde.space has PS18kh.json and AT2018zr.json but they are the same event!
This script patches this bug instead of editing the source files, we want
those to stay the same for reproducability!
'''

import os, json
from otter import Otter

def main():

    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('--otterdir', help='Directory where the otter json files will go')
    args = p.parse_args()

    otterdir = args.otterdir

    otterdir = '/home/nfranz/astro-otter/otterdb/.otter'

    # the following json files have data for the same event...
    file1 = os.path.join(otterdir, 'PS18kh.json')
    file2 = os.path.join(otterdir, 'AT2018zr.json')
    
    # read in the jsons
    with open(file1, 'r') as f:
        jsontokeep = json.load(f)
        
    with open(file2, 'r') as f:
        otherjson = json.load(f)

    # now we can just merge these by hand 

    # we want the default name from the 2018zr one
    jsontokeep['name']['default_name'] = otherjson['name']['default_name']
    
    # add the alias information from that one now one
    jsontokeep['name']['alias'] += otherjson['name']['alias']
    
    # merge the coordinates
    jsontokeep['coordinate'] += otherjson['coordinate']
    
    # classifications are identical so we don't need to copy those over
    # copy over the dates
    jsontokeep['date_reference'] += otherjson['date_reference']
    
    # and then add in the reference alias, there isn't anything else in it!
    for ref in otherjson['reference_alias']:
        if ref not in jsontokeep['reference_alias']:
            jsontokeep['reference_alias'].append(ref)

    # now write out the new file
    out = json.dumps(jsontokeep, indent=4)
    with open(file2, 'w+') as f:
        f.write(out)

    # remove file1
    os.remove(file1)

    # we need to regenerate the summary.csv table too
    db = Otter(otterdir)
    db.generate_summary_table(save=True)
    
if __name__ == '__main__':
    main()
