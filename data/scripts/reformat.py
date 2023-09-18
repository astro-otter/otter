'''
Reformats some of the data in the json files
'''

import json
import os, glob

d = '/home/nfranz/tide/data/stash'
for jFile in glob.glob(os.path.join(d, '*.json')):
    with open(jFile, 'r') as f:
        j = json.load(f)

    dataDict = j[0]
    photo = dataDict['photometry']
    del dataDict['photometry']
    
    dataDict['photometry'] = {'radio': photo}

    outpath = os.path.join(os.getcwd(), os.path.basename(jFile))
    with open(outpath, 'w') as f:
        out = json.dumps(dataDict, indent=4)
        out += ']'
        out = '[' + out
        f.write(out)
    
    
