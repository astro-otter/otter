'''
Script to generate the summary table for OTTER searches

This needs to be run everytime this repo is downloaded right now because
the summary table stores absolute paths in it which is user dependent
'''

import argparse
from otter import Otter

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--otterdir', help='Directory where the otter json files will go')
    args = p.parse_args()

    # find the otter directory and then generate the summary table
    db = Otter(args.otterdir)
    db.generate_summary_table(save=True)

if __name__ == '__main__':
    main()
