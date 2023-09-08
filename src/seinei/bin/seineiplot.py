import argparse
import sys
import os.path

import obspy
import numpy as np
import matplotlib.pyplot as plt

def cli():
    # Command line interface
    parser = argparse.ArgumentParser(
        description="""Generate plots of SAC files"""
    )
    parser.add_argument(
        "sac",
        nargs="+",
        default="",
        help="SAC data file(s) to plot",
    )
    args = parser.parse_args()

    return args


def main():
    args = cli()

    # Read SAC data files
    sts = {}
    for file in args.sac:
        sts[file] = obspy.read(file)

    # Loop over streams (corresponding to files)
    for name, st in sts.items():
        st.plot(outfile=name.replace(".sac",".png"))
        

    return 0


if __name__ == "__main__":
    main()
