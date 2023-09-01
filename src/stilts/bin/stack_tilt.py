import argparse
import configparser
import sys

import tipsy

# Steps
# 1. Take config file from command line and list of times from command line or stdin
# 1. Ingest config file with seismometer info and data downlad info
# 2. Download data from IRIS
# 3. Run tilt extraction
# 4. Make stacked plot, and individual event plots of tilt and seismograms


def cli():
    # Command line interface
    parser = argparse.ArgumentParser(
        description="Generate seismometer derived tilt plots, stacking multiple events"
    )
    parser.add_argument(
        "config", type=str, help="Configuration file with sensor and data download info"
    )
    parser.add_argument("-t", "--times", type=str, help="File with list of event times")
    args = parser.parse_args()

    return args


def readConfig(configFile):
    # Ingest config file
    cfg = configparser.ConfigParser()
    cfg.read(configFile)
    return cfg


def main():
    args = cli()
    cfg = readConfig(args.config)
    
    # Read times from file if provided, otherwise use stdin
    if args.times is None:
        times = sys.stdin.readlines()
    else:
        with open(args.times, "r") as f:
            times = f.readlines()
            
    times = tipsy.parseTimes(times)

    streams, times = tipsy.downloadStreams(times, cfg["source"])

    if len(streams) == 0:
        print("No data found at any times.", file=sys.stderr)
        return 1

    # Fix one-extra-sample in some records
    # This seems a bit awkward
    minLen = len(streams[0][0].data)  # initial
    for stream in streams:
        for trace in stream:
            if len(trace.data) < minLen:
                if abs(len(trace.data) - minLen) > 1:
                    print(
                        "Trace length differs by more than 1, exiting", file=sys.stderr
                    )
                    sys.exit(1)
                minLen = len(trace.data)

    for stream in streams:
        for trace in stream:
            trace.data = trace.data[:minLen]

    tilts = tipsy.calcTilt(streams, cfg["sensor"])
    tipsy.genFigure(streams, tilts, times, cfg["source"])
    tipsy.saveTilts(tilts, cfg["source"])
    return 0


if __name__ == "__main__":
    main()
