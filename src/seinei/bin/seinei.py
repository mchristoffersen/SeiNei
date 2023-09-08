import argparse
import sys
import os.path

import obspy
import numpy as np
import matplotlib.pyplot as plt

def cli():
    # Command line interface
    parser = argparse.ArgumentParser(
        description="""Calculate long-period tilt from broadband seismometer
        data. Input data can come from local SAC files or via IRIS download."""
    )
    parser.add_argument(
        "inventory",
        type=str,
        help="StationXML inventory file with station response information",
    )
    parser.add_argument(
        "sac",
        nargs="+",
        default="",
        help="SAC data file(s) with seismometer data in ADC counts",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="Directory to save output files (default=same as SAC files)",
    )
    parser.add_argument(
        "-pbl",
        "--passband_long_period",
        default=None,
        help="Long-period edge of passband (default=auto-set)",
    )
    parser.add_argument(
        "-pbs",
        "--passband_short_period",
        default=None,
        help="Short-period edge of passband (default=auto-set)",
    )
    parser.add_argument(
        "-g",
        "--gravity",
        default=9.81,
        help="Gravitaitonal acceleration in m/s^2 (default=9.81)",
    )
    parser.add_argument('-v', '--verbose', action='store_true', help="Enable verbose output") 
    args = parser.parse_args()

    return args


def main():
    args = cli()

    # Read StationXML inventory file, complain if obspy doesn't like it
    try:
        if(args.verbose):
            print("Reading inventory file:\n\t%s" % args.inventory)
        inv = obspy.core.inventory.inventory.read_inventory(
            args.inventory, format="STATIONXML", level="response"
        )
    except Exception as e:
        print("Failed to read StationXML inventory file.")
        print(e)
        sys.exit()

    # Read SAC data files
    if(args.verbose):
        print("Reading SAC file(s):")
    sts = {}
    for file in args.sac:
        if(args.verbose):
            print("\t%s" % file)
        sts[file] = obspy.read(file)
        # TODO: figure out how this complains if a response is missing and
        # use that to throw an error (or just skip the station)
        sts[file].attach_response(inv)

    # Loop over streams (corresponding to files)
    if(args.verbose):
        print("Estimating tilt:")
    for name, st in sts.items():
        tr = st[0]
        if(args.verbose):
            print("\t%s.%s.%s.%s - %s" % (tr.stats.network, tr.stats.station, tr.stats.location, tr.stats.channel, tr.stats.starttime))
        # Decide on low pass filter period if necessary
        if args.passband_short_period is None:
            # Generate instrument responses to translation (velocity) and tilt
            Ca, freq = tr.stats.response.get_evalresp_response(
                t_samp=1, nfft=2 ** 13, output="ACC"
            )
            Ct = Ca * -args.gravity
            Cv, freq = tr.stats.response.get_evalresp_response(
                t_samp=1, nfft=2 ** 13, output="VEL"
            )
            dbSplit = 10 * np.log10(np.abs(Ct[1:] / Cv[1:]))
            args.passband_short_period = 1.0 / np.max(freq[np.where(dbSplit >= 20)])
        else:
            args.passband_short_period = float(args.passband_short_period)
        
        # Decide on high pass filter period if necessary
        if args.passband_long_period is None:
            args.passband_long_period = tr.stats.npts / tr.stats.sampling_rate / 2
        else:
            args.passband_long_period = float(args.passband_long_period)

        if(args.verbose):
            print("\t\tPassband: %.2fs to %.2fs" % (args.passband_long_period, args.passband_short_period))
            
        # Linear detrend
        tr.detrend(type="linear")

        # Filter
        tr.filter(
            type="bandpass",
            freqmin=1.0/args.passband_long_period,
            freqmax=1.0/args.passband_short_period,
            corners=4,
            zerophase=False,
        )

        # Decimate
        fdec = 1

        # TODO: some smarter decimation frequency logic

        if(fdec > tr.stats.sampling_rate or fdec < 2*(1/args.passband_short_period)):
            print("\t\tFailed to find appropriate decimation frequency, skipping")
            continue

        tr.decimate(int(tr.stats.sampling_rate), no_filter=True)

        # Get instrument response to tilt sampled at proper frequencies
        Ca, freq = tr.stats.response.get_evalresp_response(tr.stats.delta, tr.stats.npts, output="ACC")
        Ct = -np.sqrt(args.gravity)*Ca
        DATA = np.fft.rfft(tr.data)[1:]/Ct[1:]
        DATA = np.append(0, DATA)
        tr.data = np.fft.irfft(DATA)

        if(args.output is None):
            tiltName = os.path.dirname(name) + "/tilt_" + os.path.basename(name)
        else:
            tiltName = args.output + "/tilt_" + os.path.basename(name)
            
        st.write(tiltName, format="SAC")

    return 0


if __name__ == "__main__":
    main()
