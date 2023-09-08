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
        "-lpp",
        "--low_pass_period",
        default=None,
        help="Corner period for low-pass filter in s (default=auto-set)",
    )
    parser.add_argument(
        "-hpp",
        "--high_pass_period",
        default=None,
        help="Corner period for high-pass filter in s (default=auto-set)",
    )
    parser.add_argument(
        "-g",
        "--gravity",
        default=9.81,
        help="Gravitaitonal acceleration in m/s^2 (default=9.81)",
    )
    args = parser.parse_args()

    return args


def main():
    args = cli()

    # Read StationXML inventory file, complain if obspy doesn't like it
    try:
        inv = obspy.core.inventory.inventory.read_inventory(
            args.inventory, format="STATIONXML", level="response"
        )
    except Exception as e:
        print("Failed to read StationXML inventory file.")
        print(e)
        sys.exit()

    # Read SAC data files
    sts = {}
    for file in args.sac:
        sts[file] = obspy.read(file)
        # TODO: figure out how this complains if a response is missing and
        # use that to throw an error (or just skip the station)
        sts[file].attach_response(inv)

    # Loop over streams (corresponding to files)
    for name, st in sts.items():
        tr = st[0]
        # Decide on low pass filter period if necessary
        if args.low_pass_period is None:
            # Generate instrument responses to translation (velocity) and tilt
            Ca, freq = tr.stats.response.get_evalresp_response(
                t_samp=1, nfft=2 ** 13, output="ACC"
            )
            Ct = Ca * -args.gravity
            Cv, freq = tr.stats.response.get_evalresp_response(
                t_samp=1, nfft=2 ** 13, output="VEL"
            )
            dbSplit = 10 * np.log10(np.abs(Ct / Cv))
            args.low_pass_period = 1.0 / np.max(freq[np.where(dbSplit >= 20)])
            print("Low-pass period: %.2fs" % args.low_pass_period)

        # Decide on high pass filter period if necessary
        if args.high_pass_period is None:
            args.high_pass_period = tr.stats.npts / tr.stats.sampling_rate / 2
            print("High-pass period: %.2fs" % args.high_pass_period)

        # Linear detrend
        tr.detrend(type="linear")

        # Filter
        tr.filter(
            type="bandpass",
            freqmin=1.0/args.high_pass_period,
            freqmax=1.obspy plot stream0/args.low_pass_period,
            corners=4,
            zerophase=False,
        )

        # Decimate
        fdec = 1

        # TODO: some smarter decimation frequency logic

        if(fdec > tr.stats.sampling_rate or fdec < 2*(1/args.low_pass_period)):
            print("Failed to find appropriate decimation frequency, exiting")
            sys.exit()

        tr.decimate(int(tr.stats.sampling_rate), no_filter=True)

        # Get instrument response to tilt sampled at proper frequencies
        Ca, freq = tr.stats.response.get_evalresp_response(tr.stats.delta, tr.stats.npts, output="ACC")
        Ct = -np.sqrt(args.gravity)*Ca
        DATA = np.fft.rfft(tr.data)/Ct
        DATA[0] = 0  # Set DC to 0
        tr.data = np.fft.irfft(DATA)

        tiltName = os.path.dirname(name) + "/tilt_" + os.path.basename(name)
        st.write(tiltName, format="SAC")

    return 0


if __name__ == "__main__":
    main()
