import obspy
import obspy.clients.fdsn
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import sys


def parseTimes(times):
    """Take list of strings reprenenting times (preferably ISO8601:2004 compliant) and return list of obspy.UTCDateTime objects.

    Args:
        times (list of str): list of time strings

    Returns:
        list of obspy time objects, created from the list of time strings

    """
    times = [obspy.UTCDateTime(time) for time in times]

    return times


def downloadStreams(times, source):
    """Download data from IRIS

    Args:
        times (list of obspy.UTCDateTime): list of times to search for a window of data around
        source: dict with data source information, see readme.md for format

    Returns:
        streams: list of obspy.Stream objects containing the downloaded data
        downloadedTimes: times where a window of data was downloaded
    """
    client = obspy.clients.fdsn.Client(
        "IRIS"
    )  # Could remove hardcode to use other sources?
    window = source.getfloat("window")
    streams = []
    downloadedTimes = []

    for i, time in enumerate(times):
        print("Downloading record %d/%d" % (i + 1, len(times)), end="\r")

        try:
            stream = client.get_waveforms(
                source["network"],
                source["station"],
                source["site"],
                source["component"],
                time - window,
                time + window,
                minimumlength=2*window,
                attach_response=True,
            )
            downloadedTimes.append(time)
        except obspy.clients.fdsn.header.FDSNNoDataException:
            print("No complete data found for time: %s" % time, file=sys.stderr)
            continue

        stream.filter("lowpass", freq=1.0/4, corners=3)
        stream.resample(1, no_filter=True)  # decimate to 1 Hz
        print(stream)
        streams.append(stream)

    return streams, downloadedTimes


def calcTilt(streams, sensor):
    """Integrate seismometer records and multiply by constant to calculate tilt.

    Args:
        streams (list of obspy.Stream): raw seismometer data to convert to tilt
        sensor: dict with seismometer information, see readme.md for format

    Returns:
        list of obspy.Stream objects containing tilt data
    """
    tilts = [stream.copy() for stream in streams]

    # Tilt conversion
    S = 1 / sensor.getfloat("sensitivity")
    f0 = 1 / sensor.getfloat("cornerFreq")
    w0 = 2 * np.pi * f0
    bitDepth = sensor.getfloat("bitDepth")
    g = sensor.getfloat("gravity")

    for tilt in tilts:
        tilt.filter("lowpass", freq=1.0/30, corners=3)
        tilt.detrend(type="polynomial", order=1)
        tilt.taper(0)
        #print("Hi :)")
        #tilt.filter("bandpass", freqmin=1.0/2000, freqmax=1.0/100)
        #print("CHANGE FILTER BACK")
        # Check that corner frequency is < .5 Hz
        if f0 > 0.5:
            raise RuntimeError(
                "Seismometer corner freq must be less than .5 Hz to satisfy decimation assumption"
            )
        # Assuming fs for first trace is fs of all traces here...
        if tilt[0].stats.sampling_rate != int(tilt[0].stats.sampling_rate):
            raise RuntimeError(
                "Seismometer sampling freq must be integer to satisfy decimation assumption"
            )

        tilt.integrate(method="cumtrapz")

        for trace in tilt:
            trace.data = (-S * (w0 ** 2) / g) * (trace.data * bitDepth)

    return tilts


def calcTiltVolt(streams, sensor):
    """Integrate seismometer records and multiply by constant to calculate tilt.

    Args:
        streams (list of obspy.Stream): raw seismometer data to convert to tilt
        sensor: dict with seismometer information, see readme.md for format

    Returns:
        list of obspy.Stream objects containing tilt data
    """
    tilts = [stream.copy() for stream in streams]

    # Tilt conversion
    S = 1 / sensor.getfloat("sensitivity")
    f0 = 1 / sensor.getfloat("cornerFreq")
    w0 = 2 * np.pi * f0
    bitDepth = sensor.getfloat("bitDepth")
    g = sensor.getfloat("gravity")

    for tilt in tilts:
        for trace in tilt:
            # Convert counts to volts
            # There are a lot of assumptions here I'm not checking
            # Most importantly - is the first response stage M/S -> Volts?
            del trace.stats.response.response_stages[0]
        tilt.remove_response(output="DEF", water_level=None, taper=False, zero_mean=False)
        tilt.detrend(type="constant")
        tilt.filter("lowpass", freq=f0, corners=3)
        # Check that corner frequency is < .5 Hz
        if f0 > 0.5:
            raise RuntimeError(
                "Seismometer corner freq must be less than .5 Hz to satisfy decimation assumption"
            )
        # Assuming fs for first trace is fs of all traces here...
        if tilt[0].stats.sampling_rate != int(tilt[0].stats.sampling_rate):
            raise RuntimeError(
                "Seismometer sampling freq must be integer to satisfy decimation assumption"
            )
        tilt.integrate(method="cumtrapz")

        for trace in tilt:
            trace.data = (-S * (w0 ** 2) / g) * (trace.data)

    return tilts


def calcTilt_xferFunction(streams, sensor):
    """Modify velocity transfer function to calculate tilt.

    Args:
        streams (list of obspy.Stream): raw seismometer data to convert to tilt
        sensor: dict with seismometer information, see readme.md for format

    Returns:
        list of obspy.Stream objects containing tilt data
    """
    # streams must have responses attached

    # Should be checking various assumptions here
    # 1. First response stage is from M/S -> V
    # 2. First response stage is poles and zeros form
    # Other things?

    # Work on a copy
    tilts = [st.copy() for st in streams]

    f0 = 1.0 / sensor.getfloat("cornerFreq")

    for st in tilts:
        for tr in st:
            # Add a pole at 0
            tr.stats.response.response_stages[0].poles.append(0)

            # Multiply norm factor by -9.81 (-g)
            tr.stats.response.response_stages[0].normalization_factor *= -9.81

            # Specify the correct units
            tr.stats.response.response_stages[0].input_units = "RAD"
            tr.stats.response.response_stages[
                0
            ].input_units_description = "tilt in radians"

            # Low pass filter
            tr.filter(type="lowpass", freq=f0, zerophase=False)

            # decimate to 1 Hz
            tr.resample(1, no_filter=True, window="hann")

            # Remove instrument response
            tr.remove_response(output="DEF", water_level=None, zero_mean=True)

            # Remove linear trend from tilt (equivalent to removing mean from voltage)
            # tr.detrend(type="polynomial", order=1)

    return tilts


def calcTilt_disp(streams, sensor):
    """Modify velocity transfer function to calculate tilt.

    Args:
        streams (list of obspy.Stream): raw seismometer data to convert to tilt
        sensor: dict with seismometer information, see readme.md for format

    Returns:
        list of obspy.Stream objects containing tilt data
    """
    # streams must have responses attached

    # Should be checking various assumptions here
    # 1. First response stage is from M/S -> V
    # 2. First response stage is poles and zeros form
    # Other things?

    # Work on a copy
    tilts = [st.copy() for st in streams]

    f0 = 1.0 / sensor.getfloat("cornerFreq")

    for st in tilts:
        for tr in st:
            # Remove response
            tr.remove_response(output="DEF", water_level=10)
            tr.integrate(method="cumtrapz")
            tr.filter(type="lowpass", freq=f0, zerophase=False)
            #tr.filter(type="highpass", freq=1.0/(60*5), zerophase=False)
            tr.data = tr.data * ((2*np.pi*f0)**2)/sensor.getfloat("gravity")

            # Remove linear trend from tilt (equivalent to removing mean from voltage)
            # tr.detrend(type="polynomial", order=1)

    return tilts


def saveTilts(tilts, source):
    """Save calculated tilts to a miniseed file.

    Args:
        tilts (list of obspy.Stream): calculated seismometer tilts
        source: dict with data source information, see readme.md for format

    Returns:
        0 on success
    """
    # Enforce 32 bit floating point for all streams
    for tilt in tilts:
        for trace in tilt:
            trace.data = trace.data.astype(np.float32)

    # Concatenate all streams
    sts = tilts[0]
    for tilt in tilts[1:]:
        sts += tilt

    # Save as miniseed, encoding=4 means IEEE 32 bit float
    sts.write(
        "%s_%s_tilt.mseed" % (source["network"], source["station"]),
        format="MSEED",
        encoding=4,
    )

    return 0


def genFigure(streams, tilts, times, source):
    """Generate a standard set of tilt plots as a multi-page PDF.

    Args:
        streams (list of obspy.Stream): raw seismometer data to convert to tilt
        tilts (list of obspy.Stream): calculated seismometer tilts
        times (list of obspy.UTCDateTime): list of times to search for a window of data around
        source: dict with data source information, see readme.md for format

    Returns:
        0 on success
    """

    # Stacked tilt figure
    stackE = np.zeros((len(tilts), len(tilts[0][0].data)), dtype=np.float64)
    stackN = np.zeros_like(stackE)

    try:
        for i, tilt in enumerate(tilts):
            # Add in a check to make sure streams only have one trace per component?
            stackE[i, :] = tilt.select(component="E")[0].data
            stackN[i, :] = tilt.select(component="N")[0].data
    except IndexError:
        for i, tilt in enumerate(tilts):
            stackE[i, :] = tilt.select(component="1")[0].data
            stackN[i, :] = tilt.select(component="2")[0].data

    meanE = np.mean(stackE, axis=0)
    meanN = np.mean(stackN, axis=0)

    stdE = np.std(stackE, axis=0)
    stdN = np.std(stackN, axis=0)

    t = (tilts[0][0].times() - source.getfloat("window")) / 60

    fig0, ax = plt.subplots(2, 1, figsize=(8, 4), dpi=300)
    # fig0, ax = plt.subplots(1, 1, figsize=(8, 2), dpi=300)
    # ax = [ax]

    fig0.suptitle(
        "%s %s Stacked Apparent-Tilt (%d Events)"
        % (source["network"], source["station"], len(tilts))
    )

    ax[0].plot(t, meanE * 1e9, "k", rasterized=True)
    ax[0].fill_between(
        t,
        (meanE - stdE) * 1e9,
        (meanE + stdE) * 1e9,
        color="grey",
        alpha=0.3,
        rasterized=True,
    )
    ax[0].legend(["$\mu$", "$\pm \sigma$"])
    ax[0].set_ylabel("HHE Tilt (nrad)")
    ax[0].set_xlim(t[0], t[-1])
    ax[0].grid(True, linestyle="--")

    ax[1].plot(t, meanN * 1e9, "k", rasterized=True)
    ax[1].fill_between(
        t,
        (meanN - stdN) * 1e9,
        (meanN + stdN) * 1e9,
        color="grey",
        alpha=0.3,
        rasterized=True,
    )
    ax[1].legend(["$\mu$", "$\pm \sigma$"])
    ax[1].set_ylabel("HHN Tilt (nrad)")
    ax[1].set_xlabel("Time (minutes)")
    ax[1].set_xlim(t[0], t[-1])
    ax[1].grid(True, linestyle="--")

    fig0.tight_layout()

    # All east tilts
    nrow = np.ceil(len(tilts) / 3).astype(np.int32)
    vdim = 2 * nrow
    fig1, ax = plt.subplots(nrow, 3, figsize=(8, vdim), dpi=300)
    ax = np.atleast_2d(ax)

    fig1.suptitle("East Tilt (nrad vs minutes)")

    for i in range(nrow):
        for j in range(3):
            if i * 3 + j < len(tilts):
                try:
                    tilt = tilts[i * 3 + j].select(component="E")[0].data
                except IndexError:
                    tilt = tilts[i * 3 + j].select(component="1")[0].data
                ax[i, j].plot(
                    t,
                    tilt * 1e9,
                    "k",
                    rasterized=True,
                )
                ax[i, j].set_title(
                    times[i * 3 + j].strftime("%Y-%m-%d %H:%M:%S"), fontsize=8
                )
                ax[i, j].set_xlim(t[0], t[-1])
            else:
                ax[i, j].axis("off")

    fig1.tight_layout()

    # All north tilts
    fig2, ax = plt.subplots(nrow, 3, figsize=(8, vdim), dpi=300)
    ax = np.atleast_2d(ax)

    fig2.suptitle("North Tilt (nrad vs minutes)")

    for i in range(nrow):
        for j in range(3):
            if i * 3 + j < len(tilts):
                try:
                    tilt = tilts[i * 3 + j].select(component="N")[0].data
                except IndexError:
                    tilt = tilts[i * 3 + j].select(component="2")[0].data
                ax[i, j].plot(
                    t,
                    tilt * 1e9,
                    "k",
                    rasterized=True,
                )
                ax[i, j].set_title(
                    times[i * 3 + j].strftime("%Y-%m-%d %H:%M:%S"), fontsize=8
                )
                ax[i, j].set_xlim(t[0], t[-1])
            else:
                ax[i, j].axis("off")

    fig2.tight_layout()

    # All east seismograms
    t = (streams[0][0].times() - source.getfloat("window")) / 60
    fig3, ax = plt.subplots(nrow, 3, figsize=(8, vdim), dpi=300)
    ax = np.atleast_2d(ax)

    fig3.suptitle("East Seismogram (counts vs minutes)")

    for i in range(nrow):
        for j in range(3):
            if i * 3 + j < len(streams):
                try:
                    stream = streams[i * 3 + j].select(component="E")[0].data
                except IndexError:
                    stream = streams[i * 3 + j].select(component="1")[0].data
                ax[i, j].plot(
                    t,
                    stream,
                    "k",
                    rasterized=True,
                )
                ax[i, j].set_title(
                    times[i * 3 + j].strftime("%Y-%m-%d %H:%M:%S"), fontsize=8
                )
                ax[i, j].set_xlim(t[0], t[-1])
            else:
                ax[i, j].axis("off")

    fig3.tight_layout()

    # All north seismograms
    fig4, ax = plt.subplots(nrow, 3, figsize=(8, vdim), dpi=300)
    ax = np.atleast_2d(ax)

    fig4.suptitle("North Seismogram (counts vs minutes)")

    for i in range(nrow):
        for j in range(3):
            if i * 3 + j < len(streams):
                try:
                    stream = streams[i * 3 + j].select(component="N")[0].data
                except IndexError:
                    stream = streams[i * 3 + j].select(component="2")[0].data
                ax[i, j].plot(
                    t,
                    stream,
                    "k",
                    rasterized=True,
                )
                ax[i, j].set_title(
                    times[i * 3 + j].strftime("%Y-%m-%d %H:%M:%S"), fontsize=8
                )
                ax[i, j].set_xlim(t[0], t[-1])
            else:
                ax[i, j].axis("off")

    fig4.tight_layout()

    pdf = PdfPages("%s_%s.pdf" % (source["network"], source["station"]))
    pdf.savefig(fig0)
    pdf.savefig(fig1)
    pdf.savefig(fig2)
    pdf.savefig(fig3)
    pdf.savefig(fig4)
    pdf.close()

    return 0
