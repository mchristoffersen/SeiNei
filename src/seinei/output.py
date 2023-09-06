import obspy
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import sys


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
