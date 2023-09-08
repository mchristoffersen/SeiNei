import obspy.clients.fdsn


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
                minimumlength=2 * window,
                attach_response=True,
            )
            downloadedTimes.append(time)
        except obspy.clients.fdsn.header.FDSNNoDataException:
            print("No complete data found for time: %s" % time, file=sys.stderr)
            continue

        stream.filter("lowpass", freq=1.0 / 4, corners=3)
        stream.resample(1, no_filter=True)  # decimate to 1 Hz
        print(stream)
        streams.append(stream)

    return streams, downloadedTimes
