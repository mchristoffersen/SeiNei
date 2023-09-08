import obspy


def parseTimes(times):
    """Take list of strings reprenenting times (preferably ISO8601:2004 compliant) and return list of obspy.UTCDateTime objects.

    Args:
        times (list of str): list of time strings

    Returns:
        list of obspy time objects, created from the list of time strings

    """
    times = [obspy.UTCDateTime(time) for time in times]

    return times
