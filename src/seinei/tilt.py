import obspy
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import sys


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
        tilt.filter("lowpass", freq=1.0 / 30, corners=3)
        tilt.detrend(type="polynomial", order=1)
        tilt.taper(0)
        # print("Hi :)")
        # tilt.filter("bandpass", freqmin=1.0/2000, freqmax=1.0/100)
        # print("CHANGE FILTER BACK")
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
        tilt.remove_response(
            output="DEF", water_level=None, taper=False, zero_mean=False
        )
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
            # tr.filter(type="highpass", freq=1.0/(60*5), zerophase=False)
            tr.data = tr.data * ((2 * np.pi * f0) ** 2) / sensor.getfloat("gravity")

            # Remove linear trend from tilt (equivalent to removing mean from voltage)
            # tr.detrend(type="polynomial", order=1)

    return tilts
