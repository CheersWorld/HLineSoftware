import virgo
import pandas as pd
import datetime
from moving_average import moving_average
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import glob
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u
from tqdm import tqdm
import numpy as np
import plotly.express as px
import globalVarsManager as gvars
import configurationManager as cm

"""
    This file handles observing and analysing data. Handling configuration files and console arguments
    is performed by configurationManager.py.
    
    Shared variables between the two are stored in globalVarsManager.
"""


def main():
    cm.setParameters()

    if gvars.args.baseline:
        baseline(gvars.obs)
        return

    if gvars.args.plotAll:
        files = glob.glob(
            gvars.storagePath + "/observationData/observations/" + "*.dat"
        )
        for filePath in tqdm(files):
            analyzeData(filePath.split("/")[(len(filePath.split("/")) - 1)], gvars.obs)

    if gvars.args.heatmap:
        createHeatmaps()
        return

    if gvars.args.noObservation:
        return

    if gvars.args.infinite:
        i = 0
        while True:
            print("Starting observation {0}".format(i + 1))
            observe(gvars.obs, gvars.args.noPlot)
            i += 1

    if gvars.args.runs:
        gvars.runs = int(gvars.args.runs)
    i = 0
    while i < gvars.runs:
        print("Starting observation {0} of {1}".format(i + 1, gvars.runs))
        observe(gvars.obs, gvars.args.noPlot)
        i += 1
    print("Done")


"""
    Record a baseline measurement and store it in baseline.dat
"""


def baseline(obs):
    print("\n--------Starting baseline measurement--------\n")
    virgo.observe(
        obs_parameters=gvars.obs,
        obs_file=gvars.storagePath + "observationData/baseline/baseline.dat",
        spectrometer="wola",
        start_in=0,
    )
    print("\n--------Calibration finished--------\n")


# Observe and analyze data
def observe(obs, noPlot):
    print("\n--------Starting observation--------")
    print(
        "Current time: {0}.\nEstimated completion: {1}\n".format(
            Time.now(), Time.now() + datetime.timedelta(seconds=obs["duration"])
        )
    )
    # Create file
    filename = "observation " + str(Time.now()).replace(":", "_") + ".dat"
    fp = open(gvars.storagePath + "observationData/observations/" + filename, "x")
    fp.close()
    # Begin data aqcuisition
    virgo.observe(
        obs_parameters=gvars.obs,
        obs_file=gvars.storagePath + "observationData/observations/" + filename,
        spectrometer="wola",
        start_in=0,
    )
    print("\n--------Observation finished--------\n")
    if noPlot == False:
        analyzeData(filename, gvars.obs)


def analyzeData(filename, obs):
    # Extract time of recording from filename
    recordedTime = filename.split(" ")[1] + " "
    recordedTime += filename.split(" ")[2].replace("_", ":").split(".")[0] + "."
    recordedTime += filename.split(" ")[2].replace("_", ":").split(".")[1]

    # Extract csv files from .dat file
    virgo.plot(
        obs_parameters=gvars.obs,
        dB=True,
        n=20,
        m=20,
        slope_correction=True,
        f_rest=1420.4057517667e6,
        obs_file=gvars.storagePath + "observationData/observations/" + filename,
        cal_file=gvars.storagePath + "observationData/baseline/baseline.dat",
        spectra_csv=gvars.storagePath
        + "csvFiles/spectra/spectrum "
        + recordedTime.replace(":", "_")
        + ".csv",
        power_csv=gvars.storagePath
        + "csvFiles/timeSeries/timeSeries "
        + recordedTime.replace(":", "_")
        + ".csv",
        plot_file=gvars.storagePath
        + "plots/Virgo/plot"
        + recordedTime.replace(":", "_")
        + ".png",
    )

    # Add headers to csv file
    headerList = ["hz", "spectrum", "baseline", "calibrated_spectrum"]
    file = pd.read_csv(
        gvars.storagePath
        + "csvFiles/spectra/spectrum "
        + recordedTime.replace(":", "_")
        + ".csv"
    )
    file.to_csv(gvars.storagePath + "csvFiles/gfg2.csv", header=headerList, index=False)

    # Get data with headers from csv and caculate moving average
    _temporaryData = pd.read_csv(gvars.storagePath + "csvFiles/gfg2.csv")
    calibrationSpectrum = _temporaryData[:]["calibrated_spectrum"]
    spectrumAverage = moving_average(calibrationSpectrum, 20)

    # Find peaks in spectrum
    peaks = find_peaks(
        spectrumAverage,
        prominence=gvars.peakProminence,
        width=gvars.peakWidth,
        height=gvars.peakHeight,
    )
    try:
        with open(
            gvars.storagePath
            + "observationData/observations/"
            + filename.split(".")[0]
            + "."
            + filename.split(".")[1]
            + ".header"
        ) as headerFile:
            header = headerFile.readlines()
        recordedAz = float(header[12].split("=")[1].split(",")[0])
        recordedAlt = float(header[12].split("=")[1].split(",")[1])
    except:
        print("Header file for " + filename + " not found")
        recordedAz = gvars.az
        recordedAlt = gvars.alt
    # Calculate galactic coordinates for observation
    raDec = oldTimeEquatorial(
        alt=recordedAlt,
        az=recordedAz,
        lat=gvars.lat,
        lon=gvars.long,
        height=gvars.elevation,
        time=recordedTime,
    )
    galCoord = virgo.galactic(raDec[0], raDec[1])

    # Store recorded Peaks in a Dataframe and append Dataframe to overall csv file
    if peaks[0].size > 0:
        writeFrame = pd.DataFrame(
            {"time": [], "ra": [], "dec": [], "l": [], "b": [], "hz": [], "height": []}
        )
        for i in range(0, len(peaks[0])):
            export = {
                "time": [recordedTime],
                "ra": [raDec[0]],
                "dec": [raDec[1]],
                "l": [galCoord[0]],
                "b": [galCoord[1]],
                "hz": [_temporaryData[:]["hz"][peaks[0][i]]],
                "height": [peaks[1]["peak_heights"][i]],
            }
            writeFrame = pd.concat([writeFrame, pd.DataFrame(export)])
        writeFrame.to_csv(
            gvars.storagePath + "allObservations.csv",
            mode="a",
            index=False,
            header=False,
        )

    # Generate and output peakfinding plot
    plt.figure(dpi=250)
    plt.plot(_temporaryData[:]["hz"], spectrumAverage, label="moving average n = 20")
    plt.scatter(
        _temporaryData[:]["hz"][peaks[0]],
        spectrumAverage[peaks[0]],
        marker=7,
        color="r",
    )
    plt.title("Spectrum taken at " + recordedTime)
    plt.ylabel("Signal to Noise ratio (1)")
    plt.xlabel("Frequency (hz)")
    plt.savefig(
        gvars.storagePath
        + "plots/peakIdentification/plot"
        + recordedTime.replace(":", "_")
        + ".png"
    )
    plt.close("all")
    plt.clf()


def oldTimeEquatorial(alt, az, lat, lon, time, height=0):
    """
    Takes observer's location and Alt/Az as input and returns RA/Dec as a tuple.

    Args:
            alt: float. Altitude [deg]
            az: float. Azimuth [deg]
            lat: float. Observer latitude [deg]
            lon: float. Observer longitude [deg]
            height: float. Observer elevation [m]
    """

    # Set observer location
    loc = EarthLocation(lat=lat * u.deg, lon=lon * u.deg, height=height * u.m)

    # Get current system time
    current_time = time

    # Compute Alt/Az
    AltAzcoordiantes = SkyCoord(
        alt=alt * u.deg,
        az=az * u.deg,
        obstime=current_time,
        frame="altaz",
        location=loc,
    )

    # Transform to RA/Dec
    c = AltAzcoordiantes.icrs

    ra = c.ra.hour
    dec = c.dec.deg

    # Return position as tuple
    return (ra, dec)


"""
Creates 2D histograms (heatmaps) from the data in allObservations.csv.

"""


def createHeatmaps():
    # Read overall csv for ra_dec and galactic coordiantes
    data = pd.read_csv(gvars.storagePath + "allObservations.csv")

    # Generate heatmaps for the data and store in the appropriate location
    plt.figure(dpi=250)

    heatmap(data, "l", "b", "hz").write_image(
        gvars.storagePath
        + "plots/heatmaps/frequency/"
        + str(Time.now()).replace(":", "_")
        + "heatmapGalCoord.png"
    )

    heatmap(data, "ra", "dec", "hz").write_image(
        gvars.storagePath
        + "plots/heatmaps/frequency/"
        + str(Time.now()).replace(":", "_")
        + "heatmapRaDec.png"
    )

    heatmap(data, "l", "b", "height").write_image(
        gvars.storagePath
        + "plots/heatmaps/height/"
        + str(Time.now()).replace(":", "_")
        + "heatmapGalCoord.png"
    )

    heatmap(data, "ra", "dec", "height").write_image(
        gvars.storagePath
        + "plots/heatmaps/height/"
        + str(Time.now()).replace(":", "_")
        + "heatmapRaDec.png"
    )

    plt.close("all")
    print("Heatmap export successfull")


def heatmap(data, x, y, z):
    df = pd.DataFrame({str(x): data[x][:], str(y): data[y][:], str(z): data[z][:]})
    dfPivoted = pd.pivot_table(df, values=z, index=[x, y], aggfunc=np.mean)
    return px.density_heatmap(
        dfPivoted.reset_index(),
        x=x,
        y=y,
        z=z,
        histfunc="avg",
        nbinsx=gvars.nBinsX,
        nbinsy=gvars.nBinsY,
    )


if __name__ == "__main__":
    main()
