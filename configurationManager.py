# -*- coding: utf-8 -*-
import configparser
import os
import glovalVarsManager as gvars
import pandas as pd
import argparse


def makeDefaultConfigurationFile():
    config_file = configparser.ConfigParser()

    # Generate file contents
    config_file.add_section("ObservationParameters")
    config_file.add_section("Storage")
    config_file.add_section("Peakfinding")

    config_file.set("ObservationParameters", "dev_args", "rtl,bias=1")
    config_file.set("ObservationParameters", "rf_gain", "30")
    config_file.set("ObservationParameters", "if_gain", "25")
    config_file.set("ObservationParameters", "bb_gain", "18")
    config_file.set("ObservationParameters", "frequency", "1420e6")
    config_file.set("ObservationParameters", "bandwidth", "2.4e6")
    config_file.set("ObservationParameters", "channels", "2048")
    config_file.set("ObservationParameters", "t_sample", ".2")
    config_file.set("ObservationParameters", "duration", "120")
    config_file.set("ObservationParameters", "loc", "0, 0, 0")
    config_file.set("ObservationParameters", "ra_dec", "")
    config_file.set("ObservationParameters", "az_alt", "0, 0")
    config_file.set("ObservationParameters", "runCount", "1")
    config_file.set("Storage", "storagePath", "../HLineObservations/")
    config_file.set("Peakfinding", "Prominence", "0")
    config_file.set("Peakfinding", "Width", "10")
    config_file.set("Peakfinding", "Height", "3")

    # SAVE CONFIG FILE
    with open(r"configurations.ini", "w") as configfileObj:
        config_file.write(configfileObj)
        configfileObj.flush()
        configfileObj.close()
    print("Config file 'configurations.ini' created")


def makeFileStructure():
    if not os.path.exists(gvars.storagePath):
        print("Generating File Structure")
        os.makedirs(gvars.storagePath + "observationData/baseline/")
        os.makedirs(gvars.storagePath + "observationData/observations/")
        os.makedirs(gvars.storagePath + "csvFiles/timeSeries/")
        os.makedirs(gvars.storagePath + "csvFiles/spectra/")
        os.makedirs(gvars.storagePath + "plots/Virgo/")
        os.makedirs(gvars.storagePath + "plots/PeakIdentification/")
        os.makedirs(gvars.storagePath + "plots/heatmaps/")
        df = pd.DataFrame({"time": [], "ra": [], "dec": [], "l": [], "b": [], "hz": []})
        df.to_csv(gvars.storagePath + "allObservations.csv")


def parseArguments():

    # Parse console arguments
    parser = argparse.ArgumentParser()

    parser.add_argument("-r", "--runs", help="Number of measurements performed")
    parser.add_argument("-d", "--duration", help="Time per observation (s)")
    parser.add_argument("-t", "--timeSample", help="Time per sample")
    parser.add_argument("-c", "--channels", help="Number of sample channels")
    parser.add_argument(
        "-p", "--prominence", help="Minimum peak prominence for peakfinding"
    )
    parser.add_argument("-w", "--width", help="Minimum peak width for peakfinding")
    parser.add_argument("-mH", "--height", help="Minimum peak height for peakfinding")
    parser.add_argument(
        "-nO",
        "--noObservation",
        help="Disables automatic data acquisition",
        action="store_true",
    )
    parser.add_argument(
        "-b", "--baseline", help="Measure baseline. Automatic -nO", action="store_true"
    )
    parser.add_argument(
        "-sF",
        "--setupFolders",
        help="Generates default folder structure",
        action="store_true",
    )
    parser.add_argument(
        "-sC",
        "--setupConfig",
        help="Generates default configuration file",
        action="store_true",
    )
    parser.add_argument(
        "-pA", "--plotAll", help="Plot all observation files", action="store_true"
    )
    parser.add_argument(
        "-nP", "--noPlot", help="Disables automatic plotting", action="store_true"
    )
    parser.add_argument(
        "-map",
        "--heatmap",
        help="Generate heatmap from recorded data",
        action="store_true",
    )
    parser.add_argument(
        "-i",
        "--infinite",
        help="Infinite runs. Continue until stopped",
        action="store_true",
    )
    gvars.args = parser.parse_args()
    if gvars.args.setupConfig:
        makeDefaultConfigurationFile()

    if gvars.args.setupFolders:
        makeFileStructure()

    if gvars.args.duration:
        gvars.obs["duration"] = float(gvars.args.duration)

    if gvars.args.timeSample:
        gvars.obs["t_sample"] = float(gvars.args.timeSample)

    if gvars.args.channels:
        gvars.obs["channels"] = int(gvars.args.channels)

    if gvars.args.prominence:
        gvars.peakProminence = float(gvars.args.prominence)

    if gvars.args.width:
        gvars.peakWidth = float(gvars.args.width)

    if gvars.args.height:
        gvars.peakHeight = float(gvars.args.height)


def readConfig():
    # Read file and assign parameters
    config = configparser.ConfigParser()
    config.read("configurations.ini")
    gvars.obs = {
        "dev_args": config["ObservationParameters"]["dev_args"],
        "rf_gain": float(config["ObservationParameters"]["rf_gain"]),
        "if_gain": float(config["ObservationParameters"]["if_gain"]),
        "bb_gain": float(config["ObservationParameters"]["bb_gain"]),
        "frequency": float(config["ObservationParameters"]["frequency"]),
        "bandwidth": float(config["ObservationParameters"]["bandwidth"]),
        "channels": int(config["ObservationParameters"]["channels"]),
        "t_sample": float(config["ObservationParameters"]["t_sample"]),
        "duration": float(config["ObservationParameters"]["duration"]),
        "loc": config["ObservationParameters"]["loc"],
        "ra_dec": config["ObservationParameters"]["ra_dec"],
        "az_alt": config["ObservationParameters"]["az_alt"],
    }
    print("\nConfiguration file loaded succesfully\n")

    gvars.storagePath = config["Storage"]["storagePath"]
    gvars.runs = int(config["ObservationParameters"]["runCount"])
    gvars.peakProminence = float(config["Peakfinding"]["prominence"])
    gvars.peakWidth = float(config["Peakfinding"]["Width"])
    gvars.peakHeight = float(config["Peakfinding"]["Height"])
    gvars.az = float(config["ObservationParameters"]["az_alt"].split(",")[0])
    gvars.alt = float(config["ObservationParameters"]["az_alt"].split(",")[1])
    gvars.lat = float(config["ObservationParameters"]["loc"].split(",")[0])
    gvars.long = float(config["ObservationParameters"]["loc"].split(",")[1])
    gvars.elevation = float(config["ObservationParameters"]["loc"].split(",")[2])
