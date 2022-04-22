import virgo
import pandas as pd
import os
import datetime
from moving_average import moving_average
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import glob
import argparse
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u
import configparser
from tqdm import tqdm
import numpy as np
import plotly.express as px


def main():
    global obs
    
    #Parse console arguments
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-r", "--runs", help = 'Number of measurements performed')
    parser.add_argument("-d", "--duration", help ='Time per observation (s)')
    parser.add_argument("-t", "--timeSample", help ='Time per sample')
    parser.add_argument("-c", "--channels", help ='Number of sample channels')
    parser.add_argument("-p", "--prominence", help ='Minimum peak prominence for peakfinding')
    parser.add_argument("-w", "--width", help ='Minimum peak width for peakfinding')
    parser.add_argument("-mH", "--height", help ='Minimum peak height for peakfinding')
    parser.add_argument("-nO", "--noObservation", help ='Disables automatic data acquisition', action = 'store_true')
    parser.add_argument("-b", "--baseline", help ='Measure baseline. Automatic -nO', action = 'store_true')
    parser.add_argument("-sF", "--setupFolders", help = 'Generates default folder structure', action = 'store_true')
    parser.add_argument("-sC", "--setupConfig", help = 'Generates default configuration file', action = 'store_true')
    parser.add_argument("-pA", "--plotAll", help ='Plot all observation files', action = 'store_true')
    parser.add_argument("-nP", "--noPlot", help ='Disables automatic plotting', action = 'store_true')
    parser.add_argument("-map", "--heatmap", help ='Generate heatmap from recorded data', action = 'store_true')
    parser.add_argument("-i", "--infinite", help ='Infinite runs. Continue until stopped', action = 'store_true')
    
    
    
    args = parser.parse_args()
    readConfig()    
    
    if args.setupConfig:
        makeDefaultConfigurationFile()
    
    if args.setupFolders:
        makeFileStructure()
    
    if args.duration:
        obs['duration'] = float(args.duration)
        
    if args.timeSample:
        obs['t_sample'] = float(args.timeSample)
    
    if args.channels:
        obs['channels'] = int(args.channels)
    
    if args.prominence:    
        global peakProminence
        peakProminence = float(args.prominence)
    
    if args.width:    
        global peakWidth
        peakWidth = float(args.width)
    
    if args.height:    
        global peakHeight
        peakHeight = float(args.height)
    
    if args.baseline:
        baseline(obs)
        return
    
    if args.plotAll:
        files = glob.glob(storagePath + '/observationData/observations/' + '*.dat')
        for filePath in tqdm(files):
            analyzeData(filePath.split('/')[(len(filePath.split('/')) - 1)], obs)
    
    if args.heatmap:
        heatmap()
        return
    
    if args.noObservation:
        return
    
    if args.infinite:
        i = 0
        while True:
            print("Starting observation {0}".format(i + 1))
            observe(obs, args.noPlot)
            i += 1
    
    global runs
    if args.runs:
        runs = int(args.runs)
    i = 0
    while i < runs:
        print("Starting observation {0} of {1}".format(i + 1, runs))
        observe(obs, args.noPlot)
        i += 1
    print("Done")

def readConfig():
    global storagePath
    global runs
    global peakProminence
    global peakWidth
    global peakHeight
    global alt
    global az
    global lat
    global long
    global elevation
    global obs
    
    #Read file and assign parameters
    config = configparser.ConfigParser()
    config.read('configurations.ini')
    obs = {
        'dev_args': config['ObservationParameters']['dev_args'],
        'rf_gain': float(config['ObservationParameters']['rf_gain']),
        'if_gain': float(config['ObservationParameters']['if_gain']),
        'bb_gain': float(config['ObservationParameters']['bb_gain']),
        'frequency': float(config['ObservationParameters']['frequency']),
        'bandwidth': float(config['ObservationParameters']['bandwidth']),
        'channels': int(config['ObservationParameters']['channels']),
        't_sample': float(config['ObservationParameters']['t_sample']),
        'duration': float(config['ObservationParameters']['duration']),
        'loc': config['ObservationParameters']['loc'],
        'ra_dec': config['ObservationParameters']['ra_dec'],
        'az_alt': config['ObservationParameters']['az_alt']
    }
    print('\nConfiguration file loaded succesfully\n')    
    
    storagePath = config['Storage']['storagePath']
    runs = int(config['ObservationParameters']['runCount'])
    peakProminence = float(config['Peakfinding']['prominence'])
    peakWidth = float(config['Peakfinding']['Width'])
    peakHeight = float(config['Peakfinding']['Height'])
    az = float(config['ObservationParameters']['az_alt'].split(',')[0])
    alt = float(config['ObservationParameters']['az_alt'].split(',')[1])
    lat = float(config['ObservationParameters']['loc'].split(',')[0])
    long = float(config['ObservationParameters']['loc'].split(',')[1])
    elevation = float(config['ObservationParameters']['loc'].split(',')[2])
    
def makeDefaultConfigurationFile():
    config_file = configparser.ConfigParser()
    
    #Generate file contents
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
    with open(r"configurations.ini", 'w') as configfileObj:
        config_file.write(configfileObj)
        configfileObj.flush()
        configfileObj.close()
    print("Config file 'configurations.ini' created")
    

def makeFileStructure():
        if not os.path.exists(storagePath):
            print('Generating File Structure')
            os.makedirs(storagePath + 'observationData/baseline/')
            os.makedirs(storagePath + 'observationData/observations/')
            os.makedirs(storagePath + 'csvFiles/timeSeries/')
            os.makedirs(storagePath + 'csvFiles/spectra/')
            os.makedirs(storagePath + 'plots/Virgo/')
            os.makedirs(storagePath + 'plots/PeakIdentification/')
            os.makedirs(storagePath + 'plots/heatmaps/')
            df = pd.DataFrame({'time':[],
                               'ra':[],
                               'dec':[],
                               'l':[],
                               'b':[],
                               'hz':[]})
            df.to_csv(storagePath + 'allObservations.csv')


def baseline(obs):
    #Recorded baseline and store in baseine.dat
    print('\n--------Starting baseline measurement--------\n')
    virgo.observe(obs_parameters=obs, obs_file=storagePath + 'observationData/baseline/baseline.dat', spectrometer='wola', start_in = 0)
    print('\n--------Calibration finished--------\n')
    
    
def observe(obs, noPlot):
    print('\n--------Starting observation--------')
    print('Current time: {0}.\nEstimated completion: {1}\n'.format(Time.now(), Time.now() + datetime.timedelta(seconds = obs['duration'])))
    #Create file 
    filename = 'observation ' + str(Time.now()).replace(':', '_') + obs['az_alt'].replace(' ', '_') + '.dat'
    fp = open(storagePath + 'observationData/observations/' + filename, 'x')
    fp.close()
    #Begin data aqcuisition
    virgo.observe(obs_parameters=obs, obs_file=storagePath + 'observationData/observations/' + filename, spectrometer='wola', start_in = 0)
    print('\n--------Observation finished--------\n')
    if noPlot == False:
        analyzeData(filename, obs)
    

def analyzeData(filename, obs):
    #Extract time of recording from filename
    recordedTime = filename.split(' ')[1] + ' '
    recordedTime += filename.split(' ')[2].replace('_', ':').split('.')[0] + '.'
    recordedTime += filename.split(' ')[2].replace('_', ':').split('.')[1]
    
    #Extract csv files from .dat file
    virgo.plot(obs_parameters=obs, 
               dB =  True,
               n = 20,
               m = 20,
               f_rest = 1420.4057517667e6,
               obs_file = storagePath + 'observationData/observations/' + filename,
               cal_file = storagePath + 'observationData/baseline/baseline.dat',
               spectra_csv = storagePath + 'csvFiles/spectra/spectrum ' + recordedTime.replace(':', '_') + '.csv',
               power_csv = storagePath + 'csvFiles/timeSeries/timeSeries ' + recordedTime.replace(':', '_') + '.csv',
               plot_file = storagePath + 'plots/Virgo/plot' + recordedTime.replace(':', '_') + '.png')
    
    #Add headers to csv file
    headerList = ['hz', 'spectrum', 'baseline', 'calibrated_spectrum']
    file = pd.read_csv(storagePath + 'csvFiles/spectra/spectrum ' + recordedTime.replace(':', '_') + '.csv') 
    file.to_csv(storagePath + "csvFiles/gfg2.csv", header=headerList, index=False)

    #Get data with headers from csv and caculate moving average
    _temporaryData = pd.read_csv(storagePath + "csvFiles/gfg2.csv")
    calibrationSpectrum = _temporaryData[:]['calibrated_spectrum']
    spectrumAverage = moving_average(calibrationSpectrum, 20)
    
    #Find peaks in spectrum
    peaks = find_peaks(spectrumAverage, prominence = peakProminence, width= peakWidth , height = peakHeight)
    
    try:
        #This fails if az_alt has not been stored in the filename
        filename.split(' ')[4]
        recordedAz = float(filename.split(' ')[3].replace('_', ' ').split(',')[0])
        recordedAlt = float(filename.split(' ')[3].replace('_', ' ').split(',')[1])
    except:
        print('No historic orientation data provided. Using default parameters')
        recordedAz = az
        recordedAlt = alt
        
    #Calculate galactic coordinates for observation
    raDec = oldTimeEquatorial(alt = recordedAlt, az = recordedAz, lat= lat, lon = long, height = elevation, time = recordedTime)
    galCoord = virgo.galactic(raDec[0], raDec[1])
    
    #Store recorded Peaks in a Dataframe and append Dataframe to overall csv file
    if peaks[0].size > 0:        
        writeFrame = pd.DataFrame(
            {'time': [],
              'ra': [],
              'dec': [],
              'l': [],
              'b': [],
              'hz': []})
        for i in range(0, len(peaks[0])):
            export = {'time': [recordedTime],
                      'ra': [raDec[0]],
                      'dec': [raDec[1]],
                      'l': [galCoord[0]],
                      'b': [galCoord[1]],
                      'hz': [_temporaryData[:]['hz'][peaks[0][i]]]}
            writeFrame = pd.concat([writeFrame, pd.DataFrame(export)])
        writeFrame.to_csv(storagePath + 'allObservations.csv', mode='a', index=False, header=False)
    
    #Generate and output peakfinding plot
    plt.figure(dpi = 250)
    plt.plot(_temporaryData[:]['hz'], spectrumAverage)
    plt.scatter(_temporaryData[:]['hz'][peaks[0]], spectrumAverage[peaks[0]], marker = 7, color = 'r')
    plt.title('Spectrum taken at ' + recordedTime)
    plt.ylabel('Signal to Noise ratio (1)')
    plt.xlabel('Frequency (hz)')
    plt.savefig(storagePath + 'plots/peakIdentification/plot' + recordedTime.replace(':', '_') + '.png')
    plt.close('all')
    plt.clf()
    
def oldTimeEquatorial(alt, az, lat, lon, time, height=0):
	'''
	Takes observer's location and Alt/Az as input and returns RA/Dec as a tuple.
	
	Args:
		alt: float. Altitude [deg]
		az: float. Azimuth [deg]
		lat: float. Observer latitude [deg]
		lon: float. Observer longitude [deg]
		height: float. Observer elevation [m]
	'''

	# Set observer location
	loc = EarthLocation(lat=lat * u.deg, lon=lon * u.deg, height=height * u.m)

	# Get current system time
	current_time = time

	# Compute Alt/Az
	AltAzcoordiantes = SkyCoord(alt=alt * u.deg, az=az * u.deg,
				    obstime=current_time, frame='altaz', location=loc)

	# Transform to RA/Dec
	c = AltAzcoordiantes.icrs

	ra = c.ra.hour
	dec = c.dec.deg

	# Return position as tuple
	return (ra, dec)

def heatmap():
    #Read overall csv for ra_dec and galactic coordiantes
    data = pd.read_csv(storagePath + 'allObservations.csv')
    dfGal = pd.DataFrame({'l': data['l'][:], 'b': data['b'][:], 'hz': data['hz'][:]})
    dfRa = pd.DataFrame({'ra': data['ra'][:], 'dec': data['dec'][:], 'hz': data['hz'][:]})
    
    #Pivot table
    dfGalPivoted = pd.pivot_table(dfGal, values = 'hz', index = ['l', 'b'], aggfunc = np.mean)
    dfRaPivoted = pd.pivot_table(dfRa, values = 'hz', index = ['ra', 'dec'], aggfunc = np.mean)
    
    plt.figure(dpi = 250)
    
    #Generate heatmap
    figGal = px.density_heatmap(dfGalPivoted.reset_index(), x="l", y="b", z="hz", histfunc="avg", nbinsx = 100, nbinsy = 100)
    figRa = px.density_heatmap(dfRaPivoted.reset_index(), x="ra", y="dec", z="hz", histfunc="avg", nbinsx = 100, nbinsy = 100)
    
    #Export image
    figGal.write_image(storagePath + 'plots/heatmaps/ ' + str(Time.now()).replace(':', '_') + 'heatmapGalCoord.png') 
    figRa.write_image(storagePath + 'plots/heatmaps/' + str(Time.now()).replace(':', '_') + 'heatmapRaDec.png')
    
    plt.close('all')
    print("Heatmap export successfull") 

if __name__ == "__main__":
    main()

