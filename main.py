import virgo
import pandas as pd
import os
from moving_average import moving_average
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import glob
import argparse
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import astropy.units as u
import configparser


def main():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-r", "--runs", help = 'Number of measurements performed')
    parser.add_argument("-d", "--duration", help ='Time per observation (s)')
    parser.add_argument("-t", "--timeSample", help ='Time per sample')
    parser.add_argument("-nO", "--noObservation", help ='Disables automatic data acquisition', action = 'store_true')
    parser.add_argument("-b", "--baseline", help ='Measure baseline. Automatic -nO', action = 'store_true')
    parser.add_argument("-sF", "--setupFolders", help = 'Generates default folder structure', action = 'store_true')
    parser.add_argument("-sC", "--setupConfig", help = 'Generates default configuration file', action = 'store_true')
    parser.add_argument("-pA", "--plotAll", help ='Plot all observation files. Automatic -nO', action = 'store_true')
    parser.add_argument("-nP", "--noPlot", help ='Disables automatic plotting', action = 'store_true')
    
    
    args = parser.parse_args()
    
    if args.setupConfig:
        makeDefaultConfigurationFile()

    global storagePath
    
    config = readConfig()
    obs = config[0]
    runs = int(config[1])
    storagePath = config[2]
    
    if args.setupFolders:
        makeFileStructure()
    
    if args.duration:
        obs['duration'] = float(args.duration)
    if args.timeSample:
        obs['t_sample'] = float(args.timeSample)
    
    if args.baseline:
        baseline(obs)

    if args.plotAll:
        files = glob.glob(storagePath + '/observationData/observations/' + '*.dat')
        for filePath in files:
            print('Analyzing ' + filePath)
            analyzeData(filePath.split('/')[(len(filePath.split('/')) - 1)], obs)
    
    if args.noObservation or args.baseline or args.plotAll:
        return
    
    if args.runs:
        runs = int(args.runs)
        
    i = 0
    while i < runs:
        print("Starting observation {0} of {1}".format(i + 1, runs))
        observe(obs, args.noPlot)
        i += 1
    
    print("Done")

def readConfig():
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
    return obs, config['RunParameters']['runCount'], config['StorageParameters']['storagePath']

def makeDefaultConfigurationFile():
    config_file = configparser.ConfigParser()
    
    # ADD SECTION
    config_file.add_section("ObservationParameters")
    config_file.add_section("RunParameters")
    config_file.add_section("StorageParameters")
    config_file.set("ObservationParameters", "dev_args", "rtl,bias=1")
    config_file.set("ObservationParameters", "rf_gain", "30")
    config_file.set("ObservationParameters", "if_gain", "25")
    config_file.set("ObservationParameters", "bb_gain", "18")
    config_file.set("ObservationParameters", "frequency", "1420e6")
    config_file.set("ObservationParameters", "bandwidth", "2.4e6")
    config_file.set("ObservationParameters", "channels", "2048")
    config_file.set("ObservationParameters", "t_sample", ".5")
    config_file.set("ObservationParameters", "duration", "60")
    config_file.set("ObservationParameters", "loc", "")
    config_file.set("ObservationParameters", "ra_dec", "")
    config_file.set("ObservationParameters", "az_alt", "")
    config_file.set("RunParameters", "runCount", "1")
    config_file.set("StorageParameters", "storagePath", "../HLineObservations/")
    
    # SAVE CONFIG FILE
    with open(r"configurations.ini", 'w') as configfileObj:
        config_file.write(configfileObj)
        configfileObj.flush()
        configfileObj.close()
    
    print("Config file 'configurations.ini' created")
    
    # PRINT FILE CONTENT
    read_file = open("configurations.ini", "r")
    content = read_file.read()
    print("Content of the config file are:\n")
    print(content)
    read_file.flush()
    read_file.close()

def makeFileStructure():
        if not os.path.exists(storagePath):
            print('Generating File Structure')
            os.makedirs(storagePath + 'observationData/baseline/')
            os.makedirs(storagePath + 'observationData/observations/')
            os.makedirs(storagePath + 'csvFiles/timeSeries/')
            os.makedirs(storagePath + 'csvFiles/spectra/')
            os.makedirs(storagePath + 'plots/Virgo/')
            os.makedirs(storagePath + 'plots/PeakIdentification/')
            df = pd.DataFrame({'time':[],
                               'ra':[],
                               'dec':[],
                               'l':[],
                               'b':[],
                               'hz':[]})
            print(df)
            df.to_csv(storagePath + 'allObservations.csv')


def baseline(obs):
    print('\n--------Starting baseline measurement--------\n')
    virgo.observe(obs_parameters=obs, obs_file=storagePath + 'observationData/baseline/baseline.dat', spectrometer='wola', start_in = 0)
    print('\n--------Calibration finished--------\n')
    
    
def observe(obs, noPlot):
    print('\n--------Starting observation--------\n')
    # Begin data acquisition    
    filename = 'observation ' + str(Time.now()).replace(':', '_') + '.dat'
    fp = open(storagePath + 'observationData/observations/' + filename, 'x')
    fp.close()
    virgo.observe(obs_parameters=obs, obs_file=storagePath + 'observationData/observations/' + filename, spectrometer='wola', start_in = 0)
    print('\n--------Observation finished, plotting--------\n')
    if noPlot == False:
        analyzeData(filename, obs)
    

def analyzeData(filename, obs):
    recordedTime = filename.split(' ')[1] + ' '
    recordedTime += filename.split(' ')[2].replace('_', ':').split('.')[0] + '.'
    recordedTime += filename.split(' ')[2].replace('_', ':').split('.')[1]
    
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
    
    
    headerList = ['hz', 'spectrum', 'baseline', 'calibrated_spectrum']
    file = pd.read_csv(storagePath + 'csvFiles/spectra/spectrum ' + recordedTime.replace(':', '_') + '.csv') 
    file.to_csv(storagePath + "csvFiles/gfg2.csv", header=headerList, index=False)

    
    _temporaryData = pd.read_csv(storagePath + "csvFiles/gfg2.csv")
    calibrationSpectrum = _temporaryData[:]['calibrated_spectrum']
    spectrumAverage = moving_average(calibrationSpectrum, 20)
    
    peaks = find_peaks(spectrumAverage, prominence=0, width=10 , height = 1.14)
    
    raDec = oldTimeEquatorial(alt=90, az=1, lat=1, lon=1, height=0, time = recordedTime)
    galCoord = virgo.galactic(raDec[0], raDec[1])
    

    if peaks[0].size > 0:        
        writeFrame = pd.DataFrame()
        for i in range(0, len(peaks[0]) - 1):
            export = {'time': [recordedTime],
                      'ra': [raDec[0]],
                      'dec': [raDec[1]],
                      'l': [galCoord[0]],
                      'b': [galCoord[1]],
                      'hz': [_temporaryData[:]['hz'][peaks[0][i]]]}
            writeFrame = pd.DataFrame(export)
            print('\n')
            print(writeFrame)
            print('\n')
            writeFrame.to_csv(storagePath + 'allObservations.csv', mode='a', index=False, header=False)
    
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

if __name__ == "__main__":
    main()

