import virgo
import pandas as pd
from moving_average import moving_average
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from datetime import datetime

# Define observation parameters
obs = {
    'dev_args': '',
    'rf_gain': 30,
    'if_gain': 25,
    'bb_gain': 18,
    'frequency': 1420e6,
    'bandwidth': 2.4e6,
    'channels': 2048,
    't_sample': 1,
    'duration': 60,
    'loc': '',
    'ra_dec': '',
    'az_alt': ''
}


# Check source position
#virgo.predict(lat=39.8, lon=-74.9, source='Cas A', date='2020-12-26')


# Begin data acquisition
#   virgo.observe(obs_parameters=obs, obs_file='observation.dat')


#virgo.simulate(1, 1, beamwidth=16, v_min = -800, v_max = 800, plot_file='D:/SDR/plot3')

#print(virgo.gain(0.7, 1420e6, 16, u = 'dbi'))
#print(virgo.beamwidth(1, 1420e6))
#virgo.map_hi(ra = None, dec = None, plot_file='D:/SDR/map1')



#virgo.plot(obs_parameters=obs, n=20, m=35, f_rest=1420.4057517667e6,
#           vlsr=False, meta=False, avg_ylim=(-5,15), cal_ylim=(-20,260),
#           obs_file='observation.dat', cal_file='calibration.dat', rfi=[(1419.2e6, 1419.3e6), (1420.8e6, 1420.9e6)],
#           dB=True, spectra_csv='spectrum.csv', plot_file='plot4.png')

#Following code is from PICTOR telescope
#virgo.plot(obs_parameters=obs, dB =  True , n=20, m=20, f_rest=1420.4057517667e6, obs_file='observation.dat', cal_file='calibration.dat', spectra_csv='spectrum.csv', power_csv='time_series.csv', plot_file='pictorPlot.png')
        

now = datetime.now()

headerList = ['hz', 'spectrum', 'calibration', 'calibrated_spectrum']
file = pd.read_csv(r'./csvFiles/spectrum.csv') 
file.to_csv("./csvFiles/gfg2.csv", header=headerList, index=False)


_temporaryData = pd.read_csv("./csvFiles/gfg2.csv")
calibrationSpectrum = _temporaryData[:]['calibrated_spectrum']
spectrumAverage = moving_average(calibrationSpectrum, 20)

peaks = find_peaks(spectrumAverage, prominence=0, width=10 , height = 1.14)

raDec = virgo.equatorial(alt=90, az=1, lat=1, lon=1, height=0)
galCoord = virgo.galactic(raDec[0], raDec[1])

for i in range(0, len(peaks[0])):
    export = {'time': datetime.now(),
              'ra': raDec[0],
              'dec': raDec[1],
              'l': galCoord[0],
              'b': galCoord[1],
              'hz': _temporaryData[:]['hz'][peaks[0]]}

writeFrame = pd.DataFrame(export)
writeFrame.to_csv('test.csv', mode='a', index=False, header=False)

plt.figure(dpi = 250)
plt.plot(_temporaryData[:]['hz'], spectrumAverage)
plt.plot(_temporaryData[:]['hz'][peaks[0]], spectrumAverage[peaks[0]], marker = 7, color = 'r')
plt.title('Spectrum taken at ' + str(now).replace(':', '-').replace('.', '-'))
plt.ylabel('Relative Power (1)')
plt.xlabel('Frequency (hz)')

plt.savefig('./imageExport/spectrum ' + str(now).replace(':', '-').replace('.', '-') + '.png')
