# Hydrogen Line Emission Software
This software is designed to handle data aquisicion and analysis for a hydrogen line radio telescope. 

The specific goal with this application is to create a program that allows sperad-out obvservation of hydrogen line emissions with possible interruptions. Frequency peaks in each spectrum are identified and saved with thei respective galactic coordinates in a CSV file. This allows interruption of data aquisition without loosing observed data 

Foe each observation the waterfall and spectra are saved using VIRGO, peaks in the spectrum are then identified using pandas.

# Dependencies
This project is based on [VIRGO](https://github.com/0xCoto/VIRGO), therefore it shares its dependencies.
