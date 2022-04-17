# Hydrogen Line Emission Software
This software is designed to handle data aquisicion and analysis for a hydrogen line radio telescope. 

The specific goal with this application is to create a program that allows obvservation of hydrogen line emissions at differnet times with possible interruptions. Frequency peaks in each spectrum are identified and saved with thei respective galactic coordinates in a CSV file. This allows interruption of data aquisition without loosing observed data 

Foe each observation the waterfall and spectra are saved using VIRGO, peaks in the spectrum are then identified using pandas.

# Dependencies
This project is based on [VIRGO](https://github.com/0xCoto/VIRGO), therefore it shares its dependencies. The moving average filter is taken from [motorr4ik](https://github.com/motorrr4ik/moving_average_filters)
I had to modify some parts of virgo.py, you can find the modified version in this repository. 

# Installation
Clone the github repo, and in the console run

```
python main.py -sC -nO
```

This will generate the default configuration file in the repo. The following command will set up the directories necessary. You can change the folder path in the configuration file. 
Generate the directories with

```
python main.py -sF -nO
```


# Usage
The program will read the execution parameters from the configuration file. These can be overridden with the following arguments when executing the program:

```
-r, --runs: Number of measurements performed
-d, --duration_: Time per observation (s)
-t, --timeSample: Time per sample (s)
-nO, --noObservation: Disables automatic data acquisition
-b, --baseline: Measure baseline. Automatic -nO
-sF, --setupFolders: Generates default folder structure
-sC, --setupConfig: Generates default configuration file
-pA, --plotAll: Plot all observation files. Automatic -nO
-nP, --noPlot: Disables automatic plotting
```

Plotting data requires a baseline to be captured, otherwise the program will crash. You can still record data without plotting by using the -nP flag. 
