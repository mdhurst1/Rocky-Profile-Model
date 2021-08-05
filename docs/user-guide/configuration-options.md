---
# Model Parameters
Model parameters are set in a parameter file that is passed to the model as an argument and the model reads when the software is initiated:

| *Parameter* | *Description* | *Default Value* |
| ----------- | ------------- | --------------- |
|   | *File Structure* |  |
| Folder | Path fo working directory (relative or absolute) | ./RPM_CRN/ |
| ProjectName | Descriptor for the simulations, prefix for output filenames | TestProject |
|   | *Cosmogenic Radionuclides* |  |
| CRN_Predictions | Bool controlling whether or not CRN predictions are made | true |
| Berylium | Bool controlling whether <sup>10<\sup>Be is tracked | true |
| Carbon | Bool controlling whether <sup>14<\sup>C is tracked | true |
| Aluminium | Bool controlling whether <sup>26<\sup>Al is tracked | true |
|   | *Hydrodynamics* |  |
| SeaLevelFilename | Text file containing a timeseries of relative sea level | NULL |
| SeaLevelRise | Constant rate of relative sea level change, positive values indicate sea level rise (in metres yr<sup>-1<\sup>) | 0.001 |
| TidalRange | Total spring tide range from lowest to highest tide (in metres) | 1.5 |
| TidalPeriod | Time (in hours) between consecutive peak tides | 12.25 |
| WaveHeight_Mean | Mean wave height (in metres) | 1. |
| WaveHeight_StD | Standard deviation of wave heights (in metres), for sampling about a normal distribution centred on WaveHeight_Mean |  0. |
| WavePeriod_Mean | Mean wave period (in seconds) | 6. |
| WavePeriod_StD | Standard deviation of wave period (in seconds), for sampling about a normal distribution centred on WavePeriod_Mean | 0. |
| StandingWaveCoef | Coefficient for energy expenditure of a standing wave | 0.1 |
| BreakingWaveCoef | Coefficient for energy expenditure of a breaking wave | 10. |
| BrokenWaveCoef | Coefficient for energy expenditure of a broken wave | 0.01 |
| WaveAttenuationConst | Constant describing the rate of decay of wave height/energy with distance across the intertidal zone after wave breaking | 0.01 |
|   | *Geology* |  |
| Resistance | Force required to erode an unweathered cell of rock (UNTIS!) | 10. |
| WeatheringRate | Maximum rate at which resistance is reduced for exposed cells due to intertidal weathering (UNITS) | 1. |
| SubtidalEfficacy | Constant rate of subtidal weathering for exposed submarine cells as a fraction of the maximum intertidal weathering rate `WeatheringRate` above | 0.001 |
| CliffFailureDepth | Amount of undercutting required to prompt cliff failure (i.e. the removal of all overlying cells). Provides crude approach to implementing episodic cliff failure | 0.1 |
| CliffElevation | Height of the cliff (in metres) Need some info about moving reference frame here | 30 |
|   | *Spatial Domain* |  |
| MinElevation | Minimum elevation in the model's spatial domain (in metres) | -20. |
| MaxElevation | Maximum elevation in the model's spatial domain (in metres) | 10 |
|   | *Time Control* |  |
| StartTime | Time to commence model simulation (in years) | 8000 |
| EndTime | Time at which to stop the simulation (in years) | 0 |
| TimeStep | Length of a single model timestep (in years) | -1 |
| PrintInterval | Periodicity of printing model output to file (in years) | 10 |

 ## - Example paramter file

 Parameters are parsed from a text file by the `parameters.cpp` object

 ```
 # parameter file for RPM_CRN model

# File Structure
Folder: ./RPM_CRN/
ProjectName: TestProject

# Cosmogenic Isotopes
CRN_Predictions: true
Berylium: true
Carbon: false
Aluminium: false

# Hydrodynamics
SeaLevelFilename: NULL
SeaLevelRise: 0.001
TidalRange: 1.5
TidalPeriod: 12.25
WaveHeight_Mean: 1.
WaveHeight_StD: 0.
WavePeriod_Mean: 6. 
WavePeriod_StD: 0.
StandingWaveCoef: 0.1
BreakingWaveCoef: 10.
BrokenWaveCoef: 0.01
WaveAttenuationConst: 0.01

# geology
Resistance: 10.
WeatheringRate: 1.
SubtidalEfficacy: 0.001
CliffFailureDepth: 0.1
CliffElevation: 30

# space domain
MinElevation: -20
MaxElevation: 10

# time control
StartTime: 8000
EndTime: 0
TimeStep: -1
PrintInterval: 10

```