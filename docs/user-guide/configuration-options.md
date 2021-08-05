---
# Model Parameters
Model parameters are set in a parameter file that is passed to the model as an argument and the model reads when the software is initiated.
Detailed explanation of all model parameters goes here.

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
|
|   | *Spatial Domain* |  |
|
|   | *Time Control* |  |
|

