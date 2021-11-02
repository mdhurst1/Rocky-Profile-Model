Topography
===
Example topographic profile output file: *example_ShoreProfile.xz*

Output of an across-shore topographic profile is periodically written to file (default is every 10 years controlled by the `PrintInterval` parameter in the parameter file). The coastal profile erodes from  the left to the right hand side of the array (so sea on the left, land on the right). The horizontal (`X`) dimension of the cell array is coded to grow dynamically aas the model simulation proceeds. The topographic profile is output as a vector of the horizontal position of cells at each elevation (`Z`) between the maximum and minimum elevations, spaced at the vertical resolution of the model domain `dZ` (default is 0.1 m).

The first line of the output file contains the maximum elevation, minimum elevation, and resolution of cell size (in metres).

Each row of shore profile output file contains: time (years), elevation (metres) of sea level at that time, followed by a delimited vector of X positions (in metres) of rock cells for every elevation from maximium to minimum, at the spacing specified in the header.

| Z<sub>max</sub> | Z<sub>min</sub> | Z<sub>res</sub> | | | | |
|---|---|---|---|---|---|---|
| T<sub>0</sub> | SL<sub>0</sub> | Z<sub>0</sub> |  Z<sub>1</sub> |  Z<sub>2</sub> | ... | Z<sub>n</sub> |
| ... | ... |Z<sub>0</sub> | Z<sub>1</sub> |  Z<sub>2</sub> | ... | Z<sub>n</sub> |
| T<sub>T</sub> | SL<sub>T</sub> | Z<sub>0</sub> | Z<sub>1</sub> |  Z<sub>2</sub> | ... | Z<sub>n</sub> |

CRN concentrations
===

Example CRN concentration profile output file: *example_Concentrations.xn*

The CRN concentration output corresponds to the topographic profiles and are also periodically written to file at the same intervals as the topography file above (every 100 years by default, controlled by `PrintInterval` in the parameter file). The output differs from the topography output in that a horizontal array of cells is written to file; with each cell in the `X` dimension assigned a value corresponding to the CRN concentration at the rock surface. The horizontal vector of CRN concentrations is read from left to right, with the cliff position at the right side of the vector. 

The first line of the file contains a delimited list of the nuclides being tracked (e.g. 10: <sup>10</sup>Be; 14: <sup>14</sup>C; 26: <sup>26</sup>Al) 

Each row of CRN concentration profile output file: time (years BP), the nuclide being tracked (10: <sup>10</sup>Be etc.), followed by a vector of CRN concentrations (`N`) of surface rock cells (atoms g<sup>-1</sup>) for every `X` position across shore platform, from left to right. 

| 10 | 14 | 26 | | | | |
|---|---|---|---|---|---|---|
| T<sub>0</sub> | 10 | N<sub>0</sub> | N<sub>1</sub> | N<sub>2</sub> | ... | N<sub>n</sub> |
| T<sub>0</sub> | 14 | N<sub>0</sub> | N<sub>1</sub> | N<sub>2</sub> | ... | N<sub>n</sub> |
| T<sub>0</sub> | 26 | N<sub>0</sub> | N<sub>1</sub> | N<sub>2</sub> | ... | N<sub>n</sub> |
| ... | ... | ... | ... | ... | ... | ... |
| T<sub>T</sub> | 10 | N<sub>0</sub> | N<sub>1</sub> | N<sub>2</sub> | ... | N<sub>n</sub> |
| T<sub>T</sub> | 14 | N<sub>0</sub> | N<sub>1</sub> | N<sub>2</sub> | ... | N<sub>n</sub> |
| T<sub>T</sub> | 16 | N<sub>0</sub> | N<sub>1</sub> | N<sub>2</sub> | ... | N<sub>n</sub> |