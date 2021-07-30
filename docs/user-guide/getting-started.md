---
title: 'Getting Started'
---

Compiling the code
===

Model runs using RPM-CRN are handled by **driver files** that initiate objects and run the simulations. These driver files are C++ scripts that are compiled to form an executable program. The default driver generates a standalone model that ingests a parameter file that controls the software by setting values for the various control parameters in the model.

To compile the model on your local machine, use the `make` command and point it to the makefile in the repository:
```
Rocky-Profile-Model$ make -f RPM_CRN.make
```
This will result in an executable called `RPM_CRN.exe`. 

Running the code
===

The program can then be launched at the command line, requiring two input arguments.
* The path of the folder where the model will be run
* The name of the input parameter file (which must be in the folder where the model will be run)

The following command will launch the executable in its current directory with the default parameter values:
```
Rocky-Profile-Model$ RPM_CRN.exe ./ example_parameter_file.txt
```
Alternatively you could move the parameter file to a directory for your project and customise it:
```
Rocky-Profile-Model$ RPM_CRN.exe /home/mhurst/MyFirstRPMCRN/ modified_parameter_file.txt
```
For a detailed description of the contents of a parameter file, see [configuration options](configuration-options.md).

Plotting the results
===

The model output will be written to your project directory. There are two output files, one containing the timeseries of topographic evolution of the simulated rock coast, and once containing the corresponding timeseries of CRN concentrations at the surface of the topography.

For more information about the format of the output files, see [model output](model-output.md).

Visualisation and further analysis of the model run is conducted in python (though you could write your own functions in a different language such as R or Matlab :face_vomiting: :vomiting_face:). A timeseries of the evolution of the model topography and snapshhot of the final CRN concentrations can be generated using `/plotting_functions/plot_model_output.py`.

IMAGE OF OUTPUT GOES HERE