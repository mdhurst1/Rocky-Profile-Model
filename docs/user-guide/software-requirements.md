---
title: 'Software Requirements'
---
RPM-CRN is written in C++ for efficiency. The code has been written and tested extensively in a Linux/UNIX environment, has also been compiled and run on Windows using Code::Blocks and Linux Subsystem for Windows, *and has been tested on Mac?*. 

Alternatively here is a docker environment???
Alternatively here is another online portal e.g. MyBinder???

Running the model will require working at a unix-style command line interface.

There are a number of software requirements to run the model and visualise the results:
 - A C++ compiler, preferably GCC: the GNU Compiler Collection
 - A text editor or integrated development environment (e.g. Notepad++, gedit, vims, VSCode)
 - A Python distribution with Scipy, Numpy and Matplotlib libraries

Optional for Multiobjective Optimisation:
 - Dakota 
 - QUESO

# Linux/UNIX
If you do not already work in Linux or UNIX, then the easiest way to get started would be to use some virtualisation software such as VirtualBox, and installing a Linux distribution such as Ubuntu on a virtual machine in VBox.

# Windows
Linux Subsystem for Windows

Alternatively, if you prefer to continue using Windows, it is possible to get the model working using the Code::Blocks IDE with MinGW (Minimalist GNU for Windows) compilers. The pair are available to install together here. We have not tested RPM-CRN extensively in this environment but the examples here will all compile and run correctly from Code::Blocks. 