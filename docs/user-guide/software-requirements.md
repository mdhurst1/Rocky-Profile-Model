---
title: 'Software Requirements'
---
<h1>Software Requirements</h1>

RPM-CRN is written in C++ for efficiency. The code has been written and tested extensively in a Linux/UNIX environment, has also been compiled and run on Windows using Code::Blocks and Linux Subsystem for Windows, *but has not been extensively tested on Mac*. 

Running the model will require working at a unix/linux style command line interface.

There are a number of software requirements to run the model and visualise the results:
 - A C++ compiler, preferably GCC: the GNU Compiler Collection
 - The **make** utility to compile the model with our **make files**
 - A text editor or integrated development environment (e.g. Notepad++, gedit, vims, VSCode)
 - A Python distribution with Scipy, Numpy and Matplotlib libraries

Optional for Multiobjective Optimisation:
 - Dakota 
 - QUESO

## MAC
Reconsider your life choices. Apparently you can install VirtualBox on a MAC. So there's that.

## Windows
This code has been developed both in Linux and Windows environments. When working in Windows, this usually involves setting up a Linux virtual environment, for which there are a couple of options.

### - Linux on a Virtual Machine
If you do not already work in Linux or UNIX, then the easiest way to get started would be to use some virtualisation software such as [VirtualBox](https://www.virtualbox.org/wiki/Downloads), and installing a Linux debian-based distribution such as [Ubuntu](https://ubuntu.com/) or [Linux Mint](https://linuxmint.com/) as a virtual machine in VBox. There is a nice video explaining how to do this [here](https://www.youtube.com/watch?v=x5MhydijWmc).

### - Windows Subsystem for Linux
Alternatively, you can use the Windows Subsystem for Linux (WSL) Instructions for installing Linux Subsystem for Windows can be found [here](https://docs.microsoft.com/en-us/windows/wsl/install-win10).

We find that WSL integrates really nicely within [Visual Studio Code](https://code.visualstudio.com/), our preferred integrated development environment (IDE).

### - Code::Blocks
Finally, if you prefer to continue using Windows and do not wish to dive into Linux, it is possible to get the model working using the [Code::Blocks](https://www.codeblocks.org/) IDE with MinGW (Minimalist GNU for Windows) compilers. The pair are available to install together here. We have not tested RPM-CRN extensively in this environment but the examples here will all compile and run correctly from Code::Blocks. 

## Linux/UNIX
![Terminal](img/terminal.png "Terminal")

Once you have access to a Linux machine or Linux virtual machine (WSL or VBox as above), we will work in a terminal. 

### - Download code release
In a terminal, to get the code, you can download a copy:
```
HOME$ wget https://github.com/mdhurst1/Rocky-Profile-Model/archive/refs/heads/JOSS-Paper.zip
```

### - Clone the repository
Alternatively, you can clone the repository using git (particularly if you wish to do any development work). To install git (if you do not already have it: )

```
HOME$ sudo apt install git-all
```
The download a clone of the repository with:
```
HOME$ git clone https://github.com/mdhurst1/Rocky-Profile-Model.git
```

If your using VSCode, the git integration is fantastic and you can clone the repository this way too by going to `View -> Comman Pallete` and typing `Git: Clone`. It will ask for the repository URL as above.

## Dakota with QUESO

The Dakota toolkit delivers flexible and extendable software for model optimisation and uncertainty quantification. The software package [QUESO](https://github.com/libqueso/queso) (Quantification of Uncertainty for Estimation, Simulation and Optimisation) for the solution of inverse statistical problems is wrapped within Dakota. RPM-CRN works with Dakota and QUESO to optimise model parameters to match observed shore platform topography and measured CRN concentrations, and quantify uncertainty, in order to derive a liely timeseries of rock coast development.

The following [shell script](scripts/Dakota_Queso_Install.sh) contains the terminal commands to install Dakota, Queso and associated dependencies. Thanks to Geraldo Fiorini Neto for collating the shell script.

To run the script, place it in your $home directory, then to make it executable, run:
```
HOME$ chmod +x Dakota_Queso_Install.sh
```
Run the script with:
```
HOME$ ./Dakota_Queso_Install.sh
```

```
#!/bin/bash

## Dakota 6.11.0 with QUESO 0.57.1

# install lib dependencies
apt-get install -y gcc g++ gfortran cmake libboost-all-dev libblas-dev liblapack-dev libopenmpi-dev openmpi-bin gsl-bin libgsl-dev python perl cmake-curses-gui libboost-dev openmpi-doc xorg-dev libmotif-dev build-essential gfortran autotools-dev curl unzip git vim openssh-client openssh-server libgsl0-dev python-pip doxygen texlive-latex-extra

# Install python modules
pip install numpy scipy pandas mpmath==1.0.0 sympy==1.2 pyyaml

# Remove old version of Java and install newer
apt-get remove default-jre && apt-get install openjdk-8-jre default-jdk


# Download and compile Queso 0.57.1

mkdir -p /root/queso && cd /root/queso
wget https://github.com/libqueso/queso/releases/download/v0.57.1/queso-0.57.1.tar.gz
tar -zxvf queso-0.57.1.tar.gz
cd queso-0.57.1
./configure --prefix=/opt/queso/0.57.1/
make && make install 

# Download and extract Dakota 6.11.0

mkdir -p /root/dakota && cd /root/dakota
wget https://dakota.sandia.gov/sites/default/files/distributions/public/dakota-6.11.0-release-public.src-UI.tar.gz
tar -zxvf dakota-6.11.0-release-public.src-UI.tar.gz
mkdir -p /root/dakota/build

# Set environment variables and run

set ENABLE_DAKOTA_DOCS:BOOL=TRUE

echo "export PATH=\${PATH}:/opt/dakota/6.11.0/bin" >> /etc/profile.d/dakota.sh
echo "export LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:/opt/dakota/6.11.0/lib" >> /etc/profile.d/dakota.sh
. /etc/profile.d/dakota.sh

echo "export PATH=\${PATH}:/opt/queso/0.57.1/bin" >> /etc/profile.d/queso.sh
echo "export LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:/opt/queso/0.57.1/lib" >> /etc/profile.d/queso.sh
. /etc/profile.d/queso.sh

# Configure Dakota

cd /root/dakota/build
cmake  -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_C_FLAGS:STRING="-O2" -DCMAKE_CXX_FLAGS:STRING="-O2" -DCMAKE_Fortran_FLAGS:STRING="-O2" -DDAKOTA_HAVE_MPI:BOOL=TRUE -DDAKOTA_HAVE_GSL:BOOL=TRUE -DHAVE_QUESO:BOOL=TRUE  -DBoost_NO_BOOST_CMAKE:BOOL=TRUE -DCMAKE_INSTALL_PREFIX=/opt/dakota/6.11.0/ -DENABLE_DAKOTA_DOCS:BOOL=TRUE  ../dakota-6.11.0-release-public.src-UI/

# Compile and install Dakota

make -j 4 && make install 
```