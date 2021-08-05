---
title: 'Software Requirements'
---
RPM-CRN is written in C++ for efficiency. The code has been written and tested extensively in a Linux/UNIX environment, has also been compiled and run on Windows using Code::Blocks and Linux Subsystem for Windows, *but has not been extensively tested on Mac*. 

Running the model will require working at a unix-style command line interface.

There are a number of software requirements to run the model and visualise the results:
 - A C++ compiler, preferably GCC: the GNU Compiler Collection
 - A text editor or integrated development environment (e.g. Notepad++, gedit, vims, VSCode)
 - A Python distribution with Scipy, Numpy and Matplotlib libraries

Optional for Multiobjective Optimisation:
 - Dakota 
 - QUESO

# Docker

Docker is software that allows you to run "containers" within your own operating system that have all of the required tools, libraries  and software needed to run RPM-CRN. The [Docker website](https://www.docker.com/resources/what-container) has the following definition: "A container is a standard unit of software that packages up code and all its dependencies so the application runs quickly and reliably from one computing environment to another."

Download and install Docker for your OS platform from the [Docker website](https://www.docker.com/products/docker-desktop).

# Windows
This code has been developed both in Linux and Windows environments. When working in Windows, this usually involves setting up a Linux virtual environment, for which there are a couple of options.

## Linux on a Virtual Machine
If you do not already work in Linux or UNIX, then the easiest way to get started would be to use some virtualisation software such as [VirtualBox](https://www.virtualbox.org/wiki/Downloads), and installing a Linux debian-based distribution such as [Ubuntu](https://ubuntu.com/) or [Linux Mint](https://linuxmint.com/) as a virtual machine in VBox. There is a nice video explaining how to do this [here](https://www.youtube.com/watch?v=x5MhydijWmc).

## Windows Subsystem for Linux
Alternatively, you can use the Windows Subsystem for Linux (WSL) Instructions for installing Linux Subsystem for Windows can be found [https://docs.microsoft.com/en-us/windows/wsl/install-win10](here).

We find that WSL integrates really nicely within [Visual Studio Code](https://code.visualstudio.com/), our preferred integrated development environment (IDE).

## Code::Blocks
Finally, if you prefer to continue using Windows and do not wish to dive into Linux, it is possible to get the model working using the [https://www.codeblocks.org/](Code::Blocks) IDE with MinGW (Minimalist GNU for Windows) compilers. The pair are available to install together here. We have not tested RPM-CRN extensively in this environment but the examples here will all compile and run correctly from Code::Blocks. 

# Linux/UNIX
![](./img/terminal.png =100x) Once you have access to a Linux machine or Linux virtual machine (WSL or VBox as above), we will work in a terminal. 

## Download code release
In a terminal, to get the code, you can download a copy:
`
HOME$ wget ...
`

## Clone the repository
Alternatively, you can clone the repository using git (particularly if you wish to do any development work). To install git (if you do not already have it: )

`
HOME$ sudo apt install git-all`
`
The download a clone of the repository with:
`
HOME$ git clone ...
`

If your using VSCode, the git integration is fantastic and you can clone the repository this way too by going to `View -> Comman Pallete` and typing `Git: Clone`. It will ask for the repository URL as above.

# MAC

Reconsider your life choices. Apparently you can install VirtualBox on a MAC. So there's that.