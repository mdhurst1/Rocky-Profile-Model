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