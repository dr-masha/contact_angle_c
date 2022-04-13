Author: Masa Prodanovic, all rights reserved.
Note: Original code is from 2005 and thus some dependenices are old. Minor changes made in 2022 (compilation).


Prerequisites for contact angle code are:
0. Glib 
1. The GNU Triangulated Surface Library
2. Geomview (for visualizing the result)
Both are pretty old Unix based codes, but still appear to be able to compile.  


WINDOWS OR LINUX INSTALLATION NOTES

If you have Windows 10, please activate Ubuntu Linux subsystem.
The installation instructions below assume you have opened Ubuntu Terminal and placed the code in desired
directory, e.g. in your home directory, create subdirectory local/ for all of the local code installation.

I will refer to your home directory as /home/USERNAME

cd

# Create local subdirectories
mkdir local
cd local 
mkdir src bin lib

cd src
# Now copy the contact angle code here

#Installing GLIB will happen with this:
sudo apt update
sudo apt install libglib2.0-dev

cd gts-0.7.6
# Note: /home/USERNAME/local has bin/, src/ and lib/ sub-directories
# Note that on MAC, this directory is /Users/USERNAME/local
./configure --prefix=/home/USERNAME/local
make
make install

cd contact_angle_c
# open makefile and make sure all directories are pointing to correct locations
make
make install

# Install geomview; this also requires setting up visual display

# Install VcXsrv on Windows
# Run XLaunch on Windows
# (not sure if this is a must) In Turn Windows features ON or OFF, turn on virtual Machine Platform
# I assume here 
echo "export DISPLAY=localhost:0.0" > ~/.bashrc
source ~/.bashrc
sudo apt-get install geomview

MAC INSTALLATION NOTES

Prerequities: 
Install Xcode, MacPorts and Xquartz. Otherwise, selected components below will not be properly installed or possible to use (e.g. Geomview).

# install glib2
sudo port install glib2-devel

# need pkg-config
sudo port install pkgconfig

#install Geomview
sudo port install geomview 

#follow above instructions to install gts - GNU Triangulated Surface Library (download from sourceforgenet, 
