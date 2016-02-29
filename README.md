# gfncalc - Correlation function calculator

## SHORT DESCRIPTION:

As the title suggests, gfncalc calculates normalized correlation functions from photon arrival times. In physics, these functions are normally written as g^(n)(t), hence the acronym. At present, it only computes the second order. There are plans to compute third order, but this has not been implemented yet (There are some commented parts of the code that show initial work towards this goal). Currently it also reads trigger pulses from another device to implement filtering. It features a simple command-line interface and it is cross-platform. Uses a configuration file (during first run it will create one -gfncalc.conf-) and remembers last parameters used. This program is fast at processing large datasets.

## LICENSE:

This program is released under the terms of the GNU Public License v3 (GPLv3). For more information please see:
http://www.gnu.org/licenses/gpl-3.0.html

## HOW TO USE:

At present, the program can only handle input ascii text files which contain two columns of data. First column displaying the time and the second one showing the channel (or detector) on which the photon impinged.

The program requires several inputs to run successfully. It has a hardcoded detector resolution time of 164.61 picoseconds. This can be changed but it requires recompilation. The next parameter is the path of the data (an example is provided in the CLI). The channels of all the detectors and trigger pulses must be specified (APD C is set to zero, which disables it). Binning and window size are required, they are input in multiples of the detector resolution. Then there is the internal pulse delay. This is the time, in multiples of detector resolution time, that it takes for a triggering pulse to arrive after a click on channel APD A arrives. The delay tolerance allows for some margin of error for the internal pulse delay. Finally, there is the rising edge delay. This can be input in microseconds. If our trigger pulses have a preset delay in the device, then this is where we can tell the program about this. All the pulse parameters can have random values if one only wishes to calculate second-order correlations without the trigger pulses.

It is designed to accept input of the form dataX.asc. Where X can be any positive integer. In this way, it can be run for multiple datasets.

Its output is also a text file. As a header, it displays the dataset number, total number of photon "clicks" on each channel (or detector - APD) as well number of trigger pulses and rates. The rest of the output is organized in four columns. First one is time in microseconds. Second is the number of photons in the time-bin. Third is the normalized second order correlation amplitude and finally the error on this amplitude, which is basically propagated from the poissonian error in the number of photons.

To run the program on a Unix-based system, open a terminal and then type:

./gfncalc

In Windows, in a command prompt, simply type

gfncalc.exe

## COMPILATION:

It is written in plain C++ using only the standard libraries. We have only thoroughly tested g++ (The GNU C++ compiler) as compiler, on a GNU/Linux 64 bits machine. With g++ and minGW we cross-compile the program to run on Windows (We did do an old test with MS Visual Studio 2010 for an older version of the code and it compiled without a problem). We have not tested Mac OSX, but given g++ is easily available for this platform, it should be easy to produce an executable for it.

### FOR GNU/LINUX (AND POTENTIALLY MACOSX):

In our GNU/Linux machine, to compile we simply do:

g++ -Wall gfncalc.cxx read_inputs.cxx -o gfncalc
(Only if the read_inputs.cxx and read_inputs.h files are located in the same directory)

or more generally

g++ -Wall -I/path-to-header-file gfncalc.cxx /path-to-file/read_inputs.cxx -o gfncalc

This binary has been tested on Debian GNU/Linux systems version 5/6/7 and CentOS 5. However, it should work on any Unix-based OS.

### FOR WINDOWS (USING GNU/LINUX TO COMPILE):

To cross-compile using the GNU/Linux machine to obtain a working Windows 32-bit executable we do this (g++-minGW 4.6.3):

i686-w64-mingw32-g++ -Wall -static -I/path-to-header-file gfncalc.cxx /path-to-file/read_inputs.cxx -o gfncalc.exe

The executable has been tested and works well in Windows XP/Vista/7. However, performance of the program is inferior to the GNU/Linux binary (even when it was compiled with MS Visual Studio).
