# Fortran code directory

This directory contains all the Fortran source code which eventually calculates the matrix element squared for each process.

## Directory structure

Subdirectory | Description
--- | ---
GridInstructions  |  Final instructions on how the phase space grid should be scanned for each process. Also contains the TGC definition when applicable.
Helpers  | Stand-alone helper functions used for the grid instructions.  
OMega  | The core O'Mega source code for the calculation of general matrix elements.  
Processes  | The O'Mega source code for the calculation of matrix elements for specific processes.


## Credit

The `O'Mega` source code is a copy-paste of the original `O'Mega` source code found here: https://whizard.hepforge.org/omega.html