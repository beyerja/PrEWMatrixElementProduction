# PrEWMatrixElementCalculator

Code to produce matrix element level input for the PrEW fitting framework.

## Usage

The `O'Mega` fortran code (`src/fortran`) can be compile using a prepared macro:

```bash
cd macros
chmod u+x create_executables.sh
./create_executables.sh
```

This will create the binary executables for each process (in `bin`). The binaries can be called with the given bin-indices to extract the matrix element squared at that bin-index.

## Disclaimers

- The current code is specified to use a center-of-mass energy of 250GeV. 
  This is hardcoded in parts of the O'Mega code and probably also in some ROOT scripts. 
  Sorry about that, didn't have time yet to check this.
- Most of the fortran and ROOT code has been copy-pasted blind from Robert Karls work.
  I did not clean it up and did not check most of it. 
  Good luck.