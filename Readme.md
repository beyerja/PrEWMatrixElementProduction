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

