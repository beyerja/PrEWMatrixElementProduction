#!/bin/bash

# Create the binary executables needed to calculate the matrix element 
# distributions with O'Mega


# ---------------------------------------------------------------------------- #
# Establish the important directories

# Path to this current script
dir="$( cd "$( dirname "${BASH_SOURCE[0]}"  )" && pwd  )" 

src_dir="${dir}/../src/fortran"
omega_dir="${src_dir}/OMega"
helper_dir="${src_dir}/Helpers"
processes_dir="${src_dir}/Processes"

# Output directories
bin_dir="${dir}/../bin"
lib_dir="${dir}/../lib"
if [[ ! -d ${bin_dir} ]] ; then # Create if not existing
  mkdir -p ${bin_dir}
fi
if [[ ! -d ${lib_dir} ]] ; then # Create if not existing
  mkdir -p ${lib_dir}
fi

# ---------------------------------------------------------------------------- #
# Find all the process files
process_files=$(find ${processes_dir} -iname "*.f90" -print)

# ---------------------------------------------------------------------------- #
# Fortran base files need to be listed in correct order to compile
omega_files=(
omega_kinds.f90
omega_constants.f90
omega_parameters.f90
omega_spinors.f90
omega_vectors.f90
omega_polarizations.f90
omega_tensors.f90
omega_tensor_polarizations.f90
omega_couplings.f90
omega_spinor_couplings.f90
omega_utils.f90
omega95.f90
ww_sl0eq.f90
ww_sl0muq.f90
)

helper_files=(
MathOperations.f90
)

# ---------------------------------------------------------------------------- #
# Compile all needed source files into libraries 

echo "Start compiling O'Mega."
cd ${lib_dir}
for omega_file in ${omega_files[@]};
do
  echo "Compiling into library: ${omega_file}"
  gfortran -c ${omega_dir}/${omega_file}
done
echo "Done compiling O'Mega."

echo "Start compiling Helpers."
cd ${lib_dir}
for helper_file in ${helper_files[@]};
do
  echo "Compiling into library: ${helper_file}"
  gfortran -c ${helper_dir}/${helper_file}
done
echo "Done compiling Helpers."

libaries=$(find ${lib_dir} -iname "*.o" -print)
ar -qs "${lib_dir}/libomega.a" $libaries

# ---------------------------------------------------------------------------- #

echo "Start compiling the processes."
cd ${bin_dir}

for process_file in ${process_files[@]};
do
  echo "Compiling into binary: ${process_file}"
  output_name_base=$(basename --suffix=.f90 ${process_file})
  gfortran ${process_file} -J${lib_dir} -I${lib_dir} -L${lib_dir} -lomega -o "${output_name_base}.bin"
done
echo "Done compiling processes."

# ---------------------------------------------------------------------------- #
# Clean up

echo "Cleaning up intermediate libraries."
rm ${lib_dir}/*

# ---------------------------------------------------------------------------- #

echo "Done!"