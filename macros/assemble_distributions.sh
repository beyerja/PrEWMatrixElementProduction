#!/bin/bash

#!/bin/bash

# ------------------------------------------------------------------------------
# Input parameters to be defined
process=false
output_base=false

# ------------------------------------------------------------------------------
# Read input
for i in "$@"
do
case $i in
  --process=*)
    process="${i#*=}"
    shift # past argument=value
  ;;
  --output-base=*)
    output_base="${i#*=}"
    shift # past argument=value
  ;;
  -h|--help)
    echo "usage: ./assemble_distributions.sh --process=[ww_sl0muq/...] --output-base=[...]"
    shift # past argument=value
  ;;
  *)
    # unknown option
    shift
  ;;
esac
done

# ------------------------------------------------------------------------------
# Check input arguments:

if [ "$process" = false ]; then
  >&2 echo "No process given!"; exit 1
fi

if [ "$output_base" = false ]; then
  >&2 echo "No output directory given!"; exit 1
fi

# ------------------------------------------------------------------------------
# Check if ROOT is loaded

if ! command -v root &> /dev/null
then
  >&2 echo "ROOT could not be found, please load it."; exit 1
fi

# ------------------------------------------------------------------------------
# Script environment
dir="$( cd "$( dirname "${BASH_SOURCE[0]}"  )" && pwd  )" # path of this macro

# ROOT scripts directory
scripts_dir=${dir}/../src/cpp

# Output directory
output_dir=${output_base}/${process}
if [[ ! -d ${output_dir} ]] ; then # Create if not existing
  >&2 echo "No directory for process ${process} found: ${output_dir}"; exit 1
fi

# Temporary ROOT script environment
tmp_dir=${output_dir}/tmp
if [[ ! -d ${tmp_dir} ]] ; then # Create if not existing
  mkdir -p ${tmp_dir}
fi

# Other needed output directories
if [[ ! -d ${output_base}/distributions/angular ]] ; then
  mkdir -p ${output_base}/distributions/angular
fi
if [[ ! -d ${output_base}/distributions/TGC ]] ; then
  mkdir -p ${output_base}/distributions/TGC
fi
if [[ ! -d ${output_base}/distributions/combined ]] ; then
  mkdir -p ${output_base}/distributions/combined
fi

# ------------------------------------------------------------------------------
# Core functionality: 
# - Copy the scripts into a dedicated directy
# - Replace the necessary tokens 
# - Run the scripts
# - Clean up

# Find all scripts
scripts=$(find ${scripts_dir} -iname "*.C" -print)

# Function that takes the script, searches for a given token, and replaces the 
# line of the token with the given replacement
replace_token() {
  sed -i "s|${token}|${replacement}|g" ${script}
}

for script_template in ${scripts[@]};
do
  # Make a copy of the script
  script_name=$(basename ${script_template})
  echo "Running setting up script: ${script_name}"
  script="${tmp_dir}/${script_name}"
  cp ${script_template} ${script}

  # Replace all relevant tokens
  token=PROCESS_MARKER
  replacement=$process
  replace_token
  
  token=OUTPUT_DIR_MARKER
  replacement=$output_dir
  replace_token
  
  token=N_BINS_MARKER # Number of bins is process dependent!
  if [[ ${process} == "ww_sl0muq" ]]; then
    replacement=2000
  else 
    >&2 echo "Process ${process} unknown!"; exit 1
  fi
  replace_token
done

# Run the scripts in the right order
echo "Running ROOT scripts."
echo "Convert O'Mega grid files to ROOT TTree."
root -b -q ${tmp_dir}/GridAsciiToTTree.C
echo "Calculate necessary angular TGC coefficients from ROOT file."
root -b -q ${tmp_dir}/CreateAngularDistributions.C
echo "Combine the TGC coefficients and the matrix element squared distributions."
root -b -q ${tmp_dir}/CombineDistributions.C

# Clean up
echo "Removing temporary scripts."
rm ${tmp_dir}/*

echo "Done!"

# ------------------------------------------------------------------------------