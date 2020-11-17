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
    echo "usage: ./run_process.sh --process=[ww_sl0muq/...] --output-base=[...]"
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
# Script environment
dir="$( cd "$( dirname "${BASH_SOURCE[0]}"  )" && pwd  )" # path of this macro

# Standard HTCondor environment
condor_directory=${dir}/SubmitScripts
submit_script=${condor_directory}/submit_script.submit

# Output directory
output_dir=${output_base}/${process}/omega
if [[ ! -d ${output_dir} ]] ; then # Create if not existing
mkdir -p ${output_dir}
fi

# Condor logging output 
condor_output_dir=${output_base}/CondorOutput
if [[ ! -d ${condor_output_dir} ]] ; then # Create if not existing
  mkdir -p ${condor_output_dir}
fi

# ------------------------------------------------------------------------------
# Check if the needed binary is available
bin_dir=${dir}/../bin
bin_executable=${bin_dir}/grid_${process}.bin
if ! test -f "$bin_executable"; then
  >&2 echo "Binary ${bin_executable} not found!"; exit 1
fi

# ------------------------------------------------------------------------------
# Submit the calculation to the BIRD system

echo "Start initalizing process ${process} running."

# Array to collect condor job IDs so that I can keep track of if any are still running
condor_job_IDs=()
    
# Start jobs from condor directory
cd ${condor_directory}

# Function that creates a HTCondor job
add_condor_job() {
  # Command to be run on BIRD cluster (bin set later)
  command_string="cd ${output_dir} \&\& ${bin_executable} ${BIN}"
  
  # Submit job to HTCondor using standard submitting setup
  # -> Start Marlin job and keep track of job ID to know when it's done
  condor_job_output=$(condor_submit ${submit_script} log_dir=${condor_output_dir} arguments="${command_string}")
  
  # Split output up by spaces and only read last part (which is cluster ID).
  # Details at: https://stackoverflow.com/questions/3162385/how-to-split-a-string-in-shell-and-get-the-last-field
  condor_job_IDs+="${condor_job_output##* } "
  
  # Wait a little to not overload the condor scheduler
  sleep 0.05s
}

# Binning process dependent (depends on fortran grid file)
if [[ ${process} == "ww_sl0muq" ]]; then
  for i in {1..20};
  do
    for j in {1..10};
    do	
      for k in {1..10};
      do
        BIN="${i} ${j} ${k}"
        add_condor_job
      done
    done
  done
else 
  >&2 echo "Process ${process} unknown!"; exit 1
fi
      
# ------------------------------------------------------------------------------
# Wait for jobs to finish (technically not needed here)
      
echo "Waiting for jobs to finish."
for job_ID in ${condor_job_IDs[@]}; do
  job_log_path=$(ls ${condor_output_dir}/${job_ID}*.log)
  
  # Write into variable to suppress spammy output.
  # Timeout after 15 seconds to restart command, else it gets stuck sometimes.
  status=124 # Start with error code to get into loop
  while [ $status -eq 124 ]; do
    wait_output=$(timeout 15 condor_wait ${job_log_path}) 
    status=$? # Update status -> 124 means timeout, else success
  done
done
echo "Jobs finished!"

# ------------------------------------------------------------------------------