#!/bin/bash

# Path to the Julia script
JULIA_SCRIPT="./run_MCMC_server.jl"

# Check if the script exists
if [ ! -f "$JULIA_SCRIPT" ]; then
  echo "Error: $JULIA_SCRIPT not found in current directory."
  exit 1
fi

# Launch the script 15 times using nohup with a 1-second delay
for i in {1..25}
do
  echo "Starting instance $i with nohup..."
  nohup julia --startup-file=no --project=@temp "$JULIA_SCRIPT" > "nohup_$i.out" 2>&1 &
  sleep 2
done

echo "All 30 instances started with nohup and a 2-second delay between each."

