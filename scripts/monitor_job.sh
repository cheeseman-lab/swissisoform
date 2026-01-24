#!/bin/bash
JOB_ID=$1
while squeue -j $JOB_ID &>/dev/null; do
    echo "$(date): Job $JOB_ID status:"
    squeue -j $JOB_ID -o "%.18i %.9P %.2t %.10M" | head -10
    echo ""
    sleep 300  # Check every 5 minutes
done
echo "$(date): Job $JOB_ID completed"
