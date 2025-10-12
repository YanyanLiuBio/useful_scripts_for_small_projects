#!/bin/bash
QUEUE_NAME="spot_low_priority"
REGION="us-east-2"
REASON="Terminating all jobs in queue"

# Job statuses to terminate
STATUSES=("SUBMITTED" "PENDING" "RUNNABLE" "STARTING" "RUNNING")

for STATUS in "${STATUSES[@]}"; do
  echo "Checking jobs with status: $STATUS"

  NEXT_TOKEN=""
  while :; do
    if [ -z "$NEXT_TOKEN" ]; then
      RESPONSE=$(aws batch list-jobs \
        --job-queue "$QUEUE_NAME" \
        --job-status "$STATUS" \
        --region "$REGION")
    else
      RESPONSE=$(aws batch list-jobs \
        --job-queue "$QUEUE_NAME" \
        --job-status "$STATUS" \
        --next-token "$NEXT_TOKEN" \
        --region "$REGION")
    fi

    # Extract job IDs
    echo "$RESPONSE" | jq -r '.jobSummaryList[].jobId' > job_ids.txt

    # Terminate each job
    while IFS= read -r job_id; do
      if [ -n "$job_id" ]; then
        aws batch terminate-job \
          --job-id "$job_id" \
          --reason "$REASON" \
          --region "$REGION"
        echo "Terminated job: $job_id ($STATUS)"
      fi
    done < job_ids.txt

    # Handle pagination
    NEXT_TOKEN=$(echo "$RESPONSE" | jq -r '.nextToken // empty')
    [ -z "$NEXT_TOKEN" ] && break
  done
done

echo "✅ All jobs in queue '$QUEUE_NAME' have been terminated (statuses: ${STATUSES[*]})."
