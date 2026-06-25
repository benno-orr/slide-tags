#!/usr/bin/env bash
set -u

csv="/n/data1/hms/scrb/chen/lab/bco/pipelines/tags_mapping/tags_mapping.csv"
sbatch_script="/n/data1/hms/scrb/chen/lab/bco/pipelines/tags_mapping/tags_mapping.sbatch"

printf '\n' >> "$csv"

echo "Starting submission script"
echo "CSV: $csv"
echo "SBATCH script: $sbatch_script"
echo

if [[ ! -f "$csv" ]]; then
    echo "ERROR: CSV file not found: $csv" >&2
    exit 1
fi

if [[ ! -f "$sbatch_script" ]]; then
    echo "ERROR: sbatch script not found: $sbatch_script" >&2
    exit 1
fi

line_num=1
submitted=0
skipped=0
failed=0

tail -n +2 "$csv" | while IFS=',' read -r run sample whitelist puck out_dir extra; do
    ((line_num++))

    echo "----------------------------------------"
    echo "Line $line_num"
    echo "run=[$run]"
    echo "sample=[$sample]"
    echo "whitelist=[$whitelist]"
    echo "puck=[$puck]"
    echo "out_dir=[$out_dir]"

    if [[ -n "${extra:-}" ]]; then
        echo "WARNING: extra CSV columns detected: [$extra]"
    fi

    if [[ -z "${run}${sample}${whitelist}${puck}${out_dir}" ]]; then
        echo "Skipping empty/malformed line"
        ((skipped++))
        continue
    fi

    if [[ "$run" != "TRUE" ]]; then
        echo "Skipping because run != TRUE"
        ((skipped++))
        continue
    fi

    echo "Submitting job..."

    sbatch_out=$(
        sbatch \
            -J "tagsmap_${sample}" \
            --export=ALL,sample="$sample",whitelist="$whitelist",puck="$puck",out_dir="$out_dir" \
            "$sbatch_script" 2>&1
    )
    status=$?

    echo "$sbatch_out"

    if [[ $status -eq 0 ]]; then
        echo "SUCCESS: submitted for sample=$sample"
        ((submitted++))
    else
        echo "ERROR: sbatch failed with exit code $status" >&2
        ((failed++))
    fi
done

echo
echo "Finished submissions"
echo "submitted=$submitted"
echo "skipped=$skipped"
echo "failed=$failed"
