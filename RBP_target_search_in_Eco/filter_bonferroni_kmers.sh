#!/usr/bin/env bash

### Edit the following codeblock based on the locations of your files, or where you want to deposit those files
for pattern_file in *_kmer_patterns.txt; do

    base=${pattern_file%_kmer_patterns.txt}

    kmers_file="${base}_kmers.txt"
    output_file="bonferroni_filtered_results/${base}_significant_kmers.txt"
    log_file="log_bonferroni/${base}_logs.txt"

###
    echo "Processing ${base}..."

    # Get Bonferroni threshold from count_patterns.py
    threshold=$(python count_patterns.py "$pattern_file" | \
        awk '/Threshold:/ {print $2}')

    if [[ ! -f "$kmers_file" ]]; then
        echo "Missing file: $kmers_file"
        continue
    fi

    # Total variants (excluding header)
    total=$(tail -n +2 "$kmers_file" | wc -l)

    # Count bad-chisq rows
    bad_chisq=$(awk -F'\t' '
        NR>1 && $0 ~ /bad-chisq/
    ' "$kmers_file" | wc -l)

    # Count retained variants
    retained=$(awk -F'\t' -v thresh="$threshold" '
        BEGIN {
            thresh_num = thresh + 0
        }
        NR>1 {
            p = $4 + 0

            if ($0 !~ /bad-chisq/ && p < thresh_num)
                count++
        }
        END {
            print count+0
        }
    ' "$kmers_file")

    # Count filtered variants
    filtered=$((total - retained))

    # Create filtered significant kmers file
    {
        head -1 "$kmers_file"

        awk -F'\t' -v thresh="$threshold" '
            BEGIN {
                thresh_num = thresh + 0
            }

            NR>1 {
                p = $4 + 0

                if ($0 !~ /bad-chisq/ && p < thresh_num)
                    print
            }
        ' "$kmers_file"

    } > "$output_file"

    # Write log
    {
        echo "Sample: ${base}"
        echo "Threshold: ${threshold}"
        echo "Total variants: ${total}"
        echo "Retained variants: ${retained}"
        echo "Filtered variants: ${filtered}"
        echo "Filtered due to bad-chisq: ${bad_chisq}"
    } > "$log_file"

done