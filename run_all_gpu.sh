#!/bin/sh

# Every BEAGLE benchmark under Benchmarks/ on the GPU, appending each run to
# results.csv. The mat variants are skipped: the MAT likelihood is computed in
# BEAST and does not use BEAGLE, so there is no GPU run to make.
#
# The list comes from benchmarks.py, so adding or removing an XML needs no edit
# here.
#
# --beagle_gpu is not a benchmark.py option; it is passed through to BEAST,
# whose argument matching is case-insensitive, so it reaches -beagle_GPU.
#
# Any arguments are passed through to every benchmark.py call, e.g.
#     sh run_all_gpu.sh --machine machine1
#     sh run_all_gpu.sh --machine machine1 --beast /path/to/beast

cd "$(dirname "$0")" || exit 1

python3 scripts/benchmarks.py | while read -r type number model; do
    case "$model" in
        *mat) continue ;;
    esac
    python3 scripts/benchmark.py "$type" "$number" "$model" --append --beagle_gpu "$@" </dev/null
done
