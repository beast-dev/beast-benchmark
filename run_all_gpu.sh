#!/bin/sh

# Every dataset against every BEAGLE model variant on the GPU, appending each
# run to results.csv. The mat variants are left out: the MAT likelihood is
# computed in BEAST and does not use BEAGLE, so there is no GPU run to make.
#
# --beagle_gpu is not a benchmark.py option; it is passed through to BEAST,
# whose argument matching is case-insensitive, so it reaches -beagle_GPU.
#
# Any arguments are passed through to every benchmark.py call, e.g.
#     sh run_all_gpu.sh --machine machine1
#     sh run_all_gpu.sh --machine machine1 --beast /path/to/beast
#
# A joint-density mismatch makes benchmark.py exit 2; that does not stop the
# sweep, so a stale header value does not cost the rest of the runs.

run() {
    t=$1; n=$2; m=$3; shift 3
    python3 benchmark.py "$t" "$n" "$m" --append --beagle_gpu "$@"
}

run small 1 min "$@"
run small 1 gamma "$@"
run small 1 skygrid "$@"

run small 2 min "$@"
run small 2 gamma "$@"
run small 2 skygrid "$@"

run big 1 min "$@"
run big 1 gamma "$@"
run big 1 skygrid "$@"

run big 2 min "$@"
run big 2 gamma "$@"
run big 2 skygrid "$@"

run big 3 min "$@"
run big 3 gamma "$@"
run big 3 skygrid "$@"

run tall 1 min "$@"
run tall 1 gamma "$@"
run tall 1 skygrid "$@"

run wide 1 min "$@"
run wide 1 gamma "$@"
run wide 1 skygrid "$@"
