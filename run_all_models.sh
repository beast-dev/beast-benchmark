#!/bin/sh

# Every benchmark XML under Benchmarks/, against every model variant it has,
# appending each run to results.csv. run_all.sh covers the `min` set only.
#
# The list comes from benchmarks.py, so adding or removing an XML needs no edit
# here.
#
# Any arguments are passed through to every benchmark.py call, e.g.
#     sh run_all_models.sh --machine machine1
#     sh run_all_models.sh --machine machine1 --beast /path/to/beast
#
# A joint-density mismatch makes benchmark.py exit 2; that does not stop the
# sweep, so a stale header value does not cost the rest of the runs.

cd "$(dirname "$0")" || exit 1

python3 benchmarks.py | while read -r type number model; do
    python3 benchmark.py "$type" "$number" "$model" --append "$@" </dev/null
done
