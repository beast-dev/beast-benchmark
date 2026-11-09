#!/bin/sh

# The `min` variant of every benchmark XML under Benchmarks/, appending each run
# to results.csv. run_all_models.sh covers every model variant.
#
# The list comes from benchmarks.py, so adding or removing an XML needs no edit
# here.
#
# Any arguments are passed through to every benchmark.py call, e.g.
#     sh run_all.sh --machine machine1
#     sh run_all.sh --machine machine1 --beast /path/to/beast
#
# Pass them: without --beast, benchmark.py runs whatever `beast` is on PATH,
# which may be an installed release rather than the build under test.

cd "$(dirname "$0")" || exit 1

python3 benchmarks.py | while read -r type number model; do
    [ "$model" = "min" ] || continue
    python3 benchmark.py "$type" "$number" "$model" --append "$@" </dev/null
done
