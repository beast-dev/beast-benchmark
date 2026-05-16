#!/bin/sh

python3 benchmark.py small 1 min --append --beagle_gpu
python3 benchmark.py small 2 min --append --beagle_gpu
python3 benchmark.py big 1 min --append --beagle_gpu
python3 benchmark.py big 2 min --append --beagle_gpu
python3 benchmark.py tall 1 min --append --beagle_gpu
python3 benchmark.py wide 1 min --append --beagle_gpu
