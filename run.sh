#!/bin/bash

# Usage: ./run.sh

DATASET="datasets/rgg_n_2_20_s0/rgg_n_2_20_s0.mtx"
NUM_OBJ=2
NUM_CHANGES=5
SOURCE=1

g++ -fopenmp -O2 main.cpp -o mosp_program
./mosp_program "$DATASET" $NUM_OBJ $NUM_CHANGES $SOURCE
