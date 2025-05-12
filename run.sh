#!/bin/bash

# Usage: ./run.sh

# DATASET="datasets/rgg_n_2_20_s0/rgg_n_2_20_s0.mtx"
# NUM_OBJ=2
# NUM_CHANGES=5
# SOURCE=1

# g++ -fopenmp -O2 main.cpp -o mosp_program
# ./mosp_program "$DATASET" $NUM_OBJ $NUM_CHANGES $SOURCE


#######
DATASET="datasets/road_usa/road_usa.mtx"
NUM_OBJ=8
NUM_CHANGES=1
SOURCE=1

# echo "Running tests with dataset: $DATASET, num_obj: $NUM_OBJ, num_changes: $NUM_CHANGES, source: $SOURCE"

echo "Compiling and running serial version..."
g++ -fopenmp -O2 main_serial.cpp -o mosp_program
./mosp_program "$DATASET" $NUM_OBJ $NUM_CHANGES $SOURCE

echo "Compiling and running OpenMP version..."
g++ -fopenmp -O2 main_openmp.cpp -o mosp_program
./mosp_program "$DATASET" $NUM_OBJ $NUM_CHANGES $SOURCE

# ./mosp_program datasets/road_usa/road_usa.mtx 8 50 1

# DATASET="datasets/rgg_n_2_20_s0/rgg_n_2_20_s0.mtx"
# NUM_OBJ=20
# NUM_CHANGES=50
# SOURCE=1


# echo "Compiling and running OpenMP version..."
# g++ -fopenmp -O2 main_openmp.cpp -o mosp_program
# ./mosp_program "$DATASET" $NUM_OBJ $NUM_CHANGES $SOURCE


# echo "==================================================================================="
# echo "Running for diffrent number of objectives"
# # Test different num_obj (fix others)
# for NUM_OBJ in 1 2 3 4 5 6 7 8 9; do
#   echo "------------------------------------------------------"
#   echo "Testing with num_obj: $NUM_OBJ"

#   # Set fixed variables
#   export OMP_NUM_THREADS=4  # Example, you can change this if needed
  
#   echo "Running serial version..."
# #   g++ -fopenmp -O2 main_serial.cpp -o mosp_program
#   ./mosp_program "$DATASET" $NUM_OBJ $NUM_CHANGES $SOURCE
  
#   echo "Running OpenMP version..."
# #   g++ -fopenmp -O2 main_openmp.cpp -o mosp_program
#   ./mosp_program "$DATASET" $NUM_OBJ $NUM_CHANGES $SOURCE
# done
