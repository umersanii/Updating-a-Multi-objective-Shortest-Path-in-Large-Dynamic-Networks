# #!/bin/bash

# # Usage: ./run.sh

# DATASET="datasets/rgg_n_2_20_s0/rgg_n_2_20_s0.mtx"
# NUM_OBJ=20
# NUM_CHANGES=50
# SOURCE=1

# echo "Running tests with dataset: $DATASET, num_obj: $NUM_OBJ, num_changes: $NUM_CHANGES, source: $SOURCE"

# echo "Running serial version..."
# g++ -fopenmp -O2 main_serial.cpp -o mosp_program
# ./mosp_program "$DATASET" $NUM_OBJ $NUM_CHANGES $SOURCE

# echo "Running OpenMP version..."
# g++ -fopenmp -O2 main_openmp.cpp -o mosp_program
# ./mosp_program "$DATASET" $NUM_OBJ $NUM_CHANGES $SOURCE


#!/bin/bash

# Define the values for each variable
# DATASET="datasets/rgg_n_2_20_s0/rgg_n_2_20_s0.mtx"  # Fixed dataset
DATASET="datasets/road_usa/road_usa.mtx"
NUM_OBJ=2  # Fixed num_obj
NUM_CHANGES=20  # Fixed num_changes
SOURCE=1  # Fixed source
THREADS=(1 2 4 8 16)  # Varying number of threads

echo "Compiling serial version..."
g++ -fopenmp -O2 main_serial.cpp -o mosp_program_s

echo "Compiling  OpenMP version..."
g++ -fopenmp -O2 main_openmp.cpp -o mosp_program_p

echo "Objectives, Changes, Initalizing time, dijestra time, MOSP time, serial/parallel"

echo "==================================================================================="
echo "Running for diffrent number of objectives"
# Test different num_obj (fix others)
for NUM_OBJ in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 ; do
  # echo "Testing with num_obj: $NUM_OBJ"

  # Set fixed variables
  export OMP_NUM_THREADS=4  # Example, you can change this if needed
  
  # echo "Running serial version..."
  ./mosp_program_s "$DATASET" $NUM_OBJ $NUM_CHANGES $SOURCE
  
  # echo "Running OpenMP version..."
  ./mosp_program_p "$DATASET" $NUM_OBJ $NUM_CHANGES $SOURCE
done

NUM_OBJ=15

echo "==================================================================================="
echo "Running for diffrent number of changes"
# Test different num_changes (fix others)
for NUM_CHANGES in 1 5 10 600 200; do
  # echo "Testing with num_changes: $NUM_CHANGES"
  
  # Set fixed variables
  export OMP_NUM_THREADS=4  # Example, you can change this if needed
  
  # echo "Running serial version..."
  ./mosp_program_s "$DATASET" $NUM_OBJ $NUM_CHANGES $SOURCE
  
  # echo "Running OpenMP version..."
  ./mosp_program_p "$DATASET" $NUM_OBJ $NUM_CHANGES $SOURCE
done

NUM_CHANGES=16

echo "==================================================================================="
echo "Running for diffrent number of threads"
# Test different threads (fix others)
for THREAD in 1 2 4 8 16 32; do
  echo "Testing with $THREAD threads"
  
  # Set fixed variables
  export OMP_NUM_THREADS=$THREAD  # Vary number of threads
  
  # echo "Running serial version..."
  ./mosp_program_s "$DATASET" $NUM_OBJ $NUM_CHANGES $SOURCE
  
  # echo "Running OpenMP version with $THREAD threads..."
  ./mosp_program_p "$DATASET" $NUM_OBJ $NUM_CHANGES $SOURCE
done
