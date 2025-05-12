#!/bin/bash

# Usage: ./run_mpi.sh

DATASET="datasets/road_usa/road_usa.mtx"
NUM_OBJ=8
NUM_CHANGES=80
SOURCE=1
NUM_PROCESSES_LIST=(2 4 8)  # You can modify the process counts as needed

echo "Compiling MPI version..."
mpic++ -O2 main_mpi.cpp -o mosp_mpi_program
export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

for PROCS in "${NUM_PROCESSES_LIST[@]}"; do
  echo "------------------------------------------------------"
  echo "Running MPI version with $PROCS processes..."
  mpirun -np $PROCS ./mosp_mpi_program "$DATASET" $NUM_OBJ $NUM_CHANGES $SOURCE
done

# #!/bin/bash
# echo "Compiling MPI version..."
# mpicxx -std=c++11 -o mosp_mpi_program main_mpi.cpp
# echo "------------------------------------------------------"
# echo "Running MPI version with 2 processes..."
# mpirun -np 2 ./mosp_mpi_program 2 2 0 > log.txt 2>&1
# echo "------------------------------------------------------"
# echo "Running MPI version with 4 processes..."
# mpirun -np 4 ./mosp_mpi_program 2 2 0 >> log.txt 2>&1
# echo "------------------------------------------------------"
# echo "Running MPI version with 8 processes..."
# mpirun -np 8 ./mosp_mpi_program 2 2 0 >> log.txt 2>&1