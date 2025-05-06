#!/bin/bash

# Prompt for machine information
echo "Enter the IP address of this machine:"
read THIS_IP
echo "Enter the SSH port for this machine (e.g., 2222):"
read THIS_PORT
echo "Enter the IP address of the remote machine:"
read REMOTE_IP
echo "Enter the SSH port for the remote machine (e.g., 2222):"
read REMOTE_PORT

# 1. Ensure SSH is running on both machines (if not already done)
echo "Starting SSH service if not already running..."
sudo service ssh start

# 2. Set up SSH key-based authentication (skip if keys are already set up)
echo "Setting up passwordless SSH..."
if [ ! -f "$HOME/.ssh/id_rsa" ]; then
    ssh-keygen -t rsa -N "" -f "$HOME/.ssh/id_rsa"
fi

# Copy the SSH key to the remote machine for passwordless login
ssh-copy-id -i "$HOME/.ssh/id_rsa.pub" -p $REMOTE_PORT $USER@$REMOTE_IP

# 3. Create the MPI hosts file
echo "Creating the MPI hosts file..."

HOSTS_FILE="$HOME/mpi_hosts"
echo "$THIS_IP slots=4" > $HOSTS_FILE
echo "$REMOTE_IP slots=4" >> $HOSTS_FILE

# 4. Compile a basic MPI "Hello World" program
echo "Creating the MPI Hello World program..."

MPI_PROGRAM="$HOME/hello_mpi.c"
cat > $MPI_PROGRAM << 'EOF'
#include <mpi.h>
#include <stdio.h>

int main(int argc, char** argv) {
    MPI_Init(NULL, NULL);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    printf("Hello from rank %d out of %d processes\n", world_rank, world_size);

    MPI_Finalize();
    return 0;
}
EOF

# Compile the MPI program
mpicc -o $HOME/hello_mpi $MPI_PROGRAM

# 5. Run the MPI program
echo "Running the MPI program on both machines..."
mpiexec -f $HOSTS_FILE -np 8 -H $THIS_IP:$THIS_PORT,$REMOTE_IP:$REMOTE_PORT $HOME/hello_mpi

echo "MPI job completed."
