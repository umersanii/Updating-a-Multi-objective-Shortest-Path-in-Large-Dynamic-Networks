#!/bin/bash

# Configuration
THIS_IP="192.168.8.180"
THIS_PORT=22
REMOTE_IP="192.168.8.240"
REMOTE_PORT=2222
REMOTE_USER="umer"

# 1. Start SSH locally
echo "Starting SSH service if not already running on local machine..."
sudo service ssh start

# 2. Start SSH on remote machine (if not already running)
echo "Ensuring SSH is running on the remote machine..."
ssh -p $REMOTE_PORT ${REMOTE_USER}@${REMOTE_IP} "sudo service ssh start"

# 3. Setup SSH key if needed
echo "Setting up passwordless SSH..."
if [ ! -f "$HOME/.ssh/id_rsa" ]; then
    ssh-keygen -t rsa -N "" -f "$HOME/.ssh/id_rsa"
fi

ssh-copy-id -i "$HOME/.ssh/id_rsa.pub" -p $REMOTE_PORT ${REMOTE_USER}@${REMOTE_IP}

# 4. Create the MPI hosts file
HOSTS_FILE="$HOME/mpi_hosts"
echo "Creating the MPI hosts file..."
echo "$THIS_IP slots=4" > $HOSTS_FILE
echo "$REMOTE_IP slots=4" >> $HOSTS_FILE
echo "MPI hosts file created at $HOSTS_FILE:"
cat $HOSTS_FILE

# 5. Compile MPI Hello World
MPI_PROGRAM="$HOME/hello_mpi.c"
cat > $MPI_PROGRAM << 'EOF'
#include <mpi.h>
#include <stdio.h>
#include <unistd.h>

int main(int argc, char** argv) {
    MPI_Init(NULL, NULL);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    char hostname[256];
    gethostname(hostname, sizeof(hostname));
    printf("Hello from rank %d out of %d processes on host %s\n", world_rank, world_size, hostname);

    MPI_Finalize();
    return 0;
}
EOF

mpicc -o $HOME/hello_mpi $MPI_PROGRAM

# 6. Test SSH manually before MPI
echo "Testing SSH to remote..."
ssh -p $REMOTE_PORT ${REMOTE_USER}@${REMOTE_IP} "hostname"

# 7. Run MPI
echo "Running the MPI program on both machines with verbose logging..."
mpiexec --hostfile $HOSTS_FILE -np 8 \
  --mca plm_rsh_agent "ssh -p $REMOTE_PORT" \
  --debug-daemons --oversubscribe \
  $HOME/hello_mpi

echo "MPI job completed."
