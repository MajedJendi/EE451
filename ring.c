#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char const *argv[])
{
	int npes, myrank, status;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &npes);  // total number of processes
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	
	printf("From process %d out of %d, Hello World\n",myrank, npes);



	// int MPI_Recv (void *buf, int 3, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status);
	// MPI_Send(void *buf, int 3, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);

	int msg = 451;
	if (myrank == 0) {
		// send to 1
		printf("Process 0: Initially Msg = %d\n", msg);
		msg++;
		MPI_Send( &msg, 1, MPI_INT, 1, 200, MPI_COMM_WORLD);

	} else if (myrank == 1) {
		// recv from 1
		MPI_Recv(&msg, 1, MPI_INT, 0, 200, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		// send to 2
		printf("Process 1: Msg = %d\n", msg);
		msg++;
		MPI_Send( &msg, 1, MPI_INT, 2, 200, MPI_COMM_WORLD);

	} else if (myrank == 2) {
		// recv from 2
		MPI_Recv(&msg, 1, MPI_INT, 1, 200, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		// send to 3
		printf("Process 2: Msg = %d\n", msg);
		msg++;
		MPI_Send( &msg, 1, MPI_INT, 3, 200, MPI_COMM_WORLD);

	} else if (myrank == 3) {
		// recv from 3
		MPI_Recv(&msg, 1, MPI_INT, 2, 200, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		printf("Process 3: Msg = %d\n", msg);
		// print msg
		printf("Process 0: Received Msg = %d. Done!\n", msg);
	}
	MPI_Finalize();
	return 0;
}
