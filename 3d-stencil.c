#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#define BACK 0
#define FRONT 1
#define EAST 2
#define WEST 3
#define NORTH 4
#define SOUTH 5
#define MAX_NEIGHBOR 6
#define MAX_REQUESTS 12

void create_3dstencil_dt(int bnx, int bny, int bnz, int nvar, int rank, MPI_Datatype *backfront_ptr, MPI_Datatype *eastwest_ptr, MPI_Datatype *northsouth_ptr){
	MPI_Datatype backfront;	
	MPI_Datatype eastwest;	
	MPI_Datatype northsouth;
	char errstr[256];
	int errlen;
	int mpi_errno = MPI_SUCCESS;
	/* backfront */
	mpi_errno = MPI_Type_vector(nvar, bny * bnx, bnz * bny * bnx, MPI_DOUBLE, &backfront);
	if(mpi_errno != MPI_SUCCESS)
		goto fn_fail;

	mpi_errno = MPI_Type_vector(nvar * bny * bnz, 1, bnx, MPI_DOUBLE, &eastwest);
	if(mpi_errno != MPI_SUCCESS)
		goto fn_fail;

	mpi_errno = MPI_Type_vector(nvar * bnz, bnx, bnx * bny, MPI_DOUBLE, &northsouth);
	if(mpi_errno != MPI_SUCCESS)
		goto fn_fail;

	mpi_errno = MPI_Type_commit(&backfront);
	if(mpi_errno != MPI_SUCCESS)
		goto fn_fail;

	mpi_errno = MPI_Type_commit(&eastwest);
	if(mpi_errno != MPI_SUCCESS)
		goto fn_fail;

	mpi_errno = MPI_Type_commit(&northsouth);
	if(mpi_errno != MPI_SUCCESS)
		goto fn_fail;

	*backfront_ptr = backfront;
	*eastwest_ptr = eastwest;
	*northsouth_ptr = northsouth;
	// int bfz, ewz, nsz;
	// MPI_Type_size(backfront, &bfz);
	// MPI_Type_size(eastwest, &ewz);
	// MPI_Type_size(northsouth, &nsz);
	// if(rank == 0)
	// 	printf("%d %d %d\n", bfz, ewz, nsz);

fn_exit:
	return;
fn_fail:
	MPI_Error_string(mpi_errno, errstr, &errlen);
	printf("rank %d - create datatype failed, error msg: %s\n", rank, errstr);
	exit(1);
}

void free_3dstencil_dt(MPI_Datatype backfront, MPI_Datatype eastwest, MPI_Datatype northsouth){
	MPI_Type_free(&backfront);
	MPI_Type_free(&eastwest);
	MPI_Type_free(&northsouth);
}

#define REQUIRE_PARAMS 9
#define MAX_LOCAL_ELEMENTS 201326592
void parse_params(int argc, char *argv[], int rank, int nprocs, int *px, int *py, int *pz, int *nx, int *ny, int *nz, int *nvar, int *iteration,
		int *lx, int *ly, int *lz, int *lnx, int *lny, int *lnz){
	/*
		Usage:
		mpirun -n np ./3d-stencil px py pz x y z nvar iteration
	*/
	if(argc != REQUIRE_PARAMS){
		printf("Please Input: px py pz x y z nvar iteration\n");
		exit(1);
	}

	*px = atoi(argv[1]);
	*py = atoi(argv[2]);
	*pz = atoi(argv[3]);
	if(*px * *py * *pz != nprocs){
		printf("rank %d - px * py * pz (%d) != nprocs (%d)\n", rank, *px * *py * *pz, nprocs);
		exit(1);
	}

	*nx = atoi(argv[4]);
	*ny = atoi(argv[5]);
	*nz = atoi(argv[6]);

	*nvar = atoi(argv[7]);
	*iteration = atoi(argv[8]);

	*lx = rank % (*px);
	*ly = rank / (*px) % (*py);
	*lz = rank / (*px) / (*py);

	*lnx = ((*lx) + 1) * (*nx) / (*px) - (*lx) * (*nx) / (*px);
	*lny = ((*ly) + 1) * (*ny) / (*py) - (*ly) * (*ny) / (*py);
	*lnz = ((*lz) + 1) * (*nz) / (*pz) - (*lz) * (*nz) / (*pz);
	int total_elem = (*lnx + 2) * (*lny + 2) * (*lnz + 2) * (*nvar);
	if(total_elem >= MAX_LOCAL_ELEMENTS){
		printf("rank %d - total_elem %d is too large (>= %d)\n", rank, total_elem, MAX_LOCAL_ELEMENTS);
		exit(1);
	}

	return;
}

void init_neighbor(int lx, int ly, int lz, int px, int py, int pz, int neighbors[]){
	int i;
	for(i=0;i<MAX_NEIGHBOR;++i)
		neighbors[i] = -1;
	/* back neighbor */
	if(lz - 1 >= 0)
		neighbors[BACK] = lx + px * ly + px * py * (lz - 1);

	/* front neighbor */
	if(lz + 1 < pz)
		neighbors[FRONT] = lx + px * ly + px * py * (lz + 1);

	/* east neighbor */
	if(lx - 1 >= 0)
		neighbors[EAST] = lx - 1 + px * ly + px * py * lz;

	/* west neighbor */
	if(lx + 1 < px)
		neighbors[WEST] = lx + 1 + px * ly + px * py * lz;

	/* north neighbor */
	if(ly - 1 >= 0)
		neighbors[NORTH] = lx + px * (ly - 1) + px * py * lz;

	/* south neighbor */
	if(ly + 1 < py)
		neighbors[SOUTH] = lx + px * (ly + 1) + px * py * lz;

	// int rank;
	// char buf[64] = "\0";
	// char tmp[64];
	// MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// for(i=0;i<MAX_NEIGHBOR;++i){
	// 	sprintf(tmp, "%d ", neighbors[i]);	
	// 	strcat(buf, tmp);
	// }
	// printf("rank %d - %s\n", rank, buf);
	// fflush(stdout);
		
}

void init_buf(int bnx, int bny, int bnz, int nvar, int working_buf, double *dat_buf[]){
	int i,j,k;

	dat_buf[working_buf] = (double *) malloc(bnx * bny * bnz * nvar * sizeof(double));
	dat_buf[working_buf ^ 1] = (double *) malloc(bnx * bny * bnz * nvar * sizeof(double));

	for (k = 0; k < bnz; ++k)
		for (j = 0; j < bny; ++j)
			for (i = 0; i < bnx; ++i) {
				dat_buf[working_buf][i + bnx * j + bnx * bny * k] = (double)(rand() % 10) + (double) (rand() % 10) / 100.0;
				dat_buf[working_buf ^ 1][i + bnx * j + bnx * bny * k] = 0.0;
			}

	return;
}

void halo_irecv_isend(int bnx, int bny, int bnz, int neighbors[], MPI_Datatype backfront, MPI_Datatype eastwest, MPI_Datatype northsouth, double *work_buf, int *count_ptr, MPI_Request requests[]){
	int i, count;
	count = 0;
	for(i=0;i<MAX_REQUESTS;++i){
		requests[i] = MPI_REQUEST_NULL;
	}
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	/* back */
	if (neighbors[BACK] != -1) {
		count += 2;
		// printf("rank %d - begin BACK\n", rank);
		// fflush(stdout);
		MPI_Irecv(work_buf, 1, backfront, neighbors[BACK], 0, MPI_COMM_WORLD, &requests[MAX_NEIGHBOR + BACK]);
		MPI_Isend(work_buf + bnx * bny, 1, backfront, neighbors[BACK], 0, MPI_COMM_WORLD, &requests[BACK]);

		// printf("rank %d - complete BACK\n", rank);
		// fflush(stdout);
	}
	/* front */
	if (neighbors[FRONT] != -1) {
		count += 2;

		// printf("rank %d - begin FRONT\n", rank);
		// fflush(stdout);
		MPI_Irecv(work_buf + bnx * bny * (bnz - 1), 1, backfront, neighbors[FRONT], 0, MPI_COMM_WORLD, &requests[MAX_NEIGHBOR + FRONT]);
		MPI_Isend(work_buf + bnx * bny * (bnz - 2), 1, backfront, neighbors[FRONT], 0, MPI_COMM_WORLD, &requests[FRONT]);

		// printf("rank %d - complete FRONT\n", rank);
		// fflush(stdout);
	}
	/* east */
	if (neighbors[EAST] != -1) {
		count += 2;
		// printf("rank %d - begin EAST %d\n", rank, neighbors[EAST]);
		// fflush(stdout);
		MPI_Irecv(work_buf, 1, eastwest, neighbors[EAST], 0, MPI_COMM_WORLD, &requests[MAX_NEIGHBOR + EAST]);
		MPI_Isend(work_buf + 1, 1, eastwest, neighbors[EAST], 0, MPI_COMM_WORLD, &requests[EAST]);
		// printf("rank %d - complete EAST\n", rank);
		// fflush(stdout);
	}

	/* west */
	if (neighbors[WEST] != -1) {
		count += 2;
		// printf("rank %d - begin WEST %d\n", rank, neighbors[WEST]);
		// fflush(stdout);
		MPI_Irecv(work_buf + bnx - 1, 1, eastwest, neighbors[WEST], 0, MPI_COMM_WORLD, &requests[MAX_NEIGHBOR + WEST]);
		MPI_Isend(work_buf + bnx - 2, 1, eastwest, neighbors[WEST], 0, MPI_COMM_WORLD, &requests[WEST]);
		// printf("rank %d - complete WEST\n", rank);
		// fflush(stdout);
	}
	
	/* north */
	if (neighbors[NORTH] != -1) {
		count += 2;
		MPI_Irecv(work_buf, 1, northsouth, neighbors[NORTH], 0, MPI_COMM_WORLD, &requests[MAX_NEIGHBOR + NORTH]);
		MPI_Isend(work_buf + bnx, 1, northsouth, neighbors[NORTH], 0, MPI_COMM_WORLD, &requests[NORTH]);
	}
	// printf("rank %d - complete NORTH\n", rank);
	// fflush(stdout);
	/* south */
	if(neighbors[SOUTH] != -1){
		count += 2;
		MPI_Irecv(work_buf + bnx * (bny - 1), 1, northsouth, neighbors[SOUTH], 0, MPI_COMM_WORLD, &requests[MAX_NEIGHBOR + SOUTH]);
		MPI_Isend(work_buf + bnx * (bny - 2), 1, northsouth, neighbors[SOUTH], 0, MPI_COMM_WORLD, &requests[SOUTH]);
	}
	// printf("rank %d - complete SOUTH\n", rank);
	// fflush(stdout);
	*count_ptr = count;
	return;
}

void halo_exchange(int bnx, int bny, int bnz, int neighbors[], MPI_Datatype backfront, MPI_Datatype eastwest, MPI_Datatype northsouth, double *work_buf){
	int count = 0;
	MPI_Request requests[MAX_REQUESTS];
	// int rank;
	// MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	// printf("rank %d - begin irecv isend\n", rank);
	// fflush(stdout);
	halo_irecv_isend(bnx, bny, bnz, neighbors, backfront, eastwest, northsouth, work_buf, &count, requests);
	// printf("rank %d - end irecv isend, count %d\n", rank, count);
	// fflush(stdout);
	MPI_Waitall(count, requests, MPI_STATUSES_IGNORE);

	return;
}


#define SEVEN_POINTS 7.0
void stencil(double *work_buf, double *dest_buf, int bnx, int bny, int bnz){
	int i, j, k;
	for (k = 1; k < bnz - 1; ++k)
		for (j = 1; j < bny - 1; ++j)
			for (i = 1; i < bnx - 1; ++i) {
				int origin = i + bnx * j + bnx * bny * k;
				int back_pos = i + bnx * j + bnx * bny * (k - 1);
				int front_pos = i + bnx * j + bnx * bny * (k + 1);
				int east_pos = i - 1 + bnx * j  + bnx * bny * k;
				int west_pos = i + 1 + bnx * j + bnx * bny * k;
				int north_pos = i + bnx * (j - 1) + bnx * bny * k;
				int south_pos = i + bnx * (j + 1) + bnx * bny * k;
				dest_buf[origin] = (work_buf[origin] + work_buf[back_pos] + work_buf[front_pos] + work_buf[east_pos] + work_buf[west_pos] + work_buf[north_pos] + work_buf[south_pos]) / SEVEN_POINTS;
			}
	return;
}	

int main(int argc, char *argv[]){
	MPI_Datatype backfront;	
	MPI_Datatype eastwest;	
	MPI_Datatype northsouth;	

	int i;
	int rank, nprocs;
	int nx, ny, nz;
	int px, py, pz, lx, ly, lz; /* coordinate */
	int lnx, lny, lnz;
	int nvar, iteration;
	int working_buf = 0;
	int neighbors[MAX_NEIGHBOR];
	double *dat_buf[2];
	int bnx, bny, bnz;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	parse_params(argc, argv, rank, nprocs, &px, &py, &pz, &nx, &ny, &nz, &nvar, &iteration, &lx, &ly, &lz, &lnx, &lny, &lnz);
	bnx = lnx + 2;
	bny = lny + 2;
	bnz = lnz + 2;
	// printf("rank %d - lx %d ly %d lz %d, id %d\n", rank, lx, ly, lz, getpid());
	// fflush(stdout);
	// sleep(8);
	create_3dstencil_dt(bnx, bny, bnz, nvar, rank, &backfront, &eastwest, &northsouth);
	// printf("rank %d - begin neighbors\n", rank);
	// fflush(stdout);
	init_neighbor(lx, ly, lz, px, py, pz, neighbors);
	// printf("rank %d - end neighbors, bn %d %d %d nvar %d\n", rank, bnx, bny, bnz, nvar);
	// fflush(stdout);
	init_buf(bnx, bny, bnz, nvar, working_buf, dat_buf);
	// printf("rank %d - end init buf\n", rank);
	// fflush(stdout);
	/* begin 3d 7 point stencil */
	MPI_Barrier(MPI_COMM_WORLD);
	double time = MPI_Wtime();
	for(i=0;i<iteration;++i){
		/* exchange boundary */
		// if(rank == 0){
		// 	printf("ITERATION - %d\n", i);
		// 	fflush(stdout);
		// }
		// printf("rank %d - begin exchange\n", rank);
		// fflush(stdout);
		halo_exchange(bnx, bny, bnz, neighbors, backfront, eastwest, northsouth, dat_buf[working_buf]);
		// printf("rank %d - end exchange\n", rank);
		// fflush(stdout);
		/* compute */
		stencil(dat_buf[working_buf], dat_buf[working_buf ^ 1], bnx, bny, bnz);

#ifdef AMPI_LOAD_BALANCE
		const int steps = 5;
		if(i % steps == 0){
			AMPI_Migrate(AMPI_INFO_LB_SYNC);
		}
#endif

		working_buf ^= 1;
	}
	MPI_Barrier(MPI_COMM_WORLD);

	time = MPI_Wtime() - time;
	if(rank == 0){
		printf("%d %.3lf\n", nprocs, time);
		fflush(stdout);
	}

	free_3dstencil_dt(backfront, eastwest, northsouth);
	free(dat_buf[0]);
	free(dat_buf[1]);
	MPI_Finalize();
	return 0;
}	