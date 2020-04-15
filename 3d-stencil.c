#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#define BACK 0
#define FRONT 1
#define EAST 2
#define WEST 3
#define NORTH 4
#define SOUTH 5

void create_3dstencil_dt(int nx, int ny, int nz, int nvar, int rank, MPI_Datatype *backfront_ptr, MPI_Datatype *eastwest_ptr, MPI_Datatype *northsouth_ptr){
	MPI_Datatype backfront;	
	MPI_Datatype eastwest;	
	MPI_Datatype northsouth;
	char errstr[256];
	int errlen;
	int mpi_errno = MPI_SUCCESS;
	/* backfront */
	mpi_errno = MPI_Type_vector(nvar, ny * nx, nz * ny * nx, MPI_DOUBLE, &backfront);
	if(mpi_errno != MPI_SUCCESS)
		goto fn_fail;

	mpi_errno = MPI_Type_vector(nvar * ny * nz, 1, nx, MPI_DOUBLE, &eastwest);
	if(mpi_errno != MPI_SUCCESS)
		goto fn_fail;

	mpi_errno = MPI_Type_vector(nvar * nz, nx, nx * ny, MPI_DOUBLE, &northsouth);
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

	return;
}

int main(int argc, char *argv[]){
	MPI_Datatype backfront;	
	MPI_Datatype eastwest;	
	MPI_Datatype northsouth;	

	int rank, nprocs;
	int nx, ny, nz;
	int px, py, pz, lx, ly, lz; /* coordinate */
	int lnx, lny, lnz;
	int nvar, iteration;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	parse_params(argc, argv, rank, nprocs, &px, &py, &pz, &nx, &ny, &nz, &nvar, &iteration, &lx, &ly, &lz, &lnx, &lny, &lnz);
	create_3dstencil_dt(lnx + 2, lny + 2, lnz + 2, nvar, rank, &backfront, &eastwest, &northsouth);
	/* begin 3d 7 point stencil */
	


	free_3dstencil_dt(backfront, eastwest, northsouth);
	MPI_Finalize();
	return 0;
}	