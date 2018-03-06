#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>
#include <ctype.h>
#include <wordexp.h>

#include "allvars.h"
#include "proto.h"

#ifdef ANALYTIC_GRAVITY
double*** Allocate_matrix3D(const int dim1,const int dim2, const int dim3)
{

    int i,j;
    double *** mat3D = (double ***)malloc(dim1*sizeof(double**));

        for (i = 0; i< dim1; i++) {

         mat3D[i] = (double **) malloc(dim2*sizeof(double *));

          for (j = 0; j < dim2; j++) {

              mat3D[i][j] = (double *)malloc(dim3*sizeof(double));
          }

        }
    return mat3D;
}



void allocate_acc_finegrid()
{
    hid_t   file, dataset_pot, dataspace_pot, datatype;
	hid_t   dataset_x, dataset_y, dataset_z, dataspace_x, dataspace_y, dataspace_z;
    size_t      size;                    /* size of the data element stored in file */
    hsize_t     dims_out[3];             /* dataset dimensions */
    herr_t      status=0;
	int rank, status_n;
    int Nx, Ny, Nz;

    // Open the file and the dataset //
 
    //file = H5Fopen("MWpotential.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
	file = H5Fopen("Potential_FineGrid.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
	
    dataset_pot = H5Dopen(file, "/Potentials/PotentialSormani");
    dataspace_pot = H5Dget_space(dataset_pot);

    //Get datatype, size, rank and dimensions //
    datatype  = H5Dget_type(dataset_pot);     

    size  = H5Tget_size(datatype);
    //printf(" Data size is %d \n", (int)size);

    rank      = H5Sget_simple_extent_ndims(dataspace_pot);
    status_n  = H5Sget_simple_extent_dims(dataspace_pot, dims_out, NULL);
    if (ThisTask == 0) printf("Reading file MWpotential (fine): rank %d, dimensions %lu x %lu x %lu\n", rank,
      (unsigned long)(dims_out[0]), (unsigned long)(dims_out[1]), (unsigned long)(dims_out[2]));  
	  
	Nx = dims_out[0];
	Ny = dims_out[1];
	Nz = dims_out[2];

    All.potential = malloc(Nx*Ny*Nz*sizeof(double));

    status += H5Dread(dataset_pot, datatype, dataspace_pot, dataspace_pot, H5P_DEFAULT, All.potential);

	//printf("arriva qui! %e \n", All.potential_tot[100][100][100]);
	
    dataset_x = H5Dopen(file, "/Potentials/Coord_X");
    dataspace_x = H5Dget_space(dataset_x);
    dataset_y = H5Dopen(file, "/Potentials/Coord_Y");
    dataspace_y = H5Dget_space(dataset_y);
    dataset_z = H5Dopen(file, "/Potentials/Coord_Z");
    dataspace_z = H5Dget_space(dataset_z);
	
    double *XX = malloc(Nx*sizeof(double));
	double *YY = malloc(Ny*sizeof(double));
	double *ZZ = malloc(Nz*sizeof(double));
	
    status += H5Dread(dataset_x, datatype, dataspace_x, dataspace_x, H5P_DEFAULT, XX);
    status += H5Dread(dataset_y, datatype, dataspace_y, dataspace_y, H5P_DEFAULT, YY);
	status += H5Dread(dataset_z, datatype, dataspace_z, dataspace_z, H5P_DEFAULT, ZZ);
	
	if(status != 0)
	{
		printf("ERROR IN HDF5 READING!!! \n");
	    exit(0);
	}
	// End of reading file //
	
	All.Nx = Nx;
	All.Ny = Ny;
	All.Nz = Nz;		
	All.deltax = XX[1]-XX[0];
	All.deltay = YY[1]-YY[0];
	All.deltaz = ZZ[1]-ZZ[0];
	All.xx0 = XX[0];
	All.yy0 = YY[0];
	All.zz0 = ZZ[0]; 
	
	/*Conversion in unit code*/
	
	All.deltax *= CM_PER_KPC/All.UnitLength_in_cm;
	All.deltay *= CM_PER_KPC/All.UnitLength_in_cm;
	All.deltaz *= CM_PER_KPC/All.UnitLength_in_cm;
	All.xx0 *= CM_PER_KPC/All.UnitLength_in_cm;
	All.yy0 *= CM_PER_KPC/All.UnitLength_in_cm;
	All.zz0 *= CM_PER_KPC/All.UnitLength_in_cm; 
	
	//printf("Step and first values: %e %e %e %e %e %e \n", All.deltax,All.deltay,All.deltaz,XX[Nx-1],YY[Nx-1],ZZ[Nx-1]);
	
	for(int i=0; i < Nx*Ny*Nz; i++) All.potential[i]*=1e10/(All.UnitVelocity_in_cm_per_s * All.UnitVelocity_in_cm_per_s);
		
	free(XX);
	free(YY);
	free(ZZ);
	//free(potential_tot);
}

void allocate_acc_coarsegrid()
{
    hid_t   file, dataset_pot, dataspace_pot, datatype;
	hid_t   dataset_x, dataset_y, dataset_z, dataspace_x, dataspace_y, dataspace_z;
    size_t      size;                    /* size of the data element stored in file */
    hsize_t     dims_out[3];             /* dataset dimensions */
    herr_t      status=0;
	int rank, status_n;
    int Nx, Ny, Nz;

    // Open the file and the dataset //
 
    //file = H5Fopen("MWpotential.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
	file = H5Fopen("Potential_CoarseGrid.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
	
    dataset_pot = H5Dopen(file, "/Potentials/PotentialSormani");
    dataspace_pot = H5Dget_space(dataset_pot);

    //Get datatype, size, rank and dimensions //
    datatype  = H5Dget_type(dataset_pot);     

    size  = H5Tget_size(datatype);
    //printf(" Data size is %d \n", (int)size);

    rank      = H5Sget_simple_extent_ndims(dataspace_pot);
    status_n  = H5Sget_simple_extent_dims(dataspace_pot, dims_out, NULL);
    if (ThisTask == 0) printf("Reading file MWpotential (coarse): rank %d, dimensions %lu x %lu x %lu\n", rank,
      (unsigned long)(dims_out[0]), (unsigned long)(dims_out[1]), (unsigned long)(dims_out[2]));  
	  
	Nx = dims_out[0];
	Ny = dims_out[1];
	Nz = dims_out[2];

    All.coarse_potential = malloc(Nx*Ny*Nz*sizeof(double));

    status += H5Dread(dataset_pot, datatype, dataspace_pot, dataspace_pot, H5P_DEFAULT, All.coarse_potential);

	//printf("arriva qui! %e \n", All.potential_tot[100][100][100]);
	
    dataset_x = H5Dopen(file, "/Potentials/Coord_X");
    dataspace_x = H5Dget_space(dataset_x);
    dataset_y = H5Dopen(file, "/Potentials/Coord_Y");
    dataspace_y = H5Dget_space(dataset_y);
    dataset_z = H5Dopen(file, "/Potentials/Coord_Z");
    dataspace_z = H5Dget_space(dataset_z);
	
    double *XX = malloc(Nx*sizeof(double));
	double *YY = malloc(Ny*sizeof(double));
	double *ZZ = malloc(Nz*sizeof(double));
	
    status += H5Dread(dataset_x, datatype, dataspace_x, dataspace_x, H5P_DEFAULT, XX);
    status += H5Dread(dataset_y, datatype, dataspace_y, dataspace_y, H5P_DEFAULT, YY);
	status += H5Dread(dataset_z, datatype, dataspace_z, dataspace_z, H5P_DEFAULT, ZZ);
	
	if(status != 0)
	{
		printf("ERROR IN HDF5 READING!!! \n");
	    exit(0);
	}
	// End of reading file //
	
	All.coarse_Nx = Nx;
	All.coarse_Ny = Ny;
	All.coarse_Nz = Nz;		
	All.coarse_deltax = XX[1]-XX[0];
	All.coarse_deltay = YY[1]-YY[0];
	All.coarse_deltaz = ZZ[1]-ZZ[0];
	
	/*Conversion in unit code*/
	
	All.coarse_deltax *= CM_PER_KPC/All.UnitLength_in_cm;
	All.coarse_deltay *= CM_PER_KPC/All.UnitLength_in_cm;
	All.coarse_deltaz *= CM_PER_KPC/All.UnitLength_in_cm;
	
	//printf("Step and first values: %e %e %e %e %e %e \n", All.coarse_deltax,All.coarse_deltay,All.coarse_deltaz,XX[0],YY[0],ZZ[0]);
	
	for(int i=0; i < Nx*Ny*Nz; i++) All.coarse_potential[i]*=1e10/(All.UnitVelocity_in_cm_per_s * All.UnitVelocity_in_cm_per_s);
		
	free(XX);
	free(YY);
	free(ZZ);
	//free(potential_tot);
}
#endif