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



void allocate_acc()
{
    hid_t   file, dataset_bdh, dataset_bar, dataspace_bdh, dataspace_bar, datatype;
	hid_t   dataset_x, dataset_y, dataset_z, dataspace_x, dataspace_y, dataspace_z;
    size_t      size;                    /* size of the data element stored in file */
    hsize_t     dims_out[3];             /* dataset dimensions */
    herr_t      status=0;
	int rank, status_n;
    int Nx, Ny, Nz;

    // Open the file and the dataset //
 
    file = H5Fopen("MWpotential.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
	
    dataset_bdh = H5Dopen(file, "/Potentials/BULGE_DISK_HALO");
    dataspace_bdh = H5Dget_space(dataset_bdh);
    dataset_bar = H5Dopen(file, "/Potentials/BAR");
    dataspace_bar = H5Dget_space(dataset_bar);

    //Get datatype, size, rank and dimensions //
    datatype  = H5Dget_type(dataset_bdh);     

    size  = H5Tget_size(datatype);
    printf(" Data size is %d \n", (int)size);

    rank      = H5Sget_simple_extent_ndims(dataspace_bdh);
    status_n  = H5Sget_simple_extent_dims(dataspace_bdh, dims_out, NULL);
    printf("Reading file MWpotential: rank %d, dimensions %lu x %lu x %lu\n", rank,
      (unsigned long)(dims_out[0]), (unsigned long)(dims_out[1]), (unsigned long)(dims_out[2]));  
	  
	Nx = dims_out[0];
	Ny = dims_out[1];
	Nz = dims_out[2];
	
    double (*potential_tot)[Ny][Nz] = malloc(Nx*sizeof(*potential_tot));
	double (*potential_bar)[Ny][Nz] = malloc(Nx*sizeof(*potential_bar));

    status += H5Dread(dataset_bdh, datatype, dataspace_bdh, dataspace_bdh, H5P_DEFAULT, potential_tot);
    status += H5Dread(dataset_bar, datatype, dataspace_bar, dataspace_bar, H5P_DEFAULT, potential_bar); 
	
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
	
	All.deltax = XX[1]-XX[0];
	All.deltay = YY[1]-YY[0];
	All.deltaz = ZZ[1]-ZZ[0];
	All.xx0 = XX[0];
	All.yy0 = YY[0];
	All.zz0 = ZZ[0]; 
	
	All.accx = Allocate_matrix3D(Nx, Ny, Nz);
	All.accy = Allocate_matrix3D(Nx, Ny, Nz);
	All.accz = Allocate_matrix3D(Nx, Ny, Nz);
	All.accx_bar = Allocate_matrix3D(Nx, Ny, Nz);
	All.accy_bar = Allocate_matrix3D(Nx, Ny, Nz);
	All.accz_bar = Allocate_matrix3D(Nx, Ny, Nz);	
	
	/*Conversion in unit code*/
	
	All.deltax *= CM_PER_KPC/All.UnitLength_in_cm;
	All.deltay *= CM_PER_KPC/All.UnitLength_in_cm;
	All.deltaz *= CM_PER_KPC/All.UnitLength_in_cm;
	All.xx0 *= CM_PER_KPC/All.UnitLength_in_cm;
	All.yy0 *= CM_PER_KPC/All.UnitLength_in_cm;
	All.zz0 *= CM_PER_KPC/All.UnitLength_in_cm; 
	
	printf("Step and first values: %e %e %e %e %e %e \n", All.deltax,All.deltay,All.deltaz,All.xx0,All.yy0,All.zz0);
	
	for(int k=0; k < Nz; k++)
	{
		for(int j=0; j < Ny; j++)
		{
			for(int i=0; i < Nx; i++)
			{
				potential_tot[i-1][j][k]*=1e10/(All.UnitVelocity_in_cm_per_s * All.UnitVelocity_in_cm_per_s);
				potential_bar[i-1][j][k]*=1e10/(All.UnitVelocity_in_cm_per_s * All.UnitVelocity_in_cm_per_s);
			}
		}
	}
	
	//////////////////////////
	
	for(int k=0; k < Nz; k++)
	{
		for(int j=0; j < Ny; j++)
		{
			for(int i=1; i < Nx-1; i++)
			{
				All.accx[i][j][k] = -(potential_tot[i+1][j][k]-potential_tot[i-1][j][k])/(2*All.deltax);
				All.accx_bar[i][j][k] = -(potential_bar[i+1][j][k]-potential_bar[i-1][j][k])/(2*All.deltax);
			}
		}
	}
	
	for(int k=0; k < Nz; k++)
	{
		for(int j=0; j < Ny; j++)
		{
			All.accx[0][j][k] = -(potential_tot[1][j][k]-potential_tot[0][j][k])/(All.deltax);
			All.accx_bar[0][j][k] = -(potential_bar[1][j][k]-potential_bar[0][j][k])/(All.deltax);
			All.accx[Nx-1][j][k] = -(potential_tot[Nx-1][j][k]-potential_tot[Nx-2][j][k])/(All.deltax);
			All.accx_bar[Nx-1][j][k] = -(potential_bar[Nx-1][j][k]-potential_bar[Nx-2][j][k])/(All.deltax);
		}
	}
	
	
	for(int k=0; k < Nz; k++)
	{
		for(int j=1; j < Ny-1; j++)
		{
			for(int i=0; i < Nx; i++)
			{
				All.accy[i][j][k] = -(potential_tot[i][j+1][k]-potential_tot[i][j-1][k])/(2*All.deltay);
				All.accy_bar[i][j][k] = -(potential_bar[i][j+1][k]-potential_bar[i][j-1][k])/(2*All.deltay);
			}
		}
	}
	
	for(int k=0; k < Nz; k++)
	{		
		for(int i=0; i < Nx; i++)
		{
			All.accy[i][0][k] = -(potential_tot[i][1][k]-potential_tot[i][0][k])/(All.deltay);
			All.accy_bar[i][0][k] = -(potential_bar[i][1][k]-potential_bar[i][0][k])/(All.deltay);
			All.accy[i][Ny-1][k] = -(potential_tot[i][Ny-1][k]-potential_tot[i][Ny-2][k])/(All.deltay);
			All.accy_bar[i][Ny-1][k] = -(potential_bar[i][Ny-1][k]-potential_bar[i][Ny-2][k])/(All.deltay);
		}
	}
	
	for(int k=1; k < Nz-1; k++)
	{
		for(int j=0; j < Ny; j++)
		{
			for(int i=0; i < Nx; i++)
			{
				All.accz[i][j][k] = -(potential_tot[i][j][k+1]-potential_tot[i][j][k-1])/(2*All.deltaz);
				All.accz_bar[i][j][k] = -(potential_bar[i][j][k+1]-potential_bar[i][j][k-1])/(2*All.deltaz);
			}
		}
	}
	
    for(int j=0; j < Ny; j++)
	{
		for(int i=0; i < Nx; i++)
		{
			All.accz[i][j][0] = -(potential_tot[i][j][1]-potential_tot[i][j][0])/(All.deltaz);
			All.accz_bar[i][j][0] = -(potential_bar[i][j][1]-potential_bar[i][j][0])/(All.deltaz);
			All.accz[i][j][Nz-1] = -(potential_tot[i][j][Nz-1]-potential_tot[i][j][Nz-2])/(All.deltaz);
			All.accz_bar[i][j][Nz-1] = -(potential_bar[i][j][Nz-1]-potential_bar[i][j][Nz-2])/(All.deltaz);
		}
	}
	
	//printf("Accel: %e %e %e %e %e %e \n", All.accx[0][6][6], All.accx[1][6][6], All.accx[2][6][6], potential_tot[0][6][6], potential_tot[1][6][6], potential_tot[2][6][6]);
	
		
	free(XX);
	free(YY);
	free(ZZ);
	free(potential_tot);
	free(potential_bar);
}
#endif