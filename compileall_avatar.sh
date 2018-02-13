module load intel
module load intelmpi
module load python
module load boost
module load gsl
module load fftw
module load szip
module load zlib/1.2.8--gnu--6.1.0
module load hdf5
module load profile/astro
module load cfitsio

# Compile grackle
cd /avatar/larmillo/GIZMO-CMZ/grackle 
./configure 
cd src/clib 
make machine-avatar
make clean 
make -j 24 
make install

# Compile slug
cd /avatar/larmillo/GIZMO-CMZ/slug2
make clean
make lib -j 24 MPI=ENABLE_MPI GSLVERSION=2 MACHINE=avatar FITS=ENABLE_FITS
