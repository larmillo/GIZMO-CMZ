SYSTEM=${1}


if [[ ${SYSTEM} == "marconi" ]]
then
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

    GRACKLE_DIR=${HOME}/GIZMO-CMZ/grackle
    SLUG_DIR=${HOME}/GIZMO-CMZ/slug2
    
    echo Compiling grackle and slug for ${SYSTEM}
elif  [[ ${SYSTEM} == "darwin" ]]
    
fi


#GRACKLE_DIR=${PWD}/grackle
#SLUG_DIR=${PWD}/GIZMO-CMZ/slug2

# Compile grackle
#cd ${GRACKLE_DIR}
#./configure 
#cd src/clib 
#make machine-marconi
#make clean 
#make -j 24 
#make install

# Compile slug
#cd ${HOME}/GIZMO-CMZ/slug2
#make clean
#make lib -j 24 MPI=ENABLE_MPI GSLVERSION=2 MACHINE=marconi FITS=ENABLE_FITS
