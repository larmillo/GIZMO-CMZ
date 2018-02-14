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
    
    echo "Compiling grackle and slug for ${SYSTEM}"
elif  [[ ${SYSTEM} == "avatar" ]]
then 
    echo "Compiling grackle and slug for ${SYSTEM}"
elif  [[ ${SYSTEM} == "darwin" ]]
then 
    echo "Compiling grackle and slug for ${SYSTEM}"
else
    echo -e "Usage: ./compileall.sh machinename \n\nmachinename=marconi,avatar,darwin"
fi

GIZMO_DIR=${PWD}
GRACKLE_DIR=${GIZMO_DIR}/grackle
SLUG_DIR=${GIZMO_DIR}/slug2

# Compile grackle
cd ${GRACKLE_DIR}
./configure 
cd src/clib 
make machine-${SYSTEM}
make clean 
make -j 8 
make install
make clean

# Compile slug
cd ${SLUG_DIR}
make clean
make lib -j 8 MPI=ENABLE_MPI GSLVERSION=2 MACHINE=${SYSTEM} FITS=ENABLE_FITS
mv src/libslug.* .
make clean

cd ${GIZMO_DIR}
rm -rf lib
ln -s ${SLUG_DIR}/lib lib
