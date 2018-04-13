SYSTEM=${1}

NPROC=24
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
 
elif [[ ${SYSTEM} == "raijin" ]]
then
    module unload intel-cc intel-fc
    module load intel-cc/17.0.1.132
    module load intel-fc/17.0.1.132
    #module load gcc/4.9.0
    module load openmpi/2.1.1
    module load python
    module load boost/1.66.0
    module load gsl/2.4
    module load fftw3/3.3.5 
    module load szip
    module load zlib
    module load cfitsio/3.350
    
elif  [[ ${SYSTEM} == "avatar" ]]
then 
    :
elif  [[ ${SYSTEM} == "darwin" ]]
then 
    NPROC=8
else
    echo -e "Usage: ./compileall.sh machinename \n\nmachinename=marconi,raijin,avatar,darwin"
    exit
fi

GIZMO_DIR=${PWD}
GRACKLE_DIR=${GIZMO_DIR}/grackle
SLUG_DIR=${GIZMO_DIR}/slug2

echo "Compiling grackle and slug for ${SYSTEM}"

# Compile grackle
cd ${GRACKLE_DIR}
./configure 
cd src/clib 
make machine-${SYSTEM}
make clean 
make -j ${NPROC}
make install
make clean

# Compile slug
cd ${SLUG_DIR}
make clean
make lib -j ${NPROC} MPI=ENABLE_MPI GSLVERSION=2 MACHINE=${SYSTEM} FITS=ENABLE_FITS
mv src/libslug.* .
make clean

cd ${GIZMO_DIR}
rm -rf lib
ln -s ${SLUG_DIR}/lib lib
