#
# This example writes data to the existing empty dataset created by h5_crtdat.py and then reads it back.
#
import h5py
import numpy as np
#
# Simulation parameters
#
NX = 100
NY = 100
NZ = 100
N = NX*NY*NZ
DIM = 3
xmax = 800.*3.08e18
xmin = -200.*3.08e18
deltax = (xmax - xmin)/NX
ymax = 500.*3.08e18
ymin = -500.*3.08e18
deltay = (ymax - ymin)/NY
zmax = 500.*3.08e18
zmin = -500.*3.08e18
deltaz = (zmax - zmin)/NZ
cloud_radius = 70.*3.08e18
rho1 = 0.2
rho2 = 0.001
T1 = 1.e4
T2 = 2.e6
vx = 75.e5
mu = 0.6
mh = 1.66e-24 
kb = 1.38e-16
Gamma_1 = 1.67 - 1.
#
# Open an existing file using default properties.
#
file = h5py.File('GIZMO.hdf5')
#
# Create string attribute.
#
group1 = file.create_group("Header")
#
# Create integer array attribute.
#
attr_data = 1000

attr_data1 = np.zeros((6))
for i in range(6):
    attr_data1[i] = 0
    
attr_data2 = np.zeros((6))
for i in range(6):
    attr_data2[i] = 0
attr_data2[0] = N    
    
attr_data3 = np.zeros((6))
for i in range(6):
    attr_data3[i] = 0        
#
#
#group1.attrs.create("BoxSizes", attr_data, h5py.h5t.IEEE_F64LE)
#group1.attrs.create("Flag_Cooling", 0, h5py.h5t.STD_I32LE)
#group1.attrs.create("Flag_DoublePrecision", 1, (1,), h5py.h5t.STD_I32LE)
#group1.attrs.create("Flag_Feedback", 0, (1,), h5py.h5t.STD_I32LE)
#group1.attrs.create("Flag_IC_Info", 0, (1,), h5py.h5t.STD_I32LE)
#group1.attrs.create("Flag_Metals", 0, (1,), h5py.h5t.STD_I32LE)
#group1.attrs.create("Flag_Sfr", 0, (1,), h5py.h5t.STD_I32LE)
#group1.attrs.create("Flag_StellarAge", 0, (1,), h5py.h5t.STD_I32LE)
#group1.attrs.create("HubbleParam", 1, (1,), h5py.h5t.IEEE_F64LE)
group1.attrs.create("MassTable", attr_data1, (6,), h5py.h5t.IEEE_F64LE)
#group1.attrs.create("NumFilesPerSnapshot", 1, (1,), h5py.h5t.STD_I32LE)
group1.attrs.create("NumPart_ThisFile", attr_data2, (6,), h5py.h5t.STD_I32LE)
group1.attrs.create("NumPart_Total", attr_data2, (6,), h5py.h5t.STD_U32LE)
group1.attrs.create("NumPart_Total_HighWord", attr_data3, (6,), h5py.h5t.STD_U32LE)
#group1.attrs.create("Omega0", 0, (1,), h5py.h5t.IEEE_F64LE)
#group1.attrs.create("OmegaLambda", 0, (1,), h5py.h5t.IEEE_F64LE)
#group1.attrs.create("Redshift", 0, (1,), h5py.h5t.IEEE_F64LE)
#group1.attrs.create("Time", 0, (1,), h5py.h5t.IEEE_F64LE)
group1.attrs['Time'] = 0.0  # initial time
group1.attrs['Redshift'] = 0.0 # initial redshift
group1.attrs['BoxSize'] = xmax -xmin # box size
group1.attrs['NumFilesPerSnapshot'] = 1 # number of files for multi-part snapshots
group1.attrs['Omega0'] = 1.0 # z=0 Omega_matter
group1.attrs['OmegaLambda'] = 0.0 # z=0 Omega_Lambda
group1.attrs['HubbleParam'] = 1.0 # z=0 hubble parameter (small 'h'=H/100 km/s/Mpc)
group1.attrs['Flag_Sfr'] = 0 # flag indicating whether star formation is on or off
group1.attrs['Flag_Cooling'] = 0 # flag indicating whether cooling is on or off
group1.attrs['Flag_StellarAge'] = 0 # flag indicating whether stellar ages are to be saved
group1.attrs['Flag_Metals'] = 0 # flag indicating whether metallicity are to be saved
group1.attrs['Flag_Feedback'] = 0 # flag indicating whether some parts of springel-hernquist model are active
group1.attrs['Flag_DoublePrecision'] = 0; # flag indicating whether ICs are in single/double precision
group1.attrs['Flag_IC_Info'] = 0; # flag indicating extra options for ICs
#
# Open "dset" dataset under the root group.
#
group = file.create_group("PartType0")
dataset1 = file.create_dataset("/PartType0/Coordinates", (N,DIM), dtype = h5py.h5t.IEEE_F64LE) 
dataset2 = file.create_dataset("/PartType0/Velocities", (N,DIM), dtype = h5py.h5t.IEEE_F64LE) 
dataset3 = file.create_dataset("/PartType0/Masses", (N,), dtype = h5py.h5t.IEEE_F64LE) 
dataset4 = file.create_dataset("/PartType0/Density", (N,), dtype = h5py.h5t.IEEE_F64LE) 
dataset5 = file.create_dataset("/PartType0/InternalEnergy", (N,), dtype = h5py.h5t.IEEE_F64LE) 
dataset6 = file.create_dataset("/PartType0/ParticleIDs", (N,), dtype = h5py.h5t.STD_U32LE)
#
# Initialize data object with 0.
#
Coordinates = np.zeros((N,DIM))
Velocities = np.zeros((N,DIM))
Masses = np.zeros(N)
Density = np.zeros(N)
InternalEnergy = np.zeros(N)
ParticleIDs = np.zeros(N)
#
# Assign new values
#
for z in range(NZ):
    for y in range(NY):
        for x in range(NX):
            i = x+y*NX+z*NX*NY; 
            Coordinates[i][0] = (x+1./2.)*deltax	
            Coordinates[i][1] = (y+1./2.)*deltay	
            Coordinates[i][2] = (z+1./2.)*deltaz
            ParticleIDs[i] = i; 
    
for i in range(N):
    xx = Coordinates[i][0] + xmin
    yy = Coordinates[i][1] + ymin
    zz = Coordinates[i][2] + zmin
    radius = np.sqrt(xx*xx+yy*yy+zz*zz)
    dist = (radius-cloud_radius)   
    if (dist < 0.):
        Velocities[i][0] = vx
        Velocities[i][1] = 0.
        Velocities[i][2] = 0.
        Density[i] = rho1 * mu * mh
        InternalEnergy[i] = rho1 * T1 * kb / (Gamma_1 * Density[i])
        Masses[i] = Density[i] * deltax**3
    else:
        Velocities[i][0] = 0.
        Velocities[i][1] = 0.
        Velocities[i][2] = 0.
        Density[i] = rho2 * mu * mh
        InternalEnergy[i] = rho2 * T2 * kb / (Gamma_1 * Density[i])
        Masses[i] = Density[i] * deltax**3
print rho2 * mu * mh, rho1 * mu * mh        
#
# Write data
#
print "Writing data..."
dataset1[...] = Coordinates
dataset2[...] = Velocities
dataset3[...] = Masses
dataset4[...] = Density
dataset5[...] = InternalEnergy
dataset6[...] = ParticleIDs
#
# Close the file before exiting
#
file.close()
