#!/bin/bash            # this line only there to enable syntax highlighting in this file

####################################################################################################
#  Enable/Disable compile-time options as needed: this is where you determine how the code will act
#  From the list below, please activate/deactivate the
#       options that apply to your run. If you modify any of these options,
#       make sure that you recompile the whole code by typing "make clean; make".
#
# This file was originally part of the GADGET3 code developed by
#   Volker Springel (volker.springel@h-its.org). The code has been modified
#   substantially by Phil Hopkins (phopkins@caltech.edu) for GIZMO (to add new modules and clean
#   up the naming conventions and changed many of them to match the new GIZMO conventions)
#
####################################################################################################



####################################################################################################
# --------------------------------------- Boundary Conditions & Dimensions
####################################################################################################
#PERIODIC                        # Use this if periodic boundaries are needed (otherwise open boundaries are assumed)
#BND_PARTICLES                  # particles with ID=0 are forced in place (their accelerations are set =0):
                                # use for special boundary conditions where these particles represent fixed "walls"
#LONG_X=1                    # modify box dimensions (non-square periodic box): multiply X (PERIODIC and NOGRAVITY required)
#LONG_Y=1                      # modify box dimensions (non-square periodic box): multiply Y
#LONG_Z=3                       # modify box dimensions (non-square periodic box): multiply Z
#REFLECT_BND_X                  # make the x-boundary reflecting (assumes a box 0<x<1, unless PERIODIC is set)
#REFLECT_BND_Y                  # make the y-boundary reflecting (assumes a box 0<y<1, unless PERIODIC is set)
#REFLECT_BND_Z                  # make the z-boundary reflecting (assumes a box 0<z<1, unless PERIODIC is set)
#SHEARING_BOX=1                 # shearing box boundaries: 1=r-z sheet (r,z,phi coordinates), 2=r-phi sheet (r,phi,z), 3=r-phi-z box, 4=as 3, with vertical gravity
#SHEARING_BOX_Q=(3./2.)         # shearing box q=-dlnOmega/dlnr; will default to 3/2 (Keplerian) if not set
#ONEDIM                         # Switch for 1D test problems: code only follows the x-line. requires NOGRAVITY, and all y=z=0
#TWODIMS                        # Switch for 2D test problems: code only follows the xy-plane. requires NOGRAVITY, and all z=0.
####################################################################################################



####################################################################################################
# --------------------------------------- Hydro solver method
####################################################################################################
HYDRO_MESHLESS_FINITE_MASS      # Lagrangian (constant-mass) finite-volume Godunov method
#HYDRO_MESHLESS_FINITE_VOLUME   # Moving (quasi-Lagrangian) finite-volume Godunov method
## -----------------------------------------------------------------------------------------------------
# --------------------------------------- SPH methods:
#SPHEQ_DENSITY_INDEPENDENT_SPH  # force SPH to use the 'pressure-sph' formulation ("modern" SPH)
#SPHEQ_TRADITIONAL_SPH          # force SPH to use the 'density-sph' (GADGET-2 & GASOLINE SPH)
# --------------------------------------- SPH artificial diffusion options (use with SPH; not relevant for Godunov/Mesh modes)
#SPHAV_DISABLE_CD10_ARTVISC     # Disable Cullen & Dehnen 2010 'inviscid sph' (viscosity suppression outside shocks); just use Balsara switch
#SPHAV_DISABLE_PM_CONDUCTIVITY  # Disable mixing entropy (J.Read's improved Price-Monaghan conductivity with Cullen-Dehnen switches)
## -----------------------------------------------------------------------------------------------------
# --------------------------------------- Kernel Options
#KERNEL_FUNCTION=3              # Choose the kernel function (2=quadratic peak, 3=cubic spline [default], 4=quartic spline, 5=quintic spline, 6=Wendland C2, 7=Wendland C4, 8=2-part quadratic)
####################################################################################################



####################################################################################################
# --------------------------------------- Additional Fluid Physics
####################################################################################################
EOS_GAMMA=(5.0/3.0)            # Polytropic Index of Gas (for an ideal gas law): if not set and no other (more complex) EOS set, defaults to GAMMA=5/3
#EOS_HELMHOLTZ                  # Use Timmes & Swesty 2000 EOS (for e.g. stellar or degenerate equations of state)
## -----------------------------------------------------------------------------------------------------
# --------------------------------- Magneto-Hydrodynamics
# ---------------------------------  these modules are public, but if used, the user should also cite the MHD-specific GIZMO methods paper
# ---------------------------------  (Hopkins 2015: 'Accurate, Meshless Methods for Magneto-Hydrodynamics') as well as the standard GIZMO paper
#MAGNETIC                       # master switch for MHD, regardless of which Hydro solver is used
#B_SET_IN_PARAMS                # set initial fields (Bx,By,Bz) in parameter file
#MHD_NON_IDEAL                  # enable non-ideal MHD terms: Ohmic resistivity, Hall effect, and ambipolar diffusion (solved explicitly)
#CONSTRAINED_GRADIENT_MHD=1     # use CG method (in addition to cleaning, optional!) to maintain low divB: set this value to control how aggressive the div-reduction is:
                                # 0=minimal (safest), 1=intermediate (recommended), 2=aggressive (less stable), 3+=very aggressive (less stable+more expensive)
#CONDUCTION                     # Thermal conduction solved *explicitly*: isotropic if MAGNETIC off, otherwise anisotropic
#CONDUCTION_SPITZER             # Spitzer conductivity accounting for saturation: otherwise conduction coefficient is constant
#VISCOSITY                      # Navier-stokes equations solved *explicitly*: isotropic coefficients if MAGNETIC off, otherwise anisotropic
#VISCOSITY_BRAGINSKII           # Braginskii viscosity tensor for ideal MHD
#TURB_DIFF_ENERGY               # turbulent diffusion of internal energy (conduction with effective turbulent coefficients)
#TURB_DIFF_VELOCITY             # turbulent diffusion of momentum (viscosity with effective turbulent coefficients)

#--------------------------------------- Radiative Cooling physics (mostly geared towards galactic/extragalactic cooling)
#--------------------------- These modules were originally developed for a combination of -proprietary- physics modules. they can only be used with
#--------------------------- permission from the authors. email P. Hopkins to obtain the relevant permissions for the cooling routines of interest.
COOLING                        # enables radiative cooling and heating: if GALSF, also external UV background read from file "TREECOOL"
GRACKLE                        # enable GRACKLE: cooling+chemistry package (requires COOLING above; https://grackle.readthedocs.org/en/latest )
GRACKLE_CHEMISTRY=0            # choose GRACKLE cooling chemistry: (0)=tabular, (1)=Atomic, (2)=(1)+H2+H2I+H2II, (3)=(2)+DI+DII+HD
GRACKLE_OPTS			# additional options for GRACKLE chemistry solver (Lupi)
#GRACKLE_FIX_TEMPERATURE	# if the initial temperature is provided via parameter file, the intial internal energy is computed taking into account the mmw dependence on density too (Lupi)
#GRACKLE_FULLYIMPLICIT		# fully implicit solution for Grackle (customized version of the library) (Lupi)
##-----------------------------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------------------------
####################################################################################################

## ----------------------------------------------------------------------------------------------------
# --------------------------------------- Aerodynamic Particles
# --------------------------------------- (this is developed by P. Hopkins as part of the FIRE package: the same FIRE authorship & approval policies apply, see below)
#GRAIN_FLUID                    # aerodynamically-coupled grains (particle type 3 are grains); default is Stokes drag
#GRAIN_EPSTEIN=1                # uses the cross section for molecular hydrogen (times this number) to calculate Epstein drag
#GRAIN_LORENTZFORCE             # charged grains feel Lorentz forces (requires MAGNETIC)
#GENTRY_FB
####################################################################################################



####################################################################################################
# ------------------------------------- Driven turbulence (for turbulence tests, large-eddy sims)
# ------------------------------- users of these routines should cite Bauer & Springel 2012, MNRAS, 423, 3102, and thank A. Bauer for providing the core algorithms
####################################################################################################
#TURB_DRIVING                   # turns on turbulent driving/stirring. see begrun for parameters that must be set
#POWERSPEC_GRID=128             # activates on-the-fly calculation of the turbulent velocity, vorticity, and smoothed-velocity power spectra
#ADJ_BOX_POWERSPEC              # compiles in a code module that allows via restart-flag 6 the calculation of a gas velocity power spectrum of a snapshot with an adjustable box (user defined center and size)
####################################################################################################



####################################################################################################
## ------------------------ Gravity & Cosmological Integration Options ---------------------------------
####################################################################################################
# --------------------------------------- TreePM Options (recommended for cosmological sims)
#PMGRID=512                     # COSMO enable: resolution of particle-mesh grid
#PM_PLACEHIGHRESREGION=1+2+16   # COSMO enable: particle types to place high-res PMGRID around
#PM_HIRES_REGION_CLIPPING=1000  # for stability: clips particles that escape the hires region in zoom/isolated sims
#PM_HIRES_REGION_CLIPDM         # split low-res DM particles that enter high-res region (completely surrounded by high-res)
#MULTIPLEDOMAINS=64             # Multi-Domain option for the top-tree level: iso=16,COSMO=64-128
## -----------------------------------------------------------------------------------------------------
# ---------------------------------------- Adaptive Grav. Softening (including Lagrangian conservation terms!)
#ADAPTIVE_GRAVSOFT_FORGAS       # allows variable softening length (=Hsml) for gas particles
#ADAPTIVE_GRAVSOFT_FORALL=100   # enable adaptive gravitational softening lengths for all particle types
                                # (ADAPTIVE_GRAVSOFT_FORGAS should be disabled). the softening is set to the distance
                                # enclosing a neighbor number set in the parameter file. baryons search for other baryons,
                                # dm for dm, sidm for sidm, etc. If set to numerical value, the maximum softening is this times All.ForceSoftening[for appropriate particle type]
## -----------------------------------------------------------------------------------------------------
NOGRAVITY                      # turn off self-gravity (compatible with analytic_gravity)
#GRAVITY_NOT_PERIODIC           # self-gravity is not periodic, even though the rest of the box is periodic
## -----------------------------------------------------------------------------------------------------
#ANALYTIC_GRAVITY               # Specific analytic gravitational force to use instead of/with self-gravity. If set to a numerical value
                                #  (edit "gravity/analytic_gravity.h" to actually assign the analytic gravitational forces)
                                # with a bit mask, as for PM_PLACEHIGHRESREGION above (see description)
#EOS_TRUELOVE_PRESSURE          # adds artificial pressure floor force Jeans length above resolution scale (means you will get the wrong answer, but things will look smooth)
####################################################################################################



####################################################################################################
# --------------------------------------- On the fly FOF groupfinder
# ----------------- This is originally developed as part of GADGET-3 by V. Springel
####################################################################################################
#FOF_PRIMARY_LINK_TYPES=2           # 2^type for the primary dark matter type
#FOF_SECONDARY_LINK_TYPES=1+16+32   # 2^type for the types linked to nearest primaries
#DENSITY_SPLIT_BY_TYPE=1+2+16+32    # 2^type for whch the densities should be calculated seperately
#FOF_GROUP_MIN_LEN=32               # default is 32
####################################################################################################






##BH_WIND_SPAWN
##BH_SEED_GROWTH_TESTS          # Currently testing options for BH seeding 






#RT_LYMAN_WERNER                        # specific lyman-werner [narrow H2 dissociating] band



####################################################################################################
# --------------------------------------- Multi-Threading (parallelization) options
####################################################################################################
#OPENMP=2                       # Masterswitch for explicit OpenMP implementation
#PTHREADS_NUM_THREADS=4         # custom PTHREADs implementation (don't enable with OPENMP)
####################################################################################################



####################################################################################################
# --------------------------------------- Output/Input options
####################################################################################################
HAVE_HDF5						# needed when HDF5 I/O support is desired
#OUTPUT_IN_DOUBLEPRECISION      # snapshot files will be written in double precision
#INPUT_IN_DOUBLEPRECISION       # input files assumed to be in double precision (otherwise float is assumed)
#OUTPUT_POSITIONS_IN_DOUBLE     # input/output files in single, but positions in double (used in hires, hi-dynamic range sims when positions differ by < float accuracy)
#INPUT_POSITIONS_IN_DOUBLE      # as above, but specific to the ICs file
#OUTPUTPOTENTIAL                # forces code to compute+output potentials in snapshots
#OUTPUTACCELERATION             # output physical acceleration of each particle in snapshots
#OUTPUTCHANGEOFENERGY           # outputs rate-of-change of internal energy of gas particles in snapshots
#OUTPUT_VORTICITY				# outputs the vorticity vector
#OUTPUTTIMESTEP                 # outputs timesteps for each particle
#OUTPUTCOOLRATE					# outputs cooling rate, and conduction rate if enabled
#POWERSPEC_ON_OUTPUT            # compute and output power spectra (not used)
#RECOMPUTE_POTENTIAL_ON_OUTPUT	# update potential every output even it EVALPOTENTIAL is set
#OUTPUT_ADDITIONAL_RUNINFO      # enables extended simulation output data (can slow down machines significantly in massively-parallel runs)
####################################################################################################



####################################################################################################
# -------------------------------------------- De-Bugging & special (usually test-problem only) behaviors
####################################################################################################
#DEVELOPER_MODE                 # allows you to modify various numerical parameters (courant factor, etc) at run-time
#EOS_ENFORCE_ADIABAT=(1.0)      # if set, this forces gas to lie -exactly- along the adiabat P=EOS_ENFORCE_ADIABAT*(rho^GAMMA)
#SLOPE_LIMITER_TOLERANCE=1      # sets the slope-limiters used. higher=more aggressive (less diffusive, but less stable). 1=default. 0=conservative. use on problems where sharp density contrasts in poor particle arrangement may cause errors. 2=same as AGGRESSIVE_SLOPE_LIMITERS below
#AGGRESSIVE_SLOPE_LIMITERS      # use the original GIZMO paper (more aggressive) slope-limiters. more accurate for smooth problems, but
                                # these can introduce numerical instability in problems with poorly-resolved large noise or density contrasts (e.g. multi-phase, self-gravitating flows)
#ENERGY_ENTROPY_SWITCH_IS_ACTIVE # enable energy-entropy switch as described in GIZMO methods paper. This can greatly improve performance on some problems where the
                                # the flow is very cold and highly super-sonic. it can cause problems in multi-phase flows with strong cooling, though, and is not compatible with non-barytropic equations of state
#FORCE_ENTROPIC_EOS_BELOW=(0.01) # set (manually) the alternative energy-entropy switch which is enabled by default in MFM/MFV: if relative velocities are below this threshold, it uses the entropic EOS
#TEST_FOR_IDUNIQUENESS          # explicitly check if particles have unique id numbers (only use for special behaviors)
#LONGIDS                        # use long ints for IDs (needed for super-large simulations)
#ASSIGN_NEW_IDS                 # assign IDs on startup instead of reading from ICs
#NO_CHILD_IDS_IN_ICS            # IC file does not have child IDs: do not read them (used for compatibility with snapshot restarts from old versions of the code)
#READ_HSML                      # reads hsml from IC file
#PREVENT_PARTICLE_MERGE_SPLIT   # don't allow gas particle splitting/merging operations
#PARTICLE_EXCISION              # enable dynamical excision (remove particles within some radius)
#MERGESPLIT_HARDCODE_MAX_MASS=(1.0e-6)   # manually set maximum mass for particle merge-split operations (in code units): useful for snapshot restarts and other special circumstances
#MERGESPLIT_HARDCODE_MIN_MASS=(1.0e-7)   # manually set minimum mass for particle merge-split operations (in code units): useful for snapshot restarts and other special circumstances


#USE_MPI_IN_PLACE               # MPI debugging: makes AllGatherV compatible with MPI_IN_PLACE definitions in some MPI libraries
#NO_ISEND_IRECV_IN_DOMAIN       # MPI debugging: slower, but fixes memory errors during exchange in the domain decomposition (ANY RUN with >2e9 particles MUST SET THIS OR FAIL!)
#FIX_PATHSCALE_MPI_STATUS_IGNORE_BUG # MPI debugging
#MPISENDRECV_SIZELIMIT=100      # MPI debugging
#MPISENDRECV_CHECKSUM           # MPI debugging
#DONOTUSENODELIST               # MPI debugging
#NOTYPEPREFIX_FFTW              # FFTW debugging (fftw-header/libraries accessed without type prefix, adopting whatever was
                                #   chosen as default at compile of fftw). Otherwise, the type prefix 'd' for double is used.
#DOUBLEPRECISION_FFTW           # FFTW in double precision to match libraries
#DEBUG                          # enables core-dumps and FPU exceptions
#STOP_WHEN_BELOW_MINTIMESTEP    # forces code to quit when stepsize wants to go below MinSizeTimestep specified in the parameterfile
#SEPARATE_STELLARDOMAINDECOMP   # separate stars (ptype=4) and other non-gas particles in domain decomposition (may help load-balancing)
#DISABLE_SPH_PARTICLE_WAKEUP    # don't let gas particles move to lower timesteps based on neighbor activity (use for debugging)
#EVALPOTENTIAL                  # computes gravitational potential
#MHD_ALTERNATIVE_LEAPFROG_SCHEME # use alternative leapfrog where magnetic fields are treated like potential/positions (per Federico Stasyszyn's suggestion): still testing
#FREEZE_HYDRO                   # zeros all fluxes from RP and doesn't let particles move (for testing additional physics layers)
#SUPER_TIMESTEP_DIFFUSION       # use super-timestepping to accelerate integration of diffusion operators [for testing or if there are stability concerns]
#ALLOW_IMBALANCED_GASPARTICLELOAD # increases All.MaxPartSph to All.MaxPart: can allow better load-balancing in some cases, but uses more memory. But use me if you run into errors where it can't fit the domain (where you would increase PartAllocFac, but can't for some reason)
####################################################################################################


#NOTEST_FOR_IDUNIQUENESS













