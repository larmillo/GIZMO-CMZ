/*! \file analytic_gravity.h
 *  \brief externally-specified (analytic) gravity goes here
 *
 *  This file contains supplemental code if you want to add an 
 *   -analytic- potential or gravitational force in the code, 
 *   rather than solely relying on the calculated self-gravity
 */
/*
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

void add_analytic_gravitational_forces(void);
void GravAccel_StaticPlummerSphere(void);
void GravAccel_StaticHernquist(void);
void GravAccel_StaticIsothermalSphere(void);
void GravAccel_KeplerianOrbit(void);
void GravAccel_KeplerianTestProblem(void);
void GravAccel_GrowingDiskPotential(void);
void GravAccel_StaticNFW(void);
void GravAccel_RayleighTaylorTest(void);
void GravAccel_ShearingSheet(void);
void GravAccel_PaczynskyWiita(void);
void GravAccel_CMZ(void);


/* master routine which decides which (if any) analytic gravitational forces are applied */
void add_analytic_gravitational_forces()
{
#ifdef ANALYTIC_GRAVITY
    //GravAccel_RayleighTaylorTest();     // vertical potential for RT tests
    //GravAccel_StaticPlummerSphere();    // plummer sphere
    //GravAccel_StaticHernquist();        // hernquist sphere
    //GravAccel_StaticIsothermalSphere(); // singular or cored isothermal sphere
    //GravAccel_KeplerianOrbit();         // keplerian disk
    //GravAccel_KeplerianTestProblem();   // keplerian disk with boundaries for test problem
    //GravAccel_GrowingDiskPotential();   // time-dependent (adiabatically growing) disk
    //GravAccel_StaticNFW();              // spherical NFW profile
    //GravAccel_PaczynskyWiita();         // Paczynsky-Wiita pseudo-Newtonian potential
	GravAccel_CMZ();                      // Milky-Way potential in the Central Molecular Zone
#ifdef SHEARING_BOX
    GravAccel_ShearingSheet();            // adds coriolis and centrifugal terms for shearing-sheet approximation
#endif
#endif
}

/* adds coriolis and centrifugal terms for shearing-sheet approximation */
void GravAccel_ShearingSheet()
{
#ifdef SHEARING_BOX
    int i;
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        /* centrifugal force term (depends on distance from box center) */
        P[i].GravAccel[0] += 2.*(P[i].Pos[0]-boxHalf_X) * SHEARING_BOX_Q*SHEARING_BOX_OMEGA_BOX_CENTER*SHEARING_BOX_OMEGA_BOX_CENTER;
        /* coriolis force terms */
        double vp=0;
        if(P[i].Type==0) {vp=SphP[i].VelPred[SHEARING_BOX_PHI_COORDINATE];} else {vp=P[i].Vel[SHEARING_BOX_PHI_COORDINATE];}
        P[i].GravAccel[0] += 2.*vp * SHEARING_BOX_OMEGA_BOX_CENTER;
        if(P[i].Type==0) {vp=SphP[i].VelPred[0];} else {vp=P[i].Vel[0];}
        P[i].GravAccel[SHEARING_BOX_PHI_COORDINATE] -= 2.*vp * SHEARING_BOX_OMEGA_BOX_CENTER;
#if (SHEARING_BOX==4)
        /* add vertical gravity to the force law */
        P[i].GravAccel[2] -= SHEARING_BOX_OMEGA_BOX_CENTER * SHEARING_BOX_OMEGA_BOX_CENTER * (P[i].Pos[2]-boxHalf_Z);
#endif
    }
#endif
}



/* constant vertical acceleration for Rayleigh-Taylor test problem */
void GravAccel_RayleighTaylorTest()
{
    int i;
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        /* zero out the gravity first (since this test doesn't use self-gravity) */
        P[i].GravAccel[0]=P[i].GravAccel[1]=P[i].GravAccel[2]=0;
        /* now add the constant vertical field */
        if(P[i].ID != 0) {P[i].GravAccel[1]=-0.5;}
    }
}



/* static unit Plummer Sphere (G=M=a=1) */
void GravAccel_StaticPlummerSphere()
{
    int i,k; double r, r2, dp[3];
    
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        dp[0]=P[i].Pos[0]; dp[1]=P[i].Pos[1]; dp[2]=P[i].Pos[2];
#ifdef ANALYTIC_GRAVITY_ANCHOR_TO_PARTICLE
        for(k = 0; k < 3; k++) {dp[k] = -P[i].min_xyz_to_bh[k];}
#endif
        r2 = dp[0]*dp[0] + dp[1]*dp[1] + dp[2]*dp[2];
        r = sqrt(r2);
        for(k = 0; k < 3; k++) {P[i].GravAccel[k] += -dp[k] / pow(r2 + 1, 1.5);}
        
    }
}



/* static Hernquist Profile (parameters specified in the routine below) */
void GravAccel_StaticHernquist()
{
    double HQ_M200 = 95.2401;
    double HQ_C = 9.0;
    double HQ_DARKFRACTION = 0.9;

    double r, r2, dp[3], m, a; int i, k;
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        dp[0]=P[i].Pos[0]; dp[1]=P[i].Pos[1]; dp[2]=P[i].Pos[2];
#ifdef ANALYTIC_GRAVITY_ANCHOR_TO_PARTICLE
        for(k = 0; k < 3; k++) {dp[k] = -P[i].min_xyz_to_bh[k];}
#endif
        r2 = dp[0]*dp[0] + dp[1]*dp[1] + dp[2]*dp[2]; r = sqrt(r2);
        a = pow(All.G * HQ_M200 / (100 * All.Hubble_H0_CodeUnits * All.Hubble_H0_CodeUnits), 1.0 / 3) / HQ_C * sqrt(2 * (log(1 + HQ_C) - HQ_C / (1 + HQ_C)));
        m = HQ_M200 * pow(r / (r + a), 2) * HQ_DARKFRACTION;
        if(r > 0)
        {
            for(k = 0; k < 3; k++) {P[i].GravAccel[k] += -All.G * m * dp[k] / (r * r * r);}
            
        }
    }
}



/* static Isothermal Sphere Profile (parameters specified in the routine below) */
void GravAccel_StaticIsothermalSphere()
{
    double ISO_M200=95.21;
    double ISO_R200=160.0;
    double ISO_Eps=0.1;
    double ISO_FRACTION=0.9;
    double r, r2, dp[3], m; int i, k;
    
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        dp[0]=P[i].Pos[0]; dp[1]=P[i].Pos[1]; dp[2]=P[i].Pos[2];
#ifdef ANALYTIC_GRAVITY_ANCHOR_TO_PARTICLE
        for(k = 0; k < 3; k++) {dp[k] = -P[i].min_xyz_to_bh[k];}
#endif
        r2 = dp[0]*dp[0] + dp[1]*dp[1] + dp[2]*dp[2]; r = sqrt(r2);
        m = ISO_M200 * ISO_FRACTION; if(r < ISO_R200) {m *= r/ISO_R200;}
        if(r > 0) {for(k=0;k<3;k++) {P[i].GravAccel[k] += -All.G * m * dp[k] / (r*r*r + ISO_Eps*ISO_Eps);}}
    }
}


/* time-dependent potential of an adiabatically-growing disk */
void GravAccel_GrowingDiskPotential()
{
    int n_table = 14; // number of table entries below (must match!)
    // scale factor for cosmological runs (must be in monotonic increasing order!)
    double t_disk_table[14] = {0.2, 0.250, 0.266, 0.285, 0.308, 0.333, 0.363, 0.400, 0.444, 0.500, 0.572, 0.667, 0.800, 1.000};
    // m12i parameters: from Shea's fits:
    double m_disk_table[14] = {0.0, 0.061, 0.088, 0.117, 0.153, 0.223, 0.348, 0.429, 0.581, 1.118, 2.004, 3.008, 4.403, 6.001}; // disk mass in code units
    double r_disk_table[14] = {1.0, 5.071, 7.513, 6.787, 6.162, 3.277, 4.772, 3.964, 3.418, 2.511, 2.463, 1.503, 1.005, 1.150}; // disk scale length in code units
    double z_disk_table[14] = {1.0, 4.185, 8.971, 5.089, 3.532, 3.057, 4.557, 2.117, 1.828, 0.809, 0.217, 0.148, 0.335, 0.404}; // disk scale height in code units
    /* before the particle loop, interpolate the relevant quantities to the simulation time */
    double t=All.Time, dt=0, r2, dp[3], Zterm, Rterm, Rterm2, myfacR; int i, i0=0, i1=0, k;
    if(t<=t_disk_table[0])
    {
        i0=i1=0;
    } else if(t>=t_disk_table[n_table-1]) {
        i0=i1=n_table-1;
    } else {
        for(k=1;k<n_table;k++) {if(t_disk_table[k] > t) {i1=k; break;}}
        i0=i1-1; dt=(t - t_disk_table[i0])/(t_disk_table[i1]-t_disk_table[i0]);
    }
    double m_disk = m_disk_table[i0] + dt * (m_disk_table[i1]-m_disk_table[i0]);
    double r_disk = r_disk_table[i0] + dt * (r_disk_table[i1]-r_disk_table[i0]);
    double z_disk = z_disk_table[i0] + dt * (z_disk_table[i1]-z_disk_table[i0]);
    /* ok now we can assign actual accelerations */
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        dp[0]=P[i].Pos[0]; dp[1]=P[i].Pos[1]; dp[2]=P[i].Pos[2];
#ifdef ANALYTIC_GRAVITY_ANCHOR_TO_PARTICLE
        for(k = 0; k < 3; k++) {dp[k] = -P[i].min_xyz_to_bh[k];}
#endif
        r2 = dp[0]*dp[0] + dp[1]*dp[1];
        Zterm = sqrt(z_disk*z_disk + dp[2]*dp[2]); /* sqrt((Zdisk^2 + dZ^2); appears several times  */
        Rterm = r_disk + Zterm; Rterm2 = sqrt(r2 + Rterm*Rterm); Rterm2 = Rterm2*Rterm2*Rterm2;
        myfacR = -All.G * m_disk / Rterm2; /* has units s^-2, so  multiply by length to get accel.  no sign; handle that in min_xyz_to_bh */
        /* remember, min_xyz_to_bh = x_BH - myx => positive if x_BH > myx => acceleration is in positive x if x_BH > myx, which is correct (attractive) */
        P[i].GravAccel[0] += myfacR * dp[0]; P[i].GravAccel[1] += myfacR * dp[1];
        P[i].GravAccel[2] += myfacR * dp[2] * Rterm/Zterm; // this has units of:  M*L^3*M^-1*T^-2*L^2*L^-1*L^-3 = L/T^2
    }
}


/* Keplerian forces (G=M=1): useful for orbit, MRI, planetary disk problems */
void GravAccel_KeplerianOrbit()
{
    double dp[3], r, r2; int i;
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        dp[0]=P[i].Pos[0]; dp[1]=P[i].Pos[1]; dp[2]=P[i].Pos[2];
#ifdef ANALYTIC_GRAVITY_ANCHOR_TO_PARTICLE
        int k; for(k = 0; k < 3; k++) {dp[k] = -P[i].min_xyz_to_bh[k];}
#endif
#if defined(PERIODIC)
        dp[0] -= boxHalf_X; dp[1] -= boxHalf_Y;
#endif
        r2 = dp[0]*dp[0] + dp[1]*dp[1]; r = sqrt(r2);
        P[i].GravAccel[0] = -dp[0] / (r2 * r);
        P[i].GravAccel[0] = -dp[1] / (r2 * r);
        P[i].GravAccel[2] = 0;
    }
}





/* Keplerian forces (G=M=1): this is a specific (bounded and softened) version 
 used just for the Keplerian disk test problem */
void GravAccel_KeplerianTestProblem()
{
    double x00=0;//boxHalf_X;
    double y00=0;//boxHalf_Y;
    x00=4.0;
    y00=4.0;
    int i;
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        double r = pow(pow(P[i].Pos[1]-y00,2.)+pow(P[i].Pos[0]-x00,2.),0.5);
        if((r > 0.35)&(r < 2.1))
        {
            P[i].GravAccel[0] = -(P[i].Pos[0]-x00) / pow(pow(P[i].Pos[1]-y00,2.)+pow(P[i].Pos[0]-x00,2.),1.5) ;
            P[i].GravAccel[1] = -(P[i].Pos[1]-y00) / pow(pow(P[i].Pos[1]-y00,2.)+pow(P[i].Pos[0]-x00,2.),1.5) ;
            P[i].GravAccel[2] = 0;
        }
        if(r <= 0.35)
        {
            P[i].GravAccel[0] = -(P[i].Pos[0]-x00)*pow(r/0.35,2) / pow(pow(P[i].Pos[1]-y00,2.)+pow(P[i].Pos[0]-x00,2.),1.5) ;
            P[i].GravAccel[1] = -(P[i].Pos[1]-y00)*pow(r/0.35,2) / pow(pow(P[i].Pos[1]-y00,2.)+pow(P[i].Pos[0]-x00,2.),1.5) ;
            
            P[i].GravAccel[0] += +(P[i].Pos[0]-x00)*(0.35-r)/0.35 / pow(pow(P[i].Pos[1]-y00,2.)+pow(P[i].Pos[0]-x00,2.),1.5) ;
            P[i].GravAccel[1] += +(P[i].Pos[1]-y00)*(0.35-r)/0.35 / pow(pow(P[i].Pos[1]-y00,2.)+pow(P[i].Pos[0]-x00,2.),1.5) ;
            P[i].GravAccel[2] = 0;
        }
        if(r >= 2.1)
        {
            P[i].GravAccel[0] = -(P[i].Pos[0]-x00)*(1+(r-2.1)/0.1) / pow(pow(P[i].Pos[1]-y00,2.)+pow(P[i].Pos[0]-x00,2.),1.5) ;
            P[i].GravAccel[1] = -(P[i].Pos[1]-y00)*(1+(r-2.1)/0.1) / pow(pow(P[i].Pos[1]-y00,2.)+pow(P[i].Pos[0]-x00,2.),1.5) ;
            P[i].GravAccel[2] = 0;
        }
    }
}
#ifdef ANALYTIC_GRAVITY
void GravAccel_CMZ()
{
	double dp[3], interp1_x, interp2_x, interp3_x, interp4_x, interp1_y, interp2_y, interp1_z, interp2_z;
	double abar_x, abar_y, abar_tot, theta, omega, vrot, R;
	int x,y,z,i;
	double t=All.Time;
	
	vrot = 40*(1e5/All.UnitVelocity_in_cm_per_s)*(All.UnitLength_in_cm/CM_PER_KPC); //from km s^-1 kpc^-1 to unit code
	
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
	{
		//Static potential//
		
		dp[0]=fabs(P[i].Pos[0]);
		dp[1]=fabs(P[i].Pos[0]);
		dp[2]=fabs(P[i].Pos[0]);
	
   		x = (int) ((dp[0] - All.xx0) / All.deltax);
    	x = DMAX(x, 0);
    
    	y = (int) ((dp[1] - All.yy0) / All.deltay);
    	y = DMAX(y, 0);
    
		z = (int) ((dp[2] - All.zz0) / All.deltaz);
    	z = DMAX(z, 0);
		
		/*interp1_x = All.accx[x][y][z] + (All.accx[x+1][y][z]-All.accx[x][y][z])*(dp[0] - All.xx0 - x * All.deltax)/All.deltax;
		interp2_x = All.accx[x][y+1][z] + (All.accx[x+1][y+1][z]-All.accx[x][y+1][z])*(dp[0] - All.xx0 - x * All.deltax)/All.deltax;
		interp3_x = All.accx[x][y][z+1] + (All.accx[x+1][y][z+1]-All.accx[x][y][z+1])*(dp[0] - All.xx0 - x * All.deltax)/All.deltax;
		interp4_x = All.accx[x][y+1][z+1] + (All.accx[x+1][y+1][z+1]-All.accx[x][y+1][z+1])*(dp[0] - All.xx0 - x * All.deltax)/All.deltax;
		interp1_y = interp1_x + (interp2_x-interp1_x)*(dp[1] - All.yy0 - y * All.deltay)/All.deltay;
		interp2_y = interp3_x + (interp4_x-interp3_x)*(dp[1] - All.yy0 - y * All.deltay)/All.deltay;*/
				
		
	    interp1_x = All.accx[x][y][z] + (All.accx[x+1][y][z]-All.accx[x][y][z])/All.deltax * (dp[0] - All.xx0 - x * All.deltax) + 
			(All.accx[x][y+1][z] - All.accx[x][y][z] + (All.accx[x+1][y+1][z]-All.accx[x][y+1][z])/All.deltax * (dp[0] - All.xx0 - x * All.deltax) 
				- (All.accx[x+1][y][z]-All.accx[x][y][z])/All.deltax * (dp[0] - All.xx0 - x * All.deltax))/All.deltay * (dp[1] - All.yy0 - y * All.deltay);
	
		interp2_x = All.accx[x][y][z+1] + (All.accx[x+1][y][z+1]-All.accx[x][y][z+1])/All.deltax * (dp[0] - All.xx0 - x * All.deltax) + 
			(All.accx[x][y+1][z+1] - All.accx[x][y][z+1] + (All.accx[x+1][y+1][z+1]-All.accx[x][y+1][z+1])/All.deltax * (dp[0] - All.xx0 - x * All.deltax) 
				- (All.accx[x+1][y][z+1]-All.accx[x][y][z+1])/All.deltax * (dp[0] - All.xx0 - x * All.deltax))/All.deltay * (dp[1] - All.yy0 - y * All.deltay);
		
	    interp1_y = All.accy[x][y][z] + (All.accy[x+1][y][z]-All.accy[x][y][z])/All.deltax * (dp[0] - All.xx0 - x * All.deltax) + 
			(All.accy[x][y+1][z] - All.accy[x][y][z] + (All.accy[x+1][y+1][z]-All.accy[x][y+1][z])/All.deltax * (dp[0] - All.xx0 - x * All.deltax) 
				- (All.accy[x+1][y][z]-All.accy[x][y][z])/All.deltax * (dp[0] - All.xx0 - x * All.deltax))/All.deltay * (dp[1] - All.yy0 - y * All.deltay);
	
		interp2_y = All.accy[x][y][z+1] + (All.accy[x+1][y][z+1]-All.accy[x][y][z+1])/All.deltax * (dp[0] - All.xx0 - x * All.deltax) + 
			(All.accy[x][y+1][z+1] - All.accy[x][y][z+1] + (All.accy[x+1][y+1][z+1]-All.accy[x][y+1][z+1])/All.deltax * (dp[0] - All.xx0 - x * All.deltax) 
				- (All.accy[x+1][y][z+1]-All.accy[x][y][z+1])/All.deltax * (dp[0] - All.xx0 - x * All.deltax))/All.deltay * (dp[1] - All.yy0 - y * All.deltay);
		
	    interp1_z = All.accz[x][y][z] + (All.accz[x+1][y][z]-All.accz[x][y][z])/All.deltax * (dp[0] - All.xx0 - x * All.deltax) + 
			(All.accz[x][y+1][z] - All.accz[x][y][z] + (All.accz[x+1][y+1][z]-All.accz[x][y+1][z])/All.deltax * (dp[0] - All.xx0 - x * All.deltax) 
				- (All.accz[x+1][y][z]-All.accz[x][y][z])/All.deltax * (dp[0] - All.xx0 - x * All.deltax))/All.deltay * (dp[1] - All.yy0 - y * All.deltay);
	
		interp2_z = All.accz[x][y][z+1] + (All.accz[x+1][y][z+1]-All.accz[x][y][z+1])/All.deltax * (dp[0] - All.xx0 - x * All.deltax) + 
			(All.accz[x][y+1][z+1] - All.accx[x][y][z+1] + (All.accz[x+1][y+1][z+1]-All.accz[x][y+1][z+1])/All.deltax * (dp[0] - All.xx0 - x * All.deltax) 
				- (All.accz[x+1][y][z+1]-All.accz[x][y][z+1])/All.deltax * (dp[0] - All.xx0 - x * All.deltax))/All.deltay * (dp[1] - All.yy0 - y * All.deltay);
		
		P[i].GravAccel[0] = interp1_x + (interp2_x - interp1_x)/All.deltaz * (dp[2] - All.zz0 - z * All.deltaz);
		P[i].GravAccel[1] = interp1_y + (interp2_y - interp1_y)/All.deltaz * (dp[2] - All.zz0 - z * All.deltaz);
		P[i].GravAccel[2] = interp1_z + (interp2_z - interp1_z)/All.deltaz * (dp[2] - All.zz0 - z * All.deltaz);
		
		if((P[i].Pos[0] - All.xx0) < 0) P[i].GravAccel[0] = - P[i].GravAccel[0];
		if((P[i].Pos[1] - All.yy0) < 0) P[i].GravAccel[1] = - P[i].GravAccel[1];
		if((P[i].Pos[2] - All.zz0) < 0) P[i].GravAccel[2] = - P[i].GravAccel[2];

#ifdef TWODIMS
		P[i].GravAccel[2] = 0;
#endif		
		
		//if(fabs(P[i].GravAccel[0]) > 1e-4) printf("Accel: %e %e %e %e %e %e \n", P[i].Pos[0], P[i].Pos[1],P[i].Pos[2],P[i].GravAccel[0], P[i].GravAccel[1], P[i].GravAccel[2]);
		//if(fabs(P[i].GravAccel[0]) > 1e-4) printf("Accel: %e %d %d %d %e %e %e %e \n", P[i].Pos[0], x, y, z,All.accx[x][y][z],All.accx[x+1][y][z],All.accx[x][y+1][z],All.accx[x][y][z+1]);
		//if((dp[2] - All.zz0) < 0) P[i].GravAccel[2] = - P[i].GravAccel[2];
		//printf("Sono qui");
		//Bar potential//
		
		omega = vrot*t; //rotation angle
		omega *= PI_VAL/180;//radiant	
		dp[0]=fabs(P[i].Pos[0]*cos(omega)+P[i].Pos[1]*sin(omega)); 
		dp[1]=fabs(P[i].Pos[1]*cos(omega)-P[i].Pos[0]*sin(omega)); 
		
   		x = (int) ((dp[0] - All.xx0) / All.deltax);
    	x = DMAX(x, 0);
    
    	y = (int) ((dp[1] - All.yy0) / All.deltay);
    	y = DMAX(y, 0);
		
	    interp1_x = All.accx_bar[x][y][z] + (All.accx_bar[x+1][y][z]-All.accx_bar[x][y][z])/All.deltax * (dp[0] - All.xx0 - x * All.deltax) + 
			(All.accx_bar[x][y+1][z] - All.accx_bar[x][y][z] + (All.accx_bar[x+1][y+1][z]-All.accx_bar[x][y+1][z])/All.deltax * (dp[0] - All.xx0 - x * All.deltax) 
				- (All.accx_bar[x+1][y][z]-All.accx_bar[x][y][z])/All.deltax * (dp[0] - All.xx0 - x * All.deltax))/All.deltay * (dp[1] - All.yy0 - y * All.deltay);
	
		interp2_x = All.accx_bar[x][y][z+1] + (All.accx_bar[x+1][y][z+1]-All.accx_bar[x][y][z+1])/All.deltax * (dp[0] - All.xx0 - x * All.deltax) + 
			(All.accx_bar[x][y+1][z+1] - All.accx_bar[x][y][z+1] + (All.accx_bar[x+1][y+1][z+1]-All.accx_bar[x][y+1][z+1])/All.deltax * (dp[0] - All.xx0 - x * All.deltax) 
				- (All.accx_bar[x+1][y][z+1]-All.accx_bar[x][y][z+1])/All.deltax * (dp[0] - All.xx0 - x * All.deltax))/All.deltay * (dp[1] - All.yy0 - y * All.deltay);
		
	    interp1_y = All.accy_bar[x][y][z] + (All.accy_bar[x+1][y][z]-All.accy_bar[x][y][z])/All.deltax * (dp[0] - All.xx0 - x * All.deltax) + 
			(All.accy_bar[x][y+1][z] - All.accy_bar[x][y][z] + (All.accy_bar[x+1][y+1][z]-All.accy_bar[x][y+1][z])/All.deltax * (dp[0] - All.xx0 - x * All.deltax) 
				- (All.accy_bar[x+1][y][z]-All.accy_bar[x][y][z])/All.deltax * (dp[0] - All.xx0 - x * All.deltax))/All.deltay * (dp[1] - All.yy0 - y * All.deltay);
	
		interp2_y = All.accy_bar[x][y][z+1] + (All.accy_bar[x+1][y][z+1]-All.accy_bar[x][y][z+1])/All.deltax * (dp[0] - All.xx0 - x * All.deltax) + 
			(All.accy_bar[x][y+1][z+1] - All.accy_bar[x][y][z+1] + (All.accy_bar[x+1][y+1][z+1]-All.accy_bar[x][y+1][z+1])/All.deltax * (dp[0] - All.xx0 - x * All.deltax) 
				- (All.accy_bar[x+1][y][z+1]-All.accy_bar[x][y][z+1])/All.deltax * (dp[0] - All.xx0 - x * All.deltax))/All.deltay * (dp[1] - All.yy0 - y * All.deltay);
		
	    interp1_z = All.accz_bar[x][y][z] + (All.accz_bar[x+1][y][z]-All.accz_bar[x][y][z])/All.deltax * (dp[0] - All.xx0 - x * All.deltax) + 
			(All.accz_bar[x][y+1][z] - All.accz_bar[x][y][z] + (All.accz_bar[x+1][y+1][z]-All.accz_bar[x][y+1][z])/All.deltax * (dp[0] - All.xx0 - x * All.deltax) 
				- (All.accz_bar[x+1][y][z]-All.accz_bar[x][y][z])/All.deltax * (dp[0] - All.xx0 - x * All.deltax))/All.deltay * (dp[1] - All.yy0 - y * All.deltay);
	
		interp2_z = All.accz_bar[x][y][z+1] + (All.accz_bar[x+1][y][z+1]-All.accz_bar[x][y][z+1])/All.deltax * (dp[0] - All.xx0 - x * All.deltax) + 
			(All.accz_bar[x][y+1][z+1] - All.accx_bar[x][y][z+1] + (All.accz_bar[x+1][y+1][z+1]-All.accz_bar[x][y+1][z+1])/All.deltax * (dp[0] - All.xx0 - x * All.deltax) 
				- (All.accz_bar[x+1][y][z+1]-All.accz_bar[x][y][z+1])/All.deltax * (dp[0] - All.xx0 - x * All.deltax))/All.deltay * (dp[1] - All.yy0 - y * All.deltay);
		
		abar_x = interp1_x + (interp2_x - interp1_x)/All.deltaz * (dp[2] - All.zz0 - z * All.deltaz);
		abar_y = interp1_y + (interp2_y - interp1_y)/All.deltaz * (dp[2] - All.zz0 - z * All.deltaz);
		
		if((P[i].Pos[0]*cos(omega)+P[i].Pos[1]*sin(omega) - All.xx0) < 0) abar_x = - abar_x;
		if((P[i].Pos[1]*cos(omega)-P[i].Pos[0]*sin(omega) - All.yy0) < 0) abar_y = - abar_y; 
		
		P[i].GravAccel[0] += abar_x*cos(omega) - abar_y*sin(omega);
		P[i].GravAccel[1] += abar_x*sin(omega) - abar_y*cos(omega);
		 
		
		/*abar_tot = sqrt(abar_x*abar_x+abar_y*abar_y);
		theta = fabs(atan(abar_y/abar_x));		
		
		if((dp[0] - All.xx0) < 0) P[i].GravAccel[0] -= abar_tot * cos(omega+theta); 
		if((dp[1] - All.yy0) < 0) P[i].GravAccel[1] -= abar_tot * sin(omega+theta);
		
		if((dp[0] - All.xx0) > 0) P[i].GravAccel[0] += abar_tot * cos(omega+theta); 
		if((dp[1] - All.xx0) > 0) P[i].GravAccel[1] += abar_tot * sin(omega+theta);*/
		
		if((P[i].Pos[2] - All.xx0) > 0) P[i].GravAccel[2] += interp1_z + (interp2_z - interp1_z)/All.deltaz * (dp[2] - All.zz0 - z * All.deltaz);
		if((P[i].Pos[2] - All.zz0) < 0) P[i].GravAccel[2] -= interp1_z + (interp2_z - interp1_z)/All.deltaz * (dp[2] - All.zz0 - z * All.deltaz);
#ifdef TWODIMS
		P[i].GravAccel[2] = 0;
#endif		

	}
	
}
#endif

/* static NFW potential */
void GravAccel_StaticNFW()
{
    double NFW_C=12;
    double NFW_M200=100.0;
    double NFW_Eps=0.01;
    double NFW_DARKFRACTION=0.87;
    double NFW_BOXCENTERED;
    NFW_BOXCENTERED=1;

    /* convert units */
    double R200 = pow(NFW_M200 * All.G / (100 * All.Hubble_H0_CodeUnits * All.Hubble_H0_CodeUnits), 1.0 / 3);
    double Rs = R200 / NFW_C;
    double Dc = 200.0 / 3 * NFW_C * NFW_C * NFW_C / (log(1 + NFW_C) - NFW_C / (1 + NFW_C));
    double RhoCrit = 3 * All.Hubble_H0_CodeUnits * All.Hubble_H0_CodeUnits / (8 * M_PI * All.G);
    double V200 = 10 * All.Hubble_H0_CodeUnits * R200;
    
    double r0, r, R, m, dp[3], r2, fac; int i;
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        dp[0]=P[i].Pos[0]; dp[1]=P[i].Pos[1]; dp[2]=P[i].Pos[2];
#ifdef ANALYTIC_GRAVITY_ANCHOR_TO_PARTICLE
        int k; for(k = 0; k < 3; k++) {dp[k] = -P[i].min_xyz_to_bh[k];}
#elif defined(PERIODIC)
        if(NFW_BOXCENTERED) {dp[0] -= boxHalf_X; dp[1] -= boxHalf_Y; dp[2] -= boxHalf_Z;}
#endif
        r2 = dp[0]*dp[0] + dp[1]*dp[1] + dp[2]*dp[2]; r0 = sqrt(r2);

        /* function to get enclosed mass(<r) for NFW: */
        /* Eps is in units of Rs !!!! :: use unsoftened NFW if NFW_Eps=0 */
        R = r0;
        if(NFW_Eps > 0.0)
            if(R > Rs * NFW_C)
                R = Rs * NFW_C;
        
        fac=1.0;
        if(NFW_Eps > 0.0)
        {
            m = fac * 4 * M_PI * RhoCrit * Dc * (-(Rs * Rs * Rs * (1 - NFW_Eps + log(Rs) - 2 * NFW_Eps * log(Rs) +
              NFW_Eps * NFW_Eps * log(NFW_Eps * Rs))) / ((NFW_Eps - 1) * (NFW_Eps - 1)) + (Rs * Rs * Rs * (Rs -
              NFW_Eps * Rs - (2 * NFW_Eps - 1) * (R + Rs) * log(R + Rs) + NFW_Eps * NFW_Eps * (R + Rs) * log(R + NFW_Eps * Rs)))
              / ((NFW_Eps - 1) * (NFW_Eps - 1) * (R + Rs)));
        }
        else /* analytic NFW */
        {
            m = fac * 4 * M_PI * RhoCrit * Dc *
            (-(Rs * Rs * Rs * (1 + log(Rs))) + Rs * Rs * Rs * (Rs + (R + Rs) * log(R + Rs)) / (R + Rs));
        }
        fac = V200 * V200 * V200 / (10 * All.G * All.Hubble_H0_CodeUnits) / m;
        if(NFW_Eps > 0.0)
        {
            m = fac * 4 * M_PI * RhoCrit * Dc * (-(Rs * Rs * Rs * (1 - NFW_Eps + log(Rs) - 2 * NFW_Eps * log(Rs) +
              NFW_Eps * NFW_Eps * log(NFW_Eps * Rs))) / ((NFW_Eps - 1) * (NFW_Eps - 1)) + (Rs * Rs * Rs * (Rs -
              NFW_Eps * Rs - (2 * NFW_Eps - 1) * (R + Rs) * log(R + Rs) + NFW_Eps * NFW_Eps * (R + Rs) * log(R + NFW_Eps * Rs)))
              / ((NFW_Eps - 1) * (NFW_Eps - 1) * (R + Rs)));
        }
        else /* analytic NFW */
        {
            m = fac * 4 * M_PI * RhoCrit * Dc *
            (-(Rs * Rs * Rs * (1 + log(Rs))) + Rs * Rs * Rs * (Rs + (R + Rs) * log(R + Rs)) / (R + Rs));
        }
        m *= NFW_DARKFRACTION; r=r0;

        if(r > 0)
        {
            P[i].GravAccel[0] += -All.G * m * dp[0] / (r * r * r);
            P[i].GravAccel[1] += -All.G * m * dp[1] / (r * r * r);
            P[i].GravAccel[2] += -All.G * m * dp[2] / (r * r * r);
            
        } // if(r > 0) //
    } // for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i]) //
}




/* Paczysnky Wiita pseudo-Newtonian potential, G = M_sol = c = 1 */
void GravAccel_PaczynskyWiita()
{
    double PACZYNSKY_WIITA_MASS = 1.0; // Mass to use for the Paczynksy-Wiita analytic gravity pseudo-Newtonian potential (in solar masses)
    double r_g = 2*PACZYNSKY_WIITA_MASS;
    int i, k;
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        double dp[3], r2, r;
        dp[0]=P[i].Pos[0]; dp[1]=P[i].Pos[1]; dp[2]=P[i].Pos[2];
#ifdef ANALYTIC_GRAVITY_ANCHOR_TO_PARTICLE
        for(k = 0; k < 3; k++) {dp[k] = -P[i].min_xyz_to_bh[k];}
#endif
        r2 = dp[0]*dp[0] + dp[1]*dp[1] + dp[2]*dp[2]; r = sqrt(r2);
        if(r > r_g)
        {
            double q = PACZYNSKY_WIITA_MASS/((r - r_g)*(r - r_g));
            for(k = 0; k < 3; k++) {P[i].GravAccel[k] = - q * P[i].Pos[k]/r;}
        }
    }
}

#ifdef PARTICLE_EXCISION
void apply_excision(void)
{
    double EXCISION_MASS = 0; // mass of the excised object. Used to move the excision boundary so as to capture bound objects. If zero the excision boundary will not move
    double EXCISION_INIT_RADIUS = 0; // initial excision radius
    double EXCISION_ETA = 1; // remove particles with radius < EXCISION_ETA R_excision
    double excision_radius = EXCISION_ETA * pow(EXCISION_INIT_RADIUS*EXCISION_INIT_RADIUS*EXCISION_INIT_RADIUS +
                                                3.*sqrt(2. * All.G * EXCISION_MASS) * pow(EXCISION_INIT_RADIUS, 3./2.) * All.Time +
                                                9./2. * All.G * EXCISION_MASS * All.Time*All.Time, 1./3.);
    int i;
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(P[i].Type == 0)
        {
            double dp[3], r2, r;
            dp[0]=P[i].Pos[0]; dp[1]=P[i].Pos[1]; dp[2]=P[i].Pos[2];
#ifdef ANALYTIC_GRAVITY_ANCHOR_TO_PARTICLE
            int k; for(k = 0; k < 3; k++) {dp[k] = -P[i].min_xyz_to_bh[k];}
#endif
            r2 = dp[0]*dp[0] + dp[1]*dp[1] + dp[2]*dp[2]; r = sqrt(r2);
            if(r < excision_radius) P[i].Mass = 0;
        }
    }
}
#endif

