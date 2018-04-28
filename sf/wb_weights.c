#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"
#ifdef PTHREADS_NUM_THREADS
#include <pthread.h>
#endif
#ifdef PTHREADS_NUM_THREADS
extern pthread_mutex_t mutex_nexport;
extern pthread_mutex_t mutex_partnodedrift;
#define LOCK_NEXPORT     pthread_mutex_lock(&mutex_nexport);
#define UNLOCK_NEXPORT   pthread_mutex_unlock(&mutex_nexport);
#else
#define LOCK_NEXPORT
#define UNLOCK_NEXPORT
#endif

#ifdef SN_FEEDBACK

struct wb_data_in
{
    MyDouble Pos[3];
    MyDouble Hsml;
	int Nsn;
    MyDouble Mej;
	MyDouble Mass;
	MyDouble omegab_p[3];
	MyDouble omegab_m[3];
	MyFloat NV_T[3][3];
	MyDouble Dens;

    int NodeList[NODELISTLENGTH];
}
*wbDataIn, *wbDataGet;
           
/* define properties to be injected */
void particle2in_wb(struct wb_data_in *in, int i);
void particle2in_wb(struct wb_data_in *in, int i)
{
    if((P[i].Nsn_timestep<=0)||(P[i].DensAroundStar<=0)) {in->Mej=0; return;} // trap for no sne
    int k,j;
    for(k=0; k<3; k++) {in->Pos[k] = P[i].Pos[k];}
    in->Hsml = PPP[i].Hsml;
	in->Nsn = P[i].Nsn_timestep;
	in->Mass = P[i].Mass;
	in->Dens = P[i].DensAroundStar;
	for(k=0; k<3; k++) {in->omegab_p[k] = P[i].omegab_p[k];}
	for(k=0; k<3; k++) {in->omegab_m[k] = P[i].omegab_m[k];}
    for(j=0;j<3;j++)
    {
        for(k=0;k<3;k++)
        {
            in->NV_T[j][k] = P[i].NVT[j][k];
        }
    }
	
}



struct wb_data_out
{
	MyDouble weightb_tot;
}
*wbDataResult, *wbDataOut;

void out2particle_wb(struct wb_data_out *out, int i, int mode);
void out2particle_wb(struct wb_data_out *out, int i, int mode)
{
	if (mode == 0) {P[i].wb_tot = (out->weightb_tot);}
	if (mode == 1) {P[i].wb_tot += (out->weightb_tot);}
}


struct kernel_FB
{
    double	dp[3],r,wk_i, wk_j, dwk_i, dwk_j,h_i;
};


int wb_evaluate_active_check(int i)
{
    if(P[i].Type != 4) return 0;
    if(P[i].Mass <= 0) return 0;
    if(PPP[i].Hsml <= 0) return 0;
    if(PPP[i].NumNgb <= 0) return 0;
    if(P[i].Nsn_timestep>0) {return 1;}
    return 0;
}


void *wb_evaluate_primary(void *p)
{
#define CONDITION_FOR_EVALUATION if(wb_evaluate_active_check(i)==1)
#define EVALUATION_CALL wb_evaluate(i, 0, exportflag, exportnodecount, exportindex, ngblist)
#include "../system/code_block_primary_loop_evaluation.h"
#undef CONDITION_FOR_EVALUATION
#undef EVALUATION_CALL
}
void *wb_evaluate_secondary(void *p)
{
#define EVALUATION_CALL wb_evaluate(j, 1, &dummy, &dummy, &dummy, ngblist);
#include "../system/code_block_secondary_loop_evaluation.h"
#undef EVALUATION_CALL
}


void wb_calc(void)
{
    int j, k, ngrp, ndone, ndone_flag;
    int recvTask, place;
    int save_NextParticle;
    long long n_exported = 0;
    
    /* allocate buffers to arrange communication */
    long long NTaskTimesNumPart;
    NTaskTimesNumPart = maxThreads * NumPart;
    Ngblist = (int *) mymalloc("Ngblist", NTaskTimesNumPart * sizeof(int));
    size_t MyBufferSize = All.BufferSize;
    All.BunchSize = (int) ((MyBufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
                                             sizeof(struct wb_data_in) + sizeof(struct wb_data_out) + sizemax(sizeof(struct wb_data_in),sizeof(struct wb_data_out))));
    DataIndexTable = (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
    DataNodeList = (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));
    
    NextParticle = FirstActiveParticle;	/* begin with this index */
    do
    {
        
        BufferFullFlag = 0;
        Nexport = 0;
        save_NextParticle = NextParticle;
        
        for(j = 0; j < NTask; j++)
        {
            Send_count[j] = 0;
            Exportflag[j] = -1;
        }
        
        /* do local particles and prepare export list */
#ifdef PTHREADS_NUM_THREADS
        pthread_t mythreads[PTHREADS_NUM_THREADS - 1];
        int threadid[PTHREADS_NUM_THREADS - 1];
        pthread_attr_t attr;
        
        pthread_attr_init(&attr);
        pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
        pthread_mutex_init(&mutex_nexport, NULL);
        pthread_mutex_init(&mutex_partnodedrift, NULL);
        
        TimerFlag = 0;
        
        for(j = 0; j < PTHREADS_NUM_THREADS - 1; j++)
        {
            threadid[j] = j + 1;
            pthread_create(&mythreads[j], &attr, wb_evaluate_primary, &threadid[j]);
        }
#endif
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
#ifdef _OPENMP
            int mainthreadid = omp_get_thread_num();
#else
            int mainthreadid = 0;
#endif
            wb_evaluate_primary(&mainthreadid);	/* do local particles and prepare export list */
        }
        
#ifdef PTHREADS_NUM_THREADS
        for(j = 0; j < PTHREADS_NUM_THREADS - 1; j++)
            pthread_join(mythreads[j], NULL);
#endif
        
        
        if(BufferFullFlag)
        {
            int last_nextparticle = NextParticle;
            
            NextParticle = save_NextParticle;
            
            while(NextParticle >= 0)
            {
                if(NextParticle == last_nextparticle)
                    break;
                
                if(ProcessedFlag[NextParticle] != 1)
                    break;
                
                ProcessedFlag[NextParticle] = 2;
                
                NextParticle = NextActiveParticle[NextParticle];
            }
            
            if(NextParticle == save_NextParticle)
            {
                /* in this case, the buffer is too small to process even a single particle */
                endrun(116608);
            }
            
            int new_export = 0;
            
            for(j = 0, k = 0; j < Nexport; j++)
                if(ProcessedFlag[DataIndexTable[j].Index] != 2)
                {
                    if(k < j + 1)
                        k = j + 1;
                    
                    for(; k < Nexport; k++)
                        if(ProcessedFlag[DataIndexTable[k].Index] == 2)
                        {
                            int old_index = DataIndexTable[j].Index;
                            
                            DataIndexTable[j] = DataIndexTable[k];
                            DataNodeList[j] = DataNodeList[k];
                            DataIndexTable[j].IndexGet = j;
                            new_export++;
                            
                            DataIndexTable[k].Index = old_index;
                            k++;
                            break;
                        }
                }
                else
                    new_export++;
            
            Nexport = new_export;
            
        }
        
        n_exported += Nexport;
        
        for(j = 0; j < NTask; j++)
            Send_count[j] = 0;
        for(j = 0; j < Nexport; j++)
            Send_count[DataIndexTable[j].Task]++;
        
        MYSORT_DATAINDEX(DataIndexTable, Nexport, sizeof(struct data_index), data_index_compare);
        MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);
        
        for(j = 0, Nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
        {
            Nimport += Recv_count[j];
            
            if(j > 0)
            {
                Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
                Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
            }
        }
        
        wbDataGet = (struct wb_data_in *) mymalloc("wbDataGet", Nimport * sizeof(struct wb_data_in));
        wbDataIn = (struct wb_data_in *) mymalloc("wbDataIn", Nexport * sizeof(struct wb_data_in));
        
        /* prepare particle data for export */
        
        for(j = 0; j < Nexport; j++)
        {
            place = DataIndexTable[j].Index;
            particle2in_wb(&wbDataIn[j], place);
            memcpy(wbDataIn[j].NodeList,
                   DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
        }
        
        /* exchange particle data */
        int TAG_TO_USE = TAG_FBLOOP_1A;
        for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
        {
            recvTask = ThisTask ^ ngrp;
            
            if(recvTask < NTask)
            {
                if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                {
                    /* get the particles */
                    MPI_Sendrecv(&wbDataIn[Send_offset[recvTask]],
                                 Send_count[recvTask] * sizeof(struct wb_data_in), MPI_BYTE,
                                 recvTask, TAG_TO_USE,
                                 &wbDataGet[Recv_offset[recvTask]],
                                 Recv_count[recvTask] * sizeof(struct wb_data_in), MPI_BYTE,
                                 recvTask, TAG_TO_USE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }
        
        myfree(wbDataIn);
        wbDataResult = (struct wb_data_out *) mymalloc("wbDataResult", Nimport * sizeof(struct wb_data_out));
        wbDataOut = (struct wb_data_out *) mymalloc("wbDataOut", Nexport * sizeof(struct wb_data_out));
        
        /* now do the particles that were sent to us */
        NextJ = 0;
#ifdef PTHREADS_NUM_THREADS
        for(j = 0; j < PTHREADS_NUM_THREADS - 1; j++)
            pthread_create(&mythreads[j], &attr, wb_evaluate_secondary, &threadid[j]);
#endif
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
#ifdef _OPENMP
            int mainthreadid = omp_get_thread_num();
#else
            int mainthreadid = 0;
#endif
            wb_evaluate_secondary(&mainthreadid);
        }
        
#ifdef PTHREADS_NUM_THREADS
        for(j = 0; j < PTHREADS_NUM_THREADS - 1; j++)
            pthread_join(mythreads[j], NULL);
        
        pthread_mutex_destroy(&mutex_partnodedrift);
        pthread_mutex_destroy(&mutex_nexport);
        pthread_attr_destroy(&attr);
#endif
        
        if(NextParticle < 0)
            ndone_flag = 1;
        else
            ndone_flag = 0;
        
        MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        
        /* get the result */
        TAG_TO_USE = TAG_FBLOOP_1B;
        for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
        {
            recvTask = ThisTask ^ ngrp;
            if(recvTask < NTask)
            {
                if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                {
                    /* send the results */
                    MPI_Sendrecv(&wbDataResult[Recv_offset[recvTask]],
                                 Recv_count[recvTask] * sizeof(struct wb_data_out),
                                 MPI_BYTE, recvTask, TAG_TO_USE,
                                 &wbDataOut[Send_offset[recvTask]],
                                 Send_count[recvTask] * sizeof(struct wb_data_out),
                                 MPI_BYTE, recvTask, TAG_TO_USE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }
		
		
        
        /* add the result to the local particles */
        for(j = 0; j < Nexport; j++)
        {
            place = DataIndexTable[j].Index;
            out2particle_wb(&wbDataOut[j], place, 1);
        }
        myfree(wbDataOut);
        myfree(wbDataResult);
        myfree(wbDataGet);
    }
    while(ndone < NTask);
    
    myfree(DataNodeList);
    myfree(DataIndexTable);
    myfree(Ngblist);
}


int wb_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist)
{
    int startnode, numngb_inbox, listindex = 0, j, k, n;
	double hinv, hinv3, hinv4, r2, u, hinv_j, hinv3_j, hinv4_j;
    struct kernel_FB kernel;
    struct wb_data_in local;
    struct wb_data_out out;
    memset(&out, 0, sizeof(struct wb_data_out));
    //kernel_main(0.0,1.0,1.0,&kernel_zero,&wk,-1);
	out.weightb_tot = 0;
	
    
	/* Load the data for the particle injecting feedback */
    if(mode == 0) {particle2in_wb(&local, target);} else {local = wbDataGet[target];}
    if(local.Nsn<=0) return 0; // no SNe for the master particle! nothing to do here //
    if(local.Hsml<=0) return 0; // zero-extent kernel, no particles //
    //h2 = local.Hsml*local.Hsml;
    //kernel_hinv(local.Hsml, &kernel.hinv, &kernel.hinv3, &kernel.hinv4);
    kernel.h_i = local.Hsml;
    double h2_i = kernel.h_i*kernel.h_i;
	double V_i = local.Mass/local.Dens;
		 
    /* Now start the weights computation for this particle */
    if(mode == 0)
    {
        startnode = All.MaxPart;	/* root node */
    }
    else
    {
        startnode = wbDataGet[target].NodeList[0];
        startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }
    while(startnode >= 0)
    {
        while(startnode >= 0)
        { 
			//find neighbors that can -mutually- see one another, not just single-directional searching here
            numngb_inbox = ngb_treefind_pairs_threads(local.Pos, local.Hsml, target, &startnode, mode, exportflag, exportnodecount, exportindex, ngblist);
            if(numngb_inbox < 0) return -1;
            for(n = 0; n < numngb_inbox; n++)
            {
                j = ngblist[n];
                if(P[j].Type != 0) continue; // require a gas particle //
                if(P[j].Mass <= 0) continue; // require the particle has mass //
                for(k=0; k<3; k++) {kernel.dp[k] = local.Pos[k] - P[j].Pos[k];}
#ifdef PERIODIC
                NEAREST_XYZ(kernel.dp[0],kernel.dp[1],kernel.dp[2],1); // find the closest image in the given box size  //
#endif
                r2=0; for(k=0;k<3;k++) {r2 += kernel.dp[k]*kernel.dp[k];}
                double h_j = PPP[j].Hsml;
                if(r2 <= 0) continue; // same particle //
                if((r2 >= h2_i) && (r2 >= h_j * h_j)) continue; // outside kernel //
				if(sqrt(r2) > 2 * CM_PER_KPC/All.UnitLength_in_cm) continue; //search radius not larger than 2 kpc
                kernel.r = sqrt(r2);

                if(kernel.r < kernel.h_i)
                {
				    kernel_hinv(kernel.h_i, &hinv, &hinv3, &hinv4);
                    u = kernel.r * hinv;
                    kernel_main(u, hinv3, hinv4, &kernel.wk_i, &kernel.dwk_i, 0);
                }
                else
                {
                    kernel.dwk_i = kernel.wk_i = 0;
                }
				
                if(kernel.r < h_j)
                {
                    kernel_hinv(h_j, &hinv_j, &hinv3_j, &hinv4_j);
                    u = kernel.r * hinv_j;
                    kernel_main(u, hinv3_j, hinv4_j, &kernel.wk_j, &kernel.dwk_j, 0);
                }
                else
                {
                    kernel.dwk_j = kernel.wk_j = 0;
                }
				
				double Xba[3], Xba_p[3], Xba_m[3], Xba_vers[3], fp[3], fm[3]; //gas particle position as seen from the star
				for(k=0; k<3; k++) {Xba[k] = P[j].Pos[k] - local.Pos[k];}
				for(k=0; k<3; k++) {Xba_vers[k] = Xba[k]/sqrt(r2);}
				for(k=0; k<3; k++) {Xba_p[k] = DMAX(0,Xba[k])/sqrt(r2);} //they are versors
				for(k=0; k<3; k++) {Xba_m[k] = DMIN(0,Xba[k])/sqrt(r2);} //they are versors
				
				double V_j = P[j].Mass / SphP[j].Density;
			    double wt_i,wt_j;
			    wt_i=V_i; wt_j=V_j;
				
#ifdef AGGRESSIVE_SLOPE_LIMITERS
			    wt_i=V_i; wt_j=V_j;
#else
#ifdef COOLING
			    if((fabs(V_i-V_j)/DMIN(V_i,V_j))/NUMDIMS > 1.25) {wt_i=wt_j=2.*V_i*V_j/(V_i+V_j);} else {wt_i=V_i; wt_j=V_j;}
#else
			    if((fabs(V_i-V_j)/DMIN(V_i,V_j))/NUMDIMS > 1.50) {wt_i=wt_j=(V_i*PPP[j].Hsml+V_j*local.Hsml)/(local.Hsml+PPP[j].Hsml);} else {wt_i=V_i; wt_j=V_j;}
#endif
#endif
				
				/* calculate the face area between the particles (must match what is done in the actual hydro routine!) */
				double Face_Area_Vec[3];
				double facenormal_dot_dp = 0, Face_Area_Norm = 0;
				for(k=0;k<3;k++)
                {
                    Face_Area_Vec[k] = kernel.wk_i * wt_i * (local.NV_T[k][0]*kernel.dp[0] + local.NV_T[k][1]*kernel.dp[1] + local.NV_T[k][2]*kernel.dp[2])
                                     + kernel.wk_j * wt_j * (SphP[j].NV_T[k][0]*kernel.dp[0] + SphP[j].NV_T[k][1]*kernel.dp[1] + SphP[j].NV_T[k][2]*kernel.dp[2]);
                    if(All.ComovingIntegrationOn) {Face_Area_Vec[k] *= All.cf_atime*All.cf_atime;} // Face_Area_Norm has units of area, need to convert to physical
					facenormal_dot_dp += Face_Area_Vec[k] * kernel.dp[k]; // check that face points same direction as vector normal: should be true for positive-definite (well-conditioned) NV_T 
					Face_Area_Norm += Face_Area_Vec[k] * Face_Area_Vec[k];
					//printf("AtimeX %e %e %e %e %e %e \n", kernel.wk_i, kernel.wk_j, wt_i, wt_j, V_i, V_j);
				}
			    if(facenormal_dot_dp < 0)
			    {
			        /* the effective gradient matrix is ill-conditioned (or not positive-definite!): for stability, we revert to the "RSPH" EOM */
					Face_Area_Norm = -(wt_i*V_i*kernel.dwk_i + wt_j*V_j*kernel.dwk_j) / kernel.r;
					Face_Area_Norm *= All.cf_atime*All.cf_atime; /* Face_Area_Norm has units of area, need to convert to physical */
			        Face_Area_Vec[0] = Face_Area_Norm * kernel.dp[0];
			        Face_Area_Vec[1] = Face_Area_Norm * kernel.dp[1];
			        Face_Area_Vec[2] = Face_Area_Norm * kernel.dp[2];
			        Face_Area_Norm = Face_Area_Norm * Face_Area_Norm * r2;
			    }
				
				double AtimeX = 0;				
				for(k=0; k<3; k++) {AtimeX += (-Face_Area_Vec[k])*Xba_vers[k];} //effective area has to be perpendicular to Xab
                
				SphP[j].omega_b = 0.5*(1.-1./sqrt(1.+((AtimeX)/(PI_VAL*r2))));
								
				double wbmod = 0.;
				for(k=0; k<3; k++)
				{
					//printf("omegab_m, omegab_p %e %e \n", local.omegab_m[k],local.omegab_p[k]);
					fp[k] = sqrt(0.5*(1+pow(local.omegab_m[k]/local.omegab_p[k],2)));
					fm[k] = sqrt(0.5*(1+pow(local.omegab_p[k]/local.omegab_m[k],2)));	
				    SphP[j].wb[k] = SphP[j].omega_b * (fp[k]*Xba_p[k] + fm[k]*Xba_m[k]);
					if (SphP[j].wb[k] != SphP[j].wb[k]) SphP[j].wb[k] = SphP[j].omega_b * Xba_vers[k];
					wbmod += SphP[j].wb[k]*SphP[j].wb[k];	
					
			    }
				
				out.weightb_tot += sqrt(wbmod);
				
            } // for(n = 0; n < numngb; n++)
        } // while(startnode >= 0)
        if(mode == 1)
        {
            listindex++;
            if(listindex < NODELISTLENGTH)
            {
                startnode = wbDataGet[target].NodeList[listindex];
                if(startnode >= 0) {startnode = Nodes[startnode].u.d.nextnode;}	/* open it */
            }
        } // if(mode == 1)
    } // while(startnode >= 0)
    /* Now collect the result at the right place */
    if(mode == 0) {out2particle_wb(&out, target, 0);} else {wbDataResult[target] = out;}
    return 0;
} // int wb_evaluate


#endif