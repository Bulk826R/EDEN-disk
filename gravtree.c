#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <mpi.h>
#include "halo_props.h"
#include "allvars.h"
#include "proto.h"

// Pos: comoving coordinates in units h^-1 Mpc
// physical pos: Pos*a/h = Pos*All.Time/All.HubbleParam
double AccRz(double R, double z, int type, double msz, double fstar) // comoving R,z in units Mpc/h
{
   // parameters [M]=1e10 Msun/h; [L]=Mpc/h
   // Stellar disk
   double  Rds0 = 2.5*0.001*All.HubbleParam;
   double  hds0 = 0.35*0.001*All.HubbleParam;
   double  Mds0 = 4.1*All.HubbleParam;
   // Gas disk
   double  Rdg0 = 7*0.001*All.HubbleParam;
   double  hdg0 = 0.084*0.001*All.HubbleParam;
   double  Mdg0 = 1.86*All.HubbleParam;
   // Bulge
   double  MB0  = 0.9*All.HubbleParam;
   double  rB0  = 0.5*0.001*All.HubbleParam;

   //----------------------
   double redshift = 1.0/All.Time -1.0;   
   double Msz = msz*All.HubbleParam/1e10; //UM Mstar for MW host
   //double Mdsz = Mds0/(MB0+Mds0+Mdg0)*Msz; 
   //double Mdgz = Mdg0/(MB0+Mds0+Mdg0)*Msz; 
   //double MBz  = MB0 /(MB0+Mds0+Mdg0)*Msz; 
   double Mdsz = Mds0/(MB0+Mds0)*Msz; 
   double MBz  = MB0 /(MB0+Mds0)*Msz; 
   double Mdgz = (Msz / fstar) * (1. - fstar);
   double Rdsz = Rds0*pow(1+redshift,0.28); //-0.72 -> 0.28 to take out scale factor, comoving
   double hdsz = hds0/Rds0*Rdsz;
   double Rdgz = Rdg0*pow(1+redshift,0.28); //-0.72 -> 0.28 to take out scale factor, comoving
   double hdgz = hdg0/Rdg0*Rdgz;
   double rBz  = rB0*pow(1+redshift,0.28); //-0.72 -> 0.28 to take out scale factor, comoving
   //----------------------
   double vals,valg,valb;
   if(type==0){ // minus acceleration along the z direction
      vals = Mdsz*z*pow(hdsz*hdsz + z*z, -0.5) * (Rdsz + pow(hdsz*hdsz + z*z, 0.5))*pow(R*R + pow(Rdsz + pow(hdsz*hdsz + z*z, 0.5), 2), -1.5);
      valg = Mdgz*z*pow(hdgz*hdgz + z*z, -0.5) * (Rdgz + pow(hdgz*hdgz + z*z, 0.5))*pow(R*R + pow(Rdgz + pow(hdgz*hdgz + z*z, 0.5), 2), -1.5);
      valb = MBz*z*pow(R*R + z*z, -0.5) * pow(rBz + pow(R*R + z*z, 0.5), -2);
   }
   else{ // minus acceleration along the R direction
      vals = Mdsz*R*pow(R*R + pow(Rdsz + pow(hdsz*hdsz + z*z, 0.5), 2), -1.5);
      valg = Mdgz*R*pow(R*R + pow(Rdgz + pow(hdgz*hdgz + z*z, 0.5), 2), -1.5);
      valb = MBz*R*pow(R*R + z*z,-0.5) * pow(rBz + pow(R*R + z*z, 0.5), -2);
   }
   double val = vals+valg+valb;
   return val; 

}

/*! \file gravtree.c 
 *  \brief main driver routines for gravitational (short-range) force computation
 *
 *  This file contains the code for the gravitational force computation by
 *  means of the tree algorithm. To this end, a tree force is computed for
 *  all active local particles, and particles are exported to other
 *  processors if needed, where they can receive additional force
 *  contributions. If the TreePM algorithm is enabled, the force computed
 *  will only be the short-range part.
 */

/*! This function computes the gravitational forces for all active
 *  particles.  If needed, a new tree is constructed, otherwise the
 *  dynamically updated tree is used.  Particles are only exported to other
 *  processors when really needed, thereby allowing a good use of the
 *  communication buffer.
 */
void gravity_tree(void)
{
  long long ntot;
  int numnodes, nexportsum = 0;
  int i, j, iter = 0;
  int *numnodeslist, maxnumnodes, nexport, *numlist, *nrecv, *ndonelist;
  double tstart, tend, timetree = 0, timecommsumm = 0, timeimbalance = 0, sumimbalance;
  double ewaldcount;
  double costtotal, ewaldtot, *costtreelist, *ewaldlist;
  double maxt, sumt, *timetreelist, *timecommlist;
  double fac, plb, plb_max, sumcomm;

#ifndef NOGRAVITY
  int *noffset, *nbuffer, *nsend, *nsend_local;
  long long ntotleft;
  int ndone, maxfill, ngrp;
  int k, place;
  int level, sendTask, recvTask;
  double ax, ay, az;
  MPI_Status status;
#endif

  /* set new softening lengths */
  if(All.ComovingIntegrationOn)
    set_softenings();


  /* contruct tree if needed */
  tstart = second();
  if(TreeReconstructFlag)
    {
      if(ThisTask == 0)
	printf("Tree construction.\n");

      force_treebuild(NumPart);

      TreeReconstructFlag = 0;

      if(ThisTask == 0)
	printf("Tree construction done.\n");
    }
  tend = second();
  All.CPU_TreeConstruction += timediff(tstart, tend);

  costtotal = ewaldcount = 0;

  /* Note: 'NumForceUpdate' has already been determined in find_next_sync_point_and_drift() */
  numlist = malloc(NTask * sizeof(int) * NTask);
  MPI_Allgather(&NumForceUpdate, 1, MPI_INT, numlist, 1, MPI_INT, MPI_COMM_WORLD);
  for(i = 0, ntot = 0; i < NTask; i++)
    ntot += numlist[i];
  free(numlist);


#ifndef NOGRAVITY
  if(ThisTask == 0)
    printf("Begin tree force.\n");


#ifdef SELECTIVE_NO_GRAVITY
  for(i = 0; i < NumPart; i++)
    if(((1 << P[i].Type) & (SELECTIVE_NO_GRAVITY)))
      P[i].Ti_endstep = -P[i].Ti_endstep - 1;
#endif


  noffset = malloc(sizeof(int) * NTask);	/* offsets of bunches in common list */
  nbuffer = malloc(sizeof(int) * NTask);
  nsend_local = malloc(sizeof(int) * NTask);
  nsend = malloc(sizeof(int) * NTask * NTask);
  ndonelist = malloc(sizeof(int) * NTask);

  i = 0;			/* beginn with this index */
  ntotleft = ntot;		/* particles left for all tasks together */

  while(ntotleft > 0)
    {
      iter++;

      for(j = 0; j < NTask; j++)
	nsend_local[j] = 0;

      /* do local particles and prepare export list */
      tstart = second();
      for(nexport = 0, ndone = 0; i < NumPart && nexport < All.BunchSizeForce - NTask; i++)
	if(P[i].Ti_endstep == All.Ti_Current)
	  {
	    ndone++;

	    for(j = 0; j < NTask; j++)
	      Exportflag[j] = 0;
#ifndef PMGRID
	    costtotal += force_treeevaluate(i, 0, &ewaldcount);
#else
	    costtotal += force_treeevaluate_shortrange(i, 0);
#endif
	    for(j = 0; j < NTask; j++)
	      {
		if(Exportflag[j])
		  {
		    for(k = 0; k < 3; k++)
		      GravDataGet[nexport].u.Pos[k] = P[i].Pos[k];
#ifdef UNEQUALSOFTENINGS
		    GravDataGet[nexport].Type = P[i].Type;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
		    if(P[i].Type == 0)
		      GravDataGet[nexport].Soft = SphP[i].Hsml;
#endif
#endif
		    GravDataGet[nexport].w.OldAcc = P[i].OldAcc;
		    GravDataIndexTable[nexport].Task = j;
		    GravDataIndexTable[nexport].Index = i;
		    GravDataIndexTable[nexport].SortIndex = nexport;
		    nexport++;
		    nexportsum++;
		    nsend_local[j]++;
		  }
	      }
	  }
      tend = second();
      timetree += timediff(tstart, tend);

      qsort(GravDataIndexTable, nexport, sizeof(struct gravdata_index), grav_tree_compare_key);

      for(j = 0; j < nexport; j++)
	GravDataIn[j] = GravDataGet[GravDataIndexTable[j].SortIndex];

      for(j = 1, noffset[0] = 0; j < NTask; j++)
	noffset[j] = noffset[j - 1] + nsend_local[j - 1];

      tstart = second();

      MPI_Allgather(nsend_local, NTask, MPI_INT, nsend, NTask, MPI_INT, MPI_COMM_WORLD);

      tend = second();
      timeimbalance += timediff(tstart, tend);

      /* now do the particles that need to be exported */

      for(level = 1; level < (1 << PTask); level++)
	{
	  tstart = second();
	  for(j = 0; j < NTask; j++)
	    nbuffer[j] = 0;
	  for(ngrp = level; ngrp < (1 << PTask); ngrp++)
	    {
	      maxfill = 0;
	      for(j = 0; j < NTask; j++)
		{
		  if((j ^ ngrp) < NTask)
		    if(maxfill < nbuffer[j] + nsend[(j ^ ngrp) * NTask + j])
		      maxfill = nbuffer[j] + nsend[(j ^ ngrp) * NTask + j];
		}
	      if(maxfill >= All.BunchSizeForce)
		break;

	      sendTask = ThisTask;
	      recvTask = ThisTask ^ ngrp;

	      if(recvTask < NTask)
		{
		  if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
		    {
		      /* get the particles */
		      MPI_Sendrecv(&GravDataIn[noffset[recvTask]],
				   nsend_local[recvTask] * sizeof(struct gravdata_in), MPI_BYTE,
				   recvTask, TAG_GRAV_A,
				   &GravDataGet[nbuffer[ThisTask]],
				   nsend[recvTask * NTask + ThisTask] * sizeof(struct gravdata_in), MPI_BYTE,
				   recvTask, TAG_GRAV_A, MPI_COMM_WORLD, &status);
		    }
		}

	      for(j = 0; j < NTask; j++)
		if((j ^ ngrp) < NTask)
		  nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
	    }
	  tend = second();
	  timecommsumm += timediff(tstart, tend);


	  tstart = second();
	  for(j = 0; j < nbuffer[ThisTask]; j++)
	    {
#ifndef PMGRID
	      costtotal += force_treeevaluate(j, 1, &ewaldcount);
#else
	      costtotal += force_treeevaluate_shortrange(j, 1);
              //if(j==0 && ThisTask==0) printf("Grav nbuffer[ThisTask=%d]=%d,level=%d\n",ThisTask,nbuffer[ThisTask],level);

#endif
	    }
	  tend = second();
	  timetree += timediff(tstart, tend);

	  tstart = second();
	  MPI_Barrier(MPI_COMM_WORLD);
	  tend = second();
	  timeimbalance += timediff(tstart, tend);

	  /* get the result */
	  tstart = second();
	  for(j = 0; j < NTask; j++)
	    nbuffer[j] = 0;
	  for(ngrp = level; ngrp < (1 << PTask); ngrp++)
	    {
	      maxfill = 0;
	      for(j = 0; j < NTask; j++)
		{
		  if((j ^ ngrp) < NTask)
		    if(maxfill < nbuffer[j] + nsend[(j ^ ngrp) * NTask + j])
		      maxfill = nbuffer[j] + nsend[(j ^ ngrp) * NTask + j];
		}
	      if(maxfill >= All.BunchSizeForce)
		break;

	      sendTask = ThisTask;
	      recvTask = ThisTask ^ ngrp;
	      if(recvTask < NTask)
		{
		  if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
		    {
		      /* send the results */
		      MPI_Sendrecv(&GravDataResult[nbuffer[ThisTask]],
				   nsend[recvTask * NTask + ThisTask] * sizeof(struct gravdata_in),
				   MPI_BYTE, recvTask, TAG_GRAV_B,
				   &GravDataOut[noffset[recvTask]],
				   nsend_local[recvTask] * sizeof(struct gravdata_in),
				   MPI_BYTE, recvTask, TAG_GRAV_B, MPI_COMM_WORLD, &status);

		      /* add the result to the particles */
		      for(j = 0; j < nsend_local[recvTask]; j++)
			{
			  place = GravDataIndexTable[noffset[recvTask] + j].Index;

			  for(k = 0; k < 3; k++)
			    P[place].GravAccel[k] += GravDataOut[j + noffset[recvTask]].u.Acc[k];

			  P[place].GravCost += GravDataOut[j + noffset[recvTask]].w.Ninteractions;
			}
		    }
		}

	      for(j = 0; j < NTask; j++)
		if((j ^ ngrp) < NTask)
		  nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
	    }
	  tend = second();
	  timecommsumm += timediff(tstart, tend);

	  level = ngrp - 1;
	}

      MPI_Allgather(&ndone, 1, MPI_INT, ndonelist, 1, MPI_INT, MPI_COMM_WORLD);
      for(j = 0; j < NTask; j++)
	ntotleft -= ndonelist[j];
    }

  free(ndonelist);
  free(nsend);
  free(nsend_local);
  free(nbuffer);
  free(noffset);

  /* now add things for comoving integration */

#ifndef PERIODIC
#ifndef PMGRID
  if(All.ComovingIntegrationOn)
    {
      fac = 0.5 * All.Hubble * All.Hubble * All.Omega0 / All.G;
      //printf("XXXX %g %g\n",All.Hubble,All.Omega0);

      for(i = 0; i < NumPart; i++)
	if(P[i].Ti_endstep == All.Ti_Current)
	  for(j = 0; j < 3; j++)
	    P[i].GravAccel[j] += fac * P[i].Pos[j];
    }
#endif
#endif

  for(i = 0; i < NumPart; i++)
    if(P[i].Ti_endstep == All.Ti_Current)
      {
#ifdef PMGRID
	ax = P[i].GravAccel[0] + P[i].GravPM[0] / All.G;
	ay = P[i].GravAccel[1] + P[i].GravPM[1] / All.G;
	az = P[i].GravAccel[2] + P[i].GravPM[2] / All.G;
#else
	ax = P[i].GravAccel[0];
	ay = P[i].GravAccel[1];
	az = P[i].GravAccel[2];
#endif
	P[i].OldAcc = sqrt(ax * ax + ay * ay + az * az);
      }


  if(All.TypeOfOpeningCriterion == 1)
    All.ErrTolTheta = 0;	/* This will switch to the relative opening criterion for the following force computations */

  /*  muliply by G */
  for(i = 0; i < NumPart; i++)
    if(P[i].Ti_endstep == All.Ti_Current)
      for(j = 0; j < 3; j++)
	P[i].GravAccel[j] *= All.G;


  /* Finally, the following factor allows a computation of a cosmological simulation 
     with vacuum energy in physical coordinates */
#ifndef PERIODIC
#ifndef PMGRID
  if(All.ComovingIntegrationOn == 0)
    {
      fac = All.OmegaLambda * All.Hubble * All.Hubble;

      for(i = 0; i < NumPart; i++)
	if(P[i].Ti_endstep == All.Ti_Current)
	  for(j = 0; j < 3; j++)
	    P[i].GravAccel[j] += fac * P[i].Pos[j];
    }
#endif
#endif

#ifdef SELECTIVE_NO_GRAVITY
  for(i = 0; i < NumPart; i++)
    if(P[i].Ti_endstep < 0)
      P[i].Ti_endstep = -P[i].Ti_endstep - 1;
#endif

  if(ThisTask == 0)
    printf("tree is done.\n");

#else /* gravity is switched off */

  for(i = 0; i < NumPart; i++)
    if(P[i].Ti_endstep == All.Ti_Current)
      for(j = 0; j < 3; j++)
	P[i].GravAccel[j] = 0;

#endif




  /* Now the force computation is finished */

  /*  gather some diagnostic information */

  timetreelist = malloc(sizeof(double) * NTask);
  timecommlist = malloc(sizeof(double) * NTask);
  costtreelist = malloc(sizeof(double) * NTask);
  numnodeslist = malloc(sizeof(int) * NTask);
  ewaldlist = malloc(sizeof(double) * NTask);
  nrecv = malloc(sizeof(int) * NTask);

  numnodes = Numnodestree;

  MPI_Gather(&costtotal, 1, MPI_DOUBLE, costtreelist, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(&numnodes, 1, MPI_INT, numnodeslist, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Gather(&timetree, 1, MPI_DOUBLE, timetreelist, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(&timecommsumm, 1, MPI_DOUBLE, timecommlist, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(&NumPart, 1, MPI_INT, nrecv, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Gather(&ewaldcount, 1, MPI_DOUBLE, ewaldlist, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Reduce(&nexportsum, &nexport, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&timeimbalance, &sumimbalance, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      All.TotNumOfForces += ntot;

      fprintf(FdTimings, "Step= %d  t= %g  dt= %g \n", All.NumCurrentTiStep, All.Time, All.TimeStep);
      fprintf(FdTimings, "Nf= %d%09d  total-Nf= %d%09d  ex-frac= %g  iter= %d\n",
	      (int) (ntot / 1000000000), (int) (ntot % 1000000000),
	      (int) (All.TotNumOfForces / 1000000000), (int) (All.TotNumOfForces % 1000000000),
	      nexport / ((double) ntot), iter);
      /* note: on Linux, the 8-byte integer could be printed with the format identifier "%qd", but doesn't work on AIX */

      fac = NTask / ((double) All.TotNumPart);

      for(i = 0, maxt = timetreelist[0], sumt = 0, plb_max = 0,
	  maxnumnodes = 0, costtotal = 0, sumcomm = 0, ewaldtot = 0; i < NTask; i++)
	{
	  costtotal += costtreelist[i];

	  sumcomm += timecommlist[i];

	  if(maxt < timetreelist[i])
	    maxt = timetreelist[i];
	  sumt += timetreelist[i];

	  plb = nrecv[i] * fac;

	  if(plb > plb_max)
	    plb_max = plb;

	  if(numnodeslist[i] > maxnumnodes)
	    maxnumnodes = numnodeslist[i];

	  ewaldtot += ewaldlist[i];
	}
      fprintf(FdTimings, "work-load balance: %g  max=%g avg=%g PE0=%g\n",
	      maxt / (sumt / NTask), maxt, sumt / NTask, timetreelist[0]);
      fprintf(FdTimings, "particle-load balance: %g\n", plb_max);
      fprintf(FdTimings, "max. nodes: %d, filled: %g\n", maxnumnodes,
	      maxnumnodes / (All.TreeAllocFactor * All.MaxPart));
      fprintf(FdTimings, "part/sec=%g | %g  ia/part=%g (%g)\n", ntot / (sumt + 1.0e-20),
	      ntot / (maxt * NTask), ((double) (costtotal)) / ntot, ((double) ewaldtot) / ntot);
      fprintf(FdTimings, "\n");

      fflush(FdTimings);

      All.CPU_TreeWalk += sumt / NTask;
      All.CPU_Imbalance += sumimbalance / NTask;
      All.CPU_CommSum += sumcomm / NTask;
    }


  // DY: external potential, following Volker's suggestion 
  double mstar_now = Mdisk;
  double fstar_now = Fdisk;
  double fx_disk = 0.;
  double fy_disk = 0.;
  double fz_disk = 0.;
  int itask=0;
  int otask=0;

  for(i=0; i < NumPart; i++)
     if(P[i].Ti_endstep==All.Ti_Current)
     {
        double rad = 0.;
	if(P[i].ID!=id0){
           double dx = P[i].Pos[0]-cenx;
           double dy = P[i].Pos[1]-ceny;
           double dz = P[i].Pos[2]-cenz;

           rad = sqrt(dx*dx + dy*dy + dz*dz);
           double cenjx=HALO_JX;
           double cenjy=HALO_JY;
           double cenjz=HALO_JZ;
	   double cenj = sqrt(cenjx*cenjx + cenjy*cenjy + cenjz*cenjz);
	   double costh = (cenjx*dx+cenjy*dy+cenjz*dz)/(cenj*rad);
	   double z = rad*costh; 
	   double R = rad*sqrt(1-costh*costh); 
	   if(fabs(costh)>=1){ z=0.001; printf("~~~~~~~~~~~~~~~~~~~~~~~~ costh==1\n");}
           double sforz=AccRz(R, z, 0, mstar_now, fstar_now);// -accz along the j direction
	   double forz[3] = {sforz*cenjx/cenj, sforz*cenjy/cenj, sforz*cenjz/cenj};
           double sforR=AccRz(R, z, 1, mstar_now, fstar_now);// -accR along the R direction
	   double vecR[3] = {dx-(z*cenjx/cenj), dy-(z*cenjy/cenj), dz-(z*cenjz/cenj)};
	   double forR[3] = {sforR*vecR[0]/R, sforR*vecR[1]/R, sforR*vecR[2]/R};
	   double fadia=0;
	   if(All.Time<FORCE_A && All.Time>=START_A) fadia=(All.Time-START_A)/(FORCE_A-START_A);
	   else fadia=1;
           if(All.Time>=START_A){
              P[i].GravAccel[0] -= fadia*(forz[0]+forR[0]) * All.G;
              P[i].GravAccel[1] -= fadia*(forz[1]+forR[1]) * All.G;
              P[i].GravAccel[2] -= fadia*(forz[2]+forR[2]) * All.G;
              //Recording the force back reaction from all other particles on the disk (excluding the sink)
              fx_disk += P[i].Mass*fadia*(forz[0]+forR[0]) * All.G;
              fy_disk += P[i].Mass*fadia*(forz[1]+forR[1]) * All.G;
              fz_disk += P[i].Mass*fadia*(forz[2]+forR[2]) * All.G; 
           }
	}
        else{
           itask=ThisTask;
           //printf("loop %d, task %d\n", i, ThisTask);
           //Renormalizing the sink particle's acceleration to sink-disk combination
           P[i].GravAccel[0] *= P[i].Mass / (P[i].Mass + (mstar_now / fstar_now) * All.HubbleParam/1e10);
           P[i].GravAccel[1] *= P[i].Mass / (P[i].Mass + (mstar_now / fstar_now) * All.HubbleParam/1e10);
           P[i].GravAccel[2] *= P[i].Mass / (P[i].Mass + (mstar_now / fstar_now) * All.HubbleParam/1e10);
        }
     }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allreduce(&itask, &otask, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  //Applying the disk back reaction forces and adding it to the sink's accelerations
  double Fx_disk;
  double Fy_disk;
  double Fz_disk; 
  MPI_Allreduce(&fx_disk, &Fx_disk, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&fy_disk, &Fy_disk, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&fz_disk, &Fz_disk, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  if (ThisTask == otask){
      //printf("This task = %d \n", ThisTask);
      for(i=0; i < NumPart; i++)
          if(P[i].Ti_endstep==All.Ti_Current){
            //printf("This task = %d \n", ThisTask);
            if (P[i].ID == id0){
              P[i].GravAccel[0] += Fx_disk / (P[i].Mass + (mstar_now / fstar_now) * All.HubbleParam/1e10);
              P[i].GravAccel[1] += Fy_disk / (P[i].Mass + (mstar_now / fstar_now) * All.HubbleParam/1e10);
              P[i].GravAccel[2] += Fz_disk / (P[i].Mass + (mstar_now / fstar_now) * All.HubbleParam/1e10);
            }
          }
  }

  free(nrecv);
  free(ewaldlist);
  free(numnodeslist);
  free(costtreelist);
  free(timecommlist);
  free(timetreelist);
}



/*! This function sets the (comoving) softening length of all particle
 *  types in the table All.SofteningTable[...].  We check that the physical
 *  softening length is bounded by the Softening-MaxPhys values.
 */
void set_softenings(void)
{
  int i;

  if(All.ComovingIntegrationOn)
    {
      if(All.SofteningGas * All.Time > All.SofteningGasMaxPhys)
        All.SofteningTable[0] = All.SofteningGasMaxPhys / All.Time;
      else
        All.SofteningTable[0] = All.SofteningGas;
      
      if(All.SofteningHalo * All.Time > All.SofteningHaloMaxPhys)
        All.SofteningTable[1] = All.SofteningHaloMaxPhys / All.Time;
      else
        All.SofteningTable[1] = All.SofteningHalo;
      
      if(All.SofteningDisk * All.Time > All.SofteningDiskMaxPhys)
        All.SofteningTable[2] = All.SofteningDiskMaxPhys / All.Time;
      else
        All.SofteningTable[2] = All.SofteningDisk;
      
      if(All.SofteningBulge * All.Time > All.SofteningBulgeMaxPhys)
        All.SofteningTable[3] = All.SofteningBulgeMaxPhys / All.Time;
      else
        All.SofteningTable[3] = All.SofteningBulge;
      
      if(All.SofteningStars * All.Time > All.SofteningStarsMaxPhys)
        All.SofteningTable[4] = All.SofteningStarsMaxPhys / All.Time;
      else
        All.SofteningTable[4] = All.SofteningStars;
      
      if(All.SofteningBndry * All.Time > All.SofteningBndryMaxPhys)
        All.SofteningTable[5] = All.SofteningBndryMaxPhys / All.Time;
      else
        All.SofteningTable[5] = All.SofteningBndry;
      // DY for the sink particle
      if(All.SofteningSink * All.Time > All.SofteningSinkMaxPhys)
        All.SofteningTable[6] = All.SofteningSinkMaxPhys / All.Time;
      else
        All.SofteningTable[6] = All.SofteningSink;
      
    }
  else
    {
      All.SofteningTable[0] = All.SofteningGas;
      All.SofteningTable[1] = All.SofteningHalo;
      All.SofteningTable[2] = All.SofteningDisk;
      All.SofteningTable[3] = All.SofteningBulge;
      All.SofteningTable[4] = All.SofteningStars;
      All.SofteningTable[5] = All.SofteningBndry;
      // DY
      All.SofteningTable[6] = All.SofteningSink;
    }

  for(i = 0; i < 7; i++)
    All.ForceSoftening[i] = 2.8 * All.SofteningTable[i];

  All.MinGasHsml = All.MinGasHsmlFractional * All.ForceSoftening[0];
}


/*! This function is used as a comparison kernel in a sort routine. It is
 *  used to group particles in the communication buffer that are going to
 *  be sent to the same CPU.
 */
int grav_tree_compare_key(const void *a, const void *b)
{
  if(((struct gravdata_index *) a)->Task < (((struct gravdata_index *) b)->Task))
    return -1;

  if(((struct gravdata_index *) a)->Task > (((struct gravdata_index *) b)->Task))
    return +1;

  return 0;
}
