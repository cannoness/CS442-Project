#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "nbody.h"

#ifdef GRAPHICS
#include "graphics.h"
#endif

Particle  particles[MAX_PARTICLES];   /* Particles on this node */
ParticleV pv[MAX_PARTICLES];          /* Particle velocity */
ParticleF pf[MAX_PARTICLES];          /* Particle force */
Springs springs[MAX_PARTICLES/2];
int num_springs=MAX_PARTICLES-3;

double ks = 2.25f;
float point_size = INITIAL_POINT_SIZE;

float calculateDistance(Particle a, Particle b){
	float dist=0;
	float dx = a.x-b.x;
	float dy = a.y-b.y;
	float dz = a.z-b.z;
	dist = sqrt(dx*dx+dy*dy+dz*dz);
	return dist;
}

double ComputeNewPos( Particle particles[], ParticleV pv[], ParticleF pf[], 
		      int npart, double max_f, MPI_Comm commring )
{
	
  int i;
  double      a0, a1, a2;
  static      double dt_old = 0.02, dt = 0.02;
  double      dt_est, new_dt, dt_new;
  /* integation is a0 * x^+ + a1 * x + a2 * x^- = f / m */
  a0	 = 2.0 / (dt * (dt + dt_old));
  a2	 = 2.0 / (dt_old * (dt + dt_old));
  a1	 = -(a0 + a2);      /* also -2/(dt*dt_old) */
	/* do the spring forces first */
  
  for (i=0; i<npart; i++) {
		double xi, yi, zi;
		double dx,dy,dz,len;
		dx= particles[i+1].x-particles[i].x;
		dy= particles[i+1].y-particles[i].y;
		dz= particles[i+1].z-particles[i].z;
		//len= sqrt(dx*dx+dy*dy+dz*dz);
		len= calculateDistance(particles[i],particles[i+1]);
		
		pf[i].fx+=ks*(len-springs[i].initLen)*dx / len;
		pf[i+1].fx-=ks*(len-springs[i].initLen)*dx / len;
		pf[i].fy+=ks*(len-springs[i].initLen)*dy / len;
		pf[i+1].fy-= ks*(len-springs[i].initLen)*dy / len;
		pf[i].fz+=ks*(len-springs[i].initLen)*dz / len;
		pf[i+1].fz-=ks*(len-springs[i].initLen)*dz / len;
		
    /* Very, very simple leapfrog time integration.  We use a variable 
       step version to simplify time-step control. */
     //if we're not fixed in space
		xi	           = particles[i].x;
		yi	           = particles[i].y;
		zi	           = particles[i].z;
		if(particles[i].fixed==0){	
			particles[i].x = (pf[i].fx - a1 * xi - a2 * pv[i].xold) / a0;
			particles[i].y = (pf[i].fy - a1 * yi - a2 * pv[i].yold) / a0;
			particles[i].z = (pf[i].fz - a1 * zi - a2 * pv[i].zold) / a0;
			
			pv[i].xold     = xi;
			pv[i].yold     = yi;
			pv[i].zold     = zi;
		}
		pf[i].fx       = 0;
		pf[i].fy       = 0;
		pf[i].fz       = 0;
	
  }

#if !defined(GRAPHICS)
  /* We don't use dynamic timesteps with graphics, as it's not as
   * pretty. If we're not doing graphics, we go for accuracy, and 
   * Recompute a time step. Stability criteria is roughly 
   * 2/sqrt(1/m |df/dx|) >= dt.  We leave a little room */
 /* dt_est = 1.0/sqrt(max_f);
  /* Set a minimum: */
  /*if (dt_est < 1.0e-6) dt_est = 1.0e-6;
  MPI_Allreduce( &dt_est, &dt_new, 1, MPI_DOUBLE, MPI_MIN, commring );
  /* Modify time step */
/*  if (dt_new < dt) {
    dt_old = dt;
    dt     = dt_new;
  }
  else if (dt_new > 4.0 * dt) {
    dt_old = dt;
    dt    *= 2.0;
  }*/
#endif

  return dt_old;
}

void InitParticles( Particle particles[], ParticleV pv[], ParticleF pf[], int npart )
{
	int i;
	particles[0].fixed = 1;
	particles[0].x = 0.0f;
	particles[0].y = 0.0f;
	particles[0].z = 0.0f;
	for (i = 0; i <= 31; i++){	
		float newpos = 0.0f;
		newpos +=i*0.25f;
//		particles[i].color = i%8;
		particles[i].fixed = 0; 
		particles[i].x = newpos; 
		particles[i].y = 0.0f;
		particles[i].z = 0.0f;
		particles[i].mass =  1.0f;
		pv[i].xold	  = particles[i].x;
		pv[i].yold	  = particles[i].y;
		pv[i].zold	  = particles[i].z;
		pf[i].fx	  = 0;
		pf[i].fy	  = 0;
		pf[i].fz	  = 0;
	#ifdef GRAPHICS
		particles[i].r	  = drand48();
		particles[i].g	  = drand48();
		particles[i].b	  = drand48();
	#endif
			
		particles[i].interactRange[0][0]=32;
		particles[i].interactRange[0][1]=63;
		particles[i].interactRange[1][0]=64;
		particles[i].interactRange[1][1]=95;
		particles[i].interactRange[2][0]=96;
		particles[i].interactRange[2][1]=127; 
	}
	particles[31].fixed = 1;
	
	//line number 2
	particles[32].fixed = 1;
	particles[32].x = 0.0f;
	particles[32].y = -0.25f;
	particles[32].z = 0.0f;
	for (i = 32; i <= 63; i++){	
		float newpos = 0.0f;
		newpos +=(i-32)*0.25f;
		//particles[i].color = i%8;
		particles[i].fixed = 0;
		particles[i].x = newpos; 
		particles[i].y = -0.25f;
		particles[i].z = 0.0f;
	   	particles[i].interactRange[0][0]=0;
		particles[i].interactRange[0][1]=31;
		particles[i].interactRange[1][0]=64;
		particles[i].interactRange[1][1]=95;
		particles[i].interactRange[2][0]=96;
		particles[i].interactRange[2][1]=127;
		pv[i].xold	  = particles[i].x;
		pv[i].yold	  = particles[i].y;
		pv[i].zold	  = particles[i].z;
		pf[i].fx	  = 0;
		pf[i].fy	  = 0;
		pf[i].fz	  = 0;
#ifdef GRAPHICS
	particles[i].r	  = drand48();
	particles[i].g	  = drand48();
	particles[i].b	  = drand48();
#endif
		
	}
	particles[63].fixed = 1;

	//line number 3
	particles[64].fixed = 1;
	particles[64].y = 0.0f;
	particles[64].y = 0.25f;
	particles[64].z = 0.0f;
	for (i = 65; i <= 95; i++){	
		float newpos = 0.0f;
		newpos +=(i-64)*0.25f;
//		particles[i].color = i%8;
	#ifdef GRAPHICS
	particles[i].r	  = drand48();
	particles[i].g	  = drand48();
	particles[i].b	  = drand48();
#endif
		particles[i].fixed = 0;
		particles[i].x = newpos; 
		particles[i].y = 0.25f;
		particles[i].z = 0.0f;
		particles[i].interactRange[0][0]=0;
		particles[i].interactRange[0][1]=31;
		particles[i].interactRange[1][0]=32;
		particles[i].interactRange[1][1]=63;
		particles[i].interactRange[2][0]=96;
		particles[i].interactRange[2][1]=127; 
		pv[i].xold	  = particles[i].x;
		pv[i].yold	  = particles[i].y;
		pv[i].zold	  = particles[i].z;
		pf[i].fx	  = 0;
		pf[i].fy	  = 0;
		pf[i].fz	  = 0;
#ifdef GRAPHICS
	particles[i].r	  = drand48();
	particles[i].g	  = drand48();
	particles[i].b	  = drand48();
#endif
		
	}
	particles[95].fixed = 1;
	
	//line number 4
	particles[96].fixed = 1;
	particles[96].x = 0.0f;
	particles[96].y = 0.5f;
	particles[96].z = 0.0f;
	for (i = 96; i <= 127; i++){	
		float newpos = 0.0f;
		newpos +=(i-96)*0.25f;
//		particles[i].color = i%8;
	
		particles[i].fixed = 0;
		particles[i].x = newpos; 
		particles[i].y = 0.5f-i/96;
		particles[i].z = 0.0f;
		particles[i].interactRange[0][0]=0;
		particles[i].interactRange[0][1]=31;
		particles[i].interactRange[1][0]=32;
		particles[i].interactRange[1][1]=63;
		particles[i].interactRange[2][0]=64;
		particles[i].interactRange[2][1]=95;
		pv[i].xold	  = particles[i].x;
		pv[i].yold	  = particles[i].y;
		pv[i].zold	  = particles[i].z;
		pf[i].fx	  = 0;
		pf[i].fy	  = 0;
		pf[i].fz	  = 0;
#ifdef GRAPHICS
	particles[i].r	  = drand48();
	particles[i].g	  = drand48();
	particles[i].b	  = drand48();
#endif
		
	}
	particles[127].fixed = 1;

	//find a way to save the fact that these particles interact... 
	//maybe make a structure called line 1 and line 2
	//that way I know where to start looking at point distances
	//connect every two particles with a spring
	for (i = 0; i <= 31; i++){
		springs[i].from = i;
		springs[i].to = i+1;	
		springs[i].initLen = 0.25; //calculateDistance(particles[i],particles[i+1]);
		//if (particles[i].fixed ==false)
		//particles[i].position[0]-=0.4*springs[i].initLen;
//	num_springs++;
	}
	for (i = 32; i <= 63; i++){
		springs[i].from = i;
		springs[i].to = i+1;	
		springs[i].initLen = 0.25; //calculateDistance(particles[i],particles[i+1]);	
		//if (particles[i].fixed ==false)
		//particles[i].position[0]-=0.4*springs[i].initLen;
//	num_springs++;
		
	}
	for (i = 64; i <= 95; i++){
		springs[i].from = i;
		springs[i].to = i+1;	
		springs[i].initLen = 0.25; //calculateDistance(particles[i],particles[i+1]);	
		//if (particles[i].fixed ==false)
		//particles[i].position[0]-=0.4*springs[i].initLen;
//	num_springs++;
		
	}
	for (i = 96; i <= 127; i++){
		springs[i].from = i;
		springs[i].to = i+1;	
		springs[i].initLen = 0.25; //calculateDistance(particles[i],particles[i+1]);	
		//if (particles[i].fixed ==false)
		//particles[i].position[0]-=0.4*springs[i].initLen;
//	num_springs++;
		
	}
	
}

double ComputeAllForces(Particle [], MPI_Datatype, 
			ParticleF [], MPI_Datatype,
			int, int, int,
			MPI_Comm, int, int);

void ComputeMainLoop(Particle *particles, ParticleV *pv, ParticleF *pf, int npart, 
		     int size, int rank, int cnt)
{
    MPI_Datatype particletype, forcetype;
    int         left, right, periodic;
    MPI_Comm    commring;
    double      sim_t;                      /* Simulation time */
    double      time;                       /* Computation time */

    /* Get the best ring in the topology */
    periodic = 1;
    MPI_Cart_create( MPI_COMM_WORLD, 1, &size, &periodic, 1, &commring );
    MPI_Cart_shift( commring, 0, 1, &left, &right );

    MPI_Type_contiguous( PARTICLE_SIZE, MPI_DOUBLE, &particletype );
    MPI_Type_commit( &particletype );

    MPI_Type_contiguous( FORCE_SIZE, MPI_DOUBLE, &forcetype );
    MPI_Type_commit( &forcetype );

    time = MPI_Wtime();
    sim_t = 0.0;

    while (cnt--) {
	double max_f, dt_old;
	max_f = ComputeAllForces(particles, particletype, 
				 pf, forcetype, 
				 size, npart, rank,
				 commring, right, left);

	/* Once we have the forces, we compute the changes in position */
	dt_old = ComputeNewPos( particles, pv, pf, npart, max_f, commring ); 
	sim_t += dt_old;
	/* If we're running continuously, see if we've been directed to stop */
	if (cnt < 0) {
	    int done = 0, doneout = 0;
	    MPI_Allreduce(&done, &doneout, 1, MPI_INT, MPI_LOR, commring);
	    if (doneout) {
		break;
	    }
	}
    }
    time = MPI_Wtime() - time;
    if (rank == 0) {
	printf( "Computed %d particles over %f seconds in %f seconds\n", npart * size, sim_t, time );
    }
    return;
}

#ifdef GRAPHICS
/* We use these ugly global variables if there's graphics going on 
 * since GLUT works by callback functions. */
MPI_Datatype g_particletype, g_forcetype;
static int g_npart;
static int g_size, g_rank;
static int g_left, g_right;
static MPI_Comm g_commring;
static double g_sim_t;
static int EndSimulation;
void GraphicsMainLoop(Particle *particles, ParticleV *pv, ParticleF *pf, int npart, 
		     int size, int rank)
{
    int         periodic;
    double      time;                       /* Computation time */

    g_size = size;
    g_rank = rank;
    g_npart = npart;

    /* Get the best ring in the topology */
    periodic = 1;
    MPI_Cart_create( MPI_COMM_WORLD, 1, &g_size, &periodic, 1, &g_commring );
    MPI_Cart_shift( g_commring, 0, 1, &g_left, &g_right );

    MPI_Type_contiguous( PARTICLE_SIZE, MPI_DOUBLE, &g_particletype );
    MPI_Type_commit( &g_particletype );

    MPI_Type_contiguous( FORCE_SIZE, MPI_DOUBLE, &g_forcetype );
    MPI_Type_commit( &g_forcetype );
    g_sim_t = 0.0;
    GraphicsStart();
}


void SimulationDone() 
{
	EndSimulation = 1;
} 

void SimulationDisplay() 
{
	double max_f, dt_old, endout = 0;
	GraphicsDisplayStart();
	max_f = ComputeAllForces(particles, g_particletype, 
				 pf, g_forcetype, 
				 g_size, g_npart, g_rank,
				 g_commring, g_right, g_left);

	/* Once we have the forces, we compute the changes in position */
	dt_old = ComputeNewPos( particles, pv, pf, g_npart, max_f, g_commring ); 
	g_sim_t += dt_old;
	GraphicsDisplayFinish();
	/* If we're running continuously, see if we've been directed to stop */
	MPI_Allreduce(&EndSimulation, &endout, 1, MPI_INT, MPI_LOR, g_commring);
	if (endout) {
	    GraphicsEnd();
	}
}
#endif

int main( int argc, char *argv[] )
{
    int         rank, size, npart, i, j,
	        offset;                     /* location of local particles */
    int         totpart,                    /* total number of particles */
	        cnt;                        /* number of times in loop */

    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &size );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    srand48(rank);

#ifdef GRAPHICS
    initParticleDisplay(argc, argv);
#endif

    if (argc < 2) { 
	fprintf( stderr, "Usage: %s n\n", argv[0] );
	MPI_Abort( MPI_COMM_WORLD, 1 );
    }
    totpart = atoi(argv[1]);
    
    npart = totpart/size;
    if (npart > MAX_PARTICLES) {
	fprintf( stderr, "%d is too many; max is %d\n", 
		 totpart, MAX_PARTICLES*size );
	MPI_Abort( MPI_COMM_WORLD, 1 );
    }

    if ( (totpart % size) != 0) {
	fprintf(stderr, "Warning: Number of particles truncated to %d to distribute evenly.\n",
		npart * size);
    }
    totpart = npart*size;


    /* Generate the initial values */
    InitParticles( particles, pv, pf, npart);

#ifdef GRAPHICS
    if (rank == 0) {
	GraphicsMainLoop(particles, pv, pf, npart, size, rank);
    } else {
	ComputeMainLoop(particles, pv, pf, npart, size, rank, -1);
    }
#else
    cnt = 100;
    ComputeMainLoop(particles, pv, pf, npart, size, rank, cnt);
#endif

    MPI_Finalize();

    return 0;
}

