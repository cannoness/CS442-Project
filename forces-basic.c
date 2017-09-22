#include <mpi.h>
#include <math.h>
#include <string.h>

#include "graphics.h"
#include "nbody.h"

#define G 6.67428e-11

double ComputeForces( Particle myparticles[], Particle others[], 
		      ParticleF pf[], int npart, int rank )
{
  double max_f, rmin;
  int i, j;

  max_f = 0.0;
  for (i=0; i<npart; i++) {
    double xi, yi, zi, mi, rx, ry, rz, mj, r, fx, fy, fz, f;
    rmin = 10.0;
    xi   = myparticles[i].x;
    yi   = myparticles[i].y;
    zi   = myparticles[i].z;
    fx   = 0.0;
    fy   = 0.0;
    fz   = 0.0;
    mi   = myparticles[i].mass;
#ifdef GRAPHICS
    if (rank == 0) {
        enqueueParticle(others[i].x, others[i].y, others[i].z,
		        others[i].r, others[i].g, others[i].b);
    }
#endif
    for (j=0; j<npart; j++) {
      rx = xi - others[j].x;
      ry = yi - others[j].y;
      rz = zi - others[j].z;
      mj = others[j].mass;
      r  = rx * rx + ry * ry + rz * rz;
      f = -G * mi * mj / r;
      /* ignore overlap and same particle */
      if (r < 0.1) continue;
      if (r < rmin) rmin = r;
      /* compute forces */
      r  = sqrt(r);
      fx += f * rx / r;
      fy += f * ry / r;
      fz += f * rz / r;
    }
    pf[i].fx += fx;
    pf[i].fy += fy;
    pf[i].fz += fz;
    /* Compute a rough estimate of (1/m)|df / dx| */
    fx		      = sqrt(fx*fx + fy*fy + fz*fz)/rmin;
    if (fx > max_f) max_f = fx;
  }
  return max_f;
}

double ComputeAllForces(Particle *particles, MPI_Datatype particletype,
			ParticleF *pf, MPI_Datatype forcetype, 
			int size, int npart, int rank,
			MPI_Comm commring, int right, int left)
{
    Particle  sendbuf[MAX_PARTICLES],     /* Pipeline buffers */
	    recvbuf[MAX_PARTICLES];
    MPI_Request request[2];
    MPI_Status statuses[2];
    int pipe;
    double max_f = 0.0, max_f_seg;
    /* Load the initial sendbuffer */
    memcpy( sendbuf, particles, npart * sizeof(Particle) );
#ifdef GRAPHICS
    if (rank == 0) {
    	startParticles();
    }
#endif
    max_f = 0.0;
    for (pipe=0; pipe<size; pipe++) {
        if (pipe != size-1) {
	    MPI_Isend( sendbuf, npart, particletype, right, pipe, 
		       commring, &request[0] );
	    MPI_Irecv( recvbuf, npart, particletype, left,  pipe, 
		       commring, &request[1] );
	}
	max_f_seg = ComputeForces( particles, sendbuf, pf, npart, rank );
	if (max_f_seg > max_f) max_f = max_f_seg;
	/* Push pipe */
	if (pipe != size-1) {
	    MPI_Waitall( 2, request, statuses );
	}
	memcpy( sendbuf, recvbuf, npart * sizeof(Particle) );
    }
#ifdef GRAPHICS
    if (rank == 0) {
        finishParticles();
    }
#endif
    return max_f;
}
