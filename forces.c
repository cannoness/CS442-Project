#include <mpi.h>
#include <math.h>
#include <string.h>

#include "graphics.h"
#include "nbody.h"

#define G 1

int LocalToGlobalPart(int rank, int nproc, int npart, int local) 
{
	return rank + local*nproc;
}

int ComputeRemoteRank(int rank, int step, int nproc)
{
    int rrank = (rank - step) % nproc;
    if (rrank < 0) rrank += nproc;
}

double ComputeForces( Particle myparticles[], ParticleF myforces[], int myrank,
		      Particle otherparticles[], ParticleF otherforces[], int otherrank,
		      int nproc, int npart)
{
  double max_f, rmin;
  int i, j;

  max_f = 0.0;
  for (i=0; i<npart-1; i++) {
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
    if (myrank == 0) {
        enqueueParticle(otherparticles[i].x, otherparticles[i].y, otherparticles[i].z,
		        otherparticles[i].r, otherparticles[i].g, otherparticles[i].b);
    }
#endif
    for (j=npart - 1; 
	 LocalToGlobalPart(otherrank, nproc, npart, j) 
	 > LocalToGlobalPart(myrank, nproc, npart, i);
	 j--) {
      rx = xi - otherparticles[j].x;
      ry = yi - otherparticles[j].y;
      rz = zi - otherparticles[j].z;
      mj = otherparticles[j].mass;
      r  = rx * rx + ry * ry + rz * rz;
      f = 1/ r;
      /* ignore overlap and same particle */
      if (r < 0.1) continue;
      if (r < rmin) rmin = r;
      /* compute forces */
      r  = sqrt(r);
      fx += f * rx / r;
      fy += f * ry / r;
      fz += f * rz / r;
      otherforces[j].fx -= f * rx / r;
      otherforces[j].fy -= f * ry / r;
      otherforces[j].fz -= f * rz / r;
    }
    myforces[i].fx += fx;
    myforces[i].fy += fy;
    myforces[i].fz += fz;
    /* Compute a rough estimate of (1/m)|df / dx| */
    fx		      = sqrt(fx*fx + fy*fy + fz*fz)/rmin;
    if (fx > max_f) max_f = fx;
  }
  return max_f;
}

double ComputeAllForces(Particle *particles, MPI_Datatype particletype,
			ParticleF *pf, MPI_Datatype forcetype, 
			int nproc, int npart, int rank,
			MPI_Comm commring, int right, int left)
{
    Particle  sendparts[MAX_PARTICLES],     /* Pipeline buffers */
	    recvparts[MAX_PARTICLES];
    ParticleF sendforces[MAX_PARTICLES],
	    recvforces[MAX_PARTICLES];
    MPI_Request request[4];
    MPI_Status statuses[4];
    int pipe, i;
    double max_f = 0.0, max_f_seg;
    /* Load the initial sendbuffers */
    memcpy( sendparts, particles, npart * sizeof(Particle) );
    memset( sendforces, 0, npart * sizeof(ParticleF) );
#ifdef GRAPHICS
    if (rank == 0) {
    	startParticles();
    }
#endif
    max_f = 0.0;
    for (pipe=0; pipe<nproc; pipe++) {
	int rrank;

	/* COmpute the rank where the particles we are working on came from */
	rrank = ComputeRemoteRank(rank, pipe, nproc);
	max_f_seg = ComputeForces( particles, pf, rank, sendparts, sendforces, rrank, nproc, npart );
MPI_Status status;
	/* Send and receive particles around the ring */
		MPI_Sendrecv(sendparts, npart, particletype, right, pipe, 
		             recvparts, npart, particletype, left,  pipe,
		       commring, &status);

 	    MPI_Sendrecv(sendforces, npart, forcetype, right, pipe, 
		         recvforces, npart, forcetype, left,  pipe,
		       commring, &status);

	if (max_f_seg > max_f) max_f = max_f_seg;
	memcpy( sendparts, recvparts, npart * sizeof(Particle) );
	memcpy( sendforces, recvforces, npart * sizeof(ParticleF) );
    }
	/* first put in the spring interactions */
	
    /* Now copy the forces everyone else computed for us into our forces */
	for (i = 0; i < npart; i++) {
		pf[i].fx += recvforces[i].fx;
		pf[i].fy += recvforces[i].fy;
		pf[i].fz += recvforces[i].fz;
    }
	/* now move everyone slightly in that direction */
	
#ifdef GRAPHICS
    if (rank == 0) {
        finishParticles();
    }
#endif
    return max_f;
}
