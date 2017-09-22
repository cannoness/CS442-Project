#ifndef _NBODY_H_
#define _NBODY_H_

typedef struct Particle {
	double x, y, z;
#ifdef GRAPHICS
	double r, g, b;
#endif
	double mass;
	int fixed;
	int interactRange[3][2];
} Particle;
#define PARTICLE_SIZE (sizeof(Particle)/sizeof(double))

typedef struct ParticleV {
	double xold, yold, zold;
} ParticleV;
#define VELOCITY_SIZE (sizeof(ParticleV)/sizeof(double))

typedef struct ParticleF {
	double fx, fy, fz;
} ParticleF;
#define FORCE_SIZE (sizeof(ParticleF)/sizeof(double))

typedef struct Spring {
	int to, from;
	double initLen;
} Springs;



#define INITIAL_POINT_SIZE 5.0
#define MAX_PARTICLES 128

#endif /* _NBODY_H_ */
