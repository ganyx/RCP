// InputParameters.h

#define Pi 3.14159265

#define NMAX_Iteration 10000000		// The maximum number of iterations
#define NMAX_Movement 100			// The maximum number of movements in one iteration

#define N_Particles 8421			// The number of particles in the unit cell
#define Dis_Density 0.617
#define N_Case 1
#define ContractionRate 0.02

#define Cylinder_H 4.63	// Height of the cylinder container
#define Cylinder_D 4.89	// Diameter of the cylinder container
#define Particle_D 2.3	// Diameter of the particles

#define Tolerance 1E-10

const int N_Rods = 2*N_Particles;	// The number of overlaps in the queue
const int Fresh_Freq = 10*N_Rods/NMAX_Movement;

double R_Relative_Width = 0.0;
double RelativeDensity_RMax = 0.0;
