// InputParameters.h

#define Pi 3.14159265

#define NMAX_Iteration 10000000		// The maximum number of iterations
#define NMAX_Movement 100			// The maximum number of movements in one iteration

#define N_Particles 5000			// The number of particles in the unit cell
#define N_Case 5
#define ContractionRate 0.01

#define L 1.0

#define Tolerance 0.00000001

const int N_Rods = 2*N_Particles;	// The number of overlaps in the queue
const int Fresh_Freq = 10*N_Rods/NMAX_Movement;

double R_Relative_Width = 0.0;
double RelativeDensity_RMax = 0.0;
