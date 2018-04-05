// DataStructure.cpp

struct Particle
{
	long int number;
	double x;
	double y;
	double z;
	double r;
};

struct Rod
{
	double distance;
	struct Particle *P[2];
	int N_Box;
};

struct Particle Particles[N_Particles];
struct Particle SBox[27];
struct Rod Rods[N_Rods];
struct Rod Rod0;

int Gi, Gj, Gk, Gm;
double OuterRadius, InnerRadius;
double Volume_Factor;
double upper, lower, lam, OuterRadius0;

FILE *fp, *fp0;
