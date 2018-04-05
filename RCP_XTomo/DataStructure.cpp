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
struct Particle Bottom;
struct Particle Top;
struct Particle Lateral;
//struct Particle SBox[27];
struct Rod Rods[N_Rods];
struct Rod Rod0;

int Gi, Gj, Gk, Gm;
double OuterRadius, InnerRadius;
double Volume_Factor;
double upper, lower, lam, OuterRadius0;
double L;

double Density_Roof=0.0;

double C_Volume =Pi/4.0 * Cylinder_H *Cylinder_D*Cylinder_D;
double critic_radius = pow((Dis_Density*C_Volume)/N_Particles *3.0/4.0/Pi, 0.33333333);

FILE *fp, *fp0;
