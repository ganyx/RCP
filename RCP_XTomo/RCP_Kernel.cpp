// RCP_Kernel.cpp

// Randomly distribute the intiatal packing in the cylinder container
void Initiation()
{
	int i;
	extern double SortRods(int np);
	extern bool InContainer(int i);
	Volume_Factor = 0.0;

	fprintf(fp," The Number of Particles is: %d\n",N_Particles);

	if(Cylinder_H < Cylinder_D) L = Cylinder_D;
	else L= Cylinder_H;

	srand(time(NULL));
	for(i=0;i<N_Particles;)
	{
		Particles[i].number=i;
		Particles[i].x = L*rand()/RAND_MAX - L/2.0;
		Particles[i].y = L*rand()/RAND_MAX - L/2.0;
		Particles[i].z = L*rand()/RAND_MAX - L/2.0;
		Particles[i].r = 1.0;

		// If the particle locates inside the cylinder container, generate the next one.
		if(InContainer(i)) 
		{
			Volume_Factor += pow(Particles[i].r, 3);
			i++;
		}
	}

	// The OuterRadius in unequal sizes particles' system
	// 	means the mean value of curent radius distribution.
	Volume_Factor = Volume_Factor /C_Volume;
	OuterRadius = pow(3.0/4.0/Pi/Volume_Factor, 0.33333333);
	OuterRadius0 = OuterRadius;
	lam=1.0;
	SortRods(-1);
	return;
}

// Check whether initial points is in the container 
bool InContainer(int i)
{
	double in_r = sqrt(Particles[i].x*Particles[i].x + Particles[i].y*Particles[i].y);
	if(abs(Particles[i].z) < 0.9* Cylinder_H/2.0 
		&& in_r < 0.9* Cylinder_D/2.0)
	return 1;
	else return 0;
}

// Queue the most overlaps between particles in the current system.
double SortRods(int np)
{
	int sri,srj,NP0,NP1;
	extern double NearestDist(int m, int n);
	extern bool CheckBoundary(int i);
	extern void InsertRod(int i);
	extern void DeleteRods(int nnp);

	if(np==-1) // Initial the Rods[]
	{
		for(sri=0;sri<N_Rods;sri++)
		{
			Rods[sri].distance = L;
		}
		for(sri=0;sri<N_Particles;sri++)
		{CheckBoundary(sri);}

		for(sri=0;sri<N_Particles;sri++)
		{
			for(srj=sri+1;srj<N_Particles;srj++)
			{
				NearestDist(sri, srj);
			}
		}
	}
	else if(np>=0) // Update the change by moving Rod[np]
	{
		NP0=(*Rods[np].P[0]).number;
		NP1=(*Rods[np].P[1]).number;

		DeleteRods(NP0);
		DeleteRods(NP1);

		NearestDist(NP0, NP1);

		for(sri=0;sri<N_Particles;sri++)
		{
			if(!(sri==NP0 || sri==NP1))
			{
				NearestDist(NP0, sri);
				NearestDist(NP1, sri);
			}
		}
	}

	InnerRadius = Rods[0].distance;
	return Rods[0].distance;
}

// Insert Rod0 to the queue Rods[N_Rods]
void InsertRod()
{
	int i,j;

	for(i=N_Rods-1;i>=0;i--)
	{
		if(Rod0.distance>=Rods[i].distance) break;
		continue;
	}
	if(i<N_Rods-1)
	{
		for(j=N_Rods-1;j>i+1;j--)
		{
			Rods[j].distance=Rods[j-1].distance;
			Rods[j].P[0]=Rods[j-1].P[0];
			Rods[j].P[1]=Rods[j-1].P[1];
		}
			Rods[i+1].distance=Rod0.distance;
			Rods[i+1].P[0]=Rod0.P[0];
			Rods[i+1].P[1]=Rod0.P[1];
	}
	return;
}

void DeleteRods(int np)
{
	int i,j;
	for(i=0;i<N_Rods;i++)
	{
		if((*Rods[i].P[0]).number==np || (*Rods[i].P[1]).number==np)
		{
			for(j=i;j<N_Rods-1;j++)
			{
				Rods[j].distance=Rods[j+1].distance;
				Rods[j].P[0]=Rods[j+1].P[0];
				Rods[j].P[1]=Rods[j+1].P[1];
			}
			Rods[N_Rods-1].distance = L;
		}
	}
	return;
}

bool CheckBoundary(int i) // if particle i is overlap with the boundaries, push it inside the container.
{
	extern double NearestDist(int m, int n);
	double Dgap, lat_r;
	bool move_flag=0;
	int sri;

	if(i>=0)
	{
	//Check boundaries
	//*********************************
	// Case 1: Bottom
	Dgap = Cylinder_H/2.0 + Particles[i].z;
	if(Dgap < critic_radius)
	{
		move_flag=1;
		Particles[i].z += critic_radius - Dgap;
	}

	// Case 2: Top
	Dgap= Cylinder_H/2.0 - Particles[i].z;
	if(Dgap< critic_radius)
	{
		move_flag=1;
		Particles[i].z -= critic_radius - Dgap;
	}

	// Case 3: Lateral
	Dgap = Cylinder_D/2.0 
		- sqrt(Particles[i].x*Particles[i].x+Particles[i].y*Particles[i].y);
	if(Dgap< critic_radius)
	{
		move_flag=1;
		lat_r= sqrt(Particles[i].x*Particles[i].x + Particles[i].y*Particles[i].y);
		Particles[i].x -= (critic_radius - Dgap) * Particles[i].x/lat_r;
		Particles[i].y -= (critic_radius - Dgap) * Particles[i].y/lat_r;
	}

	}
	return move_flag;
}

// Find the nearest distance between the two particles in the unit box,
//   considering the periodic boundary conditions.
double NearestDist(int m, int n)
{
	double Dist0 = L, Dist1;

	if(n>=0)
	{
		Dist0=(Particles[m].x-Particles[n].x)*(Particles[m].x-Particles[n].x)
			+ (Particles[m].y-Particles[n].y)*(Particles[m].y-Particles[n].y)
			+ (Particles[m].z-Particles[n].z)*(Particles[m].z-Particles[n].z);

	Dist0 = pow(Dist0, 0.5);
	Dist0 = Dist0 /(Particles[m].r + Particles[n].r);

	Rod0.distance=Dist0;
	Rod0.P[0] = &Particles[m];
	Rod0.P[1] = &Particles[n];

	InsertRod();
	}

	return Dist0;
}

int MoveParticles(double dist0)
{
	int i=0, Ni, NP0, NP1;
	double delta, dist1;
//	double dx, dy, dz;
	double VectorX, VectorY, VectorZ;
	double Scaler0, Scaler1;
	
	extern double SortRods(int np);

	NP0=(*Rods[i].P[0]).number;
	NP1=(*Rods[i].P[1]).number;

	delta = ((*Rods[i].P[0]).r+(*Rods[i].P[1]).r)
		*(OuterRadius - Rods[i].distance);

	if(delta<=0) goto MPEND;

		dist1=((*Rods[i].P[0]).r+(*Rods[i].P[1]).r)*Rods[i].distance;

		// Vector of two points
		VectorX = ((*Rods[i].P[0]).x-(*Rods[i].P[1]).x)/dist1;
		VectorY = ((*Rods[i].P[0]).y-(*Rods[i].P[1]).y)/dist1;
		VectorZ = ((*Rods[i].P[0]).z-(*Rods[i].P[1]).z)/dist1;

		// Scaler for the movements of two particles
		Scaler0 = (*Rods[i].P[1]).r / ((*Rods[i].P[0]).r+(*Rods[i].P[1]).r);
		Scaler1 = (*Rods[i].P[0]).r / ((*Rods[i].P[0]).r+(*Rods[i].P[1]).r);

		Particles[NP0].x += delta * VectorX * Scaler0;
		Particles[NP0].y += delta * VectorY * Scaler0;
		Particles[NP0].z += delta * VectorZ * Scaler0;
		Particles[NP1].x -= delta * VectorX * Scaler1;
		Particles[NP1].y -= delta * VectorY * Scaler1;
		Particles[NP1].z -= delta * VectorZ * Scaler1;

		CheckBoundary(NP0);
		CheckBoundary(NP1);

		// Change the list of Rods[]
		SortRods(i);
		
MPEND: return 1;
}


void ReduceOuterRadius()
{
	double D_Ratio, powerj;
	if(upper<=lower) goto EOR;
	powerj = floor(-log10(upper-lower));
	D_Ratio = pow(0.5, powerj)*ContractionRate/N_Particles;
	lam -= D_Ratio;
	OuterRadius = lam * OuterRadius0;
EOR:	return;
}

