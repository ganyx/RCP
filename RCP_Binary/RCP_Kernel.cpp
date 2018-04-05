// RCP_Kernel.cpp

// Randomly distribute the intiatal packing

void Initiation()
{
	int i;
	extern double SortRods(int np);
	Volume_Factor = 0.0;

	fprintf(fp," The Number of Particles is: %d\n",N_Particles);

	Gm=0;
	for(Gi=-1;Gi<2;Gi++){
	for(Gj=-1;Gj<2;Gj++){
		for(Gk=-1;Gk<2;Gk++){
			SBox[Gm].number=Gm;
			SBox[Gm].x=Gk*L;
			SBox[Gm].y=Gj*L;
			SBox[Gm].z=Gi*L;
			Gm++;
	}}}
	
	// Binary system
	double R_Max = 1.0 + R_Relative_Width/2.0;
	double R_Min = 1.0 - R_Relative_Width/2.0;
	double alpha = R_Max/R_Min;

	double N_RMax = N_Particles
		/(1.0+pow(alpha,3)*(1.0/RelativeDensity_RMax-1.0));
	// End

	srand(time(NULL));
	for(i=0;i<N_Particles;i++)
	{
		Particles[i].number=i;
		Particles[i].x=L*rand()/RAND_MAX;
		Particles[i].y=L*rand()/RAND_MAX;
		Particles[i].z=L*rand()/RAND_MAX;
		
		// Distribution of particles' radius, offset/relative
		// Ri = r*R0
		if(i<N_RMax)
			Particles[i].r = R_Max;
		else Particles[i].r = R_Min;

		Volume_Factor += pow(Particles[i].r, 3);
	}

	// The OuterRadius in unequal sizes particles' system
	// 	means the mean value of curent radius distribution.
	OuterRadius = L * pow(3.0/4.0/Pi/Volume_Factor, 0.33333333);
	OuterRadius0 = OuterRadius;
	lam=1.0;
	SortRods(-1);
	return;
}

// Queue the most overlaps betweeb particles in the current system.
double SortRods(int np)
{
	int sri,srj,NP0,NP1;
	extern double NearestDist(int m, int n);
	extern void InsertRod();
	extern void DeleteRods(int nnp);

	if(np<0)
	{
		// Initial the Rods[]
		for(sri=0;sri<N_Rods;sri++)
		{
			Rods[sri].distance = L;
		}

		for(sri=0;sri<N_Particles;sri++)
		{
			for(srj=sri+1;srj<N_Particles;srj++)
			{
				NearestDist(sri, srj);
				InsertRod();
			}
		}
	}
	else
	{
		NP0=(*Rods[np].P[0]).number;
		NP1=(*Rods[np].P[1]).number;
		DeleteRods(NP0);
		DeleteRods(NP1);

		NearestDist(NP0, NP1);
		InsertRod();
		for(sri=0;sri<N_Particles;sri++)
		{
			if(!(sri==NP0 || sri==NP1))
			{
				NearestDist(sri, NP0);
				InsertRod();

				NearestDist(sri, NP1);
				InsertRod();
			}
			continue;
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
			Rods[j].N_Box=Rods[j-1].N_Box;
		}
			Rods[i+1].distance=Rod0.distance;
			Rods[i+1].P[0]=Rod0.P[0];
			Rods[i+1].P[1]=Rod0.P[1];
			Rods[i+1].N_Box=Rod0.N_Box;
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
				Rods[j].N_Box=Rods[j+1].N_Box;
			}
			Rods[N_Rods-1].distance = L;
		}
	}
	return;
}


// Find the nearest distance between the two particles in the unit box,
//   considering the periodic boundary conditions.
double NearestDist(int m, int n)
{
	int ndi;
	double Dist0 = L, Dist1;
	extern bool InKernelBox(int nnp);

/************************
 * Put checking on whether one of these two points is in kernel region, 
 * if yes, skip the PBC searching, in order to reduce the CPU time.
 ************************/

	if(InKernelBox(m)||InKernelBox(n))
	{
		Dist0=(Particles[m].x-Particles[n].x)*(Particles[m].x-Particles[n].x)
			+ (Particles[m].y-Particles[n].y)*(Particles[m].y-Particles[n].y)
			+ (Particles[m].z-Particles[n].z)*(Particles[m].z-Particles[n].z);
		Rod0.N_Box=13;
	}

	else
	{
		// Periodic Boundary Conditions
		for(ndi=0;ndi<27;ndi++)
		{
		Dist1=(Particles[m].x-Particles[n].x + SBox[ndi].x)
			*(Particles[m].x-Particles[n].x + SBox[ndi].x)
			+ (Particles[m].y-Particles[n].y + SBox[ndi].y)
			*(Particles[m].y-Particles[n].y + SBox[ndi].y)
			+ (Particles[m].z-Particles[n].z + SBox[ndi].z)
			*(Particles[m].z-Particles[n].z + SBox[ndi].z);

		if(Dist1 < Dist0)
		{
			Dist0=Dist1;
			Rod0.N_Box=ndi;
		}
		}
	}

/* r0 = distance(i,j) /(f(i)+f(j)) 
 * Searching for the minimum r0 (the worst overlapping) */
	Dist0 = pow(Dist0, 0.5);
	Dist0 = Dist0 /(Particles[m].r + Particles[n].r);

	Rod0.distance=Dist0;
	Rod0.P[0] = &Particles[m];
	Rod0.P[1] = &Particles[n];

	return Dist0;
}

int MoveParticles(double dist0)
{
	int i, Ni, NP0, NP1;
	double delta, dist1;
//	double dx, dy, dz;
	double VectorX, VectorY, VectorZ;
	double Scaler0, Scaler1;
	
	extern void CheckPBCs(struct Particle *P0);
	extern double SortRods(int np);

	Ni=1;

	for(i=0;i<Ni;i++)
	{
		delta = ((*Rods[i].P[0]).r+(*Rods[i].P[1]).r)
			*(OuterRadius - Rods[i].distance);

		if(delta<=0) goto MPEND;

		dist1=((*Rods[i].P[0]).r+(*Rods[i].P[1]).r)*Rods[0].distance;

		NP0=(*Rods[i].P[0]).number;
		NP1=(*Rods[i].P[1]).number;

		// Vector of two points
		VectorX = ((*Rods[i].P[0]).x-(*Rods[i].P[1]).x 
				+ SBox[Rods[i].N_Box].x)/dist1;
		VectorY = ((*Rods[i].P[0]).y-(*Rods[i].P[1]).y 
				+ SBox[Rods[i].N_Box].y)/dist1;
		VectorZ = ((*Rods[i].P[0]).z-(*Rods[i].P[1]).z 
				+ SBox[Rods[i].N_Box].z)/dist1;

		// Scaler for the movements of two particles
		Scaler0 = (*Rods[i].P[1]).r / ((*Rods[i].P[0]).r+(*Rods[i].P[1]).r);
		Scaler1 = (*Rods[i].P[0]).r / ((*Rods[i].P[0]).r+(*Rods[i].P[1]).r);

		Particles[NP0].x += delta * VectorX * Scaler0;
		Particles[NP0].y += delta * VectorY * Scaler0;
		Particles[NP0].z += delta * VectorZ * Scaler0;
		Particles[NP1].x -= delta * VectorX * Scaler1;
		Particles[NP1].y -= delta * VectorY * Scaler1;
		Particles[NP1].z -= delta * VectorZ * Scaler1;

		CheckPBCs(&Particles[NP0]);
		CheckPBCs(&Particles[NP1]);

		// Change the list of Rods[]
		SortRods(i);
		
MPEND: continue;
	}

	return 1;
}

void CheckPBCs(struct Particle *P0)
{	
	if((*P0).x >= L) (*P0).x -=L;
	if((*P0).x <  0) (*P0).x +=L;
	if((*P0).y >= L) (*P0).y -=L;
	if((*P0).y <  0) (*P0).y +=L;
	if((*P0).z >= L) (*P0).z -=L;
	if((*P0).z <  0) (*P0).z +=L;
	return;
}


void ReduceOuterRadius()
{
	double D_Ratio, powerj;
	powerj = floor(-log10(upper-lower));
	D_Ratio = pow(0.5, powerj)*ContractionRate/N_Particles;
	lam -= D_Ratio;
	OuterRadius = lam * OuterRadius0;
	return;
}

bool InKernelBox(int np)
{
	double Box_Max = L - (2.0 + R_Relative_Width)*OuterRadius;
	double Box_Min = (2.0 + R_Relative_Width)*OuterRadius;
	bool flag;

	if(Particles[np].x <= Box_Max && Particles[np].x >= Box_Min
		&& Particles[np].y <= Box_Max && Particles[np].y >= Box_Min
		&& Particles[np].z <= Box_Max && Particles[np].z >= Box_Min) flag = 1;
	else flag = 0;

	return flag;
}
