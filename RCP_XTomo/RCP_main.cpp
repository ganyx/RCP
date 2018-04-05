// RCP_main.cpp
// 3D Random Close Packing

/* Objective of this code: (differ from original code with periodic boundary conditions)
     Packing mono-sized spherical particles into a cylinder container (D x H)
	 Compare the microstructure with the data from X-ray tomography
			Coded by Yixiang Gan, 01.2009
*/
 
#include "IncludeFiles.h"
#include "InputParameters.h"
#include "DataStructure.cpp"

#include "RCP_Kernel.cpp" 
#include "Output.cpp"

int main()  
{   
	int mi,mj,ncase,srn=1;
	double Initial_Distance; 

	if((fp=fopen("Screen.txt","wb+"))==NULL)
	{
		return 0;
	}
	if((fp0=fopen("Data_Table.txt","wb+"))==NULL)
	{
		return 0;
	}
	
	for(ncase=0;ncase<N_Case;ncase++)
	{
	Initiation(); 

	for(mj=1; mj<=NMAX_Iteration; mj++) 
	{
		Initial_Distance=Rods[0].distance;
		for(mi=1;mi<=NMAX_Movement;mi++)
		{ 
		MoveParticles(Initial_Distance); 
		}
    
		//if(mj%10==0) 
			OutputData(&mj, &ncase);
		ReduceOuterRadius();
		if((OuterRadius-InnerRadius)/OuterRadius<=Tolerance) break;

// Refresh the queue of Rods
		if((mj/Fresh_Freq) >= srn)
		{
			SortRods(-1);
			srn++;
		}
	}

//	double InnerRadius0 = InnerRadius;
//	SortRods(-1);
//	printf("*** %10.8f, %10.8f\n", InnerRadius0, InnerRadius);
//	fprintf(fp0, "%5d, %8d, %10.8f, %10.8f, %10.8f, %10.8f\n",
//		ncase+1, mj, upper, lower, InnerRadius0, InnerRadius);

	FileOutput(&ncase);
	Density_Roof=0.0;
	}

	fclose(fp0);
	fclose(fp);
	return 0;
}
