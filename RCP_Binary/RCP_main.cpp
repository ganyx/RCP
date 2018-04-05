// RCP_main.cpp
// 3D Random Close Packing
 
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

	double InnerRadius0 = InnerRadius;
	SortRods(-1);
	printf("*** %10.8f, %10.8f\n", InnerRadius0, InnerRadius);
	fprintf(fp0, "%5d, %8d, %10.8f, %10.8f, %10.8f, %10.8f\n",
		ncase+1, mj, upper, lower, InnerRadius0, InnerRadius);

	FileOutput(&ncase);
	}

	fclose(fp0);
	fclose(fp);
	return 0;
}
