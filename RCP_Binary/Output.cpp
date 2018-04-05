// Output.cpp

void OutputData(int *i, int *j)
{
	upper = Volume_Factor*4.0/3.0*Pi* pow(OuterRadius,3.0)/pow(L,3);
	lower = Volume_Factor*4.0/3.0*Pi* pow(InnerRadius,3.0)/pow(L,3);
	printf("%4d, %8d, %10.8f, %10.8f, %5d-%5d\n", 
		*j+1, *i, upper, lower, (*Rods[0].P[0]).number, (*Rods[0].P[1]).number);
//	fprintf(fp, "%4d, %8d, %10.8f, %10.8f \n", *j+1, *i, upper, lower);
	return;
}

void FileOutput(int *i)
{
	FILE *fout;
	const char base[] = "N";
	char filename[30];
	int j;
	sprintf(filename, "%s%d_%d.dat", base, N_Particles, *i);

	if((fout=fopen(filename,"wb+"))==NULL)
	{
//		cout<<"!Error!"<<endl;
//		cout<<"Cann't Open File: P_N.dat"<<endl<<endl;
		return;
	}
	fprintf(fout,"NP: %d CR: %10.8f R: %10.8f Density: %10.8f\n", 
		N_Particles, ContractionRate, InnerRadius, lower);

	for(j=0; j<N_Particles; j++)
	{
		fprintf(fout, "%10.8f %10.8f %10.8f %10.8f\n", 
			Particles[j].x, Particles[j].y, Particles[j].z, 
			Particles[j].r * InnerRadius);
	}

	fclose(fout);
	return;
}