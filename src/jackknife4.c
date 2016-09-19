#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>


//open the file
FILE *open_file(const char *path,const char *mode)
{
	FILE *file=fopen(path,mode);
	if(file==NULL)
    {
		fprintf(stderr,"errore, file %s non aperto\n",path);
		exit(1);
    }
	return file;
}  

int jackknife (int file_created, double complex *error_diag,int ncopies,int block_size, int nconf,double complex Tr_S[][3], double complex *susc_tot, double complex** mat_restot)
{
	const double V_4= 32.0*32.0*4.0;
	int n_blocks;	
	const int nflavs=3;
	const int nind_trace=9;	
	double complex Tr_M_dM_M_dM[nflavs];
	double complex Tr_M_dM[nflavs];
	double complex Tr_M_d2M[nflavs];
	double complex Tr_M_dM2[nflavs];
	double complex Tr_M_dM3[nflavs];
	double complex Tr_M_dM_Tr_M_dM_M_dM[nflavs];
	double complex Tr_M_dM_Tr_M_d2M[nflavs];
	double complex Tr_M_dM_M_d2M[nflavs];
	double complex Tr_M_dM_M_dM_M_dM[nflavs];
        double complex Tr_M_d2M_M_d2M[nflavs];
	double complex Tr_M_dM_M_dM_M_d2M[nflavs];
	double complex Tr_M_dM_M_dM_M_dM_M_dM[nflavs];
	double complex Tr_M_dM_M_dM_Tr_M_dM_Tr_M_dM[nflavs];
	double complex Tr_M_dM_Tr_M_dM_Tr_M_dM_Tr_M_dM[nflavs];
	double complex Tr_M_dM_Tr_M_dM_Tr_M_d2M[nflavs];
	double complex Tr_M_dM_M_dM_Tr_M_dM_M_dM[nflavs];
	double complex Tr_M_dM_M_dM_M_dM_Tr_M_dM[nflavs];
	double complex Tr_M_d2M_Tr_M_d2M[nflavs];
        double complex Tr_M_d2M_Tr_M_dM_M_dM[nflavs]; 
        double complex Tr_M_dM_Tr_M_dM_M_d2M[nflavs];

	double  re_err_diag[nflavs];
	double  im_err_diag[nflavs];
	
	//inizializzo i vettori a 0
	for (int iflav=0; iflav<nflavs; iflav++)
	    {
		Tr_M_dM[iflav]=Tr_M_d2M[iflav]=Tr_M_dM_M_dM[iflav]=Tr_M_dM2[iflav]=Tr_M_dM_M_d2M[iflav]=0;
                Tr_M_d2M_M_d2M[iflav]= Tr_M_dM_M_dM_M_d2M[iflav]= Tr_M_dM_M_dM_M_dM_M_dM[iflav]= Tr_M_dM_M_dM_Tr_M_dM_Tr_M_dM[iflav]=        Tr_M_dM_Tr_M_dM_Tr_M_dM_Tr_M_dM[iflav]=Tr_M_d2M_Tr_M_dM_M_dM[iflav]=Tr_M_dM_Tr_M_dM_M_d2M[iflav]=0;
		Tr_M_dM_Tr_M_dM_Tr_M_d2M[iflav]= Tr_M_dM_M_dM_M_dM_Tr_M_dM[iflav]= Tr_M_d2M_Tr_M_d2M[iflav]=0;
		Tr_M_dM_M_dM_M_dM[iflav]=Tr_M_dM3[iflav]=Tr_M_dM_Tr_M_d2M[iflav]=Tr_M_dM_Tr_M_dM_M_dM[iflav]=0;
	         re_err_diag[iflav]=im_err_diag[iflav]= 0;	
			
			
		}

	
	//calcolo il numero di blocchi 	
	n_blocks= (int)(nconf/block_size);
	
	//loop sui blocchi
	
	for(int iblock=0;iblock<n_blocks;iblock++)
	{
		int nstart= iblock*block_size*ncopies;
		//reinizializzo i vettori a 0 (perchè calcolo le tracce su un blocco diverso)
		for (int iflav=0; iflav<nflavs; iflav++)
	       {
			Tr_M_dM[iflav]=Tr_M_d2M[iflav]=Tr_M_dM_M_dM[iflav]=Tr_M_dM2[iflav]=Tr_M_dM_M_d2M[iflav]=0;
                Tr_M_d2M_M_d2M[iflav]= Tr_M_dM_M_dM_M_d2M[iflav]= Tr_M_dM_M_dM_M_dM_M_dM[iflav]= Tr_M_dM_M_dM_Tr_M_dM_Tr_M_dM[iflav]=        Tr_M_dM_Tr_M_dM_Tr_M_dM_Tr_M_dM[iflav]=Tr_M_d2M_Tr_M_dM_M_dM[iflav]=Tr_M_dM_Tr_M_dM_M_d2M[iflav]=0;
		Tr_M_dM_Tr_M_dM_Tr_M_d2M[iflav]= Tr_M_dM_M_dM_M_dM_Tr_M_dM[iflav]= Tr_M_d2M_Tr_M_d2M[iflav]=0;
		Tr_M_dM_M_dM_M_dM[iflav]=Tr_M_dM3[iflav]=Tr_M_dM_Tr_M_d2M[iflav]=Tr_M_dM_Tr_M_dM_M_dM[iflav]=0;
		    }			
		for(int i=0; i< block_size; i++)
		{ 
			//reinit to zero temp traces
			
			
			
			for(int iflav=0;iflav<nflavs;iflav++)
			{	
				//compute all mean values
				
				//compute Tr_M_dM
				for(int icopy=0;icopy<ncopies;icopy++) 
				{
					Tr_M_dM[iflav]+=     mat_restot[nstart +i*ncopies+icopy][iflav*nind_trace+0];
					
					
				 	
				}	
				//compute Tr_M_d2M
				for(int icopy=0;icopy<ncopies;icopy++) 
				{
					Tr_M_d2M[iflav]+=    mat_restot[nstart+i*ncopies+icopy][iflav*nind_trace+1];
					
				
				}
				
				
				//compute Tr_M_dM_M_dM
				for(int icopy=0;icopy<ncopies;icopy++) 
				{
					Tr_M_dM_M_dM[iflav]+=mat_restot[nstart+i*ncopies+icopy][iflav*nind_trace+2];
					
						
				}	
				//compute Tr_M_dM_M_d2M
				for(int icopy=0;icopy<ncopies;icopy++) Tr_M_dM_M_d2M[iflav]+=mat_restot[nstart+i*ncopies+icopy][iflav*nind_trace+4];	
				
				//compute Tr_M_dM_M_dM_M_dM
				for(int icopy=0;icopy<ncopies;icopy++) Tr_M_dM_M_dM_M_dM[iflav]+=mat_restot[nstart+i*ncopies+icopy][iflav*nind_trace+5];				
				
				//compute (TrM_dM)^2
				for(int icopy=0;icopy<ncopies;icopy++)
					for(int rcopy=icopy+1;rcopy<ncopies;rcopy++)
					{
						double complex complex1=mat_restot[nstart+i*ncopies+icopy][iflav*nind_trace+0];
						double complex complex2=mat_restot[nstart+i*ncopies+rcopy][iflav*nind_trace+0];
						Tr_M_dM2[iflav]+=complex1*complex2;
						
												
					}
				
				//compute (TrM_dM)(Tr(M_dM)^2)
				for(int icopy=0;icopy<ncopies;icopy++)
					for(int rcopy=icopy+1;rcopy<ncopies;rcopy++)
					{
						double complex complex1=mat_restot[nstart+i*ncopies+icopy][iflav*nind_trace+0];
						double complex complex2=mat_restot[nstart+i*ncopies+rcopy][iflav*nind_trace+2];
						Tr_M_dM_Tr_M_dM_M_dM[iflav]+=complex1*complex2;
						//compute the products for rcopy<icopy
						complex1=mat_restot[nstart+i*ncopies+rcopy][iflav*nind_trace+0];
						complex2=mat_restot[nstart+i*ncopies+icopy][iflav*nind_trace+2];
						Tr_M_dM_Tr_M_dM_M_dM[iflav]+=complex1*complex2;					
					}
				//compute (TrM_dM)(TrM_d2M)
				for(int icopy=0;icopy<ncopies;icopy++)
					for(int rcopy=icopy+1;rcopy<ncopies;rcopy++)
					{
						double complex complex1=mat_restot[nstart+i*ncopies+icopy][iflav*nind_trace+0];
						double complex complex2=mat_restot[nstart+i*ncopies+rcopy][iflav*nind_trace+1];
						Tr_M_dM_Tr_M_d2M[iflav]+=complex1*complex2;
						//compute the products for rcopy<ncopy
						complex1=mat_restot[nstart+i*ncopies+rcopy][iflav*nind_trace+0];
						complex2=mat_restot[nstart+i*ncopies+icopy][iflav*nind_trace+1];
						Tr_M_dM_Tr_M_d2M[iflav]+=complex1*complex2;						
						
					}
				
				//compute (TrM_dM)^3
				for(int icopy=0;icopy<ncopies;icopy++)
				{
					for(int rcopy=0;rcopy<ncopies;rcopy++)
					{  	
						if(rcopy==icopy);
				        else
						{	
							for(int scopy=0;scopy<ncopies;scopy++)
							{
								if(scopy==rcopy || scopy==icopy) ;
								else
								{		
									double complex complex1=mat_restot[nstart+i*ncopies+icopy][iflav*nind_trace+0];
									double complex complex2=mat_restot[nstart+i*ncopies+rcopy][iflav*nind_trace+0];
									double complex complex3=mat_restot[nstart+i*ncopies+scopy][iflav*nind_trace+0];
									Tr_M_dM3[iflav]+=complex1*complex2*complex3;
						
					           }
				             }
						   }	   
					}
				}
			


                //parte di calcolo delle suscettività del quart'ordine

				//compute Tr_M_d2M_M_d2M

				for (int icopy=0;icopy<ncopies;icopy++)
				 {
				  Tr_M_d2M_M_d2M[iflav] += mat_restot[nstart+i*ncopies+icopy][iflav*nind_trace + 6];
				 }

				//compute Tr_M_d2M_M_dM_M_dM

				for (int icopy=0;icopy<ncopies;icopy++)
				 {
				  Tr_M_dM_M_dM_M_d2M[iflav] += mat_restot[nstart+i*ncopies+icopy][iflav*nind_trace + 7];
				 }

				 //compute Tr_M_dM_M_dM_M_dM_M_dM

				for (int icopy=0;icopy<ncopies;icopy++)
				 {
				  Tr_M_dM_M_dM_M_dM_M_dM[iflav] += mat_restot[nstart+i*ncopies+icopy][iflav*nind_trace + 8];
				 }

				 //compute Tr_M_d2M_Tr_M_dM_M_dM

				 for(int icopy=0;icopy<ncopies;icopy++)
					for(int rcopy=icopy+1;rcopy<ncopies;rcopy++)
					{
						double complex complex1=mat_restot[nstart+i*ncopies+icopy][iflav*nind_trace+1];
						double complex complex2=mat_restot[nstart+i*ncopies+rcopy][iflav*nind_trace+2];
						Tr_M_d2M_Tr_M_dM_M_dM[iflav]+=complex1*complex2;
						               //compute the products for rcopy<ncopy
						               complex1=mat_restot[nstart+i*ncopies+rcopy][iflav*nind_trace+1];
						               complex2=mat_restot[nstart+i*ncopies+icopy][iflav*nind_trace+2];
						Tr_M_d2M_Tr_M_dM_M_dM[iflav]+=complex1*complex2;

					}
				 //compute Tr_M_dM_M_dM_Tr_M_dM_M_dM

				 for(int icopy=0;icopy<ncopies;icopy++)
					for(int rcopy=icopy+1;rcopy<ncopies;rcopy++)
					{
						double complex complex1=mat_restot[nstart+i*ncopies+icopy][iflav*nind_trace+2];
						double complex complex2=mat_restot[nstart+i*ncopies+rcopy][iflav*nind_trace+2];
						Tr_M_dM_M_dM_Tr_M_dM_M_dM[iflav]+=complex1*complex2;
					}

				 //compute Tr_M_dM_M_dM_M_dM_Tr_M_dM

				 for(int icopy=0;icopy<ncopies;icopy++)
					for(int rcopy=icopy+1;rcopy<ncopies;rcopy++)
					{
						double complex complex1=mat_restot[nstart+i*ncopies+icopy][iflav*nind_trace+5];
						double complex complex2=mat_restot[nstart+i*ncopies+rcopy][iflav*nind_trace+0];
						Tr_M_dM_M_dM_M_dM_Tr_M_dM[iflav]+=complex1*complex2;
						               //compute the products for rcopy<ncopy
						               complex1=mat_restot[nstart+i*ncopies+rcopy][iflav*nind_trace+5];
						               complex2=mat_restot[nstart+i*ncopies+icopy][iflav*nind_trace+0];
						Tr_M_dM_M_dM_M_dM_Tr_M_dM[iflav]+=complex1*complex2;

					}


				  //compute Tr_M_d2M_Tr_M_d2M

				  for(int icopy=0;icopy<ncopies;icopy++)
					for(int rcopy=icopy+1;rcopy<ncopies;rcopy++)
					{
						double complex complex1=mat_restot[nstart+i*ncopies+icopy][iflav*nind_trace+1];
						double complex complex2=mat_restot[nstart+i*ncopies+rcopy][iflav*nind_trace+1];
						Tr_M_d2M_Tr_M_d2M[iflav]+=complex1*complex2;

					}


				 //compute Tr_M_dM_Tr_M_dM_M_d2M

				   for(int icopy=0;icopy<ncopies;icopy++)
					for(int rcopy=icopy+1;rcopy<ncopies;rcopy++)
					{
						double complex complex1=mat_restot[nstart+i*ncopies+icopy][iflav*nind_trace+0];
						double complex complex2=mat_restot[nstart+i*ncopies+rcopy][iflav*nind_trace+4];
						Tr_M_dM_Tr_M_dM_M_d2M[iflav]+=complex1*complex2;
						               //compute the products for rcopy<ncopy
						               complex1=mat_restot[nstart+i*ncopies+rcopy][iflav*nind_trace+0];
						               complex2=mat_restot[nstart+i*ncopies+icopy][iflav*nind_trace+4];
						Tr_M_dM_Tr_M_dM_M_d2M[iflav]+=complex1*complex2;

					}


				  //compute Tr_M_dM_M_dM_Tr_M_dM_Tr_M_dM
				  for(int icopy=0;icopy<ncopies;icopy++)
				{
					for(int rcopy=0;rcopy<ncopies;rcopy++)
					{
						if(rcopy==icopy);
				        else
						   {
						     for(int scopy=0;scopy<ncopies;scopy++)
						     {
							   if(scopy==rcopy || scopy==icopy) ;
							   else
					           {
						       double complex complex1=mat_restot[nstart+i*ncopies+icopy][iflav*nind_trace+2];
						       double complex complex2=mat_restot[nstart+i*ncopies+rcopy][iflav*nind_trace+0];
						       double complex complex3=mat_restot[nstart+i*ncopies+scopy][iflav*nind_trace+0];
						       Tr_M_dM_M_dM_Tr_M_dM_Tr_M_dM[iflav]+=complex1*complex2*complex3;

					           }
				             }
						   }
					}
				}

				//compute Tr_M_d2M_Tr_M_dM_Tr_M_dM

				for(int icopy=0;icopy<ncopies;icopy++)
				{
					for(int rcopy=0;rcopy<ncopies;rcopy++)
					{
						if(rcopy==icopy);
				        else
						   {
						     for(int scopy=0;scopy<ncopies;scopy++)
						     {
							   if(scopy==rcopy || scopy==icopy) ;
							   else
					           {
						       double complex complex1=mat_restot[nstart+i*ncopies+icopy][iflav*nind_trace+1];
						       double complex complex2=mat_restot[nstart+i*ncopies+rcopy][iflav*nind_trace+0];
						       double complex complex3=mat_restot[nstart+i*ncopies+scopy][iflav*nind_trace+0];
						       Tr_M_dM_Tr_M_dM_Tr_M_d2M[iflav]+=complex1*complex2*complex3;

					           }
				             }
						   }
					}
				}

				//compute Tr_M_dM_Tr_M_dM_Tr_M_dM_Tr_M_dM

				for(int icopy=0;icopy<ncopies;icopy++)
				{
					for(int rcopy=icopy+1;rcopy<ncopies;rcopy++)
					{
						
						     for(int scopy=rcopy+1;scopy<ncopies;scopy++)
						     {
							
					             for(int tcopy=scopy+1;tcopy<ncopies;tcopy++)
					             {
					               
						       double complex complex1=mat_restot[nstart+i*ncopies+icopy][iflav*nind_trace+0];
						       double complex complex2=mat_restot[nstart+i*ncopies+rcopy][iflav*nind_trace+0];
						       double complex complex3=mat_restot[nstart+i*ncopies+scopy][iflav*nind_trace+0];
						       double complex complex4=mat_restot[nstart+i*ncopies+tcopy][iflav*nind_trace+0];
						       Tr_M_dM_Tr_M_dM_Tr_M_dM_Tr_M_dM[iflav]+=complex1*complex2*complex3*complex4;
                                  
				             }
						   }
					}
				}			

                     
                       //esco dal loop sul flavour
			}
			
				

		//esco dal loop sulle configurazioni nel blocco i_block	
		}
		
		
		
		for(int iflav=0;iflav<nflavs;iflav++)
		{	
		//calcolo la somma delle tracce sul sottocampione

		Tr_M_dM[iflav] = Tr_S[0][iflav] - Tr_M_dM[iflav];
		Tr_M_d2M[iflav]= Tr_S[1][iflav] - Tr_M_d2M[iflav]; 	
		Tr_M_dM_M_dM[iflav]= Tr_S[2][iflav] - Tr_M_dM_M_dM[iflav];
		Tr_M_dM_M_d2M[iflav]= Tr_S[7][iflav] - Tr_M_dM_M_d2M[iflav];
		Tr_M_dM_M_dM_M_dM[iflav]= Tr_S[8][iflav] - Tr_M_dM_M_dM_M_dM[iflav];
		Tr_M_dM2[iflav]= Tr_S[3][iflav] - Tr_M_dM2[iflav];
		Tr_M_dM_Tr_M_dM_M_dM[iflav]= Tr_S[5][iflav]- Tr_M_dM_Tr_M_dM_M_dM[iflav];
		Tr_M_dM_Tr_M_d2M[iflav]= Tr_S[6][iflav]- Tr_M_dM_Tr_M_d2M[iflav];
		Tr_M_dM3[iflav]=Tr_S[4][iflav]-Tr_M_dM3[iflav];
                Tr_M_d2M_M_d2M[iflav] =   Tr_S[9][iflav]- Tr_M_d2M_M_d2M[iflav];
                Tr_M_dM_M_dM_M_d2M[iflav] =    Tr_S[10][iflav]- Tr_M_dM_M_dM_M_d2M[iflav];
                Tr_M_dM_M_dM_M_dM_M_dM[iflav] =    Tr_S[11][iflav]- Tr_M_dM_M_dM_M_dM_M_dM[iflav];
                Tr_M_dM_M_dM_Tr_M_dM_Tr_M_dM[iflav] =    Tr_S[12][iflav]- Tr_M_dM_M_dM_Tr_M_dM_Tr_M_dM[iflav];
                Tr_M_dM_Tr_M_dM_Tr_M_dM_Tr_M_dM[iflav] =    Tr_S[13][iflav]- Tr_M_dM_Tr_M_dM_Tr_M_dM_Tr_M_dM[iflav];
                Tr_M_dM_Tr_M_dM_Tr_M_d2M[iflav] =    Tr_S[14][iflav]- Tr_M_dM_Tr_M_dM_Tr_M_d2M[iflav];
                Tr_M_dM_M_dM_Tr_M_dM_M_dM[iflav] =    Tr_S[15][iflav]- Tr_M_dM_M_dM_Tr_M_dM_M_dM[iflav];
                Tr_M_dM_M_dM_M_dM_Tr_M_dM[iflav] =    Tr_S[16][iflav]- Tr_M_dM_M_dM_M_dM_Tr_M_dM[iflav];
                Tr_M_d2M_Tr_M_d2M[iflav] =    Tr_S[17][iflav]- Tr_M_d2M_Tr_M_d2M[iflav];
                Tr_M_d2M_Tr_M_dM_M_dM[iflav] = Tr_S[18][iflav] - Tr_M_d2M_Tr_M_dM_M_dM[iflav];
                Tr_M_dM_Tr_M_dM_M_d2M[iflav] = Tr_S[19][iflav] - Tr_M_dM_Tr_M_dM_M_d2M[iflav];
		}
		
		//calcolo la somma del prodotto di tracce sul sottocampione
		
		
		//normalizzo le tracce
		for(int iflav=0;iflav<nflavs;iflav++)
		{	
			Tr_M_dM2[iflav] = (2.0*Tr_M_dM2[iflav])/((nconf - block_size)*(ncopies)*(ncopies -1));
			Tr_M_dM[iflav]= ((1.0)*Tr_M_dM[iflav])/((nconf - block_size)*(ncopies));
			Tr_M_dM_M_dM[iflav]= ((1.0)*Tr_M_dM_M_dM[iflav])/((nconf - block_size)*(ncopies));
			Tr_M_d2M[iflav]= ((1.0)*Tr_M_d2M[iflav])/((nconf - block_size)*(ncopies));
			Tr_M_dM_M_d2M[iflav]= ((1.0)*Tr_M_dM_M_d2M[iflav])/((nconf - block_size)*(ncopies));	
			Tr_M_dM_M_dM_M_dM[iflav]= ((1.0)*Tr_M_dM_M_dM_M_dM[iflav])/((nconf - block_size)*(ncopies));
			Tr_M_dM_Tr_M_dM_M_dM[iflav] = ((1.0)*Tr_M_dM_Tr_M_dM_M_dM[iflav])/((nconf - block_size)*(ncopies)*(ncopies -1));
			Tr_M_dM_Tr_M_d2M[iflav] = ((1.0)*Tr_M_dM_Tr_M_d2M[iflav])/((nconf - block_size)*(ncopies)*(ncopies -1));	
			Tr_M_dM3[iflav] = ((1.0)*Tr_M_dM3[iflav])/((nconf - block_size)*(ncopies)*(ncopies -1)*(ncopies -2));	
                        Tr_M_d2M_M_d2M[iflav]*= (1.0)/((nconf - block_size)*(ncopies));
                        Tr_M_dM_M_dM_M_d2M[iflav]*= (1.0)/((nconf - block_size)*(ncopies));
                        Tr_M_dM_M_dM_M_dM_M_dM[iflav]*= (1.0)/((nconf - block_size)*(ncopies));
                        Tr_M_dM_M_dM_Tr_M_dM_Tr_M_dM[iflav]*=(1.0)/((nconf - block_size)*(ncopies)*(ncopies-1)*(ncopies-2));
                        Tr_M_dM_Tr_M_dM_Tr_M_dM_Tr_M_dM[iflav]*=(24.0)/((nconf - block_size)*(ncopies)*(ncopies-1)*(ncopies-2)*(ncopies-3));
                        Tr_M_dM_Tr_M_dM_Tr_M_d2M[iflav]*= (1.0)/((nconf - block_size)*(ncopies)*(ncopies-1)*(ncopies-2));
                        Tr_M_dM_M_dM_Tr_M_dM_M_dM[iflav]*=(2.0)/((nconf - block_size)*(ncopies)*(ncopies-1));
                        Tr_M_dM_M_dM_M_dM_Tr_M_dM[iflav]*= (1.0)/((nconf - block_size)*(ncopies)*(ncopies-1));
                        Tr_M_d2M_Tr_M_d2M[iflav]*=(2.0)/((nconf - block_size)*(ncopies)*(ncopies-1));	
                        Tr_M_d2M_Tr_M_dM_M_dM[iflav] *= (1.0)/((nconf-block_size)*(ncopies)*(ncopies-1));
                        Tr_M_dM_Tr_M_dM_M_d2M[iflav] *= (1.0)/((nconf-block_size)*(ncopies)*(ncopies-1));
		
		
		
			
			
		}
		
		
		
		//calcolo le suscettivita sul sottocampione
		//4_0_0
		for(int iflav=0;iflav<nflavs;iflav++)
	    {
			
			
			double complex susc_block_diag= -((6.0)/(64))*Tr_M_dM_M_dM_Tr_M_dM_Tr_M_dM[iflav]
+((1.0)/(256))*Tr_M_dM_Tr_M_dM_Tr_M_dM_Tr_M_dM[iflav]
+((6.0)/(64))*Tr_M_dM_Tr_M_dM_Tr_M_d2M[iflav]
-((3.0)/(64))*Tr_M_dM2[iflav]*Tr_M_d2M[iflav]
+((3.0)/(64))*Tr_M_dM2[iflav]*Tr_M_dM_M_dM[iflav]
-((3.0)/(256))*Tr_M_dM2[iflav]*Tr_M_dM2[iflav]
-((6.0)/(16))*Tr_M_d2M_Tr_M_dM_M_dM[iflav]
+((3.0)/(16))*Tr_M_dM_M_dM_Tr_M_dM_M_dM[iflav]
+((1.0)/(2))*Tr_M_dM_M_dM_M_dM_Tr_M_dM[iflav]
+((5.0)/(16))*Tr_M_dM2[iflav]                
+((3.0)/(16))*Tr_M_d2M_Tr_M_d2M[iflav]
-((9.0)/(64))*Tr_M_dM_Tr_M_dM_M_d2M[iflav]
-((3.0)/(16))*Tr_M_d2M[iflav]*Tr_M_d2M[iflav]
+((3.0)/(16))*Tr_M_d2M[iflav]*Tr_M_dM_M_dM[iflav]
-((3.0)/(64))*Tr_M_d2M[iflav]*Tr_M_dM2[iflav]
-((3.0)/(4))*Tr_M_d2M_M_d2M[iflav]
+((3.0)/(1))*Tr_M_dM_M_dM_M_d2M[iflav]
+((1.0)/(1))*Tr_M_dM_M_dM[iflav]                        
-((3.0)/(2))*Tr_M_dM_M_dM_M_dM_M_dM[iflav]
-((1.0)/(4))*Tr_M_d2M[iflav]                            
+((3.0)/(16))*Tr_M_dM_M_dM[iflav]*Tr_M_d2M[iflav]                    
-((3.0)/(16))*Tr_M_dM_M_dM[iflav]*Tr_M_dM_M_dM[iflav]
+((3.0)/(64))*Tr_M_dM_M_dM[iflav]*Tr_M_dM2[iflav];
			
		susc_block_diag *= (1.0)/V_4;	
			
		//calcolo la differenza tra la media totale e la media sul sottocampione
		
     	        re_err_diag[iflav] += (creal(susc_block_diag) - creal(susc_tot[iflav]))*(creal(susc_block_diag) - creal(susc_tot[iflav]));
		im_err_diag[iflav] += (cimag(susc_block_diag) - cimag(susc_tot[iflav]))*(cimag(susc_block_diag) - cimag(susc_tot[iflav]));
		
		}	
													  
		
		
	//esco dal loop che scorre i vari blocchi	
	}
	
	// creo un file per inserire i risultati e li stampo
	FILE *file_out=open_file("jack_4_ord.txt", file_created?"a":"w");
	for(int iflav=0;iflav<nflavs;iflav++)
	{	
	  re_err_diag[iflav]= sqrt(re_err_diag[iflav]);
	  im_err_diag[iflav]= sqrt(im_err_diag[iflav]);
			
	}

	fprintf(file_out,"block_size:%d \t",block_size);
	for(int iflav=0;iflav<nflavs;iflav++)
    {
	fprintf(file_out, "%lf %lf \t",re_err_diag[iflav],im_err_diag[iflav]); 	
        if(block_size> 3) error_diag[iflav] += re_err_diag[iflav] + I*im_err_diag[iflav];
	}
	
																														  
	  fprintf(file_out,"\n");
          

 
      
																														  
	fclose(file_out);

	
	
	
	return 0;
}	









