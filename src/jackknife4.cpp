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

complex traces_calculator_jack(int num_prod, int *traces_id, double complex** matrix, int iflav,int ncopies,int nstart, int conf_id, int nind_trace)
{

double complex result=0;

double complex X[ncopies][num_prod-1];
//set X1[ncopies_eff-1][dim] to 0
for (int i=0; i<num_prod-1;i++)
 for (int icopy=0;icopy<ncopies;icopy++)   X[icopy][i]=0;

for (int j=0;j<num_prod-1;j++)
{  
   if(j==0)
   for(int icopy=ncopies-2;icopy>=0;icopy--)
   X[icopy][j] = X[icopy+1][j] + matrix[nstart + conf_id*ncopies +icopy+1][iflav*nind_trace + traces_id[num_prod-1-j]];

   else
   for (int icopy=ncopies-2-j;icopy>=0;icopy--)
    X[icopy][j] = X[icopy+1][j] + matrix[nstart + conf_id*ncopies+icopy+1][iflav*nind_trace+traces_id[num_prod-1-j]]*X[icopy+1][j-1];

}

for(int icopy=0;icopy<ncopies;icopy++) result += X[icopy][num_prod-2]*matrix[nstart + conf_id*ncopies+icopy][iflav*nind_trace + traces_id[0]];


return result;



}














int jackknife (int file_created, double complex *error_diag,int ncopies,int block_size, int nconf,double complex Tr_S[][3], double complex *susc_tot, double complex** mat_restot)
{
        const double V_4= 32.0*32.0*4.0;
	const int nflavs=3;
	const int nind_trace=9;	
        int traces_id3[3];
        int traces_id4[4];
        int traces_id2[2];
	int n_blocks;	
	 double complex Tr_M_dM[nflavs];
  double complex Tr_M_d2M[nflavs];
  double complex Tr_M_dM_M_dM[nflavs];
  double complex Tr_M_dM2[nflavs];
  double complex Tr_M_dM3[nflavs];
  double complex Tr_M_dM_Tr_M_dM_M_dM[nflavs];
  double complex Tr_M_dM_Tr_M_d2M[nflavs];
  double complex Tr_M_dM_M_d2M[nflavs];
  double complex Tr_M_dM_M_dM_M_dM[nflavs];
  double complex Tr_M_d2M_M_d2M[nflavs];
  double complex Tr_M_dM_M_dM_M_d2M[nflavs];
  double complex Tr_M_dM_M_dM_M_dM_M_dM[nflavs];
  double complex Tr_M_d2M_Tr_M_dM_M_dM[nflavs];
  double complex Tr_M_dM_Tr_M_dM_M_d2M[nflavs];
  double complex Tr_M_dM_M_dM_Tr_M_dM_Tr_M_dM[nflavs];
  double complex Tr_M_dM_Tr_M_dM_Tr_M_dM_Tr_M_dM[nflavs];
  double complex Tr_M_dM_Tr_M_dM_Tr_M_d2M[nflavs];
  double complex Tr_M_dM_M_dM_Tr_M_dM_M_dM[nflavs];
  double complex Tr_M_dM_M_dM_M_dM_Tr_M_dM[nflavs];
  double complex Tr_M_d2M_Tr_M_d2M[nflavs];

	double  re_err_diag[nflavs];
	double  im_err_diag[nflavs];
	
	//inizializzo i vettori a 0
	for (int iflav=0; iflav<nflavs; iflav++)
	    {
		Tr_M_dM[iflav]=Tr_M_d2M[iflav]=Tr_M_dM_M_dM[iflav]=Tr_M_dM2[iflav]=Tr_M_dM_M_d2M[iflav]=0;
                Tr_M_d2M_M_d2M[iflav]= Tr_M_dM_M_dM_M_d2M[iflav]= Tr_M_dM_M_dM_M_dM_M_dM[iflav]= Tr_M_dM_M_dM_Tr_M_dM_Tr_M_dM[iflav]=        Tr_M_dM_Tr_M_dM_Tr_M_dM_Tr_M_dM[iflav]=Tr_M_d2M_Tr_M_dM_M_dM[iflav]=Tr_M_dM_Tr_M_dM_M_d2M[iflav]=0;
		Tr_M_dM_Tr_M_dM_Tr_M_d2M[iflav]= Tr_M_dM_M_dM_M_dM_Tr_M_dM[iflav]= Tr_M_d2M_Tr_M_d2M[iflav]=0;
		Tr_M_dM_M_dM_M_dM[iflav]=Tr_M_dM3[iflav]=Tr_M_dM_Tr_M_d2M[iflav]=Tr_M_dM_Tr_M_dM_M_dM[iflav]=Tr_M_dM_M_dM_Tr_M_dM_M_dM[iflav]=0;
	         re_err_diag[iflav]=im_err_diag[iflav]= 0;	
			
			
		}

	
	//calcolo il numero di blocchi 	
	n_blocks= (int)(nconf/block_size);
        printf("NBLOCK: %d \n",n_blocks);
	
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
		Tr_M_dM_M_dM_M_dM[iflav]=Tr_M_dM3[iflav]=Tr_M_dM_Tr_M_d2M[iflav]=Tr_M_dM_Tr_M_dM_M_dM[iflav]=Tr_M_dM_M_dM_Tr_M_dM_M_dM[iflav]=0;
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

                                    traces_id2[0]=0;
                                    traces_id2[1]=0;
                                    Tr_M_dM2[iflav]+= traces_calculator_jack(2,traces_id2, mat_restot, iflav,ncopies,nstart,i,nind_trace);
				
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
				traces_id3[0]=0; traces_id3[1]=0; traces_id3[2]=0;
                  Tr_M_dM3[iflav]+= traces_calculator_jack(3,traces_id3, mat_restot, iflav,ncopies,nstart,i,nind_trace);


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
			traces_id3[0]=2; traces_id3[1]=0; traces_id3[2]=0;
                  Tr_M_dM_M_dM_Tr_M_dM_Tr_M_dM[iflav]+= traces_calculator_jack(3,traces_id3, mat_restot, iflav,ncopies,nstart,i, nind_trace);
                  traces_id3[0]=0; traces_id3[1]=2; traces_id3[2]=0;
                  Tr_M_dM_M_dM_Tr_M_dM_Tr_M_dM[iflav]+= traces_calculator_jack(3,traces_id3, mat_restot, iflav,ncopies,nstart,i, nind_trace);
                  traces_id3[0]=0; traces_id3[1]=0; traces_id3[2]=2;
                  Tr_M_dM_M_dM_Tr_M_dM_Tr_M_dM[iflav]+= traces_calculator_jack(3,traces_id3, mat_restot, iflav,ncopies,nstart,i, nind_trace);

				//compute Tr_M_d2M_Tr_M_dM_Tr_M_dM

				traces_id3[0]=0; traces_id3[1]=0; traces_id3[2]=1;
                  
                  Tr_M_dM_Tr_M_dM_Tr_M_d2M[iflav]+= traces_calculator_jack(3,traces_id3, mat_restot, iflav,ncopies,nstart,i, nind_trace);
                  traces_id3[0]=0; traces_id3[1]=1; traces_id3[2]=0;
                  
                  Tr_M_dM_Tr_M_dM_Tr_M_d2M[iflav]+= traces_calculator_jack(3,traces_id3, mat_restot, iflav,ncopies,nstart,i, nind_trace);
                  traces_id3[0]=1; traces_id3[1]=0; traces_id3[2]=0;
                  
                  Tr_M_dM_Tr_M_dM_Tr_M_d2M[iflav]+= traces_calculator_jack(3,traces_id3, mat_restot, iflav,ncopies,nstart,i, nind_trace);

				//compute Tr_M_dM_Tr_M_dM_Tr_M_dM_Tr_M_dM

				traces_id4[0]=0; traces_id4[1]=0; traces_id4[2]=0; traces_id4[3]=0;
                  Tr_M_dM_Tr_M_dM_Tr_M_dM_Tr_M_dM[iflav] += traces_calculator_jack(4,traces_id4, mat_restot, iflav,ncopies,nstart,i, nind_trace);		

                     
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
                        int nconf_eff=nconf-block_size;
			Tr_M_dM2[iflav] = (2.0*Tr_M_dM2[iflav])/((double)nconf_eff*(ncopies)*(ncopies -1));
			Tr_M_dM[iflav]= ((1.0)*Tr_M_dM[iflav])/((double)nconf_eff*(ncopies));
			Tr_M_dM_M_dM[iflav]= ((1.0)*Tr_M_dM_M_dM[iflav])/((double)nconf_eff*(ncopies));
			Tr_M_d2M[iflav]= ((1.0)*Tr_M_d2M[iflav])/((double)nconf_eff*(ncopies));
			Tr_M_dM_M_d2M[iflav]= ((1.0)*Tr_M_dM_M_d2M[iflav])/((double)nconf_eff*(ncopies));	
			Tr_M_dM_M_dM_M_dM[iflav]= ((1.0)*Tr_M_dM_M_dM_M_dM[iflav])/((double)nconf_eff*(ncopies));
			Tr_M_dM_Tr_M_dM_M_dM[iflav] = ((1.0)*Tr_M_dM_Tr_M_dM_M_dM[iflav])/((double)nconf_eff*(ncopies)*(ncopies -1));
			Tr_M_dM_Tr_M_d2M[iflav] = ((1.0)*Tr_M_dM_Tr_M_d2M[iflav])/((double)nconf_eff*(ncopies)*(ncopies -1));	
			Tr_M_dM3[iflav] = ((6.0)*Tr_M_dM3[iflav])/((double)nconf_eff*(ncopies)*(ncopies -1)*(ncopies -2));	
                        Tr_M_d2M_M_d2M[iflav]*= (1.0)/((double)nconf_eff*(ncopies));
                        Tr_M_dM_M_dM_M_d2M[iflav]*= (1.0)/((double)nconf_eff*(ncopies));
                        Tr_M_dM_M_dM_M_dM_M_dM[iflav]*= (1.0)/((double)nconf_eff*(ncopies));
                        Tr_M_dM_M_dM_Tr_M_dM_Tr_M_dM[iflav]*=(2.0)/((double)nconf_eff*(ncopies)*(ncopies-1)*(ncopies-2));
                        Tr_M_dM_Tr_M_dM_Tr_M_dM_Tr_M_dM[iflav]*=(24.0)/((double)nconf_eff*(ncopies)*(ncopies-1)*(ncopies-2)*(ncopies-3));
                        Tr_M_dM_Tr_M_dM_Tr_M_d2M[iflav]*= (2.0)/((double)nconf_eff*(ncopies)*(ncopies-1)*(ncopies-2));
                        Tr_M_dM_M_dM_Tr_M_dM_M_dM[iflav]*=(2.0)/((double)nconf_eff*(ncopies)*(ncopies-1));
                        Tr_M_dM_M_dM_M_dM_Tr_M_dM[iflav]*= (1.0)/((double)nconf_eff*(ncopies)*(ncopies-1));
                        Tr_M_d2M_Tr_M_d2M[iflav]*=(2.0)/((double)nconf_eff*(ncopies)*(ncopies-1));	
                        Tr_M_d2M_Tr_M_dM_M_dM[iflav] *= (1.0)/((double)nconf_eff*(ncopies)*(ncopies-1));
                        Tr_M_dM_Tr_M_dM_M_d2M[iflav] *= (1.0)/((double)nconf_eff*(ncopies)*(ncopies-1));
		
		
		
			
			
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
	+((1.0)/(16))*Tr_M_dM2[iflav]
	+((3.0)/(16))*Tr_M_d2M_Tr_M_d2M[iflav]
	-((33.0)/(64))*Tr_M_dM_Tr_M_dM_M_d2M[iflav]
	-((3.0)/(16))*Tr_M_d2M[iflav]*Tr_M_d2M[iflav]
	+((3.0)/(16))*Tr_M_d2M[iflav]*Tr_M_dM_M_dM[iflav]
	-((3.0)/(64))*Tr_M_d2M[iflav]*Tr_M_dM2[iflav]
	-((3.0)/(4))*Tr_M_d2M_M_d2M[iflav]
	+((3.0)/(1))*Tr_M_dM_M_dM_M_d2M[iflav]
	-((1.0)/(1))*Tr_M_dM_M_dM[iflav]
	-((3.0)/(2))*Tr_M_dM_M_dM_M_dM_M_dM[iflav]
	+((1.0)/(4))*Tr_M_d2M[iflav]
	+((3.0)/(16))*Tr_M_dM_M_dM[iflav]*Tr_M_d2M[iflav]
	-((3.0)/(16))*Tr_M_dM_M_dM[iflav]*Tr_M_dM_M_dM[iflav]
	+((3.0)/(64))*Tr_M_dM_M_dM[iflav]*Tr_M_dM2[iflav];

			
		susc_block_diag *= (1.0)/V_4;	
                if(iflav==0)printf("%lf %d \n", creal(susc_block_diag),nconf-block_size);
			
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









