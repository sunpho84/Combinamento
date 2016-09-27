#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "jackknife4.cpp"

typedef double complex dcomplex;
const double N_t=8;
  const int ncopies=256;
  const int nflavs=3;
  const int nind_trace=9;
  const int conf_out_eq=0;
  const int ncopies_eff=256;

double complex traces_calculator(int num_prod, int *traces_id, double complex **MAT_RES, int iflav)
{
double complex result=0;

double complex X[ncopies_eff][num_prod-1];
//set X1[ncopies_eff-1][dim] to 0
for (int i=0; i<num_prod-1;i++)
 for (int icopy=0;icopy<ncopies_eff;icopy++)   X[icopy][i]=0;

for (int i=0;i<num_prod-1;i++)
{  
   if(i==0)
   for(int icopy=ncopies_eff-2;icopy>=0;icopy--)
   X[icopy][i] = X[icopy+1][i] + MAT_RES[icopy+1][iflav*nind_trace + traces_id[num_prod-1-i]];

   else
   for (int icopy=ncopies_eff-2-i;icopy>=0;icopy--)
    X[icopy][i] = X[icopy+1][i] + MAT_RES[icopy+1][iflav*nind_trace+traces_id[num_prod-1-i]]*X[icopy+1][i-1];

}

for(int icopy=0;icopy<ncopies_eff;icopy++) result += X[icopy][num_prod-2]*MAT_RES[icopy][iflav*nind_trace + traces_id[0]];


return result;



}


int main ()
{
  int traces_id2[2];
  int traces_id3[3];
  int traces_id4[4];

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
  //definisci i vettori per le tracce delle suscettività miste

  double complex err_diag[nflavs];
  double complex susc_jack_diag[nflavs];
  double complex susc_jack_diag_2nd[nflavs];

  //init to zero

  for(int iflav=0;iflav<nflavs;iflav++)
    {
      Tr_M_dM[iflav]=Tr_M_d2M[iflav]=Tr_M_dM_M_dM[iflav]=Tr_M_dM2[iflav]=Tr_M_dM_M_d2M[iflav]=0;
      Tr_M_d2M_M_d2M[iflav]= Tr_M_dM_M_dM_M_d2M[iflav]= Tr_M_dM_M_dM_M_dM_M_dM[iflav]= Tr_M_dM_M_dM_Tr_M_dM_Tr_M_dM[iflav]=        Tr_M_dM_Tr_M_dM_Tr_M_dM_Tr_M_dM[iflav]=Tr_M_d2M_Tr_M_dM_M_dM[iflav]=Tr_M_dM_Tr_M_dM_M_d2M[iflav]=0;
      Tr_M_dM_Tr_M_dM_Tr_M_d2M[iflav]= Tr_M_dM_M_dM_M_dM_Tr_M_dM[iflav]= Tr_M_d2M_Tr_M_d2M[iflav]=Tr_M_dM_M_dM_Tr_M_dM_M_dM[iflav]=0;
      Tr_M_dM_M_dM_M_dM[iflav]=Tr_M_dM3[iflav]=Tr_M_dM_Tr_M_d2M[iflav]=Tr_M_dM_Tr_M_dM_M_dM[iflav]=err_diag[iflav]=0;
    }


  //loop until reach end of file
   double complex **MAT_RES;
  MAT_RES= (double complex **)malloc(ncopies_eff*sizeof(double complex*));
  for(int i=0;i<ncopies_eff;i++)
    {
      MAT_RES[i]= (double complex *)malloc(nflavs*nind_trace*sizeof(double complex));
    }

     //set MAT_RES to zero
    for(int icopy=0;icopy<ncopies_eff;icopy++)
       for(int k=0;k<nflavs*nind_trace;k++) MAT_RES[icopy][k]=0;
  FILE *file_in=open_file("rende_new_run1.txt","r");
  int nconfs=0;
  int eof=0;
  int puppa;
  int err_num_copies=0;
  int id_conf=0;
  while(!eof)
    {
      //init to zero temp traces



      //loop over copies
      for(int icopy=0;icopy<ncopies;icopy++)
	{
          if(err_num_copies==1) printf("%d\n",err_num_copies);
	  //skip iconf: if not read it means we reached eof

		 if(err_num_copies==0){
			//skip iconf: if not read it means we reached eof
			if(fscanf(file_in,"%d",&puppa)!=1) eof=1;
			else
			{
			if(icopy==0) id_conf=puppa;
            if(puppa!=id_conf) err_num_copies=1;

			}
			}

           else  if(icopy==0 && err_num_copies==1 )
		 {
		 err_num_copies=0;
		 id_conf=puppa;
		 }


	  if(!eof && !err_num_copies)
	    for(int j=0;j<nflavs*nind_trace;j++)
	      {
		//load real and immaginary part
		double re,im;
		if(fscanf(file_in,"%lf %lf",&re,&im)!=2)
		  {
		    fprintf(stderr,"error reading re/im for copy %d flavour/nind_trace %d\n",icopy,j);
		    exit(1);
		  }

		//put in place
		if(icopy<ncopies_eff && nconfs>=conf_out_eq) MAT_RES[icopy][j]=re+I*im;
	      	}
	      	}
      if(!eof && !err_num_copies)
	{
	  if(nconfs>=conf_out_eq)
	    {


	      //loop over flavour
	      for(int iflav=0;iflav<nflavs;iflav++)
		{
		  //compute all mean values

		  //compute TR_M_dM
		  for(int icopy=0;icopy<ncopies_eff;icopy++)
		    {
		      Tr_M_dM[iflav]+=     MAT_RES[icopy][iflav*nind_trace+0];
		      //store the trace computed over one config. Will be used when computing Tr_(dMM)_i_Tr(dMM)_j

		    }


		  //compute TR_M_d2M
		  for(int icopy=0;icopy<ncopies_eff;icopy++)
		    {
		      Tr_M_d2M[iflav]+=    MAT_RES[icopy][iflav*nind_trace+1];

		    }






		  //compute TR_M_dM_M_dM
		  for(int icopy=0;icopy<ncopies_eff;icopy++)
		    {
		      Tr_M_dM_M_dM[iflav]+=MAT_RES[icopy][iflav*nind_trace+2];


		    }
		  //compute TR_M_dM_M_d2M
		  for(int icopy=0;icopy<ncopies_eff;icopy++) Tr_M_dM_M_d2M[iflav]+=MAT_RES[icopy][iflav*nind_trace+4];

		  //compute TR_M_dM_M_dM_M_dM
		  for(int icopy=0;icopy<ncopies_eff;icopy++) Tr_M_dM_M_dM_M_dM[iflav]+=MAT_RES[icopy][iflav*nind_trace+5];

		  //compute (TrM_dM)^2
		  traces_id2[0]=0;
                  traces_id2[1]=0;
                  Tr_M_dM2[iflav]+= traces_calculator(2,traces_id2, MAT_RES, iflav);


			//store the trace computed over one config. Will be used when computing the products of traces of different flavour matrix

		  //compute (TrM_dM)(Tr(M_dM)^2)
		  for(int icopy=0;icopy<ncopies_eff;icopy++)
		    for(int rcopy=icopy+1;rcopy<ncopies_eff;rcopy++)
		      {
			double complex complex1=MAT_RES[icopy][iflav*nind_trace+0];
			double complex complex2=MAT_RES[rcopy][iflav*nind_trace+2];
			Tr_M_dM_Tr_M_dM_M_dM[iflav]+=complex1*complex2;
			//compute the products for rcopy<icopy
			complex1=MAT_RES[rcopy][iflav*nind_trace+0];
			complex2=MAT_RES[icopy][iflav*nind_trace+2];
			Tr_M_dM_Tr_M_dM_M_dM[iflav]+=complex1*complex2;
		      }
		  //compute (TrM_dM)(TrM_d2M)
		  for(int icopy=0;icopy<ncopies_eff;icopy++)
		    for(int rcopy=icopy+1;rcopy<ncopies_eff;rcopy++)
		      {
			double complex complex1=MAT_RES[icopy][iflav*nind_trace+0];
			double complex complex2=MAT_RES[rcopy][iflav*nind_trace+1];
			Tr_M_dM_Tr_M_d2M[iflav]+=complex1*complex2;
			//compute the products for rcopy<ncopy
			complex1=MAT_RES[rcopy][iflav*nind_trace+0];
			complex2=MAT_RES[icopy][iflav*nind_trace+1];
			Tr_M_dM_Tr_M_d2M[iflav]+=complex1*complex2;

		      }

		  //compute (TrM_dM)^3

                  traces_id3[0]=0; traces_id3[1]=0; traces_id3[2]=0;
                  Tr_M_dM3[iflav]+= traces_calculator(3,traces_id3, MAT_RES, iflav);


		  //Traces at order 4

		  //compute Tr_M_d2M_M_d2M

		  for (int icopy=0;icopy<ncopies_eff;icopy++)
		    {
		      Tr_M_d2M_M_d2M[iflav] += MAT_RES[icopy][iflav*nind_trace + 6];
		    }

		  //compute Tr_M_d2M_M_dM_M_dM

		  for (int icopy=0;icopy<ncopies_eff;icopy++)
		    {
		      Tr_M_dM_M_dM_M_d2M[iflav] += MAT_RES[icopy][iflav*nind_trace + 7];
		    }

		  //compute Tr_M_dM_M_dM_M_dM_M_dM

		  for (int icopy=0;icopy<ncopies_eff;icopy++)
		    {
		      Tr_M_dM_M_dM_M_dM_M_dM[iflav] += MAT_RES[icopy][iflav*nind_trace + 8];
		    }



		  //        compute Tr_M_dM_Tr_M_dM_M_d2M

		  for(int icopy=0;icopy<ncopies_eff;icopy++)
		    for(int rcopy=icopy+1;rcopy<ncopies_eff;rcopy++)
		      {
			double complex complex1=MAT_RES[icopy][iflav*nind_trace+0];
			double complex complex2=MAT_RES[rcopy][iflav*nind_trace+4];
			Tr_M_d2M_Tr_M_dM_M_dM[iflav]+=complex1*complex2;
			//compute the products for rcopy<ncopy
			complex1=MAT_RES[rcopy][iflav*nind_trace+0];
			complex2=MAT_RES[icopy][iflav*nind_trace+4];
			Tr_M_d2M_Tr_M_dM_M_dM[iflav]+=complex1*complex2;

		      }



		  //compute Tr_M_d2M_Tr_M_dM_M_dM

		  for(int icopy=0;icopy<ncopies_eff;icopy++)
		    for(int rcopy=icopy+1;rcopy<ncopies_eff;rcopy++)
		      {
			double complex complex1=MAT_RES[icopy][iflav*nind_trace+1];
			double complex complex2=MAT_RES[rcopy][iflav*nind_trace+2];
			Tr_M_d2M_Tr_M_dM_M_dM[iflav]+=complex1*complex2;
			//compute the products for rcopy<ncopy
			complex1=MAT_RES[rcopy][iflav*nind_trace+1];
			complex2=MAT_RES[icopy][iflav*nind_trace+2];
			Tr_M_d2M_Tr_M_dM_M_dM[iflav]+=complex1*complex2;

		      }
		  //compute Tr_M_dM_M_dM_Tr_M_dM_M_dM

		  for(int icopy=0;icopy<ncopies_eff;icopy++)
		    for(int rcopy=icopy+1;rcopy<ncopies_eff;rcopy++)
		      {
			double complex complex1=MAT_RES[icopy][iflav*nind_trace+2];
			double complex complex2=MAT_RES[rcopy][iflav*nind_trace+2];
			Tr_M_dM_M_dM_Tr_M_dM_M_dM[iflav]+=complex1*complex2;
		      }

		  //compute Tr_M_dM_M_dM_M_dM_Tr_M_dM

		  for(int icopy=0;icopy<ncopies_eff;icopy++)
		    for(int rcopy=icopy+1;rcopy<ncopies_eff;rcopy++)
		      {
			double complex complex1=MAT_RES[icopy][iflav*nind_trace+5];
			double complex complex2=MAT_RES[rcopy][iflav*nind_trace+0];
			Tr_M_dM_M_dM_M_dM_Tr_M_dM[iflav]+=complex1*complex2;
			//compute the products for rcopy<ncopy
			complex1=MAT_RES[rcopy][iflav*nind_trace+5];
			complex2=MAT_RES[icopy][iflav*nind_trace+0];
			Tr_M_dM_M_dM_M_dM_Tr_M_dM[iflav]+=complex1*complex2;

		      }


		  //compute Tr_M_d2M_Tr_M_d2M

		  for(int icopy=0;icopy<ncopies_eff;icopy++)
		    for(int rcopy=icopy+1;rcopy<ncopies_eff;rcopy++)
		      {
			double complex complex1=MAT_RES[icopy][iflav*nind_trace+1];
			double complex complex2=MAT_RES[rcopy][iflav*nind_trace+1];
			Tr_M_d2M_Tr_M_d2M[iflav]+=complex1*complex2;
		      }

		  //compute Tr_M_dM_M_dM_Tr_M_dM_Tr_M_dM
                  traces_id3[0]=2; traces_id3[1]=0; traces_id3[2]=0;
                  Tr_M_dM_M_dM_Tr_M_dM_Tr_M_dM[iflav]+= traces_calculator(3,traces_id3, MAT_RES, iflav);
                  traces_id3[0]=0; traces_id3[1]=2; traces_id3[2]=0;
                  Tr_M_dM_M_dM_Tr_M_dM_Tr_M_dM[iflav]+= traces_calculator(3,traces_id3, MAT_RES, iflav);
                  traces_id3[0]=0; traces_id3[1]=0; traces_id3[2]=2;
                  Tr_M_dM_M_dM_Tr_M_dM_Tr_M_dM[iflav]+= traces_calculator(3,traces_id3, MAT_RES, iflav);

		  //compute Tr_M_d2M_Tr_M_dM_Tr_M_dM


		  traces_id3[0]=0; traces_id3[1]=0; traces_id3[2]=1;

                  Tr_M_dM_Tr_M_dM_Tr_M_d2M[iflav]+= traces_calculator(3,traces_id3, MAT_RES, iflav);
                  traces_id3[0]=0; traces_id3[1]=1; traces_id3[2]=0;

                  Tr_M_dM_Tr_M_dM_Tr_M_d2M[iflav]+= traces_calculator(3,traces_id3, MAT_RES, iflav);
                  traces_id3[0]=1; traces_id3[1]=0; traces_id3[2]=0;

                  Tr_M_dM_Tr_M_dM_Tr_M_d2M[iflav]+= traces_calculator(3,traces_id3, MAT_RES, iflav);

		  //compute Tr_M_dM_Tr_M_dM_Tr_M_dM_Tr_M_dM


                  traces_id4[0]=0; traces_id4[1]=0; traces_id4[2]=0; traces_id4[3]=0;
                  Tr_M_dM_Tr_M_dM_Tr_M_dM_Tr_M_dM[iflav] += traces_calculator(4,traces_id4, MAT_RES, iflav);


		  //** FINE DELLA PARTE DI CALCOLO DELLE TRACCE **//


		  //exit from loop over flavour
		}



	      //exit from loop over copies
	    }
	  printf("%d \n",nconfs);
	  nconfs++;

	}
      //end of the reading part
    }
  //close and print the number of confs
  fclose(file_in);
  printf("Nconfs read: %d\n",nconfs);
  printf("Nconfs used: %d\n", nconfs- conf_out_eq);
  nconfs -= conf_out_eq;


  //save the traces before normalization (for jackknife)
  const int n_traces=20;
  double complex Tr_S[n_traces][nflavs];
  for(int iflav=0;iflav<nflavs;iflav++)
    {
      Tr_S[0][iflav]= Tr_M_dM[iflav];
      Tr_S[1][iflav]= Tr_M_d2M[iflav];
      Tr_S[2][iflav]= Tr_M_dM_M_dM[iflav];
      Tr_S[3][iflav]= Tr_M_dM2[iflav];
      Tr_S[4][iflav]= Tr_M_dM3[iflav];
      Tr_S[5][iflav]= Tr_M_dM_Tr_M_dM_M_dM[iflav];
      Tr_S[6][iflav]= Tr_M_dM_Tr_M_d2M[iflav];
      Tr_S[7][iflav]= Tr_M_dM_M_d2M[iflav];
      Tr_S[8][iflav]= Tr_M_dM_M_dM_M_dM[iflav];
      Tr_S[9][iflav]= Tr_M_d2M_M_d2M[iflav];
      Tr_S[10][iflav]= Tr_M_dM_M_dM_M_d2M[iflav];
      Tr_S[11][iflav]= Tr_M_dM_M_dM_M_dM_M_dM[iflav];
      Tr_S[12][iflav]= Tr_M_dM_M_dM_Tr_M_dM_Tr_M_dM[iflav];
      Tr_S[13][iflav]= Tr_M_dM_Tr_M_dM_Tr_M_dM_Tr_M_dM[iflav];
      Tr_S[14][iflav]= Tr_M_dM_Tr_M_dM_Tr_M_d2M[iflav];
      Tr_S[15][iflav]= Tr_M_dM_M_dM_Tr_M_dM_M_dM[iflav];
      Tr_S[16][iflav]= Tr_M_dM_M_dM_M_dM_Tr_M_dM[iflav];
      Tr_S[17][iflav]= Tr_M_d2M_Tr_M_d2M[iflav];
      Tr_S[18][iflav]= Tr_M_d2M_Tr_M_dM_M_dM[iflav];
      Tr_S[19][iflav]= Tr_M_dM_Tr_M_dM_M_d2M[iflav];
    }




  const double V_4=32.0*32.0*4.0;

  // normalize the traces
  for(int iflav=0;iflav<nflavs;iflav++)
    {

  Tr_M_dM2[iflav] = (2.0*Tr_M_dM2[iflav])/((double)nconfs*(ncopies_eff)*(ncopies_eff -1));
      printf("Tr_M_dM2: %lf %lf \n", creal(Tr_M_dM2[0]), cimag(Tr_M_dM2[0]));
      Tr_M_dM[iflav]= ((1.0)*Tr_M_dM[iflav])/((double)nconfs*(ncopies_eff));
      printf("Tr_M_dM: %lf %lf \n", creal(Tr_M_dM[0]), cimag(Tr_M_dM[0]));
      Tr_M_dM_M_dM[iflav]= ((1.0)*Tr_M_dM_M_dM[iflav])/((double)nconfs*(ncopies_eff));
      printf("Tr_M_dM_M_dM: %lf %lf \n", creal(Tr_M_dM_M_dM[0]), cimag(Tr_M_dM_M_dM[0]));
      Tr_M_d2M[iflav]= ((1.0)*Tr_M_d2M[iflav])/((double)nconfs*(ncopies_eff));
      printf("Tr_M_d2M: %lf %lf \n", creal(Tr_M_d2M[0]), cimag(Tr_M_d2M[0]));
      Tr_M_dM_M_d2M[iflav]= ((1.0)*Tr_M_dM_M_d2M[iflav])/((double)nconfs*(ncopies_eff));
      printf("Tr_M_dM_M_d2M: %lf %lf \n", creal(Tr_M_dM_M_d2M[0]), cimag(Tr_M_dM_M_d2M[0]));
      Tr_M_dM_M_dM_M_dM[iflav]= ((1.0)*Tr_M_dM_M_dM_M_dM[iflav])/((double)nconfs*(ncopies_eff));
      printf("Tr_M_dM_M_dM_M_dM: %lf %lf \n", creal(Tr_M_dM_M_dM_M_dM[0]), cimag(Tr_M_dM_M_dM_M_dM_M_dM[0]));
      Tr_M_dM_Tr_M_dM_M_dM[iflav] = ((1.0)*Tr_M_dM_Tr_M_dM_M_dM[iflav])/((double)nconfs*(ncopies_eff)*(ncopies_eff -1));
      printf("Tr_M_dM_Tr_M_dM_M_dM: %lf %lf \n", creal(Tr_M_dM_Tr_M_dM_M_dM[0]), cimag(Tr_M_dM_Tr_M_dM_M_dM[0]));
      Tr_M_dM_Tr_M_d2M[iflav] = ((1.0)*Tr_M_dM_Tr_M_d2M[iflav])/((double)nconfs*(ncopies_eff)*(ncopies_eff -1));
      printf("Tr_M_dM_Tr_M_d2M: %lf %lf \n", creal(Tr_M_dM_Tr_M_d2M[0]), cimag(Tr_M_dM_Tr_M_d2M[0]));
      Tr_M_dM3[iflav] = ((6.0)*Tr_M_dM3[iflav])/((double)nconfs*(ncopies_eff)*(ncopies_eff -1)*(ncopies_eff -2));
      printf("Tr_M_dM3: %lf %lf \n", creal(Tr_M_dM3[0]), cimag(Tr_M_dM3[0]));
      Tr_M_d2M_M_d2M[iflav]*= (1.0)/((double)nconfs*(ncopies_eff));
      printf("Tr_M_d2M_M_d2M: %lf %lf \n", creal(Tr_M_d2M_M_d2M[0]), cimag(Tr_M_d2M_M_d2M[0]));
      Tr_M_dM_M_dM_M_d2M[iflav]*= (1.0)/((double)nconfs*(ncopies_eff));
      printf("Tr_M_dM_M_dM_M_d2M: %lf %lf \n", creal(Tr_M_dM_M_dM_M_d2M[0]), cimag(Tr_M_dM_M_dM_M_d2M[0]));
      Tr_M_dM_M_dM_M_dM_M_dM[iflav]*= (1.0)/((double)nconfs*(ncopies_eff));
      printf("Tr_M_dM_M_dM_M_dM_M_dM: %lf %lf \n", creal(Tr_M_dM_M_dM_M_dM_M_dM[0]), cimag(Tr_M_dM_M_dM_M_dM_M_dM[0]));
      Tr_M_dM_M_dM_Tr_M_dM_Tr_M_dM[iflav]*=(2.0)/((double)nconfs*(ncopies_eff)*(ncopies_eff-1)*(ncopies_eff-2));
      printf("Tr_M_dM_M_dM_Tr_M_dM_Tr_M_dM: %lf %lf \n", creal(Tr_M_dM_M_dM_Tr_M_dM_Tr_M_dM[0]), cimag(Tr_M_dM_M_dM_Tr_M_dM_Tr_M_dM[0]));
      Tr_M_dM_Tr_M_dM_Tr_M_dM_Tr_M_dM[iflav]*=(24.0)/((double)nconfs*(ncopies_eff)*(ncopies_eff-1)*(ncopies_eff-2)*(ncopies_eff-3));
      printf("Tr_M_dM_Tr_M_dM_Tr_M_dM_Tr_M_dM: %lf %lf \n", creal(Tr_M_dM_Tr_M_dM_Tr_M_dM_Tr_M_dM[0]), cimag(Tr_M_dM_Tr_M_dM_Tr_M_dM_Tr_M_dM[0]));
      Tr_M_dM_Tr_M_dM_Tr_M_d2M[iflav]*= (2.0)/((double)nconfs*(ncopies_eff)*(ncopies_eff-1)*(ncopies_eff-2));
      printf("Tr_M_d2M_Tr_M_dM_Tr_M_dM: %lf %lf \n", creal(Tr_M_dM_Tr_M_dM_Tr_M_d2M[0]), cimag(Tr_M_dM_Tr_M_dM_Tr_M_d2M[0]));
      Tr_M_dM_M_dM_Tr_M_dM_M_dM[iflav]*=(2.0)/((double)nconfs*(ncopies_eff)*(ncopies_eff-1));
      printf("Tr_M_dM_M_dM_Tr_M_dM_M_dM: %lf %lf \n", creal(Tr_M_dM_M_dM_Tr_M_dM_M_dM[0]), cimag(Tr_M_dM_M_dM_Tr_M_dM_M_dM[0]));
      Tr_M_dM_M_dM_M_dM_Tr_M_dM[iflav]*= (1.0)/((double)nconfs*(ncopies_eff)*(ncopies_eff-1));
      Tr_M_d2M_Tr_M_d2M[iflav]*=(2.0)/((double)nconfs*(ncopies_eff)*(ncopies_eff-1));
      Tr_M_d2M_Tr_M_dM_M_dM[iflav] *= (1.0)/((double)nconfs*(ncopies_eff)*(ncopies_eff-1));
      printf("Tr_M_d2M_Tr_M_dM_M_dM: %lf %lf \n", creal(Tr_M_d2M_Tr_M_dM_M_dM[0]), cimag(Tr_M_d2M_Tr_M_dM_M_dM[0]));
      Tr_M_dM_Tr_M_dM_M_d2M[iflav] *= (1.0)/((double)nconfs*(ncopies_eff)*(ncopies_eff-1));




    }


  //// end of normalization



  //compute susceptibilities

  //compute and store diagonal susceptibilities
  for(int iflav=0;iflav<nflavs;iflav++)
    {


      double complex susc_diag= -((6.0)/(64))*Tr_M_dM_M_dM_Tr_M_dM_Tr_M_dM[iflav]
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

      susc_diag *= (1.0)/V_4;

      //store susc for jackknife
      susc_jack_diag[iflav]=susc_diag;
      printf("susc%d = %lf \n",iflav, creal(susc_jack_diag[iflav]));

    }


//compute diagonal susceptibilities
	for(int iflav=0; iflav<nflavs;iflav++)
	   {
		double complex susc_disc=(1.0/(16*V_4))*(Tr_M_dM2[iflav]-Tr_M_dM[iflav]*Tr_M_dM[iflav]);
		double complex susc_conn=(1.0/(4*V_4))*(Tr_M_d2M[iflav]-Tr_M_dM_M_dM[iflav]);
		double complex tot=susc_conn+susc_disc;
		susc_jack_diag_2nd[iflav]=tot;
                printf("%lf %d \n", creal(susc_jack_diag_2nd[iflav]),iflav);
	   }



  //Qui memorizzo di nuovo le tracce. Salvo tutto in un matricione perchè mi serve poi nel Jackknife quando devo fare le medie sui blocchi.
  FILE *file_in2= open_file("rende_new_run1.txt","r");
   double complex **mat_restot;
  mat_restot= (double complex **)malloc(nconfs*ncopies_eff*sizeof(double complex*));
  for(int i=0;i<nconfs*ncopies_eff;i++)
    {
      mat_restot[i]= (double complex *)malloc(nflavs*nind_trace*sizeof(double complex));
    }

    for(int icopy=0;icopy<ncopies_eff*nconfs;icopy++)
       for(int k=0;k<nflavs*nind_trace;k++) mat_restot[icopy][k]=0;
  int iconf=0;
  eof=0;
  while(!eof)
    {
      //loop over copies
      for(int icopy=0;icopy<ncopies;icopy++)
	{
	  //skip iconf: if not read it means we reached eof

		 if(err_num_copies==0){
			//skip iconf: if not read it means we reached eof
			if(fscanf(file_in2,"%d",&puppa)!=1) eof=1;
			else
			{
			if(icopy==0) id_conf=puppa;
                        if(puppa!=id_conf) err_num_copies=1;

			}
			}
                   else if(icopy==0 && err_num_copies==1 )
		 {
		 err_num_copies=0;
		 id_conf=puppa;
		 }

	  if(!eof && !err_num_copies)
	    for(int j=0;j<nflavs*nind_trace;j++)
	      {
		//load real and immaginary part
		double re,im;
		if(fscanf(file_in2,"%lf %lf",&re,&im)!=2)
		  {
		    fprintf(stderr,"error reading re/im for copy %d flavour/nind_trace %d\n",icopy,j);
		    exit(1);
		  }

		//put in place
		// 2nd cond is to skip non-termalized confs
		if(iconf>=conf_out_eq && icopy<ncopies_eff) mat_restot[icopy + (iconf-conf_out_eq)*ncopies_eff ][j]=re+I*im;
	      }


	}
      if(err_num_copies)
		{
		   for(int rcopy=0;rcopy<ncopies;rcopy++)
		   {
		   for(int k=0;k<nflavs*nind_trace;k++)
				{
					//load real and immaginary part
					double re,im;
					if(fscanf(file_in2,"%lf %lf",&re,&im)!=2)
					{
						fprintf(stderr,"error reading re/im for copy %d flavour/nind_trace %d\n",rcopy,k);
						exit(1);
					}

					//put in place
					// 2nd cond is to skip non-termalized confs
					if(iconf>=conf_out_eq) mat_restot[rcopy + (iconf-conf_out_eq)*ncopies ][k]=0;
				}


		   }


		}

      if(!eof && !err_num_copies) iconf++;



    }
  
  fclose(file_in2);
  printf("%d \n",iconf);
  printf("%d \n", iconf -conf_out_eq);
  //calcolo l'errore sulle suscettivita del terz'ordine
  int block_size_max= 12;
  int file_created=0;
  for(int block_size=1; block_size< block_size_max; block_size++)
    {
      int error=0;
      error=jackknife(file_created,err_diag,ncopies_eff,block_size,nconfs,Tr_S,susc_jack_diag,mat_restot);
      if (error!=0)
	{
	  printf("errore nel blocco: %d", block_size);
	  exit(1);
	}
      file_created=1;
      printf("block_size:%d \n",block_size);
    }
  //
  //compute the average of the jackknife estimation of variance (average over block size);
  for(int iflav=0; iflav<nflavs;iflav++)
    {
      err_diag[iflav] /= block_size_max -4;

    }

  //


  //stampo i risultati
  //memorizzo il potenziale chimico

  FILE *file_pot= open_file("pot.txt", "r");  //ora questa parte è inutile ma serviva prima quando c'erano tanti run a vari potenziali chimici
  double u,d,s;

  if(fscanf(file_pot,"%lf %lf %lf",&u,&d,&s)!=3)
    {
      fprintf(stderr,"errore nel leggere i potenziali chimici");
      exit(1);

    }
  fclose(file_pot);


  //create file and store the results
  //check di esistenza del file
  int b=0;
  FILE *file_check= fopen("suscettivita_4_ord.txt","r");
  if (file_check==NULL) b=1;
  else fclose(file_check);


  FILE *file_out_bis=open_file("suscettivita_4_ord.txt",b?"w":"a");

  //stampa le suscettività
  for(int iflav=0;iflav<nflavs;iflav++)  fprintf(file_out_bis,"%lf %lf %lf %lf \t",creal(susc_jack_diag[iflav])/(N_t*N_t), creal(err_diag[iflav])/(N_t*N_t),cimag(susc_jack_diag[iflav])/(N_t*N_t), cimag(err_diag[iflav])/(N_t*N_t));

  fprintf(file_out_bis,"\n");
  fclose(file_out_bis);
  //disalloco la memoria
  for(int i=0;i<nconfs*ncopies_eff;i++) free(mat_restot[i]);

  free(mat_restot);

  for(int i=0;i<ncopies_eff;i++) free(MAT_RES[i]);

  free(MAT_RES);

  return 0;
}
