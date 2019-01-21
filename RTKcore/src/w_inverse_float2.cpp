/////////////////////////////////////////////////////////////////////
//
// ŠÏ‘ªŒë·‹¤•ªUs—ñ‚ğ‹‚ß‚é•”•ª
//
////////////////////////////////////////////////////////////////////
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "global_extern.h"

int minv(double[],double,int);

void w_inverse_float2(double w1[PRN*PRN], int satn, int flag, int nsvg, int nsvb, int nsvc, int nsvr)
{
	int i,j,k,l,m,n;
	double eps;
//	double w1[48*48]={0};
	double std1 = Carrier_noise;//gps qzs
	double std2 = Code_noise;//gps qzs
	double std1b = Carrier_noise;//other
	double std2b = Code_noise;//other
	double std3 = 0.15;

	eps=1.0e-30;

	for(i=0;i<satn*4*satn*4;i++)
		w1[i]=0;
///*
	int o=0;
	l=0;m=0;n=0;
	if(flag==1){
		//satn‚ÍŠÏ‘ª‰q¯”-1
		//L1 ‚Ì‚İ‚Ìê‡A[satn*2]*[satn*2]‚Ìs—ñ
		for(i=0;i<satn*2*satn*2;i++){
			if((i+satn*2+1)%(satn*2+1)==0){
				if(i<satn*2*satn*2/2){

					w1[i]=4*OmomiC[l];
					l++;
					if(l>nsvg){
						w1[i]=4*OmomiCB[m];
						m++;
					}


//					w1[i]=4*std1*std1*OmomiP[l];
//					l++;
//					if(l>nsvg){
//						w1[i]=4*std1*std1*OmomiC[m];
//						m++;
//					}
				}
				if(i>=satn*2*satn*2/2){

					w1[i]=4*OmomiP[n];
					n++;
					if(n>nsvg){
						w1[i]=4*OmomiPB[o];
						o++;
					}

//					w1[i]=4*std2*std2*OmomiP[n];
//					n++;
//					if(n>nsvg){
//						w1[i]=4*std2*std2*OmomiC[o];
//						o++;
//					}
				}
			}
			else if(i<satn*2*satn){
				for(j=0;j<satn;j++){
					if(i<satn+j*satn*2 && i>=j*satn*2)
						w1[i]=2*std1*std1;
				}
			}
			else if(i>=satn*2*satn && i<satn*2*satn*2){
				for(j=0;j<satn;j++){
					if(i<satn+(j+satn)*satn*2+satn && i>=(j+satn)*satn*2+satn)
						w1[i]=2*std2*std2;
				}
			}
		}
	}
//	*/

//	
	l=0;m=0;
	if(flag==2){
		//satn‚ÍŠÏ‘ª‰q¯”-1
		for(i=0;i<satn*4*satn*4;i++){
			if((i+satn*4+1)%(satn*4+1)==0){
				if(i<satn*4*satn*4/2){
					w1[i]=4*std1*std1;
///*
					w1[i]=4*OmomiC[l];
					l++;
					if(l>=nsvg+1){
						w1[i]=4*OmomiCB[l-nsvg-1];
					}
					if(l==satn)
						l=0;
//*/
				}
				if(i>=satn*4*satn*4/2){
					w1[i]=4*std2*std2;
///*					
					w1[i]=4*OmomiP[m];
					m++;
					if(m>=nsvg+1){
						w1[i]=4*OmomiPB[m-nsvg-1];
					}
					if(m==satn)
						m=0;
//*/
				}
			}
			else if(i<satn*4*satn){
				for(j=0;j<satn;j++){
					if(i<satn+j*satn*4 && i>=j*satn*4)
						w1[i]=2*std1*std1;
				}
			}
			else if(i>=satn*4*satn && i<satn*4*satn*2){
				for(j=0;j<satn;j++){
					if(i<satn+(j+satn)*satn*4+satn && i>=(j+satn)*satn*4+satn)
						w1[i]=2*std1*std1;
				}
			}
			else if(i>=satn*4*satn*2 && i<satn*4*satn*3){
				for(j=0;j<satn;j++){
					if(i<satn+(j+satn*2)*satn*4+satn+satn && i>=(j+satn*2)*satn*4+satn+satn)
						w1[i]=2*std2*std2;
				}
			}
			else if(i>=satn*4*satn*3 && i<satn*4*satn*4){
				for(j=0;j<satn;j++){
					if(i<satn+(j+satn*3)*satn*4+satn+satn+satn && i>=(j+satn*3)*satn*4+satn+satn+satn)
						w1[i]=2*std2*std2;
				}
			}
		}
	}

	k=0;l=0;
	if(flag==20){
		//satn‚ÍŠÏ‘ª‰q¯”-1
		for(i=0;i<satn*4*satn*4;i++){
			if((i+satn*4+1)%(satn*4+1)==0){
				if(i<satn*4*satn*4/2){
					w1[i]=4*std1*std1;
					w1[i]=4*OmomiC[k];
					k++;
					if(k==satn)
						k=0;
				}
				if(i>=satn*4*satn*4/2){
					w1[i]=4*std2*std2;
					w1[i]=4*OmomiP[l];
					l++;
					if(l==satn)
						l=0;
				}
			}
			else if(i<satn*4*satn){
				for(j=0;j<satn;j++){
					if(i<satn+j*satn*4 && i>=j*satn*4)
						w1[i]=2*std1*std1;
				}
			}
			else if(i>=satn*4*satn && i<satn*4*satn*2){
				for(j=0;j<satn;j++){
					if(i<satn+(j+satn)*satn*4+satn && i>=(j+satn)*satn*4+satn)
						w1[i]=2*std1*std1;
				}
			}
			else if(i>=satn*4*satn*2 && i<satn*4*satn*3){
				for(j=0;j<satn;j++){
					if(i<satn+(j+satn*2)*satn*4+satn+satn && i>=(j+satn*2)*satn*4+satn+satn)
						w1[i]=2*std2*std2;
				}
			}
			else if(i>=satn*4*satn*3 && i<satn*4*satn*4){
				for(j=0;j<satn;j++){
					if(i<satn+(j+satn*3)*satn*4+satn+satn+satn && i>=(j+satn*3)*satn*4+satn+satn+satn)
						w1[i]=2*std2*std2;
				}
			}
		}
	}



	l=0;m=0;
	if(flag==22){
		for(i=0;i<satn*4*satn*4;i++){
			if((i+satn*4+1)%(satn*4+1)==0){
				if(i<satn*4*satn*4/2){
					w1[i]=4*std1*std1;
///*
					w1[i]=4*OmomiC[l];
					l++;
					if(l>=nsvg+1){
						w1[i]=4*OmomiCB[l-nsvg-1];
					}
					if(l>=nsvg+nsvc+1){
						w1[i]=4*OmomiCR[l-nsvg-nsvc-1];
					}					
					
					if(l==satn)
						l=0;
//*/
				}
				if(i>=satn*4*satn*4/2){
					w1[i]=4*std2*std2;
///*					
					w1[i]=4*OmomiP[m];
					m++;
					if(m>=nsvg+1){
						w1[i]=4*OmomiPB[m-nsvg-1];
					}
					if(m>=nsvg+nsvc+1){
						w1[i]=4*OmomiPR[m-nsvg-nsvc-1];
					}

					if(m==satn)
						m=0;
//*/
				}
			}
			else if(i<satn*4*satn){
				for(j=0;j<satn;j++){
					if(i<satn+j*satn*4 && i>=j*satn*4)
						w1[i]=2*std1*std1;
				}
			}
			else if(i>=satn*4*satn && i<satn*4*satn*2){
				for(j=0;j<satn;j++){
					if(i<satn+(j+satn)*satn*4+satn && i>=(j+satn)*satn*4+satn)
						w1[i]=2*std1*std1;
				}
			}
			else if(i>=satn*4*satn*2 && i<satn*4*satn*3){
				for(j=0;j<satn;j++){
					if(i<satn+(j+satn*2)*satn*4+satn+satn && i>=(j+satn*2)*satn*4+satn+satn)
						w1[i]=2*std2*std2;
				}
			}
			else if(i>=satn*4*satn*3 && i<satn*4*satn*4){
				for(j=0;j<satn;j++){
					if(i<satn+(j+satn*3)*satn*4+satn+satn+satn && i>=(j+satn*3)*satn*4+satn+satn+satn)
						w1[i]=2*std2*std2;
				}
			}
		}
	}


	l=0;m=0;
	if(flag==33){
		for(i=0;i<satn*4*satn*4;i++){
			if((i+satn*4+1)%(satn*4+1)==0){
				if(i<satn*4*satn*4/2){
					w1[i]=4*std1*std1;
///*
					w1[i]=4*OmomiC[l];
					l++;
					if(l>=nsvg+1){
						w1[i]=4*OmomiCB[l-nsvg-1];
					}
					if(l>=nsvg+nsvc+1){
						w1[i]=4*OmomiCR[l-nsvg-nsvc-1];
					}					
					if(l>=nsvg+nsvc+nsvb+1){
						w1[i]=4*OmomiCRR[l-nsvg-nsvc-nsvb-1];
					}	

					if(l==satn)
						l=0;
//*/
				}
				if(i>=satn*4*satn*4/2){
					w1[i]=4*std2*std2;
///*					
					w1[i]=4*OmomiP[m];
					m++;
					if(m>=nsvg+1){
						w1[i]=4*OmomiPB[m-nsvg-1];
					}
					if(m>=nsvg+nsvc+1){
						w1[i]=4*OmomiPR[m-nsvg-nsvc-1];
					}
					if(m>=nsvg+nsvc+nsvb+1){
						w1[i]=4*OmomiPRR[m-nsvg-nsvc-nsvb-1];
					}

					if(m==satn)
						m=0;
//*/
				}
			}
			else if(i<satn*4*satn){
				for(j=0;j<satn;j++){
					if(i<satn+j*satn*4 && i>=j*satn*4)
						w1[i]=2*std1*std1;
				}
			}
			else if(i>=satn*4*satn && i<satn*4*satn*2){
				for(j=0;j<satn;j++){
					if(i<satn+(j+satn)*satn*4+satn && i>=(j+satn)*satn*4+satn)
						w1[i]=2*std1*std1;
				}
			}
			else if(i>=satn*4*satn*2 && i<satn*4*satn*3){
				for(j=0;j<satn;j++){
					if(i<satn+(j+satn*2)*satn*4+satn+satn && i>=(j+satn*2)*satn*4+satn+satn)
						w1[i]=2*std2*std2;
				}
			}
			else if(i>=satn*4*satn*3 && i<satn*4*satn*4){
				for(j=0;j<satn;j++){
					if(i<satn+(j+satn*3)*satn*4+satn+satn+satn && i>=(j+satn*3)*satn*4+satn+satn+satn)
						w1[i]=2*std2*std2;
				}
			}
		}
	}


/*	
	for(i=0;i<satn*4;i++){
		for(j=0;j<satn*4;j++){
			fprintf(fp[9],"%f,",w1[j+i*satn*4]);
			if((j+1)%(satn*4)==0)
				fprintf(fp[9],"\n");
		}
	}
	exit(0);
*/

	//‹ts—ñ‚ğ‹‚ß‚é•”•ª
	int number;
	if(flag==2)
		number = satn*4;
	else if(flag==1)
		number = satn*2;
	else if(flag==3)
		number = satn*6;
	else if(flag==22)
		number = satn*4;
	else if(flag==33)
		number = satn*4;
	else if(flag==20)
		number = satn*4;
/*	
	for(i=0;i<satn*4;i++){
		for(j=0;j<satn*4;j++){
			fprintf(fp[9],"%f,",w1[j+i*satn*4]);
			if((j+1)%(satn*4)==0)
				fprintf(fp[5],"\n");
		}
	}

	exit(0);
*/


	minv(w1,eps,number);


/*	
	for(i=0;i<satn*4;i++){
		for(j=0;j<satn*4;j++){
			fprintf(fp[5],"%f,",w1[j+i*satn*4]);
			if((j+1)%(satn*4)==0)
				fprintf(fp[5],"\n");
		}
	}

	exit(0);
*/
}
