//////////////////////////////////////////////////////////////////////
//
// GPS/QZSS F9P
//
//////////////////////////////////////////////////////////////////////

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "global_extern.h"

using namespace std;

const double cs     = 299792458.0;//Œõ‘¬
const double F1     = 1575420000.0;
const double F2     = 1227600000.0;
const double myu    = 3.986005e+14;//
const double omegae = 7.2921151467e-5;//
const double pi     = 3.1415926535898;//‰~Žü—¦
const double f      = -4.442807633e-10;//

int seisu(double input);
void bias(double,double,double);
double calc_dop(int,int,int[]);
void select_sat(int[],int);
void w_inverse_float2(double[],int,int,int,int,int,int);
int minv(double[],double,int);
double calc_dop2(int,int[],int);
extern "C" {int lambda(int, int, const double *, const double *, double *, double *); }

//int lambda(int, int, const double[], const double[], double[], double[]);

void calc_rtk_GQ_F9P(int rcvn)
{
	int i,j,k,ndim;
	int prn;
	int max_prn_gps,max_prn=1;
	int kaisu=0;
	int novalid=0,miss=0;
	int l=0;
	int flag=2,gpsnum=0;

	double n_init[PRN*2]={0};
	double lrl[PRN],lrk[PRN][PRN],lrln[PRN],lrkn[PRN][PRN];
	double lrlnb[PRN];
	double dx,dy,dz;
	double y[PRN]={0};
	double eps;
	double g[PRN*PRN]={0},sg[PRN*PRN]={0},gtw[PRN*PRN]={0},g2[PRN*PRN]={0},g3[PRN*PRN]={0};
	double gt[PRN*PRN]={0},w[PRN*PRN]={0},gtwg[PRN*PRN]={0},inv_gtwg[PRN*PRN]={0},gtwg_gtw[PRN*PRN]={0};
	double q[PRN*PRN]={0};
	double ref_diff[PRN]={0},rov_diff[PRN]={0};
	double ref_diff_l2[PRN]={0},rov_diff_l2[PRN]={0};
	double ref_diff_l1[PRN]={0},rov_diff_l1[PRN]={0};
	double max_ele=0;
	double dd_code[PRN],dd_carrier[PRN];
	double dd_code_l2[PRN],dd_carrier_l2[PRN];
	double n1[PRN]={0},n2[PRN]={0};
	double n[PRN]={0};
	double dx_float,dy_float,dz_float;
	double omomi[PRN]={0},distance;

	static double lambda_pos[3]={0};
	static double s_factor=0.5,pre_ratio=0;
	double threshold_ratio = Ratio_limit;

	int nsv = SATn[rcvn]-1;//TOTAL
	int nsvg;//GPS

	dx = 0;
	dy = 0;
	dz = 0;

	for(i=0;i<PRN;i++){
		lrl[i]=0,lrln[i]=0,lrlnb[i]=0;
		for(j=0;j<PRN;j++){
			lrk[i][j]=0,lrkn[i][j]=0;
		}
	}

	//reference satellite setting
	max_ele=0;
	for(i=0;i<SATn[rcvn];i++){
		prn = SVn[rcvn][i];
		
		if(prn<=39)
			gpsnum++;

		if(max_ele < Elevation[rcvn][prn] && prn<=39){
			max_ele = Elevation[rcvn][prn];
			max_prn_gps = prn;
		}
	}

	nsvg = gpsnum-1;

	for(i=0;i<PRN*PRN;i++)
		w[i]=0;

	j=0;
	//‰ÂŽ‹‰q¯‚Ì2dˆÊ‘Š·‚ÌŽZo
	for(i=0;i<SATn[rcvn];i++){
		prn = SVn[rcvn][i];
		if(prn != max_prn_gps && prn<=39){
			//”À‘—”g‚Ì•„†‚É’ˆÓ
			dd_carrier[prn] = cs/F1*((Cp1[0][prn] - Cp1[1][prn])
				-(Cp1[0][max_prn_gps] - Cp1[1][max_prn_gps]));
			dd_carrier_l2[prn] = cs/F2*((Cp2[0][prn] - Cp2[1][prn])
				-(Cp2[0][max_prn_gps] - Cp2[1][max_prn_gps]));

			dd_code[prn] = (Pr1[0][prn] - Pr1[1][prn])
				-(Pr1[0][max_prn_gps] - Pr1[1][max_prn_gps]);
			dd_code_l2[prn] = (Pr2[0][prn] - Pr2[1][prn])
				-(Pr2[0][max_prn_gps] - Pr2[1][max_prn_gps]);

			n1[prn] = dd_carrier[prn]*F1/cs - dd_code[prn]*F1/cs;
			n2[prn] = dd_carrier_l2[prn]*F2/cs - dd_code_l2[prn]*F2/cs;

			OmomiP[j] = Code_noise+Code_noise/sin(Elevation[rcvn][prn]*pi/180.0);
			OmomiC[j] = Carrier_noise+Carrier_noise/sin(Elevation[rcvn][prn]*pi/180.0);
//			OmomiP[j] = Code_noise*1;
//			OmomiC[j] = Carrier_noise*1;
			OmomiP[j] = OmomiP[j]*OmomiP[j];
			OmomiC[j] = OmomiC[j]*OmomiC[j];
			j++;

		}//Å‘å‹ÂŠp‰q¯‚ð‚Í‚¸‚·
	}//‰q¯”

	flag=20;
	w_inverse_float2(w,SATn[rcvn]-1,flag,nsvg,0,0,0);


	j=0;//gps l1 gpsl2 (carrier) + gps l1 gpsl2 (code)
	for(i=0;i<SATn[rcvn];i++){
		prn = SVn[rcvn][i];
		if(prn != max_prn_gps && prn<=39)
			y[j++] = dd_carrier[prn];
	}
	for(i=0;i<SATn[rcvn];i++){
		prn = SVn[rcvn][i];
		if(prn != max_prn_gps && prn<=39)
			y[j++] = dd_carrier_l2[prn];
	}
	for(i=0;i<SATn[rcvn];i++){
		prn = SVn[rcvn][i];
		if(prn != max_prn_gps && prn<=39)
			y[j++] = dd_code[prn];
	}
	for(i=0;i<SATn[rcvn];i++){
		prn = SVn[rcvn][i];
		if(prn != max_prn_gps && prn<=39)
			y[j++] = dd_code_l2[prn];
	}
	distance=0;

////////////////////////////////////////////////////////////////////////////////////////////
	lrl[0] = SVx[0][max_prn_gps]-Ref_pos[0][1];
	lrl[1] = SVy[0][max_prn_gps]-Ref_pos[1][1];
	lrl[2] = SVz[0][max_prn_gps]-Ref_pos[2][1];
	lrln[0] = lrl[0]/sqrt(lrl[0]*lrl[0]+lrl[1]*lrl[1]+lrl[2]*lrl[2]);
	lrln[1] = lrl[1]/sqrt(lrl[0]*lrl[0]+lrl[1]*lrl[1]+lrl[2]*lrl[2]);
	lrln[2] = lrl[2]/sqrt(lrl[0]*lrl[0]+lrl[1]*lrl[1]+lrl[2]*lrl[2]);
	j=0;
	for(i=0;i<SATn[rcvn];i++){
		prn = SVn[rcvn][i];
		if(prn!=max_prn_gps && prn<=39){
			lrk[0][j] = SVx[0][prn]-Ref_pos[0][1];
			lrk[1][j] = SVy[0][prn]-Ref_pos[1][1];
			lrk[2][j] = SVz[0][prn]-Ref_pos[2][1];
			lrkn[0][j] = lrk[0][j]/sqrt(lrk[0][j]*lrk[0][j]+lrk[1][j]*lrk[1][j]+lrk[2][j]*lrk[2][j]);
			lrkn[1][j] = lrk[1][j]/sqrt(lrk[0][j]*lrk[0][j]+lrk[1][j]*lrk[1][j]+lrk[2][j]*lrk[2][j]);
			lrkn[2][j] = lrk[2][j]/sqrt(lrk[0][j]*lrk[0][j]+lrk[1][j]*lrk[1][j]+lrk[2][j]*lrk[2][j]);
			j++;
		}
	}
////////////////////////////////////////////////////////////////////////////////////////////


	int iter;
	for(iter=1;iter<=4;iter++){

		//Gs—ñ‚Ì¶¬@ã‚©‚çL1”À‘—”g, L2”À‘—”g, L1ƒR[ƒh, L2ƒR[ƒh
		for(i=0;i<nsv*2+3;i++){
			for(j=0;j<nsv*4;j++){
					g[i+j*(nsv*2+3)]=-1.0*(lrkn[i][j%nsv]-lrln[i]);
			}
		}
		for(i=0;i<(nsv*2+3)*nsv;i++){
			if((i-3)%(nsv*2+3+1)==0)
				g[i] = cs/F1;
		}
		for(i=(nsv*2+3)*nsv;i<(nsv*2+3)*nsv*2;i++){
			if((i-3)%(nsv*2+3+1)==0)
				g[i] = cs/F2;
		}

		//G‚Ì“]’us—ñ‚Ì¶¬
		for(i=0;i<nsv*4;i++){
			for(j=0;j<nsv*2+3;j++){
				gt[i+j*nsv*4]=g[i*(nsv*2+3)+j];
			}
		}

		for(i=0;i<48*48;i++){
			gtw[i]=0;
			gtwg[i]=0;
			gtwg_gtw[i]=0;
		}

		//G(t)*W
		for(i=0;i<nsv*2+3;i++){
			for(j=0;j<nsv*4;j++){
				for(k=0;k<nsv*4;k++){
					gtw[i*nsv*4+j]=gtw[i*nsv*4+j]+gt[k+i*nsv*4]*w[k*nsv*4+j];
				}
			}
		}

		//G(t)*W*G
		for(i=0;i<nsv*2+3;i++){
			for(j=0;j<nsv*2+3;j++){
				for(k=0;k<nsv*4;k++){
					gtwg[i*(nsv*2+3)+j]=gtwg[i*(nsv*2+3)+j]+gtw[k+i*nsv*4]*g[k*(nsv*2+3)+j];
				}
			}
		}

		//G(t)*W*G‚Ì‹ts—ñ‚ð‹‚ß‚é
		eps=1.0e-25;
		minv(gtwg,eps,nsv*2+3);

		//inv((G(t)*W*G))*G(t)W
		for(i=0;i<nsv*2+3;i++){
			for(j=0;j<nsv*4;j++){
				for(k=0;k<=nsv*2+3;k++){
					gtwg_gtw[i*nsv*4+j]=gtwg_gtw[i*nsv*4+j]+gtwg[k+i*(nsv*2+3)]*gtw[k*nsv*4+j];
				}
			}
		}

		if(iter==1){
			j=0;
			for(i=0;i<SATn[rcvn];i++){//GPS+QZS
				prn = SVn[rcvn][i];
				if(prn!=max_prn_gps && prn<=39){
					ref_diff[j] = sqrt((SVx[1][prn]-Ref_pos[0][1])*(SVx[1][prn]-Ref_pos[0][1])
						+(SVy[1][prn]-Ref_pos[1][1])*(SVy[1][prn]-Ref_pos[1][1])
						+(SVz[1][prn]-Ref_pos[2][1])*(SVz[1][prn]-Ref_pos[2][1]))
						-sqrt((SVx[1][max_prn_gps]-Ref_pos[0][1])*(SVx[1][max_prn_gps]-Ref_pos[0][1])
						+(SVy[1][max_prn_gps]-Ref_pos[1][1])*(SVy[1][max_prn_gps]-Ref_pos[1][1])
						+(SVz[1][max_prn_gps]-Ref_pos[2][1])*(SVz[1][max_prn_gps]-Ref_pos[2][1]));
					j++;
				}
			}
			for(i=0;i<j*4;i++){
				ref_diff[i] = ref_diff[i%j];//GPS+QZS+BEIDOU
			}
		}

		j=0;
		for(i=0;i<SATn[rcvn];i++){
			prn = SVn[rcvn][i];
			if(prn!=max_prn_gps && prn<=39){
				rov_diff[j] = sqrt(lrk[0][j]*lrk[0][j]+lrk[1][j]*lrk[1][j]+lrk[2][j]*lrk[2][j])
					- sqrt(lrl[0]*lrl[0]+lrl[1]*lrl[1]+lrl[2]*lrl[2]);
				j++;
			}
		}
		for(i=0;i<j*4;i++){
			rov_diff[i] = rov_diff[i%j];
		}

			for(i=0;i<nsv*4;i++){
				dx = dx+gtwg_gtw[i]*(y[i]+ref_diff[i]-rov_diff[i]);
				dy = dy+gtwg_gtw[i+nsv*4]*(y[i]+ref_diff[i]-rov_diff[i]);
				dz = dz+gtwg_gtw[i+nsv*4*2]*(y[i]+ref_diff[i]-rov_diff[i]);

				if(iter==4){
					for(j=3;j<nsv*2+3;j++)
						n[j-3] = n[j-3]+gtwg_gtw[i+nsv*4*j]*(y[i]+ref_diff[i]-rov_diff[i]);
				}
			}

/////////////////////////////////////////////////////////////////////////////////////////////
		lrl[0] = SVx[0][max_prn_gps]-(Ref_pos[0][1]+dx);
		lrl[1] = SVy[0][max_prn_gps]-(Ref_pos[1][1]+dy);
		lrl[2] = SVz[0][max_prn_gps]-(Ref_pos[2][1]+dz);
		lrln[0] = lrl[0]/sqrt(lrl[0]*lrl[0]+lrl[1]*lrl[1]+lrl[2]*lrl[2]);
		lrln[1] = lrl[1]/sqrt(lrl[0]*lrl[0]+lrl[1]*lrl[1]+lrl[2]*lrl[2]);
		lrln[2] = lrl[2]/sqrt(lrl[0]*lrl[0]+lrl[1]*lrl[1]+lrl[2]*lrl[2]);
		j=0;
		for(i=0;i<SATn[rcvn];i++){
			prn = SVn[rcvn][i];
			if(prn!=max_prn_gps && prn<=39){
				lrk[0][j] = SVx[0][prn]-(Ref_pos[0][1]+dx);
				lrk[1][j] = SVy[0][prn]-(Ref_pos[1][1]+dy);
				lrk[2][j] = SVz[0][prn]-(Ref_pos[2][1]+dz);
				lrkn[0][j] = lrk[0][j]/sqrt(lrk[0][j]*lrk[0][j]+lrk[1][j]*lrk[1][j]+lrk[2][j]*lrk[2][j]);
				lrkn[1][j] = lrk[1][j]/sqrt(lrk[0][j]*lrk[0][j]+lrk[1][j]*lrk[1][j]+lrk[2][j]*lrk[2][j]);
				lrkn[2][j] = lrk[2][j]/sqrt(lrk[0][j]*lrk[0][j]+lrk[1][j]*lrk[1][j]+lrk[2][j]*lrk[2][j]);
				j++;
			}
		}
/////////////////////////////////////////////////////////////////////////////////////////////

	}//iter

	//float‰ð‚Ì’Šo
	dx_float = dx;
	dy_float = dy;
	dz_float = dz;

	//LAMBDA–@‚Ì“K—p
	double F[PRN]={0},s[PRN]={0};

	j=0;
	for(i=0;i<=(nsv*2+3)*(nsv*2+3);i++){
		if(i>=(nsv*2+3)*3){
			q[j] = gtwg[i];
			j++;
			if(i%(nsv*2+3)==0)
				j--;
			if((i-1)%(nsv*2+3)==0)
				j--;
			if((i-2)%(nsv*2+3)==0)
				j--;
		}
	}

	k=0;
	ndim=nsv*2;
	double det=1.0,buf;
	double aa[PRN][PRN]={0};
	for(i=0;i<ndim;i++){
		for(j=0;j<nsv*2;j++){
			aa[i][j] = q[k++];
		}
	}
	//ŽOŠps—ñ‚ðì¬
	for(i=0;i<ndim;i++){
		for(j=0;j<ndim;j++){
			if(i<j){
				buf=aa[j][i]/aa[i][i];
				for(k=0;k<ndim;k++){
					aa[j][k]-=aa[i][k]*buf;
				}
			}
		}
	}
	//‘ÎŠp•”•ª‚ÌÏ
	for(i=0;i<ndim;i++){
		det*=aa[i][i];
	}
	det = pow(det,1.0/(nsv*2*2));//ADOP


	//LAMBDA–@‚ÌŒÄ‚Ño‚µ
	lambda(nsv*2,2,n,q,F,s);

///////////////////////////////////////////////////////////////////////////////////////////////////////
	//LAMBDA–@‚Å‹‚ß‚½ƒAƒ“ƒrƒMƒ…ƒCƒeƒC‚ð—˜—p‚µ‚Ä‘ªˆÊŒvŽZ‚ðs‚¤
///////////////////////////////////////////////////////////////////////////////////////////////////////


	//Še‰q¯‚Ì‹ÂŠp‚É‚æ‚éd‚Ý‚Ì¶¬
	j=0;
	for(i=0;i<SATn[rcvn];i++){
		prn = SVn[rcvn][i];
		if(prn != max_prn_gps && prn<=39)
			omomi[j++] = 1.0*1.0*sqrt(90.0/Elevation[rcvn][prn]);
	}

	dx=0;dy=0;dz=0;

	for(i=nsv*2;i<nsv*4;i++)//Žc·Å¬‚ÌŒó•âiƒAƒ“ƒrƒMƒ…ƒCƒeƒB‚»‚Ì‚à‚Ìj‚Ì‚ÝŽc‚·
		F[i]=0;

/////////////////////////////////////////////////////////////////////////////////////////////
	lrl[0] = SVx[0][max_prn_gps]-(Ref_pos[0][1]+dx);
	lrl[1] = SVy[0][max_prn_gps]-(Ref_pos[1][1]+dy);
	lrl[2] = SVz[0][max_prn_gps]-(Ref_pos[2][1]+dz);
	lrln[0] = lrl[0]/sqrt(lrl[0]*lrl[0]+lrl[1]*lrl[1]+lrl[2]*lrl[2]);
	lrln[1] = lrl[1]/sqrt(lrl[0]*lrl[0]+lrl[1]*lrl[1]+lrl[2]*lrl[2]);
	lrln[2] = lrl[2]/sqrt(lrl[0]*lrl[0]+lrl[1]*lrl[1]+lrl[2]*lrl[2]);
	j=0;
	for(i=0;i<SATn[rcvn];i++){
		prn = SVn[rcvn][i];
		if(prn!=max_prn_gps && prn<=39){
			lrk[0][j] = SVx[0][prn]-(Ref_pos[0][1]+dx);
			lrk[1][j] = SVy[0][prn]-(Ref_pos[1][1]+dy);
			lrk[2][j] = SVz[0][prn]-(Ref_pos[2][1]+dz);
			lrkn[0][j] = lrk[0][j]/sqrt(lrk[0][j]*lrk[0][j]+lrk[1][j]*lrk[1][j]+lrk[2][j]*lrk[2][j]);
			lrkn[1][j] = lrk[1][j]/sqrt(lrk[0][j]*lrk[0][j]+lrk[1][j]*lrk[1][j]+lrk[2][j]*lrk[2][j]);
			lrkn[2][j] = lrk[2][j]/sqrt(lrk[0][j]*lrk[0][j]+lrk[1][j]*lrk[1][j]+lrk[2][j]*lrk[2][j]);
			j++;
		}
	}
/////////////////////////////////////////////////////////////////////////////////////////////


	for(iter=1;iter<=3;iter++){
		//dd_carrier‚ðm‚©‚ç”g”‚É‚·‚×‚«H
		if(iter==1){
			for(i=0;i<nsv;i++){
				y[i] = y[i]*F1/cs;
			}
		}

		//Gs—ñ‚Ì¶¬
		for(i=0;i<=2;i++){
			for(j=0;j<nsv;j++){
				g[i+j*3]=-1.0*(lrkn[i][j]-lrln[i])*F1/cs;
			}
		}

		//G‚Ì“]’us—ñ‚Ì¶¬
		for(i=0;i<nsv;i++){
			for(j=0;j<3;j++){
				gt[i+j*nsv]=g[i*3+j];
			}
		}

		//G‚Ì“]’us—ñ‚ÆW‚ÌŠ|‚¯ŽZ@‚±‚êˆÈ~‚Ígt‚Ígtw‚Æ‚È‚Á‚Ä‚¢‚é
		for(i=0;i<nsv;i++){
			for(j=0;j<3;j++){
				gt[i+j*nsv]=gt[i+j*nsv]*omomi[i];
			}
		}

		for(i=0;i<=8;i++)
			g3[i]=0;

		//G‚Ì“]’us—ñ‚ÆGs—ñ‚ðŠ|‚¯‚é
		for(i=0;i<nsv;i++){
			g3[0]=g3[0]+gt[i]*g[i*3];
			g3[1]=g3[1]+gt[i]*g[i*3+1];
			g3[2]=g3[2]+gt[i]*g[i*3+2];
			g3[3]=g3[3]+gt[i+nsv]*g[i*3];
			g3[4]=g3[4]+gt[i+nsv]*g[i*3+1];
			g3[5]=g3[5]+gt[i+nsv]*g[i*3+2];
			g3[6]=g3[6]+gt[i+nsv*2]*g[i*3];
			g3[7]=g3[7]+gt[i+nsv*2]*g[i*3+1];
			g3[8]=g3[8]+gt[i+nsv*2]*g[i*3+2];
		}

		//G3s—ñ‚Ì‹ts—ñ‚ð‹‚ß‚é
		//‹ts—ñŒvŽZ//
		eps = 1.0e-25;
		minv(g3,eps,3);
	
		for(i=0;i<nsv*3;i++)
			g2[i]=0;
		//‹ts—ñ‚ÉG‚Ì“]’us—ñ‚ðŠ|‚¯‚é
		for(i=0;i<nsv;i++){
			for(j=0;j<3;j++){
				g2[i]      =g2[i]+g3[j]*gt[i+nsv*j];
				g2[nsv+i]  =g2[nsv+i]+g3[j+3]*gt[i+nsv*j];
				g2[nsv*2+i]=g2[nsv*2+i]+g3[j+6]*gt[i+nsv*j];
			}
		}

		if(iter==1){
			j=0;
			for(i=0;i<SATn[rcvn];i++){
				prn = SVn[rcvn][i];
				if(prn!=max_prn_gps && prn<=39){
					ref_diff[j] = sqrt((SVx[1][prn]-Ref_pos[0][1])*(SVx[1][prn]-Ref_pos[0][1])
						+(SVy[1][prn]-Ref_pos[1][1])*(SVy[1][prn]-Ref_pos[1][1])
						+(SVz[1][prn]-Ref_pos[2][1])*(SVz[1][prn]-Ref_pos[2][1]))
						-sqrt((SVx[1][max_prn_gps]-Ref_pos[0][1])*(SVx[1][max_prn_gps]-Ref_pos[0][1])
						+(SVy[1][max_prn_gps]-Ref_pos[1][1])*(SVy[1][max_prn_gps]-Ref_pos[1][1])
						+(SVz[1][max_prn_gps]-Ref_pos[2][1])*(SVz[1][max_prn_gps]-Ref_pos[2][1]));
					ref_diff[j] = ref_diff[j]*F1/cs;
					j++;
				}
			}
		}

		j=0;
		for(i=0;i<SATn[rcvn];i++){
			prn = SVn[rcvn][i];
			if(prn!=max_prn_gps && prn<=39){
				rov_diff[j] = sqrt(lrk[0][j]*lrk[0][j]+lrk[1][j]*lrk[1][j]+lrk[2][j]*lrk[2][j])
					- sqrt(lrl[0]*lrl[0]+lrl[1]*lrl[1]+lrl[2]*lrl[2]);
				rov_diff[j] = rov_diff[j]*F1/cs;
				j++;
			}
		}

		//dx,dy,dz‚ÌŒvŽZ
		for(i=0;i<nsv;i++){
			dx = dx+g2[i]*(y[i]+ref_diff[i]-rov_diff[i]-F[i]);
			dy = dy+g2[i+nsv]*(y[i]+ref_diff[i]-rov_diff[i]-F[i]);
			dz = dz+g2[i+nsv*2]*(y[i]+ref_diff[i]-rov_diff[i]-F[i]);
		}

		//‰q¯‚Ì•ûŒü‚Ì’PˆÊƒxƒNƒgƒ‹‚ðŽZo
		//Šî€‹Ç‚ÆŠî€‰q¯‚Ì’PˆÊƒxƒNƒgƒ‹
		lrl[0] = SVx[0][max_prn_gps]-(Ref_pos[0][1]+dx);
		lrl[1] = SVy[0][max_prn_gps]-(Ref_pos[1][1]+dy);
		lrl[2] = SVz[0][max_prn_gps]-(Ref_pos[2][1]+dz);
		lrln[0] = lrl[0]/sqrt(lrl[0]*lrl[0]+lrl[1]*lrl[1]+lrl[2]*lrl[2]);
		lrln[1] = lrl[1]/sqrt(lrl[0]*lrl[0]+lrl[1]*lrl[1]+lrl[2]*lrl[2]);
		lrln[2] = lrl[2]/sqrt(lrl[0]*lrl[0]+lrl[1]*lrl[1]+lrl[2]*lrl[2]);

		//ˆÚ“®‹Ç‚Æ‘I‘ð‰q¯‚Ì’PˆÊƒxƒNƒgƒ‹
		//‰q¯”•ª‚Ü‚í‚·
		//Å‘å‹ÂŠp‚Ì‰q¯‚Í‚·‚Å‚É‘I‘ð‚³‚ê‚Ä‚¢‚é

		j=0;
		for(i=0;i<SATn[rcvn];i++){
			prn = SVn[rcvn][i];
			if(prn!=max_prn_gps && prn<=39){
				lrk[0][j] = SVx[0][prn]-(Ref_pos[0][1]+dx);
				lrk[1][j] = SVy[0][prn]-(Ref_pos[1][1]+dy);
				lrk[2][j] = SVz[0][prn]-(Ref_pos[2][1]+dz);
				lrkn[0][j] = lrk[0][j]/sqrt(lrk[0][j]*lrk[0][j]+lrk[1][j]*lrk[1][j]+lrk[2][j]*lrk[2][j]);
				lrkn[1][j] = lrk[1][j]/sqrt(lrk[0][j]*lrk[0][j]+lrk[1][j]*lrk[1][j]+lrk[2][j]*lrk[2][j]);
				lrkn[2][j] = lrk[2][j]/sqrt(lrk[0][j]*lrk[0][j]+lrk[1][j]*lrk[1][j]+lrk[2][j]*lrk[2][j]);
				j++;
			}
		}
	}//iter

	//x,y,z -> lat,lon,hgt
	double ecef[3]={0},t[3]={0};
	double a,f,b,p,e,oldt,nn;
	double lat_diff,lon_diff,hgt_diff;

	ecef[0]=dx+Ref_pos[0][1];
	ecef[1]=dy+Ref_pos[1][1];
	ecef[2]=dz+Ref_pos[2][1];
	lambda_pos[0]=ecef[0];lambda_pos[1]=ecef[1];lambda_pos[2]=ecef[2];

	if(s[1]/s[0]>=threshold_ratio){
		Fix_pos[0] = lambda_pos[0];Fix_pos[1] = lambda_pos[1];Fix_pos[2] = lambda_pos[2];
	}

	a = 6378137;
	f = 1/298.257223563;
	b = a*(1-f);
	p = sqrt(ecef[0]*ecef[0]+ecef[1]*ecef[1]);
	e = f*(2-f);
	t[0] = atan(ecef[2]/((1-e)*p));
	oldt = 1;
	while(fabs(t[2]-oldt) > 0.000001)
	{
		oldt = t[2];
		nn = a/sqrt(1-(e*sin(t[0])*sin(t[0])));
		t[2] = p/cos(t[0]) - nn;
		t[0] = atan((ecef[2]/p)/((1-e*nn/(nn+t[2]))));
	}
	t[1] = 2*atan(ecef[1]/(ecef[0]+sqrt(ecef[0]*ecef[0]+ecef[1]*ecef[1])));
	t[0] = t[0]*180/pi;
	t[1] = t[1]*180/pi;

	lat_diff = (t[0]-POSrcvlat[1])*111319.49;
	lon_diff = (t[1]-POSrcvlon[1])*cos(t[0]*pi/180.0)*111319.49;
	hgt_diff = t[2] - POSrcvhgt[1];

	s_factor = s[1]/s[0];

	//ADOP‚Å‚Ìƒ`ƒFƒbƒN
	if(det>=1.00)
		s_factor = 1.0;

	double aaa,bbb,xxx,yyy;
	aaa=POSrcvlat[0];
	bbb=POSrcvlon[0];
	xxx = (t[0]-aaa)*111319.49;
	yyy = (t[1]-bbb)*cos(t[0]*pi/180.0)*111319.49;

	Ratio = s_factor;


	//nmea
	double nmea_mmm[2],nmea_ssss,nmea_hight,nmea_geoid;
	int stat_flag,nmea_hh,nmea_mm,nmea_ddd[2];
	////////////////////
	stat_flag=4;
	nmea_geoid=36.209;
	////////////////////
	nmea_hh=int((DGPSTIME-LeapSecond)/60/60);
	nmea_mm=int((DGPSTIME-LeapSecond-nmea_hh*60*60)/60);
	nmea_ssss=DGPSTIME-LeapSecond-nmea_hh*60*60-nmea_mm*60;
	for(i=0;i<2;i++){
		nmea_ddd[i]=int(t[i]);
		nmea_mmm[i]=(t[i]-int(t[i]))*60;
	}
	nmea_hight=t[2]-nmea_geoid;


	if(s_factor>=threshold_ratio){
		Sol_flag[2]++;

		fprintf(fp[5],"%f,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
			DGPSTIME,SATn[rcvn],yyy,xxx,t[0],t[1],t[2],s[1]/s[0],HDOP,VDOP,det);

		fprintf(fp[7],"%f,$GPGGA,%02d%02d%05.2f,%02d%012.9f,N,%03d%012.9f,E,%1d,%02d,%.1f,%.3f,M,%.3f,M,%.1f,*40\n",DGPSTIME
			,nmea_hh,nmea_mm,nmea_ssss,nmea_ddd[0],nmea_mmm[0],nmea_ddd[1],nmea_mmm[1],stat_flag,SATn[rcvn],HDOP,nmea_hight,nmea_geoid,DGPSTIME-GPSTIME);
	}

/////////////////////////////////////////////////////////////////////


}