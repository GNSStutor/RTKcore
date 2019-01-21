//////////////////////////////////////////////////////////////////////
//
// GPS/QZSS/BeiDou/Galileo F9P
//
//////////////////////////////////////////////////////////////////////

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "global_extern.h"

using namespace std;

const double cs     = 299792458.0;//����
const double F1     = 1575420000.0;
const double F2     = 1227600000.0;
const double myu    = 3.986005e+14;//
const double omegae = 7.2921151467e-5;//
const double pi     = 3.1415926535898;//�~����
const double f      = -4.442807633e-10;//

const double B1     = 1561098000.0;
const double B2     = 1207140000.0;
const double E1     = 1575420000.0;
const double E5B    = 1207140000.0;

int seisu(double input);
void bias(double,double,double);
double calc_dop(int,int,int[]);
void select_sat(int[],int);
void w_inverse_float2(double[],int,int,int,int,int,int);
int minv(double[],double,int);
double calc_dop2(int,int[],int);
extern "C" {int lambda(int, int, const double *, const double *, double *, double *); }

//int lambda(int, int, const double[], const double[], double[], double[]);

void calc_rtk_GQEB_F9P(int rcvn)
{
	int i,j,k,ndim;
	int prn;
	int max_prn_gps,max_prn_bei,max_prn_c,max_prn=1;
	int kaisu=0;
	int novalid=0,miss=0;
	int l=0;
	int flag=2,gpsnum=0,glonum=0,beinum=0,svn[PRN];

	double n_init[PRN*2]={0};
	double lrl[PRN],lrk[PRN][PRN],lrln[PRN],lrkn[PRN][PRN];
	double lrlb[PRN],lrlnb[PRN];
	double lrlc[PRN],lrlnc[PRN];
	double dx,dy,dz;
	double y[PRN]={0};
	double eps;
	double g[PRN*PRN]={0},sg[PRN*PRN]={0},gtw[PRN*PRN]={0},g2[PRN*PRN]={0},g3[PRN*PRN]={0};
	double gt[PRN*PRN]={0},w[PRN*PRN]={0},gtwg[PRN*PRN]={0},inv_gtwg[PRN*PRN]={0},gtwg_gtw[PRN*PRN]={0};
	double q[PRN*PRN]={0};
	double ref_diff[PRN]={0},rov_diff[PRN]={0},ref_diff2[PRN]={0},rov_diff2[PRN]={0};
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
	int No_dp_flag=1;

	int nsv = SATn[rcvn]-1-1-1;//TOTAL
	int nsvg;//GPS
	int nsvb;//GAL
	int nsvc;//BeiDou
	int nsvr = 0;

	dx = 0;
	dy = 0;
	dz = 0;

	for(i=0;i<PRN;i++){
		lrl[i]=0,lrln[i]=0,lrlb[i]=0,lrlnb[i]=0,lrlc[i]=0,lrlnc[i]=0;
		for(j=0;j<PRN;j++){
			lrk[i][j]=0,lrkn[i][j]=0;
		}
	}

	//��q��
	max_ele=0;//GPS+QZS
	for(i=0;i<SATn[rcvn];i++){
		prn = SVn[rcvn][i];
		if(prn<=39)
			gpsnum++;
		if(max_ele < Elevation[rcvn][prn] && prn<=39){
			max_ele = Elevation[rcvn][prn];
			max_prn_gps = prn;
		}
	}
	max_ele=0;//Galileo
	for(i=0;i<SATn[rcvn];i++){
		prn = SVn[rcvn][i];
		if(prn<=70 && prn>=41)
			glonum++;
		if(max_ele < Elevation[rcvn][prn] && prn>=41 && prn<=70){
			max_ele = Elevation[rcvn][prn];
			max_prn_bei = prn;
		}
	}
	max_ele=0;//BEIDOU
	for(i=0;i<SATn[rcvn];i++){
		prn = SVn[rcvn][i];
		if(prn<=100 && prn>=71)
			beinum++;
		if(max_ele < Elevation[rcvn][prn] && prn>=71 && prn<=100){
			max_ele = Elevation[rcvn][prn];
			max_prn_c = prn;
		}
	}

	nsvg = gpsnum-1;//GPS
	nsvc = beinum-1;//BEIDOU
	nsvb = nsv - nsvg - nsvc;//Galileo

	for(i=0;i<PRN*PRN;i++)
		w[i]=0;

	j=0;
	//���q����2�d�ʑ����̎Z�o Galileo
	for(i=0;i<SATn[rcvn];i++){
		prn = SVn[rcvn][i];
		if(prn != max_prn_bei && prn>=41 && prn<=70){
			dd_carrier[prn] = ((Cp1[0][prn] - Cp1[1][prn])
				-(Cp1[0][max_prn_bei] - Cp1[1][max_prn_bei]));
			dd_carrier_l2[prn] = ((Cp2[0][prn] - Cp2[1][prn])
				-(Cp2[0][max_prn_bei] - Cp2[1][max_prn_bei]));

			dd_code[prn] = (E1/cs)*(Pr1[0][prn] - Pr1[1][prn])
				-(E1/cs)*(Pr1[0][max_prn_bei] - Pr1[1][max_prn_bei]);
			dd_code_l2[prn] = (E5B/cs)*(Pr2[0][prn] - Pr2[1][prn])
				-(E5B/cs)*(Pr2[0][max_prn_bei] - Pr2[1][max_prn_bei]);

			OmomiPR[j] = Code_noise+1*Code_noise/sin(Elevation[rcvn][prn]*pi/180.0);
			OmomiCR[j] = Carrier_noise+1*Carrier_noise/sin(Elevation[rcvn][prn]*pi/180.0); 

//			OmomiPR[j] = Code_noise*1;
//			OmomiCR[j] = Carrier_noise*1;
			OmomiPR[j] = OmomiPR[j]*OmomiPR[j];
			OmomiCR[j] = OmomiCR[j]*OmomiCR[j];
			j++;

		}//�ő�p�q�����͂���
	}//�q����

	j=0;
	//���q����2�d�ʑ����̎Z�o BEIDOU
	for(i=0;i<SATn[rcvn];i++){
		prn = SVn[rcvn][i];
		if(prn != max_prn_c && prn>=71 && prn<=100){
			dd_carrier[prn] = ((Cp1[0][prn] - Cp1[1][prn])
				-(Cp1[0][max_prn_c] - Cp1[1][max_prn_c]));
			dd_carrier_l2[prn] = ((Cp2[0][prn] - Cp2[1][prn])
				-(Cp2[0][max_prn_c] - Cp2[1][max_prn_c]));

			dd_code[prn] = (B1/cs)*((Pr1[0][prn] - Pr1[1][prn])
				-(Pr1[0][max_prn_c] - Pr1[1][max_prn_c]));
			dd_code_l2[prn] = (B2/cs)*((Pr2[0][prn] - Pr2[1][prn])
				-(Pr2[0][max_prn_c] - Pr2[1][max_prn_c]));

			OmomiPB[j] = Code_noise+1*Code_noise/sin(Elevation[rcvn][prn]*pi/180.0);
			OmomiCB[j] = Carrier_noise+1*Carrier_noise/sin(Elevation[rcvn][prn]*pi/180.0); 
//			OmomiPB[j] = Code_noise*1;
//			OmomiCB[j] = Carrier_noise*1;
			OmomiPB[j] = OmomiPB[j]*OmomiPB[j];
			OmomiCB[j] = OmomiCB[j]*OmomiCB[j];
			j++;

		}//�ő�p�q�����͂���
	}//�q����

	j=0;
	//���q����2�d�ʑ����̎Z�o GPS
	for(i=0;i<SATn[rcvn];i++){
		prn = SVn[rcvn][i];
		if(prn != max_prn_gps && prn<=39){
			//�����g�̕����ɒ���
			dd_carrier[prn] = ((Cp1[0][prn] - Cp1[1][prn])
				-(Cp1[0][max_prn_gps] - Cp1[1][max_prn_gps]));
			dd_carrier_l2[prn] = ((Cp2[0][prn] - Cp2[1][prn])
				-(Cp2[0][max_prn_gps] - Cp2[1][max_prn_gps]));

			dd_code[prn] = F1/cs*((Pr1[0][prn] - Pr1[1][prn])
				-(Pr1[0][max_prn_gps] - Pr1[1][max_prn_gps]));
			dd_code_l2[prn] = F2/cs*((Pr2[0][prn] - Pr2[1][prn])
				-(Pr2[0][max_prn_gps] - Pr2[1][max_prn_gps]));

			OmomiP[j] = Code_noise+1*Code_noise/sin(Elevation[rcvn][prn]*pi/180.0);
			OmomiC[j] = Carrier_noise+1*Carrier_noise/sin(Elevation[rcvn][prn]*pi/180.0); 
//			OmomiP[j] = Code_noise*1;
//			OmomiC[j] = Carrier_noise*1;
			OmomiP[j] = OmomiP[j]*OmomiP[j];
			OmomiC[j] = OmomiC[j]*OmomiC[j];
			j++;

		}//�ő�p�q�����͂���
	}//�q����

	flag=22;
	w_inverse_float2(w,SATn[rcvn]-3,flag,nsvg,nsvb,nsvc,nsvr);

	j=0;//gps bei gal L1 gps bei ggal L2(carrier) + gps bei gal L1 gps bei gal L2(code)
	for(i=0;i<SATn[rcvn];i++){
		prn = SVn[rcvn][i];
		if(prn != max_prn_gps && prn<=39){
			svn[j] = prn;
			y[j++] = dd_carrier[prn];
		}
	}
	for(i=0;i<SATn[rcvn];i++){
		prn = SVn[rcvn][i];
		if(prn != max_prn_c && prn<=100 && prn>=71){
			svn[j] = prn;
			y[j++] = dd_carrier[prn];
		}
	}
	for(i=0;i<SATn[rcvn];i++){
		prn = SVn[rcvn][i];
		if(prn != max_prn_bei && prn<=70 && prn>=41){
			svn[j] = prn;
			y[j++] = dd_carrier[prn];
		}
	}
	for(i=0;i<SATn[rcvn];i++){
		prn = SVn[rcvn][i];
		if(prn != max_prn_gps && prn<=39)
			y[j++] = dd_carrier_l2[prn];
	}
	for(i=0;i<SATn[rcvn];i++){
		prn = SVn[rcvn][i];
		if(prn != max_prn_c && prn<=100 && prn>=71)
			y[j++] = dd_carrier_l2[prn];
	}
	for(i=0;i<SATn[rcvn];i++){
		prn = SVn[rcvn][i];
		if(prn != max_prn_bei && prn<=70 && prn>=41)
			y[j++] = dd_carrier_l2[prn];
	}

	for(i=0;i<SATn[rcvn];i++){
		prn = SVn[rcvn][i];
		if(prn != max_prn_gps && prn<=39)
			y[j++] = dd_code[prn];
	}
	for(i=0;i<SATn[rcvn];i++){
		prn = SVn[rcvn][i];
		if(prn != max_prn_c && prn<=100 && prn>=71)
			y[j++] = dd_code[prn];
	}
	for(i=0;i<SATn[rcvn];i++){
		prn = SVn[rcvn][i];
		if(prn != max_prn_bei && prn<=70 && prn>=41)
			y[j++] = dd_code[prn];
	}
	for(i=0;i<SATn[rcvn];i++){
		prn = SVn[rcvn][i];
		if(prn != max_prn_gps && prn<=39)
			y[j++] = dd_code_l2[prn];
	}
	for(i=0;i<SATn[rcvn];i++){
		prn = SVn[rcvn][i];
		if(prn != max_prn_c && prn<=100 && prn>=71)
			y[j++] = dd_code_l2[prn];
	}
	for(i=0;i<SATn[rcvn];i++){
		prn = SVn[rcvn][i];
		if(prn != max_prn_bei && prn<=70 && prn>=41)
			y[j++] = dd_code_l2[prn];
	}
	distance=0;



////////////////////////////////////////////////////////////////////////////////////////////
	//GPS/QZS
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

	//BEIDOU
	lrlc[0] = SVx[0][max_prn_c]-Ref_pos[0][1];
	lrlc[1] = SVy[0][max_prn_c]-Ref_pos[1][1];
	lrlc[2] = SVz[0][max_prn_c]-Ref_pos[2][1];
	lrlnc[0] = lrlc[0]/sqrt(lrlc[0]*lrlc[0]+lrlc[1]*lrlc[1]+lrlc[2]*lrlc[2]);
	lrlnc[1] = lrlc[1]/sqrt(lrlc[0]*lrlc[0]+lrlc[1]*lrlc[1]+lrlc[2]*lrlc[2]);
	lrlnc[2] = lrlc[2]/sqrt(lrlc[0]*lrlc[0]+lrlc[1]*lrlc[1]+lrlc[2]*lrlc[2]);
	for(i=0;i<SATn[rcvn];i++){
		prn = SVn[rcvn][i];
		if(prn!=max_prn_c && prn<=100 && prn>=71){
			lrk[0][j] = SVx[0][prn]-Ref_pos[0][1];
			lrk[1][j] = SVy[0][prn]-Ref_pos[1][1];
			lrk[2][j] = SVz[0][prn]-Ref_pos[2][1];
			lrkn[0][j] = lrk[0][j]/sqrt(lrk[0][j]*lrk[0][j]+lrk[1][j]*lrk[1][j]+lrk[2][j]*lrk[2][j]);
			lrkn[1][j] = lrk[1][j]/sqrt(lrk[0][j]*lrk[0][j]+lrk[1][j]*lrk[1][j]+lrk[2][j]*lrk[2][j]);
			lrkn[2][j] = lrk[2][j]/sqrt(lrk[0][j]*lrk[0][j]+lrk[1][j]*lrk[1][j]+lrk[2][j]*lrk[2][j]);
			j++;
		}
	}

	//Galileo
	lrlb[0] = SVx[0][max_prn_bei]-Ref_pos[0][1];
	lrlb[1] = SVy[0][max_prn_bei]-Ref_pos[1][1];
	lrlb[2] = SVz[0][max_prn_bei]-Ref_pos[2][1];
	lrlnb[0] = lrlb[0]/sqrt(lrlb[0]*lrlb[0]+lrlb[1]*lrlb[1]+lrlb[2]*lrlb[2]);
	lrlnb[1] = lrlb[1]/sqrt(lrlb[0]*lrlb[0]+lrlb[1]*lrlb[1]+lrlb[2]*lrlb[2]);
	lrlnb[2] = lrlb[2]/sqrt(lrlb[0]*lrlb[0]+lrlb[1]*lrlb[1]+lrlb[2]*lrlb[2]);
	for(i=0;i<SATn[rcvn];i++){
		prn = SVn[rcvn][i];
		if(prn!=max_prn_bei && prn<=70 && prn>=41){
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
		//G�s��̐����@�ォ��L1�����g(GPS+BEI+GAL), L2�����g(GPS+BEI+GAL), L1�R�[�h(GPS+BEI+GAL), L2�R�[�h(GPS+BEI+GAL)
		for(i=0;i<nsv*2+3;i++){
			for(j=0;j<nsv*4;j++){
				if(j<nsvg)
					g[i+j*(nsv*2+3)]=-1.0*(lrkn[i][j%nsv]-lrln[i])*F1/cs;
				else if(j>=nsvg && j<nsvg+nsvc)
					g[i+j*(nsv*2+3)]=-1.0*(lrkn[i][j%nsv]*B1/cs-lrlnc[i]*B1/cs);
				else if(j>=nsvg+nsvc && j<nsvg+nsvc+nsvb)
					g[i+j*(nsv*2+3)]=-1.0*(lrkn[i][j%nsv]*E1/cs-lrlnb[i]*E1/cs);
				else if(j>=nsvg+nsvc+nsvb && j<2*nsvg+nsvc+nsvb)
					g[i+j*(nsv*2+3)]=-1.0*(lrkn[i][j%nsv]-lrln[i])*F2/cs;
				else if(j>=2*nsvg+nsvc+nsvb && j<2*nsvg+2*nsvc+nsvb)
					g[i+j*(nsv*2+3)]=-1.0*(lrkn[i][j%nsv]*B2/cs-lrlnc[i]*B2/cs);
				else if(j>=2*nsvg+2*nsvc+nsvb && j<2*nsvg+2*nsvc+2*nsvb)
					g[i+j*(nsv*2+3)]=-1.0*(lrkn[i][j%nsv]*E5B/cs-lrlnb[i]*E5B/cs);
				else if(j>=2*nsvg+2*nsvc+2*nsvb && j<3*nsvg+2*nsvc+2*nsvb)
					g[i+j*(nsv*2+3)]=-1.0*(lrkn[i][j%nsv]-lrln[i])*F1/cs;
				else if(j>=3*nsvg+2*nsvc+2*nsvb && j<3*nsvg+3*nsvc+2*nsvb)
					g[i+j*(nsv*2+3)]=-1.0*(lrkn[i][j%nsv]*B1/cs-lrlnc[i]*B1/cs);
				else if(j>=3*nsvg+3*nsvc+2*nsvb && j<3*nsvg+3*nsvc+3*nsvb)
					g[i+j*(nsv*2+3)]=-1.0*(lrkn[i][j%nsv]*E1/cs-lrlnb[i]*E1/cs);
				else if(j>=3*nsvg+3*nsvc+3*nsvb && j<4*nsvg+3*nsvc+3*nsvb)
					g[i+j*(nsv*2+3)]=-1.0*(lrkn[i][j%nsv]-lrln[i])*F2/cs;
				else if(j>=4*nsvg+3*nsvc+3*nsvb && j<4*nsvg+4*nsvc+3*nsvb)
					g[i+j*(nsv*2+3)]=-1.0*(lrkn[i][j%nsv]*B2/cs-lrlnc[i]*B2/cs);
				else if(j>=4*nsvg+4*nsvc+3*nsvb && j<4*nsvg+4*nsvc+4*nsvb)
					g[i+j*(nsv*2+3)]=-1.0*(lrkn[i][j%nsv]*E5B/cs-lrlnb[i]*E5B/cs);
			}
		}
		for(i=0;i<(nsv*2+3)*nsv;i++){
			if((i-3)%(nsv*2+3+1)==0 && i<(3+nsv*2)*(gpsnum-1))
				g[i] = 1.0;
			if((i-3)%(nsv*2+3+1)==0 && i>(3+nsv*2)*(gpsnum-1))
				g[i] = 1.0;
		}
		for(i=(nsv*2+3)*nsv;i<(nsv*2+3)*nsv*2;i++){
			if((i-3)%(nsv*2+3+1)==0 && i<(nsv*2+3)*nsv+(3+nsv*2)*(gpsnum-1))
				g[i] = 1.0;
			if((i-3)%(nsv*2+3+1)==0 && i>(nsv*2+3)*nsv+(3+nsv*2)*(gpsnum-1))
				g[i] = 1.0;
		}

		/////////////////

		//G�̓]�u�s��̐���
		for(i=0;i<nsv*4;i++){
			for(j=0;j<nsv*2+3;j++){
				gt[i+j*nsv*4]=g[i*(nsv*2+3)+j];
			}
		}

		for(i=0;i<PRN*PRN;i++){
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

		//G(t)*W*G�̋t�s������߂�
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
						+(SVz[1][prn]-Ref_pos[2][1])*(SVz[1][prn]-Ref_pos[2][1]))*F1/cs
						-sqrt((SVx[1][max_prn_gps]-Ref_pos[0][1])*(SVx[1][max_prn_gps]-Ref_pos[0][1])
						+(SVy[1][max_prn_gps]-Ref_pos[1][1])*(SVy[1][max_prn_gps]-Ref_pos[1][1])
						+(SVz[1][max_prn_gps]-Ref_pos[2][1])*(SVz[1][max_prn_gps]-Ref_pos[2][1]))*F1/cs;
					ref_diff2[j] = sqrt((SVx[1][prn]-Ref_pos[0][1])*(SVx[1][prn]-Ref_pos[0][1])
						+(SVy[1][prn]-Ref_pos[1][1])*(SVy[1][prn]-Ref_pos[1][1])
						+(SVz[1][prn]-Ref_pos[2][1])*(SVz[1][prn]-Ref_pos[2][1]))*F2/cs
						-sqrt((SVx[1][max_prn_gps]-Ref_pos[0][1])*(SVx[1][max_prn_gps]-Ref_pos[0][1])
						+(SVy[1][max_prn_gps]-Ref_pos[1][1])*(SVy[1][max_prn_gps]-Ref_pos[1][1])
						+(SVz[1][max_prn_gps]-Ref_pos[2][1])*(SVz[1][max_prn_gps]-Ref_pos[2][1]))*F2/cs;
					j++;
				}
			}
			for(i=0;i<SATn[rcvn];i++){//BEIDOU
				prn = SVn[rcvn][i];
				if(prn!=max_prn_c && prn<=100 && prn>=71){
					ref_diff[j] = sqrt((SVx[1][prn]-Ref_pos[0][1])*(SVx[1][prn]-Ref_pos[0][1])
						+(SVy[1][prn]-Ref_pos[1][1])*(SVy[1][prn]-Ref_pos[1][1])
						+(SVz[1][prn]-Ref_pos[2][1])*(SVz[1][prn]-Ref_pos[2][1]))*B1/cs
						-sqrt((SVx[1][max_prn_c]-Ref_pos[0][1])*(SVx[1][max_prn_c]-Ref_pos[0][1])
						+(SVy[1][max_prn_c]-Ref_pos[1][1])*(SVy[1][max_prn_c]-Ref_pos[1][1])
						+(SVz[1][max_prn_c]-Ref_pos[2][1])*(SVz[1][max_prn_c]-Ref_pos[2][1]))*B1/cs;
					ref_diff2[j] = sqrt((SVx[1][prn]-Ref_pos[0][1])*(SVx[1][prn]-Ref_pos[0][1])
						+(SVy[1][prn]-Ref_pos[1][1])*(SVy[1][prn]-Ref_pos[1][1])
						+(SVz[1][prn]-Ref_pos[2][1])*(SVz[1][prn]-Ref_pos[2][1]))*B2/cs
						-sqrt((SVx[1][max_prn_c]-Ref_pos[0][1])*(SVx[1][max_prn_c]-Ref_pos[0][1])
						+(SVy[1][max_prn_c]-Ref_pos[1][1])*(SVy[1][max_prn_c]-Ref_pos[1][1])
						+(SVz[1][max_prn_c]-Ref_pos[2][1])*(SVz[1][max_prn_c]-Ref_pos[2][1]))*B2/cs;
					j++;
				}
			}
			for(i=0;i<SATn[rcvn];i++){//Galileo
				prn = SVn[rcvn][i];
				if(prn!=max_prn_bei && prn<=70 && prn>=41){
					ref_diff[j] = sqrt((SVx[1][prn]-Ref_pos[0][1])*(SVx[1][prn]-Ref_pos[0][1])
						+(SVy[1][prn]-Ref_pos[1][1])*(SVy[1][prn]-Ref_pos[1][1])
						+(SVz[1][prn]-Ref_pos[2][1])*(SVz[1][prn]-Ref_pos[2][1]))*E1/cs
						-sqrt((SVx[1][max_prn_bei]-Ref_pos[0][1])*(SVx[1][max_prn_bei]-Ref_pos[0][1])
						+(SVy[1][max_prn_bei]-Ref_pos[1][1])*(SVy[1][max_prn_bei]-Ref_pos[1][1])
						+(SVz[1][max_prn_bei]-Ref_pos[2][1])*(SVz[1][max_prn_bei]-Ref_pos[2][1]))*E1/cs;
					ref_diff2[j] = sqrt((SVx[1][prn]-Ref_pos[0][1])*(SVx[1][prn]-Ref_pos[0][1])
						+(SVy[1][prn]-Ref_pos[1][1])*(SVy[1][prn]-Ref_pos[1][1])
						+(SVz[1][prn]-Ref_pos[2][1])*(SVz[1][prn]-Ref_pos[2][1]))*E5B/cs
						-sqrt((SVx[1][max_prn_bei]-Ref_pos[0][1])*(SVx[1][max_prn_bei]-Ref_pos[0][1])
						+(SVy[1][max_prn_bei]-Ref_pos[1][1])*(SVy[1][max_prn_bei]-Ref_pos[1][1])
						+(SVz[1][max_prn_bei]-Ref_pos[2][1])*(SVz[1][max_prn_bei]-Ref_pos[2][1]))*E5B/cs;
					j++;
				}
			}
			for(i=0;i<j;i++)
				ref_diff[i] = ref_diff[i];
			for(i=j;i<j*2;i++)
				ref_diff[i] = ref_diff2[i-j];
			for(i=j*2;i<j*3;i++)
				ref_diff[i] = ref_diff[i-j*2];
			for(i=j*3;i<j*4;i++)
				ref_diff[i] = ref_diff2[i-j*3];
		}

		j=0;
		for(i=0;i<SATn[rcvn];i++){//GPS+QZS
			prn = SVn[rcvn][i];
			if(prn!=max_prn_gps && prn<=39){
				rov_diff[j] = sqrt(lrk[0][j]*lrk[0][j]+lrk[1][j]*lrk[1][j]+lrk[2][j]*lrk[2][j])*F1/cs
					- sqrt(lrl[0]*lrl[0]+lrl[1]*lrl[1]+lrl[2]*lrl[2])*F1/cs;
				rov_diff2[j] = sqrt(lrk[0][j]*lrk[0][j]+lrk[1][j]*lrk[1][j]+lrk[2][j]*lrk[2][j])*F2/cs
					- sqrt(lrl[0]*lrl[0]+lrl[1]*lrl[1]+lrl[2]*lrl[2])*F2/cs;
				j++;
			}
		}
		for(i=0;i<SATn[rcvn];i++){//BEIDOU
			prn = SVn[rcvn][i];
			if(prn!=max_prn_c && prn<=100 && prn>=71){
				rov_diff[j] = sqrt(lrk[0][j]*lrk[0][j]+lrk[1][j]*lrk[1][j]+lrk[2][j]*lrk[2][j])*B1/cs
					- sqrt(lrlc[0]*lrlc[0]+lrlc[1]*lrlc[1]+lrlc[2]*lrlc[2])*B1/cs;
				rov_diff2[j] = sqrt(lrk[0][j]*lrk[0][j]+lrk[1][j]*lrk[1][j]+lrk[2][j]*lrk[2][j])*B2/cs
					- sqrt(lrlc[0]*lrlc[0]+lrlc[1]*lrlc[1]+lrlc[2]*lrlc[2])*B2/cs;
				j++;
			}
		}
		for(i=0;i<SATn[rcvn];i++){//Galileo
			prn = SVn[rcvn][i];
			if(prn!=max_prn_bei && prn<=70 && prn>=41){
				rov_diff[j] = sqrt(lrk[0][j]*lrk[0][j]+lrk[1][j]*lrk[1][j]+lrk[2][j]*lrk[2][j])*E1/cs
					- sqrt(lrlb[0]*lrlb[0]+lrlb[1]*lrlb[1]+lrlb[2]*lrlb[2])*E1/cs;
				rov_diff2[j] = sqrt(lrk[0][j]*lrk[0][j]+lrk[1][j]*lrk[1][j]+lrk[2][j]*lrk[2][j])*E5B/cs
					- sqrt(lrlb[0]*lrlb[0]+lrlb[1]*lrlb[1]+lrlb[2]*lrlb[2])*E5B/cs;
				j++;
			}
		}

		for(i=0;i<j;i++)
			rov_diff[i] = rov_diff[i];
		for(i=j;i<j*2;i++)
			rov_diff[i] = rov_diff2[i-j];
		for(i=j*2;i<j*3;i++)
			rov_diff[i] = rov_diff[i-j*2];
		for(i=j*3;i<j*4;i++)
			rov_diff[i] = rov_diff2[i-j*3];


		//Doppler���𗘗p����ꍇ�̏ꍇ�킯
		if(Doppler==1 && No_dp_flag==0){
			for(i=0;i<nsv*4;i++){
				dx=dx;
				dy=dy;
				dz=dz;
				if(iter==4){
					for(j=3;j<nsv*2+3;j++)
						n[j-3] = n_init[j-3];
				}
			}
		}
		else{
			for(i=0;i<nsv*4;i++){
				dx = dx+gtwg_gtw[i]*(y[i]+ref_diff[i]-rov_diff[i]);
				dy = dy+gtwg_gtw[i+nsv*4]*(y[i]+ref_diff[i]-rov_diff[i]);
				dz = dz+gtwg_gtw[i+nsv*4*2]*(y[i]+ref_diff[i]-rov_diff[i]);

				if(iter==4){
					for(j=3;j<nsv*2+3;j++)
						n[j-3] = n[j-3]+gtwg_gtw[i+nsv*4*j]*(y[i]+ref_diff[i]-rov_diff[i]);
				}
			}
		}


/////////////////////////////////////////////////////////////////////////////////////////////
		//GPS+QZS
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

		//BEIDOU
		lrlc[0] = SVx[0][max_prn_c]-(Ref_pos[0][1]+dx);
		lrlc[1] = SVy[0][max_prn_c]-(Ref_pos[1][1]+dy);
		lrlc[2] = SVz[0][max_prn_c]-(Ref_pos[2][1]+dz);
		lrlnc[0] = lrlc[0]/sqrt(lrlc[0]*lrlc[0]+lrlc[1]*lrlc[1]+lrlc[2]*lrlc[2]);
		lrlnc[1] = lrlc[1]/sqrt(lrlc[0]*lrlc[0]+lrlc[1]*lrlc[1]+lrlc[2]*lrlc[2]);
		lrlnc[2] = lrlc[2]/sqrt(lrlc[0]*lrlc[0]+lrlc[1]*lrlc[1]+lrlc[2]*lrlc[2]);
		for(i=0;i<SATn[rcvn];i++){
			prn = SVn[rcvn][i];
			if(prn!=max_prn_c && prn<=100 && prn>=71){
				lrk[0][j] = SVx[0][prn]-(Ref_pos[0][1]+dx);
				lrk[1][j] = SVy[0][prn]-(Ref_pos[1][1]+dy);
				lrk[2][j] = SVz[0][prn]-(Ref_pos[2][1]+dz);
				lrkn[0][j] = lrk[0][j]/sqrt(lrk[0][j]*lrk[0][j]+lrk[1][j]*lrk[1][j]+lrk[2][j]*lrk[2][j]);
				lrkn[1][j] = lrk[1][j]/sqrt(lrk[0][j]*lrk[0][j]+lrk[1][j]*lrk[1][j]+lrk[2][j]*lrk[2][j]);
				lrkn[2][j] = lrk[2][j]/sqrt(lrk[0][j]*lrk[0][j]+lrk[1][j]*lrk[1][j]+lrk[2][j]*lrk[2][j]);
				j++;
			}
		}

		//Galileo
		lrlb[0] = SVx[0][max_prn_bei]-(Ref_pos[0][1]+dx);
		lrlb[1] = SVy[0][max_prn_bei]-(Ref_pos[1][1]+dy);
		lrlb[2] = SVz[0][max_prn_bei]-(Ref_pos[2][1]+dz);
		lrlnb[0] = lrlb[0]/sqrt(lrlb[0]*lrlb[0]+lrlb[1]*lrlb[1]+lrlb[2]*lrlb[2]);
		lrlnb[1] = lrlb[1]/sqrt(lrlb[0]*lrlb[0]+lrlb[1]*lrlb[1]+lrlb[2]*lrlb[2]);
		lrlnb[2] = lrlb[2]/sqrt(lrlb[0]*lrlb[0]+lrlb[1]*lrlb[1]+lrlb[2]*lrlb[2]);
	//	j=0;
		for(i=0;i<SATn[rcvn];i++){
			prn = SVn[rcvn][i];
			if(prn!=max_prn_bei && prn<=70 && prn>=41){
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

	//float���̒��o
	dx_float = dx;
	dy_float = dy;
	dz_float = dz;

	//LAMBDA�@�̓K�p
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
	//�O�p�s����쐬
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
	//�Ίp�����̐�
	for(i=0;i<ndim;i++){
		det*=aa[i][i];
	}
	det = pow(det,1.0/(nsv*2*2));//ADOP


	//LAMBDA�@�̌Ăяo��
	lambda(nsv*2,2,n,q,F,s);

///////////////////////////////////////////////////////////////////////////////////////////////////////
	//LAMBDA�@�ŋ��߂��A���r�M���C�e�C�𗘗p���đ��ʌv�Z���s��
///////////////////////////////////////////////////////////////////////////////////////////////////////


	//�e�q���̋p�ɂ��d�݂̐���
	j=0;
	for(i=0;i<SATn[rcvn];i++){
		prn = SVn[rcvn][i];
		if(prn != max_prn_gps && prn<=39)
			omomi[j++] = 1.0*1.0*sqrt(90.0/Elevation[rcvn][prn]);
	}
	for(i=0;i<SATn[rcvn];i++){
		prn = SVn[rcvn][i];
		if(prn != max_prn_c && prn>=71 && prn<=100)
			omomi[j++] = 1.0*1.0*sqrt(90.0/Elevation[rcvn][prn]);
	}
	for(i=0;i<SATn[rcvn];i++){
		prn = SVn[rcvn][i];
		if(prn != max_prn_bei && prn<=70 && prn>=41)
			omomi[j++] = 1.0*1.0*sqrt(90.0/Elevation[rcvn][prn]);
	}

	dx=0;dy=0;dz=0;

	for(i=nsv*2;i<nsv*4;i++)//�c���ŏ��̌��i�A���r�M���C�e�B���̂��́j�̂ݎc��
		F[i]=0;



/////////////////////////////////////////////////////////////////////////////////////////////
	//GPS+QZS
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

	//BEIDOU
	lrlc[0] = SVx[0][max_prn_c]-(Ref_pos[0][1]+dx);
	lrlc[1] = SVy[0][max_prn_c]-(Ref_pos[1][1]+dy);
	lrlc[2] = SVz[0][max_prn_c]-(Ref_pos[2][1]+dz);
	lrlnc[0] = lrlc[0]/sqrt(lrlc[0]*lrlc[0]+lrlc[1]*lrlc[1]+lrlc[2]*lrlc[2]);
	lrlnc[1] = lrlc[1]/sqrt(lrlc[0]*lrlc[0]+lrlc[1]*lrlc[1]+lrlc[2]*lrlc[2]);
	lrlnc[2] = lrlc[2]/sqrt(lrlc[0]*lrlc[0]+lrlc[1]*lrlc[1]+lrlc[2]*lrlc[2]);
	for(i=0;i<SATn[rcvn];i++){
		prn = SVn[rcvn][i];
		if(prn!=max_prn_c && prn<=100 && prn>=71){
			lrk[0][j] = SVx[0][prn]-(Ref_pos[0][1]+dx);
			lrk[1][j] = SVy[0][prn]-(Ref_pos[1][1]+dy);
			lrk[2][j] = SVz[0][prn]-(Ref_pos[2][1]+dz);
			lrkn[0][j] = lrk[0][j]/sqrt(lrk[0][j]*lrk[0][j]+lrk[1][j]*lrk[1][j]+lrk[2][j]*lrk[2][j]);
			lrkn[1][j] = lrk[1][j]/sqrt(lrk[0][j]*lrk[0][j]+lrk[1][j]*lrk[1][j]+lrk[2][j]*lrk[2][j]);
			lrkn[2][j] = lrk[2][j]/sqrt(lrk[0][j]*lrk[0][j]+lrk[1][j]*lrk[1][j]+lrk[2][j]*lrk[2][j]);
			j++;
		}
	}

	///Galileo
	lrlb[0] = SVx[0][max_prn_bei]-(Ref_pos[0][1]+dx);
	lrlb[1] = SVy[0][max_prn_bei]-(Ref_pos[1][1]+dy);
	lrlb[2] = SVz[0][max_prn_bei]-(Ref_pos[2][1]+dz);
	lrlnb[0] = lrlb[0]/sqrt(lrlb[0]*lrlb[0]+lrlb[1]*lrlb[1]+lrlb[2]*lrlb[2]);
	lrlnb[1] = lrlb[1]/sqrt(lrlb[0]*lrlb[0]+lrlb[1]*lrlb[1]+lrlb[2]*lrlb[2]);
	lrlnb[2] = lrlb[2]/sqrt(lrlb[0]*lrlb[0]+lrlb[1]*lrlb[1]+lrlb[2]*lrlb[2]);
	for(i=0;i<SATn[rcvn];i++){
		prn = SVn[rcvn][i];
		if(prn!=max_prn_bei && prn<=70 && prn>=41){
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
		if(iter==1){
			for(i=0;i<nsvg+nsvc+nsvb;i++){
				y[i] = y[i];
			}
//			for(i=nsvg;i<nsvg+nsvb;i++){
//				y[i] = y[i];
//			}
		}

		//G�s��̐���
		for(i=0;i<=2;i++){
			for(j=0;j<nsvg;j++){
				g[i+j*3]=-1.0*(lrkn[i][j]-lrln[i])*F1/cs;
			}
			for(j=nsvg;j<nsvg+nsvc;j++){
				g[i+j*3]=-1.0*(lrkn[i][j]*B1/cs-lrlnc[i]*B1/cs);
			}
			for(j=nsvg+nsvc;j<nsvg+nsvc+nsvb;j++){
				g[i+j*3]=-1.0*(lrkn[i][j]*E1/cs-lrlnb[i]*E1/cs);
			}
		}

		//G�̓]�u�s��̐���
		for(i=0;i<nsv;i++){
			for(j=0;j<3;j++){
				gt[i+j*nsv]=g[i*3+j];
			}
		}

		//G�̓]�u�s���W�̊|���Z�@����ȍ~��gt��gtw�ƂȂ��Ă���
		for(i=0;i<nsv;i++){
			for(j=0;j<3;j++){
				gt[i+j*nsv]=gt[i+j*nsv]*omomi[i];
			}
		}

		for(i=0;i<=8;i++)
			g3[i]=0;

		//G�̓]�u�s���G�s����|����
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

		//G3�s��̋t�s������߂�
		//�t�s��v�Z//
		eps = 1.0e-25;
		minv(g3,eps,3);

		for(i=0;i<(nsv)*3;i++)
			g2[i]=0;
		//�t�s���G�̓]�u�s����|����
		for(i=0;i<nsv;i++){
			for(j=0;j<3;j++){
				g2[i]      =g2[i]+g3[j]*gt[i+nsv*j];
				g2[nsv+i]  =g2[nsv+i]+g3[j+3]*gt[i+nsv*j];
				g2[nsv*2+i]=g2[nsv*2+i]+g3[j+6]*gt[i+nsv*j];
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
					ref_diff[j] = ref_diff[j]*F1/cs;
					j++;
				}
			}

			for(i=0;i<SATn[rcvn];i++){//BEIDOU
				prn = SVn[rcvn][i];
				if(prn!=max_prn_c && prn<=100 && prn>=71){
					ref_diff[j] = sqrt((SVx[1][prn]-Ref_pos[0][1])*(SVx[1][prn]-Ref_pos[0][1])
						+(SVy[1][prn]-Ref_pos[1][1])*(SVy[1][prn]-Ref_pos[1][1])
						+(SVz[1][prn]-Ref_pos[2][1])*(SVz[1][prn]-Ref_pos[2][1]))*B1/cs
						-sqrt((SVx[1][max_prn_c]-Ref_pos[0][1])*(SVx[1][max_prn_c]-Ref_pos[0][1])
						+(SVy[1][max_prn_c]-Ref_pos[1][1])*(SVy[1][max_prn_c]-Ref_pos[1][1])
						+(SVz[1][max_prn_c]-Ref_pos[2][1])*(SVz[1][max_prn_c]-Ref_pos[2][1]))*B1/cs;
					j++;
				}
			}

			for(i=0;i<SATn[rcvn];i++){//GLONASS
				prn = SVn[rcvn][i];
				if(prn!=max_prn_bei && prn<=70 && prn>=41){
					ref_diff[j] = sqrt((SVx[1][prn]-Ref_pos[0][1])*(SVx[1][prn]-Ref_pos[0][1])
						+(SVy[1][prn]-Ref_pos[1][1])*(SVy[1][prn]-Ref_pos[1][1])
						+(SVz[1][prn]-Ref_pos[2][1])*(SVz[1][prn]-Ref_pos[2][1]))*E1/cs
						-sqrt((SVx[1][max_prn_bei]-Ref_pos[0][1])*(SVx[1][max_prn_bei]-Ref_pos[0][1])
						+(SVy[1][max_prn_bei]-Ref_pos[1][1])*(SVy[1][max_prn_bei]-Ref_pos[1][1])
						+(SVz[1][max_prn_bei]-Ref_pos[2][1])*(SVz[1][max_prn_bei]-Ref_pos[2][1]))*E1/cs;
					j++;
				}
			}
		}

		j=0;
		for(i=0;i<SATn[rcvn];i++){//GPS+QZS
			prn = SVn[rcvn][i];
			if(prn!=max_prn_gps && prn<=39){
				rov_diff[j] = sqrt(lrk[0][j]*lrk[0][j]+lrk[1][j]*lrk[1][j]+lrk[2][j]*lrk[2][j])
					- sqrt(lrl[0]*lrl[0]+lrl[1]*lrl[1]+lrl[2]*lrl[2]);
				rov_diff[j] = rov_diff[j]*F1/cs;
				j++;
			}
		}
		for(i=0;i<SATn[rcvn];i++){//BEIDOU
			prn = SVn[rcvn][i];
			if(prn!=max_prn_c && prn<=100 && prn>=71){
				rov_diff[j] = sqrt(lrk[0][j]*lrk[0][j]+lrk[1][j]*lrk[1][j]+lrk[2][j]*lrk[2][j])
					- sqrt(lrlc[0]*lrlc[0]+lrlc[1]*lrlc[1]+lrlc[2]*lrlc[2]);
				rov_diff[j] = rov_diff[j]*B1/cs;
				j++;
			}
		}
		for(i=0;i<SATn[rcvn];i++){//Galileo
			prn = SVn[rcvn][i];
			if(prn!=max_prn_bei && prn<=70 && prn>=41){
				rov_diff[j] = sqrt(lrk[0][j]*lrk[0][j]+lrk[1][j]*lrk[1][j]+lrk[2][j]*lrk[2][j])*E1/cs
					- sqrt(lrlb[0]*lrlb[0]+lrlb[1]*lrlb[1]+lrlb[2]*lrlb[2])*E1/cs;
				j++;
			}
		}

		//dx,dy,dz
		for(i=0;i<nsv;i++){
			dx = dx+g2[i]*(y[i]+ref_diff[i]-rov_diff[i]-F[i]);
			dy = dy+g2[i+nsv]*(y[i]+ref_diff[i]-rov_diff[i]-F[i]);
			dz = dz+g2[i+nsv*2]*(y[i]+ref_diff[i]-rov_diff[i]-F[i]);
		}

/////////////////////////////////////////////////////////////////////////////////////////////
		//GPS+QZS
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

		//BEIDOU
		lrlc[0] = SVx[0][max_prn_c]-(Ref_pos[0][1]+dx);
		lrlc[1] = SVy[0][max_prn_c]-(Ref_pos[1][1]+dy);
		lrlc[2] = SVz[0][max_prn_c]-(Ref_pos[2][1]+dz);
		lrlnc[0] = lrlc[0]/sqrt(lrlc[0]*lrlc[0]+lrlc[1]*lrlc[1]+lrlc[2]*lrlc[2]);
		lrlnc[1] = lrlc[1]/sqrt(lrlc[0]*lrlc[0]+lrlc[1]*lrlc[1]+lrlc[2]*lrlc[2]);
		lrlnc[2] = lrlc[2]/sqrt(lrlc[0]*lrlc[0]+lrlc[1]*lrlc[1]+lrlc[2]*lrlc[2]);
		for(i=0;i<SATn[rcvn];i++){
			prn = SVn[rcvn][i];
			if(prn!=max_prn_c && prn<=100 && prn>=71){
				lrk[0][j] = SVx[0][prn]-(Ref_pos[0][1]+dx);
				lrk[1][j] = SVy[0][prn]-(Ref_pos[1][1]+dy);
				lrk[2][j] = SVz[0][prn]-(Ref_pos[2][1]+dz);
				lrkn[0][j] = lrk[0][j]/sqrt(lrk[0][j]*lrk[0][j]+lrk[1][j]*lrk[1][j]+lrk[2][j]*lrk[2][j]);
				lrkn[1][j] = lrk[1][j]/sqrt(lrk[0][j]*lrk[0][j]+lrk[1][j]*lrk[1][j]+lrk[2][j]*lrk[2][j]);
				lrkn[2][j] = lrk[2][j]/sqrt(lrk[0][j]*lrk[0][j]+lrk[1][j]*lrk[1][j]+lrk[2][j]*lrk[2][j]);
				j++;
			}
		}

		//GLONASS
		lrlb[0] = SVx[0][max_prn_bei]-(Ref_pos[0][1]+dx);
		lrlb[1] = SVy[0][max_prn_bei]-(Ref_pos[1][1]+dy);
		lrlb[2] = SVz[0][max_prn_bei]-(Ref_pos[2][1]+dz);
		lrlnb[0] = lrlb[0]/sqrt(lrlb[0]*lrlb[0]+lrlb[1]*lrlb[1]+lrlb[2]*lrlb[2]);
		lrlnb[1] = lrlb[1]/sqrt(lrlb[0]*lrlb[0]+lrlb[1]*lrlb[1]+lrlb[2]*lrlb[2]);
		lrlnb[2] = lrlb[2]/sqrt(lrlb[0]*lrlb[0]+lrlb[1]*lrlb[1]+lrlb[2]*lrlb[2]);
		for(i=0;i<SATn[rcvn];i++){
			prn = SVn[rcvn][i];
			if(prn!=max_prn_bei && prn<=70 && prn>=41){
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

	lat_diff = (t[0]-POSrcvlat[0])*111319.49;
	lon_diff = (t[1]-POSrcvlon[0])*cos(t[0]*pi/180.0)*111319.49;
	hgt_diff = t[2] - POSrcvhgt[0];

//	lat_diff = (t[0]-Ref_lat)*111319.49;
//	lon_diff = (t[1]-Ref_lon)*cos(t[0]*pi/180.0)*111319.49;
//	hgt_diff = t[2] - Ref_hgt;

	s_factor = s[1]/s[0];

	//ADOP�ł̃`�F�b�N
	if(det>=1.00)
		s_factor = 1.0;

	double aaa,bbb,xxx,yyy;
	aaa=Ref_lat;
	bbb=Ref_lon;
	xxx = (t[0]-POSrcvlat[0])*111319.49;
	yyy = (t[1]-POSrcvlon[0])*cos(t[0]*pi/180.0)*111319.49;

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