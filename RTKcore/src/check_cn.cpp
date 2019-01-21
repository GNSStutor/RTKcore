///////////////////////////////////////////////////////////////////
//
//  �M�����x�ɂ��q���̎g�p���f���s��
//  check��L1�т̂݁i�}���`GNSS�j
//  ���{�ł̐ݒ�ɂȂ��Ă��邱�Ƃɒ���
//
///////////////////////////////////////////////////////////////////

#include <iostream>

#include "global_extern.h"
#include <math.h>

using namespace std;
const double pi     = 3.1415926535898;//�~����
const double F1     = 1575420000.0;
const double F2     = 1227600000.0;
const double cs     = 299792458.0;//����

void check_cn(int rcvn)
{
	int i;
	int prn;

	double set1=8.0;//�p�ɑ΂���M�����x�̕��ϒl��肱�̐ݒ�l�����Ⴂ�q���͔r��
	double threshold_l1=0;

	//������
	for(i=1;i<=PRN-1;i++){
		Reject_cn[rcvn][i]=0;
	}

	for(i=0;i<SATn[rcvn];i++){
		prn = SVn[rcvn][i];

		if(Elevation[rcvn][prn]<=45 && (prn<=39||(prn>=41&&prn<=70))){//GPS,QZS,GALILEO
			threshold_l1 = 0.2286*Elevation[rcvn][prn]+39.7;
		}
		if(Elevation[rcvn][prn]<=45 &&( prn>=71&&prn<100)){//BEIDOU
			threshold_l1 = 0.2286*Elevation[rcvn][prn]+36.2;
		}
		if(Elevation[rcvn][prn]<=45 && (prn>=101)){//GLONASS
			threshold_l1 = 0.2286*Elevation[rcvn][prn]+39.7-5.0;
		}
		if(Elevation[rcvn][prn]>45 && (prn<=39||(prn>=41&&prn<=70))){
			threshold_l1 = 50.0;
		}
		if(Elevation[rcvn][prn]>45 &&( prn>=71&&prn<100)){
			threshold_l1 = 46.5;
		}
		if(Elevation[rcvn][prn]>45 && (prn>=101)){//GLONASS
			threshold_l1 = 50.0-5;
		}

		if(prn==71){//BEIDOU GEO1
			threshold_l1 = 40.6;
		}
		if(prn==72){//BEIDOU GEO2
				threshold_l1 = 34.0;
		}
		if(prn==73){//BEIDOU GEO3
				threshold_l1 = 40.2;
		}
		if(prn==74){//BEIDOU GEO4
			threshold_l1 = 40;	
		}

		if(Cn1[rcvn][prn]<=threshold_l1-set1){
				Reject_cn[rcvn][prn]=1;
		}
	}
}
