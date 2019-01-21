/////////////////////////////////////////////////////////////////////////////////
//
// UBLOX(M8P,M8T)�̎�M�@�ϑ��f�[�^�𗘗p����RTK�̉��Z
// GPS/QZSS L1, GALILEO E1, BeiDou B1, GLONASS G1
// 2018/11/30 F9P�p�ɏ������C�iGPS/QZSS L1/L2, GALILEO E1/E5B, BeiDou B1/B2, GONASS G1�j
//
/////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "global.h"

using namespace std;

int main(){
	int i,iter;//�ǂݍ��݉�
	int rcvn;//rcvn=1�i����j�@rcvn=0�i�ړ����j
	int pos=1;//�H
	int rover_kaisu=0;
	static int SATn_moto=0,SVn_moto[PRN]={0};
	double trange[PRN]={0},stock_cn1[RCVN][PRN]={0};
	End_flag = 0;

	cout.precision(7);
	Sol_flag[0]=0;Sol_flag[1]=0;Sol_flag[2]=0;//�ړ����̓ǂݍ��݉񐔁A�P�Ƒ��ʉ񐔁ARTK��FIX��

	set_initial_value();//�����ݒ�{�t�@�C���ݒ�
	cout << "GPSTIME" << " " << "Ref" << " " << "Rov" << " " << "�S��" << " " << "���ʉ�" << " " << "FIX��" << endl;

	for(iter=1;iter<=Iteration;iter++){//�ݒ�񐔕��ǂݍ���

		rcvn = 1;//1:��� 0:�ړ���
		read_data(rcvn);//�ϑ��f�[�^�A�q�@���b�Z�[�W�̓ǂݍ���

		if(End_flag==1){//�t�@�C���ǂݍ��ݍŌ�ŋ����I��
			iter = Iteration;continue;
		}

		POS=0;
		calc_satpos(rcvn);//�G�t�F�����X����̉q���ʒu�v�Z
		calc_direction(rcvn,iter);//�p�E���ʊp�v�Z�i����Ȃ̂Ń��[�U�ʒu�͌Œ�̊�ʒu�Ōv�Z�j
		calc_iono_model(rcvn);//�d���w���f���ł̓d���w�x���ʌv�Z
		calc_tropo(rcvn);//�Η����x���ʌv�Z

		Glo_Num=0;Gps_Num=0;Gal_Num=0;Bei_Num=0;//�e�q���V�X�e���̐���������
		choose_sat(rcvn,iter);//�q���I��
		
		//�P�Ƒ��ʂ��s���iGPS�q����4�@�ȏ゠�邱�Ƃ��O��j
		POS=1;//�P�Ƒ���
		if(SATn[rcvn]>=MinSvNum){
			calc_pos(rcvn,iter,pos);//�ŏ����@�ł̒P�Ƒ���
			calc_satpos(rcvn);//��M�@�N���b�N�덷�ɂ��[���������ς�邽�߁A�ēx���ˎ������Čv�Z���A�q���ʒu���Čv�Z����
			calc_pos2(rcvn,iter,pos);//�����C�������q���ʒu�ł̒P�Ƒ��ʍČv�Z
		}

///*
	//�����Ȋ�_�ʒu�𗘗p���Đ^�̋����ƕ␳�l�iDGNSS�p�j���Z�o
	for(i=0;i<SATn[rcvn];i++){
		prn = SVn[rcvn][i];

		//�^�̋����Z�o
		trange[prn] = sqrt((SVx[rcvn][prn]-Ref_pos[0][1])*(SVx[rcvn][prn]-Ref_pos[0][1])
			+ (SVy[rcvn][prn]-Ref_pos[1][1])*(SVy[rcvn][prn]-Ref_pos[1][1])
			+ (SVz[rcvn][prn]-Ref_pos[2][1])*(SVz[rcvn][prn]-Ref_pos[2][1]));

		Correction[prn] = trange[prn] - (Pr1[rcvn][prn] + SV_corrtime[rcvn][prn]);
	}
//*/

		//RTK�ɓ���ꍇ�ړ����̊ϑ��f�[�^��ǂݍ��݁A�������r�i����ƈړ����j����RTK�̉��Z�֐i��
		if(RTK == 1){
			rcvn=0;//1:��� 0:�ړ���
			if(GPSTIME>=604800)
				GPSTIME=GPSTIME-604800;
			
			//��Ǌϑ��f�[�^�̎����ƈړ��Ǌϑ��f�[�^�̎������r���ǂݍ���
			if(GPSTIME<610000){
				while((GPSTIME-DGPSTIME)>0.05 || rover_kaisu==0){
					read_data(0);//�ړ����ϑ��f�[�^�̓ǂݍ���
					rover_kaisu++;
				}
			}

			//�������f��������
//			if(fabs(GPSTIME - DGPSTIME) < 0.1 && SATn[1]>=0){//����GPS�����̏ꍇ�i1Hz���m���j
			while(DGPSTIME-GPSTIME >= -0.01 && DGPSTIME-GPSTIME <= 0.995 && DGPSTIME>=0 && DGPSTIME<=610000){//�����1Hz�A�ړ�����5Hz���̏ꍇ

				read_data(0);//�オwhile���̂Ƃ�����
				rover_kaisu++;
				Sol_flag[0]++;//�ړ����̑S�񐔃J�E���g

				calc_satpos(0);///�G�t�F�����X����̉q���ʒu�v�Z
				calc_direction(0,iter);//�p�E���ʊp�v�Z�i�ړ����Ȃ̂Ń��[�U�ʒu����j
				calc_iono_model(rcvn);//�d���w���f���ł̓d���w�x���ʌv�Z
				calc_tropo(0);//�Η����x���ʂ̌v�Z

				check_cn(rcvn);

				Glo_Num=0;Gps_Num=0;Gal_Num=0;Bei_Num=0;
				choose_sat(0,iter);//�q���I��

				fprintf(fp[9],"%f,%d,%d,%d,%d,%d\n",DGPSTIME,SATn[rcvn],Gps_Num,Gal_Num,Bei_Num,Glo_Num);

				//�ړ����̑��ʂ��s��
				POS=1;//�P�Ƒ���
				POS=2;//DGNSS
				MinSvNum=4;
				if(SATn[rcvn]>=5 && Gps_Num>=1){//DGNSS�̍Œ�q������O�̂���5�@�i�}���`GNSS�Łj�Ƃ��Ă��܂�
					calc_pos(rcvn,iter,pos);//�ŏ����@�ł̒P�Ƒ���
					calc_satpos(rcvn);//��M�@�N���b�N�덷�ɂ��[���������ς�邽�߁A�ēx���ˎ������Čv�Z���A�q���ʒu���Čv�Z����
					calc_pos2(rcvn,iter,pos);//�����C�������q���ʒu�ł̒P�Ƒ��ʍČv�Z
					Sol_flag[1]++;//���ʂ̉񐔃J�E���g
				}

				Ratio=0;

///////////////////////////////////////////////////////////////////////////////////////////////
				//�i�֐��ɓ���Ȃ��悤100�@�ȏ�Ƃ��Ă��܂��j
				//�O��RTK�R�A��RTK����1���g�̉ӏ�
				//RTK���s������(GPS/QZS/GALILEO�̏ꍇ)�@���킹��DGNSS���s��
				if(SATn[rcvn]>=100 && Bei_Num==0 && Glo_Num==0)
					calc_rtk_GQE(rcvn);

				//RTK���s������(GPS/QZS/GALILEO+BeiDou�̏ꍇ)�@���킹��DGNSS���s��
				if(SATn[rcvn]>=100 && Gps_Num+Gal_Num>=2 && Bei_Num>=2 && Gps_Num>=2)
					calc_rtk_GQEB(rcvn);

				//RTK���s������(GPS/QZS/GALILEO+GLONASS�̏ꍇ)�@���킹��DGNSS���s��
				if(SATn[rcvn]>=100 && Gps_Num+Gal_Num>=2 && Glo_Num>=2 && Gps_Num>=2)
					calc_rtk_GQER(rcvn);
///////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////////////
				//�V���ɕt������RTK����2���g�̉ӏ��ł��iGPS+QZS+GALILEO+BeiDou�܂őΉ����Ă��܂��j

				//RTK���s������(GPS/QZS+Galileo+BeiDou�̏ꍇ 2���g)
				if(SATn[rcvn]>=7 && Gps_Num>=2 && Bei_Num>=2 && Gal_Num>=2 && Glo_Num<0.1)
					calc_rtk_GQEB_F9P(rcvn);
///*
				//RTK���s������(GPS/QZS+BeiDou�̏ꍇ 2���g)
				if(Ratio<Ratio_limit){
					for(i=41;i<=70;i++){
						stock_cn1[rcvn][i] = Cn1[rcvn][i];
						Cn1[rcvn][i]=0;
					}

					for(i=0;i<SATn[rcvn];i++)
						SVn_moto[i] = SVn[rcvn][i];
					SATn_moto = SATn[rcvn];

					Glo_Num=0;Gps_Num=0;Gal_Num=0;Bei_Num=0;
					choose_sat(0,iter);//GALILEO�������ĉq���I��

					if(SATn[rcvn]>=6 && Gps_Num>=2 && Bei_Num>=2 && Gal_Num<0.1 && Glo_Num<0.1)
						calc_rtk_GQB_F9P(rcvn);
				}

				//RTK���s������(GPS/QZS+Galileo�̏ꍇ 2���g)
				if(Ratio<Ratio_limit){
					for(i=41;i<=70;i++)
						Cn1[rcvn][i] = stock_cn1[rcvn][i];

					for(i=0;i<SATn[rcvn];i++)
						SVn[rcvn][i] = SVn_moto[i];
					SATn[rcvn] = SATn_moto;

					for(i=71;i<=100;i++)
						Cn1[rcvn][i]=0;

					Glo_Num=0;Gps_Num=0;Gal_Num=0;Bei_Num=0;
					choose_sat(0,iter);//BeiDou�������ĉq���I��
			
					if(SATn[rcvn]>=6 && Gps_Num>=2 && Gal_Num>=2 && Bei_Num<0.1 && Glo_Num<0.1)
						calc_rtk_GQE_F9P(rcvn);
				}

				//RTK���s������(GPS/QZS�̏ꍇ 2���g)
				if(Ratio<Ratio_limit){
					for(i=41;i<=100;i++)
						Cn1[rcvn][i]=0;

					for(i=0;i<SATn[rcvn];i++)
						SVn[rcvn][i] = SVn_moto[i];
					SATn[rcvn] = SATn_moto;

					Glo_Num=0;Gps_Num=0;Gal_Num=0;Bei_Num=0;
					choose_sat(0,iter);//Galileo,BeiDou�������ĉq���I��
			
					if(SATn[rcvn]>=5 && Gps_Num>=2 && Gal_Num<0.1 && Bei_Num<0.1 && Glo_Num<0.1)
						calc_rtk_GQ_F9P(rcvn);
				}
//				*/
///////////////////////////////////////////////////////////////////////////////////////////////

			}//��ǎ����ɂ��킹���ړ����̃^�C�~���O�ł�if�܂���while���̏I���

		}//DGNSS�{RTK

		if(((int)iter%100)==0 && GPSTIME>=-0.1 && DGPSTIME>=-0.1)//�r���o�߂̏����o��
			cout << GPSTIME << " " << SATn[1] << " " << SATn[0] << " " << Sol_flag[0] << " " << Sol_flag[1] << " " << Sol_flag[2] << endl;
		
	}//�����܂ł��J��Ԃ��v�Z��

	file_close(RTK);//�t�@�C�������
	return(0);
}