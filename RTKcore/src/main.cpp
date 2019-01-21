/////////////////////////////////////////////////////////////////////////////////
//
// UBLOX(M8P,M8T)の受信機観測データを利用したRTKの演算
// GPS/QZSS L1, GALILEO E1, BeiDou B1, GLONASS G1
// 2018/11/30 F9P用に少し改修（GPS/QZSS L1/L2, GALILEO E1/E5B, BeiDou B1/B2, GONASS G1）
//
/////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "global.h"

using namespace std;

int main(){
	int i,iter;//読み込み回数
	int rcvn;//rcvn=1（基準側）　rcvn=0（移動側）
	int pos=1;//？
	int rover_kaisu=0;
	static int SATn_moto=0,SVn_moto[PRN]={0};
	double trange[PRN]={0},stock_cn1[RCVN][PRN]={0};
	End_flag = 0;

	cout.precision(7);
	Sol_flag[0]=0;Sol_flag[1]=0;Sol_flag[2]=0;//移動側の読み込み回数、単独測位回数、RTKのFIX回数

	set_initial_value();//初期設定＋ファイル設定
	cout << "GPSTIME" << " " << "Ref" << " " << "Rov" << " " << "全回数" << " " << "測位回数" << " " << "FIX回数" << endl;

	for(iter=1;iter<=Iteration;iter++){//設定回数分読み込む

		rcvn = 1;//1:基準側 0:移動側
		read_data(rcvn);//観測データ、航法メッセージの読み込み

		if(End_flag==1){//ファイル読み込み最後で強制終了
			iter = Iteration;continue;
		}

		POS=0;
		calc_satpos(rcvn);//エフェメリスからの衛星位置計算
		calc_direction(rcvn,iter);//仰角・方位角計算（基準側なのでユーザ位置は固定の基準位置で計算）
		calc_iono_model(rcvn);//電離層モデルでの電離層遅延量計算
		calc_tropo(rcvn);//対流圏遅延量計算

		Glo_Num=0;Gps_Num=0;Gal_Num=0;Bei_Num=0;//各衛星システムの数を初期化
		choose_sat(rcvn,iter);//衛星選択
		
		//単独測位を行う（GPS衛星が4機以上あることが前提）
		POS=1;//単独測位
		if(SATn[rcvn]>=MinSvNum){
			calc_pos(rcvn,iter,pos);//最小二乗法での単独測位
			calc_satpos(rcvn);//受信機クロック誤差により擬似距離が変わるため、再度発射時刻を再計算し、衛星位置も再計算する
			calc_pos2(rcvn,iter,pos);//少し修正した衛星位置での単独測位再計算
		}

///*
	//精密な基準点位置を利用して真の距離と補正値（DGNSS用）を算出
	for(i=0;i<SATn[rcvn];i++){
		prn = SVn[rcvn][i];

		//真の距離算出
		trange[prn] = sqrt((SVx[rcvn][prn]-Ref_pos[0][1])*(SVx[rcvn][prn]-Ref_pos[0][1])
			+ (SVy[rcvn][prn]-Ref_pos[1][1])*(SVy[rcvn][prn]-Ref_pos[1][1])
			+ (SVz[rcvn][prn]-Ref_pos[2][1])*(SVz[rcvn][prn]-Ref_pos[2][1]));

		Correction[prn] = trange[prn] - (Pr1[rcvn][prn] + SV_corrtime[rcvn][prn]);
	}
//*/

		//RTKに入る場合移動側の観測データを読み込み、時刻を比較（基準側と移動側）してRTKの演算へ進む
		if(RTK == 1){
			rcvn=0;//1:基準側 0:移動側
			if(GPSTIME>=604800)
				GPSTIME=GPSTIME-604800;
			
			//基準局観測データの時刻と移動局観測データの時刻を比較し読み込む
			if(GPSTIME<610000){
				while((GPSTIME-DGPSTIME)>0.05 || rover_kaisu==0){
					read_data(0);//移動側観測データの読み込み
					rover_kaisu++;
				}
			}

			//時刻判断処理部分
//			if(fabs(GPSTIME - DGPSTIME) < 0.1 && SATn[1]>=0){//同じGPS時刻の場合（1Hz同士等）
			while(DGPSTIME-GPSTIME >= -0.01 && DGPSTIME-GPSTIME <= 0.995 && DGPSTIME>=0 && DGPSTIME<=610000){//基準側が1Hz、移動側が5Hz等の場合

				read_data(0);//上がwhile文のときだけ
				rover_kaisu++;
				Sol_flag[0]++;//移動側の全回数カウント

				calc_satpos(0);///エフェメリスからの衛星位置計算
				calc_direction(0,iter);//仰角・方位角計算（移動側なのでユーザ位置を基準）
				calc_iono_model(rcvn);//電離層モデルでの電離層遅延量計算
				calc_tropo(0);//対流圏遅延量の計算

				check_cn(rcvn);

				Glo_Num=0;Gps_Num=0;Gal_Num=0;Bei_Num=0;
				choose_sat(0,iter);//衛星選択

				fprintf(fp[9],"%f,%d,%d,%d,%d,%d\n",DGPSTIME,SATn[rcvn],Gps_Num,Gal_Num,Bei_Num,Glo_Num);

				//移動側の測位を行う
				POS=1;//単独測位
				POS=2;//DGNSS
				MinSvNum=4;
				if(SATn[rcvn]>=5 && Gps_Num>=1){//DGNSSの最低衛星数を念のため5機（マルチGNSSで）としています
					calc_pos(rcvn,iter,pos);//最小二乗法での単独測位
					calc_satpos(rcvn);//受信機クロック誤差により擬似距離が変わるため、再度発射時刻を再計算し、衛星位置も再計算する
					calc_pos2(rcvn,iter,pos);//少し修正した衛星位置での単独測位再計算
					Sol_flag[1]++;//測位の回数カウント
				}

				Ratio=0;

///////////////////////////////////////////////////////////////////////////////////////////////
				//（関数に入らないよう100機以上としています）
				//前回RTKコアのRTK測位1周波の箇所
				//RTKを行う部分(GPS/QZS/GALILEOの場合)　あわせてDGNSSも行う
				if(SATn[rcvn]>=100 && Bei_Num==0 && Glo_Num==0)
					calc_rtk_GQE(rcvn);

				//RTKを行う部分(GPS/QZS/GALILEO+BeiDouの場合)　あわせてDGNSSも行う
				if(SATn[rcvn]>=100 && Gps_Num+Gal_Num>=2 && Bei_Num>=2 && Gps_Num>=2)
					calc_rtk_GQEB(rcvn);

				//RTKを行う部分(GPS/QZS/GALILEO+GLONASSの場合)　あわせてDGNSSも行う
				if(SATn[rcvn]>=100 && Gps_Num+Gal_Num>=2 && Glo_Num>=2 && Gps_Num>=2)
					calc_rtk_GQER(rcvn);
///////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////////////
				//新たに付加したRTK測位2周波の箇所です（GPS+QZS+GALILEO+BeiDouまで対応しています）

				//RTKを行う部分(GPS/QZS+Galileo+BeiDouの場合 2周波)
				if(SATn[rcvn]>=7 && Gps_Num>=2 && Bei_Num>=2 && Gal_Num>=2 && Glo_Num<0.1)
					calc_rtk_GQEB_F9P(rcvn);
///*
				//RTKを行う部分(GPS/QZS+BeiDouの場合 2周波)
				if(Ratio<Ratio_limit){
					for(i=41;i<=70;i++){
						stock_cn1[rcvn][i] = Cn1[rcvn][i];
						Cn1[rcvn][i]=0;
					}

					for(i=0;i<SATn[rcvn];i++)
						SVn_moto[i] = SVn[rcvn][i];
					SATn_moto = SATn[rcvn];

					Glo_Num=0;Gps_Num=0;Gal_Num=0;Bei_Num=0;
					choose_sat(0,iter);//GALILEOを除いて衛星選択

					if(SATn[rcvn]>=6 && Gps_Num>=2 && Bei_Num>=2 && Gal_Num<0.1 && Glo_Num<0.1)
						calc_rtk_GQB_F9P(rcvn);
				}

				//RTKを行う部分(GPS/QZS+Galileoの場合 2周波)
				if(Ratio<Ratio_limit){
					for(i=41;i<=70;i++)
						Cn1[rcvn][i] = stock_cn1[rcvn][i];

					for(i=0;i<SATn[rcvn];i++)
						SVn[rcvn][i] = SVn_moto[i];
					SATn[rcvn] = SATn_moto;

					for(i=71;i<=100;i++)
						Cn1[rcvn][i]=0;

					Glo_Num=0;Gps_Num=0;Gal_Num=0;Bei_Num=0;
					choose_sat(0,iter);//BeiDouを除いて衛星選択
			
					if(SATn[rcvn]>=6 && Gps_Num>=2 && Gal_Num>=2 && Bei_Num<0.1 && Glo_Num<0.1)
						calc_rtk_GQE_F9P(rcvn);
				}

				//RTKを行う部分(GPS/QZSの場合 2周波)
				if(Ratio<Ratio_limit){
					for(i=41;i<=100;i++)
						Cn1[rcvn][i]=0;

					for(i=0;i<SATn[rcvn];i++)
						SVn[rcvn][i] = SVn_moto[i];
					SATn[rcvn] = SATn_moto;

					Glo_Num=0;Gps_Num=0;Gal_Num=0;Bei_Num=0;
					choose_sat(0,iter);//Galileo,BeiDouを除いて衛星選択
			
					if(SATn[rcvn]>=5 && Gps_Num>=2 && Gal_Num<0.1 && Bei_Num<0.1 && Glo_Num<0.1)
						calc_rtk_GQ_F9P(rcvn);
				}
//				*/
///////////////////////////////////////////////////////////////////////////////////////////////

			}//基準局時刻にあわせた移動側のタイミングでのifまたはwhile文の終わり

		}//DGNSS＋RTK

		if(((int)iter%100)==0 && GPSTIME>=-0.1 && DGPSTIME>=-0.1)//途中経過の書き出し
			cout << GPSTIME << " " << SATn[1] << " " << SATn[0] << " " << Sol_flag[0] << " " << Sol_flag[1] << " " << Sol_flag[2] << endl;
		
	}//ここまでが繰り返し計算部

	file_close(RTK);//ファイルを閉じる
	return(0);
}