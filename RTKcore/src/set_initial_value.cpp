////////////////////////////////////////////////////////////////////////////////
//
// 初期設定
//
////////////////////////////////////////////////////////////////////////////////
#include <math.h>
#include <iostream>

#include "global_extern.h"

using namespace std;

void set_initial_value()
{
	int i,j,k;
	char c,buff[128];
	char file1[128],file2[128],file3[128];
	errno_t error;

	if((fopen_s(&fp[8],"initial_setting.txt","r")) != 0){
		cout << "ファイルを開けない initial_setting.txt" << endl;exit(1);
	}

	for(k=1; k<=18; k++){
		for(j=0; j<100; j++){
			i = 0;
			while((c = fgetc(fp[8])) != ',' && c != '\n' && c != EOF && c != ' ')
			{
				buff[i++] = c;
			}
			buff[i] = '\0';

			if(c=='\n'){
				switch(k){
					case 1: Elevation_mask1=atof(buff); break;
					case 2: strcpy_s(file1,128,buff); break;
					case 3: strcpy_s(file2,128,buff); break;
					case 4: strcpy_s(file3,128,buff); break;
					case 5: POSrcvlat[1]=atof(buff); break;
					case 6: POSrcvlon[1]=atof(buff); break;
					case 7: POSrcvhgt[1]=atof(buff); break;
					case 8: Code_noise=atof(buff); break;
					case 9: Carrier_noise=atof(buff); break;
					case 10: Iteration=atoi(buff); break;
					case 11: Threshold_cn=atof(buff); break;
					case 12: RTK=atoi(buff); break;
					case 13: Ratio_limit=atof(buff); break;
					case 14: Gflag = atoi(buff);  break;
					case 15: Jflag = atoi(buff);  break;
					case 16: Eflag = atoi(buff);  break;
					case 17: Cflag = atoi(buff);  break;
					case 18: Rflag = atoi(buff);  break;
					default: cout << "error in initial setting" << endl;break;
				}
				j=100;
			}
		}
	}

/////////////////////////////////////////////////////////////////////////////////
	//////////////// 緯度、経度、高度を地心座標系に変換しておく /////////////////
	int rcvn=1;//reference
    double a,b,f,n,pi;
	pi = 3.14159265359;	
	a = 6378137.0;
	f = 1/298.257223563;
	b = a*(1-f);
	n = a*a/sqrt(a*a*cos(POSrcvlat[rcvn]*pi/180)*
		cos(POSrcvlat[rcvn]*pi/180)+b*b*sin(POSrcvlat[rcvn]*pi/180)
		*sin(POSrcvlat[rcvn]*pi/180));

	Ref_pos[0][rcvn] = (n+POSrcvhgt[rcvn]+POSrcvund[rcvn])*
		cos(POSrcvlat[rcvn]*pi/180)*cos(POSrcvlon[rcvn]*pi/180);
	Ref_pos[1][rcvn] = (n+POSrcvhgt[rcvn]+POSrcvund[rcvn])*
		cos(POSrcvlat[rcvn]*pi/180)*sin(POSrcvlon[rcvn]*pi/180);
	Ref_pos[2][rcvn] = (n*b*b/a/a+POSrcvhgt[rcvn]+POSrcvund[rcvn])*
		sin(POSrcvlat[rcvn]*pi/180);

	//6/23 rover
	POSrcvlat[0]=POSrcvlat[1];
	POSrcvlon[0]=POSrcvlon[1];
	POSrcvhgt[0]=POSrcvhgt[1];

	rcvn=0;//rover
	pi = 3.14159265359;	
	a = 6378137.0;
	f = 1/298.257223563;
	b = a*(1-f);
	n = a*a/sqrt(a*a*cos(POSrcvlat[rcvn]*pi/180)*
		cos(POSrcvlat[rcvn]*pi/180)+b*b*sin(POSrcvlat[rcvn]*pi/180)
		*sin(POSrcvlat[rcvn]*pi/180));

	Ref_pos[0][rcvn] = (n+POSrcvhgt[rcvn]+POSrcvund[rcvn])*
		cos(POSrcvlat[rcvn]*pi/180)*cos(POSrcvlon[rcvn]*pi/180);
	Ref_pos[1][rcvn] = (n+POSrcvhgt[rcvn]+POSrcvund[rcvn])*
		cos(POSrcvlat[rcvn]*pi/180)*sin(POSrcvlon[rcvn]*pi/180);
	Ref_pos[2][rcvn] = (n*b*b/a/a+POSrcvhgt[rcvn]+POSrcvund[rcvn])*
		sin(POSrcvlat[rcvn]*pi/180);
/////////////////////////////////////////////////////////////////////////////////

	//電離層モデルの係数を設定する（手入力、全ての衛星はこれで計算）
	Arp[0]=1.3411E-08;
	Arp[1]=-4.8993E-09;
	Arp[2]=-3.2170E-07;
	Arp[3]=-5.7500E-07 ;
	Bet[0]=1.1419E+05;
	Bet[1]=-1.4480E+05;
	Bet[2]=-6.0313E+05;
	Bet[3]=-2.0756E+05;

	//2017/11/24
   //1.3411D-08 -4.8993D-09 -3.2170D-07 -5.7500D-07          ION ALPHA           
   //1.1419D+05 -1.4480D+05 -6.0313D+05 -2.0756D+05          ION BETA            

	User_lat = POSrcvlat[1];//移動側の初期緯度
	User_lon = POSrcvlon[1];//移動側の初期経度
	User_pos[0] = Ref_pos[0][1];//移動側の初期X
	User_pos[1] = Ref_pos[1][1];//移動側の初期Y
	User_pos[2] = Ref_pos[2][1];//移動側の初期Z

	LeapSecond = 18.0;//うるう秒の設定
	Interval = 0.2;//観測データの間隔
//	Doppler = 0;

////////////////////// 観測データの読み込み　ファイル指定 /////////////////////////
///*
	//rover
	if((fopen_s(&fp[0],file1,"r")) != 0){
		cout << "ファイルを開けない0" << endl;exit(1);
	}
	//base
	if((fopen_s(&fp[1],file2,"r")) != 0){
		cout << "ファイルを開けない1" << endl;exit(1);
	}
	//nav
	if((fopen_s(&fp[2],file3,"r")) != 0){
		cout << "ファイルを開けない2" << endl;exit(1);
	}

//*/

////////////////////// 観測データの読み込み　ファイル指定 /////////////////////////

////////////////////// 出力する結果のためのファイル設定   /////////////////////////

	//単独測位結果
	if((fopen_s(&fp[3],"pos.csv","w")) != 0){
		cout << "出力ファイルを開けない3" << endl;exit(1);
	}
	//DGNSS結果
	if((fopen_s(&fp[4],"dgnss.csv","w")) != 0){
		cout << "出力ファイルを開けない4" << endl;exit(1);
	}
	//RTK結果
	if((fopen_s(&fp[5],"rtk.csv","w")) != 0){
		cout << "出力ファイルを開けない5"<< endl;exit(1);
	}
	//テスト用
	if((fopen_s(&fp[7],"rtkplot.pos","w")) != 0){
		cout << "出力ファイルを開けない7" << endl;exit(1);
	}
	//テスト用
	if((fopen_s(&fp[9],"test.csv","w")) != 0){
		cout << "出力ファイルを開けない9" << endl;exit(1);
	}


////////////////////// 出力する結果のためのファイル設定   /////////////////////////

}
