/////////////////////////////////////////
//
// FILE�����
//
/////////////////////////////////////////
#include <stdio.h>

#include "global_extern.h"

void file_close(int RTK)
{
	if(RTK == 1){
		fclose(fp[0]);//rover obs
	}
	fclose(fp[1]);//ref obs
	fclose(fp[2]);//navigation file
	fclose(fp[3]);//�P�Ƒ��ʃt�@�C��
	fclose(fp[4]);//DGNSS�t�@�C��
	fclose(fp[5]);//RTK�t�@�C��
	fclose(fp[7]);//rtkplot�p�t�@�C��
	fclose(fp[9]);//�e�X�g�p�t�@�C��
}