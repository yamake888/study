//#include<iostream>
#define _CRT_SECURE_NO_WARNINGS
#include<direct.h>
#include"matrix_tool.h"
using namespace std;

//定数宣言

const int row_num = 1024;
const int col_num = 1024;
const int ch_num = 8;//チャンネル数
const int angle = 1;//投影方向数
const int sup = 450;
const int inf = 496;
const char offset_name[128] = "ave_PVDF_offset_";//オフセット画像のファイル名
const char gain_name[128] = "ave_PVDF_gain_";//ゲイン画像のファイル名
const char phantom_name[128] = "ave_PVDF_";//ファントム画像のファイル名
const int ch_pulse[ch_num] = { 1500,1000,500,0,2000,3500,4000,4500 };
const char outimage_name[128] = "phantom_";
const char input_folder_name[128] = "input_og4";
const char output_folder_name[128] = "output_og4";
//const int gain_ave = 3000;//gain画像の平均値

//計算
double ave(MATRIX<float> const& gain) {
	//supとinfの間の領域でgainの平均を返す関数
	double sum = 0.0;
	for (int i = sup;i <= inf;i++) {
		for (int j = 1;j <= col_num;j++) {
			sum += gain.out(i, j);
		}
	}
	double num = (double)(inf - sup + 1)*(double)col_num;
	sum = sum / num;
	return sum;
}
int main() {
	char tmp_phantom_name[128];
	char tmp_outimage_name[128];
	MATRIX<float> image[ch_num];
	MATRIX<float> offset[ch_num];
	MATRIX<float> gain[ch_num];
	MATRIX<double> outimage(row_num, col_num);
	double gain_ave[ch_num];
	_chdir("offset_gain4");
	for (int ch = 0;ch < ch_num;ch++) {
		char tmp_gain_filename[128];
		char tmp_offset_filename[128];
		sprintf(tmp_gain_filename, "%s%d%s", gain_name, ch_pulse[ch], "p_0.raw");
		sprintf(tmp_offset_filename, "%s%d%s", offset_name, ch_pulse[ch], "p_0.raw");
		gain[ch].read(tmp_gain_filename, row_num, col_num);
		offset[ch].read(tmp_offset_filename,row_num,col_num);
		gain_ave[ch] = ave(gain[ch]);
	}
	_chdir("..");
	int tmp_row_num = row_num + 1;
	int tmp_col_num = col_num + 1;
	for (int j = 0;j < ch_num;j++) {
		for (int i = 0;i < angle;i++) {
				_chdir(input_folder_name);
				sprintf_s(tmp_phantom_name, "%s%d%s%d%s", phantom_name, ch_pulse[j],"p_",i, ".raw");//複数チャンネルの場合
//				sprintf_s(tmp_phantom_name, "%s%d%s", phantom_name, i, ".raw");//シングルチャンネルの場合
				sprintf_s(tmp_outimage_name, "%s%d%s%d%s", outimage_name, j+1, "ch_", i, ".raw");
				image[j].read(tmp_phantom_name, row_num, col_num);
				double buf;
				for (int a = 1;a < tmp_row_num;a++) {
					for (int b = 1;b < tmp_col_num;b++) {
						buf = (image[j].out(a, b) - offset[j].out(a, b)) * gain_ave[j] / (gain[j].out(a, b) - offset[j].out(a, b));
						outimage.in(a, b, buf);
					}
				}
				_chdir("../output_og4");
				outimage.write(tmp_outimage_name, "");
				printf("%s\n", tmp_outimage_name);
				_chdir("..");
		}
	}
	system("pause");
	return 0;
}