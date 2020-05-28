//#include <stdio.h>
//#include <windows.h>//　matrix_toolと同じヘッダーがあると動かない場合がある。
//複数チャンネルに対応 2019/07/30
//各列の計算において最初に各エネルギーごとに鉄の減弱を考慮しない電流値を確保することで計算の高速化を行った
//またあらかじめ鉄による減弱を各鉄厚さ*エネルギー数で求めておくことでまだ早くできそう(未実装)
//LUTを受け取ってスペクトル分解を行うプログラム2019/07/30
//LUTは各chの空気層で規格化しておくこと
//サイノグラムの規格化はここで行う
//二分探索による高速化が可能
//#define _CRT_SECURE_NO_WARNINGS

#include "matrix_tool.h"
#include <iostream>
#include<direct.h>
#include<fstream>
using namespace std;

//LUT取得に必要な定数
const int row_num = 400;//縦方向のピクセル数
const int col_num = 30000;//横方向のピクセル数
const double row_t_step = 0.014;//列方向物質の厚さの刻み幅(cm)
const double col_t_step = 0.0002;//行方向物質の厚さの刻み幅(cm)
const int energy_bin = 279;//エネルギービンの数
const double bin_size = 0.5;//エネルギービンの大きさ
const char lac_row_filename[256] = "B5.txt";//質量減弱係数
const char lac_col_filename[256] = "Si14.txt";//"SUS304(Fe74%Ni8%Cr18%).txt";//質量減弱係数
const char RF_filename[256] = "RF_140kV_Cu_0,0.1,0.2,0.3_4ch.txt";
const char Spectrum_filename[30] = "TRIX140kV.txt";
const int ch_num = 4;
const char LUT_file_name[256] = "LUT_(30000,400)";
const double density_row = 0.930;//列方向物質の密度
const double density_col = 1.4;//行方向物質の密度
const int n_len = 30;//規格化する時に両端から足し合わせる長さ

								 //物質分解してスペクトルを取得するために追加で必要な定数
const char input_Sinogram_file_name[64] = "PVC_phantom_b_";//サイノグラムの入力ファイル名
const int angle = 180;
const int pix_num = 1024;

double LUT[ch_num][row_num][col_num];
double Sinogram[ch_num][angle][pix_num];
const char er_sino_name[128] = "er_sino.txt";
const char output_filename_t[128] = "thickness.txt";
const char output_row_filename[128] = "B_Sino.raw";
const char output_col_filename[128] = "Si_Sino.raw";


bool need_norm = false;//ここで規格化するかどうか.true:する,false:しない

int max2(int const& a, int const& b) {
	int ans = a;
	if (b > a)ans = b;
	return ans;
}
int min2(int const& a, int const& b) {
	int ans = a;
	if (b < a)ans = b;
	return ans;
}

int main() {
	_chdir("input_gs");
	//MATRIXを宣言&取得
	MATRIX<double> Spectrum;
	Spectrum.read(Spectrum_filename, energy_bin, 1);
	MATRIX<double> RFunction;
	RFunction.read(RF_filename, energy_bin, ch_num);
	MATRIX<double> lac_row;
	lac_row.read(lac_row_filename, energy_bin, 1);
	MATRIX<double> lac_col;
	lac_col.read(lac_col_filename, energy_bin, 1);

	//サイノグラムを取得

	char tmp_Sino_filename[128];
	MATRIX<double> Sino[ch_num];//サイノグラムを保存する配列
	for (int ch = 0;ch < ch_num;ch++) {
		sprintf_s(tmp_Sino_filename, "%s%d%s", input_Sinogram_file_name, ch + 1, ".raw");
		Sino[ch].read(tmp_Sino_filename, angle, pix_num);
		if(need_norm)Sino[ch] = Sino[ch]/Sino[ch].norm(n_len);
	}

	int tmp_row_num = row_num + 1;
	int tmp_col_num = col_num + 1;
	int tmp_energy_bin = energy_bin + 1;

	//計算

	lac_row = lac_row * density_row;
	lac_col = lac_col * density_col;

	//LUT取得

	MATRIX<double> LUT_array[ch_num];
	char output_filename[256];
	for (int ch = 0;ch < ch_num;ch++) {
		sprintf_s(output_filename, "%s%d%s", LUT_file_name, ch + 1, "ch.raw");
		LUT_array[ch].read(output_filename, row_num, col_num);
	}

	_chdir("../output_gs");//outputに移動
	double tmp_row_t = 0;
	double tmp_col_t = 0;
	double sd = 0;
	double sum = 0;
	double tmp = 0;
	ofstream output_file(er_sino_name);
	ofstream output_t(output_filename_t);

	//MATRIXから配列に入れなおす
	MATRIX<double> Material_row(angle, pix_num);
	MATRIX<double> Material_col(angle, pix_num);


	for (int ch = 0;ch < ch_num;ch++) {
		for (int i = 0;i < row_num;i++) {
			for (int j = 0;j < col_num;j++) {
				LUT[ch][i][j] = LUT_array[ch].out(i + 1, j + 1);
			}
		}
		for (int i = 0;i < angle;i++) {
			for (int j = 0;j < pix_num;j++) {
				Sinogram[ch][i][j] = Sino[ch].out(i + 1, j + 1);
			}
		}
	}

	for (int y = 0;y < angle;y++) {
		for (int x = 0;x < pix_num;x++) {
			double min = 1000;
			//LUTから最適な減弱物質の厚さのペアを求める
			//二分探索による高速化が可能
			/*
			//全探索
			for (int i = 0;i < row_num;i++) {
			for (int j = 0;j < col_num;j++) {
			sum = 0;
			for (int ch = 0;ch < ch_num;ch++) {
			tmp = (-LUT[ch][i][j] + Sinogram[ch][y][x]);
			sum += tmp * tmp;
			}
			if (min > sum) {
			tmp_row_t = i;
			tmp_col_t = j;
			sd = sum;
			min = sum;
			}
			}
			//				}
			}
			//全探索法終わり
			*/
			//二分探索法を用いて高速に行った
			for (int i = 0;i < row_num;i++) {
				int left = 0;
				int right = col_num - 1;
				int j = (left + right) / 2;
				double tmp1;
				double tmp2;
				double sum1;
				double sum2;
				while (right - left > 2) {
					sum1 = 0;
					sum2 = 0;
					j = (left + right) / 2;
					for (int ch = 0;ch < ch_num;ch++) {
						tmp1 = (-LUT[ch][i][j] + Sinogram[ch][y][x]);
						sum1 += tmp1 * tmp1;
						tmp2 = (-LUT[ch][i][j + 1] + Sinogram[ch][y][x]);
						sum2 += tmp2 * tmp2;
					}
					if (sum2 > sum1)right = j;
					else left = j;
				}
				//二分探索ここまで

				//念のため前後conf個を調べる
				int conf = 50;
				int lower = max2(0, left - conf);
				int sup = min2(right + conf, col_num - 1);
				for (int f = lower;f < sup;f++) {
					sum = 0;
					for (int ch = 0;ch < ch_num;ch++) {
						tmp = (-LUT[ch][i][f] + Sinogram[ch][y][x]);
						sum += tmp * tmp;
					}
					if (min > sum) {
						tmp_row_t = i;
						tmp_col_t = f;
						sd = sum;
						min = sum;
					}
				}
			}
			//二分探索を使った探索はここまで

			double buf;
			Material_row.in(y+1, x+1, tmp_row_t*row_t_step);
			Material_col.in(y + 1, x + 1, tmp_col_t*col_t_step);

			for (int k = 1;k <= energy_bin;k++) {
				buf = Spectrum.out(k, 1)*exp(-lac_row.out(k, 1)*tmp_row_t*row_t_step - lac_col.out(k, 1)*tmp_col_t*col_t_step);
				//				buf = 1;

				output_file << buf << "\t";
			}
			output_file << endl;
			output_t << tmp_row_t * row_t_step << "\t" << tmp_col_t * col_t_step << "\t" << sd << endl;
			printf_s("%d方向目%d番目完了\r", y, x);
		}
	}
	output_file.close();
	output_t.close();
	Material_row.write(output_row_filename,"");
	Material_col.write(output_col_filename, "");
	return 0;
}

