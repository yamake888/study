//#include <stdio.h>
//#include <windows.h>//�@matrix_tool�Ɠ����w�b�_�[������Ɠ����Ȃ��ꍇ������B
//�����`�����l���ɑΉ� 2019/07/30
//�e��̌v�Z�ɂ����čŏ��Ɋe�G�l���M�[���ƂɓS�̌�����l�����Ȃ��d���l���m�ۂ��邱�ƂŌv�Z�̍��������s����
//�܂����炩���ߓS�ɂ�錸����e�S����*�G�l���M�[���ŋ��߂Ă������Ƃł܂������ł�����(������)
//LUT���󂯎���ăX�y�N�g���������s���v���O����2019/07/30
//LUT�͊ech�̋�C�w�ŋK�i�����Ă�������
//�T�C�m�O�����̋K�i���͂����ōs��
//�񕪒T���ɂ�鍂�������\
//#define _CRT_SECURE_NO_WARNINGS

#include "matrix_tool.h"
#include <iostream>
#include<direct.h>
#include<fstream>
using namespace std;

//LUT�擾�ɕK�v�Ȓ萔
const int row_num = 400;//�c�����̃s�N�Z����
const int col_num = 30000;//�������̃s�N�Z����
const double row_t_step = 0.014;//����������̌����̍��ݕ�(cm)
const double col_t_step = 0.0002;//�s���������̌����̍��ݕ�(cm)
const int energy_bin = 279;//�G�l���M�[�r���̐�
const double bin_size = 0.5;//�G�l���M�[�r���̑傫��
const char lac_row_filename[256] = "B5.txt";//���ʌ���W��
const char lac_col_filename[256] = "Si14.txt";//"SUS304(Fe74%Ni8%Cr18%).txt";//���ʌ���W��
const char RF_filename[256] = "RF_140kV_Cu_0,0.1,0.2,0.3_4ch.txt";
const char Spectrum_filename[30] = "TRIX140kV.txt";
const int ch_num = 4;
const char LUT_file_name[256] = "LUT_(30000,400)";
const double density_row = 0.930;//����������̖��x
const double density_col = 1.4;//�s���������̖��x
const int n_len = 30;//�K�i�����鎞�ɗ��[���瑫�����킹�钷��

								 //�����������ăX�y�N�g�����擾���邽�߂ɒǉ��ŕK�v�Ȓ萔
const char input_Sinogram_file_name[64] = "PVC_phantom_b_";//�T�C�m�O�����̓��̓t�@�C����
const int angle = 180;
const int pix_num = 1024;

double LUT[ch_num][row_num][col_num];
double Sinogram[ch_num][angle][pix_num];
const char er_sino_name[128] = "er_sino.txt";
const char output_filename_t[128] = "thickness.txt";
const char output_row_filename[128] = "B_Sino.raw";
const char output_col_filename[128] = "Si_Sino.raw";


bool need_norm = false;//�����ŋK�i�����邩�ǂ���.true:����,false:���Ȃ�

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
	//MATRIX��錾&�擾
	MATRIX<double> Spectrum;
	Spectrum.read(Spectrum_filename, energy_bin, 1);
	MATRIX<double> RFunction;
	RFunction.read(RF_filename, energy_bin, ch_num);
	MATRIX<double> lac_row;
	lac_row.read(lac_row_filename, energy_bin, 1);
	MATRIX<double> lac_col;
	lac_col.read(lac_col_filename, energy_bin, 1);

	//�T�C�m�O�������擾

	char tmp_Sino_filename[128];
	MATRIX<double> Sino[ch_num];//�T�C�m�O������ۑ�����z��
	for (int ch = 0;ch < ch_num;ch++) {
		sprintf_s(tmp_Sino_filename, "%s%d%s", input_Sinogram_file_name, ch + 1, ".raw");
		Sino[ch].read(tmp_Sino_filename, angle, pix_num);
		if(need_norm)Sino[ch] = Sino[ch]/Sino[ch].norm(n_len);
	}

	int tmp_row_num = row_num + 1;
	int tmp_col_num = col_num + 1;
	int tmp_energy_bin = energy_bin + 1;

	//�v�Z

	lac_row = lac_row * density_row;
	lac_col = lac_col * density_col;

	//LUT�擾

	MATRIX<double> LUT_array[ch_num];
	char output_filename[256];
	for (int ch = 0;ch < ch_num;ch++) {
		sprintf_s(output_filename, "%s%d%s", LUT_file_name, ch + 1, "ch.raw");
		LUT_array[ch].read(output_filename, row_num, col_num);
	}

	_chdir("../output_gs");//output�Ɉړ�
	double tmp_row_t = 0;
	double tmp_col_t = 0;
	double sd = 0;
	double sum = 0;
	double tmp = 0;
	ofstream output_file(er_sino_name);
	ofstream output_t(output_filename_t);

	//MATRIX����z��ɓ���Ȃ���
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
			//LUT����œK�Ȍ��㕨���̌����̃y�A�����߂�
			//�񕪒T���ɂ�鍂�������\
			/*
			//�S�T��
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
			//�S�T���@�I���
			*/
			//�񕪒T���@��p���č����ɍs����
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
				//�񕪒T�������܂�

				//�O�̂��ߑO��conf�𒲂ׂ�
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
			//�񕪒T�����g�����T���͂����܂�

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
			printf_s("%d������%d�Ԗڊ���\r", y, x);
		}
	}
	output_file.close();
	output_t.close();
	Material_row.write(output_row_filename,"");
	Material_col.write(output_col_filename, "");
	return 0;
}

