//#include<iostream>
//�����֐��p��offset_gain
//short��float��offset��gain���󂯎��
//�����֐����v�Z����O�ɓ��e�f�[�^�̕K�v�̂Ȃ�������؂�����ǂ���������Ȃ��i�������j

#define _CRT_SECURE_NO_WARNINGS
#include<direct.h>
#include"matrix_tool.h"
#include<vector>
#include<fstream>
using namespace std;

//�萔�錾

const int row_num = 1024;
const int col_num = 1024;
const int ch_num = 8;//�`�����l����
const int up_pix = 446;//���ς����s�N�Z���̍s�ԍ��̉����i�摜��ł͏�j
const int down_pix = 492;//���ς����s�N�Z���̍s�ԍ��̏���i�摜��ł͉��j
const int ch_pulse[ch_num] = { 1500,1000,500,0,2000,3500,4000,4500 };//�ech�̃p���X��
//const int angle = 36;//���e������
const char offset_name[128] = "ave_RF_offset_";//�I�t�Z�b�g�摜�̃t�@�C����
const char gain_name[128] = "ave_RF_gain_";//�Q�C���摜�̃t�@�C����
const char phantom_name[128] = "RF_";//�t�@���g���摜�̃t�@�C����
const char outimage_name[128] = "RF_after_";
const int u_num = 1;
const int d_num = 4;
const int ite_num = 10;//�e����̌J��Ԃ���
//const int gain_ave = 3000;//gain�摜�̕��ϒl

//�e�L�X�g�t�@�C���o�͊֘A
double result[ch_num][u_num][d_num];
const char result_filename[128] = "RF_current.txt";

//�v�Z
void solve(int up, int down,int ch,vector<MATRIX<float>> const& offset,vector<MATRIX<float>> const& gain) {
	MATRIX<double> image(row_num,col_num);
	vector<vector<double>> tmp_sum(row_num,vector<double>(col_num, 0.0));
	MATRIX<short> input_image;
	char tmp_phantom_name[128];
	_chdir("input_ogRF");
	for (int i = 0;i < ite_num;i++) {
		//ite_num��B�e�̕��ς����B
		sprintf(tmp_phantom_name, "%s%d%s%d%s%d%s%d%s",phantom_name,ch_pulse[ch],"p_",up,"u_",down,"d_",i,".raw");
		input_image.read(tmp_phantom_name, row_num, col_num);
		for (int row = 0;row < row_num;row++) {
			for (int col = 0;col < col_num;col++) {
				tmp_sum[row][col] += (double)input_image.out(row + 1, col + 1);
			}
		}
	}
	for (int i = 0;i < row_num;i++) {
		for (int j = 0;j < col_num;j++) {
			double buf;
			//offset_gain�␳�v�Z
			buf = (tmp_sum[i][j]/(double)ite_num - offset[ch].out(i + 1, j + 1))/(gain[ch].out(i+1,j+1)-offset[ch].out(i+1,j+1));
			image.in(i + 1, j + 1, buf);
		}
	}
	_chdir("../output_ogRF");
	char tmp_outimage_name[128];
	double gain_ave = 0.0;
	result[ch][up][down] = 0.0;
	int num = (down_pix - up_pix + 1)*col_num;
	//gain_ave�̌v�Z��offset_gain�̌v�Z���ʂ܂ł��o��
	for (int i = up_pix;i <= down_pix;i++) {
		for (int j = 1;j <= col_num;j++) {
			gain_ave += gain[ch].out(i, j) - offset[ch].out(i,j);
			result[ch][up][down] += image.out(i, j);
		}
	}
	image = image * gain_ave / (double)num;
	sprintf(tmp_outimage_name, "%s%d%s%d%s%d%s", outimage_name, ch+1, "ch_", up, "u_", down, "d_.raw");
	image.write(tmp_outimage_name, "");
	_chdir("..");
	result[ch][up][down] = result[ch][up][down] * gain_ave/((double)num*(double)num);
}

int main() {
	vector<MATRIX<float>> offset(ch_num);
	vector<MATRIX<float>> gain(ch_num);
	//ch������offset��gain���擾
	for (int ch = 0;ch < ch_num;ch++) {
		char tmp_offset_name[128];
		char tmp_gain_name[128];
		sprintf(tmp_offset_name, "%s%d%s", offset_name, ch_pulse[ch], "p_0.raw");
		sprintf(tmp_gain_name, "%s%d%s", gain_name, ch_pulse[ch], "p_0.raw");
		_chdir("offset");
		offset[ch].read(tmp_offset_name, row_num, col_num);
		_chdir("../gain");
		gain[ch].read(tmp_gain_name, row_num, col_num);
		_chdir("..");
	}
	//�v�Z����offset_gain�␳��̉摜���o�͂���
	for (int i = 0;i < u_num;i++) {
		for (int j = 0;j < d_num;j++) {
			for (int ch = 0;ch < ch_num;ch++) {
				solve(i, j, ch,offset,gain);
				printf("%d�i�ڂ�%dch����\r", j, ch);
			}
		}
	}
	_chdir("result_ogRF");
	//�����֐��p�t�@���g�����B�e������̕��ϓd���l���o�͂���
	ofstream result_file(result_filename);
	for (int u = 0;u < u_num;u++) {
		for (int ch = 0;ch < ch_num;ch++) {
			for (int d = 0;d < d_num;d++) {
				result_file << result[ch][u][d] << endl;
			}
		}
	}
	return 0;
}