//�������l��p���ăT�C�m�O�����𒆐S�����낦�Ďw�肵���傫���ɐ؂肾���v���O����
//1ch�ŗp�������E�̒l��S�`�����l���ɓK�p����o�[�W����
//�����~�����O���ɂ���ꍇ�̂ݎg�p�\
//�����`�����l���ɑΉ�
//�ech�̓��e�������̋�C�w�ŋK�i�����s���iX���ǂ̏o�͂̎����ϓ��ɑΉ����邽�߁j
//�`�����l�����Ƃ̋K�i�����s��

#define _CRT_SECURE_NO_WARNINGS

#include<iostream>
#include<direct.h>
#include<vector>
#include"matrix_tool.h"

using namespace std;

//�萔
int air_length = 30;//��C�w�̃`�����l����//�K�i���̎��Ɏg�p����
const int abs_l = 262;//�蓮�ŗ��[�����߂�Ƃ��Ɏw�肷�鍶�[
const int abs_r = 773;//�E�[abs_r-abs_l=col_long-1�ƂȂ�悤�ɂ���
bool sel = false;//����:false,�蓮:true
const int row_num = 1024;
const int col_num = 1024;
const int row_up = 450;
const int row_down = 496;
const int row_height = 47;//�c�����͈̔�(d-u+1)
const int col_long = 512;//�T�C�m�O�����̉������̑傫��
const int ch_num = 4;//�`�����l����
const int angle = 180;//���e������
const double thres[ch_num] = { 1000,1000,500,0 };//�������l��C�w�ƃt�@���g���̈�𕪂���ꂻ���ȓK���Ȓl������B1ch�̒l�����g��Ȃ��B
const int pass = 1;//pass�������炻�������[or�E�[

				   //�t�@�C����

const char phantom_name[128] = "phantom_";//�t�@���g���摜�̃t�@�C����
const char outimage_name[128] = "Sinogram_";
const char input_folder_name[128] = "input_sa2.1";
const char output_folder_name[128] = "output_sa2.1";
//const int gain_ave = 3000;//gain�摜�̕��ϒl
int left_1ch[angle];
int right_1ch[angle];

//�v�Z

int min2(int a, int b) {
	int ans = a;
	if (a > b)ans = b;
	return ans;
}
int max2(int a, int b) {
	int ans = a;
	if (a < b)ans = b;
	return ans;
}
void get_bound(int& left, int& right, MATRIX<double>& image, int ch) {
	//�������l��p���ăT�C�m�O�����̍��E���J�b�g����B
	double buf = 0.0;
	vector<double> tmp(col_num + 1, 0.0);
	for (int i = 1;i <= col_num;i++) {
		buf = 0.0;
		for (int j = row_up;j < row_down;j++) {
			buf += image.out(j, i) / row_height;
		}
		tmp[i] = buf;
	}
	bool ravel = true;
	int buf_left = 0;
	int buf_right = 0;
	for (int i = 1;i <= col_num;i++) {
		if (tmp[i] < thres[ch] && i<col_num / 2) {
			buf_left++;
			if (buf_left == pass) {
				left = i;
				ravel = false;
			}
		}
		if (tmp[i] > thres[ch] && i>col_num / 2) {
			buf_right++;
			if (buf_right == pass) {
				right = i;
				ravel = true;
			}
		}
	}
	right = right - pass;
	left = left - pass + 1;
	int size = right - left;
	int add = (col_long - size) / 2;
	left = max2(1, left - add);
	right = min2(col_num, right + add);
	if ((col_long - size) % 2 == 0)right--;
}

int main() {
	char tmp_phantom_name[128];
	char tmp_outimage_name[128];
	int row_long = row_down - row_up;
	MATRIX<double> image;
	MATRIX<double> outimage(angle, col_long);
	for (int r = 1;r <= ch_num;r++) {
		sprintf_s(tmp_outimage_name, "%s%d%s", outimage_name, r, ".raw");
		for (int k = 0;k < angle;k++) {
			_chdir("input_sa2.1");
			sprintf_s(tmp_phantom_name, "%s%d%s%d%s", phantom_name, r, "ch_", k, ".raw");
			//			sprintf_s(tmp_phantom_name, "%s%d%s", phantom_name, k, "after.raw");
			image.read(tmp_phantom_name, row_num, col_num);
			int left, right;
			if (!sel)get_bound(left, right, image, r - 1);//�T�C�m�O�����̗��[�������ŋ��߂�
			if (sel) {
				left = abs_l;
				right = abs_r;
			}
			if (r == 1) {
				left_1ch[k] = left;
				right_1ch[k] = right;
			}
			double air_ave = 0.0;
			for (int i = 0;i < air_length;i++) {
				for (int j = row_up;j < row_down;j++) {
					air_ave += image.out(j, i + left_1ch[k]);
					air_ave += image.out(j, right_1ch[k] - i);
				}
			}
			air_ave /= (row_height*air_length * 2);
			//			image.transpose();//�摜���������Ă�Ƃ��Ɏg��
			for (int i = left_1ch[k];i <= right_1ch[k];i++) {
				double buf = 0.0;
				for (int j = row_up;j < row_down;j++) {
					buf += image.out(j, i) / row_height;
				}

				outimage.in(k + 1, i - left_1ch[k] + 1, buf / air_ave);
			}
		}
		_chdir("../output_sa2.1");
		outimage.write(tmp_outimage_name, "");
		cout << tmp_outimage_name << endl;
		_chdir("..");
	}
	system("pause");
	return 0;
}