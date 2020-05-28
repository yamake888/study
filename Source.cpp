//���m�t�@���g���̓��e�d���l���擾����v���O����
//�~���t�@���g���̂ݑΉ�
//phantom�͎��ʌ���W��*������W���̒l���G�l���M�[���Ƃɂ���͂���
//�t�@���g���ʉߌ�̃T�C�m�O�������o�͂���
//4���o�[�W����4���̂ݑΉ�
//�e�`�����l���̋�C�w�ŋK�i��

#define _CRT_SECURE_NO_WARNINGS

#include<iostream>
#include<direct.h>
#include"matrix_tool.h"

using namespace std;

//�萔

//#define pixSize 0.004167	//�P��[cm]
//const double pixSize=0.01;	//�R�����[�g����1step
const double pix_size = 0.00488;//1pixel������̑傫��(cm)
							   //#define pixSize 0.00384	//�P��[cm]//���o��100cm,�팟��80cm�̈ʒu�ɂ������Ƃ��B
const double radius_outer = 1.5;	//�~���O�a(���a)�@�P��[cm]
const double radius_inner = 1.5;	//�~�����a(���a)�@�P��[cm]
									//#define radius_outer 0.3	//�~���O�a(���a)�@�P��[cm]
									//#define radius_inner 0.27	//�~�����a(���a)�@�P��[cm]
									//#define thickness_max 0.03	//�ő����//�K�i�t�@���g���̏ꍇ[cm]//�ŏ���0�Ƃ���B
const double thickness_step = 0.01;	//�����̃X�e�b�v[cm]//�K�i�t�@���g���̏ꍇ
									//#define pixNum 30	//�����̃X�e�b�v��
const int pix_num = 1024;	//�������s�N�Z����(���o��T�C�Y)���~���t�@���g���̏ꍇ�̓R�����g�A�E�g���͂������ƁB;
							//#define EnergyBin 199	//100kV�̏ꍇ
const int angle_num = 180;//���e������
const double angle_step = 1.0;//1step�ł̊p�x�̕ω���//angle_step*angle_num=180�̎��n�[�t�X�L����
const int energy_bin = 279;	//140kV�̏ꍇ//0.5keV����
							//const int EnergyBin = 140;	//140kV�̏ꍇ//1keV����
							//const int EnergyBin = 70;	//140kV�̏ꍇ//2keV����
							//const int EnergyBin = 28;	//140kV�̏ꍇ//5keV����
const int ch_num = 4;
const int scale = 100000;


//�K�i��t�@���g���̏ꍇ�A�P�����̏ꍇ��phantom_outer���̗p�B

const char outer_file_name[20] = "PVC(C2H3Cl).txt";
const char inner_file_name[30] = "PVC(C2H3Cl).txt";

//���̒��ɓ����t�@���g���̏��//���̏��͎g���Ȃ�
const char phantom1_name[128] = "PMMA.txt";//�E
const char phantom2_name[128] = "PE.txt";//��
const char phantom3_name[128] = "POM.txt";//��
const char phantom4_name[128] = "PTFE.txt";//���̕���
//�܂Ƃ߂����
const char phantom_name[128] = "phantom_b.txt";

double phantom_pcd = 0.75;//�Ίp�̃t�@���g�����m�̋����̔���(cm)
double phantom_radius = 0.25;//�t�@���g���̔��a(cm)

const double outer_d = 0.0;//Cu:8.96;//�O���̕����̖��x
const double inner_d = 1.38;//PMMA:1.19;//�����̕����̖��x

const char RF_name[64] = "RF_140kV_Cu_0,0.1,0.2,0.3_4ch.txt";
const char Spectrum_file_name[64] = "TRIX140kV.txt";
const char output_Sinogram_file_name[64] = "PVC_phantom_b_";//�o�̓t�@�C����


 //�v�Z

int main() {
	MATRIX<double> Spectrum;
	MATRIX<double> RF;
	MATRIX<double> inner_attenuation;
	MATRIX<double> outer_attenuation;
	MATRIX<double> output_Sinogram(angle_num, pix_num);
	MATRIX<double> phantom;
	_chdir("input");

	Spectrum.read(Spectrum_file_name, energy_bin, 1);
	RF.read(RF_name, energy_bin, ch_num);
	inner_attenuation.read(inner_file_name, energy_bin, 1);
	outer_attenuation.read(outer_file_name, energy_bin, 1);
	phantom.read(phantom_name, energy_bin, 4);

	//���ʌ���W���ɖ��x�������Đ�����W���ɂ���
	inner_attenuation = inner_attenuation * inner_d;
	outer_attenuation = outer_attenuation * outer_d;

	double iso_center = pix_num * pix_size / 2;//���̂̒��S�̍��[����̈ʒu(cm)
	//double t_step = 0;//�ŏ�����i�񂾋���(cm)
	double outer_left = iso_center - radius_outer;
	double outer_right = iso_center + radius_outer;
	double inner_left = iso_center - radius_inner;
	double inner_right = iso_center + radius_inner;
	double outer_dist;
	double inner_dist;
	char tmp_output_file_name[128];
	double norm = 1;

	//���ɓ����t�@���g���̏��
	double p_center[4];
	double p_dist[4] = { 0.0,0.0,0.0,0.0 };
	double p_left[4];
	double p_right[4];
	//double p_attenuation[4];


	_chdir("../output");

	for (int n = 1;n <= ch_num;n++) {
		sprintf_s(tmp_output_file_name, "%s%d%s", output_Sinogram_file_name, n, ".raw");
		double phantom_angle[4] = { 0.0,90.0,180.0,270.0 };//�t�@���g���̒��S����̏����p�x(��)
		for (int rad = 1;rad <= angle_num;rad++) {
			double t_step = 0.0;//�ŏ�����i�񂾋���
			//�}���t�@���g���̈ʒu���v�Z
			for (int a = 0;a < 4;a++) {
				p_center[a] = iso_center + phantom_pcd * cos(3.14159265358979 * phantom_angle[a] / 180.0);
				p_left[a] = p_center[a] - phantom_radius;
				p_right[a] = p_center[a] + phantom_radius;
				//cout << phantom_angle[a] << endl;
				phantom_angle[a] += angle_step;
			}
			for (int i = 0;i < pix_num;i++) {
				t_step += pix_size;
				double buf = 0;
				if (t_step > outer_left && t_step < outer_right) {
					//�O��������ʉ߂��鋗�����v�Z
					outer_dist = 2 * sqrt(radius_outer * radius_outer - (t_step - iso_center)*(t_step - iso_center));
				}
				else outer_dist = 0;
				if (t_step > inner_left && t_step < inner_right) {
					//����������ʉ߂��鋗�����v�Z
					inner_dist = 2 * sqrt(radius_inner * radius_inner - (t_step - iso_center)*(t_step - iso_center));
				}
				else inner_dist = 0;
				//�e�t�@���g����ʉ߂��鋗�����v�Z
				for (int b = 0;b < 4;b++) {
					if (t_step > p_left[b] && t_step < p_right[b]) {
						p_dist[b] = 2 * sqrt(phantom_radius * phantom_radius - (t_step - p_center[b])*(t_step - p_center[b]));
					}
					else p_dist[b] = 0;
				}
				for (int k = 1;k <= energy_bin;k++) {
					//����v�Z
					buf += Spectrum.out(k, 1)*RF.out(k, n)*exp(-outer_attenuation.out(k, 1)*(outer_dist - inner_dist) - inner_attenuation.out(k, 1)*(inner_dist - p_dist[0] - p_dist[1] - p_dist[2] - p_dist[3]) - phantom.out(k, 1)*p_dist[0] - phantom.out(k, 2)*p_dist[1] - phantom.out(k, 3)*p_dist[2] - phantom.out(k, 4)*p_dist[3]);
				}
				output_Sinogram.in(rad, i + 1, buf);
				printf_s("%d����%d�X�e�b�v�ڂ̓d���l�̌v�Z����\r", rad, i);
			}
		}
		output_Sinogram = output_Sinogram/output_Sinogram.out(1,1);
		output_Sinogram.write(tmp_output_file_name, "");
	}
//	double test;
//	test = cos(3.14159265 / 3+6.28318);
//	cout << test << endl;
	system("pause");
	return 0;
}

