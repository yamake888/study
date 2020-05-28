//既知ファントムの投影電流値を取得するプログラム
//円筒ファントムのみ対応
//phantomは質量減弱係数*線減弱係数の値をエネルギーごとにを入力する
//ファントム通過後のサイノグラムを出力する
//4つ穴バージョン4つ穴のみ対応
//各チャンネルの空気層で規格化

#define _CRT_SECURE_NO_WARNINGS

#include<iostream>
#include<direct.h>
#include"matrix_tool.h"

using namespace std;

//定数

//#define pixSize 0.004167	//単位[cm]
//const double pixSize=0.01;	//コリメート実験1step
const double pix_size = 0.00488;//1pixelあたりの大きさ(cm)
							   //#define pixSize 0.00384	//単位[cm]//検出器100cm,被検体80cmの位置においたとき。
const double radius_outer = 1.5;	//円筒外径(半径)　単位[cm]
const double radius_inner = 1.5;	//円筒内径(半径)　単位[cm]
									//#define radius_outer 0.3	//円筒外径(半径)　単位[cm]
									//#define radius_inner 0.27	//円筒内径(半径)　単位[cm]
									//#define thickness_max 0.03	//最大厚さ//階段ファントムの場合[cm]//最小は0とする。
const double thickness_step = 0.01;	//厚さのステップ[cm]//階段ファントムの場合
									//#define pixNum 30	//厚さのステップ数
const int pix_num = 1024;	//得たいピクセル数(検出器サイズ)●円筒ファントムの場合はコメントアウトをはずすこと。;
							//#define EnergyBin 199	//100kVの場合
const int angle_num = 180;//投影方向数
const double angle_step = 1.0;//1stepでの角度の変化量//angle_step*angle_num=180の時ハーフスキャン
const int energy_bin = 279;	//140kVの場合//0.5keV刻み
							//const int EnergyBin = 140;	//140kVの場合//1keV刻み
							//const int EnergyBin = 70;	//140kVの場合//2keV刻み
							//const int EnergyBin = 28;	//140kVの場合//5keV刻み
const int ch_num = 4;
const int scale = 100000;


//階段状ファントムの場合、１物質の場合はphantom_outerを採用。

const char outer_file_name[20] = "PVC(C2H3Cl).txt";
const char inner_file_name[30] = "PVC(C2H3Cl).txt";

//穴の中に入れるファントムの情報//この情報は使われない
const char phantom1_name[128] = "PMMA.txt";//右
const char phantom2_name[128] = "PE.txt";//上
const char phantom3_name[128] = "POM.txt";//左
const char phantom4_name[128] = "PTFE.txt";//下の物質
//まとめたやつ
const char phantom_name[128] = "phantom_b.txt";

double phantom_pcd = 0.75;//対角のファントム同士の距離の半分(cm)
double phantom_radius = 0.25;//ファントムの半径(cm)

const double outer_d = 0.0;//Cu:8.96;//外側の物質の密度
const double inner_d = 1.38;//PMMA:1.19;//内側の物質の密度

const char RF_name[64] = "RF_140kV_Cu_0,0.1,0.2,0.3_4ch.txt";
const char Spectrum_file_name[64] = "TRIX140kV.txt";
const char output_Sinogram_file_name[64] = "PVC_phantom_b_";//出力ファイル名


 //計算

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

	//質量減弱係数に密度をかけて線減弱係数にする
	inner_attenuation = inner_attenuation * inner_d;
	outer_attenuation = outer_attenuation * outer_d;

	double iso_center = pix_num * pix_size / 2;//物体の中心の左端からの位置(cm)
	//double t_step = 0;//最初から進んだ距離(cm)
	double outer_left = iso_center - radius_outer;
	double outer_right = iso_center + radius_outer;
	double inner_left = iso_center - radius_inner;
	double inner_right = iso_center + radius_inner;
	double outer_dist;
	double inner_dist;
	char tmp_output_file_name[128];
	double norm = 1;

	//穴に入れるファントムの情報
	double p_center[4];
	double p_dist[4] = { 0.0,0.0,0.0,0.0 };
	double p_left[4];
	double p_right[4];
	//double p_attenuation[4];


	_chdir("../output");

	for (int n = 1;n <= ch_num;n++) {
		sprintf_s(tmp_output_file_name, "%s%d%s", output_Sinogram_file_name, n, ".raw");
		double phantom_angle[4] = { 0.0,90.0,180.0,270.0 };//ファントムの中心からの初期角度(°)
		for (int rad = 1;rad <= angle_num;rad++) {
			double t_step = 0.0;//最初から進んだ距離
			//挿入ファントムの位置を計算
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
					//外側物質を通過する距離を計算
					outer_dist = 2 * sqrt(radius_outer * radius_outer - (t_step - iso_center)*(t_step - iso_center));
				}
				else outer_dist = 0;
				if (t_step > inner_left && t_step < inner_right) {
					//内側物質を通過する距離を計算
					inner_dist = 2 * sqrt(radius_inner * radius_inner - (t_step - iso_center)*(t_step - iso_center));
				}
				else inner_dist = 0;
				//各ファントムを通過する距離を計算
				for (int b = 0;b < 4;b++) {
					if (t_step > p_left[b] && t_step < p_right[b]) {
						p_dist[b] = 2 * sqrt(phantom_radius * phantom_radius - (t_step - p_center[b])*(t_step - p_center[b]));
					}
					else p_dist[b] = 0;
				}
				for (int k = 1;k <= energy_bin;k++) {
					//減弱計算
					buf += Spectrum.out(k, 1)*RF.out(k, n)*exp(-outer_attenuation.out(k, 1)*(outer_dist - inner_dist) - inner_attenuation.out(k, 1)*(inner_dist - p_dist[0] - p_dist[1] - p_dist[2] - p_dist[3]) - phantom.out(k, 1)*p_dist[0] - phantom.out(k, 2)*p_dist[1] - phantom.out(k, 3)*p_dist[2] - phantom.out(k, 4)*p_dist[3]);
				}
				output_Sinogram.in(rad, i + 1, buf);
				printf_s("%d方向%dステップ目の電流値の計算完了\r", rad, i);
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

