
#include "matrix_tool.h"
#include <direct.h>

void LSM_for_Zeff(MATRIX<double>, MATRIX<double>, double&, double&);//最小二乗法によって、係数項と誤差を返す関数(臨時)
void est_Z_least_square_method(MATRIX<double>, MATRIX<double>, double&, double&);
void est_Z_least_square_method_LN(MATRIX<double>, MATRIX<double>, double&, double&);
void est_Z_minimize_max_relative_error(MATRIX<double>, MATRIX<double>, double&, double&);
MATRIX<double> LUT_resize_Z1_30_281_2902(MATRIX<double>, MATRIX<double>);
void data_output(char*, double, double, double, MATRIX<double>, MATRIX<double>);

void main() {
//	_chdir("work");
	const int data_num = 120;/////////////////

	int i;
	MATRIX<double> buf;

	MATRIX<double> input_data;
	input_data.read("output60to120.txt", data_num, 2);

	MATRIX<double> Zeff_table;
//	Zeff_table.read("Zeff_LUT_Z1_30(281,2901+1).raw", 281, 2902);
	Zeff_table.read("Zeff_LUT281_to_2902.raw", 281, 2902);
	MATRIX<double> LUT;
	LUT = LUT_resize_Z1_30_281_2902(Zeff_table, input_data.pullcol(1));
	LUT.delcol(1);

	double res;
	double Neff_buf;
	double res_buf;
	double Neff;
	double Zeff;
	MATRIX<double> LUT_buf;

	res = 1.797693e+308;
	for (i = 0;i<2901;i++) {
		LUT_buf = LUT.pullcol(i + 1);

		//LSM_for_Zeff(input_data.pullcol(2),LUT_buf,Neff_buf,res_buf);
		est_Z_least_square_method(input_data.pullcol(2), LUT_buf, Neff_buf, res_buf);
		//est_Z_least_square_method_LN(input_data.pullcol(2),LUT_buf,Neff_buf,res_buf);
		//est_Z_minimize_max_relative_error(input_data.pullcol(2),LUT_buf,Neff_buf,res_buf);

		if (res>res_buf) {
			Zeff = 1.0 + 0.01*i;
			Neff = Neff_buf;
			res = res_buf;
		}
	}

	data_output("result60to120.txt", Zeff, Neff, res, input_data, LUT);

//	system("pause");
}

void LSM_for_Zeff(MATRIX<double> sample, MATRIX<double> comparison, double &Neff, double &err)
//引数のNeffおよびerrは値を受け取るための数値で、関数の戻り値の代わりに使う。
//最小二乗法によって、sampleの線減弱係数に比例係数Cを掛けた値と、comparizonの微分断面積の誤差errが最小となるCを求める。
//横軸エネルギー(線形)、縦軸微分断面積(対数)のグラフ上に、sampleの線減弱係数を乗せると考えればよい。
//最小化については、対数によって判定されるようにした。
//求める条件は f=Σ{Y-(Z+X))^2 について f'=0となる　C=exp(Z)である。
//エラーerrは決定係数R^2=1-errのerrであり、0以上の最小の値が最も良い値となる。
{
	if (sample.Row() != comparison.Row()) {
		printf("Size error! (LMS_for_Zeff(row=%d,row=%d)\n", sample.Row(), comparison.Row());
		system("pause");
		exit(EXIT_FAILURE);
	}


	MATRIX<double> mbuf;

	double err_deno, err_nume;//それぞれerrの分母と分子。

	double lnC;
	double C;

	mbuf = log(comparison) - log(sample);//対数をとった上で(err)'=0を解くとこのようにlnCを表せる。
	lnC = mbuf.sum() / mbuf.Row();
	C = exp(lnC);

	mbuf = log(comparison) - log(C*sample);
	mbuf = mbuf * mbuf;
	err_nume = mbuf.sum();

	mbuf = log(comparison);//12,Dec.,2017 hamaguchi この部分は決定係数の定義に合っていない。間違いである。
	mbuf = log(comparison) - mbuf.average();
	mbuf = mbuf * mbuf;
	err_deno = mbuf.sum();

	err = err_nume / err_deno;
	Neff = 1.0 / C;
}

void est_Z_least_square_method(MATRIX<double> LAC_data, MATRIX<double> LUT_data, double &N, double &res)
//引数のNおよびresは値を受け取るための数値で、関数の戻り値の代わりに使う。
//最小二乗法によって、LAC_dataの線減弱係数と、LUT_dataの微分断面積に比例係数Nを掛けた値の誤差resが最小となるNを求める。
//横軸エネルギー(線形)、縦軸微分断面積(対数)のグラフ上に、LAC_dataの線減弱係数を乗せると考えればよい。
//最小化については重みなしの最小二乗法を利用した。
//求める条件は f=Σ{Y-(Z*X))^2 について f'=0となる　N=Zである。
//誤差resは最小化されたfであり、0以上の最小の値が最も良い値となる。
{
	{
		if (LAC_data.Row() != LUT_data.Row()) {
			printf("Size error! (est_Z_least_square_method(row=%d,row=%d)\n", LAC_data.Row(), LUT_data.Row());
			system("pause");
			exit(EXIT_FAILURE);
		}

		MATRIX<double> LAC_LUT;
		MATRIX<double> LUT_LUT;
		MATRIX<double> RES;
		double sum_LAC_LUT, sum_LUT_LUT;

		LAC_LUT = LAC_data * LUT_data;
		LUT_LUT = LUT_data * LUT_data;
		sum_LAC_LUT = LAC_LUT.sum();
		sum_LUT_LUT = LUT_LUT.sum();

		N = sum_LAC_LUT / sum_LUT_LUT;

		RES = (LAC_data - N * LUT_data)*(LAC_data - N * LUT_data);
		res = RES.sum();
	}
}

void est_Z_least_square_method_LN(MATRIX<double> LAC_data, MATRIX<double> LUT_data, double &N, double &res)
//引数のNおよびresは値を受け取るための数値で、関数の戻り値の代わりに使う。
//最小二乗法によって、LAC_dataの線減弱係数と、LUT_dataの微分断面積に比例係数Nを掛けた値の誤差resが最小となるNを求める。
//横軸エネルギー(線形)、縦軸微分断面積(対数)のグラフ上に、LAC_dataの線減弱係数を乗せると考えればよい。
//最小化については各dataの対数についての最小二乗法を利用した。
//求める条件は f=Σ{Y-(Z*X))^2 について f'=0となる　N=Zである。
//誤差resは最小化されたfであり、0以上の最小の値が最も良い値となる。
{
	{
		if (LAC_data.Row() != LUT_data.Row()) {
			printf("Size error! (est_Z_least_square_method(row=%d,row=%d)\n", LAC_data.Row(), LUT_data.Row());
			system("pause");
			exit(EXIT_FAILURE);
		}

		MATRIX<double> LN_LAC_LUT;
		MATRIX<double> RES;
		double ave_LN_LAC_LUT;

		LN_LAC_LUT = log(LAC_data / LUT_data);
		ave_LN_LAC_LUT = LN_LAC_LUT.average();

		N = exp(ave_LN_LAC_LUT);

		RES = (LAC_data - N * LUT_data)*(LAC_data - N * LUT_data);
		res = RES.sum();
	}
}

void est_Z_minimize_max_relative_error(MATRIX<double> LAC_data, MATRIX<double> LUT_data, double &N, double &res)
//引数のNおよびresは値を受け取るための数値で、関数の戻り値の代わりに使う。
//最小二乗法によって、LAC_dataの線減弱係数と、LUT_dataの微分断面積に比例係数Nを掛けた値の誤差resが最小となるNを求める。
//横軸エネルギー(線形)、縦軸微分断面積(対数)のグラフ上に、LAC_dataの線減弱係数を乗せると考えればよい。
//最小化については重みなしの最小二乗法を利用した。
//求める条件は f=Σ{Y-(Z*X))^2 について f'=0となる　N=Zである。
//誤差resは最小化されたfであり、0以上の最小の値が最も良い値となる。
{
	{
		if (LAC_data.Row() != LUT_data.Row()) {
			printf("Size error! (est_Z_least_square_method(row=%d,row=%d)\n", LAC_data.Row(), LUT_data.Row());
			system("pause");
			exit(EXIT_FAILURE);
		}

		MATRIX<double> LUT_LAC;
		MATRIX<double> RES;


		LUT_LAC = LUT_data / LAC_data;

		N = 2.0 / (LUT_LAC.max() + LUT_LAC.min());

		RES = abs((LAC_data - N * LUT_data) / LAC_data);
		res = RES.max();
	}
}

MATRIX<double> LUT_resize_Z1_30_281_2902(MATRIX<double> Zeff_table, MATRIX<double> energy_list) {
	//Zeff_LUT_Z1_30(281,2901+1).raw 専用
	int data_num = energy_list.Row();
	MATRIX<double> LUT(data_num, Zeff_table.Col());

	int i;
	int E1_num, E2_num, E2_cnt;

	E2_cnt = 0;
	for (i = 0;i<data_num;i++) {
		E1_num = int(energy_list.out(i + 1, 1) * 10 + 0.5);//エネルギーの比較に失敗しないよう整数に丸める。切捨ての結果、5を単位とした整数になる。
		do {
			E2_cnt++;
			E2_num = int(Zeff_table.out(E2_cnt, 1) * 10 + 0.5);
		} while (E1_num != E2_num);
		LUT.putrow(i + 1, Zeff_table.pullrow(E2_cnt));
	}

	if (energy_list.out(energy_list.Row(), 1) != LUT.out(LUT.Row(), 1)) {
		printf("LUT resize error!(Z1_30_281_2902) (%lf,%lf)\n", energy_list.out(energy_list.Row(), 1), LUT.out(LUT.Row(), 1));
		system("pause");
		exit(EXIT_FAILURE);
	}

	return LUT;
}

void data_output(char *filename, double Zeff, double Neff, double res, MATRIX<double> input_data, MATRIX<double> LUT) {
	MATRIX<double> data(input_data.Row(), 4);
	data.put(1, 1, input_data);
	data.putcol(3, LUT.pullcol(int(Zeff*100.0 + 0.5) - 100 + 1));

	MATRIX<double> relative_error;
	relative_error = (input_data.pullcol(2) - Neff * LUT.pullcol(int(Zeff*100.0 + 0.5) - 100 + 1)) / input_data.pullcol(2);
	data.putcol(4, relative_error * 100);

	FILE *fp;
	fp = fopen(filename, "w");
	fprintf(fp, "Zeff=\t%.2lf\n", Zeff);
	fprintf(fp, "Neff(10^23 cm^3)=\t%.4e\n", Neff*10.0);
	fprintf(fp, "Ne(10^23 cm^3)=\t%.4e\n", Zeff*Neff*10.0);
	fprintf(fp, "res=\t%.4e\n", res);
	fprintf(fp, "relative error(%%)=\t%.4e\n", abs(relative_error).max() * 100);
	fprintf(fp, "energy(keV)\tinput(1/cm)\tZeff(barn)\trelative error(%%)\n");
	data.write(fp, "%.4e");
	fclose(fp);
}