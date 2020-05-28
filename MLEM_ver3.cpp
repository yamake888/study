#include"matrix_tool.h"
#include<iostream>
#include<math.h>
#include<vector>
#include<windows.h>
#include<direct.h>
using namespace std;

const int pix_num = 1024;
const char input_row_filename[256] = "B_recon.raw";
const char input_col_filename[256] = "Si_recon.raw";
int row_atom_number = 5;
int col_atom_number = 14;
const char output_filename[256] = "result.raw";

double calc(int AZ,int AC,int BZ,int BC) {
    int num = 2;
    vector<pair<int, int>> atom(num);
    atom[0].first = AZ;
    atom[1].first = BZ;
    atom[0].second = AC;
    atom[1].second = BC;
    int sum = 0;
    float tmp = 0.0;
    for (int i = 0; i < num; i++) {
        sum += atom[i].first * atom[i].second;
    }
    for (int i = 0; i < num; i++) {
        tmp += ((float)atom[i].first * (float)atom[i].second / (float)sum) * powf((float)atom[i].first, 2.94);
    }
    float ans;

    ans = powf(tmp, 0.34013);//y”õlz1/2.94=0.34013
    return (double)ans;
}
int main() {
    MATRIX<double> Row_image, Col_image;
	_chdir("input");
    Row_image.read(input_row_filename,pix_num, pix_num);
    Col_image.read(input_col_filename,pix_num, pix_num);
	Row_image = Row_image / (double)row_atom_number;
	Col_image = Col_image / (double)col_atom_number;
    MATRIX<double> output(pix_num, pix_num);
    for (int i = 1; i <= pix_num; ++i) {
        for (int j = 1; j <= pix_num; ++j) {
            output.in(i, j, calc(row_atom_number, (int)(Row_image.out(i, j)*10000), col_atom_number, (int)(Col_image.out(i, j)*10000)));
        }
    }
	_chdir("../output");
    output.write(output_filename, "");
    return 0;
}