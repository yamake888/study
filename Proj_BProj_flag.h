#pragma once

// 最終更新 19,JUN,2018
//投影方向数の多い再構成にも対応

#define _CRT_SECURE_NO_WARNINGS //コンパイル時のfopenなどの警告を無視するため
#define _USE_MATH_DEFINES//円周率M_PIを使うため
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "matrix_tool.h"

/*
画像再構成のための関数群
①CT2.0の焼き直し
②コンボリューション再構成法(FBPまたはCT2.0)に関係したサブルーチン群
③ML-EM画像再構成法に関係したクラスとサブルーチン群
番外 少し役立つ関数群 (サイノグラムを角度方向にずらすとか。)
matrix_tool.hのインクルードが必要
*/

/* CT2.0の焼き直し プロトタイプ宣言 */

template<typename mTYPE1, typename mTYPE2>
void BackProjection_CT2p0(const MATRIX<mTYPE2>&, MATRIX<mTYPE1>&);//CT2.0の逆投影をMARIXで書き換えたもの。アクセス違反で動作しない。

template<typename mTYPE1, typename mTYPE2>
void BackProjection_CT2p0_mod(const MATRIX<mTYPE2>&, MATRIX<mTYPE1>&);//アクセス違反を避けたもの

template<typename mTYPE>
void ConvolutionIntegral_CT2p0(MATRIX<mTYPE>&, double);//サイノグラムにコンボリューションをかける関数をMATRIXで書き換えたもの

double SLConFun_CT2p0(int, double);//コンボリューションで使われる。


/* コンボリューション再構成法(FBP)に関係したサブルーチン群 プロトタイプ宣言 */

template<typename mTYPE1, typename mTYPE2>
void BackProj_simpl(const MATRIX<mTYPE1>&, MATRIX<mTYPE2>&, int);//検出確率計算が単純なもので単純逆投影を行う

template<typename mTYPE>
MATRIX<mTYPE> ConvInte(double(*conv_func)(int, double), const MATRIX<mTYPE>&, double);//サイノグラムに対して指定したフィルタで重畳積分をかける

double Shepp_Logan_cf(int, double);//指定されるフィルタのひとつ。Shepp Logan convolution filter


/* ML-EM画像再構成法に関係したクラスとサブルーチン群 プロトタイプ宣言 */

//台形で計算した検出確率を管理するクラス。使うときには始めに続く処理のために、まずはこのクラスの変数を宣言すること。
class DP_trapezoid
{
public:
  MATRIX<short> dp_Xs;//投影位置の中心がかかるステップ番号を[0]から[ステップ番-1]までの数字で管理させる。
  MATRIX<double> dp_L;
  MATRIX<double> dp_C;
  MATRIX<double> dp_R;

  int Num_angle()const { return dp_Num_angle; }
  int Step()const { return dp_Step; }
  int Height()const { return dp_Height; }
  int Width()const { return dp_Width; }

  void Clear()
  {
    dp_Num_angle = 1;
    dp_Step = 1;
    dp_Height = 1;
    dp_Width = 1;
    dp_Xs.resize(dp_Num_angle, dp_Height*dp_Width);
    dp_L.resize(dp_Num_angle, dp_Height*dp_Width);
    dp_C.resize(dp_Num_angle, dp_Height*dp_Width);
    dp_R.resize(dp_Num_angle, dp_Height*dp_Width);
  };

  DP_trapezoid(int Num_angle, int Step, int Height, int Width)
  {
    dp_Num_angle = Num_angle;
    dp_Step = Step;
    dp_Height = Height;
    dp_Width = Width;
    dp_Xs.resize(Num_angle, Height*Width);
    dp_L.resize(Num_angle, Height*Width);
    dp_C.resize(Num_angle, Height*Width);
    dp_R.resize(Num_angle, Height*Width);
  }

private:
  int dp_Num_angle;
  int dp_Step;
  int dp_Height;
  int dp_Width;
};

//台形で計算した検出確率を管理するクラス。検出確率を保持するためのメモリが大きすぎる場合に使う。
class DP_trapezoid_fragment
{
public:
  MATRIX<double> dp;

  int Num_angle()const { return dp_Num_angle; }
  int Step()const { return dp_Step; }
  int Height()const { return dp_Height; }
  int Width()const { return dp_Width; }

  DP_trapezoid_fragment(MATRIX<double> angle, int Step, int Height, int Width)//angleはdegreeで与えること。
  {
    if (angle.Col() != 1) {
      printf("angle size error! (DP_trapezoid_fragment:angle.Col()=%d", angle.Col());
      system("pause");
    }

    dp_angle = angle;
    dp_Num_angle = dp_angle.Row();
    dp_Step = Step;
    dp_Height = Height;
    dp_Width = Width;
    dp.resize(dp_Num_angle, 5);

    printf("*** calcDP_trapezoid *** :\r");

    int i;
    double th, si, co;
    double a, b, c;
    double buf;

    for (i = 1; i <= dp_Num_angle; i++)
    {
      th = M_PI*(angle.out(i, 1) / 180.0);
      si = sin(th);
      co = cos(th);

      a = fabs(si);
      b = fabs(co);
      if (a < b) {
        buf = a; a = b; b = buf;
      }
      c = 1.0 / a;

      dp.in(i, 1, a);
      dp.in(i, 2, b);
      dp.in(i, 3, c);
      dp.in(i, 4, si);
      dp.in(i, 5, co);
    }
    printf("*** calcDP_trapezoid *** : end \n");

  }

private:
  MATRIX<double> dp_angle;
  int dp_Num_angle;
  int dp_Step;
  int dp_Height;
  int dp_Width;
};

void calcDP_trapezoid(DP_trapezoid&, int scan = 1);//検出確率を計算する。

template<typename mTYPE1, typename mTYPE2>
void ML_EM(const DP_trapezoid&, const MATRIX<mTYPE2>&, MATRIX<mTYPE1>&, double, int);//投影、逆投影をセットにしてML-EM画像再構成を行う

template<typename mTYPE1, typename mTYPE2>
void ML_EM_fragment(const DP_trapezoid_fragment&, const MATRIX<mTYPE2>&, MATRIX<mTYPE1>&, double, int, double);//投影、逆投影をセットにしてML-EM画像再構成を行う

template<typename mTYPE>
void transmission_distance(const DP_trapezoid&, const MATRIX<short>&, int, double, MATRIX<mTYPE>&);//用意した画像の各要素について透過距離のサイノグラムを作る。無視するピクセルには0、計算したいピクセルには要素ごとに1から順に整数の通し番号を付けること。下手につけるとアクセスエラーが起こる。

template<typename mTYPE1, typename mTYPE2>
void Proj_trape(const DP_trapezoid&, const MATRIX<mTYPE1>&, double, MATRIX<mTYPE2>&);//画像を投影する。

template<typename mTYPE2, typename mTYPE1>
void BackProj_trape(const DP_trapezoid&, const MATRIX<mTYPE2>&, MATRIX<mTYPE1>&);//サイノグラムを逆投影する。

template<typename mTYPE1, typename mTYPE2>
void Proj_trape_frag(const DP_trapezoid_fragment&, const MATRIX<mTYPE1>&, double, MATRIX<mTYPE2>&, double);//画像を投影する。

template<typename mTYPE2, typename mTYPE1>
void BackProj_trape_frag(const DP_trapezoid_fragment&, const MATRIX<mTYPE2>&, MATRIX<mTYPE1>&);//サイノグラムを逆投影する。


/* 少し役立つ関数群 */
template<typename mTYPE>
MATRIX<mTYPE> click_angle(const MATRIX<mTYPE>&, int, int scan = 1);//サイノグラムを入力した[行番号]から始まるようにずらす。行方向に角度を取っていること。scan=1でハーフ、2でフル
template<typename mTYPE>
MATRIX<mTYPE> devide_each(const MATRIX<mTYPE>&);//投影方向ごとの最大値で各投影を割った行列を返す。



/* ここからは関数の定義 */

/* CT2.0の焼き直し */
template<typename mTYPE1, typename mTYPE2>
void BackProjection_CT2p0(const MATRIX<mTYPE2> &sinogram, MATRIX<mTYPE1> &image)
//  CT2.0のvoid bp()のみをMATRIX_toolを使って書き写したもの。
//  CT2.0ではsinogramは、対数が取られ(void readf(char string[80])にて)、コンボリューションもされている(void con()にて)
//  このサブ関数を使う場合、これらの処理は別に行うこと。
//  ※1でアクセス違反が起こる。  
{
  double PI = atan(1.0) * 4;//CT2.0では円周率にM_PIを使わず、このように計算している。
  int k, l;//画像左上端を原点としたときの画素の位置を表す変数
  int n;//投影方向の番号を表す変数
  double x, y;//画像中央を原点としたときの画素の位置を表す変数
  double s;//ステップの中心を原点として画素の位置を投影ステップに投影した位置。
  double g;//sの小数部分
  double m01;//sの整数部分
  int m0;//m01をint型に変更して扱うための変数
  double ff;//逆投影の加算計算用の変数

  int N = sinogram.Col();//Nはステップ数。サイノグラムの列数に等しい。
  int M = sinogram.Row();//Mは投影方向数。サイノグラムの行数に等しい。
  image.resize(N, N);//再構成画像を格納する行列。CT2.0ではf[N][N]に相当する。

  for (k = 0; k < N; k++)
  {
    y = k - N / 2.0;
    for (l = 0; l < N; l++)
    {
      x = l - N / 2.0;
      image.in(k + 1, l + 1, 0.0);
      ff = 0.0;
      for (n = 0; n < M; n++)
      {
        s = x*cos((n*PI) / (M - 0.0)) + y*sin((n*PI) / (M - 0.0));
        m01 = floor(s);
        g = s - m01;
        m0 = (int)m01;
        m0 = int(m0 + N / 2.0);//CT2.0にはキャストは書かれていないが、警告が邪魔なので書き足す。暗黙のキャストと同じなので影響はない。
        ff = ff + (1 - g)*sinogram.out(n + 1, m0 + 1) + g*sinogram.out(n + 1, m0 + 2);//sinogramの行列は左上端の要素が(1,1)なので、+1ずらした。※1
      }
      image.in(k + 1, l + 1, ff*(1 / (2 * PI))*(PI / (M - 0.0)));//単純に考えると2*PIではなくPIにしておけば、投影方向数での平均化と思えるのだが。
      //ここでConvolutionIntegral_CT2p0でsinogramが2倍に大きく計算されていることは関係しないだろうか。
    }
  }

}

template<typename mTYPE1, typename mTYPE2>
void BackProjection_CT2p0_mod(const MATRIX<mTYPE2> &sinogram, MATRIX<mTYPE1> &image)
//  CT2.0のvoid bp()のみをMATRIX_toolを使って書き写したもの。
//  CT2.0ではsinogramは、対数が取られ(void readf(char string[80])にて)、コンボリューションもされている(void con()にて)
//  このサブ関数を使う場合、これらの処理は別に行うこと。
//  ※1でアクセス違反が起こる。  
{
  double PI = atan(1.0) * 4;//CT2.0では円周率にM_PIを使わず、このように計算している。
  int k, l;//画像左上端を原点としたときの画素の位置を表す変数
  int n;//投影方向の番号を表す変数
  double x, y;//画像中央を原点としたときの画素の位置を表す変数
  double s;//ステップの中心を原点として画素の位置を投影ステップに投影した位置。
  double g;//sの小数部分
  double m01;//sの整数部分
  int m0;//m01をint型に変更して扱うための変数
  double ff;//逆投影の加算計算用の変数

  int N = sinogram.Col();//Nはステップ数。サイノグラムの列数に等しい。
  int M = sinogram.Row();//Mは投影方向数。サイノグラムの行数に等しい。
  image.resize(N, N);//再構成画像を格納する行列。CT2.0ではf[N][N]に相当する。

  for (k = 0; k < N; k++)
  {
    y = k - N / 2.0;
    for (l = 0; l < N; l++)
    {
      x = l - N / 2.0;
      image.in(k + 1, l + 1, 0.0);
      ff = 0.0;
      for (n = 0; n < M; n++)
      {
        s = (x + 0.5)*cos((n*PI) / (M - 0.0)) + (-y - 0.5)*sin((n*PI) / (M - 0.0)) - 0.5;
        m01 = floor(s);
        g = s - m01;
        m0 = (int)m01;
        m0 = int(m0 + N / 2.0);//CT2.0にはキャストは書かれていないが、警告が邪魔なので書き足す。暗黙のキャストと同じなので影響はない。
        if (m0 >= 0 && m0 < N - 1)
        {
          ff = ff + (1 - g)*sinogram.out(n + 1, m0 + 1) + g*sinogram.out(n + 1, m0 + 2);//sinogramの行列は左上端の要素が(1,1)なので、+1ずらした。※1
        }
      }
      image.in(k + 1, l + 1, ff*(1 / (PI))*(PI / (M - 0.0)));//単純に考えると2*PIではなくPIにしておけば、投影方向数での平均化。
      //ここでConvolutionIntegral_CT2p0でsinogramが2倍に大きく計算されていることは関係しないだろうか。
    }
  }

}

template<typename mTYPE>
void ConvolutionIntegral_CT2p0(MATRIX<mTYPE> &sinogram, double sampling_interval)
//  CT2.0のcon()をMATRIX_toolを使って書き写したもの
//  sinogramが必要なので渡している。
//  sampling_intervalはグローバル変数で定義されていたdrのことで、1投影の幅(cm)をステップ数で割った値を表す。
{
  int n, md;//それぞれ投影方向、各投影でのステップを表す。
  int m;//1つの投影で端から重畳積分を行うために使う。
  double q_temp = 0.0;//重畳積分の値を加算するための変数

  int N = sinogram.Col();//Nはステップ数。サイノグラムの列数に等しい。
  int M = sinogram.Row();//Mは投影方向数。サイノグラムの行数に等しい。
  MATRIX copy = MATRIX(sinogram);//copyを元にして重畳積分されたsinogramを計算する。
  sinogram.constant(0.0);//重畳積分の結果を格納するのですべての数値を0にしておく。この処理はCT_2.0では見あたらなかった。
  double dr = sampling_interval;

  for (n = 0; n < M; n++)
  {
    for (md = 0; md < N; md++)
    {
      q_temp = 0;
      for (m = 0; m < N; m++)
      {
        q_temp = q_temp + dr*(sinogram.out(n + 1, md + 1) + SLConFun_CT2p0(md - m, dr)*copy.out(n + 1, m + 1));
        //SLConFun_CT2p0のコメントにも書いたが、分子が2倍大きいのでsinogramの値も2倍になっている。
      }
      sinogram.in(n + 1, md + 1, mTYPE(q_temp));
    }
  }
}

double SLConFun_CT2p0(int m, double sampling_interval)
//  CT2.0のh(int m)を書き写したもの。con()で使われている。
//  sampling_intervalはグローバル変数で定義されていたdrのことで、1投影の幅(cm)をステップ数で割った値を表す。
{
  double PI = atan(1.0) * 4;//CT2.0では円周率にM_PIを使わず、このように計算している。
  double dr = sampling_interval;

  double hh = 0.0;
  hh = 4 / (PI*dr*dr*(1 - 4 * m*m));
  //ttp://fiber.shinshu-u.ac.jp/koseki/publication/2003/bio2003-3.pdf等では分子は2なのだが。

  return hh;
}

/* コンボリューション再構成法(FBP)に関係したサブ関数群 */
template<typename mTYPE1, typename mTYPE2>
void BackProj_simpl(const MATRIX<mTYPE1> &sinogram, MATRIX<mTYPE2> &image, int scan)
/*
scan=1でハーフスキャン再構成。scan=2でフルスキャン再構成。
ピクセルの投影幅を1ステップ分と見なして検出確率を計算する単純逆投影プログラム。simpleと呼ぶことにする。
sinogramは行が投影方向、列がステップ数となるように準備すること。線減弱係数の投影の和(-log(I/I0))にしておくこと。
(仕様)
　再構成画像ピクセルサイズ(縦の長さと横の長さ)とステップの幅がすべて等しい。
 　再構成画像の縦、横ピクセル数、ステップ数がすべて等しい。ただしすべて偶数でも奇数でも対応している。
  　先頭のステップ、最後尾のステップでの検出も無視しないように作ってある。
   */
{
  if (scan != 1 && scan != 2)
  {
    printf("scan value error. BackProjection_A\n");
    system("pause");
    exit(EXIT_FAILURE);
  }

  const int Num_angle = sinogram.Row();
  const int Step = sinogram.Col();
  const int Height = Step;
  const int Width = Step;
  const int Num_pixel = Height*Width;

  const double th_U = scan*M_PI / double(Num_angle);
  const double y0 = 0.5*Height - 0.5;
  const double x0 = -0.5*Width + 0.5;

  if (Height != image.Row() || Width != image.Col()) image.resize(Height, Width);

  int i, j, k, I;
  int ixx, IXX;
  double xx, XX, x, y;
  double th, si, co;
  double dp1, dp2;

  mTYPE1 *sinogram_value;
  sinogram_value = sinogram.pVALUE();

  mTYPE2 *img;
  img = image.pVALUE();

  for (i = 0; i < Num_pixel; i++) {
    img[i] = 0.0;
  }

  for (k = 0; k < Num_angle; k++)
  {
    IXX = k*Step;
    th = th_U*k;
    si = sin(th);
    co = cos(th);
    for (i = 0; i < Height; i++)
    {
      I = i*Width;
      y = y0 - i;
      XX = y*si + Step / 2.0 - 0.5;
      for (j = 0; j < Width; j++)
      {
        x = x0 + j;
        xx = x*co + XX;

        ixx = int(floor(xx));

        if (ixx < -1 || Step - 1 < ixx) continue;//検出効率がdp1、dp2ともに0になるのでcontinue    

        dp2 = xx - double(ixx);

        if (-1 < ixx && xx < Step - 1)
        {
          dp1 = 1.0 - dp2;
          img[I + j] += mTYPE2(dp2*sinogram_value[IXX + ixx + 1] + dp1*sinogram_value[IXX + ixx]);
          continue;
        }
        else if (ixx == -1)
        {
          img[I + j] += mTYPE2(dp2*sinogram_value[IXX + ixx + 1]);
          continue;
        }
        else if (ixx == Step - 1)
        {
          dp1 = 1.0 - dp2;
          img[I + j] += mTYPE2(dp1*sinogram_value[IXX + ixx]);
          continue;
        }
      }
    }
  }

  for (i = 0; i < Num_pixel; i++) {
    img[i] /= Num_angle;
  }
}

template<typename mTYPE>
MATRIX<mTYPE> ConvInte(double(*conv_func)(int, double), const MATRIX<mTYPE> &sinogram, double step_interval)
{
  /*
  コンボリューションを行う関数。
  引数にコンボリューション関数(サブ関数の名前の部分)を与えると、これを使って計算する。
  sinogramは行が投影方向、列がステップ数となるように準備すること。線減弱係数の投影の和(-log(I/I0))にしておくこと。
  step_intervalは投影ステップのステップ間距離(cm)。ここの単位がcmなので、最終的に再構成画像の単位が(1/cm)になる。
  */

  int i, j, k, K;
  int Num_angle = sinogram.Row();
  int Step = sinogram.Col();
  int Num_data = Num_angle*Step;

  MATRIX<mTYPE> conv_sinogram(Num_angle, Step);
  mTYPE *sino = sinogram.pVALUE();
  mTYPE *conv_sino = conv_sinogram.pVALUE();

  for (i = 0; i < Num_data; i++) {
    conv_sino[i] = 0;
  }

  for (k = 0; k < Num_angle; k++)
  {
    K = k*Step;
    for (i = 0; i < Step; i++)
    {
      for (j = 0; j < Step; j++)
      {
        conv_sino[K + i] += mTYPE(conv_func(i - j, step_interval)*sino[K + j]);
      }
    }
  }

  for (i = 0; i < Num_data; i++) {
    conv_sino[i] *= mTYPE(step_interval);
  }

  return conv_sinogram;
}

double Shepp_Logan_cf(int pixel_distance, double step_interval)
{
  /*
  ConvInteに使うコンボリューション関数。
  pixel_distanceはピクセル数で数えたときの、注目ステップと影響を考えるステップとの距離。
  step_intervalは投影ステップのステップ間距離(cm)。ここの単位がcmなので、最終的に再構成画像の単位が(1/cm)になる
  shepp-loganフィルタについて調べると係数の付き方が様々で混乱する。これはフーリエ変換の係数をどこにつけるかの問題で、
  画像再構成時の最後の投影方向数で割る計算と都合を合わせないといけない。
  BackProj_Aの投影方向数で割る計算で方向数をそのまま(πを掛けずに)計算する
  bp_img=bp_img/M;
  こととし、下記の通りに計算する。2016/11/11hamaguchi
  すぐ見つかるのはh=2.0/(M_PI*M_PI*step_interval*step_interval*(1.0-4.0*(pixel_distance*pixel_distance)));だが、落ち着くこと。
  */
  double h;

  h = 2.0 / (M_PI*step_interval*step_interval*(1.0 - 4.0*(pixel_distance*pixel_distance)));

  return h;
}

// Ram_Lak_cfは未完成
double Ram_Lak_cf(int pixel_distance, double step_interval)
//ConvInteに使うコンボリューション関数。
//pixel_distanceはピクセル数で数えたときの、注目ステップと影響を考えるステップとの距離。
//step_intervalは投影ステップのステップ間距離(cm)。ここの単位がcmなので、最終的に再構成画像の単位が(1/cm)になる
//Ram_Lakフィルタについても、調べると係数の付き方が様々で混乱するかもしれない。
//Shepp_Logan_cfと同様に逆投影のプログラムとの都合を考え、下記の通りに計算する。2016/11/11hamaguchi

{
  double h;

  if (pixel_distance % 2 == 1)//奇数
  {
    h = -1.0 / (pixel_distance*pixel_distance*M_PI*step_interval*step_interval);
    return h;
  }
  else if ((pixel_distance / 2) != 0)//ゼロ以外の偶数
  {
    h = 0.0;
    return h;
  }
  else
  {
    h = M_PI / (4.0*step_interval*step_interval);//ゼロ
    return h;
  }
}


/* ML-EM画像再構成法に関係したクラスとサブ関数群 */

void calcDP_trapezoid(DP_trapezoid& DP, int scan)
{
  printf("*** calcDP_trapezoid *** :\r");

  const int Num_angle = DP.Num_angle();
  const int Step = DP.Step();
  const int Height = DP.Height();
  const int Width = DP.Width();
  const int Num_pixel = Height*Width;
  const int DP_data = Num_angle*Height*Width;

  int i, j, k, J, K;
  int ixx;
  double th, si, co;
  double a, b, c, ixx_side, d;
  double dp_LCR[3];
  double x, y, xx, XX;
  double buf;

  const double th_U = scan*M_PI / double(Num_angle);
  const double y0 = 0.5*Height - 0.5;
  const double x0 = -0.5*Width + 0.5;
  const double step0 = Step / 2.0 - 0.5;

  short *DP_ixx;
  DP_ixx = DP.dp_Xs.pVALUE();

  double *DP_LCR[3];//DPをdoubleにするならここを対応させること。
  DP_LCR[0] = DP.dp_L.pVALUE();
  DP_LCR[1] = DP.dp_C.pVALUE();
  DP_LCR[2] = DP.dp_R.pVALUE();

  for (i = 0; i < DP_data; i++) {
    DP_ixx[i] = 0;
    DP_LCR[0][i] = 0.0;
    DP_LCR[1][i] = 0.0;
    DP_LCR[2][i] = 0.0;
  }

  th = 0;
  for (k = 0; k < Num_angle; k++)
  {
    K = k*Num_pixel;

    th = th_U*k;
    si = sin(th);
    co = cos(th);

    a = fabs(si);
    b = fabs(co);
    if (a < b) {
      buf = a; a = b; b = buf;
    }
    c = 1.0 / a;

    for (i = 0; i < Height; i++)
    {
      J = K + i*Width;
      y = y0 - i;
      XX = y*si + step0;
      for (j = 0; j < Width; j++)
      {
        x = x0 + j;
        xx = x*co + XX;//1つ目のステップの中心をゼロする座標で見たピクセル中心の投影位置
        ixx = int(floor(xx + 0.5));//ピクセル中心の投影位置が、いくつ目のステップの中にあるかを計算。

        if (ixx<-1 || ixx>Step) {
          DP_ixx[J + j] = -2;
          continue;
        }
        dp_LCR[0] = dp_LCR[1] = dp_LCR[2] = 0.0;

        //ひとつ前のステップに入る検出確率(dp_LCR[0])
        ixx_side = ixx - 0.5;//ひとつ前のステップとの境界の座標
        d = ixx_side - (xx - 0.5*(a - b));//境界と台形の肩との距離
        if (d >= 0) {//境界のほうが右
          dp_LCR[0] = (0.5*b + d)*c;
        }
        else if ((-d) < b) {//境界のほうが右だが台形の足の位置よりは右
          dp_LCR[0] = 0.5*c*(b + d) / b*(b + d);
        }

        //ひとつ後のステップに入る検出確率(dp_LCR[1])
        ixx_side = ixx + 0.5;//ひとつ後のステップとの境界の座標
        d = (xx + 0.5*(a - b)) - ixx_side;//境界と台形の肩との距離
        if (d >= 0) {//境界のほうが左
          dp_LCR[2] = (0.5*b + d)*c;
        }
        else if ((-d) < b) {
          dp_LCR[2] = 0.5*c*(b + d) / b*(b + d);
        }

        //投影位置の検出器に入る検出確率(dp_LCR[2])
        dp_LCR[1] = 1 - dp_LCR[0] - dp_LCR[2];

        //格納
        DP_ixx[J + j] = ixx;
        DP_LCR[0][J + j] = double(dp_LCR[0]);//DPをdoubleにするならここを対応させること。
        DP_LCR[1][J + j] = double(dp_LCR[1]);//DPをdoubleにするならここを対応させること。
        DP_LCR[2][J + j] = double(dp_LCR[2]);//DPをdoubleにするならここを対応させること。
      }
    }
  }
  printf("*** calcDP_trapezoid *** :   end\n");
}

template<typename mTYPE1, typename mTYPE2>
void ML_EM(const DP_trapezoid &DP, const MATRIX<mTYPE2> &sinogram, MATRIX<mTYPE1> &image, double pix_size, int iteration)
//  sinogramは行が投影方向、列がステップ数となるように準備すること。線減弱係数の投影の和(-log(I/I0))にしておくこと。
//  pixel_size(cm)。ここの単位がcmなので、最終的に再構成画像の単位が(1/cm)になる。
{
  if (DP.Height() != image.Row() || DP.Width() != image.Col()) {
    image.resize(DP.Height(), DP.Width());
  }
  image.constant(1);

  MATRIX<mTYPE2> sinogram_calc;
  MATRIX<mTYPE2> sinogram_ratio;
  MATRIX<mTYPE1> image_ratio;

  int i;
  for (i = 0; i < iteration; i++)
  {
    printf("*** ML-EM *** : %d/%d\r", i + 1, iteration);

    Proj_trape(DP, image, pix_size, sinogram_calc);

    sinogram_ratio = m_cast<mTYPE2>(sinogram / sinogram_calc);

    BackProj_trape(DP, sinogram_ratio, image_ratio);

    image = m_cast<mTYPE1>(image*image_ratio);
  }
  printf("*** ML-EM *** : end  \n");
}

template<typename mTYPE1, typename mTYPE2>
void ML_EM_fragment(const DP_trapezoid_fragment &DP, const MATRIX<mTYPE2> &sinogram, MATRIX<mTYPE1> &image, double pix_size, int iteration, double threshold=-1.0)
//  sinogramは行が投影方向、列がステップ数となるように準備すること。線減弱係数の投影の和(-log(I/I0))にしておくこと。
//  pixel_size(cm)。ここの単位がcmなので、最終的に再構成画像の単位が(1/cm)になる。
{
  if (DP.Height() != image.Row() || DP.Width() != image.Col()) {
    image.resize(DP.Height(), DP.Width());
  }
  image.constant(1);

  MATRIX<mTYPE2> sinogram_calc;
  MATRIX<mTYPE2> sinogram_ratio;
  MATRIX<mTYPE1> image_ratio;

  int i;
  for (i = 0; i < iteration; i++)
  {
    printf("*** ML-EM *** : %d/%d\r", i + 1, iteration);

    Proj_trape_frag(DP, image, pix_size, sinogram_calc, threshold);

    sinogram_ratio = m_cast<mTYPE2>(sinogram / sinogram_calc);

    BackProj_trape_frag(DP, sinogram_ratio, image_ratio);

    image = m_cast<mTYPE1>(image*image_ratio);
  }
  printf("*** ML-EM *** : end  \n");
}


template<typename mTYPE>
void transmission_distance(const DP_trapezoid& DP, const MATRIX<short>& image, int Num_component, double pix_size, MATRIX<mTYPE>& TD)
{
  printf("*** transmission_distance *** :\r");

  const int Num_angle = DP.Num_angle();
  const int Step = DP.Step();
  const int Height = DP.Height();
  const int Width = DP.Width();
  const int Num_data = Num_angle*Step;
  const int Num_pixel = Height*Width;

  TD.resize(Num_component, Num_angle*Step);

  int i, j, k, I, J, K;
  int Component;
  int ixx, IXX;

  short *DP_ixx;
  DP_ixx = DP.dp_Xs.pVALUE();

  double *DP_LCR[3];//DPをdoubleにするならここを対応させること。
  DP_LCR[0] = DP.dp_L.pVALUE();
  DP_LCR[1] = DP.dp_C.pVALUE();
  DP_LCR[2] = DP.dp_R.pVALUE();

  short *img = image.pVALUE();

  mTYPE **TD_value;
  TD_value = (mTYPE**)malloc(sizeof(mTYPE*)*Num_component);
  for (i = 0; i < Num_component; i++) {
    TD_value[i] = TD.pVALUE(i + 1, 1);
    for (j = 0; j < Num_data; j++) {
      TD_value[i][j] = 0.0;
    }
  }

  for (k = 0; k < Num_angle; k++) {
    K = k*Num_pixel;
    IXX = k*Step;
    for (i = 0; i < Height; i++) {
      I = i*Width;
      J = K + i*Width;
      for (j = 0; j < Width; j++) {
        ixx = DP_ixx[J + j];

        if (ixx<-1 || ixx>Step)continue;
        if (img[I + j] - 1 < 0) continue;
        Component = img[I + j] - 1;

        if (0 < ixx && ixx < Step - 1) {
          TD_value[Component][IXX + ixx - 1] += mTYPE(DP_LCR[0][J + j]);
          TD_value[Component][IXX + ixx] += mTYPE(DP_LCR[1][J + j]);
          TD_value[Component][IXX + ixx + 1] += mTYPE(DP_LCR[2][J + j]);
          continue;
        }

        switch (ixx) {
        case -1:
          TD_value[Component][IXX] += mTYPE(DP_LCR[2][J + j]);
          continue;
        case 0:
          TD_value[Component][IXX] += mTYPE(DP_LCR[1][J + j]);
          TD_value[Component][IXX + 1] += mTYPE(DP_LCR[2][J + j]);
          continue;
        }

        switch (ixx - Step) {
        case -1:
          TD_value[Component][IXX + Step - 2] += mTYPE(DP_LCR[0][J + j]);
          TD_value[Component][IXX + Step - 1] += mTYPE(DP_LCR[1][J + j]);
          continue;
        case 0:
          TD_value[Component][IXX + Step - 1] += mTYPE(DP_LCR[0][J + j]);
          continue;
        }
      }
    }
  }

  for (i = 0; i < Num_component; i++) {
    for (j = 0; j < Num_data; j++) {
      TD_value[i][j] = mTYPE(TD_value[i][j] * pix_size);
    }
  }

  TD.transpose();
  free(TD_value);

  printf("*** transmission_distance *** :   end\n");
}

template<typename mTYPE1, typename mTYPE2>
void Proj_trape(const DP_trapezoid& DP, const MATRIX<mTYPE1>& image, double pix_size, MATRIX<mTYPE2>& sinogram) {
  //    printf("*** Proj_trape *** :\r");

  const int Num_angle = DP.Num_angle();
  const int Step = DP.Step();
  const int Height = DP.Height();
  const int Width = DP.Width();
  const int Num_pixel = Height*Width;
  const int sinogram_data = Num_angle*Step;

  if (Num_angle != sinogram.Row() || Step != sinogram.Col()) sinogram.resize(Num_angle, Step);

  int i, j, k, I, K, J, IXX;
  int ixx;

  short *DP_ixx;
  DP_ixx = DP.dp_Xs.pVALUE();

  double *DP_LCR[3];//DPをdoubleにするならここを対応させること。
  DP_LCR[0] = DP.dp_L.pVALUE();
  DP_LCR[1] = DP.dp_C.pVALUE();
  DP_LCR[2] = DP.dp_R.pVALUE();

  mTYPE1 *img;
  img = image.pVALUE();

  mTYPE2 *sinogram_value;
  sinogram_value = sinogram.pVALUE();


  for (i = 0; i < sinogram_data; i++) {
    sinogram_value[i] = 0.0;
  }

  for (k = 0; k < Num_angle; k++) {
    K = k*Num_pixel;
    IXX = k*Step;
    for (i = 0; i < Height; i++) {
      I = i*Width;
      J = K + i*Width;
      for (j = 0; j < Width; j++) {
        ixx = DP_ixx[J + j];

        if (ixx<-1 || ixx>Step)continue;

        if (0 < ixx && ixx < Step - 1) {
          sinogram_value[IXX + ixx - 1] += mTYPE2(DP_LCR[0][J + j] * img[I + j]);
          sinogram_value[IXX + ixx] += mTYPE2(DP_LCR[1][J + j] * img[I + j]);
          sinogram_value[IXX + ixx + 1] += mTYPE2(DP_LCR[2][J + j] * img[I + j]);
          continue;
        }

        switch (ixx) {
        case -1:
          sinogram_value[IXX] += mTYPE2(DP_LCR[2][J + j] * img[I + j]);
          continue;
        case 0:
          sinogram_value[IXX] += mTYPE2(DP_LCR[1][J + j] * img[I + j]);
          sinogram_value[IXX + 1] += mTYPE2(DP_LCR[2][J + j] * img[I + j]);
          continue;
        }

        switch (ixx - Step) {
        case -1:
          sinogram_value[IXX + Step - 2] += mTYPE2(DP_LCR[0][J + j] * img[I + j]);
          sinogram_value[IXX + Step - 1] += mTYPE2(DP_LCR[1][J + j] * img[I + j]);
          continue;
        case 0:
          sinogram_value[IXX + Step - 1] += mTYPE2(DP_LCR[0][J + j] * img[I + j]);
          continue;
        }
      }
    }
  }

  for (i = 0; i < sinogram_data; i++) {
    sinogram_value[i] = mTYPE2(sinogram_value[i] * pix_size);
  }

  //  printf("*** Proj_trape *** :   end\n");
}

template<typename mTYPE2, typename mTYPE1>
void BackProj_trape(const DP_trapezoid& DP, const MATRIX<mTYPE2>& sinogram, MATRIX<mTYPE1>& image)
{
  //  printf("*** BackProj_trape *** :\r");

  const int Num_angle = DP.Num_angle();
  const int Step = DP.Step();
  const int Height = DP.Height();
  const int Width = DP.Width();
  const int Num_pixel = Height*Width;

  if (Height != image.Row() || Width != image.Col()) image.resize(Height, Width);

  int i, j, k, I, J, K;
  int ixx, IXX;

  short *DP_ixx;
  DP_ixx = DP.dp_Xs.pVALUE();

  double *DP_LCR[3];//DPをdoubleにするならここを対応させること。
  DP_LCR[0] = DP.dp_L.pVALUE();
  DP_LCR[1] = DP.dp_C.pVALUE();
  DP_LCR[2] = DP.dp_R.pVALUE();

  mTYPE1 *img;
  img = image.pVALUE();

  mTYPE2 *sinogram_value;
  sinogram_value = sinogram.pVALUE();

  for (i = 0; i < Num_pixel; i++) {
    img[i] = 0.0;
  }

  for (k = 0; k < Num_angle; k++) {
    K = k*Num_pixel;
    IXX = k*Step;
    for (i = 0; i < Height; i++) {
      I = i*Width;
      J = K + i*Width;
      for (j = 0; j < Width; j++) {
        ixx = DP_ixx[J + j];

        if (ixx<-1 || ixx>Step)continue;

        if (0 < ixx && ixx < Step - 1) {
          img[I + j] += mTYPE1(DP_LCR[0][J + j] * sinogram_value[IXX + ixx - 1]);
          img[I + j] += mTYPE1(DP_LCR[1][J + j] * sinogram_value[IXX + ixx]);
          img[I + j] += mTYPE1(DP_LCR[2][J + j] * sinogram_value[IXX + ixx + 1]);
          continue;
        }

        switch (ixx) {
        case -1:
          img[I + j] += mTYPE1(DP_LCR[2][J + j] * sinogram_value[IXX + ixx + 1]);
          continue;
        case 0:
          img[I + j] += mTYPE1(DP_LCR[1][J + j] * sinogram_value[IXX + ixx]);
          img[I + j] += mTYPE1(DP_LCR[2][J + j] * sinogram_value[IXX + ixx + 1]);
          continue;
        }

        switch (ixx - Step) {
        case -1:
          img[I + j] += mTYPE1(DP_LCR[0][J + j] * sinogram_value[IXX + ixx - 1]);
          img[I + j] += mTYPE1(DP_LCR[1][J + j] * sinogram_value[IXX + ixx]);
          continue;
        case 0:
          img[I + j] += mTYPE1(DP_LCR[0][J + j] * sinogram_value[IXX + ixx - 1]);
          continue;
        }
      }
    }
  }

  for (i = 0; i < Num_pixel; i++) {
    img[i] /= Num_angle;
  }

  //  printf("*** BackProj_trape *** :   end\n");
}

template<typename mTYPE1, typename mTYPE2>
void Proj_trape_frag(const DP_trapezoid_fragment& DP, const MATRIX<mTYPE1>& image, double pix_size, MATRIX<mTYPE2>& sinogram, double threshold=-1.0) {
  //    printf("*** Proj_trape *** :\r");

  int Num_angle = DP.Num_angle();
  int Step = DP.Step();
  int Height = DP.Height();
  int Width = DP.Width();
  int Num_pixel = Height*Width;
  int sinogram_data = Num_angle*Step;

  if (Num_angle != sinogram.Row() || Step != sinogram.Col()) sinogram.resize(Num_angle, Step);

  int i, j, k, I, K, J, IXX;
  int ixx;
  double x, y, xx, XX;
  double si, co, a, b, c, ixx_side, d;;

  double y0 = 0.5*Height - 0.5;
  double x0 = -0.5*Width + 0.5;
  double step0 = Step / 2.0 - 0.5;

  double dp_LCR[3];

  mTYPE1 *img;
  img = image.pVALUE();

  sinogram.constant(0.0);
  mTYPE2 *sinogram_value;
  sinogram_value = sinogram.pVALUE();

  for (k = 0; k < Num_angle; k++) {
    K = k*Num_pixel;
    IXX = k*Step;

    a = DP.dp.out(k + 1, 1);
    b = DP.dp.out(k + 1, 2);
    c = DP.dp.out(k + 1, 3);
    si = DP.dp.out(k + 1, 4);
    co = DP.dp.out(k + 1, 5);

    for (i = 0; i < Height; i++) {
      I = i*Width;

      J = K + i*Width;
      y = y0 - i;//image中心を原点とする座標
      XX = y*si + step0;
      for (j = 0; j < Width; j++) {
        if (img[I + j] < threshold) {
          img[I + j] = 0.0;
          continue;
        }//順投影するピクセル値が閾値以下のとき計算を省略する

        x = x0 + j;//image中心を原点とする座標
        xx = x*co + XX;//1つ目のステップの中心をゼロする座標で見たピクセル中心の投影位置
        ixx = int(floor(xx + 0.5));//ピクセル中心の投影位置が、いくつ目のステップの中にあるかを計算。
        if (ixx<-1 || ixx>Step) continue;//ピクセル中心の投影位置が投影の領域外である場合は以降の計算をパスする

//検出確率計算
        dp_LCR[0] = dp_LCR[1] = dp_LCR[2] = 0;

        //ひとつ前のステップに入る検出確率(dp_LCR[0])
        ixx_side = ixx - 0.5;//ひとつ前のステップとの境界の座標
        d = ixx_side - (xx - 0.5*(a - b));//境界と台形の肩との距離
        if (d >= 0) {//境界のほうが右
          dp_LCR[0] = (0.5*b + d)*c;
        }
        else if ((-d) < b) {//境界のほうが右だが台形の足の位置よりは右
          dp_LCR[0] = 0.5*c*(b + d) / b*(b + d);
        }

        //ひとつ後のステップに入る検出確率(dp_LCR[1])
        ixx_side = ixx + 0.5;//ひとつ後のステップとの境界の座標
        d = (xx + 0.5*(a - b)) - ixx_side;//境界と台形の肩との距離
        if (d >= 0) {//境界のほうが左
          dp_LCR[2] = (0.5*b + d)*c;
        }
        else if ((-d) < b) {
          dp_LCR[2] = 0.5*c*(b + d) / b*(b + d);
        }

        //投影位置の検出器に入る検出確率(dp_LCR[2])
        dp_LCR[1] = 1 - dp_LCR[0] - dp_LCR[2];

        //検出確率計算ここまで
        //投影計算

        if (0 < ixx && ixx < Step - 1) {
          sinogram_value[IXX + ixx - 1] += mTYPE2(dp_LCR[0] * img[I + j]);
          sinogram_value[IXX + ixx] += mTYPE2(dp_LCR[1] * img[I + j]);
          sinogram_value[IXX + ixx + 1] += mTYPE2(dp_LCR[2] * img[I + j]);
          continue;
        }

        switch (ixx) {
        case -1:
          sinogram_value[IXX] += mTYPE2(dp_LCR[2] * img[I + j]);
          continue;
        case 0:
          sinogram_value[IXX] += mTYPE2(dp_LCR[1] * img[I + j]);
          sinogram_value[IXX + 1] += mTYPE2(dp_LCR[2] * img[I + j]);
          continue;
        }

        switch (ixx - Step) {
        case -1:
          sinogram_value[IXX + Step - 2] += mTYPE2(dp_LCR[0] * img[I + j]);
          sinogram_value[IXX + Step - 1] += mTYPE2(dp_LCR[1] * img[I + j]);
          continue;
        case 0:
          sinogram_value[IXX + Step - 1] += mTYPE2(dp_LCR[0] * img[I + j]);
          continue;
          //投影計算ここまで
        }
      }
    }
  }

  for (i = 0; i < sinogram_data; i++) {
    sinogram_value[i] = mTYPE2(sinogram_value[i] * pix_size);
  }

  //  printf("*** Proj_trape *** :   end\n");
}

template<typename mTYPE2, typename mTYPE1>
void BackProj_trape_frag(const DP_trapezoid_fragment& DP, const MATRIX<mTYPE2>& sinogram, MATRIX<mTYPE1>& image) {
  //  printf("*** BackProj_trape *** :\r");

  const int Num_angle = DP.Num_angle();
  const int Step = DP.Step();
  const int Height = DP.Height();
  const int Width = DP.Width();
  const int Num_pixel = Height*Width;

  if (Height != image.Row() || Width != image.Col()) image.resize(Height, Width);

  int i, j, k, I, J, K;
  int ixx, IXX;
  double x, y, xx, XX;
  double si, co, a, b, c;
  double ixx_side, d;

  double y0 = 0.5*Height - 0.5;
  double x0 = -0.5*Width + 0.5;
  double step0 = Step / 2.0 - 0.5;

  double dp_LCR[3];

  image.constant(0.0);
  mTYPE1 *img;
  img = image.pVALUE();

  mTYPE2 *sinogram_value;
  sinogram_value = sinogram.pVALUE();

  for (k = 0; k < Num_angle; k++) {
    K = k*Num_pixel;
    IXX = k*Step;//sinogramのステップのオフセット部分

    a = DP.dp.out(k + 1, 1);
    b = DP.dp.out(k + 1, 2);
    c = DP.dp.out(k + 1, 3);
    si = DP.dp.out(k + 1, 4);
    co = DP.dp.out(k + 1, 5);
    for (i = 0; i < Height; i++) {
      I = i*Width;
      J = K + i*Width;

      y = y0 - i;//image中心を原点とする座標
      XX = y*si + step0;
      for (j = 0; j < Width; j++) {
        x = x0 + j;//image中心を原点とする座標
        xx = x*co + XX;//1つ目のステップの中心をゼロする座標で見たピクセル中心の投影位置
        ixx = int(floor(xx + 0.5));//ピクセル中心の投影位置が、いくつ目のステップの中にあるかを計算。
        if (ixx<-1 || ixx>Step) continue;//ピクセル中心の投影位置が投影の領域外である場合は以降の計算をパスする

 //検出確率計算
        dp_LCR[0] = dp_LCR[1] = dp_LCR[2] = 0;

        //ひとつ前のステップに入る検出確率(dp_LCR[0])
        ixx_side = ixx - 0.5;//ひとつ前のステップとの境界の座標
        d = ixx_side - (xx - 0.5*(a - b));//境界と台形の肩との距離
        if (d >= 0) {//境界のほうが右
          dp_LCR[0] = (0.5*b + d)*c;
        }
        else if ((-d) < b) {//境界のほうが右だが台形の足の位置よりは右
          dp_LCR[0] = 0.5*c*(b + d) / b*(b + d);
        }

        //ひとつ後のステップに入る検出確率(dp_LCR[1])
        ixx_side = ixx + 0.5;//ひとつ後のステップとの境界の座標
        d = (xx + 0.5*(a - b)) - ixx_side;//境界と台形の肩との距離
        if (d >= 0) {//境界のほうが左
          dp_LCR[2] = (0.5*b + d)*c;
        }
        else if ((-d) < b) {
          dp_LCR[2] = 0.5*c*(b + d) / b*(b + d);
        }

        //投影位置の検出器に入る検出確率(dp_LCR[2])
        dp_LCR[1] = 1 - dp_LCR[0] - dp_LCR[2];

//検出確率計算ここまで
//逆投影計算
        if (0 < ixx && ixx < Step - 1) {
          img[I + j] += mTYPE1(dp_LCR[0] * sinogram_value[IXX + ixx - 1]);
          img[I + j] += mTYPE1(dp_LCR[1] * sinogram_value[IXX + ixx]);
          img[I + j] += mTYPE1(dp_LCR[2] * sinogram_value[IXX + ixx + 1]);
          continue;
        }

        switch (ixx) {
        case -1:
          img[I + j] += mTYPE1(dp_LCR[2] * sinogram_value[IXX + ixx + 1]);
          continue;
        case 0:
          img[I + j] += mTYPE1(dp_LCR[1] * sinogram_value[IXX + ixx]);
          img[I + j] += mTYPE1(dp_LCR[2] * sinogram_value[IXX + ixx + 1]);
          continue;
        }

        switch (ixx - Step) {
        case -1:
          img[I + j] += mTYPE1(dp_LCR[0] * sinogram_value[IXX + ixx - 1]);
          img[I + j] += mTYPE1(dp_LCR[1] * sinogram_value[IXX + ixx]);
          continue;
        case 0:
          img[I + j] += mTYPE1(dp_LCR[0] * sinogram_value[IXX + ixx - 1]);
          continue;
        }
        //逆投影計算ここまで
      }
    }
  }

  for (i = 0; i < Num_pixel; i++) {
    img[i] /= Num_angle;
  }

  //  printf("*** BackProj_trape *** :   end\n");

}



/* 少し役立つ関数群 */
template<typename mTYPE>
MATRIX<mTYPE> click_angle(const MATRIX<mTYPE>& sinogram, int Num, int scan)//サイノグラムの投影を入力した[行番号]から始まるようにずらす。行方向に角度を取っていること。scan=1でハーフ、2でフル
{
  int Width = sinogram.Col();
  int Height = sinogram.Row();

  if (Num == 1) {
    return sinogram;
  }
  else if (Num < 2 || Height < Num) {
    printf("cannot access. click_angle(%d,scan=%d)\n", Num, scan);
    system("pause");
    exit(EXIT_FAILURE);
  }

  MATRIX<mTYPE> ans(Height, Width);

  MATRIX<mTYPE> buf;
  buf = sinogram.pull(Num, 1, Height - Num + 1, Width);

  ans.put(1, 1, buf);

  buf = sinogram.pull(1, 1, Num - 1, Width);
  if (scan == 1) buf.invert_row();

  ans.put(Height - Num + 2, 1, buf);

  return ans;
}

template<typename mTYPE>
MATRIX<mTYPE> devide_each(const MATRIX<mTYPE>& sinogram)//投影方向ごとの最大値で各投影を割った行列を返す。
{
  int buf_row = sinogram.Row();
  int buf_col = sinogram.Col();

  MATRIX<mTYPE> buf;
  MATRIX<mTYPE> ans(buf_row, buf_col);

  int i;
  for (i = 1; i <= buf_row; i++)
  {
    buf = sinogram.pullrow(i);
    buf = m_cast<mTYPE>(buf / buf.max());
    ans.putrow(i, buf);
  }

  return ans;
}