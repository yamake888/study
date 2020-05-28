#pragma once

// �ŏI�X�V 19,JUN,2018
//���e�������̑����č\���ɂ��Ή�

#define _CRT_SECURE_NO_WARNINGS //�R���p�C������fopen�Ȃǂ̌x���𖳎����邽��
#define _USE_MATH_DEFINES//�~����M_PI���g������
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "matrix_tool.h"

/*
�摜�č\���̂��߂̊֐��Q
�@CT2.0�̏Ă�����
�A�R���{�����[�V�����č\���@(FBP�܂���CT2.0)�Ɋ֌W�����T�u���[�`���Q
�BML-EM�摜�č\���@�Ɋ֌W�����N���X�ƃT�u���[�`���Q
�ԊO �����𗧂֐��Q (�T�C�m�O�������p�x�����ɂ��炷�Ƃ��B)
matrix_tool.h�̃C���N���[�h���K�v
*/

/* CT2.0�̏Ă����� �v���g�^�C�v�錾 */

template<typename mTYPE1, typename mTYPE2>
void BackProjection_CT2p0(const MATRIX<mTYPE2>&, MATRIX<mTYPE1>&);//CT2.0�̋t���e��MARIX�ŏ������������́B�A�N�Z�X�ᔽ�œ��삵�Ȃ��B

template<typename mTYPE1, typename mTYPE2>
void BackProjection_CT2p0_mod(const MATRIX<mTYPE2>&, MATRIX<mTYPE1>&);//�A�N�Z�X�ᔽ�����������

template<typename mTYPE>
void ConvolutionIntegral_CT2p0(MATRIX<mTYPE>&, double);//�T�C�m�O�����ɃR���{�����[�V������������֐���MATRIX�ŏ�������������

double SLConFun_CT2p0(int, double);//�R���{�����[�V�����Ŏg����B


/* �R���{�����[�V�����č\���@(FBP)�Ɋ֌W�����T�u���[�`���Q �v���g�^�C�v�錾 */

template<typename mTYPE1, typename mTYPE2>
void BackProj_simpl(const MATRIX<mTYPE1>&, MATRIX<mTYPE2>&, int);//���o�m���v�Z���P���Ȃ��̂ŒP���t���e���s��

template<typename mTYPE>
MATRIX<mTYPE> ConvInte(double(*conv_func)(int, double), const MATRIX<mTYPE>&, double);//�T�C�m�O�����ɑ΂��Ďw�肵���t�B���^�ŏd���ϕ���������

double Shepp_Logan_cf(int, double);//�w�肳���t�B���^�̂ЂƂBShepp Logan convolution filter


/* ML-EM�摜�č\���@�Ɋ֌W�����N���X�ƃT�u���[�`���Q �v���g�^�C�v�錾 */

//��`�Ōv�Z�������o�m�����Ǘ�����N���X�B�g���Ƃ��ɂ͎n�߂ɑ��������̂��߂ɁA�܂��͂��̃N���X�̕ϐ���錾���邱�ƁB
class DP_trapezoid
{
public:
  MATRIX<short> dp_Xs;//���e�ʒu�̒��S��������X�e�b�v�ԍ���[0]����[�X�e�b�v��-1]�܂ł̐����ŊǗ�������B
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

//��`�Ōv�Z�������o�m�����Ǘ�����N���X�B���o�m����ێ����邽�߂̃��������傫������ꍇ�Ɏg���B
class DP_trapezoid_fragment
{
public:
  MATRIX<double> dp;

  int Num_angle()const { return dp_Num_angle; }
  int Step()const { return dp_Step; }
  int Height()const { return dp_Height; }
  int Width()const { return dp_Width; }

  DP_trapezoid_fragment(MATRIX<double> angle, int Step, int Height, int Width)//angle��degree�ŗ^���邱�ƁB
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

void calcDP_trapezoid(DP_trapezoid&, int scan = 1);//���o�m�����v�Z����B

template<typename mTYPE1, typename mTYPE2>
void ML_EM(const DP_trapezoid&, const MATRIX<mTYPE2>&, MATRIX<mTYPE1>&, double, int);//���e�A�t���e���Z�b�g�ɂ���ML-EM�摜�č\�����s��

template<typename mTYPE1, typename mTYPE2>
void ML_EM_fragment(const DP_trapezoid_fragment&, const MATRIX<mTYPE2>&, MATRIX<mTYPE1>&, double, int, double);//���e�A�t���e���Z�b�g�ɂ���ML-EM�摜�č\�����s��

template<typename mTYPE>
void transmission_distance(const DP_trapezoid&, const MATRIX<short>&, int, double, MATRIX<mTYPE>&);//�p�ӂ����摜�̊e�v�f�ɂ��ē��ߋ����̃T�C�m�O���������B��������s�N�Z���ɂ�0�A�v�Z�������s�N�Z���ɂ͗v�f���Ƃ�1���珇�ɐ����̒ʂ��ԍ���t���邱�ƁB����ɂ���ƃA�N�Z�X�G���[���N����B

template<typename mTYPE1, typename mTYPE2>
void Proj_trape(const DP_trapezoid&, const MATRIX<mTYPE1>&, double, MATRIX<mTYPE2>&);//�摜�𓊉e����B

template<typename mTYPE2, typename mTYPE1>
void BackProj_trape(const DP_trapezoid&, const MATRIX<mTYPE2>&, MATRIX<mTYPE1>&);//�T�C�m�O�������t���e����B

template<typename mTYPE1, typename mTYPE2>
void Proj_trape_frag(const DP_trapezoid_fragment&, const MATRIX<mTYPE1>&, double, MATRIX<mTYPE2>&, double);//�摜�𓊉e����B

template<typename mTYPE2, typename mTYPE1>
void BackProj_trape_frag(const DP_trapezoid_fragment&, const MATRIX<mTYPE2>&, MATRIX<mTYPE1>&);//�T�C�m�O�������t���e����B


/* �����𗧂֐��Q */
template<typename mTYPE>
MATRIX<mTYPE> click_angle(const MATRIX<mTYPE>&, int, int scan = 1);//�T�C�m�O��������͂���[�s�ԍ�]����n�܂�悤�ɂ��炷�B�s�����Ɋp�x������Ă��邱�ƁBscan=1�Ńn�[�t�A2�Ńt��
template<typename mTYPE>
MATRIX<mTYPE> devide_each(const MATRIX<mTYPE>&);//���e�������Ƃ̍ő�l�Ŋe���e���������s���Ԃ��B



/* ��������͊֐��̒�` */

/* CT2.0�̏Ă����� */
template<typename mTYPE1, typename mTYPE2>
void BackProjection_CT2p0(const MATRIX<mTYPE2> &sinogram, MATRIX<mTYPE1> &image)
//  CT2.0��void bp()�݂̂�MATRIX_tool���g���ď����ʂ������́B
//  CT2.0�ł�sinogram�́A�ΐ�������(void readf(char string[80])�ɂ�)�A�R���{�����[�V����������Ă���(void con()�ɂ�)
//  ���̃T�u�֐����g���ꍇ�A�����̏����͕ʂɍs�����ƁB
//  ��1�ŃA�N�Z�X�ᔽ���N����B  
{
  double PI = atan(1.0) * 4;//CT2.0�ł͉~������M_PI���g�킸�A���̂悤�Ɍv�Z���Ă���B
  int k, l;//�摜����[�����_�Ƃ����Ƃ��̉�f�̈ʒu��\���ϐ�
  int n;//���e�����̔ԍ���\���ϐ�
  double x, y;//�摜���������_�Ƃ����Ƃ��̉�f�̈ʒu��\���ϐ�
  double s;//�X�e�b�v�̒��S�����_�Ƃ��ĉ�f�̈ʒu�𓊉e�X�e�b�v�ɓ��e�����ʒu�B
  double g;//s�̏�������
  double m01;//s�̐�������
  int m0;//m01��int�^�ɕύX���Ĉ������߂̕ϐ�
  double ff;//�t���e�̉��Z�v�Z�p�̕ϐ�

  int N = sinogram.Col();//N�̓X�e�b�v���B�T�C�m�O�����̗񐔂ɓ������B
  int M = sinogram.Row();//M�͓��e�������B�T�C�m�O�����̍s���ɓ������B
  image.resize(N, N);//�č\���摜���i�[����s��BCT2.0�ł�f[N][N]�ɑ�������B

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
        m0 = int(m0 + N / 2.0);//CT2.0�ɂ̓L���X�g�͏�����Ă��Ȃ����A�x�����ז��Ȃ̂ŏ��������B�Öق̃L���X�g�Ɠ����Ȃ̂ŉe���͂Ȃ��B
        ff = ff + (1 - g)*sinogram.out(n + 1, m0 + 1) + g*sinogram.out(n + 1, m0 + 2);//sinogram�̍s��͍���[�̗v�f��(1,1)�Ȃ̂ŁA+1���炵���B��1
      }
      image.in(k + 1, l + 1, ff*(1 / (2 * PI))*(PI / (M - 0.0)));//�P���ɍl�����2*PI�ł͂Ȃ�PI�ɂ��Ă����΁A���e�������ł̕��ω��Ǝv����̂����B
      //������ConvolutionIntegral_CT2p0��sinogram��2�{�ɑ傫���v�Z����Ă��邱�Ƃ͊֌W���Ȃ����낤���B
    }
  }

}

template<typename mTYPE1, typename mTYPE2>
void BackProjection_CT2p0_mod(const MATRIX<mTYPE2> &sinogram, MATRIX<mTYPE1> &image)
//  CT2.0��void bp()�݂̂�MATRIX_tool���g���ď����ʂ������́B
//  CT2.0�ł�sinogram�́A�ΐ�������(void readf(char string[80])�ɂ�)�A�R���{�����[�V����������Ă���(void con()�ɂ�)
//  ���̃T�u�֐����g���ꍇ�A�����̏����͕ʂɍs�����ƁB
//  ��1�ŃA�N�Z�X�ᔽ���N����B  
{
  double PI = atan(1.0) * 4;//CT2.0�ł͉~������M_PI���g�킸�A���̂悤�Ɍv�Z���Ă���B
  int k, l;//�摜����[�����_�Ƃ����Ƃ��̉�f�̈ʒu��\���ϐ�
  int n;//���e�����̔ԍ���\���ϐ�
  double x, y;//�摜���������_�Ƃ����Ƃ��̉�f�̈ʒu��\���ϐ�
  double s;//�X�e�b�v�̒��S�����_�Ƃ��ĉ�f�̈ʒu�𓊉e�X�e�b�v�ɓ��e�����ʒu�B
  double g;//s�̏�������
  double m01;//s�̐�������
  int m0;//m01��int�^�ɕύX���Ĉ������߂̕ϐ�
  double ff;//�t���e�̉��Z�v�Z�p�̕ϐ�

  int N = sinogram.Col();//N�̓X�e�b�v���B�T�C�m�O�����̗񐔂ɓ������B
  int M = sinogram.Row();//M�͓��e�������B�T�C�m�O�����̍s���ɓ������B
  image.resize(N, N);//�č\���摜���i�[����s��BCT2.0�ł�f[N][N]�ɑ�������B

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
        m0 = int(m0 + N / 2.0);//CT2.0�ɂ̓L���X�g�͏�����Ă��Ȃ����A�x�����ז��Ȃ̂ŏ��������B�Öق̃L���X�g�Ɠ����Ȃ̂ŉe���͂Ȃ��B
        if (m0 >= 0 && m0 < N - 1)
        {
          ff = ff + (1 - g)*sinogram.out(n + 1, m0 + 1) + g*sinogram.out(n + 1, m0 + 2);//sinogram�̍s��͍���[�̗v�f��(1,1)�Ȃ̂ŁA+1���炵���B��1
        }
      }
      image.in(k + 1, l + 1, ff*(1 / (PI))*(PI / (M - 0.0)));//�P���ɍl�����2*PI�ł͂Ȃ�PI�ɂ��Ă����΁A���e�������ł̕��ω��B
      //������ConvolutionIntegral_CT2p0��sinogram��2�{�ɑ傫���v�Z����Ă��邱�Ƃ͊֌W���Ȃ����낤���B
    }
  }

}

template<typename mTYPE>
void ConvolutionIntegral_CT2p0(MATRIX<mTYPE> &sinogram, double sampling_interval)
//  CT2.0��con()��MATRIX_tool���g���ď����ʂ�������
//  sinogram���K�v�Ȃ̂œn���Ă���B
//  sampling_interval�̓O���[�o���ϐ��Œ�`����Ă���dr�̂��ƂŁA1���e�̕�(cm)���X�e�b�v���Ŋ������l��\���B
{
  int n, md;//���ꂼ�꓊�e�����A�e���e�ł̃X�e�b�v��\���B
  int m;//1�̓��e�Œ[����d���ϕ����s�����߂Ɏg���B
  double q_temp = 0.0;//�d���ϕ��̒l�����Z���邽�߂̕ϐ�

  int N = sinogram.Col();//N�̓X�e�b�v���B�T�C�m�O�����̗񐔂ɓ������B
  int M = sinogram.Row();//M�͓��e�������B�T�C�m�O�����̍s���ɓ������B
  MATRIX copy = MATRIX(sinogram);//copy�����ɂ��ďd���ϕ����ꂽsinogram���v�Z����B
  sinogram.constant(0.0);//�d���ϕ��̌��ʂ��i�[����̂ł��ׂĂ̐��l��0�ɂ��Ă����B���̏�����CT_2.0�ł͌�������Ȃ������B
  double dr = sampling_interval;

  for (n = 0; n < M; n++)
  {
    for (md = 0; md < N; md++)
    {
      q_temp = 0;
      for (m = 0; m < N; m++)
      {
        q_temp = q_temp + dr*(sinogram.out(n + 1, md + 1) + SLConFun_CT2p0(md - m, dr)*copy.out(n + 1, m + 1));
        //SLConFun_CT2p0�̃R�����g�ɂ����������A���q��2�{�傫���̂�sinogram�̒l��2�{�ɂȂ��Ă���B
      }
      sinogram.in(n + 1, md + 1, mTYPE(q_temp));
    }
  }
}

double SLConFun_CT2p0(int m, double sampling_interval)
//  CT2.0��h(int m)�������ʂ������́Bcon()�Ŏg���Ă���B
//  sampling_interval�̓O���[�o���ϐ��Œ�`����Ă���dr�̂��ƂŁA1���e�̕�(cm)���X�e�b�v���Ŋ������l��\���B
{
  double PI = atan(1.0) * 4;//CT2.0�ł͉~������M_PI���g�킸�A���̂悤�Ɍv�Z���Ă���B
  double dr = sampling_interval;

  double hh = 0.0;
  hh = 4 / (PI*dr*dr*(1 - 4 * m*m));
  //ttp://fiber.shinshu-u.ac.jp/koseki/publication/2003/bio2003-3.pdf���ł͕��q��2�Ȃ̂����B

  return hh;
}

/* �R���{�����[�V�����č\���@(FBP)�Ɋ֌W�����T�u�֐��Q */
template<typename mTYPE1, typename mTYPE2>
void BackProj_simpl(const MATRIX<mTYPE1> &sinogram, MATRIX<mTYPE2> &image, int scan)
/*
scan=1�Ńn�[�t�X�L�����č\���Bscan=2�Ńt���X�L�����č\���B
�s�N�Z���̓��e����1�X�e�b�v���ƌ��Ȃ��Č��o�m�����v�Z����P���t���e�v���O�����Bsimple�ƌĂԂ��Ƃɂ���B
sinogram�͍s�����e�����A�񂪃X�e�b�v���ƂȂ�悤�ɏ������邱�ƁB������W���̓��e�̘a(-log(I/I0))�ɂ��Ă������ƁB
(�d�l)
�@�č\���摜�s�N�Z���T�C�Y(�c�̒����Ɖ��̒���)�ƃX�e�b�v�̕������ׂē������B
 �@�č\���摜�̏c�A���s�N�Z�����A�X�e�b�v�������ׂē������B���������ׂċ����ł���ł��Ή����Ă���B
  �@�擪�̃X�e�b�v�A�Ō���̃X�e�b�v�ł̌��o���������Ȃ��悤�ɍ���Ă���B
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

        if (ixx < -1 || Step - 1 < ixx) continue;//���o������dp1�Adp2�Ƃ���0�ɂȂ�̂�continue    

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
  �R���{�����[�V�������s���֐��B
  �����ɃR���{�����[�V�����֐�(�T�u�֐��̖��O�̕���)��^����ƁA������g���Čv�Z����B
  sinogram�͍s�����e�����A�񂪃X�e�b�v���ƂȂ�悤�ɏ������邱�ƁB������W���̓��e�̘a(-log(I/I0))�ɂ��Ă������ƁB
  step_interval�͓��e�X�e�b�v�̃X�e�b�v�ԋ���(cm)�B�����̒P�ʂ�cm�Ȃ̂ŁA�ŏI�I�ɍč\���摜�̒P�ʂ�(1/cm)�ɂȂ�B
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
  ConvInte�Ɏg���R���{�����[�V�����֐��B
  pixel_distance�̓s�N�Z�����Ő������Ƃ��́A���ڃX�e�b�v�Ɖe�����l����X�e�b�v�Ƃ̋����B
  step_interval�͓��e�X�e�b�v�̃X�e�b�v�ԋ���(cm)�B�����̒P�ʂ�cm�Ȃ̂ŁA�ŏI�I�ɍč\���摜�̒P�ʂ�(1/cm)�ɂȂ�
  shepp-logan�t�B���^�ɂ��Ē��ׂ�ƌW���̕t�������l�X�ō�������B����̓t�[���G�ϊ��̌W�����ǂ��ɂ��邩�̖��ŁA
  �摜�č\�����̍Ō�̓��e�������Ŋ���v�Z�Ɠs�������킹�Ȃ��Ƃ����Ȃ��B
  BackProj_A�̓��e�������Ŋ���v�Z�ŕ����������̂܂�(�΂��|������)�v�Z����
  bp_img=bp_img/M;
  ���ƂƂ��A���L�̒ʂ�Ɍv�Z����B2016/11/11hamaguchi
  ����������̂�h=2.0/(M_PI*M_PI*step_interval*step_interval*(1.0-4.0*(pixel_distance*pixel_distance)));�����A�����������ƁB
  */
  double h;

  h = 2.0 / (M_PI*step_interval*step_interval*(1.0 - 4.0*(pixel_distance*pixel_distance)));

  return h;
}

// Ram_Lak_cf�͖�����
double Ram_Lak_cf(int pixel_distance, double step_interval)
//ConvInte�Ɏg���R���{�����[�V�����֐��B
//pixel_distance�̓s�N�Z�����Ő������Ƃ��́A���ڃX�e�b�v�Ɖe�����l����X�e�b�v�Ƃ̋����B
//step_interval�͓��e�X�e�b�v�̃X�e�b�v�ԋ���(cm)�B�����̒P�ʂ�cm�Ȃ̂ŁA�ŏI�I�ɍč\���摜�̒P�ʂ�(1/cm)�ɂȂ�
//Ram_Lak�t�B���^�ɂ��Ă��A���ׂ�ƌW���̕t�������l�X�ō������邩������Ȃ��B
//Shepp_Logan_cf�Ɠ��l�ɋt���e�̃v���O�����Ƃ̓s�����l���A���L�̒ʂ�Ɍv�Z����B2016/11/11hamaguchi

{
  double h;

  if (pixel_distance % 2 == 1)//�
  {
    h = -1.0 / (pixel_distance*pixel_distance*M_PI*step_interval*step_interval);
    return h;
  }
  else if ((pixel_distance / 2) != 0)//�[���ȊO�̋���
  {
    h = 0.0;
    return h;
  }
  else
  {
    h = M_PI / (4.0*step_interval*step_interval);//�[��
    return h;
  }
}


/* ML-EM�摜�č\���@�Ɋ֌W�����N���X�ƃT�u�֐��Q */

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

  double *DP_LCR[3];//DP��double�ɂ���Ȃ炱����Ή������邱�ƁB
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
        xx = x*co + XX;//1�ڂ̃X�e�b�v�̒��S���[��������W�Ō����s�N�Z�����S�̓��e�ʒu
        ixx = int(floor(xx + 0.5));//�s�N�Z�����S�̓��e�ʒu���A�����ڂ̃X�e�b�v�̒��ɂ��邩���v�Z�B

        if (ixx<-1 || ixx>Step) {
          DP_ixx[J + j] = -2;
          continue;
        }
        dp_LCR[0] = dp_LCR[1] = dp_LCR[2] = 0.0;

        //�ЂƂO�̃X�e�b�v�ɓ��錟�o�m��(dp_LCR[0])
        ixx_side = ixx - 0.5;//�ЂƂO�̃X�e�b�v�Ƃ̋��E�̍��W
        d = ixx_side - (xx - 0.5*(a - b));//���E�Ƒ�`�̌��Ƃ̋���
        if (d >= 0) {//���E�̂ق����E
          dp_LCR[0] = (0.5*b + d)*c;
        }
        else if ((-d) < b) {//���E�̂ق����E������`�̑��̈ʒu���͉E
          dp_LCR[0] = 0.5*c*(b + d) / b*(b + d);
        }

        //�ЂƂ�̃X�e�b�v�ɓ��錟�o�m��(dp_LCR[1])
        ixx_side = ixx + 0.5;//�ЂƂ�̃X�e�b�v�Ƃ̋��E�̍��W
        d = (xx + 0.5*(a - b)) - ixx_side;//���E�Ƒ�`�̌��Ƃ̋���
        if (d >= 0) {//���E�̂ق�����
          dp_LCR[2] = (0.5*b + d)*c;
        }
        else if ((-d) < b) {
          dp_LCR[2] = 0.5*c*(b + d) / b*(b + d);
        }

        //���e�ʒu�̌��o��ɓ��錟�o�m��(dp_LCR[2])
        dp_LCR[1] = 1 - dp_LCR[0] - dp_LCR[2];

        //�i�[
        DP_ixx[J + j] = ixx;
        DP_LCR[0][J + j] = double(dp_LCR[0]);//DP��double�ɂ���Ȃ炱����Ή������邱�ƁB
        DP_LCR[1][J + j] = double(dp_LCR[1]);//DP��double�ɂ���Ȃ炱����Ή������邱�ƁB
        DP_LCR[2][J + j] = double(dp_LCR[2]);//DP��double�ɂ���Ȃ炱����Ή������邱�ƁB
      }
    }
  }
  printf("*** calcDP_trapezoid *** :   end\n");
}

template<typename mTYPE1, typename mTYPE2>
void ML_EM(const DP_trapezoid &DP, const MATRIX<mTYPE2> &sinogram, MATRIX<mTYPE1> &image, double pix_size, int iteration)
//  sinogram�͍s�����e�����A�񂪃X�e�b�v���ƂȂ�悤�ɏ������邱�ƁB������W���̓��e�̘a(-log(I/I0))�ɂ��Ă������ƁB
//  pixel_size(cm)�B�����̒P�ʂ�cm�Ȃ̂ŁA�ŏI�I�ɍč\���摜�̒P�ʂ�(1/cm)�ɂȂ�B
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
//  sinogram�͍s�����e�����A�񂪃X�e�b�v���ƂȂ�悤�ɏ������邱�ƁB������W���̓��e�̘a(-log(I/I0))�ɂ��Ă������ƁB
//  pixel_size(cm)�B�����̒P�ʂ�cm�Ȃ̂ŁA�ŏI�I�ɍč\���摜�̒P�ʂ�(1/cm)�ɂȂ�B
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

  double *DP_LCR[3];//DP��double�ɂ���Ȃ炱����Ή������邱�ƁB
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

  double *DP_LCR[3];//DP��double�ɂ���Ȃ炱����Ή������邱�ƁB
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

  double *DP_LCR[3];//DP��double�ɂ���Ȃ炱����Ή������邱�ƁB
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
      y = y0 - i;//image���S�����_�Ƃ�����W
      XX = y*si + step0;
      for (j = 0; j < Width; j++) {
        if (img[I + j] < threshold) {
          img[I + j] = 0.0;
          continue;
        }//�����e����s�N�Z���l��臒l�ȉ��̂Ƃ��v�Z���ȗ�����

        x = x0 + j;//image���S�����_�Ƃ�����W
        xx = x*co + XX;//1�ڂ̃X�e�b�v�̒��S���[��������W�Ō����s�N�Z�����S�̓��e�ʒu
        ixx = int(floor(xx + 0.5));//�s�N�Z�����S�̓��e�ʒu���A�����ڂ̃X�e�b�v�̒��ɂ��邩���v�Z�B
        if (ixx<-1 || ixx>Step) continue;//�s�N�Z�����S�̓��e�ʒu�����e�̗̈�O�ł���ꍇ�͈ȍ~�̌v�Z���p�X����

//���o�m���v�Z
        dp_LCR[0] = dp_LCR[1] = dp_LCR[2] = 0;

        //�ЂƂO�̃X�e�b�v�ɓ��錟�o�m��(dp_LCR[0])
        ixx_side = ixx - 0.5;//�ЂƂO�̃X�e�b�v�Ƃ̋��E�̍��W
        d = ixx_side - (xx - 0.5*(a - b));//���E�Ƒ�`�̌��Ƃ̋���
        if (d >= 0) {//���E�̂ق����E
          dp_LCR[0] = (0.5*b + d)*c;
        }
        else if ((-d) < b) {//���E�̂ق����E������`�̑��̈ʒu���͉E
          dp_LCR[0] = 0.5*c*(b + d) / b*(b + d);
        }

        //�ЂƂ�̃X�e�b�v�ɓ��錟�o�m��(dp_LCR[1])
        ixx_side = ixx + 0.5;//�ЂƂ�̃X�e�b�v�Ƃ̋��E�̍��W
        d = (xx + 0.5*(a - b)) - ixx_side;//���E�Ƒ�`�̌��Ƃ̋���
        if (d >= 0) {//���E�̂ق�����
          dp_LCR[2] = (0.5*b + d)*c;
        }
        else if ((-d) < b) {
          dp_LCR[2] = 0.5*c*(b + d) / b*(b + d);
        }

        //���e�ʒu�̌��o��ɓ��錟�o�m��(dp_LCR[2])
        dp_LCR[1] = 1 - dp_LCR[0] - dp_LCR[2];

        //���o�m���v�Z�����܂�
        //���e�v�Z

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
          //���e�v�Z�����܂�
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
    IXX = k*Step;//sinogram�̃X�e�b�v�̃I�t�Z�b�g����

    a = DP.dp.out(k + 1, 1);
    b = DP.dp.out(k + 1, 2);
    c = DP.dp.out(k + 1, 3);
    si = DP.dp.out(k + 1, 4);
    co = DP.dp.out(k + 1, 5);
    for (i = 0; i < Height; i++) {
      I = i*Width;
      J = K + i*Width;

      y = y0 - i;//image���S�����_�Ƃ�����W
      XX = y*si + step0;
      for (j = 0; j < Width; j++) {
        x = x0 + j;//image���S�����_�Ƃ�����W
        xx = x*co + XX;//1�ڂ̃X�e�b�v�̒��S���[��������W�Ō����s�N�Z�����S�̓��e�ʒu
        ixx = int(floor(xx + 0.5));//�s�N�Z�����S�̓��e�ʒu���A�����ڂ̃X�e�b�v�̒��ɂ��邩���v�Z�B
        if (ixx<-1 || ixx>Step) continue;//�s�N�Z�����S�̓��e�ʒu�����e�̗̈�O�ł���ꍇ�͈ȍ~�̌v�Z���p�X����

 //���o�m���v�Z
        dp_LCR[0] = dp_LCR[1] = dp_LCR[2] = 0;

        //�ЂƂO�̃X�e�b�v�ɓ��錟�o�m��(dp_LCR[0])
        ixx_side = ixx - 0.5;//�ЂƂO�̃X�e�b�v�Ƃ̋��E�̍��W
        d = ixx_side - (xx - 0.5*(a - b));//���E�Ƒ�`�̌��Ƃ̋���
        if (d >= 0) {//���E�̂ق����E
          dp_LCR[0] = (0.5*b + d)*c;
        }
        else if ((-d) < b) {//���E�̂ق����E������`�̑��̈ʒu���͉E
          dp_LCR[0] = 0.5*c*(b + d) / b*(b + d);
        }

        //�ЂƂ�̃X�e�b�v�ɓ��錟�o�m��(dp_LCR[1])
        ixx_side = ixx + 0.5;//�ЂƂ�̃X�e�b�v�Ƃ̋��E�̍��W
        d = (xx + 0.5*(a - b)) - ixx_side;//���E�Ƒ�`�̌��Ƃ̋���
        if (d >= 0) {//���E�̂ق�����
          dp_LCR[2] = (0.5*b + d)*c;
        }
        else if ((-d) < b) {
          dp_LCR[2] = 0.5*c*(b + d) / b*(b + d);
        }

        //���e�ʒu�̌��o��ɓ��錟�o�m��(dp_LCR[2])
        dp_LCR[1] = 1 - dp_LCR[0] - dp_LCR[2];

//���o�m���v�Z�����܂�
//�t���e�v�Z
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
        //�t���e�v�Z�����܂�
      }
    }
  }

  for (i = 0; i < Num_pixel; i++) {
    img[i] /= Num_angle;
  }

  //  printf("*** BackProj_trape *** :   end\n");

}



/* �����𗧂֐��Q */
template<typename mTYPE>
MATRIX<mTYPE> click_angle(const MATRIX<mTYPE>& sinogram, int Num, int scan)//�T�C�m�O�����̓��e����͂���[�s�ԍ�]����n�܂�悤�ɂ��炷�B�s�����Ɋp�x������Ă��邱�ƁBscan=1�Ńn�[�t�A2�Ńt��
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
MATRIX<mTYPE> devide_each(const MATRIX<mTYPE>& sinogram)//���e�������Ƃ̍ő�l�Ŋe���e���������s���Ԃ��B
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