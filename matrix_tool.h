/*

�e���v���[�g�Ƃ��Ă�[short]�A[int]�A[float]�A[double]�ɑΉ��B
���Z+�A-�A*�A/�Aexp�Alog�Ȃǂ͂قƂ�ǔ{���x�v�Z����A���ʂ�double�^��MATRIX�ŕԂ�̂ŗp�ӂ����֐��ŃL���X�g���邱�ƁB
�o�C�i���[�t�@�C���J�̂��߂Ƀv���W�F�N�g�̃v���p�e�B�ŕ����Z�b�g��Unicode�����Z�b�g����}���`�o�C�g�����Z�b�g�ɕύX���邱�ƁB
hamaguchi

*/

#pragma once

#define _CRT_SECURE_NO_WARNINGS //�R���p�C������fopen�Ȃǂ̌x���𖳎����邽�߁B
#define NOMINMAX //windows.h��max�Amin�֐����`�����Ȃ����߁B

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string>
#include <windows.h>
using namespace std;


/*       MATRIX�N���X�ƃ����o�֐��Q�̐錾       */
template<typename mTYPE>
class MATRIX
{
private:
//�s��̏��
  int num_of_row;//�s���B
  int num_of_column;//�񐔁B
  mTYPE *value;//�s��v�f�̃|�C���^�B


public:
//�R���X�g���N�^�ƃf�X�g���N�^
  MATRIX();//�����Ȃ��̃R���X�g���N�^�B
  MATRIX(int,int);//�w�肵��[�s��]�A[��]�̍s���p�ӂ���B�l�͂��ׂĕs��ɂȂ�B
  MATRIX(const MATRIX&);//[�s��]�̃R�s�[�R���X�g���N�^�B�ϐ��̌^�������ł���΃R�s�[�ł���B
  ~MATRIX();//�f�X�g���N�^�B

//�A�N�Z�X�ᔽ�ɑϐ��̂���v���O�����ɗ��p����Ƃ悢�B�s��̂��߂̉��Z�֐���friend�ł͂Ȃ����̃|�C���^���g���Ă���B���ꂽ�l�����B
  mTYPE* pVALUE(int Row=1,int Column=1)const{return &value[(Row-1)*num_of_column+(Column-1)];};//private�ϐ���value�|�C���^�����o���Ē��ڑ���ł���B


//�s��̑S�ʂɊւ��鑀��
  void disp(char*)const;//�s����R���\�[���ɕ\������B���Ƃ���["%d"]�Ȃ�[�����w��]���Ă��Ƃ��̒ʂ�ɕ\������B\t�A\n�͂Ȃ��Ă悢�B
  void read(FILE*,int,int);//[�t�@�C���|�C���^]��n���Ďw�肵��[�s��]�A[��]�̍s��Ƃ݂Ȃ��ēǂݍ��ށB�e�L�X�g�t�@�C���̂ݑΉ��B
  void read(const char*,int,int);//[�t�@�C����]��n���Ďw�肵��[�s��]�A[��]�̍s��Ƃ݂Ȃ��ēǂݍ��ށB�o�C�i���[��ǂݍ��ނƂ��͏����͖��������B
  void write(FILE*,char*) const;//[�t�@�C���|�C���^]��n���āA���Ƃ���["%d"]�Ȃ�[�����w��]���Ă��ƍs����t�@�C���ɏ������ށB�e�L�X�g�t�@�C���̂ݑΉ��B
  void write(const char*,char*) const;//[�t�@�C����]���w�肵�āA���Ƃ���["%d"]�ȂǏ����w�肵�Ă��ƍs����t�@�C���ɏ������ށB�o�C�i���[���������ނƂ��͏����͖��������B
  void resize(int,int);//�w�肵��[�s��]�A[��]�̍s��ɕύX����B���l�͂��ׂĕs��ɂȂ�B
  void other_row_col(int,int);//�v�f�̐��͕ς��Ȃ����A�w�肵��[�s��]�A[��]�ɕύX����BN�~M�s���(N*M)�~1�s��ɂ���Ȃǂ̗p�r�B
  int Row() const;//�s����Ԃ��B
  int Col() const;//�񐔂�Ԃ��B
  int Data() const;//�s��v�f�̐���Ԃ��B
  void transpose();//�s���]�u����B
  void invert_row();//�s�����̕��т𔽓]����B
  void invert_col();//������̕��т𔽓]����B
  void invert_endian();//�G���f�B�A���ϊ����s���B�o�C�i���[����ǂݍ��񂾏ꍇ�̓G���f�B�A���ɒ��ӂ��邱�ƁB

//�s��v�f�̕ύX
  void in(int,int,mTYPE);//�w�肵��[�s�ԍ�]�A[��ԍ�]�̗v�f��^����[���l]�ŏ���������B�A�N�Z�X�ᔽ����t���B
  void constant(mTYPE);//���ׂĂ̗v�f���w�肵��[���l]�ɂ���B
  void unit();//�s���P�s��ɂ���B�����������s��̂݁B
  void put(int,int,const MATRIX&);//�w�肵��[�s�ԍ�]�A[��ԍ�]�̈ʒu������[�̊�Ƃ��āA�^����[�s��]�ŏ���������B
  void putrow(int,const MATRIX&);//�w�肵��[�s�ԍ�]�̍s��^����[�s��]�ŏ���������B�������s�x�N�g���Ɍ���B//20161227�Ȃ������܂������Ȃ����Ƃ�����
  void putcol(int,const MATRIX&);//�w�肵��[��ԍ�]�̗��^����[�s��]�ŏ���������B��������x�N�g���Ɍ���B//20161227�Ȃ������܂������Ȃ����Ƃ�����
  void addrow(int,int Height=1,mTYPE Value=0);//�w�肵��[�s�ԍ�]�̍s�̎��Ɏw�肵��[����=�s��]�Ŏw�肵��[���l]�̍s��ǉ�����B[�s�ԍ�]�[���Ő擪�ɒǉ��B
  void addrow(int,const MATRIX&);//�w�肵��[�s�ԍ�]�̍s�̎��̍s����^����[�s��]��ǉ�����B[�s�ԍ�]�[���Ő擪����ǉ��B
  void addcol(int,int Width=1,mTYPE Value=0);////�w�肵��[��ԍ�]�̗�̎��Ɏw�肵��[��=��]�Ŏw�肵��[���l]�̍s��ǉ�����B[��ԍ�]�[���Ő擪�ɒǉ��B
  void addcol(int,const MATRIX&);//�w�肵��[��ԍ�]�̗�̎��̍s����^����[�s��]��ǉ�����B[��ԍ�]�[���Ő擪����ǉ��B
  void delrow(int,int Height=1);//�w�肵��[�s�ԍ�]�̍s����w�肵��[����=�s��]�̍s���폜���ċl�߂�B
  void delcol(int,int Width=1);//�w�肵��[��ԍ�]�̗񂩂�w�肵��[��=��]�̗���폜���ċl�߂�B

//�s��̏o��
  mTYPE out(int,int) const;//�w�肵��[�s�ԍ�]�A[��ԍ�]�̈ʒu�̒l��Ԃ��B�A�N�Z�X�ᔽ����t���B
  MATRIX pull(int,int,int,int) const;//�w�肵��[�s�ԍ�]�A[��ԍ�]�̈ʒu������[�̊�Ƃ��āA�w�肵��[����=�s��]�A[��=��]�͈̔͂̍s���Ԃ��B
  MATRIX pullrow(int)const;//�w�肵��[�s�ԍ�]�̍s��Ԃ��B
  MATRIX pullcol(int)const;//�w�肵��[��ԍ�]�̗��Ԃ��B
  mTYPE max()const;//�v�f�̍ő�l��Ԃ��B 
  mTYPE min()const;//�v�f�̍ŏ��l��Ԃ��B 
  double average()const;//�s��̑S�v�f�̕��ς�Ԃ��B//�傫�ȉ摜�S�̂ɑ΂��Ă̓���͂��s��������
  double stddev()const;//�s��̑S�v�f�̕W���΍���Ԃ��B//�傫�ȉ摜�S�̂ɑ΂��Ă̓���͂��s��������
//  MATRIX abs(const MATRIX&);//�s��̊e�v�f�̐�Βl��������s���Ԃ��B������2016/11/09hamaguchi

//�L���X�g����MATRIX��Ԃ����߂̊֐��@���ɒu���^�C�v�̂��́B
  MATRIX<short> m_short()const;//�s��̕ϐ��^��short�ɕύX�����s���Ԃ��B
  MATRIX<int> m_int()const;//�s��̕ϐ��^��int�ɕύX�����s���Ԃ��B
  MATRIX<float> m_float()const;//�s��̕ϐ��^��float�ɕύX�����s���Ԃ��B 
  MATRIX<double> m_double()const;//�s��̕ϐ��^��double�ɕύX�����s���Ԃ��B 

//������Z�q�I�[�o�[���[�h
  MATRIX& operator = (const MATRIX&);//���ӂ̍s��ɉE�ӂ̍s���������B�Öق̌^�ϊ��͂����Ȃ����ߍs��̕ϐ��^�������ꍇ�Ɍ���B�L���X�g���g�����ƁB
};



/*       MATRIX�N���X�̂��߂̊֐��Q�̃v���g�^�C�v�錾(��`�͉��̕���)       */

//�L���X�g����MATRIX��Ԃ����߂̊֐��@�O�ɒu���^�C�v�̂���
template<typename mTYPE>
MATRIX<short> m_short(const MATRIX<mTYPE>&);//[�s��]�̕ϐ��^��short�ɕύX�����s���Ԃ��B

template<typename mTYPE>
MATRIX<int> m_int(const MATRIX<mTYPE>&);//[�s��]�̕ϐ��^��int�ɕύX�����s���Ԃ��B

template<typename mTYPE>
MATRIX<float> m_float(const MATRIX<mTYPE>&);//[�s��]�̕ϐ��^��float�ɕύX�����s���Ԃ��B

template<typename mTYPE>
MATRIX<double> m_double(const MATRIX<mTYPE>&);//[�s��]�̕ϐ��^��double�ɕύX�����s���Ԃ��B

template<typename type,typename mTYPE>
MATRIX<type> m_cast(const MATRIX<mTYPE>&);//[�s��]�̕ϐ��^���w�肵��[�^]�ɕύX�����s���Ԃ��Btype�̌���͔͂���ł��Ȃ��̂Œ��ӁB


//�s��̂��߂̉��Z�֐�
template<typename mTYPE1,typename mTYPE2>
MATRIX<double> Multiplication(const MATRIX<mTYPE1>&,const MATRIX<mTYPE2>&);//[�s��]�A[�s��]�̂����Z���s���A���ʂ̍s���double�^�ŕԂ��B

template<typename mTYPE>
MATRIX<double> exp(const MATRIX<mTYPE>&);//[�s��]�̊e�v�f��exp��������s���double�^�ŕԂ��B

template<typename mTYPE>
MATRIX<double> log(const MATRIX<mTYPE>&);//[�s��]�̊e�v�f��log��������s���double�^�ŕԂ��B


//���Z�q�I�[�o�[���[�h
template<typename mTYPE>
MATRIX<mTYPE> operator -(const MATRIX<mTYPE>&);//[�s��]�̊e�v�f�̕����𔽓]�����s���^����ꂽ�s��̌^�ŕԂ��B//�߂�l��mTYPE�̗��֐��ɂ��ėv����

template<typename mTYPE>
MATRIX<double> operator +(double,const MATRIX<mTYPE>&);//[���l]��[�s��]�̊e�v�f�ɑ������s���double�^�ŕԂ��B//�߂�l��double�̗��֐��ɂ��ėv����

template<typename mTYPE>
MATRIX<double> operator +(const MATRIX<mTYPE>&,double);//[�s��]�̊e�v�f��[���l]�𑫂����s���double�^�ŕԂ��B//�߂�l��double�̗��֐��ɂ��ėv����

template<typename mTYPE1,typename mTYPE2>
MATRIX<double> operator +(const MATRIX<mTYPE1>&,const MATRIX<mTYPE2>&);//[�s��]��[�s��]�̊e�v�f�𑫂����s���double�^�ŕԂ��B//�߂�l��mTYPE�̗��֐��ɂ��ėv����

template<typename mTYPE>
MATRIX<double> operator -(double,const MATRIX<mTYPE>&);//[���l]����[�s��]�̊e�v�f���������s���double�^�ŕԂ��B//�߂�l��double�̗��֐��ɂ��ėv����

template<typename mTYPE>
MATRIX<double> operator -(const MATRIX<mTYPE>&,double);//[�s��]�̊e�v�f����[���l]���������s���double�^�ŕԂ��B//�߂�l��double�̗��֐��ɂ��ėv����

template<typename mTYPE1,typename mTYPE2>
MATRIX<double> operator -(const MATRIX<mTYPE1>&,const MATRIX<mTYPE2>&);//[�s��]����[�s��]�̊e�v�f�����������s���double�^�ŕԂ��B//�߂�l��double�̗��֐��ɂ��ėv����

template<typename mTYPE>
MATRIX<double> operator *(double,const MATRIX<mTYPE>&);//[���l]��[�s��]�̊e�v�f�Ɋ|�����s���double�^�ŕԂ��B//�߂�l��double�̗��֐��ɂ��ėv����

template<typename mTYPE>
MATRIX<double> operator *(const MATRIX<mTYPE>&,double);//[�s��]�̊e�v�f��[���l]���|�����s���double�^�ŕԂ��B//�߂�l��double�̗��֐��ɂ��ėv����

template<typename mTYPE1,typename mTYPE2>
MATRIX<double> operator *(const MATRIX<mTYPE1>&,const MATRIX<mTYPE2>&);//[�s��]��[�s��]�̊e�v�f���|�����s���double�^�ŕԂ��B//�߂�l��double�̗��֐��ɂ��ėv����

template<typename mTYPE>
MATRIX<double> operator /(double,const MATRIX<double>&);//[���l]��[�s��]�̊e�v�f�Ŋ������s���double�^�ŕԂ��B//�߂�l��double�̗��֐��ɂ��ėv����

template<typename mTYPE>
MATRIX<double> operator /(const MATRIX<double>&,double);//[�s��]�̊e�v�f��[���l]�Ŋ������s���double�^�ŕԂ��B//�߂�l��double�̗��֐��ɂ��ėv����

template<typename mTYPE1,typename mTYPE2>
MATRIX<double> operator /(const MATRIX<mTYPE1>&,const MATRIX<mTYPE2>&);//[���l]��[�s��]�̊e�v�f�Ŋ������s���double�^�ŕԂ��B//�߂�l��mTYPE�̗��֐��ɂ��ėv����



/*       MATRIX�N���X�̃����o�֐��Q�̒�`       */

//�R���X�g���N�^�ƃf�X�g���N�^
MATRIX<short>::MATRIX()
{
  num_of_row=1;
  num_of_column=1;
 
  value=(short *)malloc(sizeof(short)*1);
  if(value==NULL)
  {
    printf("MEMORY ALLOCATION ERROR.\n");
    system("pause");
  }
  
}

MATRIX<int>::MATRIX()
{
  num_of_row=1;
  num_of_column=1;
 
  value=(int *)malloc(sizeof(int)*1);
  if(value==NULL)
  {
    printf("MEMORY ALLOCATION ERROR.\n");
    system("pause");
  }
  
}

MATRIX<float>::MATRIX()
{
  num_of_row=1;
  num_of_column=1;
 
  value=(float *)malloc(sizeof(float)*1);
  if(value==NULL)
  {
    printf("MEMORY ALLOCATION ERROR.\n");
    system("pause");
  }  
}

MATRIX<double>::MATRIX()
{
  num_of_row=1;
  num_of_column=1;
 
  value=(double *)malloc(sizeof(double)*1);
  if(value==NULL)
  {
    printf("MEMORY ALLOCATION ERROR.\n");
    system("pause");
  }
  
}

MATRIX<short>::MATRIX(int Row,int Column)
{
  int buf_row=Row;
  int buf_column=Column;

  if(buf_row<1 || buf_column<1)
  {
    printf("Size error.%d*%d_matrix",buf_row,buf_column);
    system("pause");
    exit(EXIT_FAILURE);
  }

  int buf_data=Row*Column;

  num_of_row=buf_row;
  num_of_column=buf_column;


  value=(short *)malloc(sizeof(short)*buf_data);
  if(value==NULL)
  {
    printf("MEMORY ALLOCATION ERROR.\n");
    system("pause");
    exit(EXIT_FAILURE);
  }
}

MATRIX<int>::MATRIX(int Row,int Column)
{
  int buf_row=Row;
  int buf_column=Column;

  if(buf_row<1 || buf_column<1)
  {
    printf("Size error.%d*%d_matrix",buf_row,buf_column);
    system("pause");
    exit(EXIT_FAILURE);
  }

  int buf_data=Row*Column;

  num_of_row=buf_row;
  num_of_column=buf_column;


  value=(int *)malloc(sizeof(int)*buf_data);
  if(value==NULL)
  {
    printf("MEMORY ALLOCATION ERROR.\n");
    system("pause");
    exit(EXIT_FAILURE);
  }
}

MATRIX<float>::MATRIX(int Row,int Column)
{
  int buf_row=Row;
  int buf_column=Column;

  if(buf_row<1 || buf_column<1)
  {
    printf("Size error.%d*%d_matrix",buf_row,buf_column);
    system("pause");
    exit(EXIT_FAILURE);
  }

  int buf_data=Row*Column;

  num_of_row=buf_row;
  num_of_column=buf_column;


  value=(float *)malloc(sizeof(float)*buf_data);
  if(value==NULL)
  {
    printf("MEMORY ALLOCATION ERROR.\n");
    system("pause");
    exit(EXIT_FAILURE);
  }
}

MATRIX<double>::MATRIX(int Row,int Column)
{
  int buf_row=Row;
  int buf_column=Column;

  if(buf_row<1 || buf_column<1)
  {
    printf("Size error.%d*%d_matrix",buf_row,buf_column);
    system("pause");
    exit(EXIT_FAILURE);
  }

  int buf_data=Row*Column;

  num_of_row=buf_row;
  num_of_column=buf_column;


  value=(double *)malloc(sizeof(double)*buf_data);
  if(value==NULL)
  {
    printf("MEMORY ALLOCATION ERROR.\n");
    system("pause");
    exit(EXIT_FAILURE);
  }
}

template<typename mTYPE>
MATRIX<mTYPE>::MATRIX(const MATRIX<mTYPE> &matrix)
{
  int buf_row=matrix.num_of_row;
  int buf_column=matrix.num_of_column;
  int buf_data=buf_row*buf_column;

  num_of_row=buf_row;
  num_of_column=buf_column;
  value=(mTYPE *)malloc(sizeof(mTYPE)*buf_data);
  if(value==NULL)
  {
    printf("MEMORY ALLOCATION ERROR.\n");
    system("pause");
  }
  
  int i;
  mTYPE *buf_m_value=matrix.value;
  mTYPE *buf_value=value;

  for(i=0;i<buf_data;i++)
  {
    buf_value[i]=buf_m_value[i];
  }  
  
}

template<typename mTYPE>
MATRIX<mTYPE>::~MATRIX()
{
  free(value);
}


//�s��̑S�ʂɊւ��鑀��
template<typename mTYPE>
void MATRIX<mTYPE>::disp(char* format)const//�s����R���\�[���ɕ\������B
{
  int buf_row=num_of_row;
  int buf_column=num_of_column;

  mTYPE *buf_value=value;

  char *tab="\t";
  string TAB=tab;

  string FORMAT=format;
  FORMAT=TAB+FORMAT;

  int i,j;
  for(i=0;i<buf_row;i++)
  {
    for(j=0;j<buf_column;j++)
    {
      printf(FORMAT.c_str(),buf_value[i*buf_column+j]);
    }
    printf("\n");
  }
  printf("\n");
}

template<typename mTYPE>
void MATRIX<mTYPE>::read(FILE *fp,int Row,int Column)
{
  if(fp==NULL)
  {
    printf("file cannot open. read(*fp,(%d,%d))x\n",Row,Column);
    system("pause");
    exit(EXIT_FAILURE);
  }

  resize(Row,Column);

  int i;

  double buf_Value;
  int buf_data=Row*Column;
  mTYPE *buf_value=value;

  for(i=0;i<buf_data;i++)
  {
    fscanf(fp,"%lf",&buf_Value);
    buf_value[i]=mTYPE(buf_Value);        
  }  
}

template<typename mTYPE>
void MATRIX<mTYPE>::read(const char *File_name,int Row,int Column)
{
  resize(Row,Column);

  int i;
  double buf_Value;
  int buf_data=Row*Column;
  mTYPE *buf_value=value;


//.txt�g���q�̔��肱������
  char c_name_end[5];
  string s_name_end;

  bool txt=false;
  int size=0;
  while(File_name[size]!=0x00) size++;
  if(size>=5)
  {
    for(i=0;i<5;i++) c_name_end[i]=File_name[size-4+i];
    s_name_end=c_name_end;
    if(s_name_end==".txt") txt=true; 
  }
//.txt�g���q�̔��肨���

  if(txt)
  {  
    FILE *fp;
    fp=fopen(File_name,"r");
    if(fp==NULL)
    {
      printf("file cannot open. read(%s,(%d,%d))\n",File_name,Row,Column);
      system("pause");
      exit(EXIT_FAILURE);     
   }

    for(i=0;i<buf_data;i++)
    {
      if(fscanf(fp,"%lf",&buf_Value)==EOF) break;
      buf_value[i]=mTYPE(buf_Value);        
    }  

    fclose(fp);
  }
  else
  {
    HANDLE hFile=CreateFile(File_name,GENERIC_READ,0,NULL,OPEN_EXISTING,FILE_ATTRIBUTE_NORMAL,NULL);
    if(hFile==INVALID_HANDLE_VALUE)
    {
      printf("file cannot open. read(%s,(%d,%d))\n",File_name,Row,Column);
      system("pause");
      exit(EXIT_FAILURE);     
    }

    int data_size=buf_data*sizeof(mTYPE);
    unsigned long BytesDone;

    ReadFile(hFile,buf_value,data_size,&BytesDone,NULL);

    CloseHandle(hFile);   
  }
}

template<typename mTYPE>
void MATRIX<mTYPE>::write(FILE *fp,char *format) const
{
  if(fp==NULL)
  {
    printf("file cannot open. write(*fp)\n");
    system("pause");
    exit(EXIT_FAILURE);
  }

  int i,j;
  int buf_row=num_of_row;
  int buf_column=num_of_column;
  mTYPE *buf_value=value;

  char *tab="\t";
  string TAB=tab;

  string FORMAT=format;
  FORMAT=FORMAT+TAB;

  for(i=0;i<buf_row;i++)
  {
    for(j=0;j<buf_column;j++)
    {
      fprintf(fp,FORMAT.c_str(),buf_value[i*buf_column+j]);
    }
    fprintf(fp,"\n");
  }
}

template<typename mTYPE>
void MATRIX<mTYPE>::write(const char *File_name,char *format) const
{

  int i,j;
  int buf_row=num_of_row;
  int buf_column=num_of_column;
  int buf_data=buf_row*buf_column;
  mTYPE *buf_value=value;


//.txt�g���q�̔��肱������
  char c_name_end[5];
  string s_name_end;

  bool txt=false;
  int size=0;
  while(File_name[size]!=0x00) size++;
  if(size>=5)
  {
    for(i=0;i<5;i++) c_name_end[i]=File_name[size-4+i];
    s_name_end=c_name_end;
    if(s_name_end==".txt") txt=true; 
  }
//.txt�g���q�̔��肨���

  if(txt)
  {
    FILE *fp;
    fp=fopen(File_name,"w");
    if(fp==NULL)
    {
      printf("file cannot open. write(%s)\n",File_name);
      system("pause");
      exit(EXIT_FAILURE);     
    }

    char *tab="\t";
    string TAB=tab;

    string FORMAT=format;
    FORMAT=FORMAT+TAB;

    for(i=0;i<buf_row;i++)
    {
      for(j=0;j<buf_column;j++)
      {
        fprintf(fp,FORMAT.c_str(),buf_value[i*buf_column+j]);
      }
      fprintf(fp,"\n");
    }

    fclose(fp);
  }
  else
  {
    int data_size=buf_data*sizeof(mTYPE);
    unsigned long BytesDone;

    HANDLE hFile=CreateFile(File_name,GENERIC_WRITE,0,NULL,CREATE_ALWAYS,FILE_ATTRIBUTE_NORMAL,NULL);
    if(hFile!=INVALID_HANDLE_VALUE)
    {
      WriteFile(hFile, buf_value,data_size, &BytesDone, NULL);
      CloseHandle(hFile);  
    }    
  }

}

template<typename mTYPE>
void MATRIX<mTYPE>::resize(int Row,int Column)//�s��̃T�C�Y��ύX����B���l�͂��ׂĕs��ɂȂ�B
{
  int buf_row=Row;
  int buf_column=Column;

  if(buf_row<1 || buf_column<1)
  {
    printf("Size error.(%d,%d)\n",buf_row,buf_column);
    system("pause");
    exit(EXIT_FAILURE);
  }

  free(value);

  num_of_row=buf_row;
  num_of_column=buf_column;
  
  int buf_data=buf_row*buf_column;
  value=(mTYPE *)malloc(sizeof(mTYPE)*buf_data);
  if(value==NULL)
  {
    printf("MEMORY ALLOCATION ERROR.\n");
    system("pause");
    exit(EXIT_FAILURE);
  }
}

template<typename mTYPE>
void MATRIX<mTYPE>::other_row_col(int Row,int Column)
{
  if(Row*Column!=num_of_row*num_of_column)
  {
    printf("Size error.other_row_col(%d,%d)\n",Row,Column);
    system("pause");
    exit(EXIT_FAILURE);
  }
  num_of_row=Row;
  num_of_column=Column;
}

template<typename mTYPE>
int MATRIX<mTYPE>::Row()const
{
  return num_of_row;
}

template<typename mTYPE>
int MATRIX<mTYPE>::Col()const
{
  return num_of_column;
}

template<typename mTYPE>
int MATRIX<mTYPE>::Data()const
{
  return num_of_row*num_of_column;
}

template<typename mTYPE>
void MATRIX<mTYPE>::transpose()//�s���]�u����B
{
  int i,j;
  mTYPE *buf_value;
  int buf_row=num_of_row;
  int buf_column=num_of_column;
  int buf_data=buf_row*buf_column;

  buf_value=(mTYPE *)malloc(sizeof(mTYPE)*buf_data);
  if(buf_value==NULL)
  {
    printf("MEMORY ALLOCATION ERROR.transpose()\n");
    system("pause");
    exit(EXIT_FAILURE);
  }
  for(i=0;i<buf_row;i++)//�]�u��̔z��ɕ��בւ���
  {
    for(j=0;j<buf_column;j++)
    {
      buf_value[j*buf_row+i]=value[i*buf_column+j];
    }
  }

  //���בւ������̂�V�����s��v�f�Ƃ���
  num_of_row=buf_column;
  num_of_column=buf_row;
  free(value);
  value=buf_value;
}

template<typename mTYPE>
void MATRIX<mTYPE>::invert_row()
{
  int i,j;
  int buf_row=num_of_row;
  int buf_column=num_of_column;
  mTYPE *buf_value=value;

  mTYPE *buf_invert_value;
  buf_invert_value=(mTYPE *)malloc(sizeof(mTYPE)*buf_row*buf_column);
  
  for(i=0;i<buf_row;i++)
  {
    for(j=0;j<num_of_column;j++)
    {
      buf_invert_value[i*buf_column+(buf_column-1-j)]=buf_value[i*buf_column+j];
    }
  }

  //���בւ������̂�V�����s��v�f�Ƃ���
  free(value);
  value=buf_invert_value;
}

template<typename mTYPE>
void MATRIX<mTYPE>::invert_col()
{
  int i,j;
  int buf_row=num_of_row;
  int buf_column=num_of_column;
  mTYPE *buf_value=value;

  mTYPE *buf_invert_value;
  buf_invert_value=(mTYPE *)malloc(sizeof(mTYPE)*buf_row*buf_column);
  
  for(i=0;i<buf_row;i++)
  {
    for(j=0;j<num_of_column;j++)
    {
      buf_invert_value[(buf_row-1-i)*buf_column+j]=buf_value[i*buf_column+j];
    }
  }

  //���בւ������̂�V�����s��v�f�Ƃ���
  free(value);
  value=buf_invert_value;
}

template<typename mTYPE>
void MATRIX<mTYPE>::invert_endian()
{
  char *p;
  char *q;
  int buf_data=num_of_row*num_of_column;
  int data_size=sizeof(mTYPE);
  mTYPE Value;
  mTYPE *buf_value=value;

  int i,j;
  for(i=0;i<buf_data;i++)
  {
    Value=buf_value[i];
    p=(char*)&Value;
    q=(char*)&buf_value[i];
    for(j=0;j<data_size;j++)
    {
      *(q+j)=*(p+data_size-1-j);
    }    
  }
}


//�s��v�f�̕ύX
template<typename mTYPE>
void MATRIX<mTYPE>::in(int Row,int Column,mTYPE Value)//�w�肵���s�A��̗v�f��^�����l�ŏ���������B
{

  if(Row<1 || Column<1 || Row>num_of_row || Column>num_of_column)
  {
    printf("Cannot access.in(%d,%d)\n",Row,Column);
    system("pause");
    exit(EXIT_FAILURE);
  }

  int i;
  i=(Row-1)*num_of_column+(Column-1);

  mTYPE *buf_value=value;

  buf_value[i]=mTYPE(Value); 
}

template<typename mTYPE>
void MATRIX<mTYPE>::constant(mTYPE Value)//���ׂĂ̗v�f���w�肵�����l�ɂ���
{
  int i;
  int buf_data=num_of_row*num_of_column;

  mTYPE *buf_value=value;

  for(i=0;i<buf_data;i++)
  {
    buf_value[i]=mTYPE(Value);
  }
}

template<typename mTYPE>
void MATRIX<mTYPE>::unit()//�s���P�s��ɂ���B�����������s��̂݁B
{
  int buf_row=num_of_row;
  int buf_column=num_of_column;

  if(buf_row!=buf_column)
  {
    printf("Not square matrix.unit()\n");
    system("pause");
    exit(EXIT_FAILURE);
  }

  mTYPE *buf_value=value;

  int i,j;
  for(i=0;i<buf_row;i++)
  {
    for(j=0;j<buf_column;j++)
    {
      buf_value[i*buf_column+j]=mTYPE(0);
    }
  }

  for(i=0;i<buf_row;i++)
  {
    buf_value[i*buf_column+i]=mTYPE(1);
  }

}

template<typename mTYPE>
void MATRIX<mTYPE>::put(int Row,int Column,const MATRIX<mTYPE> &matrix)//�w�肵���s�A��̗v�f����Ƃ��āA�^�����s��ŏ���������B
{
  int buf_row=num_of_row;
  int buf_column=num_of_column;
  int buf_m_row=matrix.num_of_row;
  int buf_m_column=matrix.num_of_column;

  if(Row<1 || Column<1 || Row>buf_row || Column>buf_column || (Row+buf_m_row-1)>buf_row || (Column+buf_m_column-1)>buf_column)
  {
    printf("Size error.put(%d,%d,%d,%d)\n",Row,Column,buf_m_row,buf_m_column);
    system("pause");
    exit(EXIT_FAILURE);
  }

  int i,j;
  mTYPE *buf_m_value=matrix.value;
  mTYPE *buf_value=value;

  for(i=0;i<buf_m_row;i++)
  {
    for(j=0;j<buf_m_column;j++)
    {
      buf_value[(i+Row-1)*buf_column+(j+Column-1)]=buf_m_value[i*buf_m_column+j];
    }
  }
}

template<typename mTYPE>
void MATRIX<mTYPE>::putrow(int Row,const MATRIX<mTYPE> &matrix)//�w�肵���s��^�����s��ŏ���������B�������s�x�N�g���Ɍ���B
{
  int buf_row=num_of_row;
  int buf_column=num_of_column;
  int buf_m_row=matrix.num_of_row;
  int buf_m_column=matrix.num_of_column;

  if(Row<1 || Row>buf_row || buf_m_column!=buf_column || buf_m_row!=1)
  {
    printf("Size error.putrow(%d,%d,%d)\n",Row,buf_m_row,buf_m_column);
    system("pause");
    exit(EXIT_FAILURE);
  }

  int i;
  mTYPE *buf_m_value=matrix.value;
  mTYPE *buf_value=value;

  for(i=0;i<buf_m_column;i++)
  {
    buf_value[(Row-1)*buf_column+i]=buf_m_value[i];
  }
}

template<typename mTYPE>
void MATRIX<mTYPE>::putcol(int Column,const MATRIX<mTYPE> &matrix)//�w�肵�����^�����s��ŏ���������B��������x�N�g���Ɍ���B
{
  int buf_row=num_of_row;
  int buf_column=num_of_column;
  int buf_m_row=matrix.num_of_row;
  int buf_m_column=matrix.num_of_column;

  if(Column<1 || Column>buf_column || buf_m_row!=buf_row || buf_m_column!=1)
  {
    printf("Size error.putcol(%d,%d,%d)\n",Column,buf_m_row,buf_m_column);
    system("pause");
    exit(EXIT_FAILURE);
  }

  int i;
  mTYPE *buf_m_value=matrix.value;
  mTYPE *buf_value=value;

  for(i=0;i<buf_m_row;i++)
  {
    buf_value[i*buf_column+(Column-1)]=buf_m_value[i];
  }
}

template<typename mTYPE>
void MATRIX<mTYPE>::addrow(int Row,int Height,mTYPE Value)//�w�肵���s�̉���Value�̍s���w�肵�����Œǉ�����B
{
  int buf_row=num_of_row;
  int buf_column=num_of_column;

  if(Row<0 || Row>buf_row)
  {
    printf("Cannot access.addrow(%d)\n",Row);
    system("pause");
    exit(EXIT_FAILURE);
  }
  
  int i,j,buf;
  mTYPE *buf_a_value;
  buf_a_value=(mTYPE *)malloc(sizeof(mTYPE)*(buf_row+Height)*buf_column);
  if(buf_a_value==NULL)
  {
    printf("MEMORY ALLOCATION ERROR.addrow(%d,%d)\n",Row,Height);
    system("pause");
    exit(EXIT_FAILURE);
  }
  
  mTYPE *buf_value=value;

  for(i=0;i<Row;i++)
  {
    for(j=0;j<buf_column;j++)
    {
      buf_a_value[i*buf_column+j]=buf_value[i*buf_column+j];
    }
  }
  buf=Row+Height;
  for(i=Row;i<buf;i++)
  {
    for(j=0;j<buf_column;j++)
    {
      buf_a_value[i*buf_column+j]=mTYPE(Value);
    }  
  }
  buf=buf_row+Height;
  for(i=(Row+Height);i<buf;i++)
  {
    for(j=0;j<buf_column;j++)
    {
      buf_a_value[i*buf_column+j]=buf_value[(i-Height)*buf_column+j];
    }  
  }

  num_of_row=buf_row+Height;
  free(value);
  value=buf_a_value;
}

template<typename mTYPE>
void MATRIX<mTYPE>::addrow(int Row,const MATRIX<mTYPE> &matrix)
{
  if(num_of_column!=matrix.Col())
  {
    printf("Size error.addrow(%d,(%d,%d))\n",Row,matrix.Row(),matrix.Col());
    system("pause");
    exit(EXIT_FAILURE);
  }
  addrow(Row,matrix.Row());
  put(Row+1,1,matrix);
}

template<typename mTYPE>
void MATRIX<mTYPE>::addcol(int Column,int Width,mTYPE Value)//�w�肵����̉��ɒlValue�̗���w�肵�����Œǉ�����B
{
  int buf_row=num_of_row;
  int buf_column=num_of_column;

  if(Column<0 || Column>buf_column)
  {
    printf("Cannot access.addcol(%d)\n",Column);
    system("pause");
    exit(EXIT_FAILURE);
  }
  
  int i,j,buf;
  mTYPE *buf_a_value;
  buf_a_value=(mTYPE *)malloc(sizeof(mTYPE)*buf_row*(buf_column+Width));
  if(buf_a_value==NULL)
  {
    printf("MEMORY ALLOCATION ERROR.addcol(%d,%d)\n",Column,Width);
    system("pause");
    exit(EXIT_FAILURE);
  }
  
  transpose();//�]�u����addrow�𗘗p

  int Row=Column;  
  int Height=Width;
  buf=buf_row;//�]�u�ɍ��킹�čs���Ɨ񐔂�����
  buf_row=buf_column;
  buf_column=buf;

  mTYPE *buf_value=value;

  for(i=0;i<Row;i++)
  {
    for(j=0;j<buf_column;j++)
    {
      buf_a_value[i*buf_column+j]=buf_value[i*buf_column+j];
    }
  }
  buf=Row+Height;
  for(i=Row;i<buf;i++)
  {
    for(j=0;j<buf_column;j++)
    {
      buf_a_value[i*buf_column+j]=mTYPE(Value);
    }  
  }
  buf=buf_row+Height;
  for(i=(Row+Height);i<buf;i++)
  {
    for(j=0;j<buf_column;j++)
    {
      buf_a_value[i*buf_column+j]=buf_value[(i-Height)*buf_column+j];
    }  
  }

  num_of_row=buf_row+Height;
  free(value);
  value=buf_a_value;
      
  transpose();//�]�u���Č��ɖ߂�
}

template<typename mTYPE>
void MATRIX<mTYPE>::addcol(int Column, const MATRIX<mTYPE> &matrix)
{
  if(num_of_row!=matrix.Row())
  {
    printf("Size error.addcol(%d,(%d,%d))\n",Column,matrix.Row(),matrix.Col());
    system("pause");
    exit(EXIT_FAILURE);
  }

  addcol(Column,matrix.Col());
  put(1,Column+1,matrix);
}

template<typename mTYPE>
void MATRIX<mTYPE>::delrow(int Row,int Height)
{
  int buf_row=num_of_row;
  int buf_column=num_of_column;

  if(Row<1 || Row>buf_row || (Row+Height-1)>buf_row || Height>=buf_row)
  {
    printf("Cannot access.delrow(%d,%d)\n",Row,Height);
    system("pause");
    exit(EXIT_FAILURE);
  }

  int i,j,buf;
  mTYPE *buf_d_value;
  
  buf_d_value=(mTYPE *)malloc(sizeof(mTYPE)*(buf_row-Height)*buf_column);
  if(buf_d_value==NULL)
  {
    printf("MEMORY ALLOCATION ERROR.delrow(%d,%d)\n",Row,Height);
    system("pause");
    exit(EXIT_FAILURE);
  }
  
  mTYPE *buf_value=value;

  buf=Row-1;
  for(i=0;i<buf;i++)
  {
    for(j=0;j<buf_column;j++)
    {
      buf_d_value[i*buf_column+j]=buf_value[i*buf_column+j];
    }
  }
  buf=buf_row-Height;
  for(i=Row-1;i<buf;i++)
  {
    for(j=0;j<buf_column;j++)
    {
      buf_d_value[i*buf_column+j]=buf_value[(i+Height)*buf_column+j];
    }  
  }

  num_of_row=buf_row-Height;
  free(value);
  value=buf_d_value; 
 
}

template<typename mTYPE>
void MATRIX<mTYPE>::delcol(int Column,int Width)
{
  int buf_row=num_of_row;
  int buf_column=num_of_column;

  if(Column<1 || Column>buf_column || (Column+Width-1)>buf_column || Width>=buf_column)
  {
    printf("Cannot access.delcol(%d,%d)\n",Column,Width);
    system("pause");
    exit(EXIT_FAILURE);
  }

  int i,j,buf;
  mTYPE *buf_b_value;
  buf_b_value=(mTYPE *)malloc(sizeof(mTYPE)*buf_row*(buf_column-Width));
  if(buf_b_value==NULL)
  {
    printf("MEMORY ALLOCATION ERROR.delcol(%d,%d)\n",Column,Width);
    system("pause");
    exit(EXIT_FAILURE);
  }
  
  transpose();//�]�u����delrow�𗘗p
  int Row=Column;  
  int Height=Width;
  buf=buf_row;//�]�u�ɍ��킹�čs���Ɨ񐔂�����
  buf_row=buf_column;
  buf_column=buf;

  mTYPE *buf_value=value;

  buf=Row-1;
  for(i=0;i<buf;i++)
  {
    for(j=0;j<buf_column;j++)
    {
      buf_b_value[i*buf_column+j]=buf_value[i*buf_column+j];
    }
  }
  buf=buf_row-Height;
  for(i=Row-1;i<buf;i++)
  {
    for(j=0;j<buf_column;j++)
    {
      buf_b_value[i*buf_column+j]=buf_value[(i+Height)*buf_column+j];
    }  
  }

  num_of_row=buf_row-Height;
  free(value);
  value=buf_b_value;

  transpose();  
}


//�s��v�f�̏o��
template<typename mTYPE>
mTYPE MATRIX<mTYPE>::out(int Row,int Column)const//�w�肵���s�A��̗v�f��Ԃ��B
{

  if(Row<1 || Column<1 || Row>num_of_row || Column>num_of_column)
  {
    printf("Cannot access.out(%d,%d)\n",Row,Column);
    system("pause");
    exit(EXIT_FAILURE);
  }

  int i;
  i=num_of_column*(Row-1)+(Column-1);

  mTYPE *buf_value=value;

  return buf_value[i];
}

template<typename mTYPE>
MATRIX<mTYPE> MATRIX<mTYPE>::pull(int Row,int Column,int Height,int Width)const//�w�肵���s�A��̗v�f����Ƃ��Ďw�肵���s���A�񐔂͈̔͂̍s��𔲂��o���B
{
  int buf_row=num_of_row;
  int buf_column=num_of_column;

  if(Row<1 || Column<1 || Row>buf_row || Column>buf_column || (Row+Height-1)>buf_row || (Column+Width-1)>buf_column)
  {
    printf("Cannot access.pull(%d,%d,%d,%d)\n",Row,Column,Height,Width);
    system("pause");
    exit(EXIT_FAILURE);
  }

  int i,j;
  MATRIX<mTYPE> pulled(Height,Width);
  mTYPE *buf_p_value=pulled.value;

  mTYPE *buf_value=value;

  for(i=0;i<Height;i++)
  {
    for(j=0;j<Width;j++)
    {
      buf_p_value[i*Width+j]=buf_value[(i+(Row-1))*buf_column+(j+(Column-1))];
    }
  }

  pulled.value=buf_p_value;

  return pulled;  
}

template<typename mTYPE>
MATRIX<mTYPE> MATRIX<mTYPE>::pullrow(int Row)const//�w�肵���s��Ԃ��B
{
  int buf_row=num_of_row;
  int buf_column=num_of_column;

  if(Row<1 || Row>buf_row)
  {
    printf("Cannot access.pullrow(%d)\n",Row);
    system("pause");
    exit(EXIT_FAILURE);
  }

  int i;
  MATRIX<mTYPE> pulled(1,buf_column);
  mTYPE *buf_p_value=pulled.value;

  mTYPE *buf_value=value;

  for(i=0;i<buf_column;i++)
  {
    buf_p_value[i]=buf_value[(Row-1)*buf_column+i];
  }

  pulled.value=buf_p_value;

  return pulled;  
}

template<typename mTYPE>
MATRIX<mTYPE> MATRIX<mTYPE>::pullcol(int Column)const//�w�肵���s��Ԃ��B
{
  int buf_row=num_of_row;
  int buf_column=num_of_column;

  if(Column<1 || Column>buf_column)
  {
    printf("Cannot access.pullcol(%d)\n",Column);
    system("pause");
    exit(EXIT_FAILURE);
  }

  int i;
  MATRIX<mTYPE> pulled(buf_row,1);

  mTYPE *buf_p_value=pulled.value;

  mTYPE *buf_value=value;

  for(i=0;i<buf_row;i++)
  {
    buf_p_value[i]=buf_value[i*buf_column+(Column-1)];
  }

  pulled.value=buf_p_value;

  return pulled;  
}

template<typename mTYPE>
mTYPE MATRIX<mTYPE>::max() const
{
  mTYPE max_value;
  int i;
  int buf_data=num_of_row*num_of_column;

  mTYPE *buf_value=value;

  max_value=buf_value[0];
  for(i=0;i<buf_data;i++)
  {
    if(max_value<buf_value[i]) max_value=buf_value[i];
  }

  return max_value;
}

template<typename mTYPE>
mTYPE MATRIX<mTYPE>::min() const
{
  mTYPE min_value;
  int i;
  int buf_data=num_of_row*num_of_column;

  mTYPE *buf_value=value;

  min_value=buf_value[0];
  for(i=0;i<buf_data;i++)
  {
    if(min_value>buf_value[i]) min_value=buf_value[i];
  }

  return min_value;
}

template<typename mTYPE>
double MATRIX<mTYPE>::average()const{//�傫�ȉ摜�S�̂ɑ΂��Ă̓���͂��s��������
    int i;
    int buf_data=num_of_row*num_of_column;
    mTYPE *buf_value=value;

    double average=0.0;
    for(i=0;i<buf_data;i++){
      average+=buf_value[i];
    }

    average/=buf_data;
    return average;
}

template<typename mTYPE>
double MATRIX<mTYPE>::stddev()const{//�s��̑S�v�f�̕W���΍���Ԃ��B
  int i;
  int buf_data=num_of_row*num_of_column;
  mTYPE *buf_value=value;

  double sum_square=0.0;
  for(i=0;i<buf_data;i++){
    sum_square+=value[i]*value[i];
  }
  
  double ave=average();
  
  double stddev=sqrt((sum_square-buf_data*ave*ave)/(buf_data-1));

  return stddev;
}


//�L���X�g����MATRIX��Ԃ����߂̊֐��@���ɒu���^�C�v�̂��́B
template<typename mTYPE>
MATRIX<short> MATRIX<mTYPE>::m_short()const //short�^��MATRIX��Ԃ�
{
  int buf_row=num_of_row;
  int buf_column=num_of_column;
  int buf_data=buf_row*buf_column;

  MATRIX<short> ans(buf_row,buf_column);

  short *buf_ans_value=ans.pVALUE();
  mTYPE *buf_value=value;

  int i;
  for(i=0;i<buf_data;i++)
  {
    buf_ans_value[i]=short(buf_value[i]);
  }

  return ans;
}

template<typename mTYPE>
MATRIX<int> MATRIX<mTYPE>::m_int()const //int�^��MATRIX��Ԃ�
{
  int buf_row=num_of_row;
  int buf_column=num_of_column;
  int buf_data=buf_row*buf_column;

  MATRIX<int> ans(buf_row,buf_column);

  int *buf_ans_value=ans.pVALUE();
  mTYPE *buf_value=value;

  int i;
  for(i=0;i<buf_data;i++)
  {
    buf_ans_value[i]=int(buf_value[i]);
  }

  return ans;
}

template<typename mTYPE>
MATRIX<float> MATRIX<mTYPE>::m_float()const //float�^��MATRIX��Ԃ�
{
  int buf_row=num_of_row;
  int buf_column=num_of_column;
  int buf_data=buf_row*buf_column;

  MATRIX<float> ans(buf_row,buf_column);

  float *buf_ans_value=ans.pVALUE();
  mTYPE *buf_value=value;

  int i;
  for(i=0;i<buf_data;i++)
  {
    buf_ans_value[i]=float(buf_value[i]);
  }

  return ans;
}

template<typename mTYPE>
MATRIX<double> MATRIX<mTYPE>::m_double()const //double�^��MATRIX��Ԃ�
{
  int buf_row=num_of_row;
  int buf_column=num_of_column;
  int buf_data=buf_row*buf_column;

  MATRIX<double> ans(buf_row,buf_column);

  double *buf_ans_value=ans.pVALUE();
  mTYPE *buf_value=value;

  int i;
  for(i=0;i<buf_data;i++)
  {
    buf_ans_value[i]=double(buf_value[i]);
  }

  return ans;
}


//������Z�q�I�[�o�[���[�h
template<typename mTYPE>
MATRIX<mTYPE>& MATRIX<mTYPE>::operator =(const MATRIX<mTYPE> &matrix)//matrix��������B
{
  int buf_m_row=matrix.num_of_row;
  int buf_m_column=matrix.num_of_column;
  int buf_row=buf_m_row;
  int buf_column=buf_m_column;

  mTYPE *buf_value;
  buf_value=(mTYPE *)malloc(sizeof(mTYPE)*buf_m_row*buf_m_column);
  if(buf_value==NULL)
  {
    printf("MEMORY ALLOCATION ERROR.\n");
    system("pause");
    exit(EXIT_FAILURE);
  }

  mTYPE *buf_m_value;
  buf_m_value=matrix.value;

  int i,j;
  for(i=0;i<buf_row;i++)
  {
    for(j=0;j<buf_column;j++)
    {
      buf_value[i*buf_column+j]=buf_m_value[i*buf_m_column+j];
    }
  }

  num_of_row=buf_row;
  num_of_column=buf_column;
  free(value);
  value=buf_value;

  return *this;
}



/*       MATRIX�N���X�̂��߂̊֐��Q�̒�`       */

//�L���X�g����MATRIX��Ԃ����߂̊֐��@�O�ɒu���^�C�v�̂��́B
template<typename mTYPE>
MATRIX<short> m_short(const MATRIX<mTYPE>& A)//short�^��MATRIX��Ԃ�
{
  int buf_row=A.Row();
  int buf_column=A.Col();
  int buf_data=buf_row*buf_column;

  MATRIX<short> ans(buf_row,buf_column);

  short *buf_ans_value=ans.pVALUE();
  mTYPE *buf_A_value=A.pVALUE();

  int i;
  for(i=0;i<buf_data;i++)
  {
    buf_ans_value[i]=short(buf_A_value[i]);
  }

  return ans;
}

template<typename mTYPE>
MATRIX<int> m_int(const MATRIX<mTYPE>& A)//int�^��MATRIX��Ԃ�
{
  int buf_row=A.Row();
  int buf_column=A.Col();
  int buf_data=buf_row*buf_column;

  MATRIX<int> ans(buf_row,buf_column);

  int *buf_ans_value=ans.pVALUE();
  mTYPE *buf_A_value=A.pVALUE();

  int i;
  for(i=0;i<buf_data;i++)
  {
    buf_ans_value[i]=int(buf_A_value[i]);
  }

  return ans;
}

template<typename mTYPE>
MATRIX<float> m_float(const MATRIX<mTYPE>& A)//float�^��MATRIX��Ԃ�
{
  int buf_row=A.Row();
  int buf_column=A.Col();
  int buf_data=buf_row*buf_column;

  MATRIX<float> ans(buf_row,buf_column);

  float *buf_ans_value=ans.pVALUE();
  mTYPE *buf_A_value=A.pVALUE();

  int i;
  for(i=0;i<buf_data;i++)
  {
    buf_ans_value[i]=float(buf_A_value[i]);
  }

  return ans;
}

template<typename mTYPE>
MATRIX<double> m_double(const MATRIX<mTYPE>& A)//double�^��MATRIX��Ԃ�
{
  int buf_row=A.Row();
  int buf_column=A.Col();
  int buf_data=buf_row*buf_column;

  MATRIX<double> ans(buf_row,buf_column);

  double *buf_ans_value=ans.pVALUE();
  mTYPE *buf_A_value=A.pVALUE();

  int i;
  for(i=0;i<buf_data;i++)
  {
    buf_ans_value[i]=double(buf_A_value[i]);
  }

  return ans;
}

template<typename type,typename mTYPE>
MATRIX<type> m_cast(const MATRIX<mTYPE> &A)//[�s��]�̕ϐ��^���w�肵��[�^]�ɕύX�����s���Ԃ��Btype�̌���͔͂���ł��Ȃ��̂Œ��ӁB
{
  int buf_row=A.Row();
  int buf_column=A.Col();
  int buf_data=buf_row*buf_column;

  MATRIX<type> ans(buf_row,buf_column);

  type *buf_ans_value=ans.pVALUE();
  mTYPE *buf_A_value=A.pVALUE();

  int i;
  for(i=0;i<buf_data;i++)
  {
    buf_ans_value[i]=type(buf_A_value[i]);
  }

  return ans;  
}



//�s��̂��߂̉��Z�֐�
template<typename mTYPE1,typename mTYPE2>
MATRIX<double> Multiplication(const MATRIX<mTYPE1> &A,const MATRIX<mTYPE2> &B)//�s�񓯎m�̂����Z���s���B
{
  int buf_A_row=A.Row();
  int buf_A_column=A.Col();
  int buf_B_column=B.Col();

  if(buf_A_column!=B.Row())
  {
    printf("Size error.Multiplication((%d,%d),(%d,%d))\n",buf_A_row,buf_A_column,B.Row(),buf_B_column);
    system("pause");
    exit(EXIT_FAILURE);
  }

  int i,j,k;
  double buf_ans_Value;

  int buf_ans_row=buf_A_row;
  int buf_ans_column=buf_B_column;
  MATRIX<double> ans(buf_ans_row,buf_ans_column);

  double *buf_ans_value=ans.pVALUE();
  mTYPE1 *buf_A_value=A.pVALUE();
  mTYPE2 *buf_B_value=B.pVALUE();

  for(i=0;i<buf_ans_row;i++)
  {
    for(j=0;j<buf_ans_column;j++)
    {
      buf_ans_Value=0.0;
      for(k=0;k<buf_A_column;k++)
      {
        buf_ans_Value+=buf_A_value[i*buf_A_column+k]*buf_B_value[k*buf_B_column+j];
      }
      buf_ans_value[i*buf_ans_column+j]=buf_ans_Value;
    }
  }
  
  return ans; 
}

template<typename mTYPE>
MATRIX<double> exp(const MATRIX<mTYPE> &A)//�s��̊e�v�f��exp��������s���Ԃ��B
{
  int i;
  int buf_row=A.Row();
  int buf_column=A.Col();
  int buf_data=buf_row*buf_column;
 
  MATRIX<double> ans(buf_row,buf_column);

  double *buf_ans_value=ans.pVALUE();
  mTYPE *buf_A_value=A.pVALUE();

  for(i=0;i<buf_data;i++)
  {
    buf_ans_value[i]=exp(buf_A_value[i]);
  }
  
  return ans; 
}

template<typename mTYPE>
MATRIX<double> log(const MATRIX<mTYPE> &A)//�s��̊e�v�f�̎��R�ΐ���������s���Ԃ��B
{
  int i;
  int buf_row=A.Row();
  int buf_column=A.Col();
  int buf_data=buf_row*buf_column;
 
  MATRIX<double> ans(buf_row,buf_column);

  double *buf_ans_value=ans.pVALUE();
  mTYPE *buf_A_value=A.pVALUE();

  for(i=0;i<buf_data;i++)
  {
    buf_ans_value[i]=log(buf_A_value[i]);
  }
  
  return ans; 
}


//���Z�q�I�[�o�[���[�h
template<typename mTYPE>
MATRIX<mTYPE> operator -(const MATRIX<mTYPE>& A)//�v�f�̕����𔽓]�����s���Ԃ��B
{
  int buf_row=A.Row();
  int buf_column=A.Col();
  int buf_data=buf_row*buf_column;

  MATRIX<mTYPE> ans(buf_row,buf_column);

  mTYPE *buf_ans_value=ans.pVALUE();

  mTYPE *buf_A_value=A.pVALUE();

  int i;
  for(i=0;i<buf_data;i++)
  {
    buf_ans_value[i]=mTYPE((-1.0)*buf_A_value[i]);
  }
 
  return ans;
}

template<typename mTYPE>
MATRIX<double> operator +(double r,const MATRIX<mTYPE> &A)//���l���e�v�f�ɑ������s���Ԃ�
{
  int i;

  int buf_row=A.Row();
  int buf_column=A.Col();
  int buf_data=buf_row*buf_column;

  MATRIX<double> ans(buf_row,buf_column);

  double *buf_ans_value=ans.pVALUE();

  mTYPE *buf_A_value=A.pVALUE();

  for(i=0;i<buf_data;i++)
  {
    buf_ans_value[i]=double(r+buf_A_value[i]);
  }

  return ans; 
}

template<typename mTYPE>
MATRIX<double> operator +(const MATRIX<mTYPE> &A,double r)//�e�v�f�ɐ��l�𑫂����s���Ԃ�
{
  int i;

  int buf_row=A.Row();
  int buf_column=A.Col();
  int buf_data=buf_row*buf_column;

  MATRIX<double> ans(buf_row,buf_column);

  double *buf_ans_value=ans.pVALUE();

  mTYPE *buf_A_value=A.pVALUE();

  for(i=0;i<buf_data;i++)
  {
    buf_ans_value[i]=double(buf_A_value[i]+r);
  }

  return ans; 
}

template<typename mTYPE1,typename mTYPE2>
MATRIX<double> operator +(const MATRIX<mTYPE1> &A,const MATRIX<mTYPE2> &B)//�T�C�Y�̓������s��̊e�v�f�̘a�����
{
  if(A.Row()!=B.Row() || A.Col()!=B.Col())
  {
    printf("Size error.+\n");
    system("pause");
    exit(EXIT_FAILURE);
  }

  int buf_row=A.Row();
  int buf_column=A.Col();
  int buf_data=buf_row*buf_column;

  MATRIX<double> ans(buf_row,buf_column);

  mTYPE1 *buf_A_value=A.pVALUE();
  mTYPE2 *buf_B_value=B.pVALUE();
  double *buf_ans_value=ans.pVALUE();

  int i;
  for(i=0;i<buf_data;i++)
  {
    buf_ans_value[i]=double(buf_A_value[i]+buf_B_value[i]);
  } 

  return ans;
}

template<typename mTYPE>
MATRIX<double> operator -(double r,const MATRIX<mTYPE> &A)//���l����e�v�f���������s���Ԃ�
{
  int i;

  int buf_row=A.Row();
  int buf_column=A.Col();
  int buf_data=buf_row*buf_column;

  MATRIX<double> ans(buf_row,buf_column);

  double *buf_ans_value=ans.pVALUE();

  mTYPE *buf_A_value=A.pVALUE();

  for(i=0;i<buf_data;i++)
  {
    buf_ans_value[i]=double(r-buf_A_value[i]);
  }

  return ans; 
}

template<typename mTYPE>
MATRIX<double> operator -(const MATRIX<mTYPE> &A,double r)//�e�v�f���琔�l���������s���Ԃ�
{
  int i;

  int buf_row=A.Row();
  int buf_column=A.Col();
  int buf_data=buf_row*buf_column;

  MATRIX<double> ans(buf_row,buf_column);

  double *buf_ans_value=ans.pVALUE();

  mTYPE *buf_A_value=A.pVALUE();

  for(i=0;i<buf_data;i++)
  {
    buf_ans_value[i]=double(buf_A_value[i]-r);
  }

  return ans; 
}

template<typename mTYPE1,typename mTYPE2>
MATRIX<double> operator -(const MATRIX<mTYPE1> &A,const MATRIX<mTYPE2> &B)//�T�C�Y�̓������s��̊e�v�f�̍������
{
  if(A.Row()!=B.Row() || A.Col()!=B.Col())
  {
    printf("Size error.+\n");
    system("pause");
    exit(EXIT_FAILURE);
  }

  int buf_row=A.Row();
  int buf_column=A.Col();
  int buf_data=buf_row*buf_column;

  MATRIX<double> ans(buf_row,buf_column);

  mTYPE1 *buf_A_value=A.pVALUE();
  mTYPE2 *buf_B_value=B.pVALUE();
  double *buf_ans_value=ans.pVALUE();

  int i;
  for(i=0;i<buf_data;i++)
  {
    buf_ans_value[i]=double(buf_A_value[i]-buf_B_value[i]);
  } 
  
  return ans;
}

template<typename mTYPE>
MATRIX<double> operator *(double r,const MATRIX<mTYPE> &A)//���l���e�v�f�Ɋ|�����s���Ԃ�
{
  int i;

  int buf_row=A.Row();
  int buf_column=A.Col();
  int buf_data=buf_row*buf_column;

  MATRIX<double> ans(buf_row,buf_column);

  double *buf_ans_value=ans.pVALUE();

  mTYPE *buf_A_value=A.pVALUE();

  for(i=0;i<buf_data;i++)
  {
    buf_ans_value[i]=double(r*buf_A_value[i]);
  }

  return ans; 
}

template<typename mTYPE>
MATRIX<double> operator *(const MATRIX<mTYPE> &A,double r)//�e�v�f�ɐ��l���|�����s���Ԃ�
{
  int i;

  int buf_row=A.Row();
  int buf_column=A.Col();
  int buf_data=buf_row*buf_column;

  MATRIX<double> ans(buf_row,buf_column);

  double *buf_ans_value=ans.pVALUE();

  mTYPE *buf_A_value=A.pVALUE();

  for(i=0;i<buf_data;i++)
  {
    buf_ans_value[i]=double(buf_A_value[i]*r);
  }

  return ans; 
}

template<typename mTYPE1,typename mTYPE2>
MATRIX<double> operator *(const MATRIX<mTYPE1> &A,const MATRIX<mTYPE2> &B)//�T�C�Y�̓������s��̊e�v�f�̐ς����
{
  if(A.Row()!=B.Row() || A.Col()!=B.Col())
  {
    printf("Size error.+\n");
    system("pause");
    exit(EXIT_FAILURE);
  }

  int buf_row=A.Row();
  int buf_column=A.Col();
  int buf_data=buf_row*buf_column;

  MATRIX<double> ans(buf_row,buf_column);

  mTYPE1 *buf_A_value=A.pVALUE();
  mTYPE2 *buf_B_value=B.pVALUE();
  double *buf_ans_value=ans.pVALUE();

  int i;
  for(i=0;i<buf_data;i++)
  {
    buf_ans_value[i]=double(buf_A_value[i]*buf_B_value[i]);
  } 
  
  return ans;
}

template<typename mTYPE>
MATRIX<double> operator /(double r,const MATRIX<mTYPE> &A)//���l���e�v�f�Ŋ������s���Ԃ�
{
  int i;

  int buf_row=A.Row();
  int buf_column=A.Col();
  int buf_data=buf_row*buf_column;

  MATRIX<double> ans(buf_row,buf_column);

  double *buf_ans_value=ans.pVALUE();

  mTYPE *buf_A_value=A.pVALUE();

  for(i=0;i<buf_data;i++)
  {
    buf_ans_value[i]=double(r/buf_A_value[i]);
  }

  return ans; 
}

template<typename mTYPE>
MATRIX<double> operator /(const MATRIX<mTYPE> &A,double r)//�e�v�f�𐔒l�Ŋ������s���Ԃ�
{
  int i;

  int buf_row=A.Row();
  int buf_column=A.Col();
  int buf_data=buf_row*buf_column;

  MATRIX<double> ans(buf_row,buf_column);

  double *buf_ans_value=ans.pVALUE();

  mTYPE *buf_A_value=A.pVALUE();

  for(i=0;i<buf_data;i++)
  {
    buf_ans_value[i]=double(buf_A_value[i]/r);
  }

  return ans; 
}

template<typename mTYPE1,typename mTYPE2>
MATRIX<double> operator /(const MATRIX<mTYPE1> &A,const MATRIX<mTYPE2> &B)//�T�C�Y�̓������s��̊e�v�f�̏������
{
  if(A.Row()!=B.Row() || A.Col()!=B.Col())
  {
    printf("Size error.+\n");
    system("pause");
    exit(EXIT_FAILURE);
  }

  int buf_row=A.Row();
  int buf_column=A.Col();
  int buf_data=buf_row*buf_column;

  MATRIX<double> ans(buf_row,buf_column);

  mTYPE1 *buf_A_value=A.pVALUE();
  mTYPE2 *buf_B_value=B.pVALUE();
  double *buf_ans_value=ans.pVALUE();

  int i;
  for(i=0;i<buf_data;i++)
  {
    buf_ans_value[i]=double(buf_A_value[i])/buf_B_value[i];//���̏ꍇ�����������鐮���̊ۂߍ��݂��S�z�����̂ŁAdouble�ł̃L���X�g����ڂ̒l�ɂ�����B
  } 
  
  return ans;
}