/*

テンプレートとしては[short]、[int]、[float]、[double]に対応。
演算+、-、*、/、exp、logなどはほとんど倍精度計算され、結果がdouble型のMATRIXで返るので用意した関数でキャストすること。
バイナリーファイル開閉のためにプロジェクトのプロパティで文字セットをUnicode文字セットからマルチバイト文字セットに変更すること。
hamaguchi

*/

#pragma once

#define _CRT_SECURE_NO_WARNINGS //コンパイル時のfopenなどの警告を無視するため。
#define NOMINMAX //windows.hにmax、min関数を定義させないため。

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string>
#include <windows.h>
using namespace std;


/*       MATRIXクラスとメンバ関数群の宣言       */
template<typename mTYPE>
class MATRIX
{
private:
//行列の情報
  int num_of_row;//行数。
  int num_of_column;//列数。
  mTYPE *value;//行列要素のポインタ。


public:
//コンストラクタとデストラクタ
  MATRIX();//引数なしのコンストラクタ。
  MATRIX(int,int);//指定した[行数]、[列数]の行列を用意する。値はすべて不定になる。
  MATRIX(const MATRIX&);//[行列]のコピーコンストラクタ。変数の型が同じであればコピーできる。
  ~MATRIX();//デストラクタ。

//アクセス違反に耐性のあるプログラムに利用するとよい。行列のための演算関数はfriendではなくこのポインタを使っている。慣れた人向け。
  mTYPE* pVALUE(int Row=1,int Column=1)const{return &value[(Row-1)*num_of_column+(Column-1)];};//private変数のvalueポインタを取り出して直接操作できる。


//行列の全般に関する操作
  void disp(char*)const;//行列をコンソールに表示する。たとえば["%d"]など[書式指定]してやるとその通りに表示する。\t、\nはなくてよい。
  void read(FILE*,int,int);//[ファイルポインタ]を渡して指定した[行数]、[列数]の行列とみなして読み込む。テキストファイルのみ対応。
  void read(const char*,int,int);//[ファイル名]を渡して指定した[行数]、[列数]の行列とみなして読み込む。バイナリーを読み込むときは書式は無視される。
  void write(FILE*,char*) const;//[ファイルポインタ]を渡して、たとえば["%d"]など[書式指定]してやると行列をファイルに書き込む。テキストファイルのみ対応。
  void write(const char*,char*) const;//[ファイル名]を指定して、たとえば["%d"]など書式指定してやると行列をファイルに書き込む。バイナリーを書き込むときは書式は無視される。
  void resize(int,int);//指定した[行数]、[列数]の行列に変更する。数値はすべて不定になる。
  void other_row_col(int,int);//要素の数は変えないが、指定した[行数]、[列数]に変更する。N×M行列を(N*M)×1行列にするなどの用途。
  int Row() const;//行数を返す。
  int Col() const;//列数を返す。
  int Data() const;//行列要素の数を返す。
  void transpose();//行列を転置する。
  void invert_row();//行方向の並びを反転する。
  void invert_col();//列方向の並びを反転する。
  void invert_endian();//エンディアン変換を行う。バイナリーから読み込んだ場合はエンディアンに注意すること。

//行列要素の変更
  void in(int,int,mTYPE);//指定した[行番号]、[列番号]の要素を与えた[数値]で書き換える。アクセス違反判定付き。
  void constant(mTYPE);//すべての要素を指定した[数値]にする。
  void unit();//行列を単行列にする。ただし正方行列のみ。
  void put(int,int,const MATRIX&);//指定した[行番号]、[列番号]の位置を左上端の基準として、与えた[行列]で書き換える。
  void putrow(int,const MATRIX&);//指定した[行番号]の行を与えた[行列]で書き換える。ただし行ベクトルに限る。//20161227なぜかうまく動かないことがある
  void putcol(int,const MATRIX&);//指定した[列番号]の列を与えた[行列]で書き換える。ただし列ベクトルに限る。//20161227なぜかうまく動かないことがある
  void addrow(int,int Height=1,mTYPE Value=0);//指定した[行番号]の行の次に指定した[高さ=行数]で指定した[数値]の行を追加する。[行番号]ゼロで先頭に追加。
  void addrow(int,const MATRIX&);//指定した[行番号]の行の次の行から与えた[行列]を追加する。[行番号]ゼロで先頭から追加。
  void addcol(int,int Width=1,mTYPE Value=0);////指定した[列番号]の列の次に指定した[幅=列数]で指定した[数値]の行を追加する。[列番号]ゼロで先頭に追加。
  void addcol(int,const MATRIX&);//指定した[列番号]の列の次の行から与えた[行列]を追加する。[列番号]ゼロで先頭から追加。
  void delrow(int,int Height=1);//指定した[行番号]の行から指定した[高さ=行数]の行を削除して詰める。
  void delcol(int,int Width=1);//指定した[列番号]の列から指定した[幅=列数]の列を削除して詰める。

//行列の出力
  mTYPE out(int,int) const;//指定した[行番号]、[列番号]の位置の値を返す。アクセス違反判定付き。
  MATRIX pull(int,int,int,int) const;//指定した[行番号]、[列番号]の位置を左上端の基準として、指定した[高さ=行数]、[幅=列数]の範囲の行列を返す。
  MATRIX pullrow(int)const;//指定した[行番号]の行を返す。
  MATRIX pullcol(int)const;//指定した[列番号]の列を返す。
  mTYPE max()const;//要素の最大値を返す。 
  mTYPE min()const;//要素の最小値を返す。 
  double average()const;//行列の全要素の平均を返す。//大きな画像全体に対しての動作はやや不安がある
  double stddev()const;//行列の全要素の標準偏差を返す。//大きな画像全体に対しての動作はやや不安がある
//  MATRIX abs(const MATRIX&);//行列の各要素の絶対値を取った行列を返す。未実装2016/11/09hamaguchi

//キャストしたMATRIXを返すための関数　後ろに置くタイプのもの。
  MATRIX<short> m_short()const;//行列の変数型をshortに変更した行列を返す。
  MATRIX<int> m_int()const;//行列の変数型をintに変更した行列を返す。
  MATRIX<float> m_float()const;//行列の変数型をfloatに変更した行列を返す。 
  MATRIX<double> m_double()const;//行列の変数型をdoubleに変更した行列を返す。 

//代入演算子オーバーロード
  MATRIX& operator = (const MATRIX&);//左辺の行列に右辺の行列を代入する。暗黙の型変換はさせないため行列の変数型が同じ場合に限る。キャストも使うこと。
};



/*       MATRIXクラスのための関数群のプロトタイプ宣言(定義は下の方に)       */

//キャストしたMATRIXを返すための関数　前に置くタイプのもの
template<typename mTYPE>
MATRIX<short> m_short(const MATRIX<mTYPE>&);//[行列]の変数型をshortに変更した行列を返す。

template<typename mTYPE>
MATRIX<int> m_int(const MATRIX<mTYPE>&);//[行列]の変数型をintに変更した行列を返す。

template<typename mTYPE>
MATRIX<float> m_float(const MATRIX<mTYPE>&);//[行列]の変数型をfloatに変更した行列を返す。

template<typename mTYPE>
MATRIX<double> m_double(const MATRIX<mTYPE>&);//[行列]の変数型をdoubleに変更した行列を返す。

template<typename type,typename mTYPE>
MATRIX<type> m_cast(const MATRIX<mTYPE>&);//[行列]の変数型を指定した[型]に変更した行列を返す。typeの誤入力は判定できないので注意。


//行列のための演算関数
template<typename mTYPE1,typename mTYPE2>
MATRIX<double> Multiplication(const MATRIX<mTYPE1>&,const MATRIX<mTYPE2>&);//[行列]、[行列]のかけ算を行い、結果の行列をdouble型で返す。

template<typename mTYPE>
MATRIX<double> exp(const MATRIX<mTYPE>&);//[行列]の各要素のexpを取った行列をdouble型で返す。

template<typename mTYPE>
MATRIX<double> log(const MATRIX<mTYPE>&);//[行列]の各要素のlogを取った行列をdouble型で返す。


//演算子オーバーロード
template<typename mTYPE>
MATRIX<mTYPE> operator -(const MATRIX<mTYPE>&);//[行列]の各要素の符号を反転した行列を与えられた行列の型で返す。//戻り値のmTYPEの利便性について要検討

template<typename mTYPE>
MATRIX<double> operator +(double,const MATRIX<mTYPE>&);//[数値]を[行列]の各要素に足した行列をdouble型で返す。//戻り値のdoubleの利便性について要検討

template<typename mTYPE>
MATRIX<double> operator +(const MATRIX<mTYPE>&,double);//[行列]の各要素に[数値]を足した行列をdouble型で返す。//戻り値のdoubleの利便性について要検討

template<typename mTYPE1,typename mTYPE2>
MATRIX<double> operator +(const MATRIX<mTYPE1>&,const MATRIX<mTYPE2>&);//[行列]と[行列]の各要素を足した行列をdouble型で返す。//戻り値のmTYPEの利便性について要検討

template<typename mTYPE>
MATRIX<double> operator -(double,const MATRIX<mTYPE>&);//[数値]から[行列]の各要素を引いた行列をdouble型で返す。//戻り値のdoubleの利便性について要検討

template<typename mTYPE>
MATRIX<double> operator -(const MATRIX<mTYPE>&,double);//[行列]の各要素から[数値]を引いた行列をdouble型で返す。//戻り値のdoubleの利便性について要検討

template<typename mTYPE1,typename mTYPE2>
MATRIX<double> operator -(const MATRIX<mTYPE1>&,const MATRIX<mTYPE2>&);//[行列]から[行列]の各要素を引いたを行列をdouble型で返す。//戻り値のdoubleの利便性について要検討

template<typename mTYPE>
MATRIX<double> operator *(double,const MATRIX<mTYPE>&);//[数値]を[行列]の各要素に掛けた行列をdouble型で返す。//戻り値のdoubleの利便性について要検討

template<typename mTYPE>
MATRIX<double> operator *(const MATRIX<mTYPE>&,double);//[行列]の各要素に[数値]を掛けた行列をdouble型で返す。//戻り値のdoubleの利便性について要検討

template<typename mTYPE1,typename mTYPE2>
MATRIX<double> operator *(const MATRIX<mTYPE1>&,const MATRIX<mTYPE2>&);//[行列]と[行列]の各要素を掛けた行列をdouble型で返す。//戻り値のdoubleの利便性について要検討

template<typename mTYPE>
MATRIX<double> operator /(double,const MATRIX<double>&);//[数値]を[行列]の各要素で割った行列をdouble型で返す。//戻り値のdoubleの利便性について要検討

template<typename mTYPE>
MATRIX<double> operator /(const MATRIX<double>&,double);//[行列]の各要素を[数値]で割った行列をdouble型で返す。//戻り値のdoubleの利便性について要検討

template<typename mTYPE1,typename mTYPE2>
MATRIX<double> operator /(const MATRIX<mTYPE1>&,const MATRIX<mTYPE2>&);//[数値]を[行列]の各要素で割った行列をdouble型で返す。//戻り値のmTYPEの利便性について要検討



/*       MATRIXクラスのメンバ関数群の定義       */

//コンストラクタとデストラクタ
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


//行列の全般に関する操作
template<typename mTYPE>
void MATRIX<mTYPE>::disp(char* format)const//行列をコンソールに表示する。
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


//.txt拡張子の判定ここから
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
//.txt拡張子の判定おわり

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


//.txt拡張子の判定ここから
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
//.txt拡張子の判定おわり

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
void MATRIX<mTYPE>::resize(int Row,int Column)//行列のサイズを変更する。数値はすべて不定になる。
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
void MATRIX<mTYPE>::transpose()//行列を転置する。
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
  for(i=0;i<buf_row;i++)//転置後の配列に並べ替える
  {
    for(j=0;j<buf_column;j++)
    {
      buf_value[j*buf_row+i]=value[i*buf_column+j];
    }
  }

  //並べ替えたものを新しい行列要素とする
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

  //並べ替えたものを新しい行列要素とする
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

  //並べ替えたものを新しい行列要素とする
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


//行列要素の変更
template<typename mTYPE>
void MATRIX<mTYPE>::in(int Row,int Column,mTYPE Value)//指定した行、列の要素を与えた値で書き換える。
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
void MATRIX<mTYPE>::constant(mTYPE Value)//すべての要素を指定した数値にする
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
void MATRIX<mTYPE>::unit()//行列を単行列にする。ただし正方行列のみ。
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
void MATRIX<mTYPE>::put(int Row,int Column,const MATRIX<mTYPE> &matrix)//指定した行、列の要素を基準として、与えた行列で書き換える。
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
void MATRIX<mTYPE>::putrow(int Row,const MATRIX<mTYPE> &matrix)//指定した行を与えた行列で書き換える。ただし行ベクトルに限る。
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
void MATRIX<mTYPE>::putcol(int Column,const MATRIX<mTYPE> &matrix)//指定した列を与えた行列で書き換える。ただし列ベクトルに限る。
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
void MATRIX<mTYPE>::addrow(int Row,int Height,mTYPE Value)//指定した行の下にValueの行を指定した数で追加する。
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
void MATRIX<mTYPE>::addcol(int Column,int Width,mTYPE Value)//指定した列の下に値Valueの列を指定した数で追加する。
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
  
  transpose();//転置してaddrowを利用

  int Row=Column;  
  int Height=Width;
  buf=buf_row;//転置に合わせて行数と列数を交換
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
      
  transpose();//転置して元に戻す
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
  
  transpose();//転置してdelrowを利用
  int Row=Column;  
  int Height=Width;
  buf=buf_row;//転置に合わせて行数と列数を交換
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


//行列要素の出力
template<typename mTYPE>
mTYPE MATRIX<mTYPE>::out(int Row,int Column)const//指定した行、列の要素を返す。
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
MATRIX<mTYPE> MATRIX<mTYPE>::pull(int Row,int Column,int Height,int Width)const//指定した行、列の要素を基準として指定した行数、列数の範囲の行列を抜き出す。
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
MATRIX<mTYPE> MATRIX<mTYPE>::pullrow(int Row)const//指定した行を返す。
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
MATRIX<mTYPE> MATRIX<mTYPE>::pullcol(int Column)const//指定した行を返す。
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
double MATRIX<mTYPE>::average()const{//大きな画像全体に対しての動作はやや不安がある
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
double MATRIX<mTYPE>::stddev()const{//行列の全要素の標準偏差を返す。
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


//キャストしたMATRIXを返すための関数　後ろに置くタイプのもの。
template<typename mTYPE>
MATRIX<short> MATRIX<mTYPE>::m_short()const //short型のMATRIXを返す
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
MATRIX<int> MATRIX<mTYPE>::m_int()const //int型のMATRIXを返す
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
MATRIX<float> MATRIX<mTYPE>::m_float()const //float型のMATRIXを返す
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
MATRIX<double> MATRIX<mTYPE>::m_double()const //double型のMATRIXを返す
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


//代入演算子オーバーロード
template<typename mTYPE>
MATRIX<mTYPE>& MATRIX<mTYPE>::operator =(const MATRIX<mTYPE> &matrix)//matrixを代入する。
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



/*       MATRIXクラスのための関数群の定義       */

//キャストしたMATRIXを返すための関数　前に置くタイプのもの。
template<typename mTYPE>
MATRIX<short> m_short(const MATRIX<mTYPE>& A)//short型のMATRIXを返す
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
MATRIX<int> m_int(const MATRIX<mTYPE>& A)//int型のMATRIXを返す
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
MATRIX<float> m_float(const MATRIX<mTYPE>& A)//float型のMATRIXを返す
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
MATRIX<double> m_double(const MATRIX<mTYPE>& A)//double型のMATRIXを返す
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
MATRIX<type> m_cast(const MATRIX<mTYPE> &A)//[行列]の変数型を指定した[型]に変更した行列を返す。typeの誤入力は判定できないので注意。
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



//行列のための演算関数
template<typename mTYPE1,typename mTYPE2>
MATRIX<double> Multiplication(const MATRIX<mTYPE1> &A,const MATRIX<mTYPE2> &B)//行列同士のかけ算を行う。
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
MATRIX<double> exp(const MATRIX<mTYPE> &A)//行列の各要素のexpを取った行列を返す。
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
MATRIX<double> log(const MATRIX<mTYPE> &A)//行列の各要素の自然対数を取った行列を返す。
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


//演算子オーバーロード
template<typename mTYPE>
MATRIX<mTYPE> operator -(const MATRIX<mTYPE>& A)//要素の符号を反転した行列を返す。
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
MATRIX<double> operator +(double r,const MATRIX<mTYPE> &A)//数値を各要素に足した行列を返す
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
MATRIX<double> operator +(const MATRIX<mTYPE> &A,double r)//各要素に数値を足した行列を返す
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
MATRIX<double> operator +(const MATRIX<mTYPE1> &A,const MATRIX<mTYPE2> &B)//サイズの等しい行列の各要素の和を取る
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
MATRIX<double> operator -(double r,const MATRIX<mTYPE> &A)//数値から各要素を引いた行列を返す
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
MATRIX<double> operator -(const MATRIX<mTYPE> &A,double r)//各要素から数値を引いた行列を返す
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
MATRIX<double> operator -(const MATRIX<mTYPE1> &A,const MATRIX<mTYPE2> &B)//サイズの等しい行列の各要素の差を取る
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
MATRIX<double> operator *(double r,const MATRIX<mTYPE> &A)//数値を各要素に掛けた行列を返す
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
MATRIX<double> operator *(const MATRIX<mTYPE> &A,double r)//各要素に数値を掛けた行列を返す
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
MATRIX<double> operator *(const MATRIX<mTYPE1> &A,const MATRIX<mTYPE2> &B)//サイズの等しい行列の各要素の積を取る
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
MATRIX<double> operator /(double r,const MATRIX<mTYPE> &A)//数値を各要素で割った行列を返す
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
MATRIX<double> operator /(const MATRIX<mTYPE> &A,double r)//各要素を数値で割った行列を返す
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
MATRIX<double> operator /(const MATRIX<mTYPE1> &A,const MATRIX<mTYPE2> &B)//サイズの等しい行列の各要素の商を取る
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
    buf_ans_value[i]=double(buf_A_value[i])/buf_B_value[i];//この場合だけ整数割る整数の丸め込みが心配されるので、doubleでのキャストを一つ目の値にかける。
  } 
  
  return ans;
}