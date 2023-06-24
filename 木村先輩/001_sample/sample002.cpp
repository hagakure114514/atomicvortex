//動的配列の確保
//通常、配列の長さは定数でなければならないが、ある手法を使うと変数により決められる。

#include<cstdlib>
#include<cstdio>
#include<math.h>

int main(){

	int i = 1000;
	double *x;				//double型 配列 xを定義
	x = new double[i]{};	//長さ i の配列を配列 x 用に確保

	delete[] x;				//配列 x の要素を消去（メモリを解放）

	return 0;
}