//csvファイルの特定列の読み出し方

#include<cstdlib>
#include<cstdio>
#include<math.h>

FILE *file1;
FILE *file2;
double j;

int main(){

	if( ( file1 = fopen("integer2.csv","w"))==NULL ){										//integer2.csvを w (writing mode) で起動。そして、 integer.csvのアドレスをファイルポインター file1 に格納。
		printf("'integer2.csv' could not be opened in reading mode.");
		return 0;
	}

	for(int i=1; i<=100; i++){
		fprintf( file1, "%.16e,%.16e,%.16e\n", (double)i, 2.0*(double)i, 3.0*(double)i );	//integer2.csv に小数 (double)i, 2.0*(double)i, 3.0*(double)i を書き出し。
																							//(double)はキャスト演算子と呼ばれ、変数型を強制的に変更させる演算子。ここでは double 型に変更している。
	}
	fclose( file1 );


	if( ( file2 = fopen("integer2.csv","r"))==NULL ){
		printf("'integer2.csv' could not be opened in reading mode." );
		return 0;
	}
	for(int i=1; i<=100; i++){
		fscanf( file2, "%*lE,%lE,%*lE\n", &j ); 											// * を使うことで、読み飛ばしを可能にする。
		printf( "%.16e\n", j );
	}
	fclose( file2 );

	return 0;
}