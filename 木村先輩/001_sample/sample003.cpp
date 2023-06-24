//fprintf, fscanf の使い方

#include<cstdlib>
#include<cstdio>
#include<math.h>

FILE *file1;
FILE *file2;
int j;

int main(){


	//===How to use fprintf===
	if( ( file1 = fopen("integer.csv","w"))==NULL ){						//integer.csvを w (writing mode) で起動。そして、 integer.csvのアドレスをファイルポインター file1 に格納。
		printf("'integer.csv' could not be opened in reading mode.");
		return 0;
	}

	for(int i=1; i<=100; i++){
		fprintf( file1, "%d\n", i );										//integer.csv に整数 i を書き出し。
	}
	fclose( file1 );




	//===How to use fscanf===
	if( ( file2 = fopen("integer.csv","r"))==NULL ){						//integer.csvを r (reading mode) で起動。そして、 integer.csvのアドレスをファイルポインター file2 に格納。
		printf("'integer.csv' could not be opened in reading mode.");
		return 0;
	}

	for(int i=1; i<=100; i++){
		fscanf( file2, "%d\n", &j );										//integer.csv から整数 i を読み出し、int 型変数 j に格納。
		printf( "%d ", j );													//変数 j を標準出力（コマンドプロンプトに書き出し）
	}




	return 0;
}