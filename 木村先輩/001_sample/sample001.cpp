//0mW_file.csvから 300mW_file.csvまで、ファイル名冒頭のワット数を 100mW 刻みで自動生成

#include<cstdlib>
#include<cstdio>
#include<math.h>

FILE *file;

int main(){

	char file_name[50];								//ファイル名格納用　char型変数

	for(int n=0; n<=300; n+=100){
		//open file
		sprintf( file_name, "%dmW_file.csv", n );	//変数 file_name に %dmW_file.csv を格納（ %d には n が入る）
		if( ( file = fopen(file_name,"w"))==NULL ){
			printf("'%s' could not be opened in writing mode.", file_name);
			return 0;
		}

		//close file
		fclose(file);
	}

	return 0;
}