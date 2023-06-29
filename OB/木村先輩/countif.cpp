#include<cstdlib>
#include<cstdio>
#include<math.h>

FILE *fp;			//file pointer for "f_state.csv"
FILE *fp_w;			//file pointer for "guide_effic_ver.2.csv"
int state=0;		//final state of atom
int trap=0;			//how many atoms were trapped
int still=0;		//how many atoms were still dropping
int t_still=0;		//total of the atoms that are still dropping
int SAMPLE = 1000;	//number of sample atoms
int j=0;
int m=0;
int n=0;

int main(){

	char filenm[100];

	if( (fp_w = fopen("G:/graduation/temperature7/guide_effic_ver.2.csv", "w"))==NULL ){		//open guide_effic_ver.2.csv to write guiding efficiency in it
		fprintf( stderr, "guide_effic_ver.2.csv could not be booted.%s\n" );
		exit(EXIT_FAILURE);
	}

	for(int j=1; j<=100; j+=1){
		state = 0; trap = 0; still = 0;	m = j/10; n = j%10;
		sprintf( filenm, "G:/graduation/detuning4/%d.%dGHz_f_state.csv", m, n );				//store the address of "f_state.csv" into filenm

		if( (fp = fopen(filenm, "r"))==NULL ){													//open f_state.csv to read out final state(loss, trapped or still dropping)
			fprintf( stderr, "%s could not be booted.%s\n", filenm );
			exit(EXIT_FAILURE);
		}

		for(int i=1; i<=SAMPLE; i++){
			fscanf( fp, "%d\n", &state );
			printf("%d %d\n", i, state );
			if(state==1){
				trap++;
			}
			if(state==2){
				still++;
			}
		}

		t_still += still;

		fprintf( fp_w, "%.16e,%d,%d\n", (double)trap/(double)(1000-still)*100.0, trap, still );	//output "guiding efficiency(percent)", "number of trapped atoms" and "number of dropping atoms"

		fclose( fp );
	}
		fclose( fp_w );
		printf("\n%d\n", t_still);


	return 0;
}