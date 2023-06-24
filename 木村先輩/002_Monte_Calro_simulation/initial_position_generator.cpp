#include<cstdlib>
#include<cstdio>
#include<math.h>
#include<time.h>

FILE *rand_no;
double r0 = 1.0e-3;			//MOT radius[m]
double sigma3 = r0;			//triple standard deviation of normal distribution[m]
double sigma = sigma3/3.0;	//standard deviation of normal distribution[m]
int SAMPLE = 1000; 			//number of sample atoms
char var[30];

double gauss(double x);
double gauss_pos(double x);
double Uniform_1(void);
double Uniform_2(void);

int main(){

	srand((unsigned)time(NULL));	//a magic phrase to make random numbers based on time

		for(int i=1; i<=1000; i++){
			sprintf( var, "G:/graduation/rand_gauss2/rand_gauss%d.csv", i );							//CHANGE THE ADDRESS WEHERE YOU WANT TO GENERATE CSV FILES!!!
			if( (rand_no=fopen(var,"w"))==NULL ){
				printf( "'%s' could not be opened in writing mode.", var );
				return 0;
			}

			for(int i=1; i<=SAMPLE; i++){
				fprintf( rand_no,"%.16e,%.16e,%.16e\n", gauss_pos(sigma3), Uniform_1(), Uniform_2() );	//output radii, azimuths and elevations into csv files
			}

			fflush(rand_no); fclose(rand_no);
		}

	return 0;
}


//Gaussian function (3*sigma = r0)
double gauss(double x){
	return sqrt( 1.0/(2.0*M_PI*r0*r0/9.0) )*exp( -x*x/(2.0*r0*r0/9.0) );
}

double gauss_pos(double x){ //3sigma should be substituted into x
	double xi,yj;
	do{
		xi = 2.0*x*rand()/RAND_MAX-x;//generate uniform random no. xi between [-x,x] (x=n*wr/2)
		yj = 1205.0*rand()/RAND_MAX;//1205.0 is a value that is a little bit greater than the max of boltz()
	}while( gauss(xi)<yj );
	return xi;
}

double Uniform_1(void){
	double i;
		while(1){
			i = M_PI * (double)rand() / (double)RAND_MAX;//i is in (0,1)
			if( i==M_PI ) continue;
			else break;
		} 
		return i;
}

double Uniform_2(void){
	return M_PI * (double)rand() / (double)RAND_MAX;//returned number is in [0,1]
}