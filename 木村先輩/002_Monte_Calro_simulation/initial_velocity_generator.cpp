#include<cstdlib>
#include<cstdio>
#include<math.h>
#include<time.h>

FILE *rand_no;
double Rb = 86.90918053;		//atomic mass of 87Rb[g/mol]
double NA = 6.022140e+23;		//Avogadro constant[/mol]
double mass = Rb/NA/1000.0;		//mass of Rb : Rb/NA/1000[kg]
double k_b = 1.380649e-23;		//Boltzmann const.[J/K]
double temp = 10.0e-6;			//temperature of MOT[K]
int SAMPLE = 1000; 				//number of sample atoms
double *xx;
double *yy;
double *zz;
double max[21];					//a value that is little bit bigger than the maximum of Maxwell-Boltzmann distribution
double broad[21];				//broadning of Gaussian distribution
double maxx = 0.0;
double broadd  = 0.0;
char var[30];
int select;

double max_boltz(double x, double t);
double max_boltz_velo_x(double x, double m, double t);
double max_boltz_velo_y(double x, double m, double t);
double max_boltz_velo_z(double x, double m, double t);

int main(){

	srand((unsigned)time(NULL));//a magic phrase to make random numbers based on time
	printf("If you want to generate random no. for 5-100[uK], input 1.\n");
	printf("If you want to generate random no. for 10[uK], input 2.\n");
	scanf("%d",&select);

	//====generate random no. based on Maxwell-Boltzmann for each temperature (5-100[uK])====
		if( select==1 ){
			max[1] = 18.5;		broad[1] = 0.1;
			max[2] = 13.0;		broad[2] = 0.15;
			max[3] = 10.7;		broad[3] = 0.2;
			max[4] = 9.2;		broad[4] = 0.25;
			max[5] = 8.2;		broad[5] = 0.3;
			max[6] = 7.5;		broad[6] = 0.3;
			max[7] = 7.0;		broad[7] = 0.3;
			max[8] = 6.6;		broad[8] = 0.3;
			max[9] = 6.2;		broad[9] = 0.3;
			max[10] = 5.9;		broad[10] = 0.4;
			max[11] = 5.6;		broad[11] = 0.4;
			max[12] = 5.4;		broad[12] = 0.4;
			max[13] = 5.2;		broad[13] = 0.4;
			max[14] = 5.0;		broad[14] = 0.4;
			max[15] = 4.8;		broad[15] = 0.5;
			max[16] = 4.6;		broad[16] = 0.5;
			max[17] = 4.5;		broad[17] = 0.5;
			max[18] = 4.4;		broad[18] = 0.5;
			max[19] = 4.3;		broad[19] = 0.5;
			max[20] = 4.2;		broad[20] = 0.5;

			for(int i=1; i<=20; i++){
				maxx = max[i]; broadd = broad[i]; temp = (double)i*5.0e-6;

				sprintf( var, "C:/Users/PCNO1/Desktop/%duK_rand_maxbol.csv", i*5 );						//CHANGE THE ADDRESS WHERE YOU WANT TO GENERATE CSV FILES!!!
				if( (rand_no=fopen(var,"w"))==NULL ){
					printf("'%s'could not be opened in writing mode.", var);
					return 0;
				}

				for(int i=1; i<=SAMPLE; i++){
					fprintf( rand_no,"%.16e,%.16e,%.16e\n", max_boltz_velo_x(broadd,maxx,temp), max_boltz_velo_y(broadd,maxx,temp), max_boltz_velo_z(broadd,maxx,temp) );
				}

				fflush(rand_no); fclose(rand_no);
			}
		}
	//====generate random no. based on Maxwell-Boltzmann for each temperature (5-100[uK])====


	//====generate random no. based on Maxwell-Boltzmann for 10[uK]====
		if( select==2 ){
			for(int j=1; j<=100; j++){
				for(int i=1; i<=1000; i++){
					sprintf( var, "G:/graduation/rand_maxbol_10uK_2/10uK_rand_maxbol%d_%d.csv", j, i );	//CHANGE THE ADDRESS WHERE YOU WANT TO GENERATE CSV FILES!!!
					if( (rand_no=fopen(var,"w"))==NULL ){
						printf("'%s'could not be opened in writing mode.", var);
						return 0;
					}

					for(int i=1; i<=SAMPLE; i++){
						fprintf( rand_no,"%.16e,%.16e,%.16e\n", max_boltz_velo_x(0.15,13.0,10.0e-6), max_boltz_velo_y(0.15,13.0,10.0e-6), max_boltz_velo_z(0.15,13.0,10.0e-6) );
					}

					fflush(rand_no); fclose(rand_no);
				}
			}
		}
	//====generate random no. based on Maxwell-Boltzmann for 10[uK]====


	return 0;
}

	//Maxwell-Boltzmann distribution of velocity (x or y or z component)
	//This function is defined only to simplify the definition of the function max_boltz_velo_x(y,z).
	double max_boltz(double x, double t){	//temp should be substituted into t
		return sqrt( mass/(2.0*M_PI*k_b*t) )*exp( -mass*x*x/(2.0*k_b*t) );
	}

	//random velocity (x component) generator based on Maxwell-Boltzmann distribution
	double max_boltz_velo_x(double x, double m, double t){ //if x=0.15, then random no. between [-0.15[m/s],0.15[m/s]] are generated
		double xi,yj;
		do{
				xi = 2.0*x*rand()/RAND_MAX-x;//generate uniform random no. xi between [-x,x] (x=n*wr/2)
				yj = m*rand()/RAND_MAX;//12.8 is a value that is a little bit greater than the max of max_boltz()
			}while( max_boltz(xi, t)<yj );
		return xi;
	}

	//random velocity (y component) generator based on Maxwell-Boltzmann distribution
	double max_boltz_velo_y(double x, double m, double t){ //if x=0.15, then random no. between [-0.15[m/s],0.15[m/s]] are generated
		double xi,yj;
		do{
				xi = 2.0*x*rand()/RAND_MAX-x;//generate uniform random no. xi between [-x,x] (x=n*wr/2)
				yj = m*rand()/RAND_MAX;//12.8 is a value that is a little bit greater than the max of max_boltz()
			}while( max_boltz(xi, t)<yj );
		return xi;
	}

	//random velocity (z component) generator based on Maxwell-Boltzmann distribution
	double max_boltz_velo_z(double x, double m, double t){ //if x=0.15, then random no. between [-0.15[m/s],0.15[m/s]] are generated
		double xi,yj;
		do{
				xi = 2.0*x*rand()/RAND_MAX-x;//generate uniform random no. xi between [-x,x] (x=n*wr/2)
				yj = m*rand()/RAND_MAX;//12.8 is a value that is a little bit greater than the max of max_boltz()
			}while( max_boltz(xi, t)<yj );
		return xi;
	}
