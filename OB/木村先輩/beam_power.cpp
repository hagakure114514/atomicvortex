#include <iostream>
#include <math.h>
#include <cstdio>
#include <cstdlib>	//necessary when using rand(),srand(),EXIT_FAILURE etc.
#include <time.h>	//necessary when calling the current time


//===ATTENTION!====================================================================================================================================================
	//このプログラムの詳細を知りたい方は、木村の学士論文や研究室PCの木村のフォルダを参照してください。test

	//何か質問があれば、sb103_pic_sh_115og@yahoo.co.jpへメールを送ってください。

	//### 質量(Rb)や温度の値を変更した場合は、max_boltz_velo_x(y,z)のyjで使用している数値を変更する。

	//### MOTの半径を変更した場合は、gauss_pos_x(y,z)のyjの数値を変更し、ガウスランダムNo.も更新する。

	//### ビームの種類を変更した場合、損失判定の条件を変更する。初期値を変更した場合、損失判定の条件を変更する。

	//### 初期固有状態を変更する場合は、sta1またはsta2の初期値を変更する。


	// このファイルをコンパイルする際には、'gnuplot.exe', 'rand_gauss.csv', 'rand_maxbol.csv' のパスが正しいかどうか確認して下さい。

	// 注意 : "1/2"と書くと、整数の除算として計算され、答えは 0 になります。
	// 整数の計算の場合、答えの小数点以下の数字は無視されます。
	// しかし、"1.0/2.0"と書くと、double型の数値の割り算として計算されるので、答えは0.500000になります。

	// "#define" を使えば RAM 上のデータ量を節約することができます。

	// gnuplot.exe が時々何らかの理由で停止することがありますが、具体的な原因は分かりませんでした。


//===ATTENTION! <end>==============================================================================================================================================



//===SHOULD BE MODIFIED!===========================================================================================================================================
	//Is an atom really recoiled by a laser photon when a sp emission occurs?
//===SHOULD BE MODIFIED! <end>=====================================================================================================================================



//===Constants that CANNOT be changed==============================================================================================================================
	#define GNUPLOT_PATH "C:/PROGRA~2/gnuplot/bin/gnuplot.exe"
	//#define GNUPLOT_PATH "C:/Users/itolab/Desktop/gnuplot/gnuplot/binary/pgnuplot.exe"	//when executing this program on the labo PC
	#define hbar 1.0545718e-34		//Dirac constant[Js]
	#define h 6.626070e-34			//Planck constant[Js]
	#define c 2.99792458e+8			//light speed in vacuum[m/s]
	#define Rb 86.90918053			//atomic mass of 87Rb[g/mol]
	#define NA 6.022140e+23			//Avogadro constant[/mol]
	#define mass Rb/NA/1000.0		//mass of Rb : Rb/NA/1000[kg]
	#define k_b 1.380649e-23		//Boltzmann const.[J/K]
	#define G 9.80665				//standard gravity[m/s^2]
//===Constants that CANNOT be changed <end>========================================================================================================================


//===Constants that can be changed=================================================================================================================================
	//====variable====
	#define r0 1.0e-3												//MOT radius[m]
	#define waist 1.89e-3											//beam radius[m]
	#define r_out 1.50*waist										//outer radius of beam[m]
	#define r_simu 2.0*waist										//maximum radius up to which an atom is traced (2.0 is an arbitrarily chosen number.)
	#define r_judge waist/sqrt(2.0)									//radius to judge whther an atom is TRAPPED or LOSS (the radius of the peak of TEM)
	#define delta 1.0e+9*2.0*M_PI									//detuning[rad/s]
	#define branch 0.75												//branching ratio into |1>(lower hyperfine ground state) (based on "Gravitational laser trap for atoms with evanescent-wave cooling" pp.654, 661)

	//====fixed====
		#define SAMPLE 1000											//number of sample atoms[-]
		#define jloop 2000000										//how many times the time t is advanced[-]
		#define jjloop jloop/100									//the length of a uniform random no. array
		#define dt 5.0e-5											//interval time[s]
		#define temp 10.0e-6										//temperature of MOT[K]
		#define delta_hfs 6.834682610904e+9*2.0*M_PI				//frequency between the two hyperfine ground states of 87Rb[rad/s]
		#define gamma 38.1e+6										//natural linewidth of D2 line[rad/s]
		#define gamma1 gamma*branch									//spontaneous decay rate between |e> and |1>[rad/s]
		#define gamma2 gamma*(1.0-branch)							//spontaneous decay rate between |e> and |2>[rad/s]
		#define sigma3 1.0e-3										//triple standard deviation of density distribution of MOT[m]
		#define sigma sigma3/3.0									//standard deviation of density distribution of MOT[m]
	//====fixed====
	
	//############# temporary sentence for power loop ########################
		#define min_power 0.0e-3									//minimum value of when varying beam power[W]
		#define max_power 1005.0e-3									//maximum value of when varying beam power[W]
		#define d_power 10.0e-3										//increment of beam power[W]
	//############# temporary sentence for power loop <end> ##################

	//====depending on variables above====
		#define lambda c/(384.2e+12+delta/(2.0*M_PI))				//wavelength of detuned D2 line[m]
		#define lambda1 c/(384.2e+12)								//wavelength between |e> and |g1>[m]
		#define lambda2 c/(384.2e+12)								//wavelength between |e> and |g2>[m]
		#define k_wave 2.0*M_PI/lambda								//wavenumber[/m]
		#define I_s1 M_PI*h*c*gamma1/(3.0*lambda1*lambda1*lambda1)	//saturation intensity (|1>)[W/m^2]
		#define I_s2 M_PI*h*c*gamma2/(3.0*lambda2*lambda2*lambda2)	//saturation intensity (|2>)[W/m^2]
		double s1 = 0.0;											//saturation parameter between |g1> and |e>
		double s2 = 0.0;											//saturation parameter between |g2> and |e>
	//====depending on variables above====
//===Constants that can be changed <end>===========================================================================================================================


//===Variables for statements, temporary variables=================================================================================================================
	int gnu = 0;						//whether plotting the atomic positions with GNUPLOT
	int lump = 0;						//whether plotting the atomic positions in a lump
	int D_2_3 = 0;						//decide 2D plot or 3D plot from the input on command prompt
	int xyz = 0;						//whether plotting the atomic positions on x-y plane or x-z plane in 2D
	int sp = 0;							//whther a spontaneous emission is considered
	int csv_out = 0;					//whether output variables changing every moment into csv files
	int rep_sup = 0;					//the way of 2D plotting (superposition or replacemanet)
	int sta1 = 0;						//how many times an atom has been in |1>
	int sta2 = 0;						//how many times an atom has been in |2>
	int spon = 0;						//how many times spontaneous emissions have occured
	#define i_state 2					//initial energy eigenstate
	#define plots 0.1					//the size of dots plotted by GNUPLOT
	#define type 6						//the type of a marker plotted by GNUPLOT
	#define xlim_u 0.002				//upper limit of x coord. when plotting by GNUPLOT
	#define xlim_l -0.002				//lower limit of x coord. when plotting by GNUPLOT
	#define ylim_u 0.002				//upper limit of y coord. when plotting by GNUPLOT
	#define ylim_l -0.002				//lower limit of y coord. when plotting by GNUPLOT
	#define zlim_u 0.050				//upper limit of z coord. when plotting by GNUPLOT
	#define zlim_l -0.260				//lower limit of z coord. when plotting by GNUPLOT
	#define intv 100					//interval of k (=plotting)
	char var[20];						//values of variables for file names
	//############# temporary sentence for power loop ########################
	int vara = 0;						//values of variables for file names
	//############# temporary sentence for power loop <end> ##################
	double *xx;							//pointer for an array for gnuplot plotting (x coord.)
	double *yy;							//pointer for an array for gnuplot plotting (y coord.)
	double *zz;							//pointer for an array for gnuplot plotting (z coord.)
	double *tt;							//pointer for an array for gnuplot plotting (time)
	double *r_xx;						//pointer for an array for initial random position (x coord.)
	double *r_yy;						//pointer for an array for initial random position (y coord.)
	double *r_zz;						//pointer for an array for initial random position (z coord.)
	double *v_xx;						//pointer for an array for initial random velocity (x coord.)
	double *v_yy;						//pointer for an array for initial random velocity (y coord.)
	double *v_zz;						//pointer for an array for initial random velocity (z coord.)
	double *r1;							//pointer for an array for uniform random no.
	double *r2;							//pointer for an array for uniform random no.
	double *r3;							//pointer for an array for uniform random no.
	double *r4;							//pointer for an array for uniform random no.
	double *r5;							//pointer for an array for uniform random no.
	int count = 1;						//how many conditions have been executed
	char fn_rand_gauss[60];				//file name of Gaussian random no.
	char fn_rand_maxbol[60];			//file name of Maxwell-Boltzmann random no.

	//====file pointer====
		//file pointer for gnuplot.exe
		FILE *gp;						//file pointer for gnuplot.exe (The file address of gnuplot.exe is substituted into gp.)

		//file pointers for csv files
		FILE *i_velo;					//initial velocities
		FILE *i_velo_on_xy;				//initial velocities on xy plane
		FILE *i_velo_x_y_z;				//initial velocities in x, y and z direction
		FILE *i_velo_r_phi;				//initial radial and azimuthal velocities
		FILE *i_x_y_z;					//initial x, y and z coord.
		FILE *i_angmomentum;			//initial angular momentum
		FILE *i_kinetic;				//initial kinetic energy
		FILE *x_y_z_coord;				//x, y and z coordinates
		FILE *r_phi_coord;				//radii and azimuths
		FILE *velo;						//velocities
		FILE *x_y_z_velo;				//velocities in x, y and z direction
		FILE *r_max;					//maximum radii
		FILE *r_phi_velo;				//velocities in radial and azimuthal direction
		FILE *momentum;					//momentums
		FILE *angmomentum;				//angular momentums
		FILE *dipole;					//dipole forces
		FILE *passed_time;				//passed times till an atom reaches the substrate
		FILE *eigenstate;				//eigenstates
		FILE *kinetic;					//kinetic energies
		FILE *Udip;						//optical potential
		FILE *Ugrav;					//gravitational potential
		FILE *Energy;					//kinetic + potential
		FILE *Mom_x_y_z;				//momentum in x, y and z direction
		FILE *dEnergy;					//difference of kinetic + potential
		FILE *dMomentum;				//difference of momentum
		FILE *dMom_x_y_z;				//difference of momentum in x, y and z direction
		FILE *f_position;				//final x, y and z coord.
		FILE *f_state;					//final state (out of beam, substrate or still dropping)
		FILE *f_velo;					//final velocities for each direction
		FILE *f_spon;					//how many times spontaneous emission occured
		FILE *psp;						//transition probability
		FILE *rand_no1;					//initial random position
		FILE *rand_no2;					//initial random velocity
	//====file pointer====
//===Variables for statements, temporary variables <end>===========================================================================================================


//===Function Definition (prototype)===============================================================================================================================
	//if you want to know the details of the functions, go to the last part of this program.
	double Uniform_A(void);
	double Uniform_B1(void);
	double Uniform_B2(void);
	double Uniform_B3(void);
	double Uniform_B4(void);
	double Uniform_B5(void);
	double gauss(double x);
	double gauss_pos_x(double x);
	double gauss_pos_y(double x);
	double gauss_pos_z(double x);
	double max_boltz(double x, double t);
	double max_boltz_velo_x(double x);
	double max_boltz_velo_y(double x);
	double max_boltz_velo_z(double x);
	double int_TEM(double x, double y, double i, double w);
	double d_int_TEM(double x, double y, double i, double w);
	double U1_opt(double z, double d);
	double U2_opt(double z, double d);
	double dipole_f1(double z, double w, double d);
	double dipole_f2(double z, double w, double d);
	double dipole_f_x1(double x, double y, double z, double w);
	double dipole_f_y1(double x, double y, double z, double w);
	double dipole_f_x2(double x, double y, double z, double w);
	double dipole_f_y2(double x, double y, double z, double w);
//===Function Definition (prototype) <end>=========================================================================================================================







//===MAIN PROGRAM==================================================================================================================================================

int main(){

	srand((unsigned)time(NULL));	//a magic phrase to make random numbers be based on time

	//===ask...====================================================================================================================================================
		//===whether using GNUPLOT===========================================================================================
			printf("If you (don't) want to plot the atomic positions with GNUPLOT, then input 1(0).\n");
			scanf("%d", &gnu);
			printf("Please wait...\n");
		//===whether using GNUPLOT <end>=====================================================================================

		//===whether plotting in a lump======================================================================================
			if( gnu==1 ){
				printf("If you (don't) want to plot the atomic positions in a lump, input 1(0).\n");
				scanf("%d", &lump);
				printf("Please wait...\n");
			}
		//===whether plotting in a lump <end>================================================================================

		//===2D plot or 3D plot==============================================================================================
			if( gnu==1 ){
				printf("If you want to plot the atomic positions in 2D, 3D or both, input 2, 3 or 4, respectively.\n");
				scanf("%d", &D_2_3);
				printf("Please wait...\n");
			}
		//===2D plot or 3D plot <end>========================================================================================

		//===2D(x-y) plot or 2D(x-z) plot====================================================================================
			if( gnu==1 && (D_2_3==2 || D_2_3==4) ){
				printf("If you want to plot the atomic positions on x-y plane, x-z plane or both in 2D, input 1, 2 or 3, respectively.\n");
				scanf("%d", &xyz);
				printf("Please wait...\n");
			}
		//===2D(x-y) plot or 2D(x-z) plot <end>==============================================================================

		//===the way of plotting (replacement or superposition)==============================================================
			if( gnu==1 && (D_2_3==2 || D_2_3==4) && lump==0 ){
				printf("If you select replacemanet (superposition) as the way of 2D plotting, input 1(2).\n");
				scanf("%d", &rep_sup);
				printf("Please wait...\n");
			}
		//===the way of plotting (replacement or superposition) <end>========================================================

		//===whther a spontaneous emission is considered=====================================================================
			printf("If you (don't) want to take spontaneous emission into account, input 1(0).\n");
			scanf("%d", &sp);
			printf("Please wait...\n");
		//===whther a spontaneous emission is considered <end>===============================================================

		//===whther output variables changing every moment into csv files====================================================
			printf("If you (don't) want to output variables changing every moment into csv files, input 1(0).\n");
			scanf("%d", &csv_out);
			printf("Please wait...\n");
		//===whther output variables changing every moment into csv files <end>==============================================
	//===ask... <end>==============================================================================================================================================


	//===def of variables==========================================================================================================================================
		double x = 0.0, y = 0.0, z = 0.0;						//position of an atom[m]
		double r = 0.0, theta = 0.0, phi = 0.0;					//spherical coord. only to make spherical random no. for initial positions
		double theta2 = 0.0, phi2 = 0.0;						//spherical coord. only to make spherical random no. to determine the direction of the photon emitted when a spontaneous emission occurs
		double rad = 0.0;										//radial coord. in cylindrical coord. system
		double vx = 0.0, vy = 0.0, vz = 0.0, v = 0.0;			//initial velocity of an atom[m/s]
		double vr = 0.0, vphi = 0.0;							//initial cylindrical velocity of an atom[m/s]
		double a_dx = 0.0;										//acceleration by dipole force (x component)[m/s^2]
		double a_dy = 0.0;										//acceleration by dipole force (y component)[m/s^2]
		double t = 0.0;											//time[s]
		double kinetic_E = 0.0;									//kinetic energy[J]
		double mom = 0.0, ang_mom = 0.0;						//momentum[kg m/s], angular momentum[kg m^2/s]
		double rmax = 0.0;										//maximum radii[m]
		double dipole_f = 0.0;									//dipole force[N]
		double U_dip = 0.0;										//optical potential[J]
		double Ene = 0.0;										//kinetic + potential
		double Mom_x = 0.0;										//momentum in x direction
		double Mom_y = 0.0;										//momentum in y direction
		double Mom_z = 0.0;										//momentum in z direction
		double dEne = 0.0;										//difference of kinetic + potential
		double dMom = 0.0;										//difference of momentum
		double dMom_x = 0.0;									//difference of momentum in x direction
		double dMom_y = 0.0;									//difference of momentum in y direction
		double dMom_z = 0.0;									//difference of momentum in z direction
		double psp1 = 0.0;										//probability of spontaneous emission from |1>
		double psp2 = 0.0;										//probability of spontaneous emission from |2>
		double power = 0.0;										//beam power[W]
		int state = 0;											//energy eigenstate

		//temporary variables
		double x_t = 0.0, y_t = 0.0, z_t = 0.0;					//position of an atom
		double phi_t = 0.0;										//spherical coord. only to make spherical random no.
		double rad_t = 0.0;										//radial coord. in cylindrical coord. system
		double vx_t = 0.0, vy_t = 0.0, vz_t = 0.0, v_t = 0.0;	//temporary velocity of an atom[m/s]
		double vx_tt = 0.0, vy_tt = 0.0, vz_tt = 0.0, v_tt = 0.0;//temporary velocity of an atom[m/s]
		double vr_t = 0.0, vphi_t = 0.0;						//initial cylindrical velocity of an atom[m/s]
		double t_t = 0.0;										//time[s]
		int state_t = 0;										//eigenstate
		double ran1 = 0.0, ran2 = 0.0, ran3 = 0.0;				//random no.
	//===def of variables <end>====================================================================================================================================


	//===for loop (repeat changing one variable)===================================================================================================================
		for( power = min_power; power <=max_power ; power += d_power, vara += 10, count++ ){

			//############# temporary sentence for power loop ##################
				double I_0 = power/(M_PI*waist*waist/2.0);	//[W/m^2]
				sprintf( var, "%dmW", vara );
			//############# temporary sentence for power loop <end>#############

			//===generate random no. (initial positions)===========================================================================================================
				r_xx = new double[SAMPLE+1]{}; r_yy = new double[SAMPLE+1]{}; r_zz = new double[SAMPLE+1]{};
				sprintf( fn_rand_gauss, "../rand_gauss2/rand_gauss%d.csv", count );
				if( (rand_no1=fopen(fn_rand_gauss,"r"))==NULL ){
					printf("'%s' could not be opened in reading mode.", fn_rand_gauss);
					return 0;
				}
				for(int i=1; i<=SAMPLE; i++){
					fscanf( rand_no1, "%lE,%*lE,%*lE\n", &r_xx[i] ); //By using * ,one can omit unnecessary data.
				}
				fclose(rand_no1);
				
				if( (rand_no1=fopen(fn_rand_gauss,"r"))==NULL ){
					printf("'%s' could not be opened in reading mode.", fn_rand_gauss);
					return 0;
				}
				for(int i=1; i<=SAMPLE; i++){
					fscanf( rand_no1, "%*lE,%lE,%*lE\n", &r_yy[i] ); //By using * ,one can omit unnecessary data.
				}
				fclose(rand_no1);
				
				if( (rand_no1=fopen(fn_rand_gauss,"r"))==NULL ){
					printf("'%s' could not be opened in reading mode.", fn_rand_gauss);
					return 0;
				}
				for(int i=1; i<=SAMPLE; i++){
					fscanf( rand_no1, "%*lE,%*lE,%lE\n", &r_zz[i] ); //By using * ,one can omit unnecessary data.
				}
				fclose(rand_no1);
			//===generate random no. (initial positions) <end>=====================================================================================================

			//===generate random no. (initial velocities)==========================================================================================================
				v_xx = new double[SAMPLE+1]{}; v_yy = new double[SAMPLE+1]{}; v_zz = new double[SAMPLE+1]{};
				sprintf( fn_rand_maxbol, "../rand_maxbol_10uK/10uK_rand_maxbol%d.csv", count );
				if( (rand_no2=fopen(fn_rand_maxbol,"r"))==NULL ){
					printf("'%s' could not be opened in reading mode.", fn_rand_maxbol);
					return 0;
				}
				for(int i=1; i<=SAMPLE; i++){
					fscanf( rand_no2, "%lE,%*lE,%*lE\n", &v_xx[i] ); //By using * ,one can omit unnecessary data.
				}
				fclose(rand_no2);

				if( (rand_no2=fopen(fn_rand_maxbol,"r"))==NULL ){
					printf("'%s' could not be opened in reading mode.", fn_rand_maxbol);
					return 0;
				}
				for(int i=1; i<=SAMPLE; i++){
					fscanf( rand_no2, "%*lE,%lE,%*lE\n", &v_yy[i] ); //By using * ,one can omit unnecessary data.
				}
				fclose(rand_no2);
				
				if( (rand_no2=fopen(fn_rand_maxbol,"r"))==NULL ){
					printf("'%s' could not be opened in reading mode.", fn_rand_maxbol);
					return 0;
				}
				for(int i=1; i<=SAMPLE; i++){
					fscanf( rand_no2, "%*lE,%*lE,%lE\n", &v_zz[i] ); //By using * ,one can omit unnecessary data.
				}
				fclose(rand_no2);
			//===generate random no. (initial velocities) <end>====================================================================================================


			//===open files for outputting initial and final values================================================================================================
				//initial radii (on x-y plane)
				char fn_i_x_y_z[50];
				sprintf( fn_i_x_y_z, "%s_i_x_y_z.csv", var );
				if( (i_x_y_z=fopen(fn_i_x_y_z,"w"))==NULL ){
					printf("'%s' could not be opened in writing mode.", fn_i_x_y_z);
					return 0;
				}

				//initial velocities (on x-y plane)
				char fn_i_velo_on_xy[50];
				sprintf( fn_i_velo_on_xy, "%s_i_velocity_on_xy.csv", var );
				if( (i_velo_on_xy=fopen(fn_i_velo_on_xy,"w"))==NULL ){
					printf("'%s' could not be opened in writing mode.", fn_i_velo_on_xy);
					return 0;
				}

				//initial velocities in x, y and z direction
				char fn_i_velo_x_y_z[50];
				sprintf( fn_i_velo_x_y_z, "%s_i_velocity_x_y_z.csv", var );
				if( (i_velo_x_y_z=fopen(fn_i_velo_x_y_z,"w"))==NULL ){
					printf("'%s' could not be opened in writing mode.", fn_i_velo_x_y_z);
					return 0;
				}

				//initial velocities
				char fn_i_velo[50];
				sprintf( fn_i_velo, "%s_i_velocity.csv", var );
				if( (i_velo=fopen(fn_i_velo,"w"))==NULL ){
					printf("'%s' could not be opened in writing mode.", fn_i_velo);
					return 0;
				}

				//initial radial velocities
				char fn_i_velo_r_phi[50];
				sprintf( fn_i_velo_r_phi, "%s_i_velocity_r_phi.csv", var );
				if( (i_velo_r_phi=fopen(fn_i_velo_r_phi,"w"))==NULL ){
					printf("'%s' could not be opened in writing mode.", fn_i_velo_r_phi);
					return 0;
				}

				//maximum radii
				char fn_r_max[50];
				sprintf( fn_r_max, "%s_r_max.csv", var );
				if( (r_max=fopen(fn_r_max,"w"))==NULL ){
					printf("'%s' could not be opened in writing mode.", fn_r_max);
					return 0;
				}

				//final position
				char fn_f_position_rxyz[50];
				sprintf( fn_f_position_rxyz, "%s_f_position_r_x_y_z.csv", var );
				if( (f_position=fopen(fn_f_position_rxyz,"w"))==NULL ){
					printf("'%s' could not be opened in writing mode.", fn_f_position_rxyz);
					return 0;
				}

				//passed time
				char fn_passed_time[50];
				sprintf( fn_passed_time, "%s_passed_time.csv", var );
				if( (passed_time=fopen(fn_passed_time,"w"))==NULL ){
					printf("'%s' could not be opened in writing mode.", fn_passed_time);
					return 0;
				}

				//final state (out of beam, substrate or still dropping)
				char fn_f_state[50];
				sprintf( fn_f_state, "%s_f_state.csv", var );
				if( (f_state=fopen(fn_f_state,"w"))==NULL ){
					printf("'%s' could not be opened in writing mode.", fn_f_state);
					return 0;
				}

				//initial angular momentum
				char fn_i_angmomentum[50];
				sprintf( fn_i_angmomentum, "%s_i_ang_mom.csv", var );
				if( (i_angmomentum=fopen(fn_i_angmomentum,"w"))==NULL ){
					printf("'%s' could not be opened in writing mode.", fn_i_angmomentum );
					return 0;
				}

				//initial kinetic energy
				char fn_i_kinetic[50];
				sprintf( fn_i_kinetic, "%s_i_kinetic.csv", var );
				if( (i_kinetic=fopen(fn_i_kinetic,"w"))==NULL ){
					printf("'%s' could not be opened in writing mode.", fn_i_kinetic );
					return 0;
				}

				//final velocity in each direction
				char fn_f_velo[50];
				sprintf( fn_f_velo, "%s_f_velo_x_y_z.csv", var );
				if( (f_velo=fopen(fn_f_velo,"w"))==NULL ){
					printf("'%s' could not be opened in writing mode.", fn_f_velo );
					return 0;
				}

				//how many times spontaneous emission occured
				char fn_spon[50];
				sprintf( fn_spon, "%s_sp_emission.csv", var );
				if( (f_spon=fopen(fn_spon,"w"))==NULL ){
					printf("'%s' could not be opened in writing mode.", fn_spon );
					return 0;
				}
			//===open files for outputting initial and final values <end>==========================================================================================


			//===for loop (repeat the number of SAMPLES)===========================================================================================================
				for( int i=1; i<=SAMPLE; i++ ){	


				//===initialize variables and arrays with 0========================================================================================================
					state = 0;										//the current eigenstate (The author initialized "state" with "2" under the assumption of all of the atoms emitted from MOT being in |2>)
					x = 0.0, y = 0.0, z = 0.0;						//position of an atom[m]
					r = 0.0, theta = 0.0, phi = 0.0;				//spherical coord. only to make spherical random no. for initial positions
					theta2 = 0.0, phi2 = 0.0;						//spherical coord. only to make spherical random no. to determine the direction of the photon emitted when a spontaneous emission occurs
					rad = 0.0;										//cylindrical coord.
					vx = 0.0, vy = 0.0, vz = 0.0, v = 0.0;			//initial velocity of an atom[m/s]
					vr = 0.0, vphi = 0.0;							//initial cylindrical velocity of an atom[m/s]
					a_dx = 0.0;										//acceleration by dipole force (x component)[m/s^2]
					a_dy = 0.0;										//acceleration by dipole force (y component)[m/s^2]
					t = 0.0;										//time[s]
					kinetic_E = 0.0;								//kinetic energy[J]
					mom = 0.0, ang_mom = 0.0;						//momentum[kg m/s], angular momentum[kg m^2/s]
					Mom_x = 0.0; Mom_y = 0.0; Mom_z = 0.0;			//momentum in x, y and z direction[kg m/s]
					rmax = 0.0;										//maximum radii[m]
					dipole_f = 0.0;									//dipole force[N]
					U_dip = 0.0;									//optical potential[J]
					Ene = 0.0;										//kinetic + potential[J]
					s1 = 0.0;										//saturation parameter between |g1> and |e>[-]
					s2 = 0.0;										//saturation parameter between |g2> and |e>[-]
					psp1 = 0.0;										//probability of spontaneous emission from |1>
					psp2 = 0.0;										//probability of spontaneous emission from |2>
					sta1 = 0;										//how many times the atom has been in |1>
					sta2 = 0;										//how many times the atom has been in |2>
					spon = 0;										//how many times spontaneous emissions have occured

					//temporary variables
					x_t = 0.0, y_t = 0.0, z_t = 0.0;				//position of an atom
					phi_t = 0.0;									//spherical coord. only to make spherical random no.
					rad_t = 0.0;									//cylindrical coord.
					vx_t = 0.0, vy_t = 0.0, vz_t = 0.0, v_t = 0.0;	//temporary velocity of an atom[m/s]
					vx_tt = 0.0, vy_tt = 0.0, vz_tt = 0.0, v_tt = 0.0;//temporary velocity of an atom[m/s]
					vr_t = 0.0, vphi_t = 0.0;						//initial cylindrical velocity of an atom[m/s]
					t_t = 0.0;										//time[s]
					state_t = 0;									//eigenstate

					//arrays (for plotting)
					if(gnu==1){
						xx = new double[jloop+1]{};	yy = new double[jloop+1]{};	zz = new double[jloop+1]{};	tt = new double[jloop+1]{};
					}
				//===initialize variables and arrays with 0 <end>==================================================================================================

				//===generate random no. (uniform) ================================================================================================================
					//arrays (for Uniform)
					r1 = new double[jjloop+1]{}; r2 = new double[jjloop+1]{}; r3 = new double[jjloop+1]{}; r4 = new double[jjloop+1]{}; r5 = new double[jjloop+1]{};
					for(int k=1; k<=jjloop; k++){
						r1[k]=Uniform_B1(); r2[k]=Uniform_B2(); r3[k]=Uniform_B3(); r4[k]=Uniform_B4(); r5[k]=Uniform_B5();
					}
					//One can make appropriate random no. by using arrays.
					//If possible, r1-r5 should be generated jloop times, but it consumes lots of RAM and RAM may overflow.
					//So, the author generated jjloop (=jloop/100) times that is enough in the case where dt = 5.0e-5 ms.
				//===generate random no. (uniform) <end>===========================================================================================================

				//===initialize variables and arrays===============================================================================================================
					//energy eigenstate
					state = i_state;

					//position (spherical coord.)
					r = abs(r_xx[i]);
					if(r_xx[i]>=0) phi = r_yy[i];
					else if(r_xx[i]<0) phi = M_PI + r_yy[i];
					theta = r_zz[i];

					//position (Cartesian coord.)
					x = r*sin(theta)*cos(phi); y = r*sin(theta)*sin(phi); z = r*cos(theta);
					rad = sqrt(x*x+y*y);

					//velocity
					vx = v_xx[i]; vy = v_yy[i]; vz = v_zz[i];//initial speeds are given as values between [-100.0[m/s] 100.0[m/s]]
					v = sqrt( vx*vx + vy*vy + vz*vz );
					vr = vx*cos(phi)+vy*sin(phi); vphi = -vx*sin(phi)+vy*cos(phi);

					//kinetic energy
					kinetic_E = 1.0/2.0*mass*v*v;

					//momentum
					mom = mass*v; ang_mom = mass*vphi*r;

					//momentum in x, y and z direction
					Mom_x = mass*vx; Mom_y = mass*vy; Mom_z = mass*vz;

					//maximum radius
					rmax = rad;

					//dipole force
					if( state==1 ){
						dipole_f = dipole_f1(int_TEM(x,y,I_0,waist),d_int_TEM(x,y,I_0,waist),delta);
					}else if( state==2 ){
						dipole_f = dipole_f2(int_TEM(x,y,I_0,waist),d_int_TEM(x,y,I_0,waist),delta);
					}

					//optical potential
					if( state==1 ){
						U_dip = U1_opt(int_TEM(x,y,I_0,waist), delta);
					}else if( state==2 ){
						U_dip = U2_opt(int_TEM(x,y,I_0,waist), delta);
					}

					//kinetic + potential
					if( state==1 ){
						Ene = 1.0/2.0*mass*v*v + U1_opt(int_TEM(x,y,I_0,waist), delta) + mass*G*z;
					}else if( state==2 ){
						Ene = 1.0/2.0*mass*v*v + U2_opt(int_TEM(x,y,I_0,waist), delta) + mass*G*z;
					}

					//saturation parameter
					s1 = int_TEM(x,y,I_0,waist)/(I_s1*(1.0+4.0*delta*delta/(gamma1*gamma1)));
					s2 = int_TEM(x,y,I_0,waist)/(I_s2*(1.0+4.0*(delta+delta_hfs)*(delta+delta_hfs)/(gamma2*gamma2)));

					//probability of spontaneous emission
					psp1 = gamma1*s1/2.0*dt;
					psp2 = gamma2*s2/2.0*dt;

					//how many times the atom has been in |1>(|2>)
					if( state==1 ){
						sta1 = 1;
					}else if( state==2 ){
						sta2 = 1;
					}

					//====insert x, y and z coord. into arrays for plotting============================================
					if(gnu==1){
						xx[1] = x; yy[1] = y; zz[1] = z; tt[1] = t;
					}

					//====temporary variables for f_position.csv, passed_time.csv, f_state.csv, f_velo.csv.====
					//position
					rad_t = rad; phi_t = phi;
					x_t = x; y_t = y; z_t = z;

					//velocity
					vx_t = vx; vy_t = vy; vz_t = vz;
					v_t = v;
					vr_t = vr; vphi_t = vphi;

					//eigenstate
					state_t = state;
				//===initialize variables and arrays <end>=========================================================================================================
					

				//===open files for outputting values==============================================================================================================
					if( csv_out==1 ){
						//x coordinates
						char fn_x_y_z_coord[50];
						sprintf( fn_x_y_z_coord, "%s_x_y_z_%d.csv", var, i );
						if( (x_y_z_coord=fopen(fn_x_y_z_coord,"w"))==NULL ){
							printf("'%s' could not be opened in writing mode.", fn_x_y_z_coord);
							return 0;
						}

						//velocities
						char fn_velo[50];
						sprintf( fn_velo, "%s_velocity_%d.csv", var, i );
						if( (velo=fopen(fn_velo,"w"))==NULL ){
							printf("'%s' could not be opened in writing mode.", fn_velo);
							return 0;
						}

						//x, y and z velocities
						char fn_x_y_z_velo[50];
						sprintf( fn_x_y_z_velo, "%s_velocity_x_y_z_%d.csv", var, i );
						if( (x_y_z_velo=fopen(fn_x_y_z_velo,"w"))==NULL ){
							printf("'%s' could not be opened in writing mode.", fn_x_y_z_velo);
							return 0;
						}

						//radius
						char fn_r_phi_coord[50];
						sprintf( fn_r_phi_coord, "%s_r_phi_%d.csv", var, i );
						if( (r_phi_coord=fopen(fn_r_phi_coord,"w"))==NULL ){
							printf("'%s' could not be opened in writing mode.", fn_r_phi_coord);
							return 0;
						}

						//r velocity
						char fn_r_phi_velo[50];
						sprintf( fn_r_phi_velo, "%s_velocity_r_phi_%d.csv", var, i );
						if( (r_phi_velo=fopen(fn_r_phi_velo,"w"))==NULL ){
							printf("'%s' could not be opened in writing mode.", fn_r_phi_velo);
							return 0;
						}

						//kinetic energy
						char fn_kinetic[50];
						sprintf( fn_kinetic, "%s_kinetic_%d.csv", var, i );
						if( (kinetic=fopen(fn_kinetic,"w"))==NULL ){
							printf("'%s' could not be opened in writing mode.", fn_kinetic);
							return 0;
						}

						//momentum
						char fn_momentum[50];
						sprintf( fn_momentum, "%s_mom_%d.csv", var, i );
						if( (momentum=fopen(fn_momentum,"w"))==NULL ){
							printf("'%s' could not be opened in writing mode.", fn_momentum);
							return 0;
						}

						//momentum in x, y and z direction
						char fn_momentum_x_y_z[50];
						sprintf( fn_momentum_x_y_z, "%s_mom_x_y_z_%d.csv", var, i );
						if( (Mom_x_y_z=fopen(fn_momentum_x_y_z,"w"))==NULL ){
							printf("'%s' could not be opened in writing mode.", fn_momentum_x_y_z);
							return 0;
						}

						//angular momentum
						char fn_ang_momentum[50];
						sprintf( fn_ang_momentum, "%s_ang_mom_%d.csv", var, i );
						if( (angmomentum=fopen(fn_ang_momentum,"w"))==NULL ){
							printf("'%s' could not be opened in writing mode.", fn_ang_momentum);
							return 0;
						}

						//dipole force
						char fn_dipole_f[50];
						sprintf( fn_dipole_f, "%s_dipole_f_%d.csv", var, i );
						if( (dipole=fopen(fn_dipole_f,"w"))==NULL ){
							printf("'%s' could not be opened in writing mode.", fn_dipole_f);
							return 0;
						}

						//optical potential
						char fn_opt_pot[50];
						sprintf( fn_opt_pot, "%s_opt_pot_%d.csv", var, i );
						if( (Udip=fopen(fn_opt_pot,"w"))==NULL ){
							printf("'%s' could not be opened in writing mode.", fn_opt_pot);
							return 0;
						}

						//gravitational potential
						char fn_grav_pot[50];
						sprintf( fn_grav_pot, "%s_grav_pot_%d.csv", var, i );
						if( (Ugrav=fopen(fn_grav_pot,"w"))==NULL ){
							printf("'%s' could not be opened in writing mode.", fn_grav_pot);
							return 0;
						}

						//eigenstate
						char fn_eigstate[50];
						sprintf( fn_eigstate, "%s_eigstate_%d.csv", var, i );
						if( (eigenstate=fopen(fn_eigstate,"w"))==NULL ){
							printf("'%s' could not be opened in writing mode.", fn_eigstate);
							return 0;
						}

						//kinetic + potential
						char fn_Energy[50];
						sprintf( fn_Energy, "%s_Energy_%d.csv", var, i );
						if( (Energy=fopen(fn_Energy,"w"))==NULL ){
							printf("'%s' could not be opened in writing mode.", fn_Energy);
							return 0;
						}

						//transition probability
						char fn_psp[50];
						sprintf( fn_psp, "%s_psp_%d.csv", var, i );
						if( (psp=fopen(fn_psp,"w"))==NULL ){
							printf("'%s' could not be opened in writing mode.", fn_psp);
							return 0;
						}

						//difference of kinetic + potential
						char fn_dEnergy[50];
						sprintf( fn_dEnergy, "%s_dEnergy_%d.csv", var, i );
						if( (dEnergy=fopen(fn_dEnergy,"w"))==NULL ){
							printf("'%s' could not be opened in writing mode.", fn_dEnergy);
							return 0;
						}

						//difference of momentum
						char fn_dMom[50];
						sprintf( fn_dMom, "%s_dMom_%d.csv", var, i );
						if( (dMomentum=fopen(fn_dMom,"w"))==NULL ){
							printf("'%s' could not be opened in writing mode.", fn_dMom);
							return 0;
						}

						//difference of momentum in x, y and z direction
						char fn_dMomentum[50];
						sprintf( fn_dMomentum, "%s_dMom_x_y_z_%d.csv", var, i );
						if( (dMom_x_y_z=fopen(fn_dMomentum,"w"))==NULL ){
							printf("'%s' could not be opened in writing mode.", fn_dMomentum);
							return 0;
						}
					}
				//===open files for outputting values <end>========================================================================================================

					//===for loop (repeat the number of exporting times)===========================================================================================
						for( int j=1; j<=jloop; j++ ){		

							//===output initial values and variations into csv files===============================================================================
								if( j==1 ){
									fprintf( i_x_y_z, "%.16e,%.16e,%.16e\n", x, y, z );					//output initial x, y and z coord. into csv
									fprintf( i_velo_on_xy, "%.16e\n", sqrt(vx*vx+vy*vy) );				//output initial velocities (on x-y plane) into csv
									fprintf( i_velo_x_y_z, "%.16e,%.16e,%.16e\n", vx, vy, vz );			//output initial velocities in x, y and z direction into csv
									fprintf( i_velo, "%.16e\n", v );									//output initial velocities into csv
									fprintf( i_velo_r_phi, "%.16e,%.16e\n", vr, vphi );					//output initial radial and azimuthal velocities into csv
									fprintf( i_angmomentum, "%.16e\n", ang_mom );						//output initial angular momentum into csv
									fprintf( i_kinetic, "%.16e\n", kinetic_E );							//output initial kinetic energy into csv
									if( csv_out==1 ){
									fprintf( dEnergy, "%.16e\n", dEne );								//output difference of kinetic + potential energy into csv
									fprintf( dMomentum, "%.16e\n", dMom );								//output difference of momentum into csv
									fprintf( dMom_x_y_z, "%.16e,%.16e,%.16e\n", dMom_x, dMom_y, dMom_z );//output difference of momentum in x, y and z direction into csv
									}
								}
							//===output initial values and variations into csv files <end>=========================================================================


							//===output values into csv files======================================================================================================
								if( csv_out==1 ){
									fprintf( x_y_z_coord, "%.16e,%.16e,%.16e\n", x, y, z );				//outuput x coord. into csv
									fprintf( velo, "%.16e\n", v );										//outuput velocity into csv
									fprintf( x_y_z_velo, "%.16e,%.16e,%.16e\n", vx, vy, vz );			//outuput velocity in x, y and z direction into csv
									fprintf( r_phi_coord, "%.16e,%.16e\n", rad, phi );					//outuput radius and azimuth into csv
									fprintf( r_phi_velo, "%.16e,%.16e\n", vr, vphi );					//outuput radial and azimuthal velocity into csv
									fprintf( kinetic, "%.16e\n", kinetic_E );							//outuput kinetic energy into csv
									fprintf( momentum, "%.16e\n", mom );								//outuput momentum into csv
									fprintf( Mom_x_y_z, "%.16e,%.16e,%.16e\n", Mom_x, Mom_y, Mom_z );	//outuput momentum in x, y and z direction into csv
									fprintf( angmomentum, "%.16e\n", ang_mom );							//outuput angular momentum into csv
									fprintf( dipole, "%.16e\n", dipole_f );								//outuput dipole force into csv
									fprintf( Udip, "%.16e\n", U_dip );									//outuput optical potential into csv
									fprintf( Ugrav, "%.16e\n", mass*G*z );								//outuput gravitational potential into csv
									fprintf( Energy, "%.16e\n", Ene );									//outuput kinetic + potential energy into csv
									fprintf( eigenstate, "%d\n", state );								//outuput energy eigenstate into csv
									fprintf( psp, "%.16e,%.16e\n", psp1, psp2 );						//outuput transition probability into csv
								}
							//===output values into csv files <end>================================================================================================


							//===check the CURRENT STATE and branch into the CORRESPONDING STEP====================================================================
								if( state == 1 ){			//If "state" is "1", then the atom goes to the |1> step.

									//===whether spontaneous emission occurs from |1>==============================================================================
										if( sp==1 ){
											if( r1[sta1+1]<=psp1 && I_0!=0.0 ){
												//===advance velosity for dt based on absorption of laser photonand on spontaneous emission========================
													theta2 = M_PI*r3[spon+1]; phi2 = 2.0*M_PI*r4[spon+1];
													vx += hbar*k_wave*sin(theta2)*cos(phi2)/mass;
													vy += hbar*k_wave*sin(theta2)*sin(phi2)/mass;
													vz += hbar*k_wave/mass + hbar*k_wave*cos(theta2)/mass;
													vx_tt = vx; vy_tt = vy; vz_tt = vz; v_tt = sqrt(vx*vx+vy*vy+vz*vz);//substitute velocities into temporary variables after spontaneous emission
													spon++;
												//===advance velosity for dt based on absorption of laser photonand on spontaneous emission <end>==================

												//===branch into a eigenstate based on branching ratio=============================================================
													if( r5[spon+1]<=branch ){
														//sp emission into |1> (advance dt)
														//===advance position, velosity, time etc. for dt on |1>===================================================
															a_dx = dipole_f_x1(x,y,int_TEM(x,y,I_0,waist),d_int_TEM(x,y,I_0,waist))/mass;	//acceleration by dipole force (x component)
															a_dy = dipole_f_y1(x,y,int_TEM(x,y,I_0,waist),d_int_TEM(x,y,I_0,waist))/mass;	//acceleration by dipole force (y component)
															x += vx*dt + (1.0/2.0)*a_dx*dt*dt; y += vy*dt + (1.0/2.0)*a_dy*dt*dt; z += vz*dt - (1.0/2.0)*G*dt*dt;
															vx += a_dx*dt; vy += a_dy*dt; vz += - G*dt;
															v = sqrt( vx*vx + vy*vy + vz*vz );
															rad = sqrt( x*x + y*y );
															if( x==0.0 ){
																if( y==0.0 ) phi = phi;		//An azimuth need not be defined when r=0 (i.e. x=0 and y=0) but I defined it as it becomes continuous.
																else if( y>0.0 ) phi = M_PI/2.0;
																else if( y<0.0 ) phi = -M_PI/2.0;
															}else if( x>0.0 ){
																if( y>0.0 ) phi = atan( y/x );
																else if( y<0.0 ) phi = atan( y/x )+2.0*M_PI;
															}else if( x<0.0 ){
																if( y>0.0 ) phi = atan( y/x )+M_PI;
																else if( y<0.0 ) phi = atan( y/x )+M_PI;
															}
															vr = vx*cos(phi)+vy*sin(phi); vphi = -vx*sin(phi)+vy*cos(phi);
															kinetic_E = 1.0/2.0*mass*v*v;
															mom = mass*v; ang_mom = mass*vphi*r;
															Mom_x = mass*vx; Mom_y = mass*vy; Mom_z = mass*vz;
															if( rmax<rad ) rmax = rad;
															dipole_f = dipole_f1(int_TEM(x,y,I_0,waist), d_int_TEM(x,y,I_0,waist),delta);
															U_dip = U1_opt(int_TEM(x,y,I_0,waist), delta);
															Ene = 1.0/2.0*mass*v*v + U1_opt(int_TEM(x,y,I_0,waist), delta) + mass*G*z;
															s1 = int_TEM(x_t,y_t,I_0,waist)/(I_s1*(1.0+4.0*delta*delta/(gamma1*gamma1)));
															s2 = int_TEM(x_t,y_t,I_0,waist)/(I_s2*(1.0+4.0*(delta+delta_hfs)*(delta+delta_hfs)/(gamma2*gamma2)));
															psp1 = gamma1*s1/2.0*dt;
															psp2 = gamma2*s2/2.0*dt;
															t += dt;
															sta1++;

															//====update arrays====
															if(gnu==1){
																xx[j+1] = x; yy[j+1] = y; zz[j+1] = z; tt[j+1] = t;
															}
														//===advance position, velosity, time etc. for dt on |1> <end>=============================================
														state = 1;
													}else{
														//sp emission into |2> (advance dt)
														//===advance position, velosity, time etc. for dt on |2>===================================================
															a_dx = dipole_f_x2(x,y,int_TEM(x,y,I_0,waist),d_int_TEM(x,y,I_0,waist))/mass;	//acceleration by dipole force (x component)
															a_dy = dipole_f_y2(x,y,int_TEM(x,y,I_0,waist),d_int_TEM(x,y,I_0,waist))/mass;	//acceleration by dipole force (y component)
															x += vx*dt + (1.0/2.0)*a_dx*dt*dt; y += vy*dt + (1.0/2.0)*a_dy*dt*dt; z += vz*dt - (1.0/2.0)*G*dt*dt;
															vx += a_dx*dt; vy += a_dy*dt; vz += - G*dt;
															v = sqrt( vx*vx + vy*vy + vz*vz );
															rad = sqrt( x*x + y*y );
															if( x==0.0 ){
																if( y==0.0 ) phi = phi;		//An azimuth need not be defined when r=0 (i.e. x=0 and y=0) but I defined it as it becomes continuous.
																else if( y>0.0 ) phi = M_PI/2.0;
																else if( y<0.0 ) phi = -M_PI/2.0;
															}else if( x>0.0 ){
																if( y>0.0 ) phi = atan( y/x );
																else if( y<0.0 ) phi = atan( y/x )+2.0*M_PI;
															}else if( x<0.0 ){
																if( y>0.0 ) phi = atan( y/x )+M_PI;
																else if( y<0.0 ) phi = atan( y/x )+M_PI;
															}
															vr = vx*cos(phi)+vy*sin(phi); vphi = -vx*sin(phi)+vy*cos(phi);
															kinetic_E = 1.0/2.0*mass*v*v;
															mom = mass*v; ang_mom = mass*vphi*r;
															Mom_x = mass*vx; Mom_y = mass*vy; Mom_z = mass*vz;
															if( rmax<rad ) rmax = rad;
															dipole_f = dipole_f2(int_TEM(x,y,I_0,waist), d_int_TEM(x,y,I_0,waist),delta);
															U_dip = U2_opt(int_TEM(x,y,I_0,waist), delta);
															Ene = 1.0/2.0*mass*v*v + U2_opt(int_TEM(x,y,I_0,waist), delta) + mass*G*z;
															s1 = int_TEM(x_t,y_t,I_0,waist)/(I_s1*(1.0+4.0*delta*delta/(gamma1*gamma1)));
															s2 = int_TEM(x_t,y_t,I_0,waist)/(I_s2*(1.0+4.0*(delta+delta_hfs)*(delta+delta_hfs)/(gamma2*gamma2)));
															psp1 = gamma1*s1/2.0*dt;
															psp2 = gamma2*s2/2.0*dt;
															t += dt;
															sta2++;

															//====update arrays====
															if(gnu==1){
																xx[j+1] = x; yy[j+1] = y; zz[j+1] = z; tt[j+1] = t;
															}
														//===advance position, velosity, time etc. for dt on |2> <end>=============================================
														state = 2;
													}
												//===branch into a eigenstate based on branching ratio <end>=======================================================
											}else{
												//not sp emission (advance dt)
												//===advance position, velosity, time etc. for dt on |1>===========================================================
													a_dx = dipole_f_x1(x,y,int_TEM(x,y,I_0,waist),d_int_TEM(x,y,I_0,waist))/mass;	//acceleration by dipole force (x component)
													a_dy = dipole_f_y1(x,y,int_TEM(x,y,I_0,waist),d_int_TEM(x,y,I_0,waist))/mass;	//acceleration by dipole force (y component)
													x += vx*dt + (1.0/2.0)*a_dx*dt*dt; y += vy*dt + (1.0/2.0)*a_dy*dt*dt; z += vz*dt - (1.0/2.0)*G*dt*dt;
													vx += a_dx*dt; vy += a_dy*dt; vz += - G*dt;
													v = sqrt( vx*vx + vy*vy + vz*vz );
													rad = sqrt( x*x + y*y );
													if( x==0.0 ){
														if( y==0.0 ) phi = phi;		//An azimuth need not be defined when r=0 (i.e. x=0 and y=0) but I defined it as it becomes continuous.
														else if( y>0.0 ) phi = M_PI/2.0;
														else if( y<0.0 ) phi = -M_PI/2.0;
													}else if( x>0.0 ){
														if( y>0.0 ) phi = atan( y/x );
														else if( y<0.0 ) phi = atan( y/x )+2.0*M_PI;
													}else if( x<0.0 ){
														if( y>0.0 ) phi = atan( y/x )+M_PI;
														else if( y<0.0 ) phi = atan( y/x )+M_PI;
													}
													vr = vx*cos(phi)+vy*sin(phi); vphi = -vx*sin(phi)+vy*cos(phi);
													kinetic_E = 1.0/2.0*mass*v*v;
													mom = mass*v; ang_mom = mass*vphi*r;
													Mom_x = mass*vx; Mom_y = mass*vy; Mom_z = mass*vz;
													if( rmax<rad ) rmax = rad;
													dipole_f = dipole_f1(int_TEM(x,y,I_0,waist), d_int_TEM(x,y,I_0,waist),delta);
													U_dip = U1_opt(int_TEM(x,y,I_0,waist), delta);
													Ene = 1.0/2.0*mass*v*v + U1_opt(int_TEM(x,y,I_0,waist), delta) + mass*G*z;
													s1 = int_TEM(x_t,y_t,I_0,waist)/(I_s1*(1.0+4.0*delta*delta/(gamma1*gamma1)));
													s2 = int_TEM(x_t,y_t,I_0,waist)/(I_s2*(1.0+4.0*(delta+delta_hfs)*(delta+delta_hfs)/(gamma2*gamma2)));
													psp1 = gamma1*s1/2.0*dt;
													psp2 = gamma2*s2/2.0*dt;
													t += dt;
													sta1++;

													//====update arrays====
													if(gnu==1){
														xx[j+1] = x; yy[j+1] = y; zz[j+1] = z; tt[j+1] = t;
													}
												//===advance position, velosity, time etc. for dt on |1> <end>=====================================================
											}
										}else if( sp==0 ){
												//===advance position, velosity, time etc. for dt on |1>===========================================================
													a_dx = dipole_f_x1(x,y,int_TEM(x,y,I_0,waist),d_int_TEM(x,y,I_0,waist))/mass;	//acceleration by dipole force (x component)
													a_dy = dipole_f_y1(x,y,int_TEM(x,y,I_0,waist),d_int_TEM(x,y,I_0,waist))/mass;	//acceleration by dipole force (y component)
													x += vx*dt + (1.0/2.0)*a_dx*dt*dt; y += vy*dt + (1.0/2.0)*a_dy*dt*dt; z += vz*dt - (1.0/2.0)*G*dt*dt;
													vx += a_dx*dt; vy += a_dy*dt; vz += - G*dt;
													v = sqrt( vx*vx + vy*vy + vz*vz );
													rad = sqrt( x*x + y*y );
													if( x==0.0 ){
														if( y==0.0 ) phi = phi;		//An azimuth need not be defined when r=0 (i.e. x=0 and y=0) but I defined it as it becomes continuous.
														else if( y>0.0 ) phi = M_PI/2.0;
														else if( y<0.0 ) phi = -M_PI/2.0;
													}else if( x>0.0 ){
														if( y>0.0 ) phi = atan( y/x );
														else if( y<0.0 ) phi = atan( y/x )+2.0*M_PI;
													}else if( x<0.0 ){
														if( y>0.0 ) phi = atan( y/x )+M_PI;
														else if( y<0.0 ) phi = atan( y/x )+M_PI;
													}
													vr = vx*cos(phi)+vy*sin(phi); vphi = -vx*sin(phi)+vy*cos(phi);
													kinetic_E = 1.0/2.0*mass*v*v;
													mom = mass*v; ang_mom = mass*vphi*r;
													Mom_x = mass*vx; Mom_y = mass*vy; Mom_z = mass*vz;
													if( rmax<rad ) rmax = rad;
													dipole_f = dipole_f1(int_TEM(x,y,I_0,waist), d_int_TEM(x,y,I_0,waist),delta);
													U_dip = U1_opt(int_TEM(x,y,I_0,waist), delta);
													Ene = 1.0/2.0*mass*v*v + U1_opt(int_TEM(x,y,I_0,waist), delta) + mass*G*z;
													t += dt;

													//====update arrays====
													if(gnu==1){
														xx[j+1] = x; yy[j+1] = y; zz[j+1] = z; tt[j+1] = t;
													}
												//===advance position, velosity, time etc. for dt on |1> <end>=====================================================
										}
									//===whether spontaneous emission occurs from |1> <end>========================================================================

								}else if( state == 2 ){		//If "state" is "2", then the atom goes to the |2> step.

									//===whether spontaneous emission occurs from |2>==============================================================================
										if( sp==1 ){
											if( r2[sta2]<=psp2 && I_0!=0.0 ){
												//===advance velosity for dt based on absorption of laser photonand on spontaneous emission========================
													theta2 = M_PI*r3[spon+1]; phi2 = 2.0*M_PI*r4[spon+1];
													vx += hbar*k_wave*sin(theta2)*cos(phi2)/mass;
													vy += hbar*k_wave*sin(theta2)*sin(phi2)/mass;
													vz += hbar*k_wave/mass + hbar*k_wave*cos(theta2)/mass;
													vx_tt = vx; vy_tt = vy; vz_tt = vz; v_tt = sqrt(vx*vx+vy*vy+vz*vz);//substitute velocities into temporary variables after spontaneous emission
													spon++;
												//===advance velosity for dt based on absorption of laser photonand on spontaneous emission <end>==================

												//===branch into a eigenstate based on branching ratio=============================================================
													if( r5[spon+1]<=branch ){
														//sp emission into |1> (advance dt)
														//===advance position, velosity, time etc. for dt on |1>===================================================
															a_dx = dipole_f_x1(x,y,int_TEM(x,y,I_0,waist),d_int_TEM(x,y,I_0,waist))/mass;	//acceleration by dipole force (x component)
															a_dy = dipole_f_y1(x,y,int_TEM(x,y,I_0,waist),d_int_TEM(x,y,I_0,waist))/mass;	//acceleration by dipole force (y component)
															x += vx*dt + (1.0/2.0)*a_dx*dt*dt; y += vy*dt + (1.0/2.0)*a_dy*dt*dt; z += vz*dt - (1.0/2.0)*G*dt*dt;
															vx += a_dx*dt; vy += a_dy*dt; vz += - G*dt;
															v = sqrt( vx*vx + vy*vy + vz*vz );
															rad = sqrt( x*x + y*y );
															if( x==0.0 ){
																if( y==0.0 ) phi = phi;		//An azimuth need not be defined when r=0 (i.e. x=0 and y=0) but I defined it as it becomes continuous.
																else if( y>0.0 ) phi = M_PI/2.0;
																else if( y<0.0 ) phi = -M_PI/2.0;
															}else if( x>0.0 ){
																if( y>0.0 ) phi = atan( y/x );
																else if( y<0.0 ) phi = atan( y/x )+2.0*M_PI;
															}else if( x<0.0 ){
																if( y>0.0 ) phi = atan( y/x )+M_PI;
																else if( y<0.0 ) phi = atan( y/x )+M_PI;
															}
															vr = vx*cos(phi)+vy*sin(phi); vphi = -vx*sin(phi)+vy*cos(phi);
															kinetic_E = 1.0/2.0*mass*v*v;
															mom = mass*v; ang_mom = mass*vphi*r;
															Mom_x = mass*vx; Mom_y = mass*vy; Mom_z = mass*vz;
															if( rmax<rad ) rmax = rad;
															dipole_f = dipole_f1(int_TEM(x,y,I_0,waist), d_int_TEM(x,y,I_0,waist),delta);
															U_dip = U1_opt(int_TEM(x,y,I_0,waist), delta);
															Ene = 1.0/2.0*mass*v*v + U1_opt(int_TEM(x,y,I_0,waist), delta) + mass*G*z;
															s1 = int_TEM(x_t,y_t,I_0,waist)/(I_s1*(1.0+4.0*delta*delta/(gamma1*gamma1)));
															s2 = int_TEM(x_t,y_t,I_0,waist)/(I_s2*(1.0+4.0*(delta+delta_hfs)*(delta+delta_hfs)/(gamma2*gamma2)));
															psp1 = gamma1*s1/2.0*dt;
															psp2 = gamma2*s2/2.0*dt;
															t += dt;
															sta1++;

															//====update arrays====
															if(gnu==1){
																xx[j+1] = x; yy[j+1] = y; zz[j+1] = z; tt[j+1] = t;
															}
														//===advance position, velosity, time etc. for dt on |1> <end>=============================================
														state = 1;
													}else{
														//sp emission into |2> (advance dt)
														//===advance position, velosity, time etc. for dt on |2>===================================================
															a_dx = dipole_f_x2(x,y,int_TEM(x,y,I_0,waist),d_int_TEM(x,y,I_0,waist))/mass;	//acceleration by dipole force (x component)
															a_dy = dipole_f_y2(x,y,int_TEM(x,y,I_0,waist),d_int_TEM(x,y,I_0,waist))/mass;	//acceleration by dipole force (y component)
															x += vx*dt + (1.0/2.0)*a_dx*dt*dt; y += vy*dt + (1.0/2.0)*a_dy*dt*dt; z += vz*dt - (1.0/2.0)*G*dt*dt;
															vx += a_dx*dt; vy += a_dy*dt; vz += - G*dt;
															v = sqrt( vx*vx + vy*vy + vz*vz );
															rad = sqrt( x*x + y*y );
															if( x==0.0 ){
																if( y==0.0 ) phi = phi;		//An azimuth need not be defined when r=0 (i.e. x=0 and y=0) but I defined it as it becomes continuous.
																else if( y>0.0 ) phi = M_PI/2.0;
																else if( y<0.0 ) phi = -M_PI/2.0;
															}else if( x>0.0 ){
																if( y>0.0 ) phi = atan( y/x );
																else if( y<0.0 ) phi = atan( y/x )+2.0*M_PI;
															}else if( x<0.0 ){
																if( y>0.0 ) phi = atan( y/x )+M_PI;
																else if( y<0.0 ) phi = atan( y/x )+M_PI;
															}
															vr = vx*cos(phi)+vy*sin(phi); vphi = -vx*sin(phi)+vy*cos(phi);
															kinetic_E = 1.0/2.0*mass*v*v;
															mom = mass*v; ang_mom = mass*vphi*r;
															Mom_x = mass*vx; Mom_y = mass*vy; Mom_z = mass*vz;
															if( rmax<rad ) rmax = rad;
															dipole_f = dipole_f2(int_TEM(x,y,I_0,waist), d_int_TEM(x,y,I_0,waist),delta);
															U_dip = U2_opt(int_TEM(x,y,I_0,waist), delta);
															Ene = 1.0/2.0*mass*v*v + U2_opt(int_TEM(x,y,I_0,waist), delta) + mass*G*z;
															s1 = int_TEM(x_t,y_t,I_0,waist)/(I_s1*(1.0+4.0*delta*delta/(gamma1*gamma1)));
															s2 = int_TEM(x_t,y_t,I_0,waist)/(I_s2*(1.0+4.0*(delta+delta_hfs)*(delta+delta_hfs)/(gamma2*gamma2)));
															psp1 = gamma1*s1/2.0*dt;
															psp2 = gamma2*s2/2.0*dt;
															t += dt;
															sta2++;

															//====update arrays====
															if(gnu==1){
																xx[j+1] = x; yy[j+1] = y; zz[j+1] = z; tt[j+1] = t;
															}
														//===advance position, velosity, time etc. for dt on |2> <end>=============================================
														state = 2;
													}
												//===branch into a eigenstate based on branching ratio <end>=======================================================
											}else{
												//not sp emission (advance dt)
												//===advance position, velosity, time etc. for dt on |2>===========================================================
													a_dx = dipole_f_x2(x,y,int_TEM(x,y,I_0,waist),d_int_TEM(x,y,I_0,waist))/mass;	//acceleration by dipole force (x component)
													a_dy = dipole_f_y2(x,y,int_TEM(x,y,I_0,waist),d_int_TEM(x,y,I_0,waist))/mass;	//acceleration by dipole force (y component)
													x += vx*dt + (1.0/2.0)*a_dx*dt*dt; y += vy*dt + (1.0/2.0)*a_dy*dt*dt; z += vz*dt - (1.0/2.0)*G*dt*dt;
													vx += a_dx*dt; vy += a_dy*dt; vz += - G*dt;
													v = sqrt( vx*vx + vy*vy + vz*vz );
													rad = sqrt( x*x + y*y );
													if( x==0.0 ){
														if( y==0.0 ) phi = phi;		//An azimuth need not be defined when r=0 (i.e. x=0 and y=0) but I defined it as it becomes continuous.
														else if( y>0.0 ) phi = M_PI/2.0;
														else if( y<0.0 ) phi = -M_PI/2.0;
													}else if( x>0.0 ){
														if( y>0.0 ) phi = atan( y/x );
														else if( y<0.0 ) phi = atan( y/x )+2.0*M_PI;
													}else if( x<0.0 ){
														if( y>0.0 ) phi = atan( y/x )+M_PI;
														else if( y<0.0 ) phi = atan( y/x )+M_PI;
													}
													vr = vx*cos(phi)+vy*sin(phi); vphi = -vx*sin(phi)+vy*cos(phi);
													kinetic_E = 1.0/2.0*mass*v*v;
													mom = mass*v; ang_mom = mass*vphi*r;
													Mom_x = mass*vx; Mom_y = mass*vy; Mom_z = mass*vz;
													if( rmax<rad ) rmax = rad;
													dipole_f = dipole_f2(int_TEM(x,y,I_0,waist), d_int_TEM(x,y,I_0,waist),delta);
													U_dip = U2_opt(int_TEM(x,y,I_0,waist), delta);
													Ene = 1.0/2.0*mass*v*v + U2_opt(int_TEM(x,y,I_0,waist), delta) + mass*G*z;
													s1 = int_TEM(x_t,y_t,I_0,waist)/(I_s1*(1.0+4.0*delta*delta/(gamma1*gamma1)));
													s2 = int_TEM(x_t,y_t,I_0,waist)/(I_s2*(1.0+4.0*(delta+delta_hfs)*(delta+delta_hfs)/(gamma2*gamma2)));
													psp1 = gamma1*s1/2.0*dt;
													psp2 = gamma2*s2/2.0*dt;
													t += dt;
													sta2++;

													//====update arrays====
													if(gnu==1){
														xx[j+1] = x; yy[j+1] = y; zz[j+1] = z; tt[j+1] = t;
													}
												//===advance position, velosity, time etc. for dt on |2> <end>=====================================================
											}
										}else if( sp==0 ){
												//===advance position, velosity, time etc. for dt on |2>===========================================================
													a_dx = dipole_f_x2(x,y,int_TEM(x,y,I_0,waist),d_int_TEM(x,y,I_0,waist))/mass;	//acceleration by dipole force (x component)
													a_dy = dipole_f_y2(x,y,int_TEM(x,y,I_0,waist),d_int_TEM(x,y,I_0,waist))/mass;	//acceleration by dipole force (y component)
													x += vx*dt + (1.0/2.0)*a_dx*dt*dt; y += vy*dt + (1.0/2.0)*a_dy*dt*dt; z += vz*dt - (1.0/2.0)*G*dt*dt;
													vx += a_dx*dt; vy += a_dy*dt; vz += - G*dt;
													v = sqrt( vx*vx + vy*vy + vz*vz );
													rad = sqrt( x*x + y*y );
													if( x==0.0 ){
														if( y==0.0 ) phi = phi;		//An azimuth need not be defined when r=0 (i.e. x=0 and y=0) but I defined it as it becomes continuous.
														else if( y>0.0 ) phi = M_PI/2.0;
														else if( y<0.0 ) phi = -M_PI/2.0;
													}else if( x>0.0 ){
														if( y>0.0 ) phi = atan( y/x );
														else if( y<0.0 ) phi = atan( y/x )+2.0*M_PI;
													}else if( x<0.0 ){
														if( y>0.0 ) phi = atan( y/x )+M_PI;
														else if( y<0.0 ) phi = atan( y/x )+M_PI;
													}
													vr = vx*cos(phi)+vy*sin(phi); vphi = -vx*sin(phi)+vy*cos(phi);
													kinetic_E = 1.0/2.0*mass*v*v;
													mom = mass*v; ang_mom = mass*vphi*r;
													Mom_x = mass*vx; Mom_y = mass*vy; Mom_z = mass*vz;
													if( rmax<rad ) rmax = rad;
													dipole_f = dipole_f2(int_TEM(x,y,I_0,waist), d_int_TEM(x,y,I_0,waist),delta);
													U_dip = U2_opt(int_TEM(x,y,I_0,waist), delta);
													Ene = 1.0/2.0*mass*v*v + U2_opt(int_TEM(x,y,I_0,waist), delta) + mass*G*z;
													t += dt;

													//====update arrays====
													if(gnu==1){
														xx[j+1] = x; yy[j+1] = y; zz[j+1] = z; tt[j+1] = t;
													}
												//===advance position, velosity, time etc. for dt on |2> <end>=====================================================
										}
									//===whether spontaneous emission occurs from |2> <end>========================================================================
								}
							//===check the CURRENT STATE and branch into the CORRESPONDING STEP <end>==============================================================





							//===check whether the atom is "out of the beam", "on the substrate" or "still dropping" and plot the trajectory=======================
								if( sqrt(x*x+y*y) > r_simu || z < -0.26 || j==jloop ) {		//If the atom is out of beam or on the substrate or if j==jloop, this program quits to simulate the current atom.
									if( sqrt(x*x+y*y) > r_simu ){
										printf("Out of the beam\n");
									}else if( z < -0.26 ){
										printf("Substrate\n");
									}else if( j==jloop ){
										printf("Still dropping\n");
									}
									
									//===plot position=============================================================================================================
										if( gnu==1 && lump==0 ){
												
											if( ( gp = _popen(GNUPLOT_PATH, "w") ) == NULL ){							//open GNUPLOT by using "_popen" command
												fprintf( stderr, "The application could not be booted.%s\n", GNUPLOT_PATH );
												exit(EXIT_FAILURE);
											}
											
											//===plot (the way of plotting is replacement)==========================================================================
												if( rep_sup==1 ){
													if( (D_2_3 == 2 || D_2_3 == 4) && (xyz==1 || xyz==3) ){
														for( int k=1; k<=j; k+=intv ){
														fprintf( gp, "set terminal gif\n" );								//decide the exporting file type of graphs (here, the file type is .gif)
														fprintf( gp, "set grid\n" );										//draw grid in the graph
														fprintf( gp, "set title't = %f[s], z = %f[cm]'\n", tt[k], zz[k]*100.0 );//graph title
														fprintf( gp, "set xlabel 'x (m)'\n" );								//decide the x label
														fprintf( gp, "set ylabel 'y (m)'\n" );								//decide the y label
														fprintf( gp, "set xtics rotate by 45 center\n" );					//rotate scales 45 degrees
														fprintf( gp, "set xtics center offset 0,-1\n" );					//change the offset of scales on x axis
														fprintf( gp, "set xrange[%f:%f]\n", xlim_l, xlim_u );				//confine the length of x axis
														fprintf( gp, "set yrange[%f:%f]\n", ylim_l, ylim_u );				//confine the length of y axis
														fprintf( gp, "set size ratio -1\n" );								//equalize the ratio of x and y axes
														fprintf( gp, "set output '%s_pos2D_xy_%d_%d.gif'\n", var, i, k );	//decide the file name
														fprintf( gp, "plot \"-\" with points pt %d ps %f\n", type, plots );	//plot the atom position using a circle
														fprintf( gp, "%f, %f\n", xx[k], yy[k] );							//plot the atom position on x-z plane
														fprintf( gp, "e\n" );												//e means the end of inputs of coordinates
														fflush( gp );
														}
														if( D_2_3 == 2 && xyz==1 ){
															fprintf( gp, "exit\n");
															_pclose( gp );
														}
													}
													if( (D_2_3 == 2 || D_2_3 == 4) && (xyz==2 || xyz==3) ){
														for( int k=1; k<=j; k+=intv ){
														fprintf( gp, "set terminal gif\n" );								//decide the exporting file type of graphs (here, the file type is .gif)
														fprintf( gp, "set grid\n" );										//draw grid in the graph
														fprintf( gp, "set title't = %f[s], z = %f[cm]'\n", tt[k], zz[k]*100.0 );//graph title
														fprintf( gp, "set xlabel 'x (m)'\n" );								//decide the x label
														fprintf( gp, "set ylabel 'z (m)'\n" );								//decide the z label
														fprintf( gp, "set xtics rotate by 45 center\n" );					//rotate scales 45 degrees
														fprintf( gp, "set xtics center offset 0,-1\n" );					//change the offset of scales on x axis
														fprintf( gp, "set xrange[%f:%f]\n", xlim_l, xlim_u );				//confine the length of x axis
														fprintf( gp, "set yrange[%f:%f]\n", zlim_l, zlim_u );				//confine the length of z axis
														fprintf( gp, "set size noratio\n" );								//cancel of equalization of the ratio of axes
														fprintf( gp, "set output '%s_pos2D_xz_%d_%d.gif'\n", var, i, k );	//decide the file name
														fprintf( gp, "plot \"-\" with points pt %d ps %f\n", type, plots );	//plot the atom position using a circle
														fprintf( gp, "%f, %f\n", xx[k], zz[k] );							//plot the atom position on x-z plane
														fprintf( gp, "e\n" );												//e means the end of inputs of coordinates
														fflush( gp );
														}
														if( D_2_3 == 2 ){
															fprintf( gp, "exit\n");
															_pclose( gp );
														}
													}
													if( D_2_3 == 3 || D_2_3 == 4 ){
														for( int k=1; k<=j; k+=intv ){
														fprintf( gp, "set terminal gif\n" );								//decide the exporting file type of graphs (here, the file type is .gif)
														fprintf( gp, "set grid\n" );										//draw grid in the graph
														fprintf( gp, "set title't = %f[s], z = %f[cm]'\n", tt[k], zz[k]*100.0 );//graph title
														fprintf( gp, "set xlabel 'x (m)'\n" );								//decide the x label
														fprintf( gp, "set ylabel 'y (m)'\n" );								//decide the y label
														fprintf( gp, "set zlabel 'z (m)'\n" );								//decide the z label
														fprintf( gp, "set xrange[%f:%f]\n", xlim_l, xlim_u );				//confine the length of x axis
														fprintf( gp, "set yrange[%f:%f]\n", ylim_l, ylim_u );				//confine the length of y axis
														fprintf( gp, "set zrange[%f:%f]\n", zlim_l, zlim_u );				//confine the length of z axis
														//fprintf( gp, "set view equal xyz\n" );							//equalize the ratio of x and y axes
														fprintf( gp, "set ticslevel 0\n");									//set the bottom of z axis on x-y plane
														fprintf( gp, "set output '%s_pos3D_%d_%d.gif'\n", var, i, k );		//decide the file name
														fprintf( gp, "splot \"-\" with points pt %d ps %f\n", type, plots );//plot the atom position using a circle
														fprintf( gp, "%f, %f, %f\n", xx[k], yy[k], zz[k] );
														fprintf( gp, "e\n" );												//e means the end of inputs of coordinates
														fflush( gp );
														}
														fprintf( gp, "exit\n");
														_pclose( gp );
													}
												}
											//===plot (the way of plotting is replacement) <end>====================================================================

											//===plot (the way of plotting is superposition)========================================================================
												if( rep_sup==2 ){
													if( (D_2_3 == 2 || D_2_3 == 4) && (xyz==1 || xyz==3) ){
														for( int k=1; k<=j; k+=intv ){
															fprintf( gp, "set terminal gif\n" );								//decide the exporting file type of graphs (here, the file type is .gif)
															fprintf( gp, "set grid\n" );										//draw grid in the graph
															fprintf( gp, "set title't = %f[s], z = %f[cm]'\n", tt[k], zz[k]*100.0 );//graph title
															fprintf( gp, "set xlabel 'x (m)'\n" );								//decide the x label
															fprintf( gp, "set ylabel 'y (m)'\n" );								//decide the y label
															fprintf( gp, "set xtics rotate by 45 center\n" );					//rotate scales 45 degrees
															fprintf( gp, "set xtics center offset 0,-1\n" );					//change the offset of scales on x axis
															fprintf( gp, "set xrange[%f:%f]\n", xlim_l, xlim_u );				//confine the length of x axis
															fprintf( gp, "set yrange[%f:%f]\n", ylim_l, ylim_u );				//confine the length of y axis
															fprintf( gp, "set size ratio -1\n" );								//equalize the ratio of x and y axes
															fprintf( gp, "set output '%s_pos2D_xy_%d_%d.gif'\n", var, i, k );	//decide the file name
															fprintf( gp, "plot \"-\" with points pt %d ps %f\n", type, plots );	//plot the atom position using a circle
															for(int l=1; l<=k; l++){
															fprintf( gp, "%f, %f\n", xx[l], yy[l] );							//plot the atom position on x-y plane
															}
															fprintf( gp, "e\n" );												//e means the end of inputs of coordinates
															fflush( gp );

															if( k!=j && k+intv > j ){
															fprintf( gp, "set terminal gif\n" );								//decide the exporting file type of graphs (here, the file type is .gif)
															fprintf( gp, "set grid\n" );										//draw grid in the graph
															fprintf( gp, "set title't = %f[s], z = %f[cm]'\n", tt[k], zz[k]*100.0 );//graph title
															fprintf( gp, "set xlabel 'x (m)'\n" );								//decide the x label
															fprintf( gp, "set ylabel 'y (m)'\n" );								//decide the y label
															fprintf( gp, "set xtics rotate by 45 center\n" );					//rotate scales 45 degrees
															fprintf( gp, "set xtics center offset 0,-1\n" );					//change the offset of scales on x axis
															fprintf( gp, "set xrange[%f:%f]\n", xlim_l, xlim_u );				//confine the length of x axis
															fprintf( gp, "set yrange[%f:%f]\n", ylim_l, ylim_u );				//confine the length of y axis
															fprintf( gp, "set size ratio -1\n" );								//equalize the ratio of x and y axes
															fprintf( gp, "set output '%s_pos2D_xy_%d_%d.gif'\n", var, i, j );	//decide the file name
															fprintf( gp, "plot \"-\" with points pt %d ps %f\n", type, plots );	//plot the atom position using a circle
															for(int l=1; l<=j; l++){
															fprintf( gp, "%f, %f\n", xx[l], yy[l] );							//plot the atom position on x-y plane
															}
															fprintf( gp, "e\n" );												//e means the end of inputs of coordinates
															fflush( gp );
															}
														}
														if( D_2_3 == 2 && xyz==1 ){
															fprintf( gp, "exit\n");
															_pclose( gp );
														}
													}
													if( (D_2_3 == 2 || D_2_3 == 4) && (xyz==2 || xyz==3) ){
														for( int k=1; k<=j; k+=intv ){
															fprintf( gp, "set terminal gif\n" );								//decide the exporting file type of graphs (here, the file type is .gif)
															fprintf( gp, "set grid\n" );										//draw grid in the graph
															fprintf( gp, "set title't = %f[s], z = %f[cm]'\n", tt[k], zz[k]*100.0 );//graph title
															fprintf( gp, "set xlabel 'x (m)'\n" );								//decide the x label
															fprintf( gp, "set ylabel 'z (m)'\n" );								//decide the z label
															fprintf( gp, "set xtics rotate by 45 center\n" );					//rotate scales 45 degrees
															fprintf( gp, "set xtics center offset 0,-1\n" );					//change the offset of scales on x axis
															fprintf( gp, "set xrange[%f:%f]\n", xlim_l, xlim_u );				//confine the length of x axis
															fprintf( gp, "set yrange[%f:%f]\n", zlim_l, zlim_u );				//confine the length of z axis
															fprintf( gp, "set size noratio\n" );								//cancel of equalization of the ratio of axes
															fprintf( gp, "set output '%s_pos2D_xz_%d_%d.gif'\n", var, i, k );	//decide the file name
															fprintf( gp, "plot \"-\" with points pt %d ps %f\n", type, plots );	//plot the atom position using a circle
															for(int l=1; l<=k; l++){
															fprintf( gp, "%f, %f\n", xx[l], zz[l] );							//plot the atom position on x-z plane
															}
															fprintf( gp, "e\n" );												//e means the end of inputs of coordinates
															fflush( gp );

															if( k!=j && k+intv > j ){
															fprintf( gp, "set terminal gif\n" );								//decide the exporting file type of graphs (here, the file type is .gif)
															fprintf( gp, "set grid\n" );										//draw grid in the graph
															fprintf( gp, "set title't = %f[s], z = %f[cm]'\n", tt[k], zz[k]*100.0 );//graph title
															fprintf( gp, "set xlabel 'x (m)'\n" );								//decide the x label
															fprintf( gp, "set ylabel 'z (m)'\n" );								//decide the y label
															fprintf( gp, "set xtics rotate by 45 center\n" );					//rotate scales 45 degrees
															fprintf( gp, "set xtics center offset 0,-1\n" );					//change the offset of scales on x axis
															fprintf( gp, "set xrange[%f:%f]\n", xlim_l, xlim_u );				//confine the length of x axis
															fprintf( gp, "set yrange[%f:%f]\n", zlim_l, zlim_u );				//confine the length of y axis
															fprintf( gp, "set size noratio\n" );								//cancel of equalization of the ratio of axes
															fprintf( gp, "set output '%s_pos2D_xz_%d_%d.gif'\n", var, i, j );	//decide the file name
															fprintf( gp, "plot \"-\" with points pt %d ps %f\n", type, plots );	//plot the atom position using a circle
															for(int l=1; l<=j; l++){
															fprintf( gp, "%f, %f\n", xx[l], zz[l] );							//plot the atom position on x-z plane
															}
															fprintf( gp, "e\n" );												//e means the end of inputs of coordinates
															fflush( gp );
															}
														}
														if( D_2_3 == 2 ){
															fprintf( gp, "exit\n");
															_pclose( gp );
														}
													}
													if( D_2_3 == 3 || D_2_3 == 4 ){
														for( int k=1; k<=j; k+=intv ){
															fprintf( gp, "set terminal gif\n" );								//decide the exporting file type of graphs (here, the file type is .gif)
															fprintf( gp, "set grid\n" );										//draw grid in the graph
															fprintf( gp, "set title't = %f[s], z = %f[cm]'\n", tt[k], zz[k]*100.0 );//graph title
															fprintf( gp, "set xlabel 'x (m)'\n" );								//decide the x label
															fprintf( gp, "set ylabel 'y (m)'\n" );								//decide the y label
															fprintf( gp, "set zlabel 'z (m)'\n" );								//decide the z label
															fprintf( gp, "set xrange[%f:%f]\n", xlim_l, xlim_u );				//confine the length of x axis
															fprintf( gp, "set yrange[%f:%f]\n", ylim_l, ylim_u );				//confine the length of y axis
															fprintf( gp, "set zrange[%f:%f]\n", zlim_l, zlim_u );				//confine the length of z axis
															//fprintf( gp, "set view equal xyz\n" );							//equalize the ratio of x and y axes
															fprintf( gp, "set ticslevel 0\n");									//set the bottom of z axis on x-y plane
															fprintf( gp, "set output '%s_pos3D_%d_%d.gif'\n", var, i, k );		//decide the file name
															fprintf( gp, "splot \"-\" with points pt %d ps %f\n", type, plots );//plot the atom position using a circle
															for(int l=1; l<=k; l++){
															fprintf( gp, "%f, %f, %f\n", xx[l], yy[l], zz[l] );
															}
															fprintf( gp, "e\n" );												//e means the end of inputs of coordinates
															fflush( gp );

															if( k!=j && k+intv > j ){
															fprintf( gp, "set terminal gif\n" );								//decide the exporting file type of graphs (here, the file type is .gif)
															fprintf( gp, "set grid\n" );										//draw grid in the graph
															fprintf( gp, "set title't = %f[s], z = %f[cm]'\n", tt[k], zz[k]*100.0 );//graph title
															fprintf( gp, "set xlabel 'x (m)'\n" );								//decide the x label
															fprintf( gp, "set ylabel 'y (m)'\n" );								//decide the y label
															fprintf( gp, "set zlabel 'z (m)'\n" );								//decide the z label
															fprintf( gp, "set xrange[%f:%f]\n", xlim_l, xlim_u );				//confine the length of x axis
															fprintf( gp, "set yrange[%f:%f]\n", ylim_l, ylim_u );				//confine the length of y axis
															fprintf( gp, "set zrange[%f:%f]\n", zlim_l, zlim_u );				//confine the length of z axis
															//fprintf( gp, "set view equal xyz\n" );							//equalize the ratio of x and y axes
															fprintf( gp, "set ticslevel 0\n");									//set the bottom of z axis on x-y plane
															fprintf( gp, "set output '%s_pos3D_%d_%d.gif'\n", var, i, j );		//decide the file name
															fprintf( gp, "plot \"-\" with points pt %d ps %f\n", type, plots );	//plot the atom position using a circle
															for(int l=1; l<=j; l++){
															fprintf( gp, "%f, %f, %f\n", xx[l], yy[l], zz[l] );
															}
															fprintf( gp, "e\n" );												//e means the end of inputs of coordinates
															fflush( gp );
															}
														}
														fprintf( gp, "exit\n");
														_pclose( gp );
													}
												}
											//===plot (the way of plotting is superposition) <end>==================================================================

										}

										if( gnu==1 && lump==1 ){
											if( ( gp = _popen(GNUPLOT_PATH, "w") ) == NULL ){							//open GNUPLOT by using "_popen" command
												fprintf( stderr, "The application could not be booted.%s\n", GNUPLOT_PATH );
												exit(EXIT_FAILURE);
											}

											//===plot==============================================================================================================
												if( (D_2_3 == 2 || D_2_3 == 4) && (xyz==1 || xyz==3) ){
													fprintf( gp, "set terminal gif\n" );								//decide the exporting file type of graphs (here, the file type is .gif)
													fprintf( gp, "set grid\n" );										//draw grid in the graph
													fprintf( gp, "set title't = %f[s], z = %f[cm]'\n", t_t, z_t*100.0 );//graph title
													fprintf( gp, "set xlabel 'x (m)'\n" );								//decide the x label
													fprintf( gp, "set ylabel 'y (m)'\n" );								//decide the y label
													fprintf( gp, "set xtics rotate by 45 center\n" );					//rotate scales 45 degrees
													fprintf( gp, "set xtics center offset 0,-1\n" );					//change the offset of scales on x axis
													fprintf( gp, "set xrange[%f:%f]\n", xlim_l, xlim_u );				//confine the length of x axis
													fprintf( gp, "set yrange[%f:%f]\n", ylim_l, ylim_u );				//confine the length of y axis
													fprintf( gp, "set size ratio -1\n" );								//equalize the ratio of x and y axes
													fprintf( gp, "set output '%s_pos2D_xy_%d.gif'\n", var, i );			//decide the file name
													fprintf( gp, "plot \"-\" with points pt %d ps %f\n", type, plots );	//plot the atom position using a circle
													for( int k=1; k<=j; k++ ){
													fprintf( gp, "%f, %f\n", xx[k], yy[k] );							//plot the atom position on x-z plane
													}
													fprintf( gp, "e\n" );												//e means the end of inputs of coordinates
													fflush( gp );
													if( D_2_3 == 2 && xyz==1 ){
														fprintf( gp, "exit\n");
														_pclose( gp );
													}
												}
												if( (D_2_3 == 2 || D_2_3 == 4) && (xyz==2 || xyz==3) ){
													fprintf( gp, "set terminal gif\n" );								//decide the exporting file type of graphs (here, the file type is .gif)
													fprintf( gp, "set grid\n" );										//draw grid in the graph
													fprintf( gp, "set title't = %f[s], z = %f[cm]'\n", t_t, z_t*100.0 );//graph title
													fprintf( gp, "set xlabel 'x (m)'\n" );								//decide the x label
													fprintf( gp, "set ylabel 'z (m)'\n" );								//decide the z label
													fprintf( gp, "set xtics rotate by 45 center\n" );					//rotate scales 45 degrees
													fprintf( gp, "set xtics center offset 0,-1\n" );					//change the offset of scales on x axis
													fprintf( gp, "set xrange[%f:%f]\n", xlim_l, xlim_u );				//confine the length of x axis
													fprintf( gp, "set yrange[%f:%f]\n", zlim_l, zlim_u );				//confine the length of z axis
													fprintf( gp, "set size noratio\n" );								//cancel of equalization of the ratio of axes
													fprintf( gp, "set output '%s_pos2D_xz_%d.gif'\n", var, i );			//decide the file name
													fprintf( gp, "plot \"-\" with points pt %d ps %f\n", type, plots );	//plot the atom position using a circle
													for( int k=1; k<=j; k++ ){
													fprintf( gp, "%f, %f\n", xx[k], zz[k] );							//plot the atom position on x-z plane
													}
													fprintf( gp, "e\n" );												//e means the end of inputs of coordinates
													fflush( gp );
													if( D_2_3 == 2 ){
														fprintf( gp, "exit\n");
														_pclose( gp );
													}
												}
												if( D_2_3 == 3 || D_2_3 == 4 ){
													fprintf( gp, "set terminal gif\n" );								//decide the exporting file type of graphs (here, the file type is .gif)
													fprintf( gp, "set grid\n" );										//draw grid in the graph
													fprintf( gp, "set title't = %f[s], z = %f[cm]'\n", t_t, z_t*100.0 );//graph title
													fprintf( gp, "set xlabel 'x (m)'\n" );								//decide the x label
													fprintf( gp, "set ylabel 'y (m)'\n" );								//decide the y label
													fprintf( gp, "set zlabel 'z (m)'\n" );								//decide the z label
													fprintf( gp, "set xrange[%f:%f]\n", xlim_l, xlim_u );				//confine the length of x axis
													fprintf( gp, "set yrange[%f:%f]\n", ylim_l, ylim_u );				//confine the length of y axis
													fprintf( gp, "set zrange[%f:%f]\n", zlim_l, zlim_u );				//confine the length of z axis
													//fprintf( gp, "set view equal xyz\n" );							//equalize the ratio of x and y axes
													fprintf( gp, "set ticslevel 0\n");									//set the bottom of z axis on x-y plane
													fprintf( gp, "set output '%s_pos3D_%d.gif'\n", var, i );			//decide the file name
													fprintf( gp, "splot \"-\" with points pt %d ps %f\n", type, plots );//plot the atom position using a circle
													for( int k=1; k<=j; k++ ){
													fprintf( gp, "%f, %f, %f\n", xx[k], yy[k], zz[k] );
													}
													fprintf( gp, "e\n" );												//e means the end of inputs of coordinates
													fflush( gp );
													fprintf( gp, "exit\n");
													_pclose( gp );
												}
											//===plot <end>========================================================================================================
										}
									//===plot position <end>=======================================================================================================

									//====output data into csv files====
									fprintf( passed_time, "%.16e\n", t_t );
									fprintf( f_position, "%.16e,%.16e,%.16e,%.16e\n", rad_t, x_t, y_t, z_t );
									if( sqrt(x*x+y*y) > r_simu ){
										fprintf( f_state, "%d\n", 0 );			//"0" means "out of beam"
									}else if( z < -0.26 ){
										if( rad <= r_judge ){
											fprintf( f_state, "%d\n", 1 );		//"1" means "substrate" (If the radius is same as or smaller than the inner diameter of the beam, this program judges it "substrate".)
										}else if( rad > r_judge ){
											fprintf( f_state, "%d\n", 0 );		//"0" means "out of beam" (If the radius is larger than the inner diameter of the beam, this program judges it "out of beam".)
										}
									}else if( j==jloop ){
										fprintf( f_state, "%d\n", 2 );			//"2" means "still dropping"
									}
									fprintf( f_velo, "%.16e,%.16e,%.16e\n", vx_t, vy_t, vz_t );
									fprintf( f_spon, "%d\n", spon );
									//fprintf( r_max, "%d\n", rmax );

									//release memories
									if(gnu==1){
										delete[] xx; delete[] yy; delete[] zz; delete[] tt;
									}
									delete[] r1; delete[] r2; delete[] r3; delete[] r4; delete[] r5;

									//====output data by compulsion into csv files if there are unwrittend data on the RAM==================
									fflush( x_y_z_coord );		fflush( velo );			fflush( x_y_z_velo );		fflush( r_phi_coord );
									fflush( r_phi_velo );		fflush( kinetic );		fflush( momentum );			fflush( Mom_x_y_z );
									fflush( angmomentum );		fflush( dipole );		fflush( Udip );				fflush( Energy );
									fflush( dEnergy );			fflush( dMomentum );	fflush( dMom_x_y_z );		fflush( psp );
									fflush( Ugrav );

									//====close files====
									fclose( x_y_z_coord );		fclose( velo );			fclose( x_y_z_velo );		fclose( r_phi_coord );
									fclose( r_phi_velo );		fclose( kinetic );		fclose( momentum );			fclose( Mom_x_y_z );
									fclose( angmomentum );		fclose( dipole );		fclose( Udip );				fclose( Energy );
									fclose( dEnergy );			fclose( dMomentum );	fclose( dMom_x_y_z );		fclose( psp );
									fclose( Ugrav );

									if( sqrt(x*x+y*y) > r_simu || z < -0.26 ){
										break;
									}
								}else{

									if( csv_out==1 ){
										if( state_t==1&&state==1 ){
											//when a spontaneous emission has not occured
											fprintf( dEnergy, "%.16e\n", 1.0/2.0*mass*(v*v-v_t*v_t) + (U1_opt(int_TEM(x,y,I_0,waist), delta)-U1_opt(int_TEM(x_t,y_t,I_0,waist), delta)) + mass*G*(z-z_t) );
											fprintf( dMom_x_y_z, "%.16e,%.16e,%.16e\n", mass*vx - (mass*vx_t + dipole_f_x1( x_t, y_t, int_TEM(x_t,y_t,I_0,waist), d_int_TEM(x_t,y_t,I_0,waist) )*dt), mass*vy - (mass*vy_t + dipole_f_y1( x_t, y_t, int_TEM(x_t,y_t,I_0,waist), d_int_TEM(x_t,y_t,I_0,waist) )*dt), mass*vz - (mass*vz_t - mass*G*dt) );
										}else if( state_t==2&&state==2 ){
											//when a spontaneous emission has not occured
											fprintf( dEnergy, "%.16e\n", 1.0/2.0*mass*(v*v-v_t*v_t) + (U2_opt(int_TEM(x,y,I_0,waist), delta)-U2_opt(int_TEM(x_t,y_t,I_0,waist), delta)) + mass*G*(z-z_t) );
											fprintf( dMom_x_y_z, "%.16e,%.16e,%.16e\n", mass*vx - (mass*vx_t + dipole_f_x2( x_t, y_t, int_TEM(x_t,y_t,I_0,waist), d_int_TEM(x_t,y_t,I_0,waist) )*dt), mass*vy - (mass*vy_t + dipole_f_y2( x_t, y_t, int_TEM(x_t,y_t,I_0,waist), d_int_TEM(x_t,y_t,I_0,waist) )*dt), mass*vz - (mass*vz_t - mass*G*dt) );
										}else if( state_t==1&&state==2 ){
											//when a spontaneous emission has occured
											fprintf( dEnergy, "%.16e\n", hbar*2.0*M_PI*c/lambda + 1.0/2.0*mass*(v_tt*v_tt-v_t*v_t) + (U2_opt(int_TEM(x_t,y_t,I_0,waist), delta)-U1_opt(int_TEM(x_t,y_t,I_0,waist), delta)) + (U1_opt(int_TEM(x_t,y_t,I_0,waist), delta)-U2_opt(int_TEM(x_t,y_t,I_0,waist), delta)) + 1.0/2.0*mass*(v*v-v_tt*v_tt) + (U2_opt(int_TEM(x,y,I_0,waist), delta)-U2_opt(int_TEM(x_t,y_t,I_0,waist), delta)) + mass*G*(z-z_t) );
											fprintf( dMom_x_y_z, "%.16e,%.16e,%.16e\n", mass*vx - (mass*vx_t + hbar*k_wave*sin(theta2)*cos(phi2) + dipole_f_x2( x_t, y_t, int_TEM(x_t,y_t,I_0,waist), d_int_TEM(x_t,y_t,I_0,waist) )*dt), mass*vy - (mass*vy_t + hbar*k_wave*sin(theta2)*sin(phi2) +dipole_f_y2( x_t, y_t, int_TEM(x_t,y_t,I_0,waist), d_int_TEM(x_t,y_t,I_0,waist) )*dt), mass*vz - (mass*vz_t + hbar*k_wave*cos(theta2) - mass*G*dt) );
										}else if( state_t==2&&state==1 ){
											//when a spontaneous emission has occured
											fprintf( dEnergy, "%.16e\n", hbar*2.0*M_PI*c/lambda + 1.0/2.0*mass*(v_tt*v_tt-v_t*v_t) + (U1_opt(int_TEM(x_t,y_t,I_0,waist), delta)-U2_opt(int_TEM(x_t,y_t,I_0,waist), delta)) + (U2_opt(int_TEM(x_t,y_t,I_0,waist), delta)-U1_opt(int_TEM(x_t,y_t,I_0,waist), delta)) + 1.0/2.0*mass*(v*v-v_tt*v_tt) + (U1_opt(int_TEM(x,y,I_0,waist), delta)-U1_opt(int_TEM(x_t,y_t,I_0,waist), delta)) + mass*G*(z-z_t) );
											fprintf( dMom_x_y_z, "%.16e,%.16e,%.16e\n", mass*vx - (mass*vx_t + hbar*k_wave*sin(theta2)*cos(phi2) + dipole_f_x1( x_t, y_t, int_TEM(x_t,y_t,I_0,waist), d_int_TEM(x_t,y_t,I_0,waist) )*dt), mass*vy - (mass*vy_t + hbar*k_wave*sin(theta2)*sin(phi2) +dipole_f_y1( x_t, y_t, int_TEM(x_t,y_t,I_0,waist), d_int_TEM(x_t,y_t,I_0,waist) )*dt), mass*vz - (mass*vz_t + hbar*k_wave*cos(theta2) - mass*G*dt) );
										}
										fprintf( dMomentum, "%.16e\n", mass*v - mass*v_t );
									}

									//===advance dt on temporary variables===========================================================================
										//position
										phi_t = phi;
										x_t = x; y_t = y; z_t = z;
										rad_t = rad;

										//velocity
										vx_t = vx; vy_t = vy; vz_t = vz;
										v_t = v;
										vr_t = vr; vphi_t = vphi;

										//time
										t_t = t;

										//eigenstate
										state_t = state;
									//===advance dt on temporary variables <end>=====================================================================
								}
							//===check whether the atom is "out of the beam", "on the substrate" or "still dropping" and plot the trajectory <end>=================

						}
					//===for loop (repeat the number of exporting times (jloop)) <end>=============================================================================




					//===close csv files if SAMPLE is at the end===================================================================================================
						if( i==SAMPLE ){		//If SAMPLE loop (i loop) is at the end, csv files are closed.
							//====output data by compulsion into csv files if there are unwrittend data on the RAM====
							fflush( i_velo_on_xy );		fflush( i_velo_x_y_z );		fflush( i_velo );		fflush( i_velo_r_phi );
							fflush( i_x_y_z );			fflush( passed_time );		fflush( f_position );	fflush( r_max );
							fflush( f_state );			fflush( i_angmomentum );	fflush( i_kinetic );	fflush( f_velo );
							fflush( f_spon );			fflush( eigenstate );

							//====close files====
							fclose( i_velo_on_xy );		fclose( i_velo_x_y_z );		fclose( i_velo );		fclose( i_velo_r_phi );
							fclose( i_x_y_z );			fclose( passed_time );		fclose( f_position );	fclose( r_max );
							fclose( f_state );			fclose( i_angmomentum );	fclose( i_kinetic );	fclose( f_velo );
							fclose( f_spon );			fclose( eigenstate );

							//====release memories====
							delete[] r_xx; delete[] r_yy; delete[] r_zz; delete[] v_xx; delete[] v_yy; delete[] v_zz;
						}
					//===close csv files if SAMPLE is at the end <end>=============================================================================================

				}
			//===for loop (repeat the number of SAMPLES) <end>=====================================================================================================

		}
	//===for loop (repeat changing one variable) <end>=============================================================================================================

	return 0;
}

//===MAIN PROGRAM <end>============================================================================================================================================













//===Function Definiton======================================================================================================
	//generate a random number within (0,1)
	double Uniform_A(void){
		double i;
		while(1){
			i = (double)rand() / (double)RAND_MAX;//i is in (0,1)
			if( i==0.0 || i==1.0 ) continue;
			else break;
		} 
		return i;
	}

	//generate a random number within [0,1]
	double Uniform_B1(void){
		return (double)rand() / (double)RAND_MAX;//returned number is in [0,1]
	}

	//generate a random number within [0,1]
	double Uniform_B2(void){
		return (double)rand() / (double)RAND_MAX;//returned number is in [0,1]
	}

	//generate a random number within [0,1]
	double Uniform_B3(void){
		return (double)rand() / (double)RAND_MAX;//returned number is in [0,1]
	}

	//generate a random number within [0,1]
	double Uniform_B4(void){
		return (double)rand() / (double)RAND_MAX;//returned number is in [0,1]
	}

	//generate a random number within [0,1]
	double Uniform_B5(void){
		return (double)rand() / (double)RAND_MAX;//returned number is in [0,1]
	}

	//Gaussian distribution of position (x or y or z component)
	//This function is defined only to simplify the definition of the function gauss_pos_x(y,z).
	double gauss(double x){
		return sqrt( 1.0/(2.0*M_PI*1.0e-6/9.0) )*exp( -x*x/(2.0*1.0e-6/9.0) );
	}

	//random position (x component) generator based on Boltzmann distribution
	double gauss_pos_x(double x){ //sigma3 should be substituted into x
		double xi,yj;
		do{
				xi = 2.0*x*rand()/RAND_MAX-x;//generate uniform random no. xi between [-x,x] (x=n*wr/2)
				yj = 1205.0*rand()/RAND_MAX;//1205.0 is a value that is a little bit greater than the max of gauss()
			}while( gauss(xi)<yj );
		return xi;
	}

	//random position (y component) generator based on Boltzmann distribution
	double gauss_pos_y(double x){ //sigma3 should be substituted into x
		double xi,yj;
		do{
				xi = 2.0*x*rand()/RAND_MAX-x;//generate uniform random no. xi between [-x,x] (x=n*wr/2)
				yj = 1205.0*rand()/RAND_MAX;//1205.0 is a value that is a little bit greater than the max of gauss()
			}while( gauss(xi)<yj );
		return xi;
	}

	//random position (z component) generator based on Boltzmann distribution
	double gauss_pos_z(double x){ //sigma3 should be substituted into x
		double xi,yj;
		do{
				xi = 2.0*x*rand()/RAND_MAX-x;//generate uniform random no. xi between [-x,x] (x=n*wr/2)
				yj = 1205.0*rand()/RAND_MAX;//1205.0 is a value that is a little bit greater than the max of gauss()
			}while( gauss(xi)<yj );
		return xi;
	}

	//Maxwell-Boltzmann distribution of velocity (x or y or z component)
	//This function is defined only to simplify the definition of the function max_boltz_velo_x(y,z).
	double max_boltz(double x, double t){	//temp should be substituted into t
		return sqrt( mass/(2.0*M_PI*k_b*t) )*exp( -mass*x*x/(2.0*k_b*t) );
	}

	//random velocity (x component) generator based on Maxwell-Boltzmann distribution
	double max_boltz_velo_x(double x){ //if x=0.15, then random no. between [-0.15[m/s],0.15[m/s]] are generated
		double xi,yj;
		do{
				xi = 2.0*x*rand()/RAND_MAX-x;//generate uniform random no. xi between [-x,x] (x=n*wr/2)
				yj = 12.8*rand()/RAND_MAX;//12.8 is a value that is a little bit greater than the max of max_boltz()
			}while( max_boltz(xi, temp)<yj );
		return xi;
	}

	//random velocity (y component) generator based on Maxwell-Boltzmann distribution
	double max_boltz_velo_y(double x){ //if x=0.15, then random no. between [-0.15[m/s],0.15[m/s]] are generated
		double xi,yj;
		do{
				xi = 2.0*x*rand()/RAND_MAX-x;//generate uniform random no. xi between [-x,x] (x=n*wr/2)
				yj = 12.8*rand()/RAND_MAX;//12.8 is a value that is a little bit greater than the max of max_boltz()
			}while( max_boltz(xi, temp)<yj );
		return xi;
	}

	//random velocity (z component) generator based on Maxwell-Boltzmann distribution
	double max_boltz_velo_z(double x){ //if x=0.15, then random no. between [-0.15[m/s],0.15[m/s]] are generated
		double xi,yj;
		do{
				xi = 2.0*x*rand()/RAND_MAX-x;//generate uniform random no. xi between [-x,x] (x=n*wr/2)
				yj = 12.8*rand()/RAND_MAX;//12.8 is a value that is a little bit greater than the max of max_boltz()
			}while( max_boltz(xi, temp)<yj );
		return xi;
	}

	//intensity of TEM01*
	double int_TEM( double x, double y, double i, double w ){	//I_0 and waist should be substituted into i and w, respectively.
		return i*2.0*(x*x+y*y)/(w*w)*exp( -2.0*(x*x+y*y)/(w*w) );
	}

	//differentiated intensity of TEM01*
	double d_int_TEM( double x, double y, double i, double w ){		//I_0 and waist should be substituted into i and w, respectively.
		return i*4.0*sqrt(x*x+y*y)/(w*w)*( 1.0 - 2.0*(x*x+y*y)/(w*w) )*exp( -2.0*(x*x+y*y)/(w*w) );
	}

	//optical potential on |1>
	double U1_opt( double z, double d ){	//Intensity and detuning must be substituted into z and d, respectively.
		return hbar*d/2.0*log( 1.0 + 1.0/I_s1*gamma1*gamma1/(gamma1*gamma1 + 4.0*d*d)*z );
	}

	//optical potential on |2>
	double U2_opt( double z, double d ){	//Intensity and detuning must be substituted into z and d, respectively.
		return hbar*(d+delta_hfs)/2.0*log( 1.0 + 1.0/I_s2*gamma2*gamma2/(gamma2*gamma2 + 4.0*(d+delta_hfs)*(d+delta_hfs))*z );
	}

	//dipole force on |1>
	double dipole_f1( double z, double w, double d ){		//Intensity, differentiated intensity and detuning must be substituted into z, w and d, respectively.
		return -hbar*d*1.0/(2.0*I_s1)*gamma1*gamma1*w/( gamma1*gamma1 + 4.0*d*d + 1.0/I_s1*gamma1*gamma1*z );
	}

	//dipole force on |2>
	double dipole_f2( double z, double w, double d ){		//Intensity, differentiated intensity and detuning must be substituted into z, w and d, respectively.
		return -hbar*( d + delta_hfs )*1.0/(2.0*I_s2)*gamma2*gamma2*w/( gamma2*gamma2 + 4.0*(d+delta_hfs)*(d+delta_hfs) + 1.0/I_s2*gamma2*gamma2*z );
	}

	//x component of dipole force on |1>
	double dipole_f_x1( double x, double y, double z, double w ){
		return dipole_f1( z,w,delta )*x/sqrt(x*x+y*y);
	}

	//y component of dipole force on |1>
	double dipole_f_y1( double x, double y, double z, double w ){
		return dipole_f1( z,w,delta )*y/sqrt(x*x+y*y);
	}

	//x component of dipole force on |2>
	double dipole_f_x2( double x, double y, double z, double w ){
		return dipole_f2( z,w,delta )*x/sqrt(x*x+y*y);
	}

	//y component of dipole force on |2>
	double dipole_f_y2( double x, double y, double z, double w ){
		return dipole_f2( z,w,delta )*y/sqrt(x*x+y*y);
	}

//===Function Definition <end>==============================================================================================

//END OF FILE