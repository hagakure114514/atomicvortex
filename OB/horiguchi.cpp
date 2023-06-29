#include <stdafx.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <iomanip>
#include <time.h>
#include <strstream>
#include <stdio.h>
#include <complex>
#include <sstream>
#include <fstream>

#define GNUPLOT_PATH "C:/Users/itolab/Desktop/gnuplot/gnuplot/binary/pgnuplot.exe" 
#define M_PI 3.14159265359
#define G 9.8
#define MASS (87E-3/6.022E+23)
#define HBAR 1.0545E-34
#define BOLTZ 1.380658E-23
#define J_to_K 7.2464e22
#define box_number 100
#define RAMUDA 780.24e-9
#define c 2.99792458e8
#define Is 16.4//飽和強度[W/m^2]
#define sizenn_haba (6.1e6 * 2 * M_PI)//自然幅[Hz]
#define h_bar 1.054e-34//ディラック定数[J/s]
#define delta_1 6835.0e6//基底準位間幅[Hz]

#define w 6.2e-3



//いじる数字
#define MOT 9.1e7//原子数
#define n1 1.51
#define r0 1e-3//MOT原子郡半径 単位[m]
#define r1 100e-06//開口半径 単位[m]
#define r2 3.0e-3//ファネル厚さ（水平方向） 単位[m]
#define slope (tan((M_PI*(42.7))/180))//ファネル壁面傾き
#define HEIGHT (r2+2*r1+r0)//MOT中心高さz座標 単位[m]
#define TEMP 9.0E-6//MOT初期温度[K]
#define SAMPLE 1000
#define time_max 1.0//シミュレーション時間（原子の許容最大滞在時間）[s]
#define time_with 1e-2//time_max/box_number//フラックス強度評価時間[s]
#define gennsui_tyou 288e-9
#define capplar_h 0.0e-3//カップラー高さ[m] カップラー無の時
//#define capplar_h r2//0.3e-3//カップラー高さ[m]  カップラー有の時
#define w2 170e-06//ビーム中空半径[m]c+0.5e-3
//#define w2 (r1+capplar_h)//ビーム中空半径[m] カップラー有の時
#define w1 (r1+r2-1.2e-3)//ビーム半径[m]r-1.2(200u)1.4(100u),c-0
//#define w1 (r1+r2+capplar_h)//ビーム半径[m] カップラー有の時
#define beam_error_x 0.0e-6//中空ビームのx方向位置ずれ[m]
	#define P 200e-03//ビームパワー[W]　
#define Enhance 7.9*3.5//エンハンス 円偏光の時*0.5になるようにしてある
#define cappling_rate_fanel 1//ファネルのエネルギー伝達効率
#define cappling_rate_coupler 0.5//カップラーのエネルギー伝達効率
#define cappling_rate cappling_rate_fanel*cappling_rate_coupler//システム全体のエネルギー伝達効率
#define delta (1e9*2*M_PI)//レーザー周波数離調[Hz]
#define dt 1.0e-5
#define Ip 6.0e4//許容残留強度[W/m^2] 散逸力が無視できない強度

#define U1max 7.34e-28
#define U2max 0.938e-28


double para = r1;
std::string para_name = "r1";

//手計算の後、毎回入力
#define a 100//time_max / time_with
#define b 10000//time_max / dt



const double rand_max = RAND_MAX;


//関数定義
void fig_3D_xyz0(double x[], double y[], double z[], int n);
void fig_3D_xyz_at_fanel(double x[], double y[], double z[], int n);
void fig_2D_r_U(double x[], double y[], int m);
double Uniform(void);
double rand_gauss(double mu, double sigma);
void fig_BOLTZ_2D(double x[], int o);
int reflection_jadge(double xa, double ya, double za, double vxa, double vya, double vza);
double U_laser(double xa, double ya);
double reflection_vx(double xa, double ya, double za, double vxa, double vya, double vza);
double reflection_vy(double xa, double ya, double za, double vxa, double vya, double vza);
double reflection_vz(double xa, double ya, double za, double vxa, double vya, double vza);
double vr1(double xa, double ya, double za, double vxa, double vya, double vza);




double vr2(double xa, double ya, double za, double vxa, double vya, double vza);
void fig_flux_intensity_2D(double x[], int o);
void fig_laser_3D_to_2Dmap(void);
double I_laser(double xa, double ya);
double I_laser_fanel(double xa, double ya);

double Lpotential_U1(double x);
double Lpotential_U2(double x);
double Lpotential_U3(double x);

double U1_A(double x);
double U1_B(double x);
double U1_C(double x);

double U2_A(double x);
double U2_B(double x);

int main(){



	srand((unsigned)time(NULL));
	double vx0_data[SAMPLE], vy0_data[SAMPLE], vz0_data[SAMPLE], v0[SAMPLE], x0_data[SAMPLE], y0_data[SAMPLE], z0_data[SAMPLE];
	double gauss_sigma = sqrt(BOLTZ*TEMP / MASS);
	double gauss_mu = 0;
	double prev_x, prev_y, prev_z, prev_vx, prev_vy, prev_vz, prev_t, t = 0;
	double x, y, z, vx, vy, vz;
	int sample_No = 0;
	int loss_jadge;
	double dt_1 = dt;
	int time_back = 0;
	double loss_count = 0;
	double loss_count_Ip = 0;
	double gain_count = 0;
	double time_out_count = 0;
	int flux[a * 10];
	double flux_intensity[a * 10];
	int lossjadge = 0;
	int gainjadge = 0;
	int losscount = 0;
	int gaincount = 0;
	double Kinetic1 = 0;
	
	double delta_x = 0;
	double ax = 0;
	double Kinetic2 = 0;
	double H1 = 0;
	double H2 = 0;


	/*double xx[b];
	double yy[b];
	double zz[b];
	*/

	int e;
	for (e = 0; e < a * 10; e++){
		flux[e] = 0;
		flux_intensity[e] = 0.0;
		//printf("fluxe\t%d\n", flux[e]);
	}





	int j = 0;
	int i;
	for (i = 0; i < SAMPLE; i++){

		//進行度チェック
		printf("%d\n", (i + 1));




		//初期速度設定
		vx0_data[i] = rand_gauss(gauss_mu, gauss_sigma);//単位[m/s]
		vy0_data[i] = rand_gauss(gauss_mu, gauss_sigma);
		vz0_data[i] = rand_gauss(gauss_mu, gauss_sigma);
		v0[i] = sqrt(vx0_data[i] * vx0_data[i] + vy0_data[i] * vy0_data[i] + vz0_data[i] * vz0_data[i]);

		//初期位置設定
		double r = Uniform()*r0;//単位[m]
		double theta = Uniform()*M_PI;
		double fai = Uniform()*M_PI * 2;
		x0_data[i] = r*sin(theta)*cos(fai)*1e3;//単位[mm]
		y0_data[i] = r*sin(theta)*sin(fai)*1e3;
		z0_data[i] = (r*cos(theta) + HEIGHT)*1e3;
		/*double xf = rand_gauss(gauss_mu, 1e-10);//単位[m]
		double yf = rand_gauss(gauss_mu, 1e-10);
		double zf = rand_gauss(gauss_mu, 1e-10);
		x0_data[i] = xf*1e3;//単位[mm]
		y0_data[i] = yf*1e3;
		z0_data[i] = zf*1e3;*/
		x = x0_data[i] * 1e-3;//単位[m]
		y = y0_data[i] * 1e-3;
		z = z0_data[i] * 1e-3;
		vx = vx0_data[i];//単位[m/s]
		vy = vy0_data[i];
		vz = vz0_data[i];
		int u = 1;
		int ss = 0;

		for (;;){
			loss_jadge = 0;


			//自由落下
			prev_x = x; prev_y = y; prev_z = z;
			prev_vx = vx; prev_vy = vy; prev_vz = vz;
			prev_t = t;


			if (sqrt(x*x + y*y)>r1){

				x += vx*dt_1;
				y += vy*dt_1;
				z += vz*dt_1 - G*dt_1*dt_1 *0.5;
				vz += -G*dt_1;
				t += dt_1;

			}

			else {
				//漏れ励起光による散逸力

				double F, arufa;
				F = h_bar*(2 * M_PI) / RAMUDA*sizenn_haba / 2 * (I_laser(x, y) / (I_laser(x, y) + Is*(1 + (4 * (delta + delta_1) *(delta + delta_1)) / (sizenn_haba*sizenn_haba))));
				arufa = F / MASS;

				//擾乱有の場合の運動 **どちらか消すこと
				/*
				x += vx*dt_1;
				y += vy*dt_1;
				z += vz*dt_1 - (G - arufa)*dt_1*dt_1 / 2.0;
				vz += -(G - arufa)*dt_1;
				t += dt_1;
				*/

				//擾乱無の場合の運動 **どちらか消すこと
			    
				x += vx*dt_1;
				y += vy*dt_1;
				z += vz*dt_1 - G*dt_1*dt_1 *0.5;
				vz += -G*dt_1;
				t += dt_1;
				

			}



			//出射したかどうか判定
			if (z < r1*slope){
				gain_count++;
				flux[int((t / time_with))]++;
				printf("出射out\n");
				//printf("out\tx[um]%f\ty[um]%f\tz[um]%f\t%f\n", x*1e6, y*1e6, z*1e6, sqrt(x*x + y*y)*1e6);

				Kinetic1 = 0.5*MASS*(vx*vx+vy*vy);
				H1 = Kinetic1+Lpotential_U1(x);

				printf("vx[m/s]%f\tvy[m/s]%f\tvz[m/s]%f\t%f\n", vx, vy, vz);
				printf("出射時のx%f\n", x*1e6);
				printf("Kinetic1*10E+28%f\n",Kinetic1*1e28);
			
				

			if(U1max<H1){
				lossjadge++;
				//break;
				printf("U1max超えloss%d\n",lossjadge);
				printf("vx[m/s]%f\tvy[m/s]%f\tvz[m/s]%f\t%f\n", vx, vy, vz);
				printf("Kinetic3*10E+28%f\n",Kinetic1*1e28);
				
				break;

			}else{

				
				while(z < r1*slope){

					/*delta_x = vx*dt_1;

					Kinetic2 = 0.5*MASS*(vx*vx+vy*vy);

					H2 = Kinetic2 + Lpotential_U1(x);

					if(x>0 && vx>0){

						ax = -(Lpotential_U1(x+delta_x)-Lpotential_U1(x))/(sqrt(delta_x*delta_x));
					
					}else if(x>0 && vx<0){

						ax = (Lpotential_U1(x+delta_x)-Lpotential_U1(x))/(sqrt(delta_x*delta_x));

					}else if(x<0 && vx>0){

						ax = -(Lpotential_U1(x+delta_x)-Lpotential_U1(x))/(sqrt(delta_x*delta_x));

					}else if(x<0 && vx<0){

						ax = (Lpotential_U1(x+delta_x)-Lpotential_U1(x))/(sqrt(delta_x*delta_x));

					}
					*/

					//x += vx*dt_1;
					x = x;
				    y += vy*dt_1;
				    z += vz*dt_1 - G*dt_1*dt_1 * 0.5;
				    vz += -G*dt_1;
					//vx += vx+ax*dt_1;

				    t += dt_1;
					
						

					
					
					if(U1max < H1){
						lossjadge++;
						printf("U1(x)超えloss%d\n",lossjadge);
						printf("H1%f\n",H1*1e28);
						printf("H2%f\n",H2*1e28);
						printf("vx%f\n",vx);
				        printf("x%f\n", x*1e6);
						
				        break;

					}else{
						if(z<-0.026){
							gainjadge++;
				            printf("26cm落下でget%d\n",gainjadge);
				            break;

						}else{
							double vr = vr1(x, y, z, vx, vy, vz);
						    double psp;//自然放出が起こる確率
						    psp = 1 - exp((-MASS*vr*gennsui_tyou*sizenn_haba) / (h_bar*delta));
						    //printf("%f\n", psp);

						    if (Uniform() < psp){

								if(U2max<H1){
									lossjadge++;
						            printf("U2maxより高いloss%d\n",lossjadge);
									
				                    break;

								}else{
									gainjadge++;
				                    printf("U2maxより低いget%d\n",gainjadge);
									printf("E2(x+delta_x)-E1(x)/E1を検証%f\n",(Lpotential_U1(x+vx*dt_1)-Lpotential_U1(x))/Lpotential_U1(x));
									
									//printf("Lpotential1_max*10e28%f\n",Lpotential_U1(0.0045)*1e28);
									//printf("Lpotential2_max%f\n",Lpotential_U2(0.0045)*1e28);

									/*
									printf("U1_A%f\n",U1_A(0.0045)*1e25);
									printf("U1_B%f\n",U1_B(0.0045));
									printf("U1_C%f\n",U1_C(0.0045));
									printf("delta%f\n",delta);
									printf("U1_D%f\n",log(1+U1_B(0.0045*U1_C(0.0045))));
									*/

									/*
									printf("U2_A%f\n",U2_A(0.0045)*1e25);
									printf("U2_B%f\n",U2_B(0.0045)*1e4);
									*/
								

				                    break;
								}
							}
						}
					}
				}
			}





				break;
			//	break;
			}


			//原子の位置はどこ？　1.ファネル内（エバ領域ではない）？  2.ファネルの外？  3.エバ領域内？
			int s = 0;
			if ((923/1000)*sqrt(x*x + y*y) > z){
				s = 2;
			}
			else if (slope*sqrt(x*x + y*y) > z - gennsui_tyou*sqrt(slope*slope+1)){
				s = 3;
			}
			else {
				s = 1;
			}


			switch (s)
			{
			case 1:
				break;

			case 2:
				//printf("dtを1/10に\n");
				dt_1 = dt_1*0.1;
				x = prev_x, y = prev_y; z = prev_z;
				vx = prev_vx; vy = prev_vy; vz = prev_vz;
				t = prev_t;
				time_back = 1;
				break;

			case 3:
				dt_1 = dt;
				switch (u)
				{
				case 1:
					u = 1;

					//反射をするかどうか判定
					loss_jadge = reflection_jadge(x, y, z, vx, vy, vz);
					if (loss_jadge == 1){

						//printf("loss1\n");
						break;//switch(u)からのbreak
					}
					else{//弾性反射
						//printf("反射あり1\n");
						double vx2, vy2, vz2;
						vx2 = reflection_vx(x, y, z, vx, vy, vz);
						vy2 = reflection_vy(x, y, z, vx, vy, vz);
						vz2 = reflection_vz(x, y, z, vx, vy, vz);
						//printf("%f\t%f\t%f\t\t%f\t%f\t%f\n",vx,vy,vz,vx2,vy2,vz2);

						vx = vx2;
						vy = vy2;
						vz = vz2;

						vx2 = 0; vy2 = 0; vz2 = 0;//初期化


						//自然放出が起こるか否か
						double vr = vr1(x, y, z, vx, vy, vz);
						double psp;//自然放出が起こる確率
						psp = 1 - exp((-MASS*vr*gennsui_tyou*sizenn_haba) / (h_bar*delta));
						//printf("%f\n", psp);
						if (Uniform() < psp){
							//printf("自然放出発生\n");
							//反跳による加速
							double theta_2 = Uniform()*M_PI;
							double fai_2 = Uniform()*M_PI * 2;
							vx += (2 * M_PI*h_bar / (RAMUDA*MASS))*sin(theta)*cos(fai);
							vy += (2 * M_PI*h_bar / (RAMUDA*MASS))*sin(theta)*sin(fai);
							vz += (2 * M_PI*h_bar / (RAMUDA*MASS))*cos(theta);


							//Sisyphus冷却が起こるかどうか
							double q1 = 13.0 / 18.0;
							double q2 = 5.0 / 18.0;
							double qq = Uniform();
							//printf("q2は%f\n", q2);

							if (qq < q2){

								//printf("Sisyphus冷却発生\n");
								//printf("%f\n",Uniform());
								double vrx = vr2(x, y, z, vx, 0, 0);
								double vry = vr2(x, y, z, 0, vy, 0);
								double vrz = vr2(x, y, z, 0, 0, vz);
								/*
								double vrx_2 = vrx*0.8;
								double vry_2 = vry*0.8;
								double vrz_2 = vrz*0.8;
								*/
								double vrx_2 = vrx*sqrt(delta / (delta + delta_1)*q2 / q1);
								double vry_2 = vry*sqrt(delta / (delta + delta_1)*q2 / q1);
								double vrz_2 = vrz*sqrt(delta / (delta + delta_1)*q2 / q1);

								vx = vrx_2;
								vy = vry_2;
								vz = vrz_2;

							}//Sisyphus冷却終了

						}//自然放出による加熱終了
					}//弾性反射終了

					break;//switch(u)からのbreak




				case 2://case 2はいらないかも？

					//原子の位置はどこ？　1.ファネル内（エバ領域ではない）？  2.ファネルの外？  3.エバ領域内？

					/*if (sqrt(x*x + y*y) > z){
					ss= 2;
					}
					else if (sqrt(x*x + y*y) > z - gennsui_tyou*sqrt(2)){
					ss = 3;
					}
					else {
					ss = 1;
					}

					if (ss==1){
					u = 1;
					}


					u = 1;
					//反射をするかどうか判定
					loss_jadge = reflection_jadge(x, y, z, vx, vy, vz);
					if (loss_jadge == 1){
					printf("loss2\t%f\t%f\t%f\n", x, y, z);
					break;//switch(u)からのbreak
					}
					else{//弾性反射
					printf("反射あり2\n");
					double vx2, vy2, vz2;
					vx2 = reflection_vx(x, y, z, vx, vy, vz);
					vy2 = reflection_vy(x, y, z, vx, vy, vz);
					vz2 = reflection_vz(x, y, z, vx, vy, vz);
					//printf("%f\t%f\t%f\t\t%f\t%f\t%f\n",vx,vy,vz,vx2,vy2,vz2);

					vx = vx2;
					vy = vy2;
					vz = vz2;

					vx2 = 0; vy2 = 0; vz2 = 0;//初期化

					}


					break;//switch(u)からのbreak
					*/


				default:
					break;
				}//switch(u)ここまで
				break;
			default:
				break;//switch(s)からのbreak
			}//switch(s)ここまで


			if (loss_jadge == 1){
				loss_count++;
				printf("ファネルでのloss\n");
				break;//for文からのbreak
			}

			if (t > time_max){
				time_out_count++;
				printf("timeout\n");
				break;//for文からのbreak
			}


			/*
			if (time_back == 0){
			xx[j] = x*1e3;
			yy[j] = y*1e3;
			zz[j] = z*1e3;
			j++;
			}
			else time_back = 0;
			*/

			


			//新しいプログラムここまで
			//
			//
			//
			//




		}//for文ここまで


	//		gaincount=gaincount + gainjadge;
		
	//	    losscount=losscount + lossjadge;

		//fig_3D_xyz_at_fanel(xx, yy, zz, j);

		t = 0;
		/*
		for (j; j = 0; j--){
		xx[j] = 0;
		yy[j] = 0;
		zz[j] = 0;
		}
		*/


	}//for文全sampleについて(i)ここまで

	printf("loss_count\t%f\t\tgain_count\t%f\t\ttime_out_count\t%f\t\tvx\t%f\tvy\t%f\tvz\t%f\n", loss_count, gain_count, time_out_count,vx/gain_count,vy/gain_count,vz/gain_count);
	printf("loss%d\n",lossjadge);
	printf("gain%d\n",gainjadge);



	//fig_3D_xyz0(x0_data, y0_data, z0_data, SAMPLE);
	//fig_laser_3D_to_2Dmap();


	flux_intensity[0] = (flux[0] * MOT / SAMPLE) / (M_PI*(r1*1e2)*(r1*1e2)*time_with);
	double flux_intensity_max = flux_intensity[0];
	int flux_max = flux[0];
	int data;
	for (data = 1; data <a; data++){
		flux_intensity[data] = (flux[data] * MOT / SAMPLE) / (M_PI*(r1*1e2)*(r1*1e2)*time_with);

		if (flux_intensity[data]>flux_intensity_max){
			flux_intensity_max = flux_intensity[data];
		}
		if (flux[data]>flux_max){
			flux_max = flux[data];
		}
	}
	double flux_1s;
	flux_1s = (gain_count * MOT / SAMPLE) / (M_PI*(r1*1e2)*(r1*1e2)*time_max);

	printf("flux_intensity_max[10^12]\t%f\tflux_max\t%d\tflux_1s[10^12]\t%f\n", flux_intensity_max*1e-12, flux_max, flux_1s*1e-12);

	// fig_flux_intensity_2D(flux_intensity, a);

	std::ofstream ofs1;
	std::ostringstream filename;
	filename << "parameter = " << para_name << ".csv";
	ofs1.open(filename.str(), std::ios::app);
	ofs1 << "para = ," << para << ", flux_1s = ," << flux_1s << ",";	// 変わる
	ofs1 << ", loss_count = ," << loss_count << ", time_out_count = ," << time_out_count << ", gain_count = ," << gain_count << std::endl;
	ofs1.close();
	system("pause");
	return 0;

}


double Lpotential_U1(double x){
	double U1;
	U1=((h_bar*delta)/3.0)*log(1.0+((2.0*P*x*x*sizenn_haba*sizenn_haba)/(Is*M_PI*w*w*w*w*(4.0*delta*delta+sizenn_haba*sizenn_haba)))*exp((-2.0*x*x)/(w*w)));
		return U1;
}

double U1_A(double x){
	double U1_A;
	U1_A=((h_bar*delta)/3.0);
		return U1_A;
}

double U1_B(double x){
	double U1_B;
	U1_B=((2.0*P*x*x*sizenn_haba*sizenn_haba)/(Is*M_PI*w*w*w*w*(4.0*delta*delta+sizenn_haba*sizenn_haba)));
		return U1_B;
}

double U1_C(double x){
	double U1_C;
	U1_C=exp((-2.0*x*x)/(w*w));
		return U1_C;
}

double Lpotential_U2(double x){
	double U2;
	U2=((h_bar*(delta+delta_1))/3.0)*log(1.0+((2.0*P*x*x*sizenn_haba*sizenn_haba)/(Is*M_PI*pow(w,4.0)*(4.0*(delta+delta_1)*(4.0*(delta+delta_1)+sizenn_haba*sizenn_haba)))*exp(-(2.0*x*x)/(w*w))));
	return U2;
}

double U2_A(double x){
	double U2_A;
	U2_A=((h_bar*(delta+delta_1))/3.0);
		return U2_A;
}

double U2_B(double x){
	double U2_B;
	U2_B=(2.0*P*x*x*sizenn_haba*sizenn_haba)/(Is*M_PI*pow(w,4.0)*(4.0*(delta+delta_1)*(4.0*(delta+delta_1)+sizenn_haba*sizenn_haba)));
		return U2_B;
}



double Lpotential_U3(double x){
	double U3;
	U3=log(x);
	return U3;
}



void fig_laser_3D_to_2Dmap(void){
	char filename[256];
	printf("励起レーザー強度分布--");
	gets_s(filename);


	//gnuplotの起動・書き込み・終了
	FILE *gp;

	if ((gp = _popen(GNUPLOT_PATH, "w")) == NULL) {	// gnuplotをパイプで起動
		fprintf(stderr, "ファイルが見つかりません %s.", GNUPLOT_PATH);
		exit(EXIT_FAILURE);
	}

	fprintf(gp, "cd 'C:\\Users\\itolab\\Desktop\\研究\\モンテカルロシミュレーション\\関数作成練習\\ConsoleApplication1\\%s'\n", filename);





	fprintf(gp, "set xrange [-%f:%f]\n", w1, w1);
	fprintf(gp, "set yrange [-%f:%f]\n", w1, w1);
	//fprintf(gp, "set zrange [0:%f]\n",Ip*10);

	fprintf(gp, "set xl \"x (mm)\"\n");
	fprintf(gp, "set yl \"y (mm)\"\n");
	fprintf(gp, "set zl \"z (mm)\"\n");
	fprintf(gp, "set grid\n");

	fprintf(gp, "set view equal xy\n");
	fprintf(gp, "set ticslevel 0\n");
	fprintf(gp, "set isosamples 100,100\n");
	fprintf(gp, "set pm3d map\n");

	//fprintf(gp, "splot sqrt(x*x+y*y)\n");
	//fprintf(gp, "replot \"%s\" using 1:2:3\n", filename);
	fprintf(gp, "splot (2 * %f) / (%f*(%f*%f))*(exp(-(2 * ((x - %f)*(x - %f) + y*y)) / (%f*%f)) - exp(-(2 * ((x - %f)*(x - %f) + y*y)) / (%f*%f)))\n", P, M_PI, w1, w1, beam_error_x, beam_error_x, w1, w1, beam_error_x, beam_error_x, w2, w2);
	/*
	fprintf(gp, "set terminal png\n");
	fprintf(gp, "set out \"%s.png\"\n", filename);

	fprintf(gp, "splot \"%s\" using 1:2:3\n", filename);
	*/


	fflush(gp);
	system("pause");
	fprintf(gp, "exit\n");
	_pclose(gp);
}

void fig_3D_xyz0(double x[], double y[], double z[], int n){

	FILE *fp;
	errno_t error;
	char *format = "%-8f %-8f %-8f \n";
	char filename[256];



	//.datファイル（データ元）の作成

	printf("原子集団初期分布（位置）3Dのファイル名(拡張子.datまで入力)--");
	gets_s(filename);



	error = fopen_s(&fp, filename, "wb");
	if (error != 0) {
		perror("ファイルをオープンできません\n");
		return;
	}



	int N;
	for (N = 0; N < n; N++){
		if (fprintf_s(fp, format, x[N], y[N], z[N]) < 0){
			perror("datファイルへの書き込みエラー");
		}
	}


	fclose(fp);






	//gnuplotの起動・書き込み・終了
	FILE *gp;

	if ((gp = _popen(GNUPLOT_PATH, "w")) == NULL) {	// gnuplotをパイプで起動
		fprintf(stderr, "ファイルが見つかりません %s.", GNUPLOT_PATH);
		exit(EXIT_FAILURE);
	}

	fprintf(gp, "cd 'C:\\Users\\itolab\\Desktop\\研究\\モンテカルロシミュレーション\\関数作成練習\\ConsoleApplication1\\%s'\n", filename);





	fprintf(gp, "set xrange [-5:5]\n");
	fprintf(gp, "set yrange [-5:5]\n");
	fprintf(gp, "set zrange [0:5]\n");

	fprintf(gp, "set xl \"x (mm)\"\n");
	fprintf(gp, "set yl \"y (mm)\"\n");
	fprintf(gp, "set zl \"z (mm)\"\n");
	fprintf(gp, "set grid\n");

	fprintf(gp, "set view equal xyz\n");
	fprintf(gp, "set ticslevel 0\n");
	fprintf(gp, "set isosamples 100,100\n");
	//fprintf(gp, "set pm3d map\n");

	fprintf(gp, "splot sqrt(x*x+y*y)\n");
	fprintf(gp, "replot \"%s\" using 1:2:3\n", filename);
	//fprintf(gp, "splot (2 * %f) / (%f*(%f*%f))*(exp(-(2 * ((x - %f)*(x - %f) + y*y)) / (%f*%f)) - exp(-(2 * ((x - %f)*(x - %f) + y*y)) / (%f*%f)))\n", P, M_PI, w1, w1, beam_error_x, beam_error_x, w1, w1, beam_error_x, beam_error_x, w2, w2);
	/*
	fprintf(gp, "set terminal png\n");
	fprintf(gp, "set out \"%s.png\"\n", filename);

	fprintf(gp, "splot \"%s\" using 1:2:3\n", filename);
	*/


	fflush(gp);
	system("pause");
	fprintf(gp, "exit\n");
	_pclose(gp);
}

void fig_3D_xyz_at_fanel(double x[], double y[], double z[], int n){

	FILE *fp;
	errno_t error;
	char *format = "%-8f %-8f %-8f \n";
	char filename[256];



	//.datファイル（データ元）の作成

	printf("原子のファネル壁面分布（位置）3Dグラフのファイル名(拡張子.datまで入力)--");
	gets_s(filename);



	error = fopen_s(&fp, filename, "wb");
	if (error != 0) {
		perror("ファイルをオープンできません\n");
		return;
	}



	int N;
	for (N = 0; N < n; N++){
		if (fprintf_s(fp, format, x[N], y[N], z[N]) < 0){
			perror("datファイルへの書き込みエラー");
		}
	}


	fclose(fp);






	//gnuplotの起動・書き込み・終了
	FILE *gp;

	if ((gp = _popen(GNUPLOT_PATH, "w")) == NULL) {	// gnuplotをパイプで起動
		fprintf(stderr, "ファイルが見つかりません %s.", GNUPLOT_PATH);
		exit(EXIT_FAILURE);
	}

	fprintf(gp, "cd 'C:\\Users\\itolab\\Desktop\\研究\\モンテカルロシミュレーション\\関数作成練習\\ConsoleApplication1\\%s'\n", filename);





	fprintf(gp, "set xrange [-5:5]\n");
	fprintf(gp, "set yrange [-5:5]\n");
	fprintf(gp, "set zrange [-0:5]\n");

	fprintf(gp, "set xl \"x (mm)\"\n");
	fprintf(gp, "set yl \"y (mm)\"\n");
	fprintf(gp, "set zl \"z (mm)\"\n");
	fprintf(gp, "set grid\n");

	fprintf(gp, "set view equal xyz\n");
	fprintf(gp, "set ticslevel 0\n");
	fprintf(gp, "set isosamples 100,100\n");

	fprintf(gp, "splot sqrt(x*x+y*y)\n");

	fprintf(gp, "replot \"%s\" using 1:2:3 w p pointsize 0.001\n", filename);

	/*
	fprintf(gp, "set terminal png\n");
	fprintf(gp, "set out \"%s.png\"\n", filename);

	fprintf(gp, "splot \"%s\" using 1:2:3\n", filename);
	*/


	fflush(gp);
	system("pause");
	fprintf(gp, "exit\n");
	_pclose(gp);
}

void fig_2D_r_U(double x[], double y[], int m){

	FILE *fp_2D;
	errno_t error_2D;
	char *format_2D = "%-8f %-8f \n";
	char filename_2D[256];



	//.datファイル（データ元）の作成

	printf("ファネル衝突１回目のエネルギー分布2Dグラフのファイル名(拡張子.datまで入力)--");
	gets_s(filename_2D);



	error_2D = fopen_s(&fp_2D, filename_2D, "wb");
	if (error_2D != 0) {
		perror("ファイルをオープンできません\n");
		return;
	}



	int M;
	for (M = 0; M < m; M++){
		if (fprintf_s(fp_2D, format_2D, x[M], y[M]) < 0){
			perror("datファイルへの書き込みエラー");
		}
	}


	fclose(fp_2D);






	//gnuplotの起動・書き込み・終了
	FILE *gp_2D;

	if ((gp_2D = _popen(GNUPLOT_PATH, "w")) == NULL) {	// gnuplotをパイプで起動
		fprintf(stderr, "ファイルが見つかりません %s.", GNUPLOT_PATH);
		exit(EXIT_FAILURE);
	}

	fprintf(gp_2D, "cd 'C:\\Users\\itolab\\Desktop\\研究\\モンテカルロシミュレーション\\関数作成練習\\ConsoleApplication1\\%s'\n", filename_2D);





	//fprintf(gp_2D, "set xrange [-5:5]\n");
	//fprintf(gp_2D, "set yrange [-5:5]\n");

	fprintf(gp_2D, "set xl \"開口中心からの距離 (mm)\"\n");
	fprintf(gp_2D, "set yl \"運動エネルギー (mK)\"\n");
	fprintf(gp_2D, "set grid\n");

	fprintf(gp_2D, "set view equal xy\n");
	fprintf(gp_2D, "set ticslevel 0\n");




	fprintf(gp_2D, "plot \"%s\" using 1:2\n", filename_2D);

	/*
	fprintf(gp_2D, "set terminal png\n");
	fprintf(gp_2D, "set out \"%s.png\"\n", filename_2D);

	fprintf(gp_2D, "plot \"%s\" using 1:2\n", filename_2D);
	*/


	fflush(gp_2D);
	system("pause");
	fprintf(gp_2D, "exit\n");
	_pclose(gp_2D);
}

double Uniform(void){
	return ((double)rand() + 1.0) / ((double)RAND_MAX + 2.0);
}

double rand_gauss(double mu, double sigma){
	double z = sqrt(-2.0*log(Uniform())) * sin(2.0*M_PI*Uniform());
	return mu + sigma*z;
}

void fig_BOLTZ_2D(double x[], int o){
	double X_axis_max = 0.5;
	double v0_data[box_number];

	int p;
	for (p = 0; p < box_number; p++){
		v0_data[p] = 0;
	}

	int h;
	for (h = 0; h < o; h++){
		int s9 = int(floor(x[h] / (X_axis_max / box_number)));
		v0_data[s9]++;
	}



	FILE *fp_B;
	errno_t error_B;
	char *format_B = "%-8f %-8f \n";
	char filename_B[256];



	//.datファイル（データ元）の作成

	printf("初期速度v0ボルツマン分布2Dのファイル名(拡張子.datまで入力)--");
	gets_s(filename_B);



	error_B = fopen_s(&fp_B, filename_B, "wb");
	if (error_B != 0) {
		perror("ファイルをオープンできません\n");
		return;
	}



	int M;
	for (M = 0; M < box_number; M++){
		if (fprintf_s(fp_B, format_B, X_axis_max / box_number*M, v0_data[M]) < 0){
			perror("datファイルへの書き込みエラー");
		}
	}


	fclose(fp_B);






	//gnuplotの起動・書き込み・終了
	FILE *gp_B;

	if ((gp_B = _popen(GNUPLOT_PATH, "w")) == NULL) {	// gnuplotをパイプで起動
		fprintf(stderr, "ファイルが見つかりません %s.", GNUPLOT_PATH);
		exit(EXIT_FAILURE);
	}

	fprintf(gp_B, "cd 'C:\\Users\\itolab\\Desktop\\研究\\モンテカルロシミュレーション\\関数作成練習\\ConsoleApplication1\\%s'\n", filename_B);





	fprintf(gp_B, "set xrange [0:%f]\n", X_axis_max);
	//fprintf(gp_B, "set yrange [-5:5]\n");

	fprintf(gp_B, "set xl \"初期速度v0 (m/s)\"\n");
	fprintf(gp_B, "set yl \"atoms number\"\n");
	fprintf(gp_B, "set grid\n");

	fprintf(gp_B, "set view equal xy\n");
	fprintf(gp_B, "set ticslevel 0\n");




	fprintf(gp_B, "plot \"%s\" using 1:2 with boxes\n", filename_B);

	/*
	fprintf(gp_B, "set terminal png\n");
	fprintf(gp_B, "set out \"%s.png\"\n", filename_B);

	fprintf(gp_B, "plot \"%s\" using 1:2\n", filename_B);
	*/


	fflush(gp_B);
	system("pause");
	fprintf(gp_B, "exit\n");
	_pclose(gp_B);



}

void fig_flux_intensity_2D(double x[], int o){
	double X_axis_max = time_max;


	FILE *fp_B;
	errno_t error_B;
	char *format_B = "%-8f %-8f \n";
	char filename_B[256];



	//.datファイル（データ元）の作成

	printf("flux強度2Dのファイル名(拡張子.datまで入力)--");
	gets_s(filename_B);



	error_B = fopen_s(&fp_B, filename_B, "wb");
	if (error_B != 0) {
		perror("ファイルをオープンできません\n");
		return;
	}



	int M;
	for (M = 0; M < a * 10; M++){
		if (fprintf_s(fp_B, format_B, X_axis_max / box_number*M, x[M]) < 0){
			perror("datファイルへの書き込みエラー");
		}
	}


	fclose(fp_B);






	//gnuplotの起動・書き込み・終了
	FILE *gp_B;

	if ((gp_B = _popen(GNUPLOT_PATH, "w")) == NULL) {	// gnuplotをパイプで起動
		fprintf(stderr, "ファイルが見つかりません %s.", GNUPLOT_PATH);
		exit(EXIT_FAILURE);
	}



	fprintf(gp_B, "cd 'C:\\Users\\itolab\\Desktop\\研究\\モンテカルロシミュレーション\\関数作成練習\\ConsoleApplication1\\%s'\n", filename_B);





	fprintf(gp_B, "set xrange [0:%f]\n", X_axis_max);
	//fprintf(gp_B, "set yrange [-5:5]\n");

	fprintf(gp_B, "set xl \"時刻t (s)\"\n");
	fprintf(gp_B, "set yl \"flux intensity ( atoms/(cm^2・s) )\"\n");
	fprintf(gp_B, "set grid\n");

	fprintf(gp_B, "set view equal xy\n");
	fprintf(gp_B, "set ticslevel 0\n");




	fprintf(gp_B, "plot \"%s\" using 1:2 with boxes\n", filename_B);

	/*
	fprintf(gp_B, "set terminal png\n");
	fprintf(gp_B, "set out \"%s.png\"\n", filename_B);

	fprintf(gp_B, "plot \"%s\" using 1:2\n", filename_B);
	*/


	fflush(gp_B);
	system("pause");
	fprintf(gp_B, "exit\n");
	_pclose(gp_B);



}

int reflection_jadge(double x, double y, double z, double vx, double vy, double vz){
	double U_Rb;
	int loss_jadge_1;
	U_Rb = MASS / (2 * (x*x + y*y + z*z))*pow(-x*vx - y*vy + z*vz, 2);//*J_to_K;
	//printf("%f\t%f\n", U_Rb*J_to_K, U_laser(x, y)*J_to_K);
	if (U_Rb>U_laser(x, y)){
		loss_jadge_1 = 1;
	}
	else loss_jadge_1 = 0;
	//printf("UI\t%f\t\tURb\t%f\n", U_laser(x, y)*J_to_K*1e3, U_Rb*J_to_K*1e3);
	return loss_jadge_1;
}

double U_laser(double xa, double ya){

	double Ixy, Ixy1, U_laser_1, U_laser_2;
	//Ixy = (2 * P) / (M_PI*(w1*w1))*(exp(-(2 * ((xa - beam_error_x)*(xa - beam_error_x) + ya*ya)) / (w1*w1)) - exp(-(2 * ((xa - beam_error_x)*(xa - beam_error_x) + ya*ya)) / (w2*w2)));
	//Ixy = (2 * P) / (M_PI*w1_2*w1_2)*(exp(-(2 * (((xa - beam_error_x) + r2*((xa - beam_error_x) / sqrt((xa - beam_error_x)*(xa - beam_error_x) + ya*ya)))*((xa - beam_error_x) + r2*((xa - beam_error_x) / sqrt((xa - beam_error_x)*(xa - beam_error_x) + ya*ya))) + (ya + r2*(ya / sqrt((xa - beam_error_x)*(xa - beam_error_x) + ya*ya)))*(ya + r2*(ya / sqrt((xa - beam_error_x)*(xa - beam_error_x) + ya*ya))))) / (w1_2*w1_2)) - exp(-(2 * (((xa - beam_error_x) + r2*((xa - beam_error_x) / sqrt((xa - beam_error_x)*(xa - beam_error_x) + ya*ya)))*((xa - beam_error_x) + r2*((xa - beam_error_x) / sqrt((xa - beam_error_x)*(xa - beam_error_x) + ya*ya))) + (ya + r2*(ya / sqrt((xa - beam_error_x)*(xa - beam_error_x) + ya*ya)))*(ya + r2*(ya / sqrt((xa - beam_error_x)*(xa - beam_error_x) + ya*ya))))) / (w2_2*w2_2)));
	//Ixy = I_laser_fanel(xa, ya);
	Ixy = I_laser_fanel(xa, ya);

	Ixy1 = Ixy*cappling_rate*Enhance;
	//printf("%f\t%f\n", Ixy, Ixy1);
	double omega = sizenn_haba*sqrt(Ixy1 / (2 * Is));
	double s4 = (omega*omega / 2) / (delta*delta + (sizenn_haba*sizenn_haba / 4));
	U_laser_2 = (h_bar*delta / 2)*(log(1 + s4));


	/*if (xa*xa+ya*ya<r1){
	U_laser_1 = 0;
	}
	else*/ U_laser_1 = U_laser_2;

	return U_laser_1;
}

double I_laser(double xa, double ya){

	double Ixy_laser;
	//Ixy_laser =((2 * 0.5 * P) / (M_PI*(9e-6))*((0.89)*exp(-(2 * ((xa - beam_error_x)*(xa - beam_error_x) + ya*ya)) / (9e-6))*(1 - exp((-2 * ((xa - beam_error_x)*(xa - beam_error_x) + ya*ya)) / (1e-6))) + (0.11)*(2*((xa - beam_error_x)*(xa - beam_error_x) + ya*ya))/(3.6e-6)*exp((-2*((xa - beam_error_x)*(xa - beam_error_x) + ya*ya))/(3.6e-6))));
	Ixy_laser = (2 * P) / (M_PI*(w1*w1))*(2*((xa - beam_error_x)*(xa - beam_error_x) + ya*ya) / (w1*w1))*exp(-2*((xa - beam_error_x)*(xa - beam_error_x) + ya*ya) / (w1*w1));
	//printf("%f\tx\t%f\ty\t%f\n", Ixy_laser,xa,ya);
	return Ixy_laser;
}

double I_laser_fanel(double xa, double ya){
	double Ixy_fanel;
	if (sqrt(xa*xa + ya*ya)>=(r1 + capplar_h) && sqrt(xa*xa + ya*ya)<=(r2 + r1)){
		Ixy_fanel = I_laser(xa, ya);
	}
	else if (sqrt(xa*xa + ya*ya)>=(r1) && sqrt(xa*xa + ya*ya)<(r1 + capplar_h)){
		Ixy_fanel = I_laser(xa*(1 + r2 / sqrt(xa*xa + ya*ya)), ya*(1 + r2 / (sqrt(xa*xa + ya*ya))));
	}
	else Ixy_fanel = 0;
	//printf("%f\n", Ixy_fanel);
	return Ixy_fanel;
}

double reflection_vx(double xa, double ya, double za, double vxa, double vya, double vza){
	double lx, ly, lz, nx, ny, nz, kx, ky, kz, v1;

	lx = xa / sqrt(xa*xa + ya*ya + za*za);
	ly = ya / sqrt(xa*xa + ya*ya + za*za);
	lz = za / sqrt(xa*xa + ya*ya + za*za);

	nx = -lx;
	ny = -ly;
	nz = lz;

	kx = -nx*(vxa*nx + vya*ny + vza*nz);
	ky = -ny*(vxa*nx + vya*ny + vza*nz);
	kz = -nz*(vxa*nx + vya*ny + vza*nz);

	v1 = vxa + 2 * kx;
	return v1;


}

double reflection_vy(double xa, double ya, double za, double vxa, double vya, double vza){
	double lx, ly, lz, nx, ny, nz, kx, ky, kz, v1;

	lx = xa / sqrt(xa*xa + ya*ya + za*za);
	ly = ya / sqrt(xa*xa + ya*ya + za*za);
	lz = za / sqrt(xa*xa + ya*ya + za*za);

	nx = -lx;
	ny = -ly;
	nz = lz;

	kx = -nx*(vxa*nx + vya*ny + vza*nz);
	ky = -ny*(vxa*nx + vya*ny + vza*nz);
	kz = -nz*(vxa*nx + vya*ny + vza*nz);

	v1 = vya + 2 * ky;
	return v1;


}

double reflection_vz(double xa, double ya, double za, double vxa, double vya, double vza){
	double lx, ly, lz, nx, ny, nz, kx, ky, kz, v1;

	lx = xa / sqrt(xa*xa + ya*ya + za*za);
	ly = ya / sqrt(xa*xa + ya*ya + za*za);
	lz = za / sqrt(xa*xa + ya*ya + za*za);

	nx = -lx;
	ny = -ly;
	nz = lz;

	kx = -nx*(vxa*nx + vya*ny + vza*nz);
	ky = -ny*(vxa*nx + vya*ny + vza*nz);
	kz = -nz*(vxa*nx + vya*ny + vza*nz);

	v1 = v1 = vza + 2 * kz;
	return v1;


}

double vr1(double xa, double ya, double za, double vxa, double vya, double vza){
	double lx, ly, lz, nx, ny, nz, kx, ky, kz, vr1;

	lx = (xa - r1) / sqrt((xa - r1)*(xa - r1) + (ya - r1)*(ya - r1) + za*za);
	ly = (ya - r1) / sqrt((xa - r1)*(xa - r1) + (ya - r1)*(ya - r1) + za*za);
	lz = (za) / sqrt((xa - r1)*(xa - r1) + (ya - r1)*(ya - r1) + za*za);

	nx = -lx;
	ny = -ly;
	nz = lz;

	kx = -nx*(vxa*nx + vya*ny + vza*nz);
	ky = -ny*(vxa*nx + vya*ny + vza*nz);
	kz = -nz*(vxa*nx + vya*ny + vza*nz);

	vr1 = sqrt((vxa*nx + vya*ny + vza*nz)*(vxa*nx + vya*ny + vza*nz));
	return vr1;


}

double vr2(double xa, double ya, double za, double vxa, double vya, double vza){
	double lx, ly, lz, nx, ny, nz, kx, ky, kz, vr2;

	lx = (xa - r1) / sqrt((xa - r1)*(xa - r1) + (ya - r1)*(ya - r1) + za*za);
	ly = (ya - r1) / sqrt((xa - r1)*(xa - r1) + (ya - r1)*(ya - r1) + za*za);
	lz = (za) / sqrt((xa - r1)*(xa - r1) + (ya - r1)*(ya - r1) + za*za);

	nx = -lx;
	ny = -ly;
	nz = lz;

	kx = -nx*(vxa*nx + vya*ny + vza*nz);
	ky = -ny*(vxa*nx + vya*ny + vza*nz);
	kz = -nz*(vxa*nx + vya*ny + vza*nz);

	vr2 = (vxa*nx + vya*ny + vza*nz);
	return vr2;


}

