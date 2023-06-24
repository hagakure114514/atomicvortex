// random_kikyaku_method.cpp : このファイルには 'main' 関数が含まれています。プログラム実行の開始と終了がそこで行われます。
//

#include <stdio.h>
#include <cstdio>
#include <corecrt_math_defines.h>
#include <math.h>
#include <random>
#include <iostream>
#include <fstream>

using std::endl;
using std::ofstream;

const double mass = 1.40999e-25;		//mass of Rb : Rb/NA/1000[kg]
const double kb = 1.380649e-23;		//Boltzmann const.[J/K]
const int N = 100;		// csvファイルの生成する数
const double t = 10.0e-6;			// 冷却原子団 templature[K]
const double r0 = 1.0e-3;			//MOT radius[m]
const double sigma = r0 / 3.0;		//standard deviation of normal distribution[m]
const double v_max = 0.15;
const int SAMPLE = 1000; 			//number of sample atoms

void gauss(double* x);
void max_boltz(double* v);


int main()
{
    printf("Random generation start.\n");

    char fname[30];
	sprintf_s(fname, "random_init_pos_vel.csv");
	ofstream ofs(fname);        // ファイルパスを指定する

	for (int j = 1; j <= SAMPLE; j++) {

		double x, vx;
		gauss(&x);
		max_boltz(&vx);

		ofs << x << ", "<< vx << endl;
			
	}

	printf("Random generation finish.\n");
}


void gauss(double* x) {

	std::mt19937 rand_src(std::random_device{}());
	std::uniform_real_distribution<double> dist(0.0, 1.0);

	double p_uni, giji_ransu;
	do {
		p_uni = dist(rand_src);
		giji_ransu = 2.0 * 3.0 * sigma * dist(rand_src) - 3.0 * sigma;

	} while (p_uni > exp(-giji_ransu * giji_ransu / (2.0 * sigma * sigma)));

	*x = giji_ransu;
}


void max_boltz(double* vx) {

	std::mt19937 rand_src(std::random_device{}());
	std::uniform_real_distribution<double> dist(0.0, 1.0);

	double p_uni, giji_ransu;
	do {
		p_uni = dist(rand_src);
		giji_ransu = v_max * dist(rand_src);

	} while (p_uni > mass * giji_ransu * giji_ransu / (2.0 * kb * t) *exp( 1.0 -mass*giji_ransu * giji_ransu / (2.0 *kb*t)));

	*vx = giji_ransu;

}