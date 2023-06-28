#pragma once

#include <stdio.h>
#include <corecrt_math_defines.h>

//===’è”=================================================================================================================================
const double hbar = 1.0545e-34;			//Dirac constant[Js]
const double h = 6.626070e-34;			//Planck constant[Js]
const double e0 = 8.85418e-12;			//—U“d—¦
const double c = 2.99792458e+8;			//light speed in vacuum[m/s]
const double mass = 1.40999e-25;		//mass of Rb[kg]
const double k_b = 1.380649e-23;		//Boltzmann const.[J/K]
const double G = 9.80665;				//standard gravity[m/s^2]
const double delta_hfs = 6.834682610904e+9 * 2.0 * M_PI;				//frequency between the two hyperfine ground states of 87Rb[rad/s]
const double branch = 0.75;												//branching ratio into |1>(lower hyperfine ground state) (based on "Gravitational laser trap for atoms with evanescent-wave cooling" pp.654, 661)
const double gamma = 38.1e+6;										//natural linewidth of D2 line[rad/s]

//===Œõ‰Q ƒpƒ‰[ƒ[ƒ^[=================================================================================================================================
const int l = 1;				//•ûˆÊŠpƒ‚[ƒhŽw”
const int prop = -1;			//ƒr[ƒ‚Ìis•ûŒü -1:-z•ûŒü, +1:+z•ûŒü
const double detuning0 = 2.0 * M_PI * 1e9;			//—£’² ƒ¢=ƒÖ-ƒÖ_0 [rad/s]
const double lambda = c / (384.26e12 + detuning0 / (2.0 * M_PI));				//wavelength of detuned D2 line[m]
const double w0 = 2e-3;			//ƒr[ƒ•[m]
const double beam_power = 100e-3;			//ƒr[ƒo—Í[W]


//===ƒŠƒ|ƒ“ƒvŒõ ƒpƒ‰[ƒ[ƒ^[=================================================================================================================================
const double detuning_pm = 2.0 * M_PI * 1e9;			//—£’² [rad/s]
const double lambda_pm = lambda;						//wavelength of detuned D2 line[m]
const double w0_pm = 5e-3;								//ƒr[ƒ•[m]
const double beam_power_pm = 10e-3;						//ƒr[ƒo—Í[W]


//===—â‹pŒ´Žq’c ƒpƒ‰ƒ[ƒ^[=================================================================================================================================
const double r0 = 1.0e-3;												//MOT radius[m]
const double temp = 10.0e-6;										//temperature of MOT[K]
const double v_max = 0.15;

//===ƒVƒ~ƒ…ƒŒ[ƒVƒ‡ƒ“ ƒpƒ‰ƒ[ƒ^[=================================================================================================================================
const int SAMPLE = 1000;											//number of sample atoms[-]
const int jloop = 20000;										//how many times the time t is advanced[-]
const double dt = 2.0e-5;											//interval time[s]


//==\‘¢‘Ì=================================================================================================================================
enum class state { d1, d2, d3 };			// |d1>=|g1,n>, |d2>=|g2,n>, |d3>=|e, n-1>

typedef struct {
	double x; double y; double z;
}position;

typedef struct {
	double vx; double vy; double vz;
}velocity;


//==kikyakuhou function=================================================================================================================================
void gauss(double* x);
void max_boltz(double* v);
void rm_position(double* x, double* y, double* z);
void rm_velocity(double* vx, double* vy, double* vz);
