// F_dip_freq.cpp : このファイルには 'main' 関数が含まれています。プログラム実行の開始と終了がそこで行われます。
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

// const
const double hbar = 1.0545e-34;   // ディラック定数
const double e0 = 8.85418e-12;  // 誘電率
const double c = 2.998e8;   // 光速[m / s]

// atom parameter
const double m = 1.443160e-25;   // 原子質量[kg
const double lambda = 780.2e-9;				//wavelength of detuned D2 line[m]
const double gamma = 38.1e+6;    // 自然幅[Hz]
const double delta_hfs = 2.0 * M_PI * 6.834682610904e+9;				//frequency between the two hyperfine ground states of 87Rb[rad/s]
const double q = 0.75; // 87Rbの励起状態からF = 2への分岐比
const double gamma1 = q * gamma;  // 自然幅 1[Hz]
const double gamma2 = (1 - q) * gamma;    // 自然幅 2[Hz]
const double I_s = 2 * M_PI * M_PI * c * hbar * gamma / (3 * lambda * lambda * lambda);
const double I_s1 = q * I_s;  // 飽和強度[W / m ^ 2]
const double I_s2 = (1 - q) * I_s;  // 飽和強度[W / m ^ 2]

//===Optical Vortex=================================================================================================================================
const int l = 1;				//OAM
const double detuning0 = 2 * M_PI * 1e9;        // 離調 a
const double w0 = 2e-3;			//beam waist[m]
const double beam_power = 100e-3;			//beam power[W]
const double I0 = 2.0 * beam_power / (M_PI * w0 * w0);	//beam intensity coefficient [V/m]
const double k = 2 * M_PI / lambda + detuning0 / c;             //波数
double I(double rh);                             //ビーム強度[W/m^2]

//===Repump beam(Optical Vortex)=================================================================================================================================
const int l_pm = -l;				//OAM
const double detuning0_pm = -2.0 * M_PI * 10e6;			//detuning [rad/s]
const double w0_pm = 1e-3;								//beam waist[m]
const double beam_power_pm = 10e-3;						//beam power[W]
const double I0_pm = 2.0 * beam_power_pm / (M_PI * w0_pm * w0_pm);                      //ビームパワー[W/m^2]
const double k_pm = -2 * M_PI / lambda + (detuning0_pm - delta_hfs) / c;                //波数
double I_pm(double rh);

// sパラメータ
double detuning = 0.0;
double detuning_pm = 0.0;
double s1(double rh);
double s2(double rh);
double grad_s1(double rh);
double grad_s2(double rh);
double s1_pm(double rh);
double s2_pm(double rh);
double grad_s1_pm(double rh);
double grad_s2_pm(double rh);


int main()
{
    
    

    char fname[30];
    sprintf_s(fname, "f_dip_500um_1mm.csv");
    ofstream ofs(fname);        // ファイルパスを指定する

    for (int i = 0; i < 2000; i++) {

        detuning = 2*M_PI*(5e6 * i - 8e9);
        detuning_pm = 2 * M_PI * (5e6 * i - 2e9);

        double r_c = 500e-6;
        // 離調計算
        double detuning1 = detuning;                        //LGビームa 離調 |e>-|g1>
        double detuning2 = detuning + delta_hfs;           //LGビームa 離調 |e>-|g2>
        double detuning1_pm = detuning_pm - delta_hfs;           //LGビームb 離調 |e>-|g1>
        double detuning2_pm = detuning_pm;                        //LGビームb 離調 |e>-|g2>


        // Rabi周波数計算
        double r1 = sqrt(s1(r_c) * (2.0 * detuning1 * detuning1 + 0.5 * gamma1 * gamma1));
        double r2 = sqrt(s2(r_c) * (2.0 * detuning2 * detuning2 + 0.5 * gamma2 * gamma2));
        double gr1 = grad_s1(r_c) * (detuning1 * detuning1 + gamma1 * gamma1 / 4.0) / r1;
        double gr2 = grad_s2(r_c) * (detuning2 * detuning2 + gamma2 * gamma2 / 4.0) / r2;
        double r1_pm = sqrt(s1_pm(r_c) * (2.0 * detuning1_pm * detuning1_pm + 0.5 * gamma1 * gamma1));
        double r2_pm = sqrt(s2_pm(r_c) * (2.0 * detuning2_pm * detuning2_pm + 0.5 * gamma2 * gamma2));
        double gr1_pm = grad_s1_pm(r_c) * (detuning1_pm * detuning1_pm + gamma1 * gamma1 / 4.0) / r1;
        double gr2_pm = grad_s2_pm(r_c) * (detuning2_pm * detuning2_pm + gamma2 * gamma2 / 4.0) / r2;


        //双極子力
 //       double f_dip_500 = 2 * hbar * (gr1 * (detuning1 * r1 * r2 * r2 * gamma2 * gamma2 - (gamma * gamma + detuning2 * detuning2 + (1 + gamma2 / gamma) * r2 * r2) * detuning1 * r1 * gamma1 * gamma)
 //           + gr2 * (detuning2 * r1 * r1 * r2 * gamma1 * (gamma + gamma2) - (gamma * gamma + detuning1 * detuning1 + (1 + gamma1 / gamma) * r1 * r1) * detuning2 * r2 * gamma2 * gamma))
 //           / (gamma * gamma * (gamma * gamma + detuning1 * detuning1 + (1 + gamma1 / gamma) * r1 * r1) * (gamma * gamma + detuning2 * detuning2 + (1 + gamma2 / gamma) * r2 * r2) - r1 * r1 * r2 * r2 * (gamma + gamma1) * (gamma + gamma2));
        //双極子力
        double f_dip_500 = 2 * hbar * (gr1_pm * (detuning1_pm * r1_pm * r2_pm * r2_pm * gamma2 * gamma2 - (gamma * gamma + detuning2_pm * detuning2_pm + (1 + gamma2 / gamma) * r2_pm * r2_pm) * detuning1_pm * r1_pm * gamma1 * gamma)
            + gr2_pm * (detuning2_pm * r1_pm * r1_pm * r2_pm * gamma1 * (gamma + gamma2) - (gamma * gamma + detuning1_pm * detuning1_pm + (1 + gamma1 / gamma) * r1_pm * r1_pm) * detuning2_pm * r2_pm * gamma1 * gamma))
           / (gamma * gamma * (gamma * gamma + detuning1_pm * detuning1_pm + (1 + gamma1 / gamma) * r1_pm * r1_pm) * (gamma * gamma + detuning2_pm * detuning2_pm + (1 + gamma2 / gamma) * r2_pm * r2_pm) - r1_pm * r1_pm * r2_pm * r2_pm * (gamma1 + gamma) * (gamma + gamma2));




        r_c = 1e-3;

        // Rabi周波数計算
         r1 = sqrt(s1(r_c) * (2.0 * detuning1 * detuning1 + 0.5 * gamma1 * gamma1));
         r2 = sqrt(s2(r_c) * (2.0 * detuning2 * detuning2 + 0.5 * gamma2 * gamma2));
         gr1 = grad_s1(r_c) * (detuning1 * detuning1 + gamma1 * gamma1 / 4.0) / r1;
         gr2 = grad_s2(r_c) * (detuning2 * detuning2 + gamma2 * gamma2 / 4.0) / r2;
         r1_pm = sqrt(s1_pm(r_c) * (2.0 * detuning1_pm * detuning1_pm + 0.5 * gamma1 * gamma1));
         r2_pm = sqrt(s2_pm(r_c) * (2.0 * detuning2_pm * detuning2_pm + 0.5 * gamma2 * gamma2));
         gr1_pm = grad_s1_pm(r_c) * (detuning1_pm * detuning1_pm + gamma1 * gamma1 / 4.0) / r1;
         gr2_pm = grad_s2_pm(r_c) * (detuning2_pm * detuning2_pm + gamma2 * gamma2 / 4.0) / r2;



        //双極子力
  //       double f_dip_1000 = 2 * hbar * (gr1 * (detuning1 * r1 * r2 * r2 * gamma2 * gamma2 - (gamma * gamma + detuning2 * detuning2 + (1 + gamma2 / gamma) * r2 * r2) * detuning1 * r1 * gamma1 * gamma)
  //          + gr2 * (detuning2 * r1 * r1 * r2 * gamma1 * (gamma + gamma2) - (gamma * gamma + detuning1 * detuning1 + (1 + gamma1 / gamma) * r1 * r1) * detuning2 * r2 * gamma2 * gamma))
  //          / (gamma * gamma * (gamma * gamma + detuning1 * detuning1 + (1 + gamma1 / gamma) * r1 * r1) * (gamma * gamma + detuning2 * detuning2 + (1 + gamma2 / gamma) * r2 * r2) - r1 * r1 * r2 * r2 * (gamma + gamma1) * (gamma + gamma2));

        double f_dip_1000 = 2 * hbar * (gr1_pm * (detuning1_pm * r1_pm * r2_pm * r2_pm * gamma2 * gamma2 - (gamma * gamma + detuning2_pm * detuning2_pm + (1 + gamma2 / gamma) * r2_pm * r2_pm) * detuning1_pm * r1_pm * gamma1 * gamma)
            + gr2_pm * (detuning2_pm * r1_pm * r1_pm * r2_pm * gamma1 * (gamma + gamma2) - (gamma * gamma + detuning1_pm * detuning1_pm + (1 + gamma1 / gamma) * r1_pm * r1_pm) * detuning2_pm * r2_pm * gamma1 * gamma))
            / (gamma * gamma * (gamma * gamma + detuning1_pm * detuning1_pm + (1 + gamma1 / gamma) * r1_pm * r1_pm) * (gamma * gamma + detuning2_pm * detuning2_pm + (1 + gamma2 / gamma) * r2_pm * r2_pm) - r1_pm * r1_pm * r2_pm * r2_pm * (gamma1 + gamma) * (gamma + gamma2));

        ofs << detuning_pm/ (2 * M_PI) << ", " << f_dip_500 << ", " << f_dip_1000 << endl;

    }

    return 0;

}


// 関数
double I(double rh) {                              //ビーム強度[W/m^2]
    return 2.0 * I0 * rh * rh / (w0 * w0) * exp(-2.0 * rh * rh / (w0 * w0));			// intensity [V/m]
}

double I_pm(double rh) {
    return 2.0 * I0_pm * rh * rh / (w0_pm * w0_pm) * exp(-2.0 * rh * rh / (w0_pm * w0_pm));			// intensity [V/m]
}

// saturation parameter between |g1> and |e>
double s1(double rh)
{
    double s = I(rh) / (I_s1 * (1.0 + 4.0 * detuning * detuning / (gamma1 * gamma1)));

    if (s > 0.1) {
        printf("s parameter1 is error %e, %e\n", s ,detuning / ( 2 * M_PI ) );
    }

    return s;
}

double s2(double rh)
{
    double s = I(rh) / (I_s2 * (1.0 + 4.0 * (detuning + delta_hfs) * (detuning + delta_hfs) / (gamma2 * gamma2)));

    if (s > 0.1) {
        printf("s parameter2 is error %e\n", s);
    }

    return s;
}

double grad_s1(double rh)
{
    double grad_int = 4.0 * I0 / (w0 * w0) * (rh - 2.0 * rh * rh * rh / (w0 * w0)) * exp(-2 * rh * rh / (w0 * w0));
    double gs = grad_int / (I_s1 * (1 + 4 * detuning * detuning / (gamma1 * gamma1)));
    if (gs < 0.0) {
        printf("grad s parameter1 pm is error\n");
    }
    return gs;
}

double grad_s2(double rh)
{
    double grad_int = 4.0 * I0 / (w0 * w0) * (rh - 2.0 * rh * rh * rh / (w0 * w0)) * exp(-2 * rh * rh / (w0 * w0));
    double gs = grad_int / (I_s2 * (1 + 4 * (detuning + delta_hfs) * (detuning + delta_hfs) / (gamma2 * gamma2)));
    if (gs < 0.0) {
        printf("grad s parameter2 pm is error\n");
    }
    return gs;
    return gs;
}

double s1_pm(double rh)
{
    double s = I_pm(rh) / (I_s1 * (1.0 + 4.0 * (detuning_pm - delta_hfs) * (detuning_pm - delta_hfs) / (gamma1 * gamma1)));

    if (s > 0.1) {
        printf("s parameter1 pm is error %e\n", s);
    }

    return s;
}

double s2_pm(double rh)
{
    double s = I_pm(rh) / (I_s2 * (1 + 4 * detuning_pm * detuning_pm / (gamma2 * gamma2)));

    if (s > 0.1) {
        printf("s parameter2 pm is error %e\n", s);
    }

    return s;
}

double grad_s1_pm(double rh)
{
    double grad_int = 4.0 * I0_pm / (w0_pm * w0_pm) * (rh - 2.0 * rh * rh * rh / (w0_pm * w0_pm)) * exp(-2 * rh * rh / (w0_pm * w0_pm));
    double gs = grad_int / (I_s1 * (1 + 4 * (detuning_pm - delta_hfs) * (detuning_pm - delta_hfs) / (gamma1 * gamma1)));
    return gs;
}

double grad_s2_pm(double rh)
{
    double grad_int = 4.0 * I0_pm / (w0_pm * w0_pm) * (rh - 2.0 * rh * rh * rh / (w0_pm * w0_pm)) * exp(-2 * rh * rh / (w0_pm * w0_pm));
    double gs = grad_int / (I_s2 * (1 + 4 * detuning_pm * detuning_pm / (gamma2 * gamma2)));
    return gs;
}