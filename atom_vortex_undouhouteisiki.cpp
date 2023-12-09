// LG_cp_flux.cpp : このファイルには 'main' 関数が含まれています。プログラム実行の開始と終了がそこで行われます。
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
const double k_b = 1.380649e-23;
const double e0 = 8.85418e-12;  // 誘電率
const double g = 9.83065;   // 重力加速度[m / s ^ 2]
const double c = 2.998e8;   // 光速[m / s]
const double dt = 5.0e-5; // 時間ステップ[s]
const double temp = 10.0e-6;										//temperature of MOT[K]
const double r0 = 1.0e-3;           //MOT radius[m]
const double sigma = r0 / 3.0;      //standard deviation of normal distribution[m]
const int SAMPLE = 1000;           //number of sample atoms

// atom parameter
const double m = 1.443160e-25;   // 原子質量[kg
const double lambda = 780.2e-9;				//wavelength of detuned D2 line[m]
const double gamma = 38.1e+6;    // 自然幅[Hz]
const double delta_hfs = 2.0 * M_PI* 6.834682610904e+9;				//frequency between the two hyperfine ground states of 87Rb[rad/s]
const double q = 0.75; // 87Rbの励起状態からF = 2への分岐比
const double gamma1 = q * gamma;  // 自然幅 1[Hz]
const double gamma2 = (1 - q) * gamma;    // 自然幅 2[Hz]
const double I_s = 2 * M_PI * M_PI * c * hbar * gamma / (3 * lambda * lambda * lambda);
const double I_s1 = q*I_s;  // 飽和強度[W / m ^ 2]
const double I_s2 = (1-q)*I_s;  // 飽和強度[W / m ^ 2]

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


// 原子の初期位置、初速度　
const double v_max = 0.15;
void gauss(double* x);
void max_boltz(double* v);
void rm_position(double* x, double* y, double* z);
void rm_velocity(double* vx, double* vy, double* vz);



int main()
{    

    int flag_mode=0;

    if (flag_mode == 1) {
  
    }
    else {

        // 宣言
        //double x[10001], y[10001], z[10001];
        int count_trap = 0, count_cw = 0;
        double x0[SAMPLE + 1] = {}, y0[SAMPLE + 1] = {}, z0[SAMPLE + 1] = {};
        double vx0[SAMPLE + 1] = {}, vy0[SAMPLE + 1] = {}, vz0[SAMPLE + 1] = {};

        rm_position(x0, y0, z0);
        rm_velocity(vx0, vy0, vz0);

        
        for (int j = 0; j < SAMPLE; j++) {
            double x = x0[j];   double y = y0[j];   double z = z0[j];         // 初期位置[m]
            double Vx = vx0[j]; double Vy = vy0[j]; double Vz = vz0[j];        // 初期速度[m / s]
            double a_x = 0.0, a_y = 0.0, a_z = 0.0;
            double v_phi;

            for (int i = 0; i < 100000; i++) {
                double r_c = sqrt(x * x + y * y);
                double phi = atan2(y, x);
                v_phi = -Vx * sin(phi) + Vy * cos(phi);

                    // 離調計算
                detuning = detuning0 - (k * Vz + l / r_c * v_phi);
                detuning_pm = detuning0_pm - (k_pm * Vz - l_pm / r_c * v_phi);
                double detuning1 = detuning;                        //LGビームa 離調 |e>-|g1>
                double detuning2 = detuning + delta_hfs;           //LGビームa 離調 |e>-|g2>
                double detuning1_pm = detuning_pm - delta_hfs;           //LGビームb 離調 |e>-|g1>
                double detuning2_pm = detuning_pm;                        //LGビームb                      離調 |e>-

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
                double f_dip =2 * hbar * (gr1 * (detuning1 * r1 * r2 * r2 * gamma2 * gamma2 - (gamma * gamma + detuning2 * detuning2 + (1 + gamma2 / gamma) * r2 * r2) * detuning1 * r1 * gamma1 * gamma)
                                        + gr2 * (detuning2 * r1 * r1 * r2 * gamma1 * (gamma + gamma2) - (gamma * gamma + detuning1 * detuning1 + (1 + gamma1 / gamma) * r1 * r1) * detuning2 * r2 * gamma2 * gamma))
                                    / (gamma * gamma * (gamma * gamma + detuning1 * detuning1 + (1 + gamma1 / gamma) * r1 * r1) * (gamma * gamma + detuning2 * detuning2 + (1 + gamma2 / gamma) * r2 * r2) - r1 * r1 * r2 * r2 * (gamma + gamma1) * (gamma + gamma2))
                                +2 * hbar * (gr1_pm * (detuning1_pm * r1_pm * r2_pm * r2_pm * gamma2 * gamma2 - (gamma * gamma + detuning2_pm * detuning2_pm + (1 + gamma2 / gamma) * r2_pm * r2_pm) * detuning1_pm * r1_pm * gamma1 * gamma)
                                        + gr2_pm * (detuning2_pm * r1_pm * r1_pm * r2_pm * gamma1 * (gamma + gamma2) - (gamma * gamma + detuning1_pm * detuning1_pm + (1 + gamma1 / gamma) * r1_pm * r1_pm) * detuning2_pm * r2_pm * gamma1 * gamma))
                                    / (gamma * gamma * (gamma * gamma + detuning1_pm * detuning1_pm + (1 + gamma1 / gamma) * r1_pm * r1_pm) * (gamma * gamma + detuning2_pm * detuning2_pm + (1 + gamma2 / gamma) * r2_pm * r2_pm) - r1_pm * r1_pm * r2_pm * r2_pm * (gamma1 + gamma) * (gamma + gamma2));

                double f_dip_x = f_dip * cos(phi);
                double f_dip_y = f_dip * sin(phi);


                //自発力

                double f_sp = 2 * hbar * k * gamma * (-r1 * r1 * r2 * r2 * (gamma * gamma1 + gamma * gamma2 + gamma2 * gamma2) + (gamma*gamma+detuning2*detuning2+(1+gamma2/gamma)*gamma2*gamma2)*r1*r1*gamma1*gamma + (gamma*gamma+detuning1*detuning1+(1+gamma1/gamma)*r1*r1)*r2*r2*gamma2*gamma)
                                    / (gamma * gamma * (gamma * gamma + detuning1 * detuning1 + (1 + gamma1 / gamma) * r1 * r1) * (gamma * gamma + detuning2 * detuning2 + (1 + gamma2 / gamma) * r2 * r2) - r1 * r1 * r2 * r2 * (gamma + gamma1) * (gamma + gamma2))
                                + 2 * hbar * k_pm * gamma * (-r1_pm * r1_pm * r2_pm * r2_pm * (gamma * gamma2 + gamma * gamma1 + gamma1 * gamma1) + (gamma*gamma+detuning1_pm*detuning1_pm+(1+gamma1/gamma)*gamma1*gamma1)*r2_pm*r2_pm*gamma2*gamma + (gamma*gamma+detuning2_pm*detuning2_pm+(1+gamma2/gamma)*r2_pm*r2_pm)*r1_pm*r1_pm*gamma1*gamma)
                                    / (gamma * gamma * (gamma * gamma + detuning1_pm * detuning1_pm + (1 + gamma1 / gamma) * r1_pm * r1_pm) * (gamma * gamma + detuning2_pm * detuning2_pm + (1 + gamma2 / gamma) * r2_pm * r2_pm) - r1_pm * r1_pm * r2_pm * r2_pm * (gamma1 + gamma) * (gamma + gamma2));
                
                //トルク
                double f_trq = 2 * hbar * (double)l * gamma / r_c * (-r1 * r1 * r2 * r2 * (gamma * gamma1 + gamma * gamma2 + gamma2 * gamma2) + (gamma*gamma+detuning2*detuning2+(1+gamma2/gamma)*gamma2*gamma2)*r1*r1*gamma1*gamma + (gamma*gamma+detuning1*detuning1+(1+gamma1/gamma)*r1*r1)*r2*r2*gamma2*gamma)
                                    / (gamma * gamma * (gamma * gamma + detuning1 * detuning1 + (1 + gamma1 / gamma) * r1 * r1) * (gamma * gamma + detuning2 * detuning2 + (1 + gamma2 / gamma) * r2 * r2) - r1 * r1 * r2 * r2 * (gamma + gamma1) * (gamma + gamma2))
                                - 2 * hbar * (double)l_pm * gamma / r_c * (-r1_pm * r1_pm * r2_pm * r2_pm * (gamma * gamma2 + gamma * gamma1 + gamma1 * gamma1) + (gamma * gamma + detuning1_pm * detuning1_pm + (1 + gamma1 / gamma) * gamma1 * gamma1) * r2_pm * r2_pm * gamma2 * gamma + (gamma * gamma + detuning2_pm * detuning2_pm + (1 + gamma2 / gamma) * r2_pm * r2_pm) * r1_pm * r1_pm * gamma1 * gamma)
                                    / (gamma * gamma * (gamma * gamma + detuning1_pm * detuning1_pm + (1 + gamma1 / gamma) * r1_pm * r1_pm) * (gamma * gamma + detuning2_pm * detuning2_pm + (1 + gamma2 / gamma) * r2_pm * r2_pm) - r1_pm * r1_pm * r2_pm * r2_pm * (gamma1 + gamma) * (gamma + gamma2));

                double f_trq_x = -f_trq * sin(phi);
                double f_trq_y = f_trq * cos(phi);



                //　離散化運動方程式
                a_x = 1 / m * (f_dip_x + f_trq_x); a_y = 1 / m * (f_dip_y + f_trq_y); a_z = 1 / m * (f_sp)-g;
                x += Vx * dt + 0.5 * a_x * dt * dt; y += Vy * dt + 0.5 * a_y * dt * dt; z += Vz * dt + 0.5 * a_z * dt * dt;
                Vx += a_x * dt; Vy += a_y * dt; Vz += a_z * dt;


                    //      printf("%e, %e, %e\n", x[i], y[i], z[i]);

                if (r_c > w0 / sqrt(2.0))   break;

                if (z < -2.6e-1) {
                    count_trap++;
                    break;
                }
            }

            if (v_phi > 0)  count_cw++;
        }

        printf("trap efficiency[%%]:%e\n", (double)count_trap / (double)SAMPLE * 1e2);
        printf("rotation direction[%%]:%e\n", (2*(double)count_cw / (double)count_trap -1.0 )* 1e2);

        return 0;


    }
}


// 関数
double I(double rh) {                              //ビーム強度[W/m^2]
    return 2.0 * I0 * rh * rh / (w0 * w0) * exp(-2.0 * rh * rh / (w0 * w0));            // intensity [V/m]
}

double I_pm(double rh) {
    return 2.0 * I0_pm * rh * rh / (w0_pm * w0_pm) * exp(-2.0 * rh * rh / (w0_pm * w0_pm));         // intensity [V/m]
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


// functions for the Reject method
void gauss(double* x) {
    std::mt19937 rand_src(std::random_device{}());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    double p_uni, giji_ransu;
    do {
        p_uni = dist(rand_src);
        giji_ransu = 2.0 * 3.0 * r0 / 3.0 * dist(rand_src) - 3.0 * r0 / 3.0;
    } while (p_uni > exp(-giji_ransu * giji_ransu / (2.0 * r0 / 3.0 * r0 / 3.0)));
    *x = giji_ransu;
}

void max_boltz(double* vx) {
    std::mt19937 rand_src(std::random_device{}());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    double p_uni, giji_ransu;
    do {
        p_uni = dist(rand_src);
        giji_ransu = v_max * dist(rand_src);
    } while (p_uni > m * giji_ransu * giji_ransu / (2.0 * k_b * temp) * exp(1.0 - m * giji_ransu * giji_ransu / (2.0 * k_b * temp)));
    *vx = giji_ransu;
}

void rm_position(double* x, double* y, double* z) {
    std::mt19937 rand_src(std::random_device{}());
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    for (int j = 0; j < SAMPLE; j++) {
        double r0;
        gauss(&r0);

        double phi = 2.0 * M_PI * dist(rand_src);
        double psi = M_PI * dist(rand_src);

        x[j] = r0 * cos(phi) * sin(psi);
        y[j] = r0 * sin(phi) * sin(psi);
        z[j] = r0 * cos(psi);
    }
}

void rm_velocity(double* vx, double* vy, double* vz) {
    std::mt19937 rand_src(std::random_device{}());
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    for (int j = 0; j < SAMPLE; j++) {
        double v0;
        max_boltz(&v0);

        double phi = 2.0 * M_PI * dist(rand_src);
        double psi = M_PI * dist(rand_src);

        vx[j] = v0 * cos(phi) * sin(psi);
        vy[j] = v0 * sin(phi) * sin(psi);
        vz[j] = -v0 * abs(cos(psi));
    }
}