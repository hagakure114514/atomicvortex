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
const double hb = 1.0545e-34;   // ディラック定数
const double k_b = 1.380649e-23;
const double e0 = 8.85418e-12;  // 誘電率
const double g = 9.83065;   // 重力加速度[m / s ^ 2]
const double c = 2.998e8;   // 光速[m / s]
const double dt = 2.0e-5; // 時間ステップ[s]
const double t = 10.0e-6;			// 冷却原子団 templature[K]
const double r0 = 1.0e-3;			//MOT radius[m]
const double sigma = r0 / 3.0;		//standard deviation of normal distribution[m]
const double v_max = 0.15;
const int SAMPLE = 10000; 			//number of sample atoms

// atom parameter
const double m = 1.40999e-25;   // 原子質量[kg]
const double lambda = 780.241e-9;   // 波長[m]
const double nw = 2 * M_PI * 6.0666e6;    // 自然幅[Hz]
const double q = 0.741; // 85Rbの励起状態からF = 2への分岐比
const double nw1 = q * nw;  // 自然幅 1[Hz]
const double nw2 = (1 - q) * nw;    // 自然幅 2[Hz]
const double Isat = 2 * M_PI * M_PI * c * hb * nw / (3 * lambda * lambda * lambda);
const double Isat1 = 2 * M_PI * M_PI * c * hb * nw1 / (3 * lambda * lambda * lambda);  // 飽和強度[W / m ^ 2]
const double Isat2 = 2 * M_PI * M_PI * c * hb * nw2 / (3 * lambda * lambda * lambda);  // 飽和強度[W / m ^ 2]

// 光渦
const int l = -1;       //方位角方向
const int p = 0;        //動径方向
const double w0 = 2e-3;    //ビーム幅[m]
const double P = 300.0e-3;            // 出力[W]

// 光ビームa(+z方向)
const double det0_a = 2 * M_PI * 1e9;        // 離調 a
const double I0_a = 2.0 * P / (M_PI * w0 * w0);       //ビームパワー[W/m^2]
const double ka = 2 * M_PI / lambda + det0_a / c;             //波数
double E0a = sqrt((2 * I0_a) / (e0 * c));
double Ea(double rh);                                                     //電界[V/m]
double Ia(double rh);                             //ビーム強度[W/m^2]

// 光ビームb(-z方向)
const double det0_b = 2 * M_PI * (1e9 + 3.035e9);        // 離調 b
const double I0_b = 2.0 * P / (M_PI * w0 * w0);                      //ビームパワー[W/m^2]
const double kb = -2 * M_PI / lambda + det0_b / c;             //波数
const double E0b = sqrt((2 * I0_b) / (e0 * c));
double Eb(double rh);
double Ib(double rh);

// Rabi周波数
double rabi1_a(double rh);
double rabi2_a(double rh);
double rabi1_b(double rh);
double rabi2_b(double rh);


// Rabi周波数 grad
double grad_r1a(double rh);
double grad_r2a(double rh);
double grad_r1b(double rh);
double grad_r2b(double rh);


// 原子の初期位置、初速度　生成
void gauss(double* x, double* y);
void max_boltz(double* vx, double* vy, double* vz);
void generation_init_state(double* x0, double* y0, double* vx0, double* vy0, double* vz0);


int main()
{
    // Optical force on/off
    printf("Optical Force calculate? 0:yes, 1:no \n");
    char c = std::cin.get();

    // 宣言
    double x[2], y[2], z[2];
    double x0[SAMPLE], y0[SAMPLE], vx0[SAMPLE], vy0[SAMPLE], vz0[SAMPLE];
    double x1[SAMPLE], z1[SAMPLE], x2[SAMPLE], z2[SAMPLE], x3[SAMPLE], z3[SAMPLE];
    int count_trap = 0, count_cw = 0;
    double sum_vel2 = 0.0;

    generation_init_state(x0, y0, vx0, vy0, vz0);

    for (int j = 0; j < SAMPLE; j++) {
        x[0] = x0[j];   y[0] = y0[j];   z[0] = 0.0;         // 初期位置[m]
        double Vx = vx0[j]; double Vy = vy0[j]; double Vz = vz0[j];        // 初期速度[m / s]

        double a_x = 0.0, a_y = 0.0, a_z = 0.0;
        double r1a = 0.0, r1b = 0.0, r2a = 0.0, r2b = 0.0, gr1a = 0.0, gr2a = 0.0, gr1b = 0.0, gr2b = 0.0;
        double v_phi;


        for (int i = 0; i < 10000; i++) {
            double r_c = sqrt(x[i] * x[i] + y[i] * y[i]);
            double phi = atan2(y[i], x[i]);
            v_phi = -Vx * sin(phi) + Vy * cos(phi);

            // 離調計算
            double det1a = det0_a - (ka * Vz + l / r_c * v_phi);                        //LGビームa 離調 |e>-|g1>
            double det2a = det0_a - (ka * Vz + l / r_c * v_phi) + 2 * M_PI * 3.035e9;           //LGビームa 離調 |e>-|g2>
            double det1b = det0_b - (kb * Vz + l / r_c * v_phi) - 2 * M_PI * 3.035e9;           //LGビームb 離調 |e>-|g1>
            double det2b = det0_b - (kb * Vz + l / r_c * v_phi);                        //LGビームb 離調 |e>-|g2>

            // Rabi周波数計算
            r1a = rabi1_a(r_c); r1b = rabi1_b(r_c); r2a = rabi2_a(r_c); r2b = rabi2_b(r_c);
            gr1a = grad_r1a(r_c); gr2a = grad_r2a(r_c); gr1b = grad_r1b(r_c); gr2b = grad_r2b(r_c);


            //双極子力
            double f_dip = -2 * hb * (((-det1a * r1a * r2a * r2a * nw2 * nw2 + (nw * nw + det2a * det2a + (2 - q) * r2a * r2a) * det1a * r1a * nw1 * nw)
                / (nw * nw * (nw * nw + det1a * det1a + (1 + q) * r1a * r1a) * (nw * nw + det2a * det2a + (2 - q) * r2a * r2a) - r1a * r1a * r2a * r2a * (2 - q) * nw * nw2) * gr1a
                + ((-det2a * r2a * r1a * r1a * nw1 * (2 - q) * nw + (nw * nw + det1a * det1a + (1 + q) * r1a * r1a) * det2a * r2a * nw2 * nw)
                    / (nw * nw * (nw * nw + det1a * det1a + (1 + q) * r1a * r1a) * (nw * nw + det2a * det2a + (2 - q) * r2a * r2a) - r1a * r1a * r2a * r2a * (2 - q) * nw * nw2)) * gr2a)
                + ((-det1b * r1b * r2b * r2b * nw2 * nw2 + (nw * nw + det2b * det2b + (2 - q) * r2b * r2b) * det1b * r1b * nw1 * nw)
                    / (nw * nw * (nw * nw + det1b * det1b + (1 + q) * r1b * r1b) * (nw * nw + det2b * det2b + (2 - q) * r2b * r2b) - r1b * r1b * r2b * r2b * (2 - q) * nw * nw2) * gr1b
                    + ((-det2b * r2b * r1b * r1b * nw1 * (2 - q) * nw + (nw * nw + det1b * det1b + (1 + q) * r1b * r1b) * det2b * r2b * nw2 * nw)
                        / (nw * nw * (nw * nw + det1b * det1b + (1 + q) * r1b * r1b) * (nw * nw + det2b * det2b + (2 - q) * r2b * r2b) - r1b * r1b * r2b * r2b * (2 - q) * nw * nw2)) * gr2b));

            double f_dip_x = f_dip * cos(phi);
            double f_dip_y = f_dip * sin(phi);


            //自発力

            double f_sp = 2 * hb * (((-(r1a * r1a * r2a * r2a * (nw1 * nw + nw1 * nw2 + nw2 * nw2)) / nw + (nw * nw + det2a * det2a + (2 - q) * r2a * r2a) * r1a * r1a * nw1 + (nw * nw + det1a * det1a + (1 + q) * r1a * r1a) * r2a * r2a * nw2)
                / (nw * (nw * nw + det1a * det1a + (1 + q) * r1a * r1a) * (nw * nw + det2a * det2a + (2 - q) * r2a * r2a) - (r2a * r2a * r1a * r1a * (nw + nw2) * nw2) / nw)) * nw * ka
                + ((-(r1b * r1b * r2b * r2b * (nw1 * nw + nw1 * nw2 + nw2 * nw2)) / nw + (nw * nw + det2b * det2b + (2 - q) * r2b * r2b) * r1b * r1b * nw1 + (nw * nw + det1b * det1b + (1 + q) * r1b * r1b) * r2b * r2b * nw2)
                    / (nw * (nw * nw + det1b * det1b + (1 + q) * r1b * r1b) * (nw * nw + det2b * det2b + (2 - q) * r2b * r2b) - (r2b * r2b * r1b * r1b * (nw + nw2) * nw2) / nw)) * nw * kb);


            //トルク

            double f_trq = 2 * hb * (((-(r1a * r1a * r2a * r2a * (nw1 * nw + nw1 * nw2 + nw2 * nw2)) / nw + (nw * nw + det2a * det2a + (2 - q) * r2a * r2a) * r1a * r1a * nw1 + (nw * nw + det1a * det1a + (1 + q) * r1a * r1a) * r2a * r2a * nw2)
                / (nw * (nw * nw + det1a * det1a + (1 + q) * r1a * r1a) * (nw * nw + det2a * det2a + (2 - q) * r2a * r2a) - (r2a * r2a * r1a * r1a * (nw + nw2) * nw2) / nw))
                + ((-(r1b * r1b * r2b * r2b * (nw1 * nw + nw1 * nw2 + nw2 * nw2)) / nw + (nw * nw + det2b * det2b + (2 - q) * r2b * r2b) * r1b * r1b * nw1 + (nw * nw + det1b * det1b + (1 + q) * r1b * r1b) * r2b * r2b * nw2)
                    / (nw * (nw * nw + det1b * det1b + (1 + q) * r1b * r1b) * (nw * nw + det2b * det2b + (2 - q) * r2b * r2b) - (r2b * r2b * r1b * r1b * (nw + nw2) * nw2) / nw))) * (double)l * nw / r_c;

            double f_trq_x = -f_trq * sin(phi);
            double f_trq_y = f_trq * cos(phi);

            if (c == '1') {
                f_dip_x = 0; f_dip_y = 0; f_sp = 0; f_trq_x = 0; f_trq_y = 0;
            }

            //　離散化運動方程式
            a_x = 1 / m * (f_dip_x + f_trq_x); a_y = 1 / m * (f_dip_y + f_trq_y); a_z = 1 / m * (f_sp)-g;
            x[1] = Vx * dt + 0.5 * a_x * dt * dt + x[0]; y[1] = Vy * dt + 0.5 * a_y * dt * dt + y[0]; z[1] = Vz * dt + 0.5 * a_z * dt * dt + z[0];
            Vx += a_x * dt; Vy += a_y * dt; Vz += a_z * dt;

            x[0] = x[1];  y[0] = y[1]; z[0] = z[1];

            if (i == 3000) {
                x1[j] = x[1]; z1[j] = z[1];
            }

            if (i == 6000) {
                x2[j] = x[1]; z2[j] = z[1];
            }

            if (i == 9000) {
                x3[j] = x[1]; z3[j] = z[1];
            }
        }
    }

    char fname[30];
    sprintf_s(fname, "traj_1.csv");
    ofstream ofs(fname);        // ファイルパスを指定する
    for (int j = 0; j < SAMPLE; j++) {
        ofs << x1[j] << ", " << z1[j] << ", " << x2[j] << ", " << z2[j] << ", " << x3[j] << ", " << z3[j] << endl;
    }

    printf("Output position finish\n");
    return 0;
}

double Ea(double rh) {                                                      //電界[V/m]
    return E0a * exp(-rh * rh / (w0 * w0)) * sqrt(2) * rh / w0;
}

double Ia(double rh) {                              //ビーム強度[W/m^2]
    return 0.5 * e0 * c * Ea(rh) * Ea(rh);
}


double Eb(double rh) {
    return E0b * exp(-rh * rh / (w0 * w0)) * sqrt(2) * rh / w0;
}
double Ib(double rh) {
    return 0.5 * e0 * c * Eb(rh) * Eb(rh);
}


double rabi1_a(double rh) {
    return nw1 * sqrt(Ia(rh) / (2.0 * Isat1));
}
double rabi2_a(double rh) {
    return nw2 * sqrt(Ia(rh) / (2.0 * Isat2));
}
double rabi1_b(double rh) {
    return nw1 * sqrt(Ib(rh) / (2.0 * Isat1));
}
double rabi2_b(double rh) {
    return nw2 * sqrt(Ib(rh) / (2.0 * Isat2));
}

double grad_r1a(double rh) {
    return nw1 * sqrt(I0_a / Isat1) / w0 * (1 - 2 * rh * rh / (w0 * w0)) * exp(-rh * rh / (w0 * w0));
}
double grad_r2a(double rh) {
    return nw2 * sqrt(I0_a / Isat2) / w0 * (1 - 2 * rh * rh / (w0 * w0)) * exp(-rh * rh / (w0 * w0));
}
double grad_r1b(double rh) {
    return nw1 * sqrt(I0_b / Isat1) / w0 * (1 - 2 * rh * rh / (w0 * w0)) * exp(-rh * rh / (w0 * w0));
}
double grad_r2b(double rh) {
    return nw2 * sqrt(I0_b / Isat2) / w0 * (1 - 2 * rh * rh / (w0 * w0)) * exp(-rh * rh / (w0 * w0));
}


void gauss(double* x, double* y) {

    std::mt19937 rand_src(std::random_device{}());
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    double p_uni, giji_ransu;
    do {
        p_uni = dist(rand_src);
        giji_ransu = 2.0 * 3.0 * sigma * dist(rand_src) - 3.0 * sigma;

    } while (p_uni > exp(-giji_ransu * giji_ransu / (2.0 * sigma * sigma)));

    double phi = 2.0 * M_PI * dist(rand_src);

    *x = giji_ransu * cos(phi); *y = giji_ransu * sin(phi);
}

void max_boltz(double* vx, double* vy, double* vz) {

    std::mt19937 rand_src(std::random_device{}());
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    double p_uni, giji_ransu;
    do {
        p_uni = dist(rand_src);
        giji_ransu = v_max * dist(rand_src);

    } while (p_uni > m * giji_ransu * giji_ransu / (2.0 * k_b * t) * exp(1.0 - m * giji_ransu * giji_ransu / (2.0 * k_b * t)));

    double phi = 2.0 * M_PI * dist(rand_src);
    double psi = M_PI * dist(rand_src);

    *vx = giji_ransu * cos(phi) * sin(psi); *vy = giji_ransu * sin(phi) * sin(psi); *vz = -giji_ransu * cos(psi);

}

void generation_init_state(double* x0, double* y0, double* vx0, double* vy0, double* vz0) {
    for (int j = 0; j < SAMPLE; j++) {
        double x, y, vx, vy, vz;

        gauss(&x, &y);
        max_boltz(&vx, &vy, &vz);

        x0[j] = x, y0[j] = y, vx0[j] = vx, vy0[j] = vy, vz0[j] = vz;
    }
}