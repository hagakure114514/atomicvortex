// atomic_vortex.cpp : このファイルには 'main' 関数が含まれています。プログラム実行の開始と終了がそこで行われます。
//

#include <iostream>
#include "const_param.h"
#include "atom.h"
#include "DressedAtom.h"
#include <fstream>

using std::endl;
using std::ofstream;

int main()
{
	int flag_sp_tmp = 3;

	// spontaneous emission mode  
	printf("select spontaneous emission mode?\n 1:no OAM, 2:dipole-radiation, 3:all direction, 4:all direction(coherent OAM rad.)\n");
	std::cin >> flag_sp_tmp;

	
	   	// 変数の定義
	   	int sum_sp=0;
	   	int count_guide = 0;
	   	int count_vphi = 0;

	    // 配列の定義
	    double x0[SAMPLE + 1] = {}, y0[SAMPLE + 1] = {}, z0[SAMPLE + 1] = {};
	    double vx0[SAMPLE + 1] = {}, vy0[SAMPLE + 1] = {}, vz0[SAMPLE + 1] = {};

	    rm_position(x0, y0, z0);
		rm_velocity(vx0, vy0, vz0);

		printf("simulation execution\n");

		char fname[30];
		sprintf_s(fname, "stat_azimuthal_vel.csv");
		ofstream ofs(fname);        // ファイルパスを指定する

		for (int ii = 0; ii < SAMPLE; ii++) {

			position r0 = { x0[ii], y0[ii], z0[ii] };
			velocity v0 = { vx0[ii], vy0[ii], vz0[ii] };

			// オブジェクトのコンストラクタ
			atom Rb87(r0, v0, state::d1);		// 原子オブジェクト
			atom* rb87 = &Rb87;
			DressedAtom OV1;			// dressed-atom状態オブジェクト
			OV1.flag_sp = flag_sp_tmp;

			// 時間ステップごとの運動
			for (int i = 0; i <= jloop; i++) {
				OV1.process_repump(rb87);
				OV1.process_dipole(rb87);
				OV1.process_diss(rb87);
				OV1.step_motion(rb87);

				if (rb87->r.z < -0.26) {
					count_guide++;
					sum_sp += OV1.count_sp;
					double angVf = ( -(rb87->v.vx) * sin(rb87->phi) + (rb87->v.vy) * cos(rb87->phi))* rb87->radius + rb87->l_rot / (mass);
					ofs << angVf << endl;
					if (angVf > 0) count_vphi++;
					break;
				}

				if (rb87->radius > w0 / sqrt(2.0)) {
					break;
				}
			}
		}
			
		printf("avarage spontaneous emission %d/%d times, ", sum_sp, count_guide);
		printf("guide efficiency %lf, ", (double)count_guide/(double)SAMPLE *100);
		printf("unity of azimuthal direction %lf \n", (2.0*(double)count_vphi-(double)count_guide) / (double)count_guide *100);			// R-L/R+L [%]

	 	ofs << (double)count_guide/(double)SAMPLE << ", " <<  (double)sum_sp/(double)count_guide << ", " << (2.0*(double)count_vphi-(double)count_guide) / (double)count_guide << endl;
	    return 0;

}


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
	} while (p_uni > mass * giji_ransu * giji_ransu / (2.0 * k_b * temp) *exp( 1.0 -mass*giji_ransu * giji_ransu / (2.0 *k_b*temp)));
	*vx = giji_ransu;
}

void rm_position(double* x, double* y, double* z){
	std::mt19937 rand_src(std::random_device{}());
    std::uniform_real_distribution<double> dist(0.0, 1.0);

	for (int j = 0; j < SAMPLE; j++) {
		double r0;
		gauss(&r0);

		double phi = 2.0 * M_PI * dist(rand_src);
    	double psi = M_PI * dist(rand_src);

		x[j] = r0* cos(phi) * sin(psi);
		y[j] = r0* sin(phi) * sin(psi);
		z[j] = r0* cos(psi);
	}
}

void rm_velocity(double* vx, double* vy, double* vz){
	std::mt19937 rand_src(std::random_device{}());
    std::uniform_real_distribution<double> dist(0.0, 1.0);

	for (int j = 0; j < SAMPLE; j++) {
		double v0;
		max_boltz(&v0);

		double phi = 2.0 * M_PI * dist(rand_src);
    	double psi = M_PI * dist(rand_src);

		vx[j] = v0* cos(phi) * sin(psi);
		vy[j] = v0* sin(phi) * sin(psi);
		vz[j] = - v0 * abs(cos(psi));
	}
}