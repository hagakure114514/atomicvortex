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
	// spontaneous emission mode  
    printf("select spontaneous emission mode?\n 0:no-sp, 1:z-axis-sp, 2:dipole-radiation, 3:all direction \n");
    char c = std::cin.get();
	

	char fname[30];
	sprintf_s(fname, "traj_pos50_vel0.csv");
	ofstream ofs(fname);        // ファイルパスを指定する

	// オブジェクトのコンストラクタ
	position r0={ 0.0, 50.0e-6, 0.0};
	velocity v0={ 1.0e-2, 0.0, 0.0};
    atom Rb87(r0, v0, state::d2);		// 原子オブジェクト
    atom* rb87 = &Rb87;
    DressedAtom OV1;			// dressed-atom状態オブジェクト
	OV1.flag_sp = (int)c;


    // 配列の定義
    double x[jloop + 1] = {}, y[jloop + 1] = {}, z[jloop+1] = {};

    // 時間ステップごとの運動
	int lim=0;
    for (int i = 0; i <= jloop; i++) {
		lim = i;
        OV1.process_dipole(rb87);
        OV1.process_diss(rb87);
        OV1.step_motion(rb87);

        x[i] = rb87->r.x;  y[i] = rb87->r.y;  z[i] = rb87->r.z;

		if (z[i] < -0.30) {
			lim = i - 1;
			printf("z potision under -30 cm with %d processes (%e s)\n", lim, 5.0e-5 *i);
			break;
		}

		if (rb87->radius > w0/sqrt(2.0)) {
			printf("conf radius out at z=%e cm\n", z[i] * 1e2);
			break;
		}

		//printf("%e, %e, %e\n", x[i], y[i], z[i]);
		ofs << x[i] << ", " << y[i] << ", " << z[i] << endl;
    }

	printf("spontaneous_emission_%d_times\n", OV1.count_sp);

    return 0;
}