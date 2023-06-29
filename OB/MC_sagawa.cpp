//
//  main.cpp
//  funnel_MC
//
//  Created by KEN on 2014/01/20.
//  Copyright (c) 2014年 KEN. All rights reserved.
//
//	コンパイルは"g++ -O2 Source.cpp"

// 2014-12-17 ファンデルワールス力(vdw_force)を抜いた式で計算


#include <iostream>
#include <Eigen/Dense>
#include <math.h>
#include <stdlib.h>
#include <time.h>	
#include <fstream>
#include <sstream>
#include <iomanip>
using namespace Eigen;
using namespace std;

static const int  SAMPLE_ALL = 1000;
static const int  MOT_ATOMS = 10 * 1000 * 1000;
static const double ANGLE_SLOPE = 44.6;
static const int NUM_SLOPE = 360;     //斜面数360：円錐, 4：四角錐, 3：三角錐

static const double ANGLE_ONEFACE = 360 / NUM_SLOPE;
static const double MASS = 87.0e-3 / 6.022e23;
static const double K_BOLTZ = 1.380648E-23;  // ボルツマン定数
static const double TEMP_MOT = 10E-6;
static const double SIZE_MOT = 1.0E-3 / 2;   // MOT半径
static const double HEIGHT_MOT = 3.0E-3;
static const double D_HFS = 2.0 * M_PI * 6834.0E6;
static const double DELTA1 = 2.0 * M_PI * 1000.0E6;     //離調1GHz
static const double DELTA2 = DELTA1 + D_HFS;
static const double N_FUN = 1.45;
static const double HBAR = 1.0545E-34;
static const double LAMDA = 780E-9;
static const double DECAY = LAMDA / (2 * M_PI * sqrt(pow(N_FUN, 2) * pow(sin(ANGLE_SLOPE / 180 * M_PI), 2) - 1));
//static const double DECAY = 406E-9;
static const double GRAVITY = 9.80665;
static const int MAX_BOUND = 100000;
static const int MAX_TIME = 1000;
static const int t_plot_hist = 1000;


static const double EVANESCENT_AREA = LAMDA;
static const double JUDGE_AREA = 1.0E-9;
static const double GAMMA0 = 6.0E6;
static const double Q1 = 0.75;
static const double I_SAT1 = Q1 * 1.58E1;
static const double I_SAT2 = (1 - Q1) * 1.58E1;
static const double GMM1 = 2 * M_PI * Q1 * GAMMA0;
static const double GMM2 = 2 * M_PI * (1 - Q1) * GAMMA0;
static const double KJ1 = 2.0 * M_PI / 795.0E-9;
static const double KJ2 = 2.0 * M_PI / 780.2E-9;

static const double PK = 1 / MASS * HBAR * (2 * M_PI / LAMDA);  //# 自然放出光子の反跳
static const double EPK = PK * N_FUN * sin(ANGLE_SLOPE);  //# エバネッセント光から受ける反跳

double f1_x(double v, double r);
double f1_y(double v, double r);
double f1_z(double v, double r);
double f2(double v, double r);
double vdw_force(double r);
double intensity_dash(double r, double r_xy, double CENT, double WAIST, double POWER);
double fn(double r, double r_xy, int state, double CENT, double WAIST, double POWER, int area, int pol_num);

double gn(double vr);

double radians(double(degree))
{
	return degree / 180 * M_PI;
}
double degrees(double(radian))
{
	return radian / M_PI * 180;
}

////////////////////////////
//  どの面上にいるかチェック
////////////////////////////
int area_check(Vector3d vec_position)
{
	int area_num = 0;
	double x = vec_position(0);
	double y = vec_position(1);
	double angle = degrees(atan(y / x));
	if (x < 0)
	{
		angle += 180;
	}
	else if (x > 0 && y < 0)
	{
		angle += 360;
	}
	if (angle > 360 - ANGLE_ONEFACE / 2)
	{
		angle = angle - 360;
	}

	for (int i = 0; i < NUM_SLOPE; i++)
	{
		double phy_min = 0.0 - ANGLE_ONEFACE / 2 + i * ANGLE_ONEFACE;
		double phy_max = 0.0 + ANGLE_ONEFACE / 2 + i * ANGLE_ONEFACE;

		if (phy_min <= angle && angle < phy_max)
		{
			area_num = i;
			break;
		}
	}
	return area_num;
}

////////////////////////////
//  面と点との距離
////////////////////////////

//todo nvecを何度も呼び出してたら重い？
Vector3d check_n_vec(MatrixXd n_vec_matrix, int area_num)
{
	Vector3d n_vec;
	n_vec = n_vec_matrix.row(area_num);
	return n_vec;
}
double measure_height(Vector3d vector_p, Vector3d vector_n)
{
	//Vector3d OP = vector_p;
	//Vector3d n_vec = vector_n;
	return vector_p.dot(vector_n);
}

////////////////////////////
//  バウンド
////////////////////////////
Vector3d bound(Vector3d vector_v, Vector3d vector_n)
{
	//Vector3d n_vec = vector_n;
	Vector3d velocity_bound = vector_v - 2.0 * vector_v.dot(vector_n) * vector_n;
	return velocity_bound;
}

Vector3d check_h_vec(MatrixXd h_vec_matrix, int area_num)
{
	Vector3d h_vec;
	h_vec = h_vec_matrix.row(area_num);
	return h_vec;
}

Vector3d bound_sisyphus(Vector3d vector_v, Vector3d vector_n, Vector3d vector_h)
{
	//Vector3d n_vec = vector_n;
	//Vector3d h_vec = vector_h;
	double theta = rand() / RAND_MAX * M_PI;
	double phi = rand() / RAND_MAX * 2 * M_PI;
	Vector3d random_vec(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
	Vector3d bound_vec;

	bound_vec = vector_v - (1.0 + sqrt(DELTA1 / DELTA2)) * vector_v.dot(vector_n) * vector_n + EPK * vector_h + PK * random_vec;
	return bound_vec;
}

Vector3d bound_spontaneous(Vector3d vector_v, Vector3d vector_n, Vector3d vector_h)
{
	//vector3d n_vec = vector_n;
	//vector3d h_vec = vector_h;
	double theta = rand() / RAND_MAX * M_PI;
	double phi = rand() / RAND_MAX * 2 * M_PI;
	Vector3d bound_vec;
	Vector3d random_vec(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
	bound_vec = vector_v - (1.0 + 1.0) * vector_v.dot(vector_n) * vector_n + EPK * vector_h + PK * random_vec;
	return bound_vec;
}

Vector3d pumping(Vector3d vector_v)
{
	double theta = rand() / RAND_MAX * M_PI;
	double phi = rand() / RAND_MAX * 2 * M_PI;
	Vector3d pump_vec;
	Vector3d random_vec(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
	return vector_v + PK * random_vec;
}

////////////////////////////
//  斜面への入射速度成分
////////////////////////////
double projection(Vector3d vector_v, Vector3d vector_n)
{
	return  vector_v.dot(vector_n);
}


int main(int argc, char* argv[])
{
	if (argc != 4)
	{
		std::cout << "please input valuables" << endl;
		std::cout << "val1 = exit, val2 = power, val3 = polarization(0=c-pol, 1=p-pol, 2=SPRなし)" << endl;
		exit(1);
	}
	int input_exit;
	input_exit = atoi(argv[1]);
	//std::cout << input_exit << endl;
	int input_power;
	input_power = atoi(argv[2]);
	//std::cout << input_power << endl;
	int pol_num;
	pol_num = atoi(argv[3]);
	//std::cout << pol_num << endl;
	static const int POWER_filename = input_power;
	static const int EXIT_filename = input_exit;

	string pol = "circle";
	int ENHANCE = 7;
	if (pol_num == 1)
	{
		ENHANCE = 14;
		pol = "linear";
	}
	else if (pol_num == 2)
	{
		ENHANCE = 1;
		pol_num = 0;
		pol = "no_pol";
	}

	//std::cout << pol << "-pol, " << ENHANCE << endl;

	static const double POWER = ENHANCE * POWER_filename * 1.0E-3;
	static const double EXIT = EXIT_filename * 1.0E-6;
	static const double CENT = EXIT * 0.5;
	static const double WAIST = 3.0E-3;
	static const double HEIGHT_EXIT = EXIT / 2 * tan(ANGLE_SLOPE / 180 * M_PI);

	////////////////////////////
	//  計算時間計測
	////////////////////////////
	clock_t start, end;
	start = clock();


	////////////////////////////
	//  ファネルの面の頂点を決める
	////////////////////////////
	double phy_zero;
	double angle_phy_zero;
	double height_vertex;
	double phy, angle_phy;
	double xx, yy, zz;
	Vector3d p0, p1;
	MatrixXd p2, p3;
	p0 = Vector3d::Zero(3);
	p1 = Vector3d::Zero(3);
	p2 = MatrixXd::Zero(NUM_SLOPE, 3);
	p3 = MatrixXd::Zero(NUM_SLOPE, 3);

	phy_zero = -ANGLE_ONEFACE / 2;
	angle_phy_zero = radians(phy_zero);
	height_vertex = cos(angle_phy_zero) * tan(radians(ANGLE_SLOPE));

	for (int i = 0; i < NUM_SLOPE; i++)
	{
		phy = ANGLE_ONEFACE * i - ANGLE_ONEFACE / 2;
		angle_phy = radians(phy);
		xx = cos(angle_phy);
		yy = sin(angle_phy);
		zz = height_vertex;
		p2.row(i) << xx / zz, yy / zz, zz / zz;

		phy = ANGLE_ONEFACE * (i + 1) - ANGLE_ONEFACE / 2;
		angle_phy = radians(phy);
		xx = cos(angle_phy);
		yy = sin(angle_phy);
		zz = height_vertex;
		p3.row(i) << xx / zz, yy / zz, zz / zz;
	}
	//ofstream ofs;
	//ofs.open("vertex.dat");
	//ofs << p2 << endl;
	//ofs.close();

	//cout << p3 << endl;

	////////////////////////////
	//  面の単位法線ベクトル
	////////////////////////////
	MatrixXd n_vec_matrix;
	Vector3d n_vec, AB, AC;
	n_vec_matrix = MatrixXd::Zero(NUM_SLOPE, 3);
	n_vec = Vector3d::Zero(3);
	//cout << p2 << endl << endl;
	//cout << p3 << endl << endl;
	for (int i = 0; i < NUM_SLOPE; i++)
	{
		AB = p2.row(i);
		AC = p3.row(i);
		//cout << AB << " " << AC << endl;

		n_vec = AB.cross(AC);
		double vec_len = sqrt(n_vec.dot(n_vec));
		n_vec /= vec_len;
		n_vec_matrix.row(i) = n_vec;
	}
	//cout << n_vec_matrix << endl;

	////////////////////////////
	//  面平行ベクトル
	////////////////////////////
	MatrixXd h_vec_matrix;
	Vector3d h_vec;
	h_vec_matrix = MatrixXd::Zero(NUM_SLOPE, 3);
	h_vec = Vector3d::Zero(3);
	for (int i = 0; i < NUM_SLOPE; i++)
	{
		AB = p2.row(i);
		AC = p3.row(i);
		h_vec = AB + AC;
		double vec_len = sqrt(h_vec.dot(h_vec));
		h_vec /= vec_len;
		h_vec_matrix.row(i) = h_vec;
	}
	//cout << h_vec_matrix << endl;

	////////////////////////////
	//  速度の初期値
	////////////////////////////
	//乱数で正規分布を与える(Box - Muller法)
	srand((unsigned)time(NULL));
	double v_xyz[3 * SAMPLE_ALL] = {};
	for (int i = 0; i < 3 * SAMPLE_ALL; i++)
	{
		double rnd1 = rand() / (double)RAND_MAX;
		double rnd2 = rand() / (double)RAND_MAX;
		v_xyz[i] = sqrt(-2 * log(rnd1)) * cos(2 * M_PI * rnd2);
	}

	//各成分の和
	double v_sum[3] = {};
	for (int i = 0; i < SAMPLE_ALL * 3; i = i + 3)
	{
		for (int j = 0; j < 3; j++)
		{
			v_sum[j] += v_xyz[i + j];
		}
	}

	//成分の和をゼロにする(中心からのズレを修正する)
	double v_rev[3 * SAMPLE_ALL] = {};
	double v_average[3] = {};
	for (int i = 0; i < 3; i++)
	{
		v_average[i] = v_sum[i] / SAMPLE_ALL;
	}
	for (int i = 0; i < SAMPLE_ALL * 3; i = i + 3)
	{
		for (int j = 0; j < 3; j++)
		{
			v_rev[i + j] = v_xyz[i + j] - v_average[j];
		}
	}

	//温度修正
	double temp_atom = 0.0;
	double v_xyz_for_check = 0.0;
	for (int i = 0; i < SAMPLE_ALL; i++)
	{
		v_xyz_for_check += pow(v_rev[i * 3], 2) + pow(v_rev[i * 3 + 1], 2) + pow(v_rev[i * 3 + 2], 2);
	}
	//temp_atom = v_xyz_for_check / (3 * (SAMPLE_ALL - 1)) * MASS / K_BOLTZ;
	temp_atom = v_xyz_for_check / (3 * (SAMPLE_ALL)) * MASS / K_BOLTZ;
	
	//設定温度になるように速度に係数をかける
	double coefficient_temp = sqrt(TEMP_MOT / temp_atom);
	double v_rev_rev[3 * SAMPLE_ALL] = {};
	for (int i = 0; i < SAMPLE_ALL * 3; i = i + 3)
	{
		for (int j = 0; j < 3; j++)
		{
			v_rev_rev[i + j] = v_rev[i + j] * coefficient_temp;
		}
	}

	//速度を３次元配列に格納
	MatrixXd v_initialized;
	v_initialized = MatrixXd::Zero(SAMPLE_ALL, 3);
	for (int row = 0; row < SAMPLE_ALL; row++)
	{
		for (int column = 0; column < 3; column++)
		{
			v_initialized(row, column) = v_rev_rev[row * 3 + column];
		}
	}

	/////////////////////////////////////////
	//    位置の初期値
	/////////////////////////////////////////
	//乱数で正規分布を与える(Box - Muller法)
	srand((unsigned)time(NULL));
	double p_xyz[3 * SAMPLE_ALL] = {};
	for (int i = 0; i < 3 * SAMPLE_ALL; i++)
	{
		double rnd1 = rand() / (double)RAND_MAX;
		double rnd2 = rand() / (double)RAND_MAX;
		p_xyz[i] = SIZE_MOT * sqrt(-2 * log(rnd1)) * cos(2 * M_PI * rnd2);
	}

	//各成分の和
	double p_sum[3] = {};
	for (int i = 0; i < SAMPLE_ALL * 3; i = i + 3)
	{
		for (int j = 0; j < 3; j++)
		{
			p_sum[j] += p_xyz[i + j];
		}
	}

	//成分の和をゼロにする(中心からのズレを修正する)
	double p_rev[3 * SAMPLE_ALL] = {};
	double p_average[3] = {};
	for (int i = 0; i < 3; i++)
	{
		p_average[i] = p_sum[i] / SAMPLE_ALL;
	}
	for (int i = 0; i < SAMPLE_ALL * 3; i = i + 3)
	{
		for (int j = 0; j < 3; j++)
		{
			p_rev[i + j] = p_xyz[i + j] - p_average[j];
		}
	}

	//MOTの高さを加える
	for (int i = 0; i < 3 * SAMPLE_ALL; i = i + 3)
	{
		p_rev[i + 2] += HEIGHT_MOT;
	}

	//位置を３次元配列に格納
	MatrixXd p_initialized;
	p_initialized = MatrixXd::Zero(SAMPLE_ALL, 3);
	for (int row = 0; row < SAMPLE_ALL; row++)
	{
		for (int column = 0; column < 3; column++)
		{
			p_initialized(row, column) = p_rev[row * 3 + column];
		}
	}
	//ofstream ofs_init_p;
	//ofs_init_p.open("init_p.dat");
	//ofs_init_p << p_initialized << endl;
	//ofs_init_p.close();
	//ofstream ofs_init_v;
	//ofs_init_v.open("init_v.dat");
	//ofs_init_v << v_initialized << endl;
	//ofs_init_v.close();


	/////////////////////////////////
	//  file open
	/////////////////////////////////

	//ofstream ofs_output;
	//ofs_output.open("output.dat");

	//ofstream ofs_trace;
	//ofs_trace.open("trace.dat");



	ofstream ofs_time;
	ostringstream osst;
	osst << "time_" << pol << "_" << EXIT_filename << "um_" << POWER_filename << "mw_" << NUM_SLOPE << "face_" << ANGLE_SLOPE << ".dat";
	ofs_time.open(osst.str().c_str());
	//cout << osst.str() << endl;


	ofstream loss_log;
	ostringstream losslog;
	losslog << "losslog_" << pol << "_" << EXIT_filename << "um_" << POWER_filename << "mw_" << NUM_SLOPE << "face_" << ANGLE_SLOPE << ".csv";
	loss_log.open(losslog.str().c_str());




	// for (int i = 1; i <= t_plot_hist; i++)
	// {
	// ostringstream oss1;
	// ostringstream oss2;
	// oss1 << "time_log_" << i;
	// ofstream file;
	// oss2 << "timelog_" << pol << "_" << EXIT_filename << "um_" << POWER_filename << "mw_" << NUM_SLOPE << "face_" << ANGLE_SLOPE << "_" << std::setw( 3 ) << std::setfill( '0' ) << i << ".csv";
	// file.open(oss2.str(), ios::trunc);
	// file.close();
	// //file.open(oss2.str(), ios::app);
	// }

/*
	ofstream time_log;
	ostringstream timelog;
	timelog << "timelog_" << pol << "_" << EXIT_filename << "um_" << POWER_filename << "mw_" << NUM_SLOPE << "face_" << ANGLE_SLOPE << "_01.csv";
	time_log.open(timelog.str(), ios::trunc);
	time_log.close();
	time_log.open(timelog.str(), ios::app);
*/
	/////////////////////////////////
	//ここまで準備
	/////////////////////////////////
	double x, y, z;
	double vx, vy, vz;
	double r, r_xy;
	double vr, pre_vr, v_in;
	double p;
	double k1x, k1y, k1z;
	double l1x, l1y, l1z;
	double k2x, k2y, k2z;
	double l2x, l2y, l2z;
	double k3x, k3y, k3z;
	double l3x, l3y, l3z;
	double k4x, k4y, k4z;
	double l4x, l4y, l4z;
	double k1, k2, k3, k4;
	double l1, l2, l3, l4;
	double t, dt;
	double pre_t;
	double output_time;
	int output_atom;
	int timeout_count;
	double height;
	int area;
	int quitloop;
	int count_bound;
	int state;
	int flag_plot;
	int flag_plot_pre;
	double temperature_output_atom;
	double rnd1, rnd2;
	double t_sum;
	double temperature_sum;
	double flux_intensity;
	double average_t;
	double average_temp;
	double gain;
	double t_plot;
	stringstream plot_file;
	static const int n_hist = MAX_TIME;
	int hist[n_hist] = { 0 };
	int hist_total[n_hist] = { 0 };
	Vector3d output_v;
	Vector3d init_p, init_v;
	Vector3d after_p, after_v;
	Vector3d vec_n, vec_h;
	Vector3d pre_p, pre_v;
	Vector3d incident_p, incident_v;

	/////////////////////////////////
	//  初期状態
	/////////////////////////////////
	output_atom = 0;
	t_sum = 0.0;
	temperature_sum = 0.0;
	timeout_count = 0;
	t_plot = 0.0;

	for (int sample_num = 0; sample_num < SAMPLE_ALL; sample_num++)
	{
		if ((sample_num * 100) % SAMPLE_ALL == 0)
		{
			std::cout << (sample_num * 100) / SAMPLE_ALL << "%  " << pol << "_" << EXIT_filename << "um_" << POWER_filename << "mw_" << NUM_SLOPE << "face" << endl;
		}


		init_p = p_initialized.row(sample_num);
		init_v = v_initialized.row(sample_num);

		after_p = init_p;
		after_v = init_v;

		t = 0.0;
		area = 0;
		area = area_check(after_p);

		vec_n = check_n_vec(n_vec_matrix, area);
		height = measure_height(after_p, vec_n);

		flag_plot = 0;
		flag_plot_pre = 0;

		/////////////////////////////////
		//  スタート
		/////////////////////////////////
		quitloop = 0;
		count_bound = 0;
		state = 1;


		/////////////////////////////////
		//  原子解放追跡 "free space loop" and "evanescent loop"
		/////////////////////////////////

		while (1)
		{
			dt = 1.0E-4;

			/////////////////////////////////
			//  free space loop
			/////////////////////////////////
			while (1)// free
			{
				pre_p = after_p;
				pre_v = after_v;
				pre_t = t;
				t += dt;

				if (t >= MAX_TIME)
				{
					std::cout << "Timeout!" << endl;
					timeout_count++;
					if (timeout_count > 1000)
					{
						std::cout << "timeout>10 exit_program";
						exit(1);
					}
					quitloop = 1;
					break;
				}

				// for (int i = 1; i <= t_plot_hist; i++)
				// {
				// 	t_plot = i * 0.001;
				// 	if (t >= t_plot && t < t_plot + 0.001)
				// 	{
				// 		flag_plot_pre = flag_plot;
				// 		flag_plot = i;
				// 		if (flag_plot - flag_plot_pre == 1)
				// 		{
				// 			//plot_file << "file";
				// 			//std::string aaa = plot_file.str();
				// 			//std::cout << t << endl;
				// 			ofstream time_log;
				// 			ostringstream ossss;
				// 			ossss << "timelog_" << pol << "_" << EXIT_filename << "um_" << POWER_filename << "mw_" << NUM_SLOPE << "face_" << ANGLE_SLOPE << "_" << std::setw( 3 ) << std::setfill( '0' ) << i << ".csv";
				// 			time_log.open(ossss.str(), ios::app);
				// 			time_log << x * 1.0E3 << ", " << y * 1.0E3 << ", " << z * 1.0E3 << endl;
				// 			time_log.close();

				// 		}
				// 	}
				// }
				/*
				if (t >= 0.1 && t < 0.2)
				{
					flag_plot_pre = flag_plot;
					flag_plot = 1;
					if (flag_plot - flag_plot_pre == 1)
					{
						std::cout << t << endl;
						time_log << x * 1.0E3 << ", " << y * 1.0E3 << ", " << z * 1.0E3 << endl;
					}
				}
				else if (t >= 0.2 && t < 0.3)
				{
					flag_plot_pre = flag_plot;
					flag_plot = 2;
					if (flag_plot - flag_plot_pre == 1)
					{
						std::cout << t << endl;
						time_log2 << x * 1.0E3 << ", " << y * 1.0E3 << ", " << z * 1.0E3 << endl;
					}
				}
				else if (t >= 0.3 && t < 0.4)
				{
					flag_plot_pre = flag_plot;
					flag_plot = 3;
					if (flag_plot - flag_plot_pre == 1)
					{
						std::cout << t << endl;
						time_log3 << x * 1.0E3 << ", " << y * 1.0E3 << ", " << z * 1.0E3 << endl;
					}
				}

				else if (t >= 0.4 && t < 0.5)
				{
					flag_plot_pre = flag_plot;
					flag_plot = 4;
					if (flag_plot - flag_plot_pre == 1)
					{
						std::cout << t << endl;
						time_log4 << x * 1.0E3 << ", " << y * 1.0E3 << ", " << z * 1.0E3 << endl;
					}
				}

				else if (t >= 0.5 && t < 0.6)
				{
					flag_plot_pre = flag_plot;
					flag_plot = 5;
					if (flag_plot - flag_plot_pre == 1)
					{
						std::cout << t << endl;
						time_log5 << x * 1.0E3 << ", " << y * 1.0E3 << ", " << z * 1.0E3 << endl;
					}
				}

				else if (t >= 0.6 && t < 0.7)
				{
					flag_plot_pre = flag_plot;
					flag_plot = 6;
					if (flag_plot - flag_plot_pre == 1)
					{
						std::cout << t << endl;
						time_log6 << x * 1.0E3 << ", " << y * 1.0E3 << ", " << z * 1.0E3 << endl;
					}
				}

				else if (t >= 0.7 && t < 0.8)
				{
					flag_plot_pre = flag_plot;
					flag_plot = 7;
					if (flag_plot - flag_plot_pre == 1)
					{
						std::cout << t << endl;
						time_log7 << x * 1.0E3 << ", " << y * 1.0E3 << ", " << z * 1.0E3 << endl;
					}
				}

				else if (t >= 0.8 && t < 0.9)
				{
					flag_plot_pre = flag_plot;
					flag_plot = 8;
					if (flag_plot - flag_plot_pre == 1)
					{
						std::cout << t << endl;
						time_log8 << x * 1.0E3 << ", " << y * 1.0E3 << ", " << z * 1.0E3 << endl;
					}
				}

				else if (t >= 0.9 && t < 1.0)
				{
					flag_plot_pre = flag_plot;
					flag_plot = 9;
					if (flag_plot - flag_plot_pre == 1)
					{
						std::cout << t << endl;
						time_log9 << x * 1.0E3 << ", " << y * 1.0E3 << ", " << z * 1.0E3 << endl;
					}
				}

				else if (t > 1.0)
				{
					flag_plot_pre = flag_plot;
					flag_plot = 10;
					if (flag_plot - flag_plot_pre == 1)
					{
						std::cout << t << endl;
						time_log10 << x * 1.0E3 << ", " << y * 1.0E3 << ", " << z * 1.0E3 << endl;
					}
				}
*/

				x = pre_p(0);
				y = pre_p(1);
				z = pre_p(2);
				vx = pre_v(0);
				vy = pre_v(1);
				vz = pre_v(2);

				k1x = dt * f1_x(vx, x);
				k1y = dt * f1_y(vy, y);
				k1z = dt * f1_z(vz, z);
				l1x = dt * f2(vx, x);
				l1y = dt * f2(vy, y);
				l1z = dt * f2(vz, z);

				k2x = dt * f1_x(vx + k1x / 2, x + l1x / 2);
				k2y = dt * f1_y(vy + k1y / 2, y + l1y / 2);
				k2z = dt * f1_z(vz + k1z / 2, z + l1z / 2);
				l2x = dt * f2(vx + k1x / 2, x + l1x / 2);
				l2y = dt * f2(vy + k1y / 2, y + l1y / 2);
				l2z = dt * f2(vz + k1z / 2, z + l1z / 2);

				k3x = dt * f1_x(vx + k2x / 2, x + l2x / 2);
				k3y = dt * f1_y(vy + k2y / 2, y + l2y / 2);
				k3z = dt * f1_z(vz + k2z / 2, z + l2z / 2);
				l3x = dt * f2(vx + k2x / 2, x + l2x / 2);
				l3y = dt * f2(vy + k2y / 2, y + l2y / 2);
				l3z = dt * f2(vz + k2z / 2, z + l2z / 2);

				k4x = dt * f1_x(vx + k3x, x + l3x);
				k4y = dt * f1_y(vy + k3y, y + l3y);
				k4z = dt * f1_z(vz + k3z, z + l3z);
				l4x = dt * f2(vx + k3x, x + l3x);
				l4y = dt * f2(vy + k3y, y + l3y);
				l4z = dt * f2(vz + k3z, z + l3z);

				after_v(0) = vx + (k1x + 2 * k2x + 2 * k3x + k4x) / 6;
				after_v(1) = vy + (k1y + 2 * k2y + 2 * k3y + k4y) / 6;
				after_v(2) = vz + (k1z + 2 * k2z + 2 * k3z + k4z) / 6;
				after_p(0) = x + (l1x + 2 * l2x + 2 * l3x + l4x) / 6;
				after_p(1) = y + (l1y + 2 * l2y + 2 * l3y + l4y) / 6;
				after_p(2) = z + (l1z + 2 * l2z + 2 * l3z + l4z) / 6;

				if (state == 2)
				{
					after_v = pumping(after_v);
					state = 1;
				}

				area = area_check(after_p);
				vec_n = check_n_vec(n_vec_matrix, area);
				height = measure_height(after_p, vec_n);
				//cout << after_p(2) << endl;
				if (pow(after_p(0), 2) + pow(after_p(1), 2) < pow(EXIT / 2, 2) && after_p(2) < HEIGHT_EXIT)
				{
					output_atom++;
					t_sum += t;
					std::cout << "      " << output_atom << " atoms output! ";
					output_time = t;
					output_v = after_v;
					temperature_output_atom = (pow(output_v(0), 2) + pow(output_v(1), 2) + pow(output_v(2), 2)) / 3 * MASS / K_BOLTZ;
					temperature_sum += temperature_output_atom;

					ofs_time << output_time << endl;
					//cout << "output time = " << output_time << endl;
					//cout << "output temp = " << temperature_output_atom << endl;

					//ofs_output << "output time = " << output_time << endl;
					//ofs_output << "output temp = " << temperature_output_atom << endl;
					//time_log << after_p(0) * 1.0E3 << ", " << after_p(1) * 1.0E3 << ", " << after_p(2) * 1.0E3 << endl;

					quitloop = 1;

					for (int i = 0; i < n_hist; i++)
					{
						if (output_time < i + 1)
						{
							hist[i]++;
							break;
						}
					}

					break;
				}

				if (height < EVANESCENT_AREA)
				{
					after_p = pre_p;
					after_v = pre_v;
					t = pre_t;
					dt /= 10;
					continue;
				}
				else if (height < EVANESCENT_AREA + JUDGE_AREA)
				{
					incident_p = after_p;
					incident_v = after_v;
					break;
				}
				else
				{
					x = after_p(0);
					y = after_p(1);
					z = after_p(2);
					//ofs_trace << t << ", " << x << ", " << y << ", " << z << endl;
					continue;
				}

			}//while free
			if (quitloop == 1)
			{
				std::cout << "Finish at free area" << endl;
				break;
			}

			/////////////////////////////////
			//  evanescent loop
			/////////////////////////////////
			r = height;
			vr = projection(incident_v, vec_n);
			v_in = vr;
			r_xy = sqrt(pow(incident_p(0), 2) + pow(incident_p(1), 2));

			dt = 1.0E-9;

			while (1)
			{
				pre_vr = vr;
				t += dt;

				k1 = dt * fn(r, r_xy, state, CENT, WAIST, POWER, area, pol_num);
				l1 = dt * gn(vr);

				k2 = dt * fn(r, r_xy, state, CENT, WAIST, POWER, area, pol_num);
				l2 = dt * gn(vr + k1 / 2.0);

				k3 = dt * fn(r, r_xy, state, CENT, WAIST, POWER, area, pol_num);
				l3 = dt * gn(vr + k2 / 2.0);

				k4 = dt * fn(r, r_xy, state, CENT, WAIST, POWER, area, pol_num);
				l4 = dt * gn(vr + k3);

				vr += (k1 + 2 * k2 + 2 * k3 + k4) / 6;
				r += (l1 + 2 * l2 + 2 * l3 + l4) / 6;

				if (r < 0.0)
				{
					std::cout << "      Loss" << endl;
					loss_log << after_p(0) * 1.0E3 << ", " << after_p(1) * 1.0E3 << ", " << after_p(2) * 1.0E3 << endl;
					quitloop = 1;
					break;
				}
				else if (vr * pre_vr < 0)
				{
					count_bound++;
					if (count_bound > MAX_BOUND)
					{
						std::cout << "      bound max" << endl;
						quitloop = 1;
						break;
					}
					break;
				}
				else if (r > EVANESCENT_AREA + JUDGE_AREA)
				{
					break;
				}
				else
				{
					x = after_p(0);
					y = after_p(1);
					z = after_p(2);
					continue;
				}

			}//while evanescent
			if (quitloop == 1)
			{
				break;
			}

			if (state == 1)
			{
				rnd1 = rand() / double(RAND_MAX);
				rnd2 = rand() / double(RAND_MAX);

				p = 1.0 - exp(v_in * MASS * DECAY * GAMMA0 * M_PI / (HBAR * DELTA1));
				//p = 0.0;  //without spontaneous emittion
				//cout << rnd1 << ", " << p << endl;
				if (rnd1 < p)
				{
					vec_h = check_h_vec(h_vec_matrix, area);
					if (rnd2 < 0.25)
					{
						after_v = bound_sisyphus(incident_v, vec_n, vec_h);
						state = 2;
					}
					else
					{
						after_v = bound_spontaneous(incident_v, vec_n, vec_h);
						state = 1;
					}
				}
				else
				{
					after_v = bound(incident_v, vec_n);
					state = 1;
				}
			}
			else
			{
				after_v = bound(incident_v, vec_n);
				state = 2;
			}
			after_p = incident_p;


		}//while free + evanescent
	}

	average_t = t_sum / output_atom;
	average_temp = temperature_sum / output_atom * 1.0E6;
	flux_intensity = output_atom / double(SAMPLE_ALL - timeout_count) * double(MOT_ATOMS) / (M_PI * pow(EXIT / 2, 2)) / 1.0 * 1.0E-4;
	gain = output_atom / double(SAMPLE_ALL - timeout_count) * 100;


	for (int i = 0; i < n_hist; i++)
	{
		for (int j = 0; j <= i; j++)
		{
			hist_total[i] += hist[j];
		}

	}

	double total_flux = hist_total[n_hist - 1] / double(SAMPLE_ALL - timeout_count) * double(MOT_ATOMS) / (M_PI * pow(EXIT / 2, 2)) / 1.0 * 1.0E-4;


	std::cout << "average t " << average_t << endl;
	std::cout << "average temp (uk)" << average_temp << endl;
	std::cout << "flux intensity " << flux_intensity << endl;
	std::cout << "fluxtotal = " << total_flux << endl;
	std::cout << "gain " << gain << endl;


	ofstream ofs_flux;
	ostringstream oss;
	oss << "flux_" << pol << "_" << EXIT_filename << "um_" << POWER_filename << "mw_" << NUM_SLOPE << "face_" << ANGLE_SLOPE << ".dat";
	ofs_flux.open(oss.str().c_str());
	//std::cout << oss.str() << endl;

	ofstream power_vs_flux;
	ostringstream pvsf;
	pvsf << "plot_" << pol << "_" << NUM_SLOPE << "face_" << ANGLE_SLOPE << ".csv";
	//pvsf << "plot_" << pol << "_" << EXIT_filename << "um_" << NUM_SLOPE << "face_" << ANGLE_SLOPE << ".csv";
	power_vs_flux.open(pvsf.str().c_str(), ios::app);

	ofstream hist_time;
	ostringstream hist_t;
	hist_t << "hist_" << pol << "_" << EXIT_filename << "um_" << POWER_filename << "mw_" << NUM_SLOPE << "face_" << ANGLE_SLOPE << ".csv";
	hist_time.open(hist_t.str().c_str());


	ofs_flux << "#flux_intensity, " << "gain, " << "temp, " << "t, " << "beam_waist, " << "delta "<< endl;
	ofs_flux << flux_intensity << ", " << gain << ", " << average_temp << ", " << average_t << ", "<<  WAIST << ", "<< DELTA1 << endl;
	std::cout << oss.str() << endl;

	power_vs_flux << POWER_filename << ", " << double(POWER_filename) / 1000 / 3.14159 / (WAIST*WAIST*1E4) << ", " << EXIT_filename << ", " << total_flux << ", " << gain << ", " << average_temp << ", " << average_t << endl;


	for (int i = 0; i < 100; i++)
	{
		hist_time << i << ", " << hist_total[i] / double(SAMPLE_ALL) * double(MOT_ATOMS) / (M_PI * pow(EXIT / 2, 2)) / 1.0 * 1.0E-4 << ", " << average_t << endl;
	}


	/////////////////////////////////
	//  file close
	/////////////////////////////////
	//ofs_output.close();
	//ofs_trace.close();
	ofs_flux.close();
	ofs_time.close();
	power_vs_flux.close();
	hist_time.close();
	loss_log.close();
	/*
	time_log2.close();
	time_log3.close();
	time_log4.close();
	time_log5.close();
	time_log6.close();
	time_log7.close();
	time_log8.close();
	time_log9.close();
	time_log10.close();
	*/

	end = clock();
	std::cout << "計算時間は" << (end - start) / CLOCKS_PER_SEC << "秒" << endl;

	start = clock();
	while (1)
	{
		end = clock();
		if ((end - start) / CLOCKS_PER_SEC > 10)
		{
			break;
		}
	}



}

double f1_x(double v, double r)
{
	return 0;
}
double f1_y(double v, double r)
{
	return 0;
}
double f1_z(double v, double r)
{
	return -GRAVITY;
}
double f2(double v, double r)
{
	return v;
}

double vdw_force(double r)
{
	double sigma = HBAR * (GMM1 / (KJ1 * KJ1 * KJ1) + GMM2 / (KJ2 * KJ2 * KJ2));
	double nn = (N_FUN * N_FUN - 1) / (N_FUN * N_FUN + 1);
	return -nn * sigma * 3 / (16 * r * r * r * r);
}
double intensity_dash(double r, double r_xy, double CENT, double WAIST, double POWER)
{
	double gauss1 = exp(-2.0 * pow(r_xy, 2) / (WAIST * WAIST));
	double gauss2 = exp(-2.0 * pow(r_xy, 2) / (CENT * CENT));
	return 2 * POWER / (M_PI * WAIST * WAIST) * (gauss1 - gauss2) * (2.0 / DECAY) * exp(-(2.0 / DECAY) * r);
}
double fn(double r, double r_xy, int state, double CENT, double WAIST, double POWER, int area, int pol_num)
{
	if (state == 1)
	{
		double delta = DELTA1;
		double gamma = GMM1;
		double i_sat = I_SAT1;
		double s_dash;
		if (pol_num == 0)
		{
			s_dash = intensity_dash(r, r_xy, CENT, WAIST, POWER) * gamma * gamma / (4 * i_sat * delta * delta);  //円偏光
		}
		else
		{
			s_dash = (cos(radians(area * ANGLE_ONEFACE)) * cos(radians(area * ANGLE_ONEFACE))) * intensity_dash(r, r_xy, CENT, WAIST, POWER) * gamma * gamma / (4 * i_sat * delta * delta);  //直線偏光
		}
		return HBAR * delta * s_dash / 2.0 / MASS;
		//return HBAR * delta * s_dash / 2.0 / MASS + vdw_force(r) / MASS;
	}
	else
	{
		double delta = DELTA2;
		double gamma = GMM2;
		double i_sat = I_SAT2;
		double s_dash;
		if (pol_num == 0)
		{
			s_dash = intensity_dash(r, r_xy, CENT, WAIST, POWER) * gamma * gamma / (4 * i_sat * delta * delta);  //円偏光
		}
		else
		{
			
			s_dash = (cos(radians(area * ANGLE_ONEFACE)) * cos(radians(area * ANGLE_ONEFACE))) * intensity_dash(r, r_xy, CENT, WAIST, POWER) * gamma * gamma / (4 * i_sat * delta * delta);  //直線偏光
		}
		return HBAR * delta * s_dash / 2.0 / MASS;
		//return HBAR * delta * s_dash / 2.0 / MASS + vdw_force(r) / MASS;
	}
}
double gn(double vr)
{
	return vr;
}



