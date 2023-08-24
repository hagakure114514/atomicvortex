#pragma once

#include <math.h>
#include "const_param.h"

class atom
{
public:
	position r;
	velocity v;			//並進運動の速度
	velocity v_pre;		//並進運動の速度
	state s;
	double radius;
	double phi;
	double acc_x, acc_y;
	double l_rot = 0.0;		// 重心運動の角運動量 [kg*m^2/s]
	double E_kin = 0.0;

	atom();
	atom(position r0, velocity v0, state s0);
	atom(position r0, velocity v0, double l_rot0, state s0);
};

