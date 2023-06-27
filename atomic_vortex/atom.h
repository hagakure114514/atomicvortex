#pragma once

#include <math.h>
#include "const_param.h"

class atom
{
public:
	position r;
	velocity v;
	state s;
	double radius;
	double radius_pre;
	double phi;
	double acc_x, acc_y, acc_z;
	double E_kin = 0.0;

	atom();
	atom(position r0, velocity v0, state s0);
	
};

