#include "DressedAtom.h"

DressedAtom::DressedAtom()
{

	detuning = detuning0;			//detuning between |e> and |g1>
	srand((unsigned)time(NULL));	//a magic phrase to make random numbers be based on time

}

// Optcal potential
void DressedAtom::process_dipole(atom* obj)
{
	// 
	detuning_doppler(obj);
	force_dip(obj);
}

// dissipative force
void DressedAtom::process_diss(atom* obj)
{
	if(flag_sp !=0){

		// spontaneous emission occur
		if (spontaneous_emission(obj) == 1) {
			detuning_doppler(obj);
			recoil_diss(obj);
			count_sp++;
		}
	}
}
	



// motion within a time step dt
void DressedAtom::step_motion(atom* obj)
{
	obj->r.x += 1.0 / 2.0 * obj->acc_x * dt * dt + obj->v.vx * dt;
	obj->r.y += 1.0 / 2.0 * obj->acc_y * dt * dt + obj->v.vy * dt;
	obj->r.z += 1.0 / 2.0 * (obj->acc_z - G) * dt * dt + obj->v.vz * dt;

	obj->v.vx += obj->acc_x * dt;
	obj->v.vy += obj->acc_y * dt;
	obj->v.vz += (obj->acc_z - G) * dt;

	obj->radius_pre = obj->radius;
	obj->radius = sqrt(obj->r.x * obj->r.x + obj->r.y * obj->r.y);
	obj->phi = atan2(obj->r.y, obj->r.x);

	obj->acc_x = 0.0;
	obj->acc_y = 0.0;
	obj->acc_z = 0.0;

}


// factorial function x!
// double DressedAtom::factorial(int x) {
// 	double y = 1.0;
// 	if (x == 0) {
// 		return y;
// 	}
// 	else {
// 		for (int i = 1; i <= x; i++) {
// 			y = y * (double)i;
// 		}
// 		return y;
// 	}
// }

// intensity function of Optical Vortex
double DressedAtom::intensity(double x)
{
	// Optical Vortex Electric Field: E0*sqrt(2)*r/w0 *exp(-r^2/w0^2)
	double I0 = 2.0 * beam_power / (M_PI * w0 * w0);	//beam intensity coefficient [V/m]
	double ints = 2.0*I0*x*x/(w0*w0)*exp(-2.0*x * x / (w0 * w0));			// intensity [V/m]
	return ints;
}

// detuning doppler shift
void DressedAtom::detuning_doppler(atom* obj)
{
	double vphi = -(obj->v.vx) * sin(obj->phi) + (obj->v.vy) * cos(obj->phi);
	detuning = detuning0 - k_wave * (obj->v.vz) - l * prop * vphi / (obj->radius);
}

// saturation parameter between |g1> and |e>
double DressedAtom::s1(double x)
{
	double s = intensity(x) / (I_s1 * (1.0 + 4.0 * detuning * detuning / (gamma1 * gamma1)));
	
	if (s > 0.1) {
		printf("s parameter1 is error %e\n", s);
	}
	
	return s;
}

double DressedAtom::s2(double x)
{
	double s = intensity(x) / (I_s2 * (1.0 + 4.0 * (detuning + delta_hfs) * (detuning + delta_hfs) / (gamma2 * gamma2)));

	if (s > 0.1) {
		printf("s parameter2 is error %e\n", s);
	}

	return s;
}

double DressedAtom::grad_s1(double x)
{
	double I0 = 2.0 * beam_power / (M_PI * w0 * w0);	//最大ビーム強度*1/e [V/m]
	double grad_int = 4.0 * I0/(w0*w0) *(x - 2.0*x * x*x / (w0 * w0)) * exp(-2 * x * x / (w0 * w0));
	double gs = grad_int /(I_s1 *(1 + 4 * detuning * detuning / (gamma1*gamma1)));
	return gs;
}

double DressedAtom::grad_s2(double x)
{
	double I0 = 2.0 * beam_power / (M_PI * w0 * w0);	//最大ビーム強度*1/e [V/m]
	double grad_int = 4.0 * I0 / (w0 * w0) * (x - 2.0 * x * x * x / (w0 * w0)) * exp(-2 * x * x / (w0 * w0));
	double gs = grad_int /(I_s2* (1 + 4 * (detuning + delta_hfs) * (detuning + delta_hfs) / (gamma2*gamma2)));
	return gs;
}

// 双極子ポテンシャルから受ける力
void DressedAtom::force_dip(atom* obj)
{
	double dUopt = 0.0;
	if (obj->s == state::d1) {
		dUopt = 2.0/3.0 * hbar * detuning / 2.0 * ( log(1.0+s1(obj->radius)) - log(1.0+s1(obj->radius_pre)) );
	}
	else if (obj->s == state::d2) {
		dUopt = 2.0/3.0 * hbar * (detuning + delta_hfs) / 2.0 * ( log(1.0+s2(obj->radius)) - log(1.0+s2(obj->radius_pre)) );
	}
	else {
		
		std::mt19937 rand_src(std::random_device{}());
		std::uniform_real_distribution<double> dist(0.0, 1.0);
		double p_diss = dist(rand_src);

		if (p_diss <= branch) {
			dUopt = 2.0/3.0 * hbar * detuning / 2.0 * ( log(1.0+s1(obj->radius)) - log(1.0+s1(obj->radius_pre)) );
			obj->s = state::d1;
		}
		else {
			dUopt = 2.0/3.0 * hbar * (detuning + delta_hfs) / 2.0 * ( log(1.0+s2(obj->radius)) - log(1.0+s2(obj->radius_pre)) );
			obj->s = state::d2;
		}
	}

	if (obj->radius != obj->radius_pre)
	{
		double f_dip = -dUopt / (obj->radius - obj->radius_pre);
		obj->acc_x += f_dip * cos(obj->phi) / mass;
		obj->acc_y += f_dip * sin(obj->phi) / mass;
	}
}


// void DressedAtom::force_diss(atom* obj)
// {
// 	double F_diss = 0.0;
// 	if (obj->s == state::d1) {
// 		F_diss = hbar * gamma1 / 2.0 * s1(obj->radius) / (1.0+s1(obj->radius));
// 	}
// 	else if (obj->s == state::d2) {
// 		F_diss = hbar * gamma2 / 2.0 * s2(obj->radius) / (1.0+s2(obj->radius));
// 	}
// 	else {
		
// 		std::mt19937 rand_src(std::random_device{}());
// 		std::uniform_real_distribution<double> dist(0.0, 1.0);
// 		double p_diss = dist(rand_src);

// 		if (p_diss <= branch) {
// 			F_diss = hbar * gamma1 / 2.0 * s1(obj->radius) / (1.0+s1(obj->radius));
// 			obj->s = state::d1;
// 		}
// 		else {
// 			F_diss = hbar * gamma2 / 2.0 * s2(obj->radius) / (1.0+s2(obj->radius));
// 			obj->s = state::d2;
// 		}
// 	}

// 		obj->acc_x += -F_diss * sin(obj->phi) * l/(obj->radius) / mass;
// 		obj->acc_y += F_diss * cos(obj->phi) * l/(obj->radius) / mass;
// 		obj->acc_z += F_diss * abs(k_wave) / mass;

// }

// 自然放出
bool DressedAtom::spontaneous_emission(atom* obj)
{
	std::mt19937 rand_src(std::random_device{}());
	std::uniform_real_distribution<double> dist(0.0, 1.0);

	bool p;
	double psp = dist(rand_src);			//自然放出確率のダイス

	if (obj->s == state::d1) {
		
		//double p1 = gamma * s1(obj->radius) * dt * (double)life_sp / 2.0;
		double p1 = 1.0 - exp(-gamma * s1(obj->radius) * dt / 2.0);

		p = psp<=p1 ? 1 : 0;

	}
	else if (obj->s == state::d2) {
		double p2 = 1.0 - exp(-gamma * s2(obj->radius) * dt / 2.0);

		p = psp<=p2 ? 1 : 0;

	}
	else {
		// |3> は自然放出を起こさない。
		p = 0;

	}

	return p;

}

// 自然放出の反跳
void DressedAtom::recoil_diss(atom* obj)
{
	std::mt19937 rand_src(std::random_device{}());
	std::uniform_real_distribution<double> dist(0.0, 1.0);
	double p_recoil = dist(rand_src);

	switch(flag_sp){
		case 0:		//no momentum change
			obj->v.vx += 0.0;
			obj->v.vy += 0.0;
			obj->v.vz += 0.0;

		case 1:		// wavevector direction
			double k_recoil = (dist(rand_src)>0.5)? k_wave*1.0: k_wave*-1.0;
			obj->v.vx += 0.0;
			obj->v.vy += 0.0;
			obj->v.vz += hbar * k_wave / mass + hbar * k_recoil;

		case 2:		// dipole radiation direction
			double sp_psi = 2.0 * M_PI * dist(rand_src);	//psi around polarization axis
			double sp_theata = M_PI * dist(rand_src);

			obj->v.vx += 0.0;
			obj->v.vy += 0.0;
			obj->v.vz += hbar * k_wave / mass + hbar * k_wave;

		default:
			double sp_psi = 2.0 * M_PI * dist(rand_src);			//自然放出の方位角
			double sp_theata = M_PI * dist(rand_src);			//自然放出の仰角

			obj->v.vx += hbar * k_wave * sin(sp_theata) * cos(sp_psi) / mass - hbar * l / (obj->radius) * sin(obj->radius) / mass;
			obj->v.vy += hbar * k_wave * sin(sp_theata) * sin(sp_psi) / mass + hbar * l / (obj->radius) * cos(obj->radius) / mass;
			obj->v.vz += hbar * k_wave / mass + hbar * k_wave  * cos(sp_theata) / mass;

	}

	if (obj->s == state::d1) {
		if (p_recoil <= branch) {
			// |1,n-1>への自然放出
			obj->s = state::d1;
		}else{
			// |2,n-1>への自然放出
			// Sisyphus冷却
			obj->s = state::d2;

			double uopt1 = 2.0/3.0 * hbar * detuning / 2.0 * log(1.0+s1(obj->radius));
			double uopt2 = 2.0/3.0 * hbar * (detuning + delta_hfs) / 2.0 * log(1.0+s2(obj->radius));
			double Kin_r = ;

			// 運動エネルギーの変化
			obj->v.vx = 0.0;
			obj->v.vy = 0.0;
		}
	}else{
		if (p_recoil <= branch) {
			// |1,n-1>への自然放出
			obj->s = state::d1;
		}else{
			// |2,n-1>への自然放出
			obj->s = state::d2;
		}
	}
	

}