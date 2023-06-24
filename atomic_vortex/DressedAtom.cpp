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
	// dissipative force
	force_diss(obj);

	//if(flag_sp !=1){

		// spontaneous emission occur
		if (spontaneous_emission(obj) == 1) {
			detuning_doppler(obj);
			recoil_diss(obj);
			count_sp++;
			//sisyphus_cooling(obj);
		}
	//}
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
double DressedAtom::factorial(int x) {
	double y = 1.0;
	if (x == 0) {
		return y;
	}
	else {
		for (int i = 1; i <= x; i++) {
			y = y * (double)i;
		}
		return y;
	}
}

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
		dUopt = hbar * detuning * log(1.0+s1(obj->radius)) / 2.0 - hbar * detuning * log(1.0+s1(obj->radius_pre)) / 2.0;
	}
	else if (obj->s == state::d2) {
		dUopt = hbar * (detuning + delta_hfs) * log(1.0+s2(obj->radius)) / 2.0 - hbar * (detuning + delta_hfs) * log(1.0+s2(obj->radius_pre)) / 2.0;
	}
	else {
		
		std::mt19937 rand_src(std::random_device{}());
		std::uniform_real_distribution<double> dist(0.0, 1.0);
		double p_diss = dist(rand_src);

		if (p_diss <= branch) {
			dUopt = hbar * detuning * log(1.0+s1(obj->radius)) / 2.0 - hbar * detuning * log(1.0+s1(obj->radius_pre)) / 2.0;
		}
		else {
			dUopt = hbar * (detuning + delta_hfs) * log(1.0+s2(obj->radius)) / 2.0 - hbar * (detuning + delta_hfs) * log( 1.0 +s2(obj->radius_pre)) / 2.0;
		}
	}

	if (obj->radius != obj->radius_pre)
	{
		double f_dip = -dUopt / (obj->radius - obj->radius_pre);
		obj->acc_x += f_dip * cos(obj->phi) / mass;
		obj->acc_y += f_dip * sin(obj->phi) / mass;
	}
}


void DressedAtom::force_diss(atom* obj)
{
	double F_diss = 0.0;
	if (obj->s == state::d1) {
		F_diss = hbar * gamma1 / 2.0 * s1(obj->radius) / (1.0+s1(obj->radius));
	}
	else if (obj->s == state::d2) {
		F_diss = hbar * gamma2 / 2.0 * s2(obj->radius) / (1.0+s2(obj->radius));
	}
	else {
		
		std::mt19937 rand_src(std::random_device{}());
		std::uniform_real_distribution<double> dist(0.0, 1.0);
		double p_diss = dist(rand_src);

		if (p_diss <= branch) {
			F_diss = hbar * gamma1 / 2.0 * s1(obj->radius) / (1.0+s1(obj->radius));
		}
		else {
			F_diss = hbar * gamma2 / 2.0 * s2(obj->radius) / (1.0+s2(obj->radius));
		}
	}

		obj->acc_x += -F_diss * sin(obj->phi) * l/(obj->radius) / mass;
		obj->acc_y += F_diss * cos(obj->phi) * l/(obj->radius) / mass;
		obj->acc_z += F_diss * abs(k_wave) / mass;

}

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

		if(psp<=p1){
			p = 1;
		}
		else {
			p= 0;
		}
	}
	else if (obj->s == state::d2) {
		double p2 = 1.0 - exp(-gamma * s2(obj->radius) * dt / 2.0);

		if (psp <= p2) {
			p = 1;
		}
		else {
			p = 0;
		}
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

	double sp_psi = 2.0 * M_PI * dist(rand_src);			//自然放出の方位角
	double sp_theata = M_PI * dist(rand_src);			//自然放出の仰角

	double p_recoil = dist(rand_src);
	if (p_recoil <= branch) {
		obj->s = state::d1;
	}
	else {
		obj->s = state::d2;
	}

	obj->v.vx += hbar * k_wave * sin(sp_theata) * cos(sp_psi) / mass - hbar * l / (obj->radius) * sin(obj->radius) / mass;
	obj->v.vy += hbar * k_wave * sin(sp_theata) * sin(sp_psi) / mass + hbar * l / (obj->radius) * cos(obj->radius) / mass;
	obj->v.vz += hbar * k_wave  * cos(sp_theata) / mass;

}

// Sisyphus冷却
void DressedAtom::sisyphus_cooling(atom* obj)
{
	std::mt19937 rand_src(std::random_device{}());
	std::uniform_real_distribution<double> dist(0.0, 1.0);

	double psisy = dist(rand_src);			//sisyphus冷却が起こるかのダイス

	double uopt_kin = 0.0;
	double uopt1 = hbar * detuning * s1(obj->radius) / 2.0;
	double uopt2 = hbar * (detuning + delta_hfs) * s2(obj->radius) / 2.0;
	double uopt3 = -uopt1 - uopt2;

	double kinetic_r = 1.0 / 2.0 * mass * ((obj->v.vx) * cos(obj->phi) + (obj->v.vy) * sin(obj->phi)) * ((obj->v.vx) * cos(obj->phi) + (obj->v.vy) * sin(obj->phi));

	if (psisy <=branch) {

		// |1,n-1>への自然放出
		if (obj->s == state:: d1) {
			//運動エネルギーの変化なし
			uopt_kin = 0.0;
		}
		else if(obj->s == state::d2){
			//加熱
			uopt_kin = uopt1 - uopt2;
		}
		else {
			//加熱
			uopt_kin = uopt1 - uopt3;
		}
		obj->s = state::d1;
	}
	else{

		// |2,n-1>への自然放出
		if (obj->s == state::d1) {
			//冷却
			uopt_kin = uopt2 - uopt1;
		}
		else if (obj->s == state::d2) {
			//変化なし
			uopt_kin = 0.0;
		}
		else {
			//加熱
			uopt_kin = uopt2 - uopt3;
		}
		obj->s = state::d2;
	}

	// 運動エネルギーの変化
	obj->v.vx = 0.0;
	obj->v.vy = 0.0;

}
