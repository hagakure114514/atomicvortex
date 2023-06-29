//Atomic Funnle simulation program
#include <iostream.h>
#include <fstream.h>
#include <math.h>
#include <stdlib.h>
#include <iomanip.h>
#include <time.h>
#include <strstrea.h>
#include <stdio.h>

// doughnut , Note "intensity" and "intensity dash" 

#define EPS 1.0E-10
#define G 9.8
#define EVANE 0.780E-6
#define LAMDA 0.780E-6

//#define DECAY (0.25*LAMDA)
//#define DECAY (0.4*LAMDA)
//#define DECAY (0.55*LAMDA)
#define DECAY (0.7*LAMDA)

#define MASS (85E-3/6.022E+23)
#define HBAR 1.0545E-34
#define BOLTZ 1.380658E-23

//#define HEIGHT 2.0E-3
//#define HEIGHT 3.0E-3
#define HEIGHT 2.0E-3

#define TEMP 10.0E-6
//#define TEMP 50.0E-6
//#define TEMP 100.0E-6
//#define TEMP 150.0E-6

#define TIME_MAX 10
#define SAMPLE 1000

#define EXIT 2.0E-6  /* radius  */
//#define EXIT 20.0E-6
//#define EXIT 200.0E-6

#define OUT (1.224744871*EXIT)

#define POWER 1000.0E-3
//#define POWER 800.0E-3
//#define POWER 600.0E-3
//#define POWER 400.0E-3

#define CENT (EXIT/2.0)    /*  radius? diameter ?*/
//#define CENT (1800.0E-6/2.0)

//#define WAIST 5.0E-3   /*  maybe radius  */
//#define WAIST 4.0E-3
#define WAIST 500.0E-6
//#define WAIST 2.0E-3    
//#define WAIST 1.0E-3

#define M_PI 3.14
#define Q1 0.75
#define Q2 0.61
#define GAMMA0 6.0E6
#define GAMMA1 (2*M_PI*Q1*GAMMA0)
#define GAMMA2 (2*M_PI*(1-Q1)*GAMMA0)
#define I_SAT1 (Q1*1.58E1)
#define I_SAT2 ((1-Q1)*1.58E1)
#define D_HFS (2.0*M_PI*3036.0E6)

#define DELTA1 (2.0*M_PI*1000.0E6)
//#define DELTA1 (2.0*M_PI*500.0E6)

#define DELTA2 (DELTA1+D_HFS)

#define GMM1 (2.0*M_PI*5.58E6)
#define GMM2 (2.0*M_PI*GAMMA0)
#define KJ1 (2.0*M_PI/795.0E-9)
#define KJ2 (2.0*M_PI/780.2E-9)

//#define N_FUN 1.45  //glass
//#define N_FUN 1.51
#define N_FUN 2.5 //TiO2

#define PK (HBAR*2*M_PI/LAMDA/MASS)
#define EPK (PK*N_FUN*0.81614)



double nx[4],ny[4],nz[4];
const double r_max=RAND_MAX;
double fn(double r,double rxys,int state);
double gn(double vr);
void runge_kutta(double &r,double rxys,double &vr,
double &t,double dt,int state);
double intensity(double r,double rxys);
double intensity_dash(double r,double rxys);
double vdw_potential(double r);
double vdw_force(double r);




main()
{
int n,i,j,/*k,*/seed,quitloop/*,check*/;
int nl[10];
int lcount,s1count,s2count,testcount;
int p_num,state;



double gauss;
//double p0,r1,r3;
double x,y,z,vx,vy,vz,t,dt,dm,vv2,vh2,advv2=0.,advh2=0.,Tv,Th,rh,ch;
double x0,y0,z0,atomn=0.,thr=0.;
double r,/*xs,ys,*/zs,rxys,vr/*,vtx,vty,vtz*/;
double *a,*b,*c;
double *va,*vb,*vc;
time_t now;
double theta,phi;
double p,rn,rt;
double tx,ty,tz;
double prev_x,prev_y,prev_z,prev_vx,prev_vy,prev_vz,prev_t;
int prev_state;
char file[80];
//char *temp;
for(i=0;i<=9;i++){
	nl[i]=0;
}



ofstream fout1;
fout1.open("test.dat");
if(!fout1)
{
cerr << "Cannot open the output file.";
exit(-1);
}
fout1.setf(ios::scientific);
fout1<<setprecision(3);

//ofstream fout2;
//fout2.setf(ios::scientific);
//fout2<<setprecision(3);

//ostrstream ostr(file, sizeof file);

x0=HEIGHT/sqrt(3.0);y0=HEIGHT/sqrt(3.0);z0=HEIGHT/sqrt(3.0);

seed=time(&now)%37;
//seed = 1;
srand(seed);

fout1<<"# DECAY="<<DECAY<<endl;
fout1<<"# HEIGHT="<<HEIGHT<<endl;
fout1<<"# TEMP="<<TEMP<<endl;
fout1<<"# EXIT="<<EXIT<<endl;
fout1<<"# POWER="<<POWER<<endl;
fout1<<"# WAIST="<<WAIST<<endl;
fout1<<"# HOLLOW RADIUS"<<CENT<<endl;


for (n=1; n<=SAMPLE; n++) {

ostrstream ostr(file, sizeof file);
ostr << "atom"<< n << ".dat";
ostr << ends;

/* fout2.open(file,ios::out);
if(!fout2)
{
cerr << "Cannot open the output file.";
exit(-1);
}*/

//cout <<n << "st atom" <<endl;
//fout2 << "seed="<<seed<<endl;
quitloop=0;
lcount=0;
s1count=0;
s2count=0;
state=1;
t=0.0;
x=x0+(0.5-rand()/r_max)*pow(10,-3);
y=y0+(0.5-rand()/r_max)*pow(10,-3);
z=z0+(0.5-rand()/r_max)*pow(10,-3);
gauss=0.0;
for (j=1;j<=12;j++) gauss=gauss+rand()/r_max;
gauss=gauss-6.0; vx=sqrt(BOLTZ*TEMP/MASS)*gauss;
gauss=0.0;
for (j=1;j<=12;j++) gauss=gauss+rand()/r_max;
gauss=gauss-6.0; vy=sqrt(BOLTZ*TEMP/MASS)*gauss;
gauss=0.0;
for (j=1;j<=12;j++) gauss=gauss+rand()/r_max;
gauss=gauss-6.0; vz=sqrt(BOLTZ*TEMP/MASS)*gauss;

//cout <<endl;
for(;;){
dt=1.0E-4;//1.0E-5
//cout << "In free space" <<endl;
//fout2 << "_In_free_space" <<endl;
//cout << t << " ";
for(;;){
prev_x=x;prev_y=y;prev_z=z;
prev_vx=vx;prev_vy=vy;prev_vz=vz;
prev_t=t;prev_state=state;
x+=vx*dt-G/sqrt(3.0)*dt*dt/2.0;
y+=vy*dt-G/sqrt(3.0)*dt*dt/2.0;
z+=vz*dt-G/sqrt(3.0)*dt*dt/2.0;
vx+= -G/sqrt(3.0)*dt;
vy+= -G/sqrt(3.0)*dt;
vz+= -G/sqrt(3.0)*dt;
t+=dt;
if ((state==2) &&
(vx/sqrt(3.0)+vy/sqrt(3.0)+vz/sqrt(3.0)>0.0)) {
theta=M_PI*rand()/r_max;
phi=2.0*M_PI*rand()/r_max;

vx+=PK*(sin(theta)*cos(phi)-1.0/sqrt(3.0));

vy+=PK*(sin(theta)*sin(phi)-1.0/sqrt(3.0));
vz+=PK*(cos(theta)-1.0/sqrt(3.0));
state=1;
}
if (t>TIME_MAX){
//cout << "time out" <<endl;
fout1 << "2 "
<< t  << " "
<< lcount << " " << s1count << " "
<< s2count << " "
<< x  << " " << y  << " " << z  << " "
<< vx << " " << vy << " " << vz << " "
<< state <<endl;
quitloop = 1;
break;
}
if ((x-EVANE< 0.0) || (y-EVANE< 0.0) ||
(z-EVANE< 0.0)){
x=prev_x;y=prev_y;z=prev_z;
vx=prev_vx;vy=prev_vy;vz=prev_vz;
t=prev_t;state=prev_state;
dt=dt/5.0;
continue;
}
else if ((x-EVANE < EPS) || (y-EVANE < EPS)
|| (z-EVANE <
EPS)){ p_num=3;
if(y-EVANE < EPS)p_num=2;
if(x-EVANE < EPS)p_num=1;
break;
}
if (x/sqrt(3.0)+y/sqrt(3.0)+z/sqrt(3.0)<=OUT) {

atomn++;
if(lcount==0) thr++;

fout1 << "1 "
<< t  << " "
<< lcount << " " << s1count << " "
<< s2count << " "
<< x  << " " << y  << " " << z  << " "
<< vx << " " << vy << " " << vz << " "
<< state <<endl;

vv2=pow( (fabs(vx) + fabs(vy) + fabs(vz)) , 2 ) / 3.;

vh2=vx*vx+vy*vy+vz*vz-vv2;

advv2+=vv2;

advh2+=vh2; 

//cout << "gain !!" <<endl;
quitloop = 1;
break;
}
//fout2 <<x << " " << y << " " << z <<endl;
}
if(quitloop == 1)break;

//set data for evanesent calc
dt=1.0E-8;//1.0E-10
if(p_num==1){
x=EVANE+2*EPS;
r=x;
vr=vx;
a=&x,b=&y,c=&z;
va=&vx,vb=&vy,vc=&vz;
rxys=sqrt(y*y+z*z);
zs=y/sqrt(3.0)+z/sqrt(3.0);
tx=0.0;
ty=1.0/sqrt(2.0);
tz=1.0/sqrt(2.0);
}
if(p_num==2){
y=EVANE+2*EPS;
r=y;
vr=vy;
a=&y,b=&z,c=&x;
va=&vy,vb=&vz,vc=&vx;
rxys=sqrt(x*x+z*z);
zs=x/sqrt(3.0)+z/sqrt(3.0);
tx=1.0/sqrt(2.0);
ty=0.0;
tz=1.0/sqrt(2.0);
}
if(p_num==3){
z=EVANE+2*EPS;
r=z;
vr=vz;
a=&z,b=&x,c=&y;
va=&vz,vb=&vx,vc=&vy;
rxys=sqrt(x*x+y*y);
zs=x/sqrt(3.0)+y/sqrt(3.0);
tx=1.0/sqrt(2.0);
ty=1.0/sqrt(2.0);
tz=0.0;
}

if(zs < OUT){

atomn++;
if(lcount==0) thr++;

fout1 << "1 "
<< t  << " "
<< lcount << " " << s1count << " " <<
s2count << " "
<< x  << " " << y  << " " << z  << " "
<< vx << " " << vy << " " << vz << " "
<< state <<endl;

vv2=pow( (fabs(vx) + fabs(vy) + fabs(vz)) , 2 ) / 3.;

vh2=vx*vx+vy*vy+vz*vz-vv2;

advv2+=vv2;

advh2+=vh2;

//cout << "gain !!" << endl;
break;
}

//in evanesent field
lcount++;

//cout<<"In "<< lcount <<"th evanesent space"<<endl;
//fout2<<"_In_"<< lcount <<"_th_evanesent_space"<<endl;

//cout<<"time ="<< t <<"sec"<<endl;
testcount=0;
dm=0.0;
for(;;){
//testcount++;
runge_kutta(r,rxys,vr,t,dt,state);
//cout <<testcount<<"th loop"<<" r="<<r<<endl;;
dm+=dt;
if (t>TIME_MAX){
//cout << "time out" <<endl;
fout1 << "2 "
<< t  << " "
<< lcount << " " << s1count << " "
<< s2count << " "
<< x  << " " << y  << " " << z  << " "
<< vx << " " << vy << " " << vz << " "
<< state <<endl;
quitloop = 1;
break;
}
if(r<0.0){
fout1 << "0 "
<< t  << " "
<< lcount << " " << s1count << " "
<< s2count << " "
<< x  << " " << y  << " " << z  << " "
<< vx << " " << vy << " " << vz << " "
<< state <<endl;

//cout << "loss !!" << endl;

rh=sqrt(x*x+y*y+z*z);
ch=0.5E-3;




if(rh<0.5E-3) nl[1]++;
for(i=1;i<=7;i++){

	 if((rh>=ch*i)&&(rh<ch*(i+1))) nl[i+1]++;
}
if(rh>=4.0E-3) nl[9]++;




quitloop = 1;
break;
}
if(vr>0.0)break;
}
if(quitloop == 1)break;

t=t+dm;

//dm=2*dm;
//*b+=vz*dt-G/sqrt(3.0)*dm*dm/2.0;
//*c+=vz*dt-G/sqrt(3.0)*dm*dm/2.0;
//*vb+= -G/sqrt(3.0)*dm;
//*vc+= -G/sqrt(3.0)*dm;

if(state ==1){
rn=rand()/r_max;
p=1.0-exp(MASS*
*va*DECAY*GAMMA0*2.0*M_PI/(HBAR*DELTA1));
//cout << "p=" << p <<endl;
//p=0.0;
//p=1.0;
if (rn<=p) {
//cout << "spontaneous emission"<<endl;
//fout2 <<"_spontaneous_emission"<<endl;
rt=rand()/r_max;
if (rt<=0.25) {
s1count++;
*va=-*va*sqrt(DELTA1/DELTA2);
state=2;
theta=M_PI*rand()/r_max;
phi=2.0*M_PI*rand()/r_max;

vx+=EPK*tx+PK*sin(theta)*cos(phi);

vy+=EPK*ty+PK*sin(theta)*sin(phi);
vz+=EPK*tz+PK*cos(theta);
}
else {
s2count++;
*va=-*va;
theta=M_PI*rand()/r_max;
phi=2.0*M_PI*rand()/r_max;

vx+=EPK*tx+PK*sin(theta)*cos(phi);

vy+=EPK*ty+PK*sin(theta)*sin(phi);
vz+=EPK*tz+PK*cos(theta);
}
}
else {
*va=-*va;
}
}else {
*va=-*va;
}
//fout2 <<x << " " << y << " " << z << " "
// <<vx << " " <<vy << " " <<vz <<endl;
}
//fout2.close();
}
Tv=MASS*advv2/BOLTZ/atomn;
Th=MASS*advh2/BOLTZ/atomn/2.;
fout1 << "gain % " << atomn/SAMPLE*100. <<endl;
fout1 << "through % " << thr/atomn*100. <<endl;
fout1 << "temp vertical " << Tv <<endl;
fout1 << "temp horizontal " << Th <<endl;
fout1 << "loss at r<0.5mm      " << nl[1]<<endl;
fout1 << "loss at 0.5mm<r<1mm  " << nl[2]<<endl;
fout1 << "loss at 1mm<r<1.5mm  " << nl[3]<<endl;
fout1 << "loss at 1.5mm<r<2mm  " << nl[4]<<endl;
fout1 << "loss at 2mm<r<2.5mm  " << nl[5]<<endl;
fout1 << "loss at 2.5mm<r<3mm  " << nl[6]<<endl;
fout1 << "loss at 3mm<r<3.5mm  " << nl[7]<<endl;
fout1 << "loss at 3.5mm<r<4mm  " << nl[8]<<endl;
fout1 << "loss at r>4mm        " << nl[9]<<endl;


fout1.close();
return 0;
}


double fn(double r,double rxys,int state)
{
double delta = (state == 1 ? DELTA1:DELTA2);
double i_sat = (state == 1 ? I_SAT1:I_SAT2);
double gamma = (state == 1 ? GAMMA1:GAMMA2);

double s_dash=intensity_dash(r,rxys)*gamma*gamma/(4*i_sat*delta*delta);
return HBAR*delta*s_dash/2.0/MASS + vdw_force(r)/MASS ;
}
double gn(double vr)
{ return vr;}

void runge_kutta(double &r,double rxys,double &vr,double &t,double dt,int
state)
{
double kn[5],ln[5];
kn[1]=dt*fn(r,rxys,state);
ln[1]=dt*gn(vr);

kn[2]=dt*fn(r+ln[1]/2.0,rxys,state);
ln[2]=dt*gn(vr+kn[1]/2.0);

kn[3]=dt*fn(r+ln[2]/2.0,rxys,state);
ln[3]=dt*gn(vr+kn[2]/2.0);

kn[4]=dt*fn(r+ln[3]/2.0,rxys,state);
ln[4]=dt*gn(vr+kn[3]);

double kkn=(kn[1]+2.0*kn[2]+2.0*kn[3]+kn[4])/6.0;
double lln=(ln[1]+2.0*ln[2]+2.0*ln[3]+ln[4])/6.0;
r+=lln; vr+=kkn; t+=dt;
}

double intensity(double r,double rxys)
{
	double gauss1=exp(-2.0*rxys*rxys/(WAIST*WAIST));
double gauss2=exp(-2.0*rxys*rxys/(CENT*CENT));

return
2*POWER/(M_PI*WAIST*WAIST)*(gauss1-gauss2)*exp(-(2.0/DECAY)*r);
}

double intensity_dash(double r,double rxys)
{
	double gauss1=exp(-2.0*rxys*rxys/(WAIST*WAIST));
double gauss2=exp(-2.0*rxys*rxys/(CENT*CENT));

return
2*POWER/(M_PI*WAIST*WAIST)*(gauss1-gauss2)*exp(-(2.0/DECAY)*r);
}
double vdw_potential(double r)
{
double sigma=HBAR*(GMM1/(KJ1*KJ1*KJ1)+GMM2/(KJ2*KJ2*KJ2));
double nn=(N_FUN*N_FUN-1)/(N_FUN*N_FUN+1);
return -nn*sigma/(16*r*r*r);
}

double vdw_force(double r)
{
double sigma=HBAR*(GMM1/(KJ1*KJ1*KJ1)+GMM2/(KJ2*KJ2*KJ2));
double nn=(N_FUN*N_FUN-1)/(N_FUN*N_FUN+1);
return -nn*sigma*3/(16*r*r*r*r);
}
