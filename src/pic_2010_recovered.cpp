// Recovered research code from Appendix B of Ting-Yi Lin's 2010 master's thesis.
// Recovery/portability edits only: repaired PDF line wraps, added M_PI fallback,
// and disabled the Windows-only system("PAUSE") call. Scientific logic is unchanged.
// See docs/PROVENANCE.md and docs/RECOVERY_NOTES.md.

#include<iostream>
using namespace std;
#include<math.h>
#include<fstream>
#include <stdio.h>
#include <time.h> /*引述 time() 函數宣告*/
#include <sys/timeb.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif /*引述 struct timeb 結構定義與 ftime() 函數宣告*/
double L,Dx,Dt,total_times,Nd;//read parameters
int AP,NSTEPS,KE,WT,P_distribution,P_collision,P_weighting,re_loop;
double PPC;
//-------------------------------------------plasma space
double *NT;
double *Ex,*potential;
double *Exh,*q;
double *N,*ni,*ne;
//-------------------------------------------plasma space
//-------------------------------------------particle data
double *pEx,*pEy,*pEz;
double *vx,*vy,*vz,*v_perpendicular,*total_v;
double *x,*hlx,*nu;
//-------------------------------------------particle data
//-------------------------------------------datat 分析
double *totalka,*totalkh,*totalE,*totalEF,*totalP,*FTC,*FTS,*FT,*data_r;
double *Ha4,*Ha6,*Ha8,*Ha10;
double *w_time,*log_Ha4,*log_Ha6,*log_Ha8,*log_Ha10;
//-------------------------------------------datat 分析
double *cell_x,*cell_t;
#define window_x 200
#define window_y 1000
#define Thermal_wind_x 200
double hermal_velocity,a_deby_lenth,a_plasma_frequency,thermal_velocity;
const double
epaz(8.854187717e-12),e(-1.602176487e-19),u0(1.256637061e-6),em(9.10938188e-31),c0(

299792458.0),tk(1.3806503e-23);
//M_PI=pi(3.141592654)
double dx,dt,wave_length,n0,A0,Kp2,J,b;
int n,T,i,Ei,a,y,k,kc;
double a_wave_cells;
double Vmax=0,Vmin=0;
double window[window_y][window_x],window2[window_y][window_x];
double Twindow[window_y][Thermal_wind_x];
int walx,walt,wst,wsx,jumpx,jumpt,wt=0;
double Amplitude,Temperature1,k0,plasma_lenght;
double B,k2;
double acoo_velocity1,acoo_velocity2;
double acoo_k1,acoo_k2;
int acoo_time1,acoo_time2;
//-----------------------------------------------collision
double Collision_gama,Collision_nw_gama;
double c_arfa; //PAS constan
double RMK;//-----------Monte carlo constan =(1-exp(-gamat))
//-------------------------------------
double gg[Thermal_wind_x],ggx[Thermal_wind_x];
double thermal_v[Thermal_wind_x],thermal_x[Thermal_wind_x];
//--------------------------------------
//-------------------------------------------------------BL collision
double BLeta_v,BLeta_x,BLcrv,BLrdt;
double BLc0,BLc1,BLc2,BLc3;
//-------------------------------------------------------BL collision
//-------------------------------------------datat 分析
double rem4[100][1000],rm_a4[1000],rm_s4[1000];
double rem6[100][1000],rm_a6[1000],rm_s6[1000];
double rem8[100][1000],rm_a8[1000],rm_s8[1000];
double rem10[100][1000],rm_a10[1000],rm_s10[1000];

//-------------------------------------------datat 分析
void output(string D,int end,double EC[]);
void Input_parameter();
inline void first();//----------------------give initial parmeter
//----------------------------------------------------give a random number
const int IM1=2147483563, IM2=2147483399;
const int IA1=40014, IA2=40692, IQ1=53668, IQ2=52774;
const int IR1=12211, IR2=3791, NTAB=32, IMM1=IM1-1;
const int NDIV=1+IMM1/NTAB;
const double EPS=3.0e-16, RNMX=(1.0-EPS), AM=1.0/double(IM1);
static int idum2=123456789;
static int iy=0;
double n_IQ1=1.0/double(IQ1),n_IQ2=1.0/double(IQ2),n_NDIV=1.0/double(NDIV);
double idum=-time(NULL);
inline double ran2()
{
static int iv[NTAB];
int j, k;
double temp;
if (idum <=0) {
idum=(idum==0 ? 1: -idum);
idum2=idum;
for (j=NTAB+7; j >= 0; j--){
k=idum*n_IQ1;
idum=IA1*(idum-k*IQ1)-k*IR1;
if (idum < 0) idum += IM1;
if (j < NTAB) iv[j]=idum;
}
iy=iv[0];
}
k=idum*n_IQ1;
idum=IA1*(idum-k*IQ1)-k*IR1;
if (idum < 0) idum += IM1;

k=idum2*n_IQ2;
idum2=IA2*(idum2-k*IQ2)-k*IR2;
if(idum2 < 0) idum2 +=IM2;
j=iy*n_NDIV;
iy=iv[j]-idum2;
iv[j]=idum;
if (iy < 1) iy +=IMM1;
if ((temp=AM*iy) > RNMX) return RNMX;
else return temp;
}
//----------------------------------------------------give a random number
//--------------------------------------give different velocity distribution
inline void give_temperature_Linux(double v[],double temperature)
{
int randan[3],aa,dd;
double r[12],d,sum;
d=thermal_velocity;
for(i=0;i<Ei;i++)
{
sum =0.0;
for(k=0;k<12;k++)
{
sum=ran2()+sum;
}
v[i]=d*(sum-6.0);
}
}
inline void give_temperature_U(double v[])
{
int randan[3],aa,dd;
double r[12],d,sum;
d=thermal_velocity;

for(i=0;i<Ei;i++)
{
v[i]=sqrt(3.0)*2.0*d* (ran2()-0.5);
}
}
inline void ball_distribution()
{
int randan[3],aa,dd;
double r[12],d,sum;
d=thermal_velocity;
for(i=0;i<Ei;i++)
{
vx[i]=2.0*(ran2()-0.5)*sqrt(5.0);
vy[i]=2.0*(ran2()-0.5)*sqrt(5.0);
vz[i]=2.0*(ran2()-0.5)*sqrt(5.0);
total_v[i]=sqrt(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]);
for(;sqrt(5.0)<total_v[i];)
{
vx[i]=2.0*(ran2()-0.5)*sqrt(5.0);
vy[i]=2.0*(ran2()-0.5)*sqrt(5.0);
vz[i]=2.0*(ran2()-0.5)*sqrt(5.0);
total_v[i]=sqrt(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]);
}
}
for(i=0;i<Ei;i++)
{
vx[i]=d*vx[i];
vy[i]=d*vy[i];
vz[i]=d*vz[i];
}
}

//--------------------------------------give different velocity distribution
//-------------------------------------統計速度空間的粒子數目(固定上下限)
inline void fix_temperature_view(double vt[],double number_of_T[],double velocity_local[],double vmax,double vmin)
{
double vm,vd;
vd=vmax-vmin;
for(k=0;k<Thermal_wind_x;k++)
{
number_of_T[k]=0.0;
velocity_local[k]=0.0;
}
for(i=0;i<Ei;i++)
{if(vt[i]<vmax && vt[i]>vmin)
number_of_T[int((vt[i]-vmin)*Thermal_wind_x/vd)]++;
}
for(k=0;k<Thermal_wind_x;k++)
{
velocity_local[k]=k*vd/Thermal_wind_x+vmin;
}
}
//-------------------------------------統計速度空間的粒子數目(固定上下限)
//-------------------------------------統計速度空間的粒子數目(不固定上下限)
inline void temperature_view(double v[],double number_of_T[],double velocity_local[])
{
double vmax=0,vmin=0,vm,vd;
for(k=0;k<Thermal_wind_x;k++)
{
number_of_T[k]=0.0;
velocity_local[k]=0.0;
}
for(i=0;i<AP;i++)
{
if(vmax<v[i])

vmax=v[i];
if(vmin>v[i])
vmin=v[i];
}
vd=vmax-vmin;
for(i=0;i<Ei;i++)
{
number_of_T[int((v[i]-vmin)*Thermal_wind_x/vd)]++;
}
for(k=0;k<Thermal_wind_x;k++)
{
velocity_local[k]=k*vd/Thermal_wind_x+vmin;
}
}
//-------------------------------------統計速度空間的粒子數目(不固定上下限)
//-------------------------------------給速度和粒子大小的初始值(歸零)
inline void ivelocity()
{
for(i=0;i<AP;i++)
{
vx[i]=vy[i]=0;
nu[i]=1.0;
}
}
//-------------------------------------給速度和粒子大小的初始值 (歸零)
//-------------------------------------給粒子位子的初始值 (歸零)
inline void particle_give()
{
Ei=AP;
for(i=0;i<Ei;i++)
{
x[i]=double(i)*1.0/PPC+1.0;
}

}
//-------------------------------------給粒子位子的初始值 (歸零)
//----------------------------計算隔點上應有的電荷密度(使用 first weighting)
inline void first_calculate_density()
{
for(k=0;k<KE;k++)
{
N[k]=0;
}
for(i=0;i<Ei;i++)
{
a=int(x[i]/1.0);
b=x[i]-double(a);
N[a]=(1.0-b)*nu[i]+N[a];
N[a+1]=b*nu[i]+N[a+1];
}
N[1]=N[1]+N[KE-2];
N[KE-2]=N[1];
for(k=0;k<KE;k++)
{
ne[k]=e*n0*N[k]/PPC;
q[k]=ni[k]+ne[k];
}
}
//----------------------------計算隔點上應有的電荷密度(使用 first weighting)
//----------------------------計算隔點上應有的電荷密度(使用 second weighting)
inline void second_calculate_density()
{
for(k=0;k<KE;k++)
{
N[k]=0;
}
for(i=0;i<Ei;i++)
{
a=int(x[i]);
b=x[i]-double(a);

if(b>0.5)
{
N[a]=N[a]+(0.5*(1.5-b)*(1.5-b))*nu[i];
a++;
b=1.0-b;
N[a]=N[a]+(0.75-b*b)*nu[i];
a++;
b++;
N[a]=N[a]+(0.5*(1.5-b)*(1.5-b))*nu[i];
}
else
{
a--;
b++;
N[a]=N[a]+(0.5*(1.5-b)*(1.5-b))*nu[i];
a++;
b--;
N[a]=N[a]+(0.75-b*b)*nu[i];
a++;
b=1.0-b;
N[a]=N[a]+(0.5*(1.5-b)*(1.5-b))*nu[i];
}
}
//----------------------------------per_BC
N[0]=N[0]+N[KE-3];
N[KE-3]=N[0];
N[1]=N[1]+N[KE-2];
N[KE-2]=N[1];
N[2]=N[2]+N[KE-1];
N[KE-1]=N[2];
//-----------------------------------
for(k=0;k<KE;k++)
{
ne[k]=e*n0*N[k]/PPC;
q[k]=ni[k]+ne[k];

}
}
//----------------------------計算隔點上應有的電荷密度(使用 second weighting)
inline void calculate_potential()
{
double sum;
sum=0.0;
for(k=1;(k<KE-2);k++)
{
sum=q[k]*k+sum;
}
potential[1]=-sum*dx*dx/epaz/(KE-3.0);
potential[2]=-q[1]*dx*dx/epaz+2.0*potential[1];
for(k=2;k<(KE-1);k++)
{
potential[(k+1)]=-q[k]*dx*dx/epaz+2.0*potential[k]-potential[(k-1)];
}
potential[0]= potential[(KE-3)];
potential[(KE-2)]= potential[1];
potential[(KE-1)]= potential[2];
sum=0.0;
for(k=1;(k<KE-2);k++)
{
sum=q[k]+sum;
}
}
inline void calculate_Ex()
{
calculate_potential();
for(k=1;k<(KE-1);k++)

{
Ex[k]=-0.5*(potential[k+1]-potential[k-1])/dx;
}
Ex[0]=Ex[(KE-3)];
Ex[(KE-2)]=Ex[1];
Ex[(KE-1)]=Ex[2];
}
//-----------------將格點(cells)上的電場 weighting 到粒子上(使用 first weighting)
inline void first_wating_particle_f()
{
for(i=0;i<Ei;i++)
{
a=int(x[i]);
b=x[i]-a;
pEx[i]=b*Ex[a+1]+(1.0-b)*Ex[a];
}
}
//-----------------將格點(cells)上的電場 weighting 到粒子上(使用 first weighting)
//-----------------將格點(cells)上的電場 weighting 到粒子上(使用 second weighting)
inline void second_wating_particle_f()
{
double d1,d2,d3;
for(i=0;i<Ei;i++)
{
a=int(x[i]);
b=x[i]-double(a);
if(b>0.5)
{
d1=b;
d2=1.0-d1;
d3=d2+1.0;
}
else

{
a--;
d1=b+1.0;
d2=d1-1.0;
d3=1.0-d2;
}
pEx[i]=0.5*(1.5-d3)*(1.5-d3)*Ex[(a+2)]+(0.75-d2*d2)*Ex[(a+1)]+0.5*(1.5-d1)*(1.5-d1)*Ex[a];
}
}
//-----------------將格點(cells)上的電場 weighting 到粒子上(使用 second weighting)
//----------------使用 Virtamo 提出的方法進行資料分析
inline void Hermite_polynomials()
{
double Hu;
double p2,p4,p6,p8,p10;
double Hp[12];
Ha4[wt]=0.0;
Ha6[wt]=0.0;
Ha8[wt]=0.0;
Ha10[wt]=0.0;
Hp[0]=1.0;
for(i=0;i<Ei;i++)
{
Hu=vx[i]/thermal_velocity;
Hp[1]=Hu;
for(k=1;k<10;k++)
{
Hp[k+1]=Hu*Hp[k]-k*Hp[k-1];
}
Ha4[wt]+=Hp[4];
Ha6[wt]+=Hp[6];
Ha8[wt]+=Hp[8];
Ha10[wt]+=Hp[10];

}
Ha4[wt]=fabs(Ha4[wt])/AP;
Ha6[wt]=fabs(Ha6[wt])/AP;
Ha8[wt]=fabs(Ha8[wt])/AP;
Ha10[wt]=fabs(Ha10[wt])/AP;
}
//----------------使用 Virtamo 提出的方法進行資料分析
//----------------Turner 所使用的碰撞方法
inline void MK_collision()
{
double co2f;
co2f=2.0*RMK;
for(i=0;i<Ei;i++)
{
if(co2f>ran2())
{
vx[i]=-vx[i];
}
}
}
//----------------Turner 所使用的碰撞方法
//----------------小角度散射模型(small-angle-collision model)
inline void PAS_collision()
{
double jj;
double c_lumda_old,c_lumda_new,c_gama;
double c_gama_max(0),c_gama_min(0),c_lumda_old_max(0),c_lumda_old_min(0),c_lumda_new_max(0),c_lumda_new_min(0);
jj=1.0/(sqrt(3.0)*thermal_velocity);
for(i=0;i<Ei;i++)
{
total_v[i]=sqrt(vx[i]*vx[i]+v_perpendicular[i]*v_perpendicular[i]);
c_gama=Collision_gama*pow((tanh(c_arfa*jj*total_v[i])/(jj*total_v[i])),3.0);
c_lumda_old=vx[i]/total_v[i];

if(ran2()>0.5)
{
c_lumda_new=c_lumda_old*(1.0-c_gama*dt)+sqrt(fabs((1.0-pow(c_lumda_old,2.0)))*c_gama*dt);
}
else
{
c_lumda_new=c_lumda_old*(1.0-c_gama*dt)-sqrt(fabs((1.0-pow(c_lumda_old,2.0)))*c_gama*dt);
}
vx[i]=total_v[i]*c_lumda_new;
v_perpendicular[i]=total_v[i]*sqrt(1.0-c_lumda_new*c_lumda_new);
}
}
//----------------小角度散射模型(small-angle-collision model)
//----------------邊界的處理(使用週期性邊界)
inline void take_out_particle()
{
for(i=0;i<Ei;i++)
{
if(x[i]<1)
{x[i]=x[i]+KE-3.0;}
if(x[i]>=(KE-2))
{x[i]=x[i]-KE+3.0;}
}
}
//----------------邊界的處理(使用週期性邊界)
//---------------------------eq of motion(算出粒子受力、速度和位子)
inline void eq_of_motion()
{
take_out_particle();
if(P_weighting==1)
{

first_calculate_density();
}
if(P_weighting==2)
{
second_calculate_density();
}
calculate_Ex();
if(P_weighting==1)
{
first_wating_particle_f();
}
if(P_weighting==2)
{
second_wating_particle_f();
}
for(i=0;i<Ei;i++)
{
vx[i]=vx[i]+pEx[i]*dt*e/em;
x[i]=x[i]+dt*vx[i]/dx;
}
if(P_collision==1)
{
MK_collision();
}
if(P_collision==2)
{
PAS_collision();
}
}
//---------------------------eq of motion(算出粒子受力、速度和位子)

//----------------有 BL collision 模型的 eq of motion(算出粒子受力、速度和位子)
inline void BL_mode_motion()
{
double BLrandan_x,BLrandan_v,BLrandan1,BLrandan2,sum,s1,s2,r1,r2;
double BLfast1,BLfast2,BLadt,BL_vh;
BLfast1=BLc0*BLc2/BLc1;
BLfast2=(1.0-BLc0/BLc1)/BLrdt;
for(i=0;i<Ei;i++)
{
do{
s1 = ran2();
s2 = ran2();
}while(s1 == 0.0);
r1 = sqrt(-2.0 * log(s1)) * cos(2.0 * M_PI * s2);
r2 = sqrt(-2.0 * log(s1)) * sin(2.0 * M_PI * s2);
BLrandan1=r1;
BLrandan2=r2;
BLrandan_x=BLeta_x*BLrandan1;
BLrandan_v=BLeta_v*(BLcrv*BLrandan1+sqrt(1.0-BLcrv*BLcrv))*BLrandan2;
BLadt=pEx[i]*e*dt/em;
x[i]=x[i]+(BLc1*vx[i]*dt+BLc2*BLadt*dt+BLrandan_x)/dx;
vx[i]=BLc0*vx[i]+(BLc1-BLc2)*BLadt+BLrandan_v;
}
//******************************************
take_out_particle();
second_calculate_density();
calculate_Ex();
second_wating_particle_f();
//***************************************
for(i=0;i<Ei;i++)
{
BLadt=pEx[i]*e*dt/em;
vx[i]=vx[i]+BLc2*BLadt;

}
}
//----------------有 BL collision 模型的 eq of motion(算出粒子受力、速度和位子)
//----------------取 data(可用可不用)
void watching_window(double NE[])
{
int GG;
for(i=0;i<window_x;i++)
{
GG=wsx+i*jumpx;
window[wt][i]=NE[GG];
}
}
void watching_Twindow()
{
int GG;
fix_temperature_view(vx,gg,ggx,(3.0*thermal_velocity),(-3.0*thermal_velocity));
for(i=0;i<Thermal_wind_x;i++)
{
Twindow[wt][i]=gg[i];
gg[i]=0;
}
//cout<<"\n"<<T;
}
//----------------取 data(可用可不用)
void output2D (string D,int endy,int endx,double AC[window_y][window_x])
{
ofstream A2(D.c_str());
for(i=0;i<endy;i++)
{
for(k=0;k<endx;k++)
{A2 << AC[i][k]<<" ";}
A2 <<"\n";
}

A2.close();
}
void output_thermal (string D,int endy,int endx,double AC[window_y][Thermal_wind_x])
{
ofstream A2(D.c_str());
for(i=0;i<endy;i++)
{
for(k=0;k<endx;k++)
{A2 << AC[i][k]<<" ";}
A2 <<"\n";
}
A2.close();
}
//-------------做線性迴歸
void Lin_REG(double d_x[],double d_y[],double s_y[],int Number,double *a1,double *a2,double *sd_L)
{
double sum,sum_x,sum_y,sum_x2,sum_xy,sum_sd,average_x,average_y;
double r,ru,rdl,rdr;
sum=0.0;
sum_x=0.0;
sum_y=0.0;
sum_x2=0.0;
sum_xy=0.0;
sum_sd=0.0;
ru=0.0;
rdl=0.0;
rdr=0.0;
for(k=0;k<Number;k++)
{
sum+=1.0/(s_y[k]*s_y[k]);
sum_x+=d_x[k]/(s_y[k]*s_y[k]);
sum_y+=d_y[k]/(s_y[k]*s_y[k]);
sum_x2+=(d_x[k]*d_x[k])/(s_y[k]*s_y[k]);
sum_xy+=(d_x[k]*d_y[k])/(s_y[k]*s_y[k]);
}

*a1=(sum_y*sum_x2-sum_x*sum_xy)/(sum*sum_x2-sum_x*sum_x);//截距
*a2=(sum*sum_xy-sum_y*sum_x)/(sum*sum_x2-sum_x*sum_x);//斜率
average_x=0.0;
average_y=0.0;
for(k=0;k<Number;k++)
{
average_x+=d_x[k];
average_y+=d_y[k];
}
average_x=average_x/Number;
average_y=average_y/Number;
for(k=0;k<Number;k++)
{
ru+=(d_x[k]-average_x)*(d_y[k]-average_y);
rdl+=(d_x[k]-average_x)*(d_x[k]-average_x);
rdr+=(d_y[k]-average_y)*(d_y[k]-average_y);
}
r=ru/sqrt(rdl*rdr);//相關係數
*sd_L=sqrt(sum/(sum*sum_x2-sum_x*sum_x));//斜率標準差
}
//-------------做線性迴歸
//-------------算標準差
void Ave_St(double d_x[],int Number,double *a1,double *a2)
{
double sum_x=0.0,st=0.0;
for(k=0;k<Number;k++)
{
sum_x+=d_x[k];
}

*a1=sum_x/Number;
for(k=0;k<Number;k++)
{
st+=(d_x[k]-*a1)*(d_x[k]-*a1);
}
*a2=sqrt(st/(double(Number)-1.0));
}
//-------------算標準差
//-------------------------------------------------------main function
int main()
{
double Ex0;
int ss;
first();
ivelocity();
struct timeb tb;
int IT,final_time;
ftime(&tb);
IT= tb.time;
double a_4[100],a_6[100],a_8[100],a_10[100],Q;
double geta4,geta6,geta8,geta10;
int a4_end,a6_end,a8_end,a10_end;
int pass;
//---------- pass 總共要跑的次數
for(pass=0;pass<re_loop;pass++)
{
particle_give();
//////////////////////////////////////////////////give phas distribution
if(Temperature1!=0)
{
if (P_distribution==0)
{

give_temperature_Linux(vx,Temperature1);
give_temperature_Linux(vy,Temperature1);
give_temperature_Linux(vz,Temperature1);
}
if (P_distribution==1)
{
give_temperature_U(vx); //----------------------
give_temperature_U(vy);
give_temperature_U(vz);
}
if (P_distribution==2)
{
ball_distribution();
}
}
//////////////////////////////////////////////////give phas distribution
for(i=0;i<Ei;i++)
{
v_perpendicular[i]=sqrt(vy[i]*vy[i]+vz[i]*vz[i]);
}
//---------------------ion fix*/
for(k=0;k<KE;k++)
{
ni[k]=-e*n0;
}
//---------------------ion fix*/
second_calculate_density();
calculate_Ex();
second_wating_particle_f();
cout<<"thermal_velocity"<<thermal_velocity<<"\n";
wt=0;
T=0;
//--------------------一個 case 要跑幾個 time step
for(n=0;n<NSTEPS;n++)
{
T++;

cout<<"loop="<<(pass+1)<<" "<<"time_step="<<T<<"\n";
eq_of_motion();
// BL_mode_motion();
if(n%jumpt==0&&n<(window_y*jumpt-1))
{
Hermite_polynomials();
watching_window(q);
if(Temperature1!=0)
watching_Twindow();
wt++;
}
//-----------------------------------------K
for(i=0;i<AP;i++)
{
if(Vmax<vx[i])
Vmax=vx[i];
if(Vmin>vx[i])
Vmin=vx[i];
totalka[n]=0.5*em*nu[i]*n0*dx/PPC*(vx[i]*vx[i]+vy[i]*vy[i])+totalka[n];
totalP[n]=nu[i]*vx[i]+totalP[n];
}
//----------算總電場能量
for(k=1;k<KE-2;k++)
{
totalEF[n]=0.5*epaz*Ex[k]*Ex[k]*dx+totalEF[n];
}
//----------算總電場能量
//----------算總能(因為算完 V.x.後再算電場所以 Ex 會比速度快半個時間隔點)
if(n>0)
totalE[n]=0.5*totalEF[n]+0.5*totalEF[n-1]+totalka[n];
//----------算總能(因為算完 V.x.後再算電場所以 Ex 會比速度快半個時間隔點)
//----------對系統中最長的波長做傅立葉變換（Fourier Transform)

for(k=0;k<KE-3;k++)
{
FTS[n]=Ex[k]*sin(2.0*M_PI*double(k)*1.0/(KE-3))+FTS[n];
FTC[n]=Ex[k]*cos(2.0*M_PI*double(k)*1.0/(KE-3))+FTC[n];
FT[n]=FTS[n]*FTS[n]+FTC[n]*FTC[n];
}
//----------對系統中最長的波長做傅立葉變換（Fourier Transform)
}
//--------------------一個 case 要跑幾個 time step
cout<<"thermal_velocity"<<thermal_velocity<<"\n";
cout<<"\n"<<"Vmax="<<Vmax<<" "<<"Vmin="<<Vmin;
for(i=0;i<(NSTEPS/jumpt);i++)
{
w_time[i]=(double(i)*jumpt/a_plasma_frequency);
rem4[pass][i]=log(Ha4[i]);
rem6[pass][i]=log(Ha6[i]);
rem8[pass][i]=log(Ha8[i]);
rem10[pass][i]=log(Ha10[i]);
}
}
//---------- pass 總共要跑的次數
//------------------------------------資料分析，做回歸和算標準差
double NNO[100],aa1,aa2;
for(i=0;i<(NSTEPS/jumpt);i++)
{
for(k=0;k<pass;k++)
{
NNO[k]=rem4[k][i];
}
Ave_St(NNO,pass,&aa1,&aa2);

rm_a4[i]=aa1;
rm_s4[i]=aa2;
for(k=0;k<pass;k++)
{
NNO[k]=rem6[k][i];
}
Ave_St(NNO,pass,&aa1,&aa2);
rm_a6[i]=aa1;
rm_s6[i]=aa2;
for(k=0;k<pass;k++)
{
NNO[k]=rem8[k][i];
}
Ave_St(NNO,pass,&aa1,&aa2);
rm_a8[i]=aa1;
rm_s8[i]=aa2;
for(k=0;k<pass;k++)
{
NNO[k]=rem10[k][i];
}
Ave_St(NNO,pass,&aa1,&aa2);
rm_a10[i]=aa1;
rm_s10[i]=aa2;
}
double sd_a4,sd_a6,sd_a8,sd_a10;
Lin_REG(w_time,rm_a4,rm_s4,(NSTEPS/jumpt),&Q,&geta4,&sd_a4);
Lin_REG(w_time,rm_a6,rm_s6,(NSTEPS/jumpt),&Q,&geta6,&sd_a6);
Lin_REG(w_time,rm_a8,rm_s8,(NSTEPS/jumpt),&Q,&geta8,&sd_a8);
Lin_REG(w_time,rm_a10,rm_s10,(NSTEPS/jumpt),&Q,&geta10,&sd_a10);
//------------------------------------資料分析，做回歸和算標準差
//--------------計時器(計算總共花了多少秒)
ftime(&tb);
final_time= tb.time;
//--------------計時器(計算總共花了多少秒)
//---------------------輸出 data

ofstream G1("space.txt");
G1<<"cell_x"<<" "<<"ni"<<" "<<"ne"<<" "<<"q"<<" ";
G1<<"Ex"<<" "<<"potential"<<"\n";
for(i=0;i<KE;i++)
{
G1<<cell_x[i]<<" "<<ni[i]<<" "<<ne[i]<<" "<<q[i]<<" ";
G1<<Ex[i]<<" "<<potential[i]<<"\n";
}
G1.close();
ofstream G2("time.txt");
G2<<"cell_t"<<" "<<"FTC"<<" "<<"FTS"<<" "<<"FT"<<" ";
G2<<"totalEF"<<" "<<"totalE"<<" "<<"totalka"<<" "<<"totalP"<<"\n";
for(i=0;i<NSTEPS;i++)
{
G2<<cell_t[i]<<" "<<FTC[i]<<" "<<FTS[i]<<" "<<FT[i]<<" ";
G2<<totalEF[i]<<" "<<totalE[i]<<" "<<totalka[i]<<" "<<totalP[i]<<"\n";
}
G2.close();
/* ofstream G3("Hp.txt");
G3<<"cell_t"<<" "<<"Ha4"<<" "<<"Ha6"<<" "<<"Ha8"<<" ";
G3<<"Ha10"<<"\n";
for(i=0;i<(NSTEPS/jumpt);i++)
{
G3<<(double(i)*jumpt/a_plasma_frequency)<<" "<<Ha4[i]<<" "<<Ha6[i]<<"
"<<Ha8[i]<<" ";
G3<<Ha10[i]<<"\n";
}
G3.close();*/
ofstream G4("Hermite_polynomials.txt");
G4<<"w_time"<<" "<<"ln(H4)"<<" "<<"standard_ln(H4)"<<" "<<"ln(H6)"<<" "<<"standard_ln(H6)"<<" "<<"ln(H8)"<<" "<<"standard_ln(H8)"<<" "<<"ln(H10)"<<" "<<"standard_ln(H10)"<<"\n";
for(i=0;i<(NSTEPS/jumpt);i++)
{
G4<<w_time[i]<<" "<<rm_a4[i]<<" "<<rm_s4[i]<<" "<<rm_a6[i]<<" "<<rm_s6[i]<<" "<<rm_a8[i]<<" "<<rm_s8[i]<<" "<<rm_a10[i]<<" "<<rm_s10[i]<<"\n";
}
G4.close();
ofstream all("fitting_number(omega_p*tau_R).txt",ios::app);
all<<"H4"<<" "<<"standard(H4)"<<" "<<"H6"<<" "<<"standard(H6)"<<" "<<"H8"<<" "<<"standard(H8)"<<" "<<"H10"<<" "<<"standard(H10)"<<"\n";
all<<(-1.0/(geta4/4.0))<<" "<<(1.0/(geta4/4.0)*sd_a4/geta4)<<" ";
all<<(-1.0/(geta6/6.0))<<" "<<(1.0/(geta6/6.0)*sd_a6/geta6)<<" ";
all<<(-1.0/(geta8/8.0))<<" "<<(1.0/(geta8/8.0)*sd_a8/geta8)<<" ";
all<<(-1.0/(geta10/10.0))<<" "<<(1.0/(geta10/10.0)*sd_a10/geta10)<<"\n";
all.close();
ofstream G331("output.txt");
G331<<"------------------input file---------"<<"\n";
G331<<"L="<<L<<"\n";
G331<<"Dx="<<Dx<<"\n";
G331<<"Dt="<<Dt<<"\n";
G331<<"total_times="<<total_times<<"\n";
G331<<"Nd="<<Nd<<"\n";
G331<<"P_distribution="<<P_distribution<<"\n";
G331<<"P_collision_mode="<<P_collision<<"\n";
G331<<"collision_rate="<<Collision_nw_gama<<"\n";
G331<<"collision_alpha="<<c_arfa<<"\n";
G331<<"P_weighting="<<P_weighting<<"\n";
G331<<"pass_loop="<<pass<<"\n";
G331<<"------------------input file---------"<<"\n";
G331<<"------------------used time---------"<<"\n";
G331<<"used time="<<(final_time-IT)<<" second"<<"\n";
G331<<"------------------used time---------"<<"\n";
G331.close();
//---------------------輸出 data
cout<<"\n"<<(final_time-IT);

/* system("PAUSE") */;//時間暫停
}
void output(string D,int end,double EC[])
{
ofstream Amplitude(D.c_str());
for(k=0;k<end;k++)
{Amplitude << EC[k]<<"\n";}
Amplitude.close();
}
void Input_parameter()
{
//----------------------------------------------------------------------只讀數字的部分
char c[1000],r[1000];
double number,parameter[50];
int st,i,dd(0),number_star,number_end,k;
ifstream read("input.txt");
while(!read.eof() && dd<40)
{
number_star=0;
for(i=0;i<100;i++)
{
c[i]=0;
r[i]=0;
}
number_star=0;
number_end=0;
i=0;
k=0;
read.getline(c,100);
for(i=0;i<100;i++)
{
if(c[i]=='=')
{

number_star=i;
}
if(c[i]=='#')
{
number_end=i;
}
}
for(i=number_star+1;i<number_end;i++)
{
r[k]=c[i];
k++;
}
number=atof(r);
parameter[dd]=number;
if(number_star==0)
{
dd--;
}
dd++;
}
read.close();
cout<<"parameter[0]="<<parameter[7]<<"\n";
cout<<"parameter[1]="<<parameter[8]<<"\n";
cout<<"parameter[2]="<<parameter[9]<<"\n";
//----------------------------------------------------------------------只讀數字的部分
//------------------------------------------------------------------把讀到的值載入變數
L=parameter[0];
Dx=parameter[1]; //-------one cell's particle numeber
Dt=parameter[2];
total_times=parameter[3];
Nd=parameter[4];
P_distribution=int(parameter[5]);

P_collision=int(parameter[6]);
Collision_nw_gama=parameter[7];
c_arfa=parameter[8];
P_weighting=int(parameter[9]);
re_loop=int(parameter[10]);
KE=int(L/Dx)+3;
AP=int(L*Nd);
NSTEPS=int(total_times/Dt);
WT=NSTEPS+1;
//----------------------------------------------------------------------
//-------------------------------------------plasma space
NT=new double[NSTEPS];
cell_t=new double[NSTEPS];
//-------------------------------------------plasma space
Ex=new double[KE];
potential=new double[KE];
Exh=new double[KE];
q=new double[KE];
N=new double[KE];
ni=new double[KE];
ne=new double[KE];
cell_x=new double[KE];
//-------------------------------------------particle data
pEx=new double[AP];
pEy=new double[AP];
pEz=new double[AP];
vx=new double[AP];
vy=new double[AP];
vz=new double[AP];
v_perpendicular=new double[AP];
total_v=new double[AP];
x=new double[AP];
hlx=new double[AP];
nu=new double[AP];
//-------------------------------------------particle data

//-------------------------------------------datat 分析
totalka=new double[WT];
totalkh=new double[WT];
totalE=new double[WT];
totalEF=new double[WT];
totalP=new double[WT];
FTC=new double[WT];
FTS=new double[WT];
FT=new double[WT];
data_r=new double[WT];
Ha4=new double[WT];
Ha6=new double[WT];
Ha8=new double[WT];
Ha10=new double[WT];
w_time=new double[WT];
log_Ha4=new double[WT];
log_Ha6=new double[WT];
log_Ha8=new double[WT];
log_Ha10=new double[WT];
//------------------------------------------------------------------把讀到的值載入變數
}
inline void first()
{
double deby_lenth,plasma_frequency;
Input_parameter();
n0=1.0e15;
Temperature1=1000;
deby_lenth=sqrt(epaz*tk*Temperature1/(n0*e*e));
plasma_frequency=sqrt(n0*e*e/(epaz*em));
cout<<"deby_lenth="<<deby_lenth<<"\n";
cout<<"plasma_frequency="<<plasma_frequency<<"\n";
a_deby_lenth=1/Dx;/////////////////////////////////////////////////////////////
PPC=double(Nd)/a_deby_lenth;

dx=deby_lenth/a_deby_lenth;
a_plasma_frequency=1/Dt;//////////////////////////////////////////////////////
dt=1.0/(plasma_frequency*a_plasma_frequency);
thermal_velocity=sqrt(tk*Temperature1/em);
k0=2*M_PI/a_wave_cells;
T=0;
cout<<"="<<NSTEPS;
//---------------------------------------------
walt=NSTEPS;
walx=KE-1;
wst=0;
wsx=0;
jumpt=NSTEPS/500+1;
jumpx=1;
cout<<"jumpt="<<jumpt<<" "<<jumpx;
//-------------------------------------------
for(k=0;k<KE;k++)
{
cell_x[k]=k/a_deby_lenth;
}
for(k=0;k<NSTEPS;k++)
{
cell_t[k]=k/a_plasma_frequency;
}
Collision_gama=Collision_nw_gama*plasma_frequency;//一秒鐘碰撞的次數
//------------------------------------------------Monte carlo collision constan
RMK=1.0-exp(-Collision_gama*dt) ;
if(Collision_gama*dt<1.0e-6)

{
RMK=Collision_gama*dt;
}
//------------------------------------------------Monte carlo collision constan
//---------------------------------------------------BL collision 會使用到的參數
BLrdt=Collision_gama*dt;
BLeta_v=sqrt((tk*Temperature1/em)*(1.0-exp(-2.0*BLrdt)));
BLeta_x=sqrt(fabs(dt*dt*tk*Temperature1/(em*BLrdt)*(2.0-(3.0-4.0*exp(-BLrdt)+exp(-2.0*BLrdt))/BLrdt)));
//BLeta_x=dt*dt*tk*Temperature1/(em*BLrdt)*(2.0-(3.0-4.0*exp(-BLrdt)+exp(-2.0*BLrdt))/BLrdt);
// cout<<" aa "<<2*exp(-1.0*BLrdt)<<" aa"<<exp(-2.0*BLrdt)<<" aa";
BLcrv=dt*tk*Temperature1*pow((1.0-exp(-BLrdt)),2.0)/(em*BLeta_v*BLeta_x*BLrdt);
BLc0=exp(-BLrdt);
BLc1=(1.0-BLc0)/BLrdt;
BLc2=(1.0-BLc1)/BLrdt;
BLc3=(0.5-BLc2)/BLrdt;
//---------------------------------------------------BL collision 會使用到的參數
}
