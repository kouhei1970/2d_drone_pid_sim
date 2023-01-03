/******************************
 *                            
 * ルンゲクッタ法による 
 * モータシミュレーション                   
 *                            
 ******************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h> //可変長引数関数を使うために必要
#include <math.h>
#include <time.h>
#include "pid.hpp"

#define FRONT 0
#define REAR 1

//定数ノミナル
const double Rm = 1.16;//Resistance[Ohm]
const double Km = 8.73e-4;//Torque constant[Nm/A]
const double Jr = 2.75e-8;//Moment of inertia[kg m^2]
const double Jp = 3.50e-8;//Moment of inertia[kg m^2]
const double Jm = 2*Jr + Jp;
const double Cq = 4.49713785e-10;//Cofficient of torque (Propeller)
const double Dm = 1.02432352e-07;   //Cofficient of viscous damping [Nm s]
const double End_time =2.0;//Time [s]
long cpu_time;

//モータ状態構造体
typedef struct 
{
  double omega;
  double u;
  double omega_;
  double u_;
} motor_t;

//ドローン構造体
typedef struct 
{
  double q;
  double theta;
  double q_;
  double theta_;
} drone_t;

//Motor Equation of motion
//TL = Cq omega^2
//Jm domega/dt + (Dm + K^2/Rm) omega + TL = Km u/R
//omega:angular velocity
//u:value[0]=u
//t:time
double omega_dot(double omega, double t, double *value)
{
  double u =value[0];
  double TL=Cq * omega * omega;
  return (Km*u/Rm - (2*Dm + Km*Km/Rm)*omega - TL)/Jm;
}

//Drone equation of motion
// Iy dq/dt = armr 2 Ct(front_omega^2 - rear_omega^2)
//value[0]:front_omega
//value[1]:rear_omega
double qdot(double q, double t, double *value)
{
  double Iy = 2.19e-5;
  double armr=0.033;
  double Ct=0.57e-7;
  double front_omega = value[0];
  double rear_omega  = value[1];

  return 2*armr*Ct*(front_omega*front_omega - rear_omega*rear_omega)/Iy;
}

double thetadot(double theta, double t, double *value)
{
  double q = value[0];

  return q;
}



//Runge Kutta method
//dxdy:derivative
//x:state
//t:time
//h:step size
//n:number of arguments
double rk4(double (*dxdt)(double, double, double*), double x, double t, double h, int n, ...)
{
  va_list args;
  double *value;
  double k1,k2,k3,k4;

  value=(double*)malloc(sizeof(double) * n);
  va_start(args , n);
  for(int i=0;i<n;i++)
  {
    value[i]=va_arg(args, double);
  }
  va_end(args);
  
  k1 = h * dxdt(x, t, value);
  k2 = h * dxdt(x+0.5*h*k1, t+0.5*h, value);
  k3 = h * dxdt(x+0.5*h*k2, t+0.5*h, value);
  k4 = h * dxdt(x+h*k3, t+h, value);

  free(value);
  
  return x+(k1 + k2*2.0 + k3*2.0 + k4)/6;
}

void save_state(double t, motor_t* front_motor, motor_t* rear_motor, drone_t* drone)
{
  front_motor->omega_ = front_motor->omega;
  front_motor->u_ = front_motor->u;
  rear_motor->omega_ = rear_motor->omega;
  rear_motor->u_ = rear_motor->u;
  drone->q_ = drone->q;
  drone->theta_ = drone->theta;
}

void print_state(double t, motor_t* front_motor, motor_t* rear_motor, drone_t* drone)
{
  printf("%11.8f %11.8f %11.8f %11.8f %11.8f\n",
    t, 
    front_motor->omega,
    rear_motor->omega,
    drone->q,
    drone->theta
  );
}

//
//Motor simulator main
//
void drone_sim(void)
{
  motor_t motor[2];
  drone_t drone;
  

  //state init
  motor[FRONT].omega = 0;
  motor[REAR].omega = 0;
  motor[FRONT].u = 0;
  motor[REAR].u = 0;
  drone.q=0.0;
  drone.theta = 1.0*M_PI/180.0;
  
  double t       = 0.0; //time
  double step    = 0.0001;//step size

  //initial state output
  print_state(t, &motor[FRONT], &motor[REAR], &drone);
  
  //Simulation loop
  cpu_time = clock();
  while(t < End_time )    
  {
    //Control
    double kp1=0;
    double ti1=1000;
    double td1=0;
    double kp2=0;
    double ti3=1000;
    double td4=0;

    //Save state
    save_state(t, &motor[FRONT], &motor[REAR], &drone);

    //Update(Runge-Kutta method)
    motor[FRONT].omega = rk4(omega_dot, motor[FRONT].omega_, t, step, 1, motor[FRONT].u_);
    motor[REAR].omega  = rk4(omega_dot, motor[REAR].omega_,  t, step, 1, motor[REAR].u_ );
    drone.q = rk4(qdot, drone.q_, t, step, 2, motor[FRONT].omega_, motor[REAR].omega_);
    drone.theta = rk4(thetadot, drone.q_, t, step, 1, drone.q_); 

    t = t + step;
    
    //Output
    print_state(t, &motor[FRONT], &motor[REAR], &drone);
  }
  cpu_time = clock() -cpu_time;
  printf("#elapsed time %f\n", (double)cpu_time/CLOCKS_PER_SEC);
}

int main(void)
{
  drone_sim();
}
