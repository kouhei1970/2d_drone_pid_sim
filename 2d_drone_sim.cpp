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

void save_state(motor_t* motor)
{
  motor->omega_ = motor->omega;
  motor->u_ = motor->u;
}

void print_state(double t, motor_t* motor)
{
  printf("%11.8f %11.8f %11.8f\n",
    t, 
    motor->u,
    motor->omega
  );
}

//
//Motor simulator main
//
void motor_sim(void)
{
  motor_t motor;
  
  //state init
  motor.omega = 0;
  motor.u = 3.5;
  
  double t       = 0.0; //time
  double step    = 0.0001;//step size

  //initial state output
  print_state(t, &motor);
  
  //Simulation loop
  cpu_time = clock();
  while(t < End_time )    
  {
    //Save state
    save_state(&motor);

    //Update(Runge-Kutta method)
    motor.omega = rk4(omega_dot, motor.omega_, t, step, 1, motor.u_);
    t = t + step;
    
    //Output
    print_state(t, &motor);
  }
  cpu_time = clock() -cpu_time;
  printf("#elapsed time %f\n", (double)cpu_time/CLOCKS_PER_SEC);
}

int main(void)
{
  motor_sim();
}
