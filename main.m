clc; clear;
clear all;

g=9.81;

%% mesh parameters
m = 21;         %along x and i
n = 21;         %along z and j

mesh_parameters=[m n];
%% geometry - diameter, length, radial clearence
D = 0.04;   %shaft diameter;
R = D/2;
L = 0.04;      %bearing length
h0 = 75e-6;     %radial clearence

nu = 1/h0; % scaling factor

geometry_parameters=[L R h0];

%% operational parameters - RPM, viscosity, load
n_rpm = 2800;    %rotation speed
omega = n_rpm*2*pi/60;
mu = 1.11e-3;     %dynamic viscosity
PresCond=0;     %ambient pressure condition

operational_parameters=[n_rpm, mu, PresCond];
%% force parameters
mass=3.970;         %rotor full mass, kg 3.793 (REAL MASS = 3.970 4.4% different)

Jp=0.00075; %polar moment of inertia, kg m^2  3693 0.00075
Jd=0.047331604; % diametral moment of inertia, kg m^2 0.047331604;

%% geometry of two bearings system
a = 87.548e-3; % distance from centre of mass to bearing#1 87.548e-3;
b = 132.452e-3; % distance from centre of mass to bearing#2 THIS ONE IS ACTIVE 132.452e-3;
c = a + 130e-3;  % distance from centre of mass to the coupling (b#1 direction);
d = 22.452e-3; %distance from centre of mass to the point where fd is applied to 22.452e-3

s = 30e-3; %distance from bearing centre to displacement sensors

F2=a*mass*g/(a+b); % static force on bearing#2
F1=mass*g-F2;      % static force on bearing#1
delta=4e-6;       % unbalance 4e-6
force_parameters=[mass];

%% dynamics calculation parameters
Tcalc=0.15; % seconds
% equilibrium position explicit determination
FdX=0;FdY=0;
[t, x] = ode45 (@(t, x) EOM(t, x,FdX,FdY,F2, mesh_parameters, geometry_parameters, operational_parameters), [0 Tcalc], [0 0 0 0]);
xeq=x(end,1);
yeq=x(end,3);

%% dynamics training set calculation
% [XC1, FC1] = dynamics(F1, mesh_parameters, geometry_parameters, operational_parameters,T_calc);

%% parameters of CS elements
% proximity transducer properties
Vmax=10;  % max and min
Vmin=0;  % output voltage
Hmax=6e-3; % max and min
Hmin=0e-3; % measured distance in the linear range
Kps=(Vmax-Vmin)/(Hmax-Hmin); %proximity sensor voltage2gap coefficient
ps=1e-6; % accuracy
taus=0.00125; %time constant

% band limited white noise
n_pow=1e-10; %n_pow=1e-8;

% ADC resolution and sampling freq
fs=1000; % sampling freq in Hz
Ts=1/fs; % sampling time in sec
adc_acc=0.0063; % datasheet absolute accuracy in V

% Servovalve characteristics
Vmax_sv=10; % max and min
Vmin_sv=0; % input voltage
Ksv = 30; % gain coefficient
tau = 0.12; % time constant
svs = 0.1*Vmax_sv; %voltage sensitivity
%% linearised model

%[FX0,FY0]= CalculateLoadCapacity(mesh_parameters,geometry_parameters, operational_parameters,[xeq 0 yeq 0]);

mm=F2/g;

deltax=1e-6;
deltav=1e-5;

[FXX1,FYY1]= CalculateLoadCapacity(mesh_parameters,geometry_parameters, operational_parameters,[xeq+deltax 0 yeq 0]);
[FXX2,FYY2]= CalculateLoadCapacity(mesh_parameters,geometry_parameters, operational_parameters,[xeq-deltax 0 yeq 0]);

Kxx=(FXX1-FXX2)/(2*deltax);
Kyx=(FYY1-FYY2)/(2*deltax);

[FXX3,FYY3]= CalculateLoadCapacity(mesh_parameters,geometry_parameters, operational_parameters,[xeq 0 yeq+deltax 0]);
[FXX4,FYY4]= CalculateLoadCapacity(mesh_parameters,geometry_parameters, operational_parameters,[xeq 0 yeq-deltax 0]);

Kxy=(FXX3-FXX4)/(2*deltax);
Kyy=(FYY3-FYY4)/(2*deltax);

[FXX1,FYY1]= CalculateLoadCapacity(mesh_parameters,geometry_parameters, operational_parameters,[xeq 0+deltav yeq 0]);
[FXX2,FYY2]= CalculateLoadCapacity(mesh_parameters,geometry_parameters, operational_parameters,[xeq 0-deltav yeq 0]);

Bxx=(FXX1-FXX2)/(2*deltav);
Byx=(FYY1-FYY2)/(2*deltav);

[FXX3,FYY3]= CalculateLoadCapacity(mesh_parameters,geometry_parameters, operational_parameters,[xeq 0 yeq 0+deltav]);
[FXX4,FYY4]= CalculateLoadCapacity(mesh_parameters,geometry_parameters, operational_parameters,[xeq 0 yeq 0-deltav]);

Bxy=(FXX3-FXX4)/(2*deltav);
Byy=(FYY3-FYY4)/(2*deltav);

A=[0 1 0 0; Kxx/mm Bxx/mm Kxy/mm Bxy/mm; 0 0 0 1; Kyx/mm Byx/mm Kyy/mm Byy/mm];

B=[0 0; Ksv/mm 0; 0 0; 0 Ksv/mm];

C=[1 0 0 0; 0 0 0 0; 0 0 1 0; 0 0 0 0];

D=0;

b2lin=ss(A,B,C,D);
