%% PID implementation
clear;clc;close all

% Model Set: Car and Rob
%%% x_M = x; x_m = x + l*sin(q)
%%% 【Newton Formula】 
%%% (M+m)*ddx + m*L*cos(q)*ddq - m*L*sin(q)*dq^2 = F - B*dx 
%%% 【Euler Formula】 
%%% (I + mL^2)*ddq = mgLsin(q) - mLcos(q)*ddx - b*dq
M = 1;
m = 0.1;
B = 0.1;
b = 0.1;
I = 1;
L = 1; % 质心到铰点距离
g = 9.8;
%%% dX = AX+BU
%%% U is the force applied on the car
%%% X = [x,dx,theta,dtheta]
p = I*(M+m) + M*m*L^2;
A = [0,1,0,0;
     0,-(I+m*L^2)*B/p,m^2*g*L^2/p,0;
     0,0,0,1;
     0,-m*L*b*p,m*g*L*(M+m)/p,0];
B = [0;(I+m*L^2)/p;0;m*L/p];

% PID Set
Kp = 100;
Ki = 10;
Kd = 10;
% Control Set
err = [];
dt = 0.1;
tspan = 0:dt:10;
theta_ref = 0;
theta0 = 5*pi/180;
dtheta0 = 0;
x0 = 0;
dx0 = 0;
int_e0 = 0;

% solve
z0 = [x0;dx0;theta0;dtheta0;int_e0];
odefun = @(t,z)Inverse_Pendulum_PID_ODE(t,z,A,B,Kp,Ki,Kd,theta_ref);
[t,z] = ode45(odefun,tspan,z0);

% visualization
x = z(:, 1);
dx = z(:, 2);
theta = z(:, 3);
dtheta = z(:, 4);
int_e = z(:, 5);

figure;
plot(t, theta*180/pi, 'LineWidth', 1.5);
grid on;
xlabel('Time / s');
ylabel('\theta / deg');
title('Inverted Pendulum Angle Response under PID Control');

figure;
plot(t, x, 'LineWidth', 1.5);
grid on;
xlabel('Time / s');
ylabel('Cart Position x / m');
title('Cart Position Response');

% figure;
% plot(t, u, 'LineWidth', 1.5);
% grid on;
% xlabel('Time / s');
% ylabel('Control Force u / N');
% title('PID Control Input');

%% Function for PID
function dz = Inverse_Pendulum_PID_ODE(~,z,A,B,Kp,Ki,Kd,theta_ref)
    X = z(1:4);
    theta = z(3);
    dtheta =z(4);
    int_e = z(5);

    e = theta_ref - theta;
    de = -dtheta;
    u = Kp*e + Ki*int_e + Kd*de;

    dX = A*X+B*u;

    dini_e = e;
    dz = [dX;dini_e];
end