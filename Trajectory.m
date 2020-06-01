close all; clear; clc;

%% READ ME

% This script contains formulas from 'ROCKET
% PROPULSION', by Stephen Heister, William
% E. Anderson, Timothee L. Porpoint, R. Joeseph 
% Cassidy,'ATMOSPHERIC and SPACE FLIGHT DYNAMICS'
% by Ashish Tewari, ISBN-10: 0-8176-4437-7,
% Birka¨user Boston and NASA's U.S. 1979 
% Standard Atmosphere formulas,
% http://www.braeunig.us/space/atmmodel.htm.
%
% This script is only intended for civilian
% purposes.
%
% Code by: Preetham Akula (preetham.akula@rmit.edu.au) & Reza Aliakbari
% (reza.aliakbari@rmit.edu.au)
%% Standard Atmosphere
%for T_0 = 288:10:318
T_0 = 288.15;  % [K]
P_0 = 101325;  % [Pa]
h_0 = 0;       % [m]

[T,P,h,rho,a] = StandardAtmosphereV2(T_0,P_0,h_0);%,3);

%% Altitude / Air Density Relation

% Fourier 5

Alt_AirDen = rho_h_CurveFit(h, rho);%,2);
F5 = coeffvalues(Alt_AirDen);

% Altitude /  Air Density Expression

syms h
f_rho_alt = F5(1)...
    + F5(2)*cos(1*h*F5(12)) + F5(3)*sin(1*h*F5(12))...
    + F5(4)*cos(2*h*F5(12)) + F5(5)*sin(2*h*F5(12))...
    + F5(6)*cos(3*h*F5(12)) + F5(7)*sin(3*h*F5(12))...
    + F5(8)*cos(4*h*F5(12)) + F5(9)*sin(4*h*F5(12))...
    + F5(10)*cos(5*h*F5(12)) + F5(11)*sin(5*h*F5(12));

%% Thrust 

[thrust_data] = xlsread('O3400');
thrust_data_time = thrust_data(:,1);
thrust_data_force = thrust_data(:,2);

% Fourier 6

Thrust_Time = Thrust_Curve(thrust_data_time,thrust_data_force);
F6 = coeffvalues(Thrust_Time);

%% 2D & 3D Trajectory

% time

tt = 100;   % [s]
dt = 0.001; % [s]
t(1) = 0;   % [s]

% velocity

v(1) = 0;      % [m/s]
v_y(1) = 0;
v_r(1) = 0;    % [m/s]
v_w = 0/3.6;   % [m/s] with [km/h] input
v_w_y = 0/3.6;
w_dir = -1;    % wind direction
w_dir_y = -1;
r_e = 6378100;
theta(1) = 0;

% angle

phi(1) = 89.999*(pi/180);   % [rad] flight path angle from x-axis
psi(1) = phi(1);            % [rad] force direction from x-axis
alpha(1) = 89.999*(pi/180);

% displacement

x(1) = 0;    % [m]
x_actual(1) = 0;
z(1) = 0;    % [m]
y(1) = 0;
Dist(1) = 0; % [m]

% force

D(1) = 0; % [N]
L(1) = 0; % [N]

% propulsion
f_F_t = @(Tt) F6(1)...
    + F6(2)*cos(1*Tt*F6(14)) + F6(3)*sin(1*Tt*F6(14))...
    + F6(4)*cos(2*Tt*F6(14)) + F6(5)*sin(2*Tt*F6(14))...
    + F6(6)*cos(3*Tt*F6(14)) + F6(7)*sin(3*Tt*F6(14))...
    + F6(8)*cos(4*Tt*F6(14)) + F6(9)*sin(4*Tt*F6(14))...
    + F6(10)*cos(5*Tt*F6(14)) + F6(11)*sin(5*Tt*F6(14))...
    + F6(12)*cos(6*Tt*F6(14)) + F6(13)*sin(6*Tt*F6(14));
F(1) = f_F_t(0.001);
Isp = 196.5; % [s]
t_b = 6.16;  % [s]

% acceleration

g_e = 9.81;     % [rad/s]
dv_dt(1) = 0;   % [rad/s]
dphi_dt(1) = 0; % [rad/s]

% mass

m_0 = 34.675; % [kg]
m_dot = @(Tt) abs((1/(g_e*Isp))*(F6(1)...
    + F6(2)*cos(1*Tt*F6(14)) + F6(3)*sin(1*Tt*F6(14))...
    + F6(4)*cos(2*Tt*F6(14)) + F6(5)*sin(2*Tt*F6(14))...
    + F6(6)*cos(3*Tt*F6(14)) + F6(7)*sin(3*Tt*F6(14))...
    + F6(8)*cos(4*Tt*F6(14)) + F6(9)*sin(4*Tt*F6(14))...
    + F6(10)*cos(5*Tt*F6(14)) + F6(11)*sin(5*Tt*F6(14))...
    + F6(12)*cos(6*Tt*F6(14)) + F6(13)*sin(6*Tt*F6(14))));
m_p = integral(m_dot,0,t_b); % [kg]
m(1) = m_0;                  % [kg]

% aerodynamics 

number_fin = 3;
t_fin = 0.006;
root_tip_fin = 0.152;
d_ref_airf = 0.131;
length = 3.04;
A_ref_airf = pi*(d_ref_airf/2)^2; 
A_ref_fin = (t_fin*root_tip_fin)*number_fin;
A_ref = A_ref_airf + A_ref_fin;
A_ref_long = length*d_ref_airf;
A_body = pi*d_ref_airf*length;
Rho(1) = 1.225; 
Cd_long = 0.27;
Cd = 0.2;
Cl = 0;

% wind drag

D_w(1) = 0.5*Rho(1)*(v_w^2)*A_ref_long*Cd_long;

D_w_y(1) = 0.5*Rho(1)*(v_w_y^2)*A_ref_long*Cd_long;

% dynamic pressure

P_d(1) = (1/2)*Rho(1)*v(1)^2;

% flow regimes

mew = 1.802E-5;
Re_x(1) = 0;

% #1 - Blasius Solution (Laminar)

Cf_L(1) = 0;

% #2 - Prandtl's One-Seventh-Power Law (Turbulent)

Cf_T(1) = 0;

% skin friction

D_skfr_L(1) = 0;
D_skfr_T(1) = 0;
D_skfr(1) = 0;

for i = 1:(tt/dt)
    
    if phi(i) == 90
        phi(i) = 89.999999;
        psi(i)= 89.999999;
    end
   
    dv_dt(i) = (w_dir*D_w(i)/m(i))*cos(phi(i)) + ((F(i)-D(i)-D_skfr(i))*cos(psi(i)-phi(i))/m(i)...
        -(L(i)*sin(psi(i)-phi(i)))/m(i)...
        -g_e*sin(phi(i)));
    
    dv_dt_y(i) = (w_dir_y*D_w(i)/m(i))*cos(alpha(i));
    
    if t(i) == 0
        dphi_dt(i) = 0;
        
        dalpha_dt(i) = 0;
    else
        dphi_dt(i) = (w_dir*D_w(i)/(m(i)*v(i)))*sin(phi(i)) + ((F(i)-D(i)-D_skfr(i))*sin(psi(i)-phi(i))/(m(i)*v(i))...
            +(L(i)*cos(psi(i)-phi(i)))/(m(i)*v(i))...
            -g_e*cos(phi(i)))/(v(i));
        
        dalpha_dt(i) = (w_dir_y*D_w(i)/(m(i)*v(i)))*sin(alpha(i));
    end
    
    t(i+1) = t(i) + dt;
    
    if t(i) <= t_b 
        F(i+1) = F6(1)...
    + F6(2)*cos(1*t(i+1)*F6(14)) + F6(3)*sin(1*t(i+1)*F6(14))...
    + F6(4)*cos(2*t(i+1)*F6(14)) + F6(5)*sin(2*t(i+1)*F6(14))...
    + F6(6)*cos(3*t(i+1)*F6(14)) + F6(7)*sin(3*t(i+1)*F6(14))...
    + F6(8)*cos(4*t(i+1)*F6(14)) + F6(9)*sin(4*t(i+1)*F6(14))...
    + F6(10)*cos(5*t(i+1)*F6(14)) + F6(11)*sin(5*t(i+1)*F6(14))...
    + F6(12)*cos(6*t(i+1)*F6(14)) + F6(13)*sin(6*t(i+1)*F6(14));
        m(i+1) = m_0 - integral(m_dot,0,t(i+1));   
    else
        F(i+1) = 0;
        m(i+1) = (m_0 - m_p);
    end
    
    % * equations

    v_star = v(i)+dv_dt(i)*dt;
    v_y_star = v_y(i)+dv_dt_y(i)*dt;
    F_star = F(i+1);
    D_star = 0.5*Rho(i)*(v_r(i)^2)*A_ref*Cd;
    D_w_star = 0.5*Rho(i)*(v_w^2)*A_ref_long*Cd_long;
    D_w_y_star = 0.5*Rho(i)*(v_w_y^2)*A_ref_long*Cd_long;
    L_star = 0.5*Rho(i)*(v_r(i)^2)*A_ref*Cl;
    phi_star = phi(i)+dphi_dt(i)*dt;
    psi_star = atan((v_star*sin(phi_star))/(v_star*cos(phi_star)));
    m_star = m(i+1);
    
    dv_dt_star = (w_dir*D_w_star./m_star)*cos(phi(i)) + ((F_star-D_star-D_skfr(i))*cos(psi_star-phi_star)-L_star*sin(psi_star-phi_star)...
        -m_star*g_e*sin(phi_star))./m_star;
    dv_dt_y_star = (w_dir_y*D_w_y_star/m_star)*sin(alpha(i));
    dphi_dt_star = (w_dir*D_w_star./(m_star*v_star))*sin(phi(i)) + ((F_star-D_star-D_skfr(i))*sin(psi_star-phi_star)+L_star*cos(psi_star-phi_star)...
        -m_star*g_e*cos(phi_star))./(m_star*v_star);
    dalpha_dt_star = (w_dir_y*D_w_y_star/(m_star*v_star))*sin(alpha(i));
    
    v(i+1) = v(i)+(dt/2)*(dv_dt(i)+dv_dt_star);
    v_y(i+1) = v_y(i)+(dt/2)*(dv_dt_y(i)+dv_dt_y_star);
    phi(i+1) = phi(i)+(dt/2)*(dphi_dt(i)+dphi_dt_star);
    alpha(i+1) = alpha(i)+(dt/2)*(dalpha_dt(i)+dalpha_dt_star);
    
    if t(i) <= t_b
        psi(i+1) = phi(i+1);
    else
        psi(i+1) = atan((v(i+1)*sin(phi(i+1)))/(v(i+1)*cos(phi(i+1))));
    end
    
    v_r(i+1) = (v(i+1)*sin(phi(i+1)))/sin(psi(i+1));
    
    x(i+1) = x(i)+(dt/2)*(v(i)*cos(phi(i))+v(i+1)*cos(phi(i+1)));% + v_w*dt;
    
    theta(i+1) = 2*atan(x(i+1)/(2*r_e));
    
    x_actual(i+1) = theta(i+1)*r_e; 
    
    y(i+1) = y(i)+(dt/2)*(v_y(i)*cos(alpha(i))+v_y(i+1)*cos(alpha(i+1)));% + v_w*dt;
    
    if z(i) >= 0
        z(i+1) = z(i)+(dt/2)*(v(i)*sin(phi(i))+v(i+1)*sin(phi(i+1)));
    elseif z(i) < 0 && t(i+1) >= t_b
        z(i+1) = 0;
    else
        z(i+1) = 0;
        fprintf('*ERROR* Rapid Unscheduled Disassembly *ERROR*\n\n')
    end
    
    Dist(i+1) = Dist(i)+sqrt((x(i+1)-x(i))^2+(z(i+1)-z(i))^2);
    
    Rho(i+1) = F5(1)...
    + F5(2)*cos(1*(z(i)/1000)*F5(12)) + F5(3)*sin(1*(z(i)/1000)*F5(12))...
    + F5(4)*cos(2*(z(i)/1000)*F5(12)) + F5(5)*sin(2*(z(i)/1000)*F5(12))...
    + F5(6)*cos(3*(z(i)/1000)*F5(12)) + F5(7)*sin(3*(z(i)/1000)*F5(12))...
    + F5(8)*cos(4*(z(i)/1000)*F5(12)) + F5(9)*sin(4*(z(i)/1000)*F5(12))...
    + F5(10)*cos(5*(z(i)/1000)*F5(12)) + F5(11)*sin(5*(z(i)/1000)*F5(12));
    
    % parachute deployment

    if z(i+1) <= (0/3.28) && phi(i+1) < 0
        Cd = 0.8;
        Cl = 0;
        d_parachute = 1;
        A_ref = (1/4)*pi*d_parachute^2;
        %phi(i) = -89.999*(pi/180);
        psi(i+1) = phi(i);
    end
       
    D(i+1) = 0.5*Rho(i)*(v_r(i+1)^2)*A_ref*Cd;
    D_w(i+1) = 0.5*Rho(i)*(v_w^2)*A_ref_long*Cd_long;
    L(i+1) = 0.5*Rho(i)*(v_r(i+1)^2)*A_ref*Cl;
    
    % dynamic pressure
    
    P_d(i+1) = (1/2)*Rho(i)*v_r(i)^2;
    
    % flow regimes and skin friction
    
    Re_x(i+1) = (Rho(i+1)*v_r(i+1)*length)/mew;
    
    % #1 - Blasius Solution (Laminar)
    
    Cf_L(i+1) = 0.664/sqrt(Re_x(i+1));
    
    % #2 - Prandtl's One-Seventh-Power Law (Turbulent)
    
    Cf_T(i+1) = 0.027/(Re_x(i+1))^(1/7);
    
    % skin friction
    
    D_skfr_L(i+1) = (Cf_L(i+1)*Rho(i+1)*v_r(i+1)^2*A_body)/2;
    D_skfr_T(i+1) = (Cf_T(i+1)*Rho(i+1)*v_r(i+1)^2*A_body)/2;
    
    if Re_x(i+1) <= 4000
        D_skfr(i+1) = D_skfr_L(i+1);
    else
        D_skfr(i+1) = D_skfr_T(i+1);
    end
     
    if z(i+1) < 0
        break
    end
    
end

% % ballistics
% 
% v_bo = max(v);
% phi_bo = max(phi);
% r_bo = max(z) + r_e;
% mew_e = 3.986E14;
% v_c = sqrt(mew_e/r_e);
% 
% Q_bo = (v_bo/v_c)^2;
% beta = 2*acos((1-Q_bo*(cos(phi_bo))^2)/(sqrt(1+Q_bo*(Q_bo-2)*(cos(phi_bo))^2)));
% range = r_e*beta
% alt_max = (r_bo/(2-Q_bo))*(2+Q_bo*(Q_bo-2)*(cos(phi_bo))^2) - r_e

figure(6);

subplot(4,2,1);
hold on;
plot(t,v);
xlabel('Time (s)');
ylabel('Velocity (m/s)');

subplot(4,2,2);
hold on;
plot(t,phi*180/pi);
xlabel('Time (s)');
ylabel('Phi (deg)');

subplot(4,2,3);
hold on;
plot(t,psi*180/pi);
xlabel('Time (s)');
ylabel('Psi (deg)');

subplot(4,2,4);
hold on;
plot(x,z);
xlabel('Range (m)');
ylabel('Height (m)');

subplot(4,2,5);
hold on;
plot(t,z);
xlabel('Time (s)');
ylabel('Height (m)');

subplot(4,2,6);
hold on;
plot(t,D); hold on;
plot(t,D_w); hold on;
xlabel('Time (s)');
ylabel('Drag (N)');

subplot(4,2,7);
hold on;
plot(t,F);
xlabel('Time (s)');
ylabel('Thrust (N)');
xlim([0 t_b])

subplot(4,2,8);
hold on;
plot(t,m);
xlabel('Time (s)');
ylabel('Mass (kg)');
xlim([0 t_b])

sgtitle('Flight Parameters');

figure(7)
xlabel('Range (m)');
ylabel('Altitude (m)');
hold on
xlim([0 max(x)])
ylim([0 max(z)])
comet(x,z)
hold off

figure(8)
plot3(x,y,z)
xlabel('Range (m)');
ylabel('Deviation (m)');
zlabel('Altitude (m)');
grid on
hold on

figure(9)
comet3(x,y,z)
xlabel('Range (m)');
ylabel('Deviation (m)');
zlabel('Altitude (m)');
grid on
hold on
comet3(x,y,z)

figure(10)
comet(x_actual,z)

display(max(z));
display(max(z)*3.28);
