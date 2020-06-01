function [T,P,h,rho,a] = StandardAtmosphereV2(Tin,Pin,hin)%,xxx) 

%% READ ME

% This is a metric Standard Atmosphere 
% calculator for 
% an altitude range of 0 to 120 km
% or 0 to 393,701 ft.
%
% Inputs are sea level values of 
% Tin = 288.15 K, Pin = 101325 Pa,
% hin = 0 m, rhoin = 1.225l kg/m^3, these 
% can also be altered depending on 
% 'on the day' conditions.
%
% This script contains formulas from
% 'ATMOSPHERIC and SPACE FLIGHT DYNAMICS'
% by Ashish Tewari, ISBN-10: 0-8176-4437-7,
% Birka¨user Boston and NASA's U.S. 1979 
% Standard Atmosphere formulas,
% http://www.braeunig.us/space/atmmodel.htm.
%
% This script is only intended for civilian
% purposes.
%
% Code by: Reza Aliakbari.

%% Standard Atmosphere

% Initial Conditions

T(1) = Tin;       % [K]
P(1) = Pin;       % [Pa]
h(1) = hin;       % [m]
M = 0.0289644;    % [kg/mol]
R = 8.3144598;    % [J/(mol.K)] (constant up till 100 km)
R_kg = 287;       % [J/(kg.K)] (constant up till 100 km)  
g = 9.81;         % [m/s^2]
beta = 2/6371000; % [1/m]
gamma = 1.4;     

rho(1) = P(1)/(R_kg*T(1)); % [kg/m^3]
a(1) = sqrt(gamma*R_kg*T(1)); % [m/s] Speed of Sound

% Loop Conditions

h_scalar = 120; % [m]
dh = 0.01;      % [m]

% Loop

for i = 1:(h_scalar/dh)
    
    if h(i) >= 0 && h(i) < 11 % Troposphere
        
        T(i+1) = T(i) - 6.5*dh;
        P(i+1) = P(1)*(T(1)/(T(i) - 6.5*dh))^(34.1632/-6.5);
        rho(i+1) = P(i)/(R_kg*T(i));
        a(i+1) = sqrt(gamma*R_kg*T(i));
        
    elseif h(i) >= 11 && h(i) < 20 % Tropopause
        
        T(i+1) = T(i);
        P(i+1) = P((11/dh))*exp(-34.1632*(h(i) - 11)/T((11/dh)));
        rho(i+1) = P(i)/(R_kg*T(i));
        a(i+1) = sqrt(gamma*R_kg*T(i));
        
    elseif h(i) >= 20 && h(i) < 32 % Stratosphere
       
        T(i+1) = T(i) + dh;
        P(i+1) = P((20/dh))*(T((20/dh))/(T((20/dh)) + (h(i) - 20)))^(34.1632);
        rho(i+1) = P(i)/(R_kg*T(i));
        a(i+1) = sqrt(gamma*R_kg*T(i));
        
    elseif h(i) >= 32 && h(i) < 47 % Stratosphere
        
        T(i+1) = T(i) + 2.8*dh;
        P(i+1) = P((32/dh))*(T((32/dh))/(T(i) + 2.8*(h(i) - 32)))^(34.1632/2.8);
        rho(i+1) = P(i)/(R_kg*T(i));
        a(i+1) = sqrt(gamma*R_kg*T(i));
        
    elseif h(i) >= 47 && h(i) < 51 % Stratopause
        
        T(i+1) = T(i);
        P(i+1) = P((47/dh))*exp(-34.1632*(h(i) - 47)/270.65);
        rho(i+1) = P(i)/(R_kg*T(i));
        a(i+1) = sqrt(gamma*R_kg*T(i));
        
    elseif h(i) > 51 && h(i) <= 71 % Mesosphere
        
        T(i+1) = T(i) - 2.8*dh;
        P(i+1) = P((51/dh))*(T((51/dh))/(T((51/dh)) - 2.8*(h(i) - 51)))^(34.1632/-2.8);
        rho(i+1) = P(i)/(R_kg*T(i));
        a(i+1) = sqrt(gamma*R_kg*T(i));
        
    elseif h(i) > 71 && h(i) <= 84.86 % Mesosphere
     
        T(i+1) = T(i) - 2*dh;
        P(i+1) = P((71/dh))*(T((71/dh))/(T((71/dh)) - 2*(h(i) - 71)))^(34.1632/-2);
        rho(i+1) = P(i)/(R_kg*T(i));
        a(i+1) = sqrt(gamma*R_kg*T(i));
        
    elseif h(i) > 84.86 && h(i) <= 91 % Mesopause
        
        T(i+1) = T(i);
        P(i+1) = exp(-0.0000000422012*h(i)^5 + 0.0000213489*h(i)^4 - 0.00426388*h(i)^3 + 0.421404*h(i)^2 - 20.8270*h(i) + 416.225);
        rho(i+1) = P(i)/(R_kg*T(i));
        a(i+1) = sqrt(gamma*R_kg*T(i));
        
    elseif h(i) > 91 && h(i) <= 110 % Thermosphere
        
        T(i+1) = 263.1905 - (263.1905 - T((91/dh)))*sqrt(1 - ((h(i) - 91) / -19.9429)^2);
        P(i+1) = exp(-0.0000000422012*h(i)^5 + 0.0000213489*h(i)^4 - 0.00426388*h(i)^3 + 0.421404*h(i)^2 - 20.8270*h(i) + 416.225);
        rho(i+1) = exp( 0.000000075691*h(i)^5 - 0.0000376113*h(i)^4 + 0.0074765*h(i)^3 - 0.743012*h(i)^2 + 36.7280*h(i) - 729.346);
        a(i+1) = sqrt(gamma*R_kg*T(i));
        
    elseif h(i) > 110 && h(i) <= 120 % Thermosphere
        
        T(i+1) = T(i) + 12*dh;
        P(i+1) = exp(-0.0000000422012*h(i)^5 + 0.0000213489*h(i)^4 - 0.00426388*h(i)^3 + 0.421404*h(i)^2 - 20.8270*h(i) + 416.225);
        rho(i+1) = exp( 0.000000075691*h(i)^5 - 0.0000376113*h(i)^4 + 0.0074765*h(i)^3 - 0.743012*h(i)^2 + 36.7280*h(i) - 729.346);
        a(i+1) = sqrt(gamma*R_kg*T(i));
        
    else
        
        fprintf('Error!')
        break
          
    end
    
    h(i+1) = h(i) + dh;
    
end

h_ft = h.*3.28084; % [ft]

% Plots

figure;%(xxx);

subplot(2,3,1); hold on; plot(T,h);
xlabel('Temperature (K)');
ylabel('Altitude (km)');
ylim([0 h_scalar]); grid on;

subplot(2,3,2); hold on; plot(P*10^-3,h);
xlabel('Air Pressure (kPa)');
ylabel('Altitude (km)');
ylim([0 h_scalar]); grid on;

subplot(2,3,3); hold on; plot(rho,h);
xlabel('Air Density (kg/m^3)');
ylabel('Altitude (km)');
ylim([0 h_scalar]); grid on;

subplot(2,3,4); hold on; plot(T,h_ft);
xlabel('Temperature (K)');
ylabel('Altitude (kft)');
ylim([0 h_scalar*3.28084]); grid on;

subplot(2,3,5); hold on; plot(P*10^-3,h_ft);
xlabel('Air Pressure (kPa)');
ylabel('Altitude (kft)');
ylim([0 h_scalar*3.28084]); grid on;

subplot(2,3,6); hold on; plot(rho,h_ft);
xlabel('Air Density (kg/m^3)');
ylabel('Altitude (kft)');
ylim([0 h_scalar*3.28084]); grid on;

sgtitle('Temperature, Air Pressure, and Air Density with Altitude');

figure;%(xxx+1);

subplot(1,2,1); hold on; plot(a,h)
xlabel('Speed of Sound (m/s)');
ylabel('Altitude (km)');
ylim([0 h_scalar]); grid on;

subplot(1,2,2); hold on; plot(a,h_ft)
xlabel('Speed of Sound (m/s)');
ylabel('Altitude (kft)');
ylim([0 h_scalar*3.28084]); grid on;

sgtitle('Speed of Sound with Altitude');

end
