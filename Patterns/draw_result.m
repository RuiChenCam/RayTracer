%This script tries to draw the result from the Starlab test result and
%compare with the simulation
close all

%% Draw 3D Total Gain
%Draw starlab result
data=readtable('MT242021.txt');
phi=data.Phi/2/pi*360;
theta=data.Theta/2/pi*360;
gain=data.Gain_DB;
figure()
title('Total Gain in dB - Starlab')
patternCustom(gain, theta, phi);

% %Draw simulation result
% simu=readtable('simulation.txt');
% theta_simu=simu.Theta;
% phi_simu=simu.Phi;%transform the simulation coordinate to starlab one
% gain_simu=simu.Gain_Total_;
% gain_simu_lin=10.^(gain_simu/10);
% figure()
% title('Total Gain in dB - Simulation')
% patternCustom(gain_simu, theta_simu, phi_simu);


%% Draw polar plot for phi=0 degree
%Draw the starlab result

figure(3)
%index=find(phi==90.000000000292400);
%polarplot(deg2rad(theta(index)),gain(index));
%set(gca,'fontsize', 18);
patternCustom(gain, theta, phi,'CoordinateSystem','polar','Slice','phi','SliceValue',90.000000000292400);
legend('Starlab \phi=90')
%Draw simulation Result
figure(3)
hold on

patternCustom(gain_simu, theta_simu, phi_simu,'CoordinateSystem','polar','Slice','phi','SliceValue',90);
legend('Simulation \phi=90')

%% Draw Axial Ratio
%axial=readtable('Axial_ratio_dB_against_theta_phi.txt');
axial_ratio=data.AxialRatio_dB_;
figure(4)

title('Axial Ratio in dB - Starlab')
patternCustom(axial_ratio, theta, phi,'CoordinateSystem','rectangular','Slice','phi','SliceValue',90.000000000292400);
AR_simu=readtable('AR_simu.txt');
AR_simu=AR_simu.FarField1_dB_;
AR_simu=[AR_simu(38:end);AR_simu(1:37)];
figure(4)
hold on
plot(-180:5:180,AR_simu)
legend("Measured","Simulation")
title('Axial Ratio')
xlabel('\theta - degree')
ylabel('Magnitude - dB')


%% Draw Gain and Efficiency against frequency
% gain_freq=readtable('Gain_dB_against_freq.txt');
% eff_freq=readtable('Efficiency_dB_against_freq.txt');
% figure()
% 
% plot(gain_freq.Frequency/1e6,gain_freq.Gain_DB);
% hold on
% plot(eff_freq.Frequency/1e6,eff_freq.Efficiency_DB);
% xlabel('Frequency: MHz')
% ylabel('dB')
% title('Antenna Gain and Efficiency Against Frequency')
% legend('Gain','Efficiency')


%% Attempt to draw Gx, Gy, Gz
% gain_linear=10.^(gain/10);%Convert to linear
% 
% Gx = (gain_linear.*cos(phi/360*2*pi).*sin(theta/360*2*pi));
% 
% Gy = (gain_linear.*sin(theta/360*2*pi).*sin(phi/360*2*pi));
% 
% Gz = (gain_linear.*cos(theta/360*2*pi));
% 
% figure()
% patternCustom(Gx, theta, phi);
% figure()
% patternCustom(Gy, theta, phi);
% figure()
% patternCustom(Gz, theta, phi);

% %% Draq Directivity against frequency
% dire_freq=readtable('directivity_dB_against_frequency.txt');



% %% Attempt to draw Gx Gy Gz from G_theta G_phi
% simu_g_theta=simu.Gain_Theta_;
% simu_g_phi=simu.Gain_Phi_;
% %convert to linear
% simu_g_theta=10.^(simu_g_theta/10);
% simu_g_phi=10.^(simu_g_phi/10);
% %calculation
% simu_g_z=-simu_g_theta./sin(theta_simu/360*2*pi);
% simu_g_x=simu_g_theta./sin(pi/2-theta_simu/360*2*pi).*cos(phi_simu/360*2*pi)+simu_g_phi.*sin(phi_simu/360*2*pi);
% simu_g_y=simu_g_theta./sin(pi/2-theta_simu/360*2*pi).*sin(phi_simu/360*2*pi)+simu_g_phi.*cos(phi_simu/360*2*pi);
% 
% %convert to dB
% simu_g_z=10*log10(abs(simu_g_z));
% simu_g_x=10*log10(abs(simu_g_x));
% simu_g_y=10*log10(abs(simu_g_y));
% figure()
% simu_g_z(find(simu_g_z==inf))=-90;
% patternCustom(simu_g_z, theta_simu, phi_simu);
