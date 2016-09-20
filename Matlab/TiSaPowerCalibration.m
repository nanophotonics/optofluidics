clc
clear
close all

file_directory = 'R:\aa938\NanoPhotonics\Laboratory\2016.09.19 - TiSa power calibration 800 nm\';
file_name = '2016.09.19 - TiSa power calibration 800 nm.txt';

header_rows = 1;

data = dlmread([file_directory, file_name], '\t', header_rows, 0);
% data(:,1) = angle of the waveplate (deg)
% data(:,2) = power of the Ti:Sa laser after the polariser (W)
angle_degrees = data(:,1);
% angle_radians = angle_degrees*pi/180;
power_watts = data(:,2);

figure('Units','normalized','Position',[0.2 0.1 0.7 0.7]);
plot(angle_degrees, power_watts, '.k', 'MarkerSize', 16); hold on
set(gca, 'FontSize', 14)
legend('measured')

fitting_function = fittype('a*(sin(b*x*pi/180+c))^2');
starting_coefficients = [4,2,0.09];
[fit_coefficients, goodness_of_fit] = fit(angle_degrees, power_watts, ...
    fitting_function, 'StartPoint', starting_coefficients)
plot(fit_coefficients)

xlabel('Angle (degrees)')
ylabel('Ti:Sa power (W)') 
title(file_name)