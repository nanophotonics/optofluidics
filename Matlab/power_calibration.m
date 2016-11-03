function power = power_calibration(reference_angle, reference_power, angle)
% clc
% clear
% close all

file_directory = 'R:\aa938\NanoPhotonics\Laboratory\2016.10.12 - TiSa power calibration 767 nm\';
file_name = '2016.10.12 - TiSa power calibration 767 nm - fit.txt';
% [file_name, file_directory, ~] = uigetfile('.txt',...
%                                            'Select a calibration file to read',...
%                                            [file_directory file_name],...
%                                            'MultiSelect','off');
header = 4;
file_id = fopen([file_directory file_name], 'r');
for i = 1:1:header
    fgets(file_id);
end
% power(x) = a*(sin(b*x*pi/180+c))^2
% x ---> angle of the waveplate (degrees)
% power ---> power of the Ti:Sa laser (W)
coefficients = zeros(3,3);
coefficients(:,1) = fscanf(file_id, 'a =  %f  (%f, %f)  W');
fgets(file_id);
coefficients(:,2) = fscanf(file_id, 'b =  %f  (%f, %f)  1/rad');
fgets(file_id);
coefficients(:,3) = fscanf(file_id, 'c =  %f  (%f, %f)  rad');
fgets(file_id);
fclose(file_id);

disp(reference_angle)
disp(reference_power)
disp(angle)
power = 1; % same units as the reference_power

end