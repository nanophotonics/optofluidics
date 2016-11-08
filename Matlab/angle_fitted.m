% Created by Ana Andres-Arroyo (aa938)
% Calculates the angle required to achieve a certain power based on the 
% half-waveplate calibration fitting parameters for the Ti:Sa laser.

function required_angle = angle_fitted(reference_angle, reference_power, power)

% reference_angle in degrees
% reference_power in W or mW.
% power in W or mW: can be a single number or an array of numbers

% clc
% clear
% close all
% 
% power = 0:50:1500;
% reference_angle = 30;
% reference_power = 16;
% input_title = 'Parameters'; 
% input_data = {'Reference Ange (deg):',...
%               'Reference Power (mW):', ...
%               };
% default_values = {num2str(reference_angle),...
%                   num2str(reference_power),...
%                   };
% dlg_options.WindowStyle = 'normal'; dlg_options.Resize = 'on'; dim = [1 60];
% answer = inputdlg(input_data, input_title, dim, default_values, dlg_options);
% reference_angle = str2double(answer{1});   
% reference_power = str2double(answer{2});   


file_directory = '';
file_name = 'power_calibration_fit.txt';
% [file_name, file_directory, ~] = uigetfile('.txt',...
%                                            'Select a calibration file to read',...
%                                            [file_directory file_name],...
%                                            'MultiSelect','off');
header_rows = 5;
file_id = fopen([file_directory file_name], 'r');
for i = 1:1:header_rows
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

b = coefficients(1,2);
c = coefficients(1,3);

possible_angles = 0:0.1:90;
calibrated_power = reference_power * ...
    (sin(pi/180*b*possible_angles+c).^2) / ...
    (sin(pi/180*b*reference_angle+c).^2); % same units as the reference_power

[~,i_min] = min(calibrated_power);
[~,i_max] = max(calibrated_power);
required_angle = zeros(size(power));
for i = 1:1:max(size(power))
    [~,index_angle] = min(abs(calibrated_power(i_min:i_max) - power(i)));
    required_angle(i) = possible_angles(i_min + index_angle - 1);
end

figure
plot(possible_angles,calibrated_power), hold all
plot(required_angle, power, '.', 'MarkerSize', 20)
grid on
xlabel('Waveplate angle (degrees)')
ylabel('Power')
title(['b = ' num2str(b) ' 1/rad, c = ' num2str(c) ' rad'])
legend('Calibrated power', 'Required power', 'Location', 'NW')

end