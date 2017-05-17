% Created by Ana Andres-Arroyo (aa938)
% Calculates the power at a certain angle based on the half-waveplate
% calibration data for the Ti:Sa laser.

% RESULTS ARE INACCURATE FOR SMALL REFERENCE POWERS
% use power_fitted instead

function power = power_spline(reference_angle, reference_power, angle, file_path)

% reference_angle in degrees
% reference_power in W or mW.
% angle in degrees: can be a single number or an array of numbers

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

% file_directory = '';
% file_name = 'power_calibration_data.txt';
% [file_name, file_directory, ~] = uigetfile('.txt',...
%                                            'Select a calibration file to read',...
%                                            [file_directory file_name],...
%                                            'MultiSelect','off');
% file_path = [file_directory, file_name];

header_rows = 2;
data = dlmread(file_path, '\t', header_rows, 0);
% data(:,1) = angle of the waveplate (deg)
% data(:,2) = power of the Ti:Sa laser after the polariser (W)
calibration_angle = data(:,1); % degrees
calibration_power = data(:,2); % watts

power = reference_power / ...
    spline(calibration_angle, calibration_power, reference_angle) * ...
    spline(calibration_angle, calibration_power, angle); % same units as the reference_power

% power = reference_power * ...
%     (sin(pi/180*b*angle+c).^2) / ...
%     (sin(pi/180*b*reference_angle+c).^2); % same units as the reference_power

end