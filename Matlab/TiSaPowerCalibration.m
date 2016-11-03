clc
clear
close all

file_directory = 'R:\aa938\NanoPhotonics\Laboratory\2016.10.12 - TiSa power calibration 767 nm\';
file_name = '2016.10.12 - TiSa power calibration 767 nm.txt';
[file_name, file_directory, ~] = uigetfile('.txt',...
                                           'Select a calibration file to read',...
                                           [file_directory file_name],...
                                           'MultiSelect','off');
file_directory_save = file_directory;
                                       
                                       header_rows = 1;
data = dlmread([file_directory, file_name], '\t', header_rows, 0);
% data(:,1) = angle of the waveplate (deg)
% data(:,2) = power of the Ti:Sa laser after the polariser (W)
angle_degrees = data(:,1);
angle_radians = angle_degrees*pi/180;
power_watts = data(:,2);

figure('Units','normalized','Position',[0.2 0.1 0.7 0.7]);
plot(angle_degrees, power_watts, '.k', 'MarkerSize', 16); hold on
set(gca, 'FontSize', 14)
legend('measured', 'Location', 'SE')

fitting_function = fittype('a*(sin(b*x*pi/180+c))^2');
starting_coefficients = [4,2,1];
[fit_object, goodness_of_fit] = fit(angle_degrees, power_watts, ...
    fitting_function, 'StartPoint', starting_coefficients);
plot(fit_object)
legend('measured', 'fitted curve', 'Location', 'SE')
disp(fit_object)
disp(goodness_of_fit)
coefficient_values = coeffvalues(fit_object);
coefficient_confidence = confint(fit_object);

xlabel('Angle (degrees)')
ylabel('Ti:Sa power (W)') 
title(file_name(1:end-4))

%% ---
text_cell{1} = 'power(x) = a*(sin(b*x*pi/180+c))^2';
text_cell{2} = 'x ---> angle of the waveplate (degrees)';
text_cell{3} = 'power ---> power of the Ti:Sa laser (W)';
text_cell{4} = '';
text_cell{5} = ['a =  ' num2str(coefficient_values(1,1), '%.3f') ...
                '  (' num2str(coefficient_confidence(1,1), '%.3f'), ...
                ', ' num2str(coefficient_confidence(2,1), '%.3f'), ...
                ')  W'];
text_cell{6} = ['b =  ' num2str(coefficient_values(1,2), '%.3f') ...
                '  (' num2str(coefficient_confidence(1,2), '%.3f'), ...
                ', ' num2str(coefficient_confidence(2,2), '%.3f'), ...
                ')  1/rad'];
text_cell{7} = ['c =  ' num2str(coefficient_values(1,3), '%.3f') ...
                '  (' num2str(coefficient_confidence(1,3), '%.3f'), ...
                ', ' num2str(coefficient_confidence(2,3), '%.3f'), ...
                ')  rad'];
text('Units','normalized','Position',[0.08,0.92], ...
    'VerticalAlignment', 'top', 'String' , text_cell)

%% ---
file_name_save = [file_name(1:end-4) ' - fit.txt'];
[file_name_save,file_directory_save,~] = uiputfile('.txt',...
                                                   'Select a file to save the fitting parameters',...
                                                   [file_directory_save file_name_save]); 
file_id = fopen([file_directory_save file_name_save], 'wt');
for i = 1:1:size(text_cell,2)
    fprintf(file_id, [strrep(text_cell{i}, '%', 'percent') '\n']);
end
fclose(file_id);