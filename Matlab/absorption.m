function abs = absorption(medium, wavelength)
% clc
% clear
% close all

file_directory = 'R:\aa938\NanoPhotonics\Matlab\';

% figure
% legend_absorption = {};

if strfind(medium, 'H2O')
    file_name_H2O = 'H2O_loss.txt';
    data_H2O = dlmread([file_directory file_name_H2O], '\t', 1, 0);
%     semilogy(data_H2O(:,1), data_H2O(:,2), 'LineWidth', 2), hold all
%     legend_absorption{end+1} = 'H2O';
    data = data_H2O;

elseif strfind(medium, 'D2O')
    file_name_D2O = 'D2O_loss.txt';
    data_D2O = dlmread([file_directory file_name_D2O], '\t', 1, 0);
%     semilogy(data_D2O(:,1), data_D2O(:,2), 'LineWidth', 2), hold all
%     legend_absorption{end+1} = 'D2O';
    data = data_D2O;
end

% legend(legend_absorption, 'Location', 'NW')
% grid on
% xlabel('Wavelength (nm)')
% ylabel('Absorption (dB/m)')
% xlim([200,2000])

% data(:,1) = wavelength (nm)
% data(:,2) = loss (dB/m)
abs = spline(data(:,1), data(:,2), wavelength); % dB/m

end