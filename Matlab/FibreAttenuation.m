% Created by Ana Andres-Arroyo (aa938)
% Reads txt files from the spectrum analysier and calculates the
% attenuation coefficient of the fibre

clc
clear
% close all

%% READING DATA
% *************************************************************************

% specify default path
file_path = 'R:\aa938\NanoPhotonics\Laboratory\';
% file_path = 'R:\3-Temporary\aa938\';

% pop up window to choose the file(s) to read
[file_name, file_path, ~] = uigetfile('.txt',...
                                      'Spectrum Analyser Files to Read (use CTRL to select multiple files)',...
                                      file_path,...
                                      'MultiSelect','on');
file_path_save = file_path;
if strfind(file_path, 'Laboratory')
    title_cell{1} = file_path(strfind(file_path, 'Laboratory')+11:end-1);
elseif strfind(file_path, 'aa938')
    title_cell{1} = file_path(strfind(file_path, 'aa938')+6:end-1);
else
    title_cell{1} = file_path;
end
    
% if just a single file is selected, change the variable file_names
% from a string to a cell so the loops indexing will work
if isa(file_name,'char')
    temporary = file_name;
    clear file_name
    file_name{1} = temporary;
    clear temporary;
end

number_of_files = size(file_name,2);

% reading the data from the txt files
measurement_numbers = zeros(size(file_name));
raw_data = zeros(1001,2,number_of_files);
for i = 1:1:number_of_files
    measurement_numbers(i) = str2double(file_name{i}(2:5));
    raw_data(:,:,i) = dlmread([file_path file_name{i}], ',', [3, 0, 1003, 1]);
    % raw_data(:,1,i) = wavelength (nm)
    % raw_data(:,2,i) = spectral density (dBm/nm)
end

%% READING INFO FILE
% *************************************************************************

[file_name_info, file_path, ~] = uigetfile('.txt',...
                                          'Select one info file to read', ...
                                          file_path, ...
                                          'MultiSelect','off');
file_path_save = file_path;
info_data = dlmread([file_path file_name_info], '\t', 1, 0);
% info_data(:,1) = measurement number
% info_data(:,2) = fibre length (cm)

[~,measurement_indices,~] = intersect(info_data(:,1), measurement_numbers);
fibre_length = info_data(measurement_indices,2);

%% PLOTTING THE SPECTRA
% *************************************************************************
% figure('Units','normalized','Position',[0.2 0.15 0.7 0.7])
% legend_cell = {};
% for i = 1:1:number_of_files
%     h(i) = plot(raw_data(:,1,i), raw_data(:,2,i), ...
%                 'LineWidth', 2); hold all    
%     legend_cell{end+1} = [num2str(measurement_numbers(i), '%.0f') ...
%         ' // L = ' num2str(fibre_length(i), '%07.3f') ' cm'];
% end
% grid on
% xlabel('Wavelength (nm)')
% ylabel('Spectral Density (dBm/nm)')
% title(title_cell, 'interpreter', 'none')
% legend(legend_cell, 'Location', 'SE')
% set(gca, 'FontSize', 12)

%% CALCULATING THE ATTENUATION
% *************************************************************************

wavelengths = max(min(raw_data(:,1,:))) : 0.1 : min(max(raw_data(:,1,:)));
% wavelengths = 780:10:810; % nm
power = zeros(number_of_files, max(size(wavelengths)));
attenuation = zeros(size(wavelengths));
for j = 1:1:max(size(wavelengths))
    clc, disp([num2str(j/max(size(wavelengths))*100, '%.1f') '%'])
    for i = 1:1:number_of_files
        power(i,j) = spline(raw_data(:,1,i), raw_data(:,2,i), wavelengths(j));
    end
    polyfit_coefficients{j} = polyfit(fibre_length, power(:,j), 1);
    attenuation(j) = -polyfit_coefficients{j}(1)*100;
end

%% PLOTTING THE POWER vs. LENGTH
% *************************************************************************
figure('Units','normalized','Position',[0.05 0.1 0.85 0.7])
legend_wavelength = {};
for j = 1:1:max(size(wavelengths))
    h = plot(fibre_length, power(:,j), ...
             '.', 'MarkerSize', 12, ...
             'LineWidth', 1); hold all
    legend_wavelength{end+1} = [num2str(wavelengths(j), '%.0f') ' nm - data'];
    text(fibre_length + 2, ...
         power(:,j), ...
         strread(num2str(measurement_numbers),'%s'), ...
         'HorizontalAlignment', 'center',...
         'Color', 'k')
    f = plot(fibre_length, ...
             fibre_length*polyfit_coefficients{j}(1) + polyfit_coefficients{j}(2), ...
             'LineWidth', 1, ...
             'Color', h.Color); hold all  
    legend_wavelength{end+1} = [num2str(wavelengths(j), '%.0f') ' nm' ...
        ' - fit: \alpha = ' num2str(attenuation(j), '%.2f') ' dB/m'];
end

grid on
xlabel('Fibre Length (cm)')
ylabel('Spectral Density (dBm/nm)')
title(title_cell, 'interpreter', 'none')
legend(legend_wavelength, 'Location', 'NE')
set(gca, 'FontSize', 12)

%% PLOTTING THE ATTENUATION vs. WAVENLENGTH
% *************************************************************************

figure('Units','normalized','Position',[0.05 0.1 0.7 0.7])
plot(wavelengths, attenuation, 'LineWidth', 1), hold all
plot(wavelengths, absorption('D2O', wavelengths), '-.', 'LineWidth', 1), hold all
plot(wavelengths, absorption('H2O', wavelengths), '-.', 'LineWidth', 1), hold all
grid on
legend('Fibre Attenuation', 'D2O bulk absorption')
xlabel('Wavelength (nm)')
ylabel('Attenuation \alpha (dB/m)')
title(title_cell, 'interpreter', 'none')
set(gca, 'FontSize', 12)
% ylim([-1,3.5])

%% CONTOUR PLOT
% *************************************************************************
% % wavelengths = max(min(raw_data(:,1,:))) : 0.1 : min(max(raw_data(:,1,:)));
% wavelengths = 770 : 0.1 : min(max(raw_data(:,1,:)));
% power_spectral_density = zeros(max(size(wavelengths)), number_of_files);
% for i = 1:1:number_of_files
%     power_spectral_density(:,i) = spline(raw_data(:,1,i), raw_data(:,2,i), wavelengths);
% end
% figure('Units','normalized','Position',[0.02 0.065 0.7 0.7])
% contour_values = power_spectral_density';
% contour_levels = linspace(min(min(contour_values)), max(max(contour_values)), 100);
% contourf(wavelengths, ...
%          measurement_numbers,...
%          contour_values, ...
%          'LineStyle', 'none',...
%          'LevelListMode', 'manual', ...
%          'LevelList', contour_levels);
% colormap(jet)
% colorbar
% xlabel('Wavelength (nm)')
% ylabel('Measurement Number')
% title({'Spectral Density (dBm/nm)', title_cell{1}}, 'interpreter', 'none')
% set(gca, 'FontSize', 12)


