clc
clear
close all
figure_handles = {};

% directory = 'R:\3-Temporary\aa938\';
directory = 'R:\aa938\NanoPhotonics\Laboratory\';
folder_name = '2016.11.17 - PBG Au NR attenuation\spectra';
directory = [directory folder_name];
file_names = '.txt';
delimiter = ',«x20»';

%% READ wavelength
% *************************************************************************
[file_names, directory, ~] = uigetfile('.txt',...
                                      'WAVELENGTH data', ...
                                      [directory file_names], ...
                                      'MultiSelect','off');
                                  
directory_save = directory;
if strfind(directory, 'Laboratory')
    folder_name = directory(strfind(directory, 'Laboratory')+11:end-1);
elseif strfind(directory, 'aa938')
    folder_name = directory(strfind(directory, 'aa938')+6:end-1);
else
    folder_name = directory;
end
                                  
fileID = fopen([directory file_names], 'r');
raw_data = textscan(fileID, '%f','Delimiter', delimiter);
fclose(fileID);
indices = 7:6:size(raw_data{1},1);
wavelengths = raw_data{1}(indices);
clear raw_data

%% READ dBm
% *************************************************************************
file_names = '.txt';
[file_names, directory, ~] = uigetfile('.txt',...
    'INTENSITY data: multiple files possible', ...
    [directory file_names], ...
    'MultiSelect','on');

if isa(file_names,'char')
    temporary = file_names;
    clear file_name
    file_names{1} = temporary;
    clear temporary;
end
number_of_files = size(file_names,2);

intensity = zeros(number_of_files, max(size(wavelengths)));
delimiter = ',«x20»';
for i = 1:1:number_of_files
    fileID = fopen([directory file_names{i}], 'r');
    raw_data = textscan(fileID, '%f','Delimiter', delimiter);
    fclose(fileID);
    intensity(i,:) = raw_data{1}(indices);
    clear raw_data
end
%% PLOTTING spectra
% *************************************************************************
figure_handles{end+1} = figure('Units','normalized','Position',[0.05 0.1 0.7 0.7]);
legend_cell = {};
for i = 1:1:number_of_files
    plot(wavelengths, intensity(i,:), ...
        'LineWidth', 1), hold all
    legend_cell{end+1} = file_names{i}(1:end-4);
end
grid on
legend(legend_cell, 'Location', 'SE')
xlabel('Wavelength (nm)')
ylabel('Spectral Density (dBm/nm)')
title(folder_name, 'interpreter', 'none')

%% CALCULATE attenuation
% *************************************************************************
fibre_length = 17; % cm
[ref_index, sample_index] = dialog_two_lists('Select data:', ...
                                      'Reference:', file_names, ...
                                      'Sample:', file_names);
sample_attenuation = -1/(fibre_length/100) * ...
    (intensity(sample_index,:) - intensity(ref_index,:)); % dB/m


% PLOTTING attenuation
% *************************************************************************
% figure_handles{end+1} = figure('Units','normalized','Position',[0.1 0.2 0.7 0.7]);
% legend_cell = {};
plot(wavelengths, sample_attenuation, ...
    'LineWidth', 1), hold all
legend_cell{end+1} = ['ref: ' file_names{ref_index}(1:end-4) ', ' ...
    'sample: ' file_names{sample_index}(1:end-4)];
grid on
xlabel('Wavelength (nm)')
ylabel('Attenuation (dB/m)')
title(folder_name, 'interpreter', 'none')
legend(legend_cell, 'Location', 'SE')

%% CALCULATE optical density
% *************************************************************************

sample_OD = 10 * sample_attenuation / 100; % 1/cm

figure_handles{end+1} = figure('Units','normalized','Position',[0.1 0.2 0.7 0.7]);
legend_cell = {};
plot(wavelengths, sample_OD, ...
    'LineWidth', 1), hold all
legend_cell{end+1} = ['ref: ' file_names{ref_index}(1:end-4) ', ' ...
    'sample: ' file_names{sample_index}(1:end-4)];

% read manufacturer optical density
directory = 'R:\aa938\NanoPhotonics\Laboratory\2016.11.16 - PBG HCF prep for Au NR\';
OD_file_name = 'nanocomposix-spectra.csv';
% [OD_file_name, directory, ~] = uigetfile('.csv',...
%                                       'OPTICAL DENSITY data', ...
%                                       [directory OD_file_name], ...
%                                       'MultiSelect','off');
                                  
manufacturer_OD = dlmread([directory OD_file_name], ',', 0, 0);
plot(manufacturer_OD(:,1), manufacturer_OD(:,2))
legend_cell{end+1} = OD_file_name(1:end-4);

grid on
xlabel('Wavelength (nm)')
ylabel('Optical Density (cm^{-1})')
title(folder_name, 'interpreter', 'none')
legend(legend_cell, 'Location', 'SE')


%% SAVING FIGURES
% *************************************************************************
menu_save_figures = 1;
menu_save_figures = menu('Save Figures?', 'NO', 'YES');
if menu_save_figures == 2    
    for i = 1:1:max(size(figure_handles))
        if findobj(figure_handles{i}) ~= 0
            figure_save = figure_handles{i};
            file_name_save = file_names{1}(1:end-4);
            
            figure(figure_save)
            pause(0.1)
            [file_name_save,directory_save,~] = uiputfile(['.' 'png'],...
                'File to Save the Figure',[directory_save file_name_save]);
            hgexport(figure_save, [directory_save file_name_save], hgexport('factorystyle'), 'Format', 'png')
            file_name_save = strrep(file_name_save, 'png', 'fig');    
            saveas(figure_save, [directory_save file_name_save], 'fig');

        end
    end
end
