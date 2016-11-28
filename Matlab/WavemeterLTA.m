% Author: Ana Andres-Arroyo (aa938)
% Reads long term array (.lta) files from the HighFinesse wavemeter
% Plots them and fits a region to a linear fit
% Calculates the speed with which the SolsTiS laser can change wavelength

clc
clear
close all
figure_handles = {};

% directory = 'R:\3-Temporary\aa938\';
directory = 'R:\aa938\NanoPhotonics\Laboratory\';
folder_name = '2016.11.28 - TiSa wavemeter Python\';
directory = [directory folder_name];
file_names = '.lta';

%% Select files
% *************************************************************************
[file_names, directory, ~] = uigetfile('.lta',...
                                      'Wavemeter long term array files', ...
                                      [directory file_names], ...
                                      'MultiSelect','on');
                                  
directory_save = directory;
if strfind(directory, 'Laboratory')
    folder_name = directory(strfind(directory, 'Laboratory')+11:end-1);
elseif strfind(directory, 'aa938')
    folder_name = directory(strfind(directory, 'aa938')+6:end-1);
else
    folder_name = directory;
end

% if just a single file is selected, change the variable file_names
% from a string to a cell so the loops indexing will work
if isa(file_names,'char')
    temporary = file_names;
    clear file_name
    file_names{1} = temporary;
    clear temporary;
end
number_of_files = size(file_names,2);

%% Select files
% *************************************************************************
raw_data = cell(size(file_names));
header = 38;
for i = 1:1:number_of_files
    raw_data{i} = dlmread([directory file_names{i}], '\t', header, 0);
    % raw_data{i}(:,1) = time (ms)
    % raw_data{i}(:,2) = wavelength (nm)
end

%% Plot
% *************************************************************************
close all
figure_handles = {};

figure_handles{end+1} = figure('Units','normalized','Position',[0.05 0.1 0.7 0.7]);
for i = 1:1:number_of_files
    plot(raw_data{i}(:,1)/1000, raw_data{i}(:,2), ...
        'LineWidth', 1), hold all
end
grid on
xlabel('Time (s)')
ylabel('Wavelength (nm)')
title('Wavelength meter data')
legend(file_names, 'Location', 'SE')
ylim([700,1000])

%% Linear fit
% *************************************************************************

t_min = 500; % ms
t_max = 5500; % ms
files_to_fit = 1:4;

input_title = 'Time selection for linear fit'; 
input_data = {'Start t (s):','End t (s):', 'Files to fit (0 = all):'};
default_values = {num2str(t_min/1000), num2str(t_max/1000), '0'};
dlg_options.WindowStyle = 'normal'; dlg_options.Resize = 'on'; dim = [1 80];
answer = inputdlg(input_data,input_title,dim,default_values,dlg_options);
t_min = str2double(answer{1})*1000; %ms
t_max = str2double(answer{2})*1000; % ms
if str2double(answer{3}) > 0
    files_to_fit = round(str2double(answer{3}));
end

linear_data = cell(size(raw_data));
fit_coefficients = cell(size(raw_data));
fit_structre = cell(size(raw_data));
for i = files_to_fit
    [i_min, ~] = find(raw_data{i}(:,1) >= t_min);
    [i_max, ~] = find(raw_data{i}(:,1) <= t_max);
    i_intersect = intersect(i_min, i_max);
    linear_data{i} = raw_data{i}(i_intersect,:);
    
    [fit_coefficients{i}, fit_structure{i}] = polyfit(linear_data{i}(:,1), linear_data{i}(:,2), 1);
    
end
    
% figure_handles{end+1} = figure('Units','normalized','Position',[0.05 0.1 0.7 0.7]);
% legend_fit = {};
figure(figure_handles{end});
for i = files_to_fit
    plot(linear_data{i}(:,1)/1000, linear_data{i}(:,2), ...
        'LineWidth', 1), hold all
    legend_fit{end+1} = file_names{i}(1:end-4);
    
    plot(linear_data{i}(:,1)/1000, polyval(fit_coefficients{i}, linear_data{i}(:,1)), ...
        'k--', 'LineWidth', 1), hold all
    legend_fit{end+1} = ['Slope = ' num2str(fit_coefficients{i}(1)*1000, '%.2f') ' nm/s'];
end
grid on
xlabel('Time (s)')
ylabel('Wavelength (nm)')
title('Wavelength meter data')
legend(legend_fit, 'Location', 'SE')


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