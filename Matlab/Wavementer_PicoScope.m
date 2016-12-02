% Author: Ana Andres-Arroyo (aa938)
% Reads long term array (.lta) files from the HighFinesse wavemeter
% Plots them and fits a region to a linear fit
% Calculates the speed with which the SolsTiS laser can change wavelength

clc
clear
close all
figure_handles = {};

directory = 'R:\3-Temporary\aa938\';
% directory = 'R:\aa938\NanoPhotonics\Laboratory\';
folder_name = '';
directory = [directory folder_name];
file_names = '.txt';

%% Select files
% *************************************************************************
[file_names, directory, ~] = uigetfile('.txt',...
                                      'Wavemeter & PicoScope files', ...
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
    clear file_names
    file_names{1} = temporary;
    clear temporary;
end
number_of_files = size(file_names,2);

%% Select files
% *************************************************************************
raw_data = cell(size(file_names));
wavelengths = cell(size(file_names));
header = 8;
for i = 1:1:number_of_files
    raw_data{i} = dlmread([directory file_names{i}], '\t', header, 0);
    raw_data{i}(:,end) = [];
    % raw_data{i}(1,:) = wavelength (nm)
    % raw_data{i}(2:end,:) = intensity (V)
    wavelengths{i} = raw_data{i}(1,:);
    raw_data{i}(1,:) = [];
    % raw_data{i}(:,:) = intensity (V)
    
end

%% Plot
% *************************************************************************
close all
figure_handles = {};

figure_handles{end+1} = figure('Units','normalized','Position',[0.05 0.1 0.7 0.7]);
for i = 1:1:number_of_files
    subplot(number_of_files, 1, i)
    for j = 1:1:size(raw_data{i},2)
        plot(1:1:size(raw_data{i}(:,j)), raw_data{i}(:,j), ...
            'LineWidth', 1), hold all
    end
    grid on
%     xlabel('PicoScope measurement number')
    ylabel('Intensity (V)')
    title('PicoScope data')
    if size(raw_data{i},2) < 10
        legend(cellstr(num2str(wavelengths{i}')), 'Location', 'SE')
    end
    
end


%% Average
% *************************************************************************
averaged_data = cell(size(raw_data));
for i = 1:1:number_of_files
    averaged_data{i} = mean(raw_data{i},1);
end

figure_handles{end+1} = figure('Units','normalized','Position',[0.25 0.15 0.7 0.7]);
for i = 1:1:number_of_files
    plot(wavelengths{i}, averaged_data{i}, ...
        'LineWidth', 1), hold all
end

grid on
xlabel('Wavelength (nm)')
ylabel('Mean Intensity (V)')
title('Wavemeter and PicoScope data')
if number_of_files < 15
    legend(file_names, 'Location', 'SE', 'interpreter', 'none')
end





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