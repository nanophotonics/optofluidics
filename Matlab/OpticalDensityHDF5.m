warning('off', 'MATLAB:plot:IgnoreImaginaryXYPart')

clc
clear
close all
figure_handles = {};

colour_type = {'DEFAULT', ...
               'parula', 'jet', 'hsv', 'cool', ...
               'spring', 'summer', 'autumn', 'autumn reversed', 'winter', ...
               'gray', 'copper',...
               'red', 'green', 'aqua', 'blue', 'purple',...
               };
selected_colour = 1;

%% default directory
DirectoryRead = 'R:\aa938\NanoPhotonics\Laboratory\';
FolderRead = '2017.01.18 - Au NR spectra\';
FolderPathRead = [DirectoryRead FolderRead];
FileNameRead = '';

%% prompt to choose file
[FileNameRead, FolderPathRead, ~] = uigetfile('.h5',...
    'H5 file to read:',[FolderPathRead FileNameRead],'MultiSelect','off');
%     'H5 file to read:','MultiSelect','off');
FolderPathSave = FolderPathRead;
slash_index = strfind(FolderPathRead, '\');
FolderRead = FolderPathRead(slash_index(end-1)+1:end-1);
FilePathRead = [FolderPathRead FileNameRead];

%% selecting an h5 group
GroupInfo = h5info(FilePathRead);
while max(size(GroupInfo.Groups)) > 0
    menu_group = menu('Choose group:', GroupInfo.Groups.Name);
    GroupInfo = h5info(FilePathRead, GroupInfo.Groups(menu_group).Name);
end


%% reading the data
number_of_spectra = size(GroupInfo.Datasets,1);
number_of_wavelengths = GroupInfo.Datasets(1).Dataspace.Size;
DatasetPaths = cell(number_of_spectra,1);
DatasetNames = cell(number_of_spectra,1);
description = cell(number_of_spectra,1);
counts = zeros(number_of_spectra, number_of_wavelengths);
wavelengths = zeros(number_of_spectra, number_of_wavelengths);
background = zeros(number_of_spectra, number_of_wavelengths);
reference = ones(number_of_spectra, number_of_wavelengths);
for i = 1:1:number_of_spectra
    DatasetPaths{i} = [GroupInfo.Name '/' GroupInfo.Datasets(i).Name];
    slash_index = strfind(DatasetPaths{i}, '/');
    DatasetNames{i} = DatasetPaths{i}(slash_index(end)+1:end);
    try
        description(i) = h5readatt(GroupInfo.Filename, DatasetPaths{i}, 'description');
    catch
    end
    counts(i,:) = h5read(GroupInfo.Filename, DatasetPaths{i});
    wavelengths(i,:) = h5readatt(GroupInfo.Filename, DatasetPaths{i}, 'wavelengths');
    try
        background(i,:) = h5readatt(GroupInfo.Filename, DatasetPaths{i}, 'background');
    catch
    end
    try
        reference(i,:) = h5readatt(GroupInfo.Filename, DatasetPaths{i}, 'reference');
    catch
    end
end



%% sort the data array by spectrum number
index_sort = zeros(number_of_spectra,1);
for i = 1:1:number_of_spectra
    underscore_index = strfind(DatasetNames{i}, '_');
    index_sort(i) = str2double(DatasetNames{i}(underscore_index(end)+1:end))+1;
end
DatasetNames(index_sort,:) = DatasetNames;
counts(index_sort,:) = counts;
wavelengths(index_sort,:) = wavelengths;
background(index_sort,:) = background;
reference(index_sort,:) = reference;

%% plot raw data -----
figure_handles{end+1} = figure('Units','normalized','Position',[0.01 0.085 0.75 0.8]);
for i = 1:1:number_of_spectra
    h(i) = plot(wavelengths(i,:), counts(i,:), ...
        'LineWidth', 2); hold on
end
grid on
set(gca, 'FontSize', 18)
% legend(description, 'Location', 'NEO', 'interpreter', 'none')
legend(DatasetNames, 'Location', 'NEO', 'interpreter', 'none')
title({FolderRead, 'Raw Data'}, 'interpreter', 'none')
xlabel('Wavelength (nm)')
ylabel('Intensity (counts)')
xlim([300,1050])

%% plot styling
[selected_colour, ~] = listdlg('PromptString', 'Colour scheme:',...
                           'SelectionMode', 'single', ...
                           'ListString', colour_type,...
                           'InitialValue', selected_colour);

for i = 1:1:number_of_spectra
    if selected_colour > 1 
        colour_RGB = colour_gradient(i, number_of_spectra, colour_type(selected_colour));
        h(i).Color = colour_RGB;  
    end
    h(i).MarkerSize = 1;
    h(i).LineStyle = '-';
    h(i).LineWidth = 2;
end

%% data analysis options
analysis_options = {'Background', 'Reference', 'Optical Density'};
selected_analysis = [1,2,3];
[selected_analysis, ~] = listdlg(...
    'PromptString', 'Data Analysis:', ...
    'SelectionMode', 'multiple', ...
    'ListString', analysis_options, ...
    'InitialValue', selected_analysis, ...
    'OKString', 'OK', ...
    'CancelString', 'NONE');

%% apply background and reference corrections
data = counts;
reference_bg = reference;
ytext = 'Intensity (counts)';
selected_spectra = 1:1:number_of_spectra;
title_analysis = 'Correction: ';
if any(strcmp(analysis_options(selected_analysis), 'Background'))
    title_analysis = [title_analysis ' Background,'];
end
if any(strcmp(analysis_options(selected_analysis), 'Reference'))
    title_analysis = [title_analysis ' Reference,'];
    ytext = 'Transmission';
end
for i = number_of_spectra:-1:1
    if or(...
            and(any(strcmp(analysis_options(selected_analysis), 'Background')), ...
            isempty(find(background(i,:),1))),...
            and(any(strcmp(analysis_options(selected_analysis), 'Reference')), ...
            isempty(find(reference(i,:) - 1,1))))
        disp(selected_spectra(i)-1)
        selected_spectra(i) = [];
    end
    if any(strcmp(analysis_options(selected_analysis), 'Background'))
        data(i,:) = counts(i,:) - background(i,:);
        if find(reference(i,:) - 1)
            reference_bg(i,:) = reference(i,:) - background(i,:);
        end
    end
    if any(strcmp(analysis_options(selected_analysis), 'Reference'))
        data(i,:) = data(i,:) ./ reference_bg(i,:);
    end
end

%% optical density calculation
if any(strcmp(analysis_options(selected_analysis), 'Optical Density'))
    sample_length = 1; % cm
    sample_dilution = 1; % none
    input_title = 'Optical Density Parameters:';
    input_data = {'Sample Length (cm):', 'Sample Dilution (times):'};
    default_values = {num2str(sample_length), num2str(sample_dilution)};
    dlg_options.WindowStyle = 'normal'; dlg_options.Resize = 'on'; dim = [1 80];
    answer = inputdlg(input_data, input_title, dim, default_values, dlg_options);
    sample_length = str2double(answer{1});
    sample_dilution = str2double(answer{2});  
    
    data = -log10(data);
    title_analysis = 'Optical Density';
    ytext = 'Optical Density (cm^{-1})';
%     ytext = 'Optical Density';
    
end

%% select spectra to plot
selected_index = 1:1:max(size(selected_spectra));
[selected_index, ~] = listdlg('PromptString', 'Spectra to plot',...
                               'SelectionMode', 'multiple', ...
                               'ListString', DatasetNames(selected_spectra),...
                               'InitialValue', selected_index);
selected_spectra = selected_spectra(selected_index);

%% plot selected corrected data -----
if or(...
        max(size(selected_analysis)) > 0, ...
        max(size(selected_spectra)) < number_of_spectra)
    figure_handles{end+1} = figure('Units','normalized','Position',[0.23 0.085 0.75 0.8]);
    for i = selected_spectra
        h(i) = plot(wavelengths(i,:), data(i,:), ...
            'LineWidth', 2); hold on
    end
    grid on
    set(gca, 'FontSize', 18)
%     legend(description(selected_spectra), 'Location', 'NEO', 'interpreter', 'none')
    legend(DatasetNames(selected_spectra), 'Location', 'NEO', 'interpreter', 'none')
    title({FolderRead, title_analysis}, 'interpreter', 'none')
    xlabel('Wavelength (nm)')
    ylabel(ytext)
    xlim([300,1050])
%     ylim([-0.3,1])

    %% plot styling
    [selected_colour, ~] = listdlg('PromptString', 'Colour scheme:',...
                               'SelectionMode', 'single', ...
                               'ListString', colour_type,...
                               'InitialValue', selected_colour);

    for i = selected_spectra
        if selected_colour > 1 
            colour_RGB = colour_gradient(i, number_of_spectra, colour_type(selected_colour));
            h(i).Color = colour_RGB;  
        end
        h(i).MarkerSize = 1;
        h(i).LineStyle = '-';
        h(i).LineWidth = 2;
    end
    
    %% read manufacturer optical density
    directory = 'R:\aa938\NanoPhotonics\Laboratory\2016.11.16 - PBG HCF prep for Au NR\';
    OD_file_name = 'nanocomposix-spectra.csv';
    % [OD_file_name, directory, ~] = uigetfile('.csv',...
    %                                       'OPTICAL DENSITY data', ...
    %                                       [directory OD_file_name], ...
    %                                       'MultiSelect','off');

    manufacturer_OD = dlmread([directory OD_file_name], ',', 0, 0);
    plot(manufacturer_OD(:,1), manufacturer_OD(:,2), 'LineWidth', 1)
    legend_cell = DatasetNames(selected_spectra);
    legend_cell{end+1} = OD_file_name(1:end-4);
    legend(legend_cell, 'Location', 'NEO', 'interpreter', 'none')
end


%% saving figures
% *************************************************************************
menu_save_figures = 1;
menu_save_figures = menu('Save Figures?', 'NO', 'YES');

if menu_save_figures == 2    
    for i = 1:1:max(size(figure_handles))
        if findobj(figure_handles{i}) ~= 0
            figure_save = figure_handles{i};
            FileNameSave = FileNameRead(1:end-3);
            
            figure(figure_save)
            pause(0.1)
            [FileNameSave,FolderPathSave,~] = uiputfile(['.' 'png'],...
                'File to Save the Figure',[FolderPathSave FileNameSave]);
            hgexport(figure_save, [FolderPathSave FileNameSave], hgexport('factorystyle'), 'Format', 'png')
            FileNameSave = strrep(FileNameSave, 'png', 'fig');    
            saveas(figure_save, [FolderPathSave FileNameSave], 'fig');

        end
    end
end