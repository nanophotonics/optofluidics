warning('off', 'MATLAB:plot:IgnoreImaginaryXYPart')

new_figures = 2;
new_figures = menu('New Figures?', 'NO', 'YES');

if new_figures == 2
    clc
    clear
    close all
    figure_handles = {};
    new_figures = 2;
end

colour_type = {'DEFAULT', ...
               'parula', 'jet', 'hsv', 'cool', ...
               'spring', 'summer', 'autumn', 'autumn reversed', 'winter', ...
               'gray', 'copper',...
               'red', 'green', 'aqua', 'blue', 'purple',...
               };
selected_colour = 1;
legend_location = 'NEO';
legend_font_size = 8;
plot_font_size = 18;

%% default directory
% *************************************************************************
DirectoryRead = 'R:\aa938\NanoPhotonics\Laboratory\';
% DirectoryRead = 'R:\3-Temporary\aa938\';
% DirectoryRead = 'C:\Users\Ana Andres\Desktop\';
% FolderRead = '2017.01.18 - Au NR spectra cuvette\';
FolderRead = '2017.01.27 - Au NR spectra fibre\';
FolderPathRead = [DirectoryRead FolderRead];
FileNameRead = '';

%% prompt to choose file
% *************************************************************************
[FileNameRead, FolderPathRead, ~] = uigetfile('.h5',...
    'H5 file to read:',[FolderPathRead FileNameRead],'MultiSelect','off');
%     'H5 file to read:','MultiSelect','off');
FolderPathSave = FolderPathRead;
slash_index = strfind(FolderPathRead, '\');
FolderRead = FolderPathRead(slash_index(end-1)+1:end-1);
FilePathRead = [FolderPathRead FileNameRead];

%% selecting an h5 group
% *************************************************************************
GroupInfo = h5info(FilePathRead);
while max(size(GroupInfo.Groups)) > 0
    menu_group = menu('Choose group:', GroupInfo.Groups.Name);
    GroupInfo = h5info(FilePathRead, GroupInfo.Groups(menu_group).Name);
end

%% reading the data
% *************************************************************************
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
% *************************************************************************
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
description(index_sort,:) = description;

%% plot raw data -----
% *************************************************************************
new_raw_figure = new_figures;
if exist('figure_raw','var') == 0
    new_raw_figure = 2;
else
    if ishandle(figure_raw) == 0
        new_raw_figure = 2;
    end
end

if new_raw_figure == 2    
    figure_raw = figure('Units','normalized','Position',[0.01 0.085 0.75 0.8]);
    figure_handles{end+1} = figure_raw;
    legend_raw = {};
    h_raw = {};
end
figure(figure_raw);


for i = 1:1:number_of_spectra
    h_raw{end+1} = plot(wavelengths(i,:), counts(i,:), ...
        'LineWidth', 2); hold on
    legend_raw{end+1} = [DatasetNames{i} ' - ' ...
            description{i}];
end
grid on
set(gca, 'FontSize', plot_font_size)
legend(legend_raw, 'Location', legend_location, ...
    'FontSize', legend_font_size, 'interpreter', 'none')
title({FolderRead, 'Raw Data'}, 'interpreter', 'none')
xlabel('Wavelength (nm)')
ylabel('Intensity (counts)')
xlim([300,1050])

% plot styling
% *************************************************************************
[selected_colour, ~] = listdlg('PromptString', 'Colour scheme:',...
                           'SelectionMode', 'single', ...
                           'ListString', colour_type,...
                           'InitialValue', selected_colour);

for i = 1:1:max(size(h_raw))
    if selected_colour > 1 
        colour_RGB = colour_gradient(i, max(size(h_raw)), colour_type(selected_colour));
        h_raw{i}.Color = colour_RGB;  
    end
    h_raw{i}.MarkerSize = 1;
    h_raw{i}.LineStyle = '-';
    h_raw{i}.LineWidth = 2;
end

%% data analysis options
% *************************************************************************
analysis_options = {'Background', 'Reference', ...
    'Optical Density', 'Savitzky-Golay Filtering', 'Normalisation'};
selected_analysis = [1,2,3];
[selected_analysis, ~] = listdlg(...
    'PromptString', 'Data Analysis:', ...
    'SelectionMode', 'multiple', ...
    'ListString', analysis_options, ...
    'InitialValue', selected_analysis, ...
    'OKString', 'OK', ...
    'CancelString', 'NONE');

%% apply background and reference corrections
% *************************************************************************
data_bg_ref = counts;
reference_bg = reference;
ytext = 'Intensity (counts)';
corrected_spectra = 1:1:number_of_spectra;
title_analysis = 'Correction: ';
if any(strcmp(analysis_options(selected_analysis), 'Background'))
    title_analysis = [title_analysis 'Background. '];
end
if any(strcmp(analysis_options(selected_analysis), 'Reference'))
    title_analysis = [title_analysis 'Reference. '];
    ytext = 'Transmission';
end
for i = number_of_spectra:-1:1
    if or(...
            and(any(strcmp(analysis_options(selected_analysis), 'Background')), ...
            isempty(find(background(i,:),1))),...
            and(any(strcmp(analysis_options(selected_analysis), 'Reference')), ...
            isempty(find(reference(i,:) - 1,1))))
        disp(corrected_spectra(i)-1)
        corrected_spectra(i) = [];
    end
    if any(strcmp(analysis_options(selected_analysis), 'Background'))
        data_bg_ref(i,:) = counts(i,:) - background(i,:);
        if find(reference(i,:) - 1)
            reference_bg(i,:) = reference(i,:) - background(i,:);
        end
    end
    if any(strcmp(analysis_options(selected_analysis), 'Reference'))
        data_bg_ref(i,:) = data_bg_ref(i,:) ./ reference_bg(i,:);
    end
end
data_plot = data_bg_ref;

%% optical density calculation
% *************************************************************************
data_OD = data_bg_ref;
if any(strcmp(analysis_options(selected_analysis), 'Optical Density'))
    sample_length = 31; % cm
    sample_dilution = 20; % times
    input_title = 'Optical Density Parameters:';
    input_data = {'Sample Length (cm):', 'Sample Dilution in D2O (times):'};
    default_values = {num2str(sample_length), num2str(sample_dilution)};
    dlg_options.WindowStyle = 'normal'; dlg_options.Resize = 'on'; dim = [1 80];
    answer = inputdlg(input_data, input_title, dim, default_values, dlg_options);
    sample_length = str2double(answer{1});
    sample_dilution = str2double(answer{2});  
    
    data_OD = -log10(data_bg_ref);
    data_OD = data_OD / sample_length * sample_dilution;
    title_analysis = ['OD: L = ' ...
        num2str(sample_length) ' cm, ' num2str(sample_dilution) ...
        'x dilution in D2O. '];
    ytext = 'Optical Density (cm^{-1})';
%     ytext = 'Optical Density';
    
end
data_plot = data_OD;

%% Savitzky-Golay filtering
% *************************************************************************
data_filter = data_OD;
if find(strcmp(analysis_options(selected_analysis), 'Savitzky-Golay Filtering'))
    % if polynomial_order = 1 then this becomes a moving average
    polynomial_order = 1; % must be less than the frame size
    frame_size = 21; % must be odd
    
    input_title = 'Savitzky-Golay Filtering Parameters'; 
    input_data = {'Polynomial Order (odd) (1 = moving average):',...
        'Frame Size (odd):'};
    dlg_options.WindowStyle = 'normal'; dlg_options.Resize = 'on'; dim = [1 80];
    default_values = {num2str(polynomial_order),num2str(frame_size)};
    answer = inputdlg(input_data, input_title, dim, default_values, dlg_options);
    polynomial_order = str2double(answer{1});
    frame_size = str2double(answer{2});          
    
    title_analysis = [title_analysis 'SG filter: '...
        'order = ' num2str(polynomial_order) ', '...
        num2str(frame_size) ' points.'];
    
    for i = 1:1:number_of_spectra
        data_filter(i,:) = sgolayfilt(data_OD(i,:),polynomial_order,frame_size);
    end
end
data_plot = data_filter;

%% normalisation
% *************************************************************************
data_norm = data_filter;
if any(strcmp(analysis_options(selected_analysis), 'Normalisation'))    
    for i = 1:1:number_of_spectra
        data_norm(i,:) = data_filter(i,:) - min(data_filter(i,50:end-50));
        data_norm(i,:) = data_norm(i,:) ./ max(data_norm(i,:));
    end
    if any(strcmp(analysis_options(selected_analysis), 'Optical Density'))
        ytext = 'Normalised Optical Density (a.u.)';
    else
        ytext = 'Normalised Transmission (a.u.)';
    end    
end
data_plot = data_norm;

%% plot selected corrected data -----
% *************************************************************************

selected_index = 1:1:max(size(corrected_spectra));
[selected_index, ~] = listdlg('PromptString', 'Spectra to plot',...
                               'SelectionMode', 'multiple', ...
                               'ListString', DatasetNames(corrected_spectra),...
                               'InitialValue', selected_index);
selected_spectra = corrected_spectra(selected_index);

if or(...
        max(size(selected_analysis)) > 0, ...
        max(size(selected_spectra)) < number_of_spectra)
    
    new_selected_figure = new_figures;
    if exist('figure_selected','var') == 0
        new_selected_figure = 2;
    else
        if ishandle(figure_selected) == 0
            new_selected_figure = 2;
        end
    end
    
    if new_selected_figure == 2
        figure_selected = figure('Units','normalized','Position',[0.23 0.085 0.75 0.8]);
        figure_handles{end+1} = figure_selected;
        legend_selected = {};
        h_selected = {};
    end
    figure(figure_selected)
        
%     for i = selected_spectra
    for i = selected_spectra(end:-1:1)
        h_selected{end+1} = plot(wavelengths(i,:), data_plot(i,:), ...
            'LineWidth', 2); hold on
        legend_selected{end+1} = [DatasetNames{i} ' - ' ...
            description{i}];
    end

    grid on
    set(gca, 'FontSize', plot_font_size)
    legend(legend_selected, 'Location', legend_location, ...
        'FontSize', legend_font_size, 'interpreter', 'none')
    title({FolderRead, title_analysis}, 'interpreter', 'none')
    xlabel('Wavelength (nm)')
    ylabel(ytext)
    xlim([300,1050])
    if any(strcmp(analysis_options(selected_analysis), 'Optical Density'))
        ylim([-0.3,1.2]) % optical density
    end

% plot styling
% *************************************************************************
    [selected_colour, ~] = listdlg('PromptString', 'Colour scheme:',...
                               'SelectionMode', 'single', ...
                               'ListString', colour_type,...
                               'InitialValue', selected_colour);

    for i = 1:1:max(size(h_selected))
        if selected_colour > 1 
            colour_RGB = colour_gradient(i, max(size(h_selected)), colour_type(selected_colour));
            h_selected{i}.Color = colour_RGB;  
        end
        h_selected{i}.MarkerSize = 1;
        h_selected{i}.LineStyle = '-';
        h_selected{i}.LineWidth = 2;
    end
    
% read manufacturer optical density
% *************************************************************************
    if any(strcmp(analysis_options(selected_analysis), 'Optical Density'))
    %     directory = 'R:\aa938\NanoPhotonics\Laboratory\2016.11.16 - PBG HCF prep for Au NR\';
        directory = FolderPathRead;
        OD_file_name = 'nanocomposix-spectra.csv';
        % [OD_file_name, directory, ~] = uigetfile('.csv',...
        %                                       'OPTICAL DENSITY data', ...
        %                                       [directory OD_file_name], ...
        %                                       'MultiSelect','off');

        manufacturer_OD = dlmread([directory OD_file_name], ',', 0, 0);
        plot(manufacturer_OD(:,1), manufacturer_OD(:,2), 'k', 'LineWidth', 2)
        legend_selected{end+1} = OD_file_name(1:end-4);
        legend(legend_selected, 'Location', legend_location, ...
            'FontSize', legend_font_size, 'interpreter', 'none')
    end
end

      
%% saving figures
% *************************************************************************
menu_save_figures = 1;
menu_save_figures = menu('Save Figures?', 'NO', 'YES');

if menu_save_figures == 2    
    for i = 1:1:max(size(figure_handles))
        if ishandle(figure_handles{i})
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