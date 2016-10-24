clc
clear
close all

% default directory
DirectoryRead = 'C:\Users\Ana Andres\Documents\NanoPhotonics\Laboratory\';
FolderRead = '';
FolderPathRead = [DirectoryRead FolderRead];
FileNameRead = '';

% prompt to choose file
[FileNameRead, FolderPathRead, ~] = uigetfile('.h5',...
    'H5 file to read:',[FolderPathRead FileNameRead],'MultiSelect','off');
%     'H5 file to read:','MultiSelect','off');
slash_index = strfind(FolderPathRead, '\');
FolderRead = FolderPathRead(slash_index(end-1)+1:end-1);
FilePathRead = [FolderPathRead FileNameRead];

% selecting a timelapse (custom function)
TimelapseInfo = h5_TimelapseInfo(FilePathRead);
number_of_spectra = size(TimelapseInfo.Datasets,1);
number_of_wavelengths = TimelapseInfo.Datasets(1).Dataspace.Size;

% reading the data
data = zeros(number_of_spectra, number_of_wavelengths);
for i = 1:1:number_of_spectra
    DatasetName = [TimelapseInfo.Name '/' TimelapseInfo.Datasets(i).Name];
    data(i,:) = h5read(TimelapseInfo.Filename, DatasetName);
end

% sort the data array by time
index_sort = zeros(number_of_spectra,1);
for i = 1:1:number_of_spectra
    index_sort(i) = str2double(TimelapseInfo.Datasets(i).Name(10:end))+1;
end
data_sort(index_sort,:) = data;

% information (use the last dataset, assumed they all have the same parameters)
info_string = h5readatt(TimelapseInfo.Filename, DatasetName, 'information');

wavelengths = h5readatt(TimelapseInfo.Filename, DatasetName, 'wavelengths');
wavelength_start = round(wavelengths(1)); % nm
wavelength_end = round(wavelengths(end)); % nm


%% Timestamp

integration_time = h5readatt(TimelapseInfo.Filename, DatasetName, 'integration_time');
% time_interval = h5readatt(TimelapseInfo.Filename, DatasetName, 'time_interval');
% time_interval = 1;
% info_string = 'Info';
% time = 0:1:number_of_spectra-1;
% time = time * time_interval;
% time_start = time(1); % seconds
% time_end = time(end); % seconds
% time_interval_plot = time_interval; % seconds


creation_timestamp = cell(number_of_spectra,1);
% time_datetime = zeros(number_of_spectra,1);
time_unsorted = zeros(number_of_spectra,1);
time = zeros(number_of_spectra,1);
for i = 1:1:number_of_spectra
    DatasetName = [TimelapseInfo.Name '/' TimelapseInfo.Datasets(i).Name];
    creation_timestamp{i} = h5readatt(TimelapseInfo.Filename, DatasetName, 'creation_timestamp');
    time_datetime(i) = datetime(creation_timestamp{i}, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSSSSS');
    time_unsorted(i) = etime(datevec(time_datetime(i)),datevec(time_datetime(1)));
end
time(index_sort) = time_unsorted;
time = time';
time_interval = time(2);
time_interval_plot = time(2);
time_start = time(1);
time_end = time(end);


%% python analysis
% *************************************************************************
menu_python_analysis = 1; 
try
    reference = h5readatt(TimelapseInfo.Filename, DatasetName, 'reference');
    background = h5readatt(TimelapseInfo.Filename, DatasetName, 'background');
    menu_python_analysis = menu('Python Recorded Reference & Background', ...
        'IGNORE', 'APPLY just BG', 'APPLY just REF', 'APPLY BOTH'); 
catch
end

data_python = data_sort;
if menu_python_analysis == 2 % just BG
    for i = 1:1:number_of_spectra
        data_python(i,:) = data_sort(i,:) - background';
    end
    title_colourmap{5} = 'Python Background Correction';
elseif menu_python_analysis == 3 % just REF
    for i = 1:1:number_of_spectra
        data_python(i,:) = data_sort(i,:) ./ reference';
    end
    title_colourmap{5} = 'Python Reference Correction';
elseif menu_python_analysis == 4 % BG & REF
    for i = 1:1:number_of_spectra
        data_python(i,:) = (data_sort(i,:) - background') ./ (reference' - background');
    end
    title_colourmap{5} = 'Python Background and Reference Correction';
else
    title_colourmap{5} = 'No Correction in Python';
end



%% wavelength and time selection
% *************************************************************************

input_title = 'Plotting Range Selection'; 
input_data = {['Start Wavelength (nm) (spectrometer starts at ' num2str(wavelengths(1),'%.0f') ' nm) :'],...
    ['End Wavelength (nm) (spectrometer ends at ' num2str(wavelengths(end),'%.0f') ' nm) :'],...
    'Start Time (s):',...
    ['End Time (s) (data recorded for ' num2str(time(end)) ' s) :'], ...
    ['Time Interval for Plotting (s) (multiple of ' num2str(time_interval) ' s):']};
resize = 'on'; dim = [1 80];
valdef = {num2str(wavelength_start),num2str(wavelength_end)...
    num2str(time_start),num2str(time_end),num2str(time_interval_plot)};
answer = inputdlg(input_data,input_title,dim,valdef,resize);
wavelength_start = str2double(answer{1});
wavelength_end = str2double(answer{2});
time_start = str2double(answer{3});
time_end = str2double(answer{4});
time_interval_plot = str2double(answer{5});

[index_wav_start,~] = find(wavelengths >= wavelength_start);
[index_wav_end,~] = find(wavelengths <= wavelength_end);
index_wavelengths = intersect(index_wav_start,index_wav_end);

[~,index_time_start] = find(time >= time_start);
[~,index_time_end] = find(time <= time_end);
index_time_all = intersect(index_time_start,index_time_end);
time_interval_index = round(time_interval_plot/time_interval);
index_time = index_time_all(1):time_interval_index:index_time_all(end);

% %% system response reference
% *************************************************************************

menu_ref = 1;
menu_ref = menu('Correct for the System Response?', 'NO', 'YES (data/reference)', ...
    'YES (data/1st spectra)', 'YES ((data-background)/(reference-background))');

data_corrected = data_python;

if menu_ref == 2 || menu_ref == 4

    % prompt to choose reference file
    [FileNameReadRef, FolderPathReadRef, ~] = uigetfile('.h5',...
        'H5 file to read: REFERENCE:',[FolderPathRead FileNameRead],'MultiSelect','off');
    slash_index = strfind(FolderPathReadRef, '\');
    FilePathReadRef = [FolderPathReadRef FileNameReadRef];

    % selecting a reference timelapse (custom function)
    TimelapseInfoRef = h5_TimelapseInfo(FilePathReadRef);
    number_of_spectra_ref = size(TimelapseInfoRef.Datasets,1);
    number_of_wavelengths_ref = TimelapseInfoRef.Datasets(1).Dataspace.Size;

    % reading the reference data
    data_ref = zeros(number_of_spectra_ref, number_of_wavelengths_ref);
    for i = 1:1:number_of_spectra_ref
        DatasetNameRef = [TimelapseInfoRef.Name '/' TimelapseInfoRef.Datasets(i).Name];
        data_ref(i,:) = h5read(TimelapseInfoRef.Filename, DatasetNameRef);
    end

    % reference information 
    % (use the last dataset, assume they all have the same parameters)
    % assume same wavelengths as data file
    integration_time_ref = h5readatt(TimelapseInfoRef.Filename, DatasetNameRef, 'integration_time');
    info_string_ref = h5readatt(TimelapseInfoRef.Filename, DatasetNameRef, 'information');
    time_interval_ref = h5readatt(TimelapseInfoRef.Filename, DatasetNameRef, 'time_interval');
    time_ref = 0:1:number_of_spectra_ref-1;
    time_ref = time_ref * time_interval_ref;

    % sort the reference data array by time
    index_sort_ref = zeros(number_of_spectra_ref,1);
    for i = 1:1:number_of_spectra_ref
        index_sort_ref(i) = str2double(TimelapseInfoRef.Datasets(i).Name(10:end))+1;
    end
    data_sort_ref(index_sort_ref,:) = data_ref;
    
    time_start_ref = 0; % seconds
    time_end_ref = 900; % seconds
    time_interval_plot_ref = time_interval_ref; % seconds

    input_title = 'Reference Selection'; 
    input_data = {'Start Time (s):',...
        ['End Time (s) (data recorded for ' num2str(time_ref(end)) ' s) :'], ...
        ['Time Interval for Plotting (s) (multiple of ' num2str(time_interval_ref) ' s):']};
    resize = 'on'; dim = [1 80];
    valdef = {num2str(time_start_ref),num2str(time_end_ref),num2str(time_interval_plot_ref)};
    answer = inputdlg(input_data,input_title,dim,valdef,resize);
    time_start_ref = str2double(answer{1});
    time_end_ref = str2double(answer{2});
    time_interval_plot_ref = str2double(answer{3});

    [~,index_time_start_ref] = find(time_ref >= time_start_ref);
    [~,index_time_end_ref] = find(time_ref <= time_end_ref);
    index_time_all_ref = intersect(index_time_start_ref,index_time_end_ref);
    time_interval_index_ref = round(time_interval_plot_ref/time_interval_ref);
    index_time_ref = index_time_all_ref(1):time_interval_index_ref:index_time_all_ref(end);

    data_ref_mean = mean(data_sort_ref(index_time_ref,:),1);
    
    if menu_ref == 2
        for i = 1:1:number_of_spectra
            data_corrected(i,:) = data_sort(i,:)./data_ref_mean;
        end
        
    elseif menu_ref == 4
        % prompt to choose background file
        [FileNameReadBG, FolderPathReadBG, ~] = uigetfile('.h5',...
            'H5 file to read: BACKGROUND:',[FolderPathReadRef FileNameReadRef],'MultiSelect','off');
        slash_index = strfind(FolderPathReadBG, '\');
        FilePathReadBG = [FolderPathReadBG FileNameReadBG];

        % selecting a background timelapse (custom function)
        TimelapseInfoBG = h5_TimelapseInfo(FilePathReadBG);
        number_of_spectra_bg = size(TimelapseInfoBG.Datasets,1);
        number_of_wavelengths_bg = TimelapseInfoBG.Datasets(1).Dataspace.Size;

        % reading the background data
        data_bg = zeros(number_of_spectra_bg, number_of_wavelengths_bg);
        for i = 1:1:number_of_spectra_bg
            DatasetNameBG = [TimelapseInfoBG.Name '/' TimelapseInfoBG.Datasets(i).Name];
            data_bg(i,:) = h5read(TimelapseInfoBG.Filename, DatasetNameBG);
        end

        % background information 
        % (use the last dataset, assume they all have the same parameters)
        % assume same wavelengths as data file
        integration_time_bg = h5readatt(TimelapseInfoBG.Filename, DatasetNameBG, 'integration_time');
        info_string_bg = h5readatt(TimelapseInfoBG.Filename, DatasetNameBG, 'information');
        time_interval_bg = h5readatt(TimelapseInfoBG.Filename, DatasetNameBG, 'time_interval');
        time_bg = 0:1:number_of_spectra_bg-1;
        time_bg = time_bg * time_interval_bg;

        % sort the background data array by time
        index_sort_bg = zeros(number_of_spectra_bg,1);
        for i = 1:1:number_of_spectra_bg
            index_sort_bg(i) = str2double(TimelapseInfoBG.Datasets(i).Name(10:end))+1;
        end
        data_sort_bg(index_sort_bg,:) = data_bg;

        time_start_bg = 0; % seconds
        time_end_bg = 900; % seconds
        time_interval_plot_bg = time_interval_bg; % seconds

        input_title = 'Background Selection'; 
        input_data = {'Start Time (s):',...
            ['End Time (s) (data recorded for ' num2str(time_bg(end)) ' s) :'], ...
            ['Time Interval for Plotting (s) (multiple of ' num2str(time_interval_bg) ' s):']};
        resize = 'on'; dim = [1 80];
        valdef = {num2str(time_start_bg),num2str(time_end_bg),num2str(time_interval_plot_bg)};
        answer = inputdlg(input_data,input_title,dim,valdef,resize);
        time_start_bg = str2double(answer{1});
        time_end_bg = str2double(answer{2});
        time_interval_plot_bg = str2double(answer{3});

        [~,index_time_start_bg] = find(time_bg >= time_start_bg);
        [~,index_time_end_bg] = find(time_bg <= time_end_bg);
        index_time_all_bg = intersect(index_time_start_bg,index_time_end_bg);
        time_interval_index_bg = round(time_interval_plot_bg/time_interval_bg);
        index_time_bg = index_time_all_bg(1):time_interval_index_bg:index_time_all_bg(end);

        data_bg_mean = mean(data_sort_bg(index_time_bg,:),1);
        
        for i = 1:1:number_of_spectra
            data_corrected(i,:) = (data_sort(i,:)-data_bg_mean)./(data_ref_mean-data_bg_mean);
        end
    end
    
    title_colourmap{6} = ['REFERENCE: ' FileNameReadRef ' // ' ...
        strrep(TimelapseInfo.Name(2:end), '_', '\_') ' // ' info_string_ref];
    title_colourmap{7} = ['Start t = ' num2str(time_ref(index_time_ref(1)))...
        ' s // End t = ' num2str(time_ref(index_time_ref(end)))...
        ' s // Delta t = ' num2str(time_interval_index_ref*time_interval_ref) ' s.'];
    if menu_ref == 4
        title_colourmap{8} = ['BACKGROUND: ' FileNameReadBG ' // ' ...
            strrep(TimelapseInfo.Name(2:end), '_', '\_') ' // ' info_string_bg];
        title_colourmap{9} = ['Start t = ' num2str(time_bg(index_time_bg(1)))...
            ' s // End t = ' num2str(time_bg(index_time_bg(end)))...
            ' s // Delta t = ' num2str(time_interval_index_bg*time_interval_bg) ' s.'];
    end

elseif menu_ref == 3
    for i = 1:1:number_of_spectra
        data_corrected(i,:) = data_python(i,:)./data_python(1,:);
    end
    title_colourmap{6} = 'REFERENCE: same timelapse at t = 0s';
end


%% transmission / absorbance / molar attenuation coefficient
% *************************************************************************
data_spectra = data_corrected;
string_spectra = 'Transmission (a.u.)';

if menu_ref || menu_python == 2   % reference corrected
    menu_spectra = 1;
    menu_spectra = menu('Spectra', 'Transmission', 'Absorbance (-log(data))', 'Absorbance/Sample Length', 'Molar Attenuation Coefficient');
    if menu_spectra == 2 % absorbance
        for i = 1:1:number_of_spectra
            data_spectra(i,:) = -log(data_corrected(i,:));
        end
        string_spectra = 'Absorbance';
    end
     if menu_spectra == 3 % Absorbance/Sample Length
      input_title = 'Choose Sample length'; 
    input_data = {'Sample length (cm)'}; 
    resize = 'on'; dim = [1 80];
    valdef = {num2str(4.8)};
    answer = inputdlg(input_data,input_title,dim,valdef,resize);
        length = str2double(answer{1});
        for i = 1:1:number_of_spectra
            data_spectra(i,:) = (-log(data_corrected(i,:))./ length);
        end
        string_spectra = 'Absorbance / Sample Length (cm^{-1})';
    end
    if menu_spectra == 4 % molar attenuation coefficient
        input_title = 'Choose Concentration and Sample length'; 
        input_data = {'Concentration of MV (mol/L)', 'Sample length (cm)'}; 
        resize = 'on'; dim = [1 80];
        valdef = {num2str(40E-6), num2str(4.8)};
        answer = inputdlg(input_data,input_title,dim,valdef,resize);
        concentration = str2double(answer{1});
        length = str2double(answer{2});
        for i = 1:1:number_of_spectra
            data_spectra(i,:) = (-log(data_corrected(i,:))./(concentration * length));
        end
        string_spectra = 'Molar Attenuation Coefficient (M^{-1}cm^{-1})';
    end
end

%% linear / log
% *************************************************************************

menu_scale = 1;
% menu_scale = menu('Intensity Scale', 'Linear', 'Logarithmic');
data_scale = data_spectra;
string_scale = 'Linear Scale';
if menu_scale == 2
    data_scale = log(data_spectra);
    string_scale = 'Logarithmic Scale';
end

%% normalise spectra
% *************************************************************************

menu_norm = 1;
menu_norm = menu('Normalise Timelapse?', 'NO', 'YES (all spectra)', 'YES (each spectra)');
data_norm = data_scale;
% offset = 2500;
offset = 0;
if menu_norm == 2
    for i = 1:1:number_of_spectra
        data_norm(i,:) = (data_scale(i,:)-offset) / (max(max(data_scale))-offset);
    end
elseif menu_norm == 3
    for i = 1:1:number_of_spectra
        data_norm(i,:) = (data_scale(i,:)-offset) / (max(data_scale(i,:))-offset);
    end
end
data_plot = data_norm;


%% plotting 2D colourmap
% *************************************************************************

menu_plot2D = 1;
% menu_plot2D = menu('Plot 2D Colourmap?', 'NO', 'YES');

title_colourmap{1} = FolderRead;
title_colourmap{2} = ['DATA: ' FileNameRead ' // ' strrep(TimelapseInfo.Name(2:end), '_', '\_')];
% title_colourmap{3} = [info_string ' // ' string_spectra ' // ' string_scale];
title_colourmap{3} = [info_string ' // ' string_spectra];
title_colourmap{4} = ['Start t = ' num2str(time(index_time(1)))...
    ' s // End t = ' num2str(time(index_time(end)))...
    ' s // Delta t = ' num2str(time_interval_index*time_interval) ' s.'];

if menu_plot2D == 2
    fig_timelapse = figure('Units','normalized','Position',[0.01 0.07 0.75 0.7], 'tag', 'fig_timelapse');
    nlevels = linspace(min(min(data_plot(index_time,index_wavelengths))),...
        max(max(data_plot(index_time,index_wavelengths))), 200);
    % nlevels = linspace(0, max(max(data_plot(index_time,index_wavelengths))), 200);
    contourf(wavelengths(index_wavelengths), ...
        time(index_time),...
        data_plot(index_time,index_wavelengths), ...
        'LineStyle', 'none');%,...
%         'LevelStepMode', 'manual', 'LevelStep', 10)
%         'LevelListMode', 'manual', 'LevelList', nlevels)
    colormap(jet)
    colorbar
    xlabel('Wavelength (nm)')
    ylabel('Time (s)')
    set(gca, 'FontSize', 12)
    set(gca, 'FontSize', 12)
    title(title_colourmap)
end


%% plotting individual spectra
% *************************************************************************

menu_evolution = 1;
menu_evolution = menu('Plot Individual Spectra?', 'NO', 'YES');

if menu_evolution == 2
menu_irregular_interval = menu('Plot spectra with irregular interval? ','NO','YES');


if menu_irregular_interval == 1
    gradient_type = {'DEFAULT','yellow' 'red', 'green', 'aqua', 'blue', 'purple', 'gray'};
    menu_gradient = 1; 
    menu_gradient = menu('Which colour scheme to use?', gradient_type);

    fig_spectra = figure('Units','normalized','Position',[0.23 0.17 0.75 0.7], 'tag', 'fig_spectra');
    legn_individual = {};
    total_time_plot = size(index_time,2);
    for t = 1:1:total_time_plot
        if menu_gradient == 1 % default
            plot(wavelengths, data_plot(index_time(t),:), 'LineWidth', 1); hold on
        else % colour gradient
            colour_RGB = colour_gradient(t,total_time_plot, gradient_type(menu_gradient));
            colour_RGB = colour_gradient(total_time_plot-t,total_time_plot, gradient_type(menu_gradient));
            plot(wavelengths, data_plot(index_time(t),:), 'LineWidth', 1, ...
                'color', colour_RGB); hold on
        end
        legn_individual{end+1} = ['t = ' num2str(time(index_time(t))) ' s'];
    end
    xlabel('Wavelength (nm)')
    ylabel(string_spectra)
    xlim([wavelengths(index_wavelengths(1)), wavelengths(index_wavelengths(end))])
    if size(legn_individual,2)<30
        legend(legn_individual)
    end
    set(gca, 'FontSize', 12)
    title(title_colourmap)
end

if menu_irregular_interval == 2
    input_title = 'Choose how many spectra to plot';
    input_data = {'Number of spectra'};
    resize = 'on'; dim = [1 80];
    valdef = {'1'};
    answer = inputdlg(input_data,input_title,dim,valdef,resize);
    num_of_spectra = str2double(answer{1});
    clear valdef
    clear input_data
    input_title = 'Select times for spectra plotting'; 
    valdef = cell(1, num_of_spectra);
    input_data = cell(1, num_of_spectra);
    for i = 1:num_of_spectra
      input_data{i} = ['Time for Spectrum ' num2str(i) ' (s)'];  
      valdef{i} = '1';
    end
    valdef{1} = '0';
    valdef{2} = '60';
    valdef{3} = '300';
    valdef{4} = '600';
    valdef{5} = '1200';
    valdef{6} = '1800';
    
    resize = 'on'; dim = [1 80];
    answer = inputdlg(input_data,input_title,dim, valdef,resize);
    
    select_time = zeros(1, num_of_spectra);
    select_index = zeros(1, num_of_spectra);
    select_legend = cell(1, num_of_spectra);
    for i = 1:num_of_spectra
        select_time(i) = str2double(answer{i});
        [~, select_index(i)] = min(abs(time - select_time(i)));
        
        plot(wavelengths, data_plot(select_index(i),:),'linewidth',2);
        hold all
        select_legend{i} = [num2str(time(select_index(i))) ' s'];
    end
    xlabel('Wavelength (nm)');
    xlim([450 800])
    ylabel(string_spectra);
    legend(select_legend);
    set(gca, 'FontSize', 16);
    title(title_colourmap);
    
end

end



%% plotting trace
% *************************************************************************

menu_trace = 1;
menu_trace = menu('Plot Trace?', 'NO', 'YES');

if menu_trace == 2
    wavelength_start_trace = 603; % nm
    wavelength_end_trace = 604; % nm
    time_start_trace = time_start; % seconds
    time_end_trace = time_end; % seconds
    time_interval_plot_trace = time_interval_plot; % seconds

    input_title = 'Trace Range Selection'; 
    input_data = {['Start Wavelength (nm) (spectrometer starts at ' num2str(wavelengths(1),'%.0f') ' nm) :'],...
        ['End Wavelength (nm) (spectrometer ends at ' num2str(wavelengths(end),'%.0f') ' nm) :'],...
        'Start Time (s):',...
        ['End Time (s) (data recorded for ' num2str(time(end)) ' s) :'], ...
        ['Time Interval for Trace (s) (multiple of ' num2str(time_interval) ' s):']};
    resize = 'on'; dim = [1 80];
    valdef = {num2str(wavelength_start_trace),num2str(wavelength_end_trace)...
        num2str(time_start_trace),num2str(time_end_trace),num2str(time_interval_plot_trace)};
    answer = inputdlg(input_data,input_title,dim,valdef,resize);
    wavelength_start_trace = str2double(answer{1});
    wavelength_end_trace = str2double(answer{2});
    time_start_trace = str2double(answer{3});
    time_end_trace = str2double(answer{4});
    time_interval_plot_trace = str2double(answer{5});

    [index_wav_start_trace,~] = find(wavelengths >= wavelength_start_trace);
    [index_wav_end_trace,~] = find(wavelengths <= wavelength_end_trace);
    index_wavelengths_trace = intersect(index_wav_start_trace,index_wav_end_trace);

    [~,index_time_start_trace] = find(time >= time_start_trace);
    [~,index_time_end_trace] = find(time <= time_end_trace);
    index_time_all_trace = intersect(index_time_start_trace,index_time_end_trace);
    time_interval_index_trace = round(time_interval_plot_trace/time_interval);
    index_time_trace = index_time_all_trace(1):time_interval_index_trace:index_time_all_trace(end);
    
    data_plot_mean = mean(data_plot(:,index_wavelengths_trace),2);
    
    % smoothing the trace
    menu_smoothing = 1;
    menu_smoothing = menu('Smooth Data?', 'NO', 'YES: Savitzky-Golay filtering');

    polynomial_order = 1; % must be less than the frame size
    % if polynomial_order = 1 then this becomes a moving average
    frame_size = 51; % must be odd
    
    if menu_smoothing == 2 % Savitzky-Golay filtering
        input_data = {'Polynomial order (must be less than the frame size)', 'Frame size (must be odd)'};
            input_title = 'Parameters for Savitzky-Golay filtering'; 
            input_resize = 'on'; input_dim = [1 60];
            input_valdef = {num2str(polynomial_order), num2str(frame_size)};
            input_answer = inputdlg(input_data,input_title,input_dim,input_valdef,input_resize);
            polynomial_order = str2double(input_answer{1});
            frame_size = str2double(input_answer{2});

        title_cell{2} = ['Savitzky-Golay filtering,'...
            ' Polynomial order = ' num2str(polynomial_order),','...
            ' Frame size =  ' num2str(frame_size) ' points.'];

            data_plot_mean = sgolayfilt(data_plot_mean,polynomial_order,frame_size);
    end
    
    title_trace = title_colourmap;
%     title_trace = title_colourmap(1:3);
%     title_trace{4} = ['Start t = ' num2str(time(index_time_trace(1)))...
%         ' s // End t = ' num2str(time(index_time_trace(end)))...
%         ' s // Delta t = ' num2str(time_interval_index_trace*time_interval) ' s.'];
    
    menu_new_trace = 1;
    menu_new_trace = menu('Plot Trace In:', 'New Plot', 'Existing Plot');
    if menu_new_trace == 1
        fig_trace = figure('Units','normalized','Position',[0.12 0.12 0.75 0.7], 'tag', 'fig_trace');
        legn_trace = {};
    elseif menu_new_trace == 2
        figure(fig_trace)
    end
    plot(time(index_time_trace), data_plot_mean(index_time_trace), ...
        'LineWidth', 2), hold all
    legn_trace{end+1} = ['\lambda = ' num2str(wavelength_start_trace) ' to ' ...
        num2str(wavelength_end_trace) ' nm'];
    if menu_smoothing == 2 % smoothing
        legn_trace{end} = [legn_trace{end},', ' title_cell{2}]; %' // Polynomial order = ' num2str(polynomial_order), ' // Frame size = ' num2str(frame_size),' points'];
    end
    
    legend(legn_trace, 'Location', 'SE')
    xlabel('Time (s)')
    ylabel(['Averaged ' string_spectra])
    set(gca, 'FontSize', 16)
    title(title_trace)
    y = ylim;
end



%% fitting trace, saturation function
% *************************************************************************
menu_fit = 1;
if menu_trace == 2
    menu_fit = menu('Fit Trace, Saturation?', 'NO', 'YES');
end

if menu_fit == 2
    time_start_fit = time_start; % seconds
    time_end_fit = time_end; % seconds

    input_title = 'Time Selection for Curve Fitting'; 
    input_data = {'Start Time (s):',...
        ['End Time (s) (data recorded for ' num2str(time(end)) ' s) :']};
    resize = 'on'; dim = [1 80];
    valdef = {num2str(time_start_fit),num2str(time_end_fit)};
    answer = inputdlg(input_data,input_title,dim,valdef,resize);
    time_start_fit = str2double(answer{1});
    time_end_fit = str2double(answer{2});


    [~,index_time_start_fit] = find(time >= time_start_fit);
    [~,index_time_end_fit] = find(time <= time_end_fit);
    index_time_all_fit = intersect(index_time_start_fit,index_time_end_fit);
    time_interval_index_trace = round(time_interval_plot_trace/time_interval);
    index_time_trace = index_time_all_fit(1):time_interval_index_trace:index_time_all_fit(end);
    
%     data_plot_mean = mean(data_plot(:,index_wavelengths_trace),2);
    
    plot(time(index_time_trace), data_plot_mean(index_time_trace), ...
        'LineWidth', 2), hold all
    legn_trace{end+1} = ['\lambda = ' num2str(wavelength_start_trace) ' to ' ...
        num2str(wavelength_end_trace) ' nm, t = ' num2str(time_start_fit) ' to ' ...
        num2str(time_end_fit) ' s'];
    
    legend(legn_trace, 'Location', 'SE')
    xlabel('Time (s)')
    ylabel(['Averaged ' string_spectra])
    set(gca, 'FontSize', 12)
    title(title_trace)
    y = ylim;   
    
    value_a = 2.5;
%     value_b = 1.4;
    value_tau = 100;    
    value_to = 100;
    
    input_title = 'Choose Starting Parameters for a*(1 - exp(-(x-t_0)/tau))';
    input_data = {'a','tau (s)', 't_0 (s)'}; % need to be in alphabetical order!!
    resize = 'on'; dim = [1 120];
    valdef = {num2str(value_a), num2str(value_tau), num2str(value_to)};
    answer = inputdlg(input_data,input_title,dim,valdef,resize);
    value_a = str2double(answer{1});
    value_tau = str2double(answer{2});
    value_to = str2double(answer{3});
    
    Exponential = fittype('a*(1 - exp(-((x-to)/tau)))');
    Exponential_Start = [value_a, value_tau, value_to]; % need to be in alphabetical order!!
    [fit_trace,gof_trace] = fit(time(index_time_trace)', data_plot_mean(index_time_trace),...
        Exponential, 'Startpoint', Exponential_Start);
    p = plot(fit_trace, '--k');
    p.LineWidth = 1.5;
    
    ylim(y)
    xlabel('Time (s)')
    ylabel(['Averaged ' string_spectra])
    legend('Location', 'SE')
    
    line([time(index_time_trace(1)), time(index_time_trace(1))], y, ...
        'Color', 'k', 'LineWidth', 1.5)
    line([time(index_time_trace(end)), time(index_time_trace(end))], y, ...
        'Color', 'k', 'LineWidth', 1.5)

    fit_trace_confint = confint(fit_trace);
    fit_trace_delta_a = abs(fit_trace_confint(1,1) - fit_trace_confint(2,1))/2;
    fit_trace_delta_tau = abs(fit_trace_confint(1,2) - fit_trace_confint(2,2))/2;
    fit_trace_delta_to = abs(fit_trace_confint(1,3) - fit_trace_confint(2,3))/2;

    text_fit{1} = 'Exponential Fit: a*(1 - exp(-(x-t_0)/\tau))';
    text_fit{2} = ['a = ' num2str(fit_trace.a, '%.2f') ' \pm ' ...
        num2str(fit_trace_delta_a, '%.2f')];
    text_fit{3} = ['\tau = ' num2str(fit_trace.tau, '%.0f') ' \pm ' ...
        num2str(fit_trace_delta_tau, '%.0f') ' s'];
    text_fit{4} = ['t_0 = ' num2str(fit_trace.to, '%.0f') ' \pm ' ...
        num2str(fit_trace_delta_to, '%.0f') ' s'];
    text('Units','normalized','Position',[0.1,0.9],'VerticalAlignment','top','String',text_fit) 

    nameSave = 'Exp_1_TauAndT0.dat';
    pause(0.1)
    [nameSave,dirSave,~] = uiputfile(['.' 'dat'],...
        'New File to save values of the exponential fit',...
        [FolderPathRead nameSave]); % choosing the file name
    fid3 = fopen([dirSave nameSave], 'wt');
    fprintf(fid3, 'Tau = %.2f +- %.2f (s)\n',...
        fit_trace.tau, fit_trace_delta_tau);
    fprintf(fid3, 't_0 = %.2f +- %.2f (s)\n',...
       fit_trace.to, fit_trace_delta_to);
    fprintf(fid3, 'a = %.2f +- %.2f\n',...
        fit_trace.a, fit_trace_delta_a);
	fclose(fid3);
end
   
% 
%% fitting trace, decay function
% *************************************************************************
menu_fit = 1;
if menu_trace == 2
    menu_fit = menu('Fit Trace, Decay?', 'NO', 'YES');
end
if menu_fit == 2
    time_start_fit = 300; % seconds
    time_end_fit = 900; % seconds

    input_title = 'Time Selection for Curve Fitting'; 
    input_data = {'Start Time (s):',...
        ['End Time (s) (data recorded for ' num2str(time(end)) ' s) :']};
    resize = 'on'; dim = [1 80];
    valdef = {num2str(time_start_fit),num2str(time_end_fit)};
    answer = inputdlg(input_data,input_title,dim,valdef,resize);
    time_start_fit = str2double(answer{1});
    time_end_fit = str2double(answer{2});


    [~,index_time_start_fit] = find(time >= time_start_fit);
    [~,index_time_end_fit] = find(time <= time_end_fit);
    index_time_all_fit = intersect(index_time_start_fit,index_time_end_fit);
    time_interval_index_trace = round(time_interval_plot_trace/time_interval);
    index_time_trace = index_time_all_fit(1):time_interval_index_trace:index_time_all_fit(end);
    
%     data_plot_mean = mean(data_plot(:,index_wavelengths_trace),2);
    
    plot(time(index_time_trace), data_plot_mean(index_time_trace), ...
        'LineWidth', 2), hold all
    legn_trace{end+1} = ['\lambda = ' num2str(wavelength_start_trace) ' to ' ...
        num2str(wavelength_end_trace) ' nm, t = ' num2str(time_start_fit) ' to ' ...
        num2str(time_end_fit) ' s'];
    
    legend(legn_trace, 'Location', 'SE')
    xlabel('Time (s)')
    ylabel(['Averaged ' string_spectra])
    set(gca, 'FontSize', 12)
    title(title_trace)
    y = ylim;   
    
    value_a = 1.5;
    value_tau = 100;    
    value_to = 300;
    value_b = 2.5;
    
    input_title = 'Choose Starting Parameters for a*(exp(-(x-t_0)/tau))';
    input_data = {'a','tau (s)', 't_0 (s)', 'b'};
    resize = 'on'; dim = [1 120];
    valdef = {num2str(value_a), num2str(value_tau), num2str(value_to), num2str(value_b)};
    answer = inputdlg(input_data,input_title,dim,valdef,resize);
    value_a = str2double(answer{1});
    value_tau = str2double(answer{2});
    value_to = str2double(answer{3});
    value_b = str2double(answer{4});
    
    Exponential = fittype('-a*(1-exp(-(x-to)/tau))+ b');
    Exponential_Start = [value_a, value_b, value_tau, value_to];
    [fit_trace,gof_trace] = fit(time(index_time_trace)', data_plot_mean(index_time_trace),...
        Exponential, 'Startpoint', Exponential_Start, 'Lower', [0, 0, 0, time_start_fit], 'Upper', [Inf, Inf, Inf, time_start_fit])
    p = plot(fit_trace, '--k');
    p.LineWidth = 1.5;
    
    ylim(y)
    xlabel('Time (s)')
    ylabel(['Averaged ' string_spectra])
    legend('Location', 'SE')
    
    line([time(index_time_trace(1)), time(index_time_trace(1))], y, ...
        'Color', 'k', 'LineWidth', 1.5)
    line([time(index_time_trace(end)), time(index_time_trace(end))], y, ...
        'Color', 'k', 'LineWidth', 1.5)

    fit_trace_confint = confint(fit_trace);
    fit_trace_delta_a = abs(fit_trace_confint(1,1) - fit_trace_confint(2,1))/2;
    fit_trace_delta_tau = abs(fit_trace_confint(1,3) - fit_trace_confint(2,3))/2;
    fit_trace_delta_b = abs(fit_trace_confint(1,2) - fit_trace_confint(2,2))/2;

    text_fit{1} = 'Exponential Fit: -a*(1-exp(-(x-to)/tau))+ b';
    text_fit{2} = ['a = ' num2str(fit_trace.a, '%.2f') ' \pm ' ...
        num2str(fit_trace_delta_a, '%.2f')];
    text_fit{3} = ['\tau = ' num2str(fit_trace.tau, '%.0f') ' \pm ' ...
        num2str(fit_trace_delta_tau, '%.0f') ' s'];
    text_fit{4} = ['t_0 = ' num2str(time_start_fit, '%.0f'),' s'];
    text_fit{5} = ['b = ' num2str(fit_trace.b, '%.1f') ' \pm ' ...
        num2str(fit_trace_delta_b, '%.1f') ' s'];
    text('Units','normalized','Position',[0.1,0.9],'VerticalAlignment','top','String',text_fit) 

        nameSave = 'Exp_1_TauAndT0.dat';
        pause(0.1)
        [nameSave,dirSave,~] = uiputfile(['.' 'dat'],...
            'New File to save values of the exponential fit',[FolderPathRead nameSave]); % choosing the file name
        
        fid3 = fopen([dirSave nameSave], 'wt');
        fprintf(fid3, 'Tau = %.2f +- %.2f (s)\n',...
            fit_trace.tau, fit_trace_delta_tau);
        fprintf(fid3, 't_0 = %.2f +- %.2f (s)\n',...
           fit_trace.to, fit_trace_delta_to);
        fprintf(fid3, 'a = %.2f +- %.2f\n',...
            fit_trace.a, fit_trace_delta_a);
	fclose(fid3);
end
   

%% saving the figures
% *************************************************************************

menu_save_fig = 1;
menu_save_fig = menu('Save Figures?', 'NO', 'YES'); 
if menu_save_fig == 2    
    pathSave = FolderPathRead;
    nameSave = [strrep(FolderRead(1:10), '.', '-') '-' ...
            TimelapseInfo.Name(2:14) '-' TimelapseInfo.Name(end) '-timelapse'];

    if findobj('tag','fig_timelapse') ~= 0
        figure(fig_timelapse)
        pause(0.1)
        
        [nameSave,pathSave,~] = uiputfile(['.' 'png'],...
            'New File to Save the Figure',[pathSave nameSave]);
        hgexport(fig_timelapse, [pathSave nameSave], hgexport('factorystyle'), 'Format', 'png')
        nameSave = strrep(nameSave, 'png', 'fig');    
        saveas(fig_timelapse, [pathSave nameSave], 'fig');
        
    end
    
    if findobj('tag','fig_spectra') ~= 0
        figure(fig_spectra)
        pause(0.1)
        
        nameSave = strrep(nameSave, 'fig', 'png');
        nameSave = strrep(nameSave, 'timelapse', 'spectra');
        [nameSave,pathSave,~] = uiputfile(['.' 'png'],...
            'New File to Save the Figure',[pathSave nameSave]);
        hgexport(fig_spectra, [pathSave nameSave], hgexport('factorystyle'), 'Format', 'png')
        nameSave = strrep(nameSave, 'png', 'fig');    
        saveas(fig_spectra, [pathSave nameSave], 'fig');
        
    end
    
    if findobj('tag','fig_trace') ~= 0
        figure(fig_trace)
        pause(0.1)
        
        nameSave = strrep(nameSave, 'fig', 'png');
        nameSave = strrep(nameSave, 'timelapse', 'trace');
        nameSave = strrep(nameSave, 'spectra', 'trace');
        [nameSave,pathSave,~] = uiputfile(['.' 'png'],...
            'New File to Save the Figure',[pathSave nameSave]);
        hgexport(fig_trace, [pathSave nameSave], hgexport('factorystyle'), 'Format', 'png')
        nameSave = strrep(nameSave, 'png', 'fig');    
        saveas(fig_trace, [pathSave nameSave], 'fig');
        
    end
end
