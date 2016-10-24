% Created by Ana (aa938)
% Reads txt files created with the PicoScope software
% Savitzky-Golay filtering
% Plots the channels

clc
clear
close all

%% READING DATA
% *************************************************************************

% specify default path
file_path = 'R:\aa938\NanoPhotonics\Laboratory\';
% file_path = 'R:\3-Temporary\aa938\';

% pop up window to choose the file(s) to read
[file_name, file_path, ~] = uigetfile('.txt',...
    'PicoScope Files to Read (use CTRL to select multiple files)',file_path,'MultiSelect','on');
file_path_save = file_path;
if strfind(file_path, 'Laboratory')
    title_cell_channels{1} = file_path(strfind(file_path, 'Laboratory')+11:end-1);
elseif strfind(file_path, 'aa938')
    title_cell_channels{1} = file_path(strfind(file_path, 'aa938')+6:end-1);
else
    title_cell_channels{1} = file_path;
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
% it takes a while because the files are large
raw_data = cell(size(file_name));
for i = 1:1:number_of_files
    disp(['Reading File ' num2str(i) '/' num2str(number_of_files)])
    raw_data{i} = dlmread([file_path file_name{i}], '\t', 2, 0);
    
    % raw_data{i}(:,1) = time
    % raw_data{i}(:,2) = channel A
    % raw_data{i}(:,3) = channel B
    
    % read channel names and units
    file_identifier = fopen([file_path file_name{i}], 'r');
    line1 = fgets(file_identifier);
    line1 = strrep(line1, 'ChannelA', 'Channel A');
    line1 = strrep(line1, 'ChannelB', 'Channel B');
    line2 = fgets(file_identifier);
    fclose(file_identifier);
    channel_name(i,:) = strsplit(line1(1:end-1),'\t');
    channel_unit(i,:) = strsplit(line2(1:end-1),'\t');
    
    % change all channel units to Volts
    for j = 2:1:size(raw_data{i},2) % channels
        if strcmp(channel_unit{i,j}, '(V)')
        elseif strcmp(channel_unit{i,j}, '(mV)')
            channel_unit{i,j} = '(V)';
            raw_data{i}(:,j) = raw_data{i}(:,j)/1000;
        else
            disp(['Check Channel (' num2str(i) ',' num2str(j) ') Units!'])
        end
    end
end
disp('Finished reading all files!')

%% OPTIONS
% *************************************************************************
options = {'Savitzky-Golay Filtering', ...
           'Normalisation', ...
           'Average & Standard Dev.',...
           'Period Selection',...
           'Channel B / Channel A'};
selected_options = 5;
[selected_options, ~] = listdlg('PromptString', 'Select options:',...
                                'SelectionMode', 'multiple', ...
                                'ListString', options, ...
                                'InitialValue', selected_options, ...
                                'OKString', 'OK', ...
                                'CancelString', 'NONE');

%% SMOOTHING
% *************************************************************************

smoothed_data = raw_data;   
% title_cell_channels{2} = 'Raw unfiltered data';
if find(strcmp(options(selected_options), 'Savitzky-Golay Filtering'))
    % if polynomial_order = 1 then this becomes a moving average
    polynomial_order = 1; % must be less than the frame size
    frame_size = 21; % must be odd
    
    input_title = 'Savitzky-Golay Filtering Parameters'; 
    input_data = {'Polynomial Order (odd) (1 = moving average):',...
        'Frame Size (odd):'};
    resize = 'on'; dimensions = [1 80];
    default_values = {num2str(polynomial_order),num2str(frame_size)};
    answer = inputdlg(input_data, input_title, dimensions, default_values, resize);
    polynomial_order = str2double(answer{1});
    frame_size = str2double(answer{2});    
    
    
    title_cell_channels{2} = ['Savitzky-Golay filtering. '...
        'Polynomial order = ' num2str(polynomial_order) '. '...
        'Frame size =  ' num2str(frame_size) ' points.'];
    
    for i = 1:1:number_of_files % files
%         disp(['Smoothing File ' num2str(i) '/' num2str(number_of_files)])
        for j = 2:1:size(smoothed_data{i},2) % channels
            smoothed_data{i}(:,j) = sgolayfilt(raw_data{i}(:,j),polynomial_order,frame_size);
        end
    end
    disp('Finished smoothing all files!')
end

%% NORMALISATION
% *************************************************************************
normalised_data = smoothed_data;
if find(strcmp(options(selected_options), 'Normalisation'))
%     t = 200; % ms
    t = -20; % ms
%     t = -100; % ms
    for i = 1:1:number_of_files
        [~, t_index] = min(abs(smoothed_data{i}(:,1)-t));
        for j = 2:1:size(smoothed_data{i},2)
            normalised_data{i}(:,j) = smoothed_data{i}(:,j) / smoothed_data{i}(t_index,j);
        end
    end
end

%% STANDARD DEVIATION
% *************************************************************************
average = zeros(size(channel_name));
stdev = zeros(size(channel_name));
if find(strcmp(options(selected_options), 'Average & Standard Dev.'))
    for i = 1:1:number_of_files
        threshold_ref = max(period_data{i}(:,2))/2; % V
        [t_ref,~] = find(period_data{i}(:,2) > threshold_ref);
        t_ref(end-9:end) = [];
        t_ref(1:10) = [];    
        for j = 1:1:size(channel_name,2)
            average(i,j) = mean(period_data{i}(t_ref,1));
            stdev(i,j) = std(period_data{i}(t_ref,1));
        end
    end
end

%% PERIOD SELECTION
% *************************************************************************
period_data = normalised_data;
threshold_time = 5; % ms

if find(strcmp(options(selected_options), 'Period Selection'))
%     T = [18,40,18,40];
    T = ones(size(file_name));
    T = T * 10;
    for i = 1:1:number_of_files
        [~, T_index] = min(abs(period_data{i}(:,1)-T(i)));
        period_data{i}(T_index:end,:) = [];
        
%         threshold_ref = max(period_data{i}(:,2))/2; % V
%         [t_index,~] = find(period_data{i}(:,1) > threshold_time);
%         [~, T_index] = min(abs(period_data{i}(t_index,2)-threshold_ref));
%         T_index = t_index(1) + T_index - 40;
%         period_data{i}(T_index:end,:) = [];
%         T(i) = round(period_data{i}(end,1),-1);
    end
end

%% CALCULATING THE POWER
% *************************************************************************
waveplate_angle = zeros(size(file_name));
% waveplate_angle = [60, 55, 50, 50, 50, 45, 45, 40];
% waveplate_angle = [60, 55, 50, 50, 50, 45, 45, 40, 40];
waveplate_angle_unique = unique(waveplate_angle);
% Calibration parameters from 22/09/2016 @ 800 nm
a = 5.353; % W
b = 1.994; % 1/rad
c = 0.9775; % rad
laser_power = a*(sin(b*waveplate_angle*pi/180+c)).^2*1000; % mW


%% CHANNEL FACTORS
% *************************************************************************
channel_factor = ones(size(file_name));
measurement_factor = ones(size(file_name));
% measurement_factor = [1,1,1,1,1,1,1,1,1,1,1];
% measurement_factor = [1/1, 1/1, 1/0.9, 1/0.8, 1/0.7, 1/0.6, 1/0.5, 1/0.4, 1/0.3, 1/0.1,1,1];
% measurement_factor = [1.86,1,1,1,1,1,1,1,1,1,1];
% measurement_factor = [26,1,26,28,22,1,1,1,1,1,1];

ND_A = zeros(size(file_name));
ND_B = zeros(size(file_name));
% ND_A = [0.4, 0.4, 0.4, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6];
% ND_B = [0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2];
% ND_A = [4, 4, 4, 4, 4, 4, 4, 4, 4];
% ND_B = [2, 2, 2, 4, 4, 4, 4, 4, 4];

ND_A = ND_A + 0.01;
ND_B = ND_B + 0.01;

if find(strcmp(channel_name(1,:), 'Channel A'))
    input_title = 'ND Channel A'; 
    input_data = file_name;
    resize = 'on'; dimensions = [1 60];
    default_values = cellstr(num2str((ND_A')));
    answer = inputdlg(input_data, input_title, dimensions, default_values, resize);
    ND_A = cell2mat(answer);
end
if find(strcmp(channel_name(1,:), 'Channel B'))
    input_title = 'ND Channel B'; 
    input_data = file_name;
    resize = 'on'; dimensions = [1 60];
    default_values = cellstr(num2str((ND_B')));
    answer = inputdlg(input_data, input_title, dimensions, default_values, resize);
    ND_B = cell2mat(answer);
end
ND = [ND_A, ND_B];

channel_data = period_data;
for i = 1:1:number_of_files % files
    for j = 2:1:size(channel_name,2)-1 % channels
        channel_data{i}(:,j) = channel_data{i}(:,j)*channel_factor(j)*measurement_factor(i)*10^ND(i,j-1);
    end
end
plot_data = channel_data;

%% PLOTTING FIGURES
% *************************************************************************
close all
colour_type = {'DEFAULT', ...
    'parula', 'jet', 'hsv', 'cool', ...
    'spring', 'summer', 'autumn', 'winter', ...
    'gray', 'copper',...
    'red', 'green', 'aqua', 'blue', 'purple',...
    };

centre_wavelength = zeros(size(file_name));
% centre_wavelength = [459, 500, 550, 600, 650, 700, 750, 790];
% centre_wavelength = [592, 592, 676.2, 744.5, 831.0];
% centre_wavelength = [698.5, 698.5, 698.5, 698.5, 698.5];

files_to_plot = 1:1:number_of_files;
% files_to_plot = 1:1:6;
% files_to_plot = 2:1:number_of_files;
% files_to_plot = 3:4;

channels_to_plot = [2,3];
% [channels_to_plot, ~] = listdlg('PromptString', 'Channels to plot',...
%                                 'SelectionMode', 'multiple', ...
%                                 'ListString', channel_name(1,2:end),...
%                                 'InitialValue', channels_to_plot - 1);
% channels_to_plot = channels_to_plot + 1;

menu_subplots = 2;
if size(channels_to_plot,2) > 1
%     menu_subplots = menu('Plot each channel in a different subplot?', 'NO', 'YES (horizontal)', 'YES (vertical)');
end
if menu_subplots == 2 % horizontal
    layout = [1,2];
elseif menu_subplots == 3 % vertical
    layout = [2,1];
end

figure_picoscope = figure('Units','normalized','Position',[0.01 0.07 0.95 0.8],'tag','figure_picoscope');
% plotting raw or smoothed channels
legend_A = {}; % initialising the legend cell
legend_B = {}; % initialising the legend cell
% for i = files_to_plot % files
for i = files_to_plot(end:-1:1) % files
%     for j = 2:1:size(plot_data{i},2) % channels
%     for j = 2:1:3 % channels
    for j = channels_to_plot % channels
        if strcmp(channel_name{i,j}, 'Channel A')
            if menu_subplots == 1
                if size(channels_to_plot,2)>1
                    yyaxis left
                end
            elseif menu_subplots == 2
                subplot(layout(1),layout(2),1)
            end
            
            
            ylabel('Channel A (V)')
            if strcmp(options(selected_options), 'Normalisation')
                ylabel('Channel A normalised')
            end
            grid on
            legend_A{end+1} = [file_name{i}(1:end-4) ...
%                 ' // ' channel_name{i,j} ...
%                 ' // ' num2str(laser_power(i), '%.0f') ' mW'...
%                 ' // ND = ' num2str(ND_A(i), '%.0f') ...
%                 ' x ' num2str(channel_factor(j)*measurement_factor(i), '%01.2f')...
%                 ' // centre = ' num2str(centre_wavelength(i), '%03.1f') ' nm'...
                ];
%             ylim([-1,4])
%             pause(0.1)
        elseif strcmp(channel_name{i,j}, 'Channel B')
            if menu_subplots == 1
                if size(channels_to_plot,2)>1
                    yyaxis right
                end
            elseif menu_subplots == 2
                subplot(layout(1),layout(2),2)
            end
            
            ylabel('Channel B (V)')
            if strcmp(options(selected_options), 'Normalisation')
                ylabel('Channel B normalised')
            end
            grid on
            legend_B{end+1} = [file_name{i}(1:end-4) ...
%                 ' // ' channel_name{i,j} ...
%                 ' // ' num2str(laser_power(i), '%.0f') ' mW'...
%                 ' // ND = ' num2str(ND_B(i), '%.0f') ...     
%                 ' // \lambda = ' num2str(centre_wavelength(i), '%03.0f') ' nm'...
%                 ' x ' num2str(channel_factor(j)*measurement_factor(i), '%01.2f')...
%                 ' // centre = ' num2str(centre_wavelength(i), '%03.1f') ' nm'...   
                ];
%             ylim([-2,8])
        end
        
        h(i,j) = plot(plot_data{i}(:,1), ...
            plot_data{i}(:,j),...
            'LineWidth', 2); hold all
        
%         pause(0.5)
    end
end

%% colour scheme
figure(figure_picoscope)
menu_colour = 2;
% menu_colour = menu([channel_name{1,j} ' colour scheme'], colour_type);
[menu_colour, ~] = listdlg('PromptString', 'Colour scheme:',...
                           'SelectionMode', 'single', ...
                           'ListString', colour_type,...
                           'InitialValue', menu_colour);
for j = channels_to_plot
    for i = files_to_plot    
        if menu_colour > 1 
            colour_RGB = colour_gradient(i, number_of_files, colour_type(menu_colour));
%             colour_RGB = colour_gradient(find(waveplate_angle_unique == waveplate_angle(i)), ...
%                 size(waveplate_angle_unique,2), ...
%                 colour_type(menu_colour));
            h(i,j).Color = colour_RGB;  
        end
        h(i,j).MarkerSize = 1;
        h(i,j).LineStyle = '-';
        h(i,j).LineWidth = 1;
    end
end

% yyaxis left
% ylabel('Channel A (V)')
% grid on
% yyaxis right
% ylabel('Channel B (V)')
% grid on

%% axis limits
figure(figure_picoscope)
legend_location = 'NE';
% x_limits = xlim;
x_limits = [-5,25];

input_title = 'Plot formatting';
input_data = {'Legend Location', 'X Axis Min', 'X Axis Max'};
resize = 'on'; dim = [1 80];
valdef = {legend_location, num2str(x_limits(1)), num2str(x_limits(2))};
answer = inputdlg(input_data,input_title,dim,valdef,resize);
legend_location = answer{1};
x_limits = [str2double(answer{2}),str2double(answer{3})];

if menu_subplots == 1
    title(title_cell_channels, 'interpreter', 'none')
    xlabel('Time (ms)')
    legend_channels = [legend_A, legend_B];
    legend(legend_channels, 'Location', legend_location)
    axis auto
    xlim(x_limits)

    if size(channels_to_plot,2)>1
        yyaxis left
        ylim([0,7])
        yyaxis right
        ylim([0,2.5])
    end
    
elseif menu_subplots == 2
    subplot(layout(1),layout(2),1)
    legend(legend_A, 'Location', legend_location, 'interpreter', 'none')
    title(title_cell_channels, 'interpreter', 'none')
    xlabel('Time (ms)')
    xlim(x_limits)
    
    subplot(layout(1),layout(2),2)
    legend(legend_B, 'Location', legend_location, 'interpreter', 'none')
    title(title_cell_channels, 'interpreter', 'none')
    xlabel('Time (ms)')
    xlim(x_limits)
end

%% DIVIDING B/A
% *************************************************************************
if find(strcmp(options(selected_options), 'Channel B / Channel A'))
    divided_data = cell(size(plot_data));
    for i = 1:1:number_of_files % files
        divided_data{i} = plot_data{i}(:,3) ./ plot_data{i}(:,2);
    end
    title_cell_divided{1} = title_cell_channels{1};
end

%% plotting ratio B/A
if find(strcmp(options(selected_options), 'Channel B / Channel A'))
%     subplot(subplot_division)
    figure_division = figure('Units','normalized','Position',[0.01 0.07 0.95 0.8],'tag','figure_division');
    legend_division = {}; % initialising the legend cell
%     for i = files_to_plot % files
    for i = files_to_plot(end:-1:1) % files
        h_div(i) = plot(plot_data{i}(:,1), ...
            divided_data{i},...
            'LineWidth', 2); hold all
        legend_division{end+1} = [file_name{i}(1:end-4) ' // B/A'...
%             ' // centre = ' num2str(centre_wavelength(i), '%03.1f') ' nm'...
            ];
    end
    legend(legend_division, 'Location', 'best', 'interpreter', 'none')
    title(title_cell_divided, 'interpreter', 'none')
    ylabel('B/A')
    xlabel('Time (ms)')
    xlim([0,8.5])
%     xlim([-50,250])
%     xlim([-310,-90])
%     ylim([1,5])
    grid on
end

%% colour scheme
figure(figure_division)
% menu_colour = 2;
% menu_colour = menu([channel_name{1,j} ' colour scheme'], colour_type);
[menu_colour, ~] = listdlg('PromptString', 'Colour scheme:',...
                           'SelectionMode', 'single', ...
                           'ListString', colour_type,...
                           'InitialValue', menu_colour);
for i = files_to_plot    
    if menu_colour > 1 
        colour_RGB = colour_gradient(i, number_of_files, colour_type(menu_colour));
%             colour_RGB = colour_gradient(find(waveplate_angle_unique == waveplate_angle(i)), ...
%                 size(waveplate_angle_unique,2), ...
%                 colour_type(menu_colour));
        h_div(i).Color = colour_RGB;  
    end
    h_div(i).MarkerSize = 1;
    h_div(i).LineStyle = '-';
    h_div(i).LineWidth = 1;
end

%% SAVING FIGURES
% *************************************************************************
menu_save_figures = 1;
menu_save_figures = menu('Save Figures?', 'NO', 'YES');
if menu_save_figures == 2    
    if findobj('tag','figure_picoscope') ~= 0
        figure_save = figure_picoscope;
        file_name_save = file_name{1}(1:end-4);
        
        figure(figure_save)
        pause(0.1)
        [file_name_save,file_path_save,~] = uiputfile(['.' 'png'],...
            'File to Save the Figure',[file_path_save file_name_save]);
        hgexport(figure_save, [file_path_save file_name_save], hgexport('factorystyle'), 'Format', 'png')
        file_name_save = strrep(file_name_save, 'png', 'fig');    
        saveas(figure_save, [file_path_save file_name_save], 'fig');
        
    end
    
    if findobj('tag','figure_division') ~= 0
        figure_save = figure_division;
        file_name_save = file_name{1}(1:end-4);
        
        figure(figure_save)
        pause(0.1)
        [file_name_save,file_path_save,~] = uiputfile(['.' 'png'],...
            'File to Save the Figure',[file_path_save file_name_save]);
        hgexport(figure_save, [file_path_save file_name_save], hgexport('factorystyle'), 'Format', 'png')
        file_name_save = strrep(file_name_save, 'png', 'fig');    
        saveas(figure_save, [file_path_save file_name_save], 'fig');
        
    end
end

