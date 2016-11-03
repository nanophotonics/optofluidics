% Created by Ana (aa938)
% Reads txt files created with the PicoScope software
% Savitzky-Golay filtering
% Plots the channels
%#ok<*ST2NM>

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
                                      'PicoScope Files to Read (use CTRL to select multiple files)',...
                                      file_path,...
                                      'MultiSelect','on');
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
measurement_numbers = zeros(size(file_name));
raw_data = cell(size(file_name));
for i = 1:1:number_of_files
    disp(['Reading File ' num2str(i) '/' num2str(number_of_files)])
    measurement_numbers(i) = str2double(file_name{i}(10:13));
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
           'Channel B / Channel A',...
           'Read Info File',...
           'Calculate Power',...
           };
selected_options = [5,6,7];
[selected_options, ~] = listdlg('PromptString', 'Select options:',...
                                'SelectionMode', 'multiple', ...
                                'ListString', options, ...
                                'InitialValue', selected_options, ...
                                'OKString', 'OK', ...
                                'CancelString', 'NONE');
                            
colour_type = {'DEFAULT', ...
               'parula', 'jet', 'hsv', 'cool', ...
               'spring', 'summer', 'autumn', 'autumn reversed', 'winter', ...
               'gray', 'copper',...
               'red', 'green', 'aqua', 'blue', 'purple',...
               };
selected_colour = 2;

%% smoothing
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
    dlg_options.WindowStyle = 'normal'; dlg_options.Resize = 'on'; dim = [1 80];
    default_values = {num2str(polynomial_order),num2str(frame_size)};
    answer = inputdlg(input_data, input_title, dim, default_values, dlg_options);
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

%% normalisation
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

%% standard deviation
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

%% period selection
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

%% READING INFO FILE
% *************************************************************************
if find(strcmp(options(selected_options), 'Read Info File'))
    [file_name_info, file_path, ~] = uigetfile('.xlsx',...
                                          'Select one info file to read', ...
                                          file_path, ...
                                          'MultiSelect','off');
    file_path_save = file_path;
    info_data = readtable([file_path file_name_info]);
    info_numbers = info_data{:,'Number'};
    [~,measurement_indices,~] = intersect(info_numbers, measurement_numbers);
    ND_A = info_data{measurement_indices,'ND_A'};
    ND_B = info_data{measurement_indices,'ND_B'};
    waveplate_angle = info_data{measurement_indices,'Angle_deg'}; % degrees
else
    if find(strcmp(channel_name(1,:), 'Channel A'))
        ND_A = zeros(number_of_files,1);
    end
    if find(strcmp(channel_name(1,:), 'Channel B'))
        ND_B = zeros(number_of_files,1);
    end
    waveplate_angle = zeros(number_of_files,1);
end

%% ND filters
% *************************************************************************

ND = [];
ND_title = {};
if find(strcmp(channel_name(1,:), 'Channel A'))
    ND = [ND, ND_A];
    ND_title{end+1} = 'ND A';
end
if find(strcmp(channel_name(1,:), 'Channel B'))
    ND = [ND, ND_B];
    ND_title{end+1} = 'ND B';
end
ND = dialog_table(file_name, ND_title, ND);

channel_data = period_data;
for i = 1:1:number_of_files % files
    for j = 2:1:size(channel_name,2)-1 % channels
        channel_data{i}(:,j) = channel_data{i}(:,j)*10^ND(i,j-1);
    end
end
plot_data = channel_data;

%% calculating the power
% *************************************************************************
if find(strcmp(options(selected_options), 'Calculate Power'))
    reference_angle = 30; % degrees
    reference_power = 16; % mW
    input_title = 'Reference'; 
    input_data = {'Reference Angle','Reference Power (mW)'};
    dlg_options.WindowStyle = 'normal'; dlg_options.Resize = 'on'; dim = [1 80];
    default_values = {num2str(reference_angle),num2str(reference_power)};
    answer = inputdlg(input_data, input_title, dim, default_values, dlg_options);
    reference_angle = str2double(answer{1});
    reference_power = str2double(answer{2}); 
    
    reference_index = find(waveplate_angle == reference_angle);
    reference_voltage = max(plot_data{reference_index}(:,2));
    laser_power = zeros(number_of_files,1);
    for i = 1:1:number_of_files % files
        laser_power(i) = reference_power * max(plot_data{i}(:,2))/reference_voltage;
        % Use power_fitted and compare the results.
    end
    
end


%% PLOTTING FIGURES
% *************************************************************************
close all
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
            legend_A{end+1} = [num2str(measurement_numbers(i), '%02.f')...
                               ' // ' num2str(laser_power(i), '%03.0f') ' mW'...
                               ];
%                 file_name{i}(1:end-4) ...
%                 ' // ' channel_name{i,j} ...

%                 ' // ND = ' num2str(ND_A(i), '%.0f') ...
%                 ' x ' num2str(channel_factor(j)*measurement_factor(i), '%01.2f')...
%                 ' // centre = ' num2str(centre_wavelength(i), '%03.1f') ' nm'...

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
            legend_B{end+1} = [num2str(measurement_numbers(i), '%02.0f')...
                               ' // ' num2str(laser_power(i), '%03.0f') ' mW'...
                               ];
%             legend_B{end+1} = [file_name{i}(1:end-4) ...
%                 ' // ' channel_name{i,j} ...
%                 ' // ' num2str(laser_power(i), '%.0f') ' mW'...
%                 ' // ND = ' num2str(ND_B(i), '%.0f') ...     
%                 ' // \lambda = ' num2str(centre_wavelength(i), '%03.0f') ' nm'...
%                 ' x ' num2str(channel_factor(j)*measurement_factor(i), '%01.2f')...
%                 ' // centre = ' num2str(centre_wavelength(i), '%03.1f') ' nm'...   
%                 ];
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
% selected_colour = 2;
% selected_colour = menu([channel_name{1,j} ' colour scheme'], colour_type);
[selected_colour, ~] = listdlg('PromptString', 'Colour scheme:',...
                           'SelectionMode', 'single', ...
                           'ListString', colour_type,...
                           'InitialValue', selected_colour);
for j = channels_to_plot
    for i = files_to_plot    
        if selected_colour > 1 
            colour_RGB = colour_gradient(i, number_of_files, colour_type(selected_colour));
%             colour_RGB = colour_gradient(find(waveplate_angle_unique == waveplate_angle(i)), ...
%                 size(waveplate_angle_unique,2), ...
%                 colour_type(selected_colour));
            h(i,j).Color = colour_RGB;  
        end
        h(i,j).MarkerSize = 1;
        h(i,j).LineStyle = '-';
        h(i,j).LineWidth = 2;
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
default_values = {legend_location, num2str(x_limits(1)), num2str(x_limits(2))};
answer = inputdlg(input_data, input_title, dim, default_values, dlg_options);
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
    legend(legend_division, 'Location', 'NEO', 'interpreter', 'none')
    title(title_cell_divided, 'interpreter', 'none')
    ylabel('B/A')
    xlabel('Time (ms)')
    xlim([0,8.2])
%     xlim([-50,250])
%     xlim([-310,-90])
%     ylim([1,5])
    grid on
end

%% colour scheme
figure(figure_division)
% selected_colour = 2;
% selected_colour = menu([channel_name{1,j} ' colour scheme'], colour_type);
[selected_colour, ~] = listdlg('PromptString', 'Colour scheme:',...
                           'SelectionMode', 'single', ...
                           'ListString', colour_type,...
                           'InitialValue', selected_colour);
for i = files_to_plot    
    if selected_colour > 1 
        colour_RGB = colour_gradient(i, number_of_files, colour_type(selected_colour));
%             colour_RGB = colour_gradient(find(waveplate_angle_unique == waveplate_angle(i)), ...
%                 size(waveplate_angle_unique,2), ...
%                 colour_type(selected_colour));
        h_div(i).Color = colour_RGB;  
    end
    h_div(i).MarkerSize = 1;
    h_div(i).LineStyle = '-';
    h_div(i).LineWidth = 2;
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

