% Created by Ana (aa938)
% Reads txt files created with the Thorlabs Powermeter

clc
clear
close all
figures = {};

%% CHOOSING FILES
% *************************************************************************

% specify default path
% folder_path = 'R:\aa938\NanoPhotonics\Laboratory\';
folder_path = 'R:\3-Temporary\aa938\';
% folder_path = 'R:\3-Temporary\os354\';

number_of_folders = 2;
% number_of_folders = menu('Where are the files located?', 'SINGLE folder', 'MULTIPLE folders');

if number_of_folders == 1
    % pop up window to choose the file(s) to read from a SINGLE FOLDER
    [file_name, folder_path, ~] = uigetfile('.txt',...
                                          'PicoScope Files to Read (use CTRL to select multiple files)',...
                                          folder_path,...
                                          'MultiSelect','on');
    file_name = cellstr(file_name); % convert to cell array of strings
    file_path = cell(size(file_name));
    for i = 1:1:size(file_name,2)
        file_path{i} = [folder_path file_name{i}];
    end
else
    % pop up window to choose the file(s) to read from MULTIPLE FOLDERS
    paths = uipickfiles('FilterSpec', [folder_path '*.txt']);
    file_path = {};
    for i = 1:1:size(paths,2)
        if strfind(paths{i}, 'txt')
            file_path{end+1} = paths{i};
        else
            directory_contents = dir(paths{i});
            for j = 3:1:size(directory_contents,1)
%             for j = 3
                file_path{end+1} = [paths{i} '\' directory_contents(j).name];
                if strfind(file_path{end}, 'txt')
                else % delete all non-txt files
                    file_path(end) = [];
                end
            end
        end
    end
    folder_path = file_path{1};
    file_name = cell(size(file_path));
    for i = 1:1:size(file_path,2)
        slash = strfind(folder_path, '\');
        file_name{i} = file_path{i}(slash(end)+1:end);
    end
end
number_of_files = size(file_name,2);


folder_path_save = folder_path;
if strfind(folder_path, 'Laboratory')
    slash = strfind(folder_path, '\');
    slash_index = find(slash > strfind(folder_path, 'Laboratory')+11);
%     title_cell_channels{1} = folder_path(strfind(folder_path, 'Laboratory')+11:end-1);
    title_cell_channels{1} = folder_path(strfind(folder_path, 'Laboratory')+11:slash(slash_index(1))-1);
elseif strfind(folder_path, 'aa938')
    slash = strfind(folder_path, '\');
    slash_index = find(slash > strfind(folder_path, 'aa938')+6);
%     title_cell_channels{1} = folder_path(strfind(folder_path, 'aa938')+6:end-1);
    title_cell_channels{1} = folder_path(strfind(folder_path, 'aa938')+6:slash(slash_index(1))-1);
else
    title_cell_channels{1} = folder_path;
end    


%% READING DATA
% *************************************************************************

% reading the data from the txt files
% it takes a while because the files are large
measurement_numbers = zeros(size(file_name));
raw_data = cell(size(file_name));
plot_data = cell(size(file_name));
for i = 1:1:number_of_files
    disp(['Reading File ' num2str(i) '/' num2str(number_of_files)])
%     measurement_numbers(i) = str2double(file_name{i}(10:13));
    measurement_numbers(i) = 0;
    raw_data{i} = readtable(file_path{i});
    
    % raw_data{i}(:,1) = time stamp
    % raw_data{i}(:,2) = power 
    % raw_data{i}(:,3) = units
    
    xaxis_units = '(s)';
    channel_name{i,1} = 'Time';
    channel_unit{i,1} = '(s)';
    channel_name{i,2} = 'Power Meter';
    channel_unit{i,2} = ['(' raw_data{1}.Var3{1} ')'];
%     yaxis_units = raw_data{1}.Var3(1);
    % TO TO: change all powers to Watts
    
    % convert time stamp to milliseconds
    date_vector = datevec(datestr(raw_data{i}.Var1(:),'dd-mm-yyyy HH:MM:SS.FFF'));
    for j = 1:1:size(date_vector,1)
        time_data{i}(j,1) = etime(date_vector(j,:),date_vector(1,:));
    end
    time_data{i}(:,2) = raw_data{i}.Var2(:);
    
end
disp('Finished reading all files!')
% figure
% plot(time_data{i}(:,1),time_data{i}(:,2))

%% OPTIONS
% *************************************************************************
options = {'Savitzky-Golay Filtering', ...
           'Normalisation', ...
           'Average & Standard Dev.',...
           'Period Selection',...
           'Read Info File',...
           'Calculate Power',...
           'Channel B / Channel A',...
           };
selected_options = [3];
% selected_options = 2;
% selected_options = 5;
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
selected_colour = 1;

%% smoothing
% *************************************************************************

smoothed_data = time_data;   
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
%             normalised_data{i}(:,j) = smoothed_data{i}(:,j) / smoothed_data{i}(t_index,j);
            normalised_data{i}(:,j) = smoothed_data{i}(:,j) / max(smoothed_data{i}(:,j));
        end
    end
end

%% standard deviation
% *************************************************************************
average = zeros(size(channel_name));
stdev = zeros(size(channel_name));
if find(strcmp(options(selected_options), 'Average & Standard Dev.'))
    for i = 1:1:number_of_files
%         threshold_ref = max(period_data{i}(:,2))/2; % V
%         [t_ref,~] = find(period_data{i}(:,2) > threshold_ref);
%         t_ref(end-9:end) = [];
%         t_ref(1:10) = [];    
        t_ref = 1:1:size(normalised_data{i}(:,2),1);
        for j = 2:1:size(channel_name,2)
            average(i,j) = mean(normalised_data{i}(t_ref,j));
            stdev(i,j) = std(normalised_data{i}(t_ref,j));
%             normalised_data{i}(:,j) = normalised_data{i}(:,j) - average(i,j);
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

if find(cell2mat(strfind(channel_name(1,:), 'Channel A')))
    ND_A = zeros(size(file_path));
end
if find(cell2mat(strfind(channel_name(1,:), 'Channel B')))
    ND_B = zeros(size(file_path));
end
if find(cell2mat(strfind(channel_name(1,:), 'Power Meter')))
    ND_PM = zeros(size(file_path));
end
waveplate_angle = zeros(size(file_path));
% wavelength = zeros(size(file_path));
sample = cell(size(file_path));
    
if find(strcmp(options(selected_options), 'Read Info File'))
    [file_name_info, folder_path, ~] = uigetfile('.xlsx',...
                                          'Select one info file to read', ...
                                          folder_path, ...
                                          'MultiSelect','off');
    folder_path_save = folder_path;
    info_data = readtable([folder_path file_name_info]);
    info_numbers = info_data{:,'Number'};
    [~,measurement_indices,~] = intersect(info_numbers, measurement_numbers);
%     ND_A = info_data{measurement_indices,'ND_A'};
%     ND_B = info_data{measurement_indices,'ND_B'};
%     waveplate_angle = info_data{measurement_indices,'Angle_deg'}; % degrees
%     wavelength = info_data{measurement_indices,'Wavelength_nm'}; % nm
    waveplate_angle = zeros(size(file_path));
    wavelength = zeros(size(file_path));
    for i = 1:1:number_of_files
        waveplate_angle(i) = info_data{info_data.Number == measurement_numbers(i), 'Angle_deg'}; % degrees
        wavelength(i) = info_data{info_data.Number == measurement_numbers(i), 'Wavelength_nm'}; % nm
        sample(i) = info_data{info_data.Number == measurement_numbers(i), 'Sample'};
    end
end

%% ND filters
% *************************************************************************

default_values = [];
dialog_title = {};
dialog_title{end+1} = 'Waveplate angle (deg)';
default_values = [default_values; waveplate_angle];
if find(cell2mat(strfind(channel_name(1,:), 'Channel A')))
    default_values = [default_values; ND_A];
    dialog_title{end+1} = 'ND A';
end
if find(cell2mat(strfind(channel_name(1,:), 'Channel B')))
    default_values = [default_values; ND_B];
    dialog_title{end+1} = 'ND B';
end
if find(cell2mat(strfind(channel_name(1,:), 'Power Meter')))
    default_values = [default_values; ND_PM];
    dialog_title{end+1} = 'ND Power Meter';
end


% ND = default_values;
ND = default_values(2:end,:);
% answer_values = dialog_table(file_name, dialog_title, default_values);
% waveplate_angle = answer_values(:,1);
% ND = answer_values(:,2:end);

channel_data = period_data;
for i = 1:1:number_of_files % files
    for j = 2:1:size(channel_name,2) % channels
        channel_data{i}(:,j) = channel_data{i}(:,j)*10^ND(j-1,i);
    end
end
plot_data = channel_data;

%% calculating the power
% *************************************************************************
power_data = channel_data;
laser_power = zeros(number_of_files,1);
if find(strcmp(options(selected_options), 'Calculate Power'))
    reference_angle = 35; % degrees
    reference_power = 88; % mW
    input_title = 'Reference'; 
    input_data = {'Reference Angle (deg)','Reference Power (mW)'};
    dlg_options.WindowStyle = 'normal'; dlg_options.Resize = 'on'; dim = [1 80];
    default_values = {num2str(reference_angle),num2str(reference_power)};
    answer = inputdlg(input_data, input_title, dim, default_values, dlg_options);
    reference_angle = str2double(answer{1});
    reference_power = str2double(answer{2}); 
    
    reference_index = find(waveplate_angle == reference_angle);
    reference_voltage = max(plot_data{reference_index(end)});
    for i = 1:1:number_of_files % files
        laser_power(i) = reference_power * max(channel_data{i}(:,2))/reference_voltage(2);
%         for j = 2:1:size(channel_name,2) % channels
%             power_data{i}(:,j) = reference_power * channel_data{i}(:,j)/reference_voltage(j);
%         end
%         yaxis_units = 'mW';
    end
    laser_power_fitted = power_fitted(reference_angle, reference_power, waveplate_angle);
    laser_power_spline = power_spline(reference_angle, reference_power, waveplate_angle);
    
%     close all
    figures{end+1} = figure('Units','normalized','Position',[0.01 0.07 0.95 0.8]);
    
    subplot(1,2,1)
%     yyaxis left
    plot(measurement_numbers, laser_power, '.', 'MarkerSize', 16), hold all
    plot(measurement_numbers, laser_power_fitted, 'o', 'MarkerSize', 4, 'LineWidth', 1), hold all
    plot(measurement_numbers, laser_power_spline, 'x', 'MarkerSize', 8, 'LineWidth', 1.5), hold all
    title(['Ref angle = ' num2str(reference_angle) ...
        '\circ, Ref power = ' num2str(reference_power) ' mW, ' ...
        'Ref # = ' num2str(measurement_numbers(reference_index(end)))])
    xlabel('Measurement Number')
    ylabel('Power (mW)')
    legend('Picoscope reference power', ...
           'Ti:Sa calibration coefficients', ...
           'Ti:Sa calibration interpolation', ...
           'Location', 'NW')
    grid on
    yyaxis right
    plot(measurement_numbers, waveplate_angle, '<', 'MarkerSize', 6, 'LineWidth', 1.5), hold all
    ylabel('Waveplate Angle (degrees)')
    
    subplot(1,2,2)
    plot(waveplate_angle, laser_power, '.', 'MarkerSize', 16), hold all
%     plot(waveplate_angle, laser_power_fitted, 'o', 'MarkerSize', 4, 'LineWidth', 1), hold all
    plot(min(waveplate_angle):max(waveplate_angle), ...
         power_fitted(reference_angle, reference_power, min(waveplate_angle):max(waveplate_angle)), ...
         '-', 'MarkerSize', 4, 'LineWidth', 1), hold all
    plot(waveplate_angle, laser_power_spline, 'x', 'MarkerSize', 8, 'LineWidth', 1.5), hold all
    title(['Ref angle = ' num2str(reference_angle) ...
        '\circ, Ref power = ' num2str(reference_power) ' mW, ' ...
        'Ref # = ' num2str(measurement_numbers(reference_index(end)))])
    xlabel('Waveplate Angle (degrees)')
    ylabel('Power (mW)')
    legend('Picoscope reference power', ...
           'Ti:Sa calibration coefficients', ...
           'Ti:Sa calibration interpolation', ...
           'Location', 'NW')
    grid on
    
    text(waveplate_angle + 0.1, ...
         laser_power, ...
         strread(num2str(measurement_numbers),'%s'), ...
         'HorizontalAlignment', 'center',...
         'Color', 'k')
    
    
end
plot_data = power_data;

%% PLOTTING FIGURES
% *************************************************************************
% close all
% wavelength = zeros(size(file_name));
% wavelength = [459, 500, 550, 600, 650, 700, 750, 790];
% wavelength = [592, 592, 676.2, 744.5, 831.0];
% wavelength = [698.5, 698.5, 698.5, 698.5, 698.5];

% files_to_plot = 1:1:number_of_files;
files_to_plot = number_of_files:-1:1;
% files_to_plot = 1:1:6;
% files_to_plot = 2:1:number_of_files;
% files_to_plot = 3:4;

channels_to_plot = [2];
% [channels_to_plot, ~] = listdlg('PromptString', 'Channels to plot',...
%                                 'SelectionMode', 'multiple', ...
%                                 'ListString', channel_name(1,2:end),...
%                                 'InitialValue', channels_to_plot - 1);
% channels_to_plot = channels_to_plot + 1;

menu_subplots = 1;
if size(channels_to_plot,2) > 1
    menu_subplots = 2;
    menu_subplots = menu('Plot each channel in a different subplot?', 'NO', 'YES (horizontal)', 'YES (vertical)');
end
if menu_subplots == 2 % horizontal
    layout = [1,2];
elseif menu_subplots == 3 % vertical
    layout = [2,1];
end

figures{end+1} = figure('Units','normalized','Position',[0.01 0.07 0.95 0.8]);
figure_picoscope = figures{end};
% plotting raw or smoothed channels
legend_A = {}; % initialising the legend cell
legend_B = {}; % initialising the legend cell
for i = files_to_plot % files
%     disp([num2str(i) ' / ' num2str(max(size(files_to_plot))) ' files'])
%     for j = 2:1:size(plot_data{i},2) % channels
%     for j = 2:1:3 % channels
    for j = channels_to_plot % channels
%         disp([num2str(j) ' / ' num2str(max(size(channels_to_plot))) ' channels'])
        if strfind(channel_name{i,j}, 'Power Meter')
            if menu_subplots == 1
                if size(channels_to_plot,2)>1
                    yyaxis left
                end
            elseif menu_subplots == 2
                subplot(layout(1),layout(2),1)
            end
            
            if strcmp(options(selected_options), 'Normalisation')
                ylabel('Power Meter normalised')
            end
            grid on
            legend_A{end+1} = '';
%             legend_A{end+1} = file_name{i};
%             legend_A{end+1} = num2str(measurement_numbers(i), '%02.f');
            if find(strcmp(options(selected_options), 'Read Info File'))
                legend_A{end} = [legend_A{end} ' // ' sample{i}];
                legend_A{end} = [legend_A{end} ' // ' num2str(wavelength(i), '%03.0f') ' nm'];
                legend_A{end} = [legend_A{end} ' // ' num2str(waveplate_angle(i), '%02.0f') ' deg'];
            end
            if find(strcmp(options(selected_options), 'Calculate Power'))
                legend_A{end} = [legend_A{end} ' // ' num2str(laser_power(i), '%03.0f') ' mW'];
            end
            if find(strcmp(options(selected_options), 'Average & Standard Dev.'))
                legend_A{end} = [legend_A{end} 'Average = ' num2str(average(i,2), '%03.3f') ' W'];
                legend_A{end} = [legend_A{end} ' //  Std = ' num2str(stdev(i,2), '%03.6f') ' W'];
                legend_A{end} = [legend_A{end} ' = ' num2str(stdev(i,2)/average(i,2)*100, '%03.2f') ' %'];
            end
            legend_A{end} = [legend_A{end} ' // ' file_name{i}];
        end

        
        h(i,j) = plot(plot_data{i}(:,1), ...
            plot_data{i}(:,j),...
            'LineWidth', 2); hold all
        xlabel([channel_name{1,1} ' ' channel_unit{1,1}])
        ylabel([channel_name{1,2} ' ' channel_unit{1,2}])
%         pause(0.5)
    end
end

%% colour scheme
figure(figure_picoscope)
[selected_colour, ~] = listdlg('PromptString', 'Colour scheme:',...
                           'SelectionMode', 'single', ...
                           'ListString', colour_type,...
                           'InitialValue', selected_colour);
for j = channels_to_plot
    for i = files_to_plot    
        if selected_colour > 1 
            if find(strcmp(options(selected_options), 'Calculate Power'))
                colour_RGB = colour_gradient(find(unique(laser_power) == laser_power(i)), ...
                    size(unique(laser_power),1), ...
                    colour_type(selected_colour));
            elseif find(strcmp(options(selected_options), 'Read Info File'))
                colour_RGB = colour_gradient(find(unique(waveplate_angle) == waveplate_angle(i)), ...
                    size(unique(waveplate_angle),1), ...
                    colour_type(selected_colour));
            else
                colour_RGB = colour_gradient(i, number_of_files, colour_type(selected_colour));
            end

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

%% axis limits and legend location
figure(figure_picoscope)
legend_location = 'best';
x_limits = xlim;

input_title = 'Plot formatting';
input_data = {'Legend Location', 'X Axis Min', 'X Axis Max'};
default_values = {legend_location, num2str(x_limits(1)), num2str(x_limits(2))};
dlg_options.WindowStyle = 'normal'; dlg_options.Resize = 'on'; dim = [1 80];
answer = inputdlg(input_data, input_title, dim, default_values, dlg_options);
legend_location = answer{1};
x_limits = [str2double(answer{2}),str2double(answer{3})];
set(gca,'FontSize', 14)

if menu_subplots == 1
    title(title_cell_channels, 'interpreter', 'none')
    xlabel(['Time ' channel_unit{1,1}])
    legend_channels = [legend_A, legend_B];
    if strfind(legend_location, 'none')
    elseif strfind(legend_location, 'off')
    else
        legend(legend_channels, 'Location', legend_location, 'interpreter', 'none')
    end
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
    if strfind(legend_location, 'none')
    elseif strfind(legend_location, 'off')
    else
        legend(legend_A, 'Location', legend_location, 'interpreter', 'none')
    end
    title(title_cell_channels, 'interpreter', 'none')
    xlabel(['Time ' channel_unit{1,1}])
    xlim(x_limits)
    [~,measurement_indices] = unique(measurement_numbers);
    for i = 1:1:size(measurement_indices,1)
        text(mean(xlim), ...
            max(plot_data{measurement_indices(i)}(:,2))*0.98, ...
            legend_A{files_to_plot(measurement_indices(i))}, ...
            'HorizontalAlignment', 'center',...
            'VerticalAlignment', 'top', ...
            'Color', 'k',...
            'FontSize', 8)
    end
     
    
    subplot(layout(1),layout(2),2)
    if strfind(legend_location, 'none')
    elseif strfind(legend_location, 'off')
    else
        legend(legend_B, 'Location', legend_location, 'interpreter', 'none')
    end
    title(title_cell_channels, 'interpreter', 'none')
    xlabel(['Time ' channel_unit{1,1}])
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
    figures{end+1} = figure('Units','normalized','Position',[0.01 0.07 0.95 0.8]);
    figure_division = figures{end};
    legend_division = {}; % initialising the legend cell
    for i = files_to_plot % files
%     for i = files_to_plot(end:-1:1) % files
        h_div(i) = plot(plot_data{i}(:,1), ...
            divided_data{i},...
            'LineWidth', 1); hold all
%         legend_division{end+1} = file_name{i};
        legend_division{end+1} = num2str(measurement_numbers(i), '%02.f');
        if find(strcmp(options(selected_options), 'Read Info File'))
            legend_division{end} = [legend_division{end} ' // ' sample{i}];
            legend_division{end} = [legend_division{end} ' // ' num2str(wavelength(i), '%03.0f') ' nm'];
            legend_division{end} = [legend_division{end} ' // ' num2str(waveplate_angle(i), '%02.0f') ' deg'];
        end
        if find(strcmp(options(selected_options), 'Calculate Power'))
            legend_division{end} = [legend_division{end} ' // ' num2str(laser_power(i), '%03.0f') ' mW'];
        end
    end
    title(title_cell_divided, 'interpreter', 'none')
    ylabel('B/A')
    xlabel(['Time ' channel_unit{1,1}])
%     xlim([0,8.2])
%     xlim([-50,250])
%     xlim([-310,-90])
%     ylim([1,5])
%     ylim([4,5])
    grid on


    %% colour scheme
    figure(figure_division)
    [selected_colour, ~] = listdlg('PromptString', 'Colour scheme:',...
                               'SelectionMode', 'single', ...
                               'ListString', colour_type,...
                               'InitialValue', selected_colour);
    for i = files_to_plot    
        if selected_colour > 1 
            if find(strcmp(options(selected_options), 'Calculate Power'))
                colour_RGB = colour_gradient(find(unique(laser_power) == laser_power(i)), ...
                    size(unique(laser_power),1), ...
                    colour_type(selected_colour));
            elseif find(strcmp(options(selected_options), 'Read Info File'))
                colour_RGB = colour_gradient(find(unique(waveplate_angle) == waveplate_angle(i)), ...
                    size(unique(waveplate_angle),1), ...
                    colour_type(selected_colour));
            else
                colour_RGB = colour_gradient(i, number_of_files, colour_type(selected_colour));
            end
            h_div(i).Color = colour_RGB;  
        end
        h_div(i).MarkerSize = 1;
        h_div(i).LineStyle = '-';
        h_div(i).LineWidth = 1;
    end
    
    %% axis limits and legend location
    figure(figure_division)
    legend_location = 'NEO';
%     x_limits = xlim;
    x_limits = [0,8];
%     y_limits = ylim;
    y_limits = [0,1];

    input_title = 'Plot formatting';
    input_data = {'Legend Location', ...
        'X Axis Min', 'X Axis Max', ...
        'Y Axis Min', 'Y Axis Max'};
    default_values = {legend_location, ...
        num2str(x_limits(1)), num2str(x_limits(2)), ...
        num2str(y_limits(1)), num2str(y_limits(2))};
    dlg_options.WindowStyle = 'normal'; dlg_options.Resize = 'on'; dim = [1 80];
    answer = inputdlg(input_data, input_title, dim, default_values, dlg_options);
    legend_location = answer{1};
    x_limits = [str2double(answer{2}),str2double(answer{3})];
    y_limits = [str2double(answer{4}),str2double(answer{5})];
    
    if strfind(legend_location, 'none')
    elseif strfind(legend_location, 'off')
    else
        legend(legend_division, 'Location', 'NEO', 'interpreter', 'none')
    end
    
    xlim(x_limits)
    ylim(y_limits)

end

%% SAVING FIGURES
% *************************************************************************
% menu_save_figures = 1;
% menu_save_figures = menu('Save Figures?', 'NO', 'YES');
% if menu_save_figures == 2    
%     if findobj('tag','figure_picoscope') ~= 0
%         figure_save = figure_picoscope;
%         file_name_save = file_name{1}(1:end-4);
%         
%         figure(figure_save)
%         pause(0.1)
%         [file_name_save,folder_path_save,~] = uiputfile(['.' 'png'],...
%             'File to Save the Figure',[folder_path_save file_name_save]);
%         hgexport(figure_save, [folder_path_save file_name_save], hgexport('factorystyle'), 'Format', 'png')
%         file_name_save = strrep(file_name_save, 'png', 'fig');    
%         saveas(figure_save, [folder_path_save file_name_save], 'fig');
%         
%     end
%     
%     if findobj('tag','figure_division') ~= 0
%         figure_save = figure_division;
%         file_name_save = file_name{1}(1:end-4);
%         
%         figure(figure_save)
%         pause(0.1)
%         [file_name_save,folder_path_save,~] = uiputfile(['.' 'png'],...
%             'File to Save the Figure',[folder_path_save file_name_save]);
%         hgexport(figure_save, [folder_path_save file_name_save], hgexport('factorystyle'), 'Format', 'png')
%         file_name_save = strrep(file_name_save, 'png', 'fig');    
%         saveas(figure_save, [folder_path_save file_name_save], 'fig');
%         
%     end
% end

menu_save_figures = 1;
menu_save_figures = menu('Save Figures?', 'NO', 'YES');
if menu_save_figures == 2    
    for i = 1:1:max(size(figures))
        if findobj(figures{i}) ~= 0
            figure_save = figures{i};
            file_name_save = '';
            
            figure(figure_save)
            pause(0.1)
            [file_name_save,folder_path_save,~] = uiputfile(['.' 'png'],...
                'File to Save the Figure',[folder_path_save file_name_save]);
            hgexport(figure_save, [folder_path_save file_name_save], hgexport('factorystyle'), 'Format', 'png')
            file_name_save = strrep(file_name_save, 'png', 'fig');    
            saveas(figure_save, [folder_path_save file_name_save], 'fig');

        end
    end
end

