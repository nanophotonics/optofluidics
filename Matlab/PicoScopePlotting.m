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

% pop up window to choose the file(s) to read
[file_name, file_path, ~] = uigetfile('.txt',...
    'PicoScope Fil to Read (use CTRL to select multiple files)',file_path,'MultiSelect','on');
if strfind(file_path, 'Laboratory')
    title_cell_channels{1} = file_path(strfind(file_path, 'Laboratory')+11:end);
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
    raw_data{i} = dlmread([file_path file_name{i}], '\t',2,0);
    
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
        if strfind(channel_unit{i,j}, '(V)')
        elseif strfind(channel_unit{i,j}, '(mV)')
            channel_unit{i,j} = '(V)';
            raw_data{i}(:,j) = raw_data{i}(:,j)/1000;
        else
            disp(['Check Channel (' num2str(i) ',' num2str(j) ') Units!'])
        end
    end
end
disp('Finished reading all files!')

%% SMOOTHING
% *************************************************************************

menu_smoothing = 3;
menu_smoothing = menu('Smooth Data?', 'NO', 'YES: Savitzky-Golay filtering');

smoothed_data = raw_data;   
title_cell_channels{2} = 'Raw unfiltered data';
if menu_smoothing == 2 % Savitzky-Golay filtering
    % if polynomial_order = 1 then this becomes a moving average
    polynomial_order = 5; % must be less than the frame size
    frame_size = 2001; % must be odd
    
    input_title = 'Savitzky-Golay Filtering Parameters'; 
    input_data = {'Polynomial Order (odd) (1 = moving average):',...
        'Frame Size (odd):'};
    resize = 'on'; dimensions = [1 80];
    devault_values = {num2str(polynomial_order),num2str(frame_size)};
    answer = inputdlg(input_data, input_title, dimensions, devault_values, resize);
    polynomial_order = str2double(answer{1});
    frame_size = str2double(answer{2});    
    
    
    title_cell_channels{2} = ['Savitzky-Golay filtering. '...
        'Polynomial order = ' num2str(polynomial_order) '. '...
        'Frame size =  ' num2str(frame_size) ' points.'];
    
    for i = 1:1:number_of_files % files
        disp(['Smoothing File ' num2str(i) '/' num2str(number_of_files)])
        for j = 2:1:size(smoothed_data{i},2) % channels
            smoothed_data{i}(:,j) = sgolayfilt(raw_data{i}(:,j),polynomial_order,frame_size);
        end
    end
    disp('Finished smoothing all files!')
end

%% DIVIDING B/A
% *************************************************************************
menu_divide = 1;
menu_divide = menu('Divide Data?', 'NO', 'YES: B/A');
if menu_divide == 2
    divided_data = cell(size(smoothed_data));
    for i = 1:1:number_of_files % files
        divided_data{i} = smoothed_data{i}(:,3) ./ smoothed_data{i}(:,2);
    end
    title_cell_divided{1} = title_cell_channels{1};
end


%% PLOTTING FIGURES
% *************************************************************************
close all
figure_picoscope = figure('Units','normalized','Position',[0.01 0.07 0.75 0.8],'tag','figure_picoscope');
subplot_channels = subplot(1,1,1);
if menu_divide == 2
    subplot_channels = subplot(2,1,1);
    subplot_division = subplot(2,1,2);
end
channel_factor = [1,1,1,1];
measurement_factor = ones(size(file_name));
% measurement_factor = [1,1,1,1,1,1,1,1,1,1,1];
% measurement_factor = [1.86,1,1,1,1,1,1,1,1,1,1];
% measurement_factor = [26,1,26,28,22,1,1,1,1,1,1];
centre_wavelength = zeros(size(file_name));
% centre_wavelength = [592, 592, 676.2, 744.5, 831.0];
% centre_wavelength = [698.5, 698.5, 698.5, 698.5, 698.5];

% plotting raw or smoothed channels
subplot(subplot_channels)
legend_A = {}; % initialising the legend cell
legend_B = {}; % initialising the legend cell
for i = 1:1:number_of_files % files
    for j = 2:1:size(smoothed_data{i},2) % channels
%     for j = 2:1:3 % channels
%     for j = 3 % channels
        if strfind(channel_name{i,j}, 'Channel A')
            yyaxis left
            ylabel('Channel A (V)')
            grid on
            legend_A{end+1} = [file_name{i}(1:end-4) ' // ' ...
                channel_name{i,j} ...
%                 ' x ' num2str(channel_factor(j)*measurement_factor(i), '%01.2f')...
%                 ' // centre = ' num2str(centre_wavelength(i), '%03.1f') ' nm'...
                ];
%             ylim([-1,4])
            pause(0.1)
        elseif strfind(channel_name{i,j}, 'Channel B')
            subplot(subplot_channels)
            yyaxis right
            ylabel('Channel B (V)')
            grid on
            legend_B{end+1} = [file_name{i}(1:end-4) ' // ' ...
                channel_name{i,j} ...
%                 ' x ' num2str(channel_factor(j)*measurement_factor(i), '%01.2f')...
%                 ' // centre = ' num2str(centre_wavelength(i), '%03.1f') ' nm'...
                ];
%             ylim([-2,8])
        end
        
        plot(smoothed_data{i}(:,1), ...
            smoothed_data{i}(:,j)*channel_factor(j)*measurement_factor(i),...
            'LineWidth', 1), hold all
    end
end
yyaxis left
ylabel('Channel A (V)')
grid on
yyaxis right
ylabel('Channel B (V)')
grid on

legend_channels = [legend_A, legend_B];
legend(legend_channels, 'Location', 'S', 'interpreter', 'none')
title(title_cell_channels, 'interpreter', 'none')
xlabel('Time (ms)')
xlim([-50,250])

% plotting ratio B/A
if menu_divide == 2
    subplot(subplot_division)
    legend_division = {}; % initialising the legend cell
    for i = 1:1:number_of_files % files
        plot(smoothed_data{i}(:,1), ...
            divided_data{i},...
            'LineWidth', 2), hold all
        legend_division{end+1} = [file_name{i}(1:end-4) ' // B/A'...
%             ' // centre = ' num2str(centre_wavelength(i), '%03.1f') ' nm'...
            ];
    end
    legend(legend_division, 'Location', 'S', 'interpreter', 'none')
%     title(title_cell_divided, 'interpreter', 'none')
    ylabel('B/A')
    xlabel('Time (ms)')
    xlim([-50,250])
    ylim([1,2])
    grid on
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
            'File to Save the Figure',[file_path file_name_save]);
        hgexport(figure_save, [file_path_save file_name_save], hgexport('factorystyle'), 'Format', 'png')
        file_name_save = strrep(file_name_save, 'png', 'fig');    
        saveas(figure_save, [file_path_save file_name_save], 'fig');
        
    end
end

