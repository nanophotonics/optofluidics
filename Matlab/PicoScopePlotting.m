clc
clear
close all

%% READING DATA
% *************************************************************************

% select default path
file_path = 'C:\Users\Ana Andres\Documents\NanoPhotonics\Laboratory\2016.06.17 - 60 nm Au\';

% pop up window to choose the file(s) to read
[file_name, file_path, ~] = uigetfile('.txt',...
    'PicoScope Fil to Read (use CTRL to select multiple files)',file_path,'MultiSelect','on');
if strfind(file_path, 'Laboratory')
    title_cell{1} = file_path(strfind(file_path, 'Laboratory')+11:end);
else
    title_cell{1} = file_path;
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
end

%% SMOOTHING DATA
% *************************************************************************

menu_smoothing = 3;
menu_smoothing = menu('Smooth Data?', 'NO', 'YES: Savitzky-Golay filtering');

smoothed_data = raw_data;   
title_cell{2} = 'Raw unfiltered data';
if menu_smoothing == 2 % Savitzky-Golay filtering
    % if polynomial_order = 1 then this becomes a moving average
    polynomial_order = 4; % must be less than the frame size
    frame_size = 2001; % must be odd
    title_cell{2} = ['Savitzky-Golay filtering.'...
        ' Polynomial order = ' num2str(polynomial_order)...
        ' Frame size =  ' num2str(frame_size) ' points.'];
    
    for i = 1:1:number_of_files % files
        disp(['Smoothing File ' num2str(i) '/' num2str(number_of_files)])
        for j = 2:1:size(smoothed_data{i},2) % channels
            smoothed_data{i}(:,j) = sgolayfilt(raw_data{i}(:,j),polynomial_order,frame_size);
        end
    end
end


%% PLOTTING DATA
% *************************************************************************
% close all
figure_picoscope = figure('Units','normalized','Position',[0.01 0.1 0.97 0.64],'tag','figure_picoscope');
hold all
channel_name = {'Time', 'Channel A', 'Channel B', 'Channel C', 'Channel D'};
channel_factor = [1,1,1];
measurement_factor = [26,1,26,28,22,1,1,1,1,1,1];
centre_wavelength = [592, 592, 676.2, 744.5, 831.0];

legend_picoscope = {}; % initialising the legend cell
for i = 1:1:number_of_files % files
%     for j = 2:1:size(smoothed_data{i},2) % channels
    for j = 3 % channels
        plot(smoothed_data{i}(:,1), ...
            smoothed_data{i}(:,j)*channel_factor(j)*measurement_factor(i),...
            'LineWidth', 2)
        legend_picoscope{end+1} = [file_name{i}(1:end-4) ' // ' ...
            channel_name{j} ' x ' num2str(channel_factor(j)*measurement_factor(i), '%02.0f')...
            ' // centre = ' num2str(centre_wavelength(i), '%03.1f') ' nm'];
    end
end
grid on
legend(legend_picoscope, 'Location', 'best', 'interpreter', 'none')
title(title_cell, 'interpreter', 'none')
xlabel('Time (ms)')
ylabel('Signal (V)')
xlim([-50,250])


