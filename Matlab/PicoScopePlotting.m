clc
clear
close all

%% READING DATA
% *************************************************************************

% select default path
file_path = 'C:\Users\Ana Andres\Documents\NanoPhotonics\Laboratory\';

% pop up window to choose the file(s) to read
[file_name, file_path, ~] = uigetfile('.txt',...
    'PicoScope Fil to Read (use CTRL to select multiple files)',file_path,'MultiSelect','on');

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
    raw_data{i} = dlmread([file_path file_name{i}], '\t',2,0);
    % raw_data{i}(:,1) = time
    % raw_data{i}(:,2) = channel A
    % raw_data{i}(:,3) = channel B
end

% plotting the data
figure_picoscope = figure('Units','normalized','Position',[0.01 0.1 0.97 0.64],'tag','figure_picoscope');
hold all
channel_name = {'Time', 'Channel A', 'Channel B', 'Channel C', 'Channel D'};
legend_picoscope = {}; % initialising the legend cell
for i = 1:1:number_of_files % files
    for j = 2:1:size(raw_data{i},2) % channels
        plot(raw_data{i}(:,1), raw_data{i}(:,j))
        legend_picoscope{end+1} = [file_name{i}(1:end-4) ' // ' channel_name{j}];
    end
end
grid on
legend(legend_picoscope, 'Location', 'best')
title(file_path, 'interpreter', 'none')
xlabel('Time (ms)')
ylabel('Signal (mV)')


