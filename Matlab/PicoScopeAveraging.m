% Created by Ana (aa938)
% Reads txt files created with the PicoScope software
% Averages the files

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
    raw_data{i} = dlmread([file_path file_name{i}], '\t',2,0);
    
    % raw_data{i}(:,1) = time
    % raw_data{i}(:,2) = channel A
    % raw_data{i}(:,3) = channel B
    
    % read channel names and units
    file_identifier = fopen([file_path file_name{i}], 'rt');
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

plot_data = raw_data;


%% PLOTTING FIGURES
% *************************************************************************
close all
figure_picoscope = figure('Units','normalized','Position',[0.01 0.07 0.75 0.8],'tag','figure_picoscope');
subplot_channels = subplot(1,1,1);

files_to_plot = 1:1:number_of_files;

menu_channels = 1;
% menu_channels = menu('Channels to plot', 'A & B', 'Just A (ref)', 'Just B (signal)');
if menu_channels == 1
    channels_to_plot = 2:3; % A and B
elseif menu_channels == 2
    channels_to_plot = 2; % A
elseif menu_channels == 3
    channels_to_plot = 3; % B
end

menu_subplots = 2;
if size(channels_to_plot,2) > 1
%     menu_subplots = menu('Plot each channel in a different subplot?', 'NO', 'YES (horizontal)', 'YES (vertical)');
end
if menu_subplots == 2 % horizontal
    layout = [1,2];
elseif menu_subplots == 3 % vertical
    layout = [2,1];
end

% plotting raw or smoothed channels
subplot(subplot_channels)
legend_A = {}; % initialising the legend cell
legend_B = {}; % initialising the legend cell
for i = files_to_plot % files
    for j = channels_to_plot % channels
        if strfind(channel_name{i,j}, 'Channel A')
            if menu_subplots == 1
                subplot(subplot_channels)
                if size(channels_to_plot,2)>1
                    yyaxis left
                end
            elseif menu_subplots == 2
                subplot(layout(1),layout(2),1)
            end
            
            ylabel('Channel A (V)')
            grid on
            legend_A{end+1} = [file_name{i}(1:end-4) ...
                ' // ' channel_name{i,j} ...
                ];
%             pause(0.1)

        elseif strfind(channel_name{i,j}, 'Channel B')
            if menu_subplots == 1
                subplot(subplot_channels)
                if size(channels_to_plot,2)>1
                    yyaxis right
                end
            elseif menu_subplots == 2
                subplot(layout(1),layout(2),2)
            end
            
            ylabel('Channel B (V)')
            grid on
            legend_B{end+1} = [file_name{i}(1:end-4) ...
                ' // ' channel_name{i,j} ...
                ];
        end
        
        h(i,j) = plot(plot_data{i}(:,1), ...
            plot_data{i}(:,j),...
            'LineWidth', 1); hold all
        
%         pause(0.5)
    end
end

%% colour scheme
colour_type = {'DEFAULT', 'yellow' 'red', 'green', 'aqua', 'blue', 'purple', 'gray'};
for j = channels_to_plot
    menu_colour = 1;
%     menu_colour = menu([channel_name{1,j} ' colour scheme'], colour_type);
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
        h(i,j).LineWidth = 0.5;
    end
end


%% axis limits
legend_location = 'NEO';
x_limits = [-5,15];

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

%% SELECTING BAD FILES
% *************************************************************************

menu_ok = menu('Ready to select bad files?', 'OK');

number_of_bad_files = 0;
input_title = 'Selecting bad files';
input_data = {'Number of bad files'};
resize = 'on'; dim = [1 80];
valdef = {num2str(number_of_bad_files)};
answer = inputdlg(input_data,input_title,dim,valdef,resize);
number_of_bad_files = str2double(answer{1});

files_to_average = files_to_plot;
if number_of_bad_files > 0
    
    bad_files = cell(number_of_bad_files,1);
    input_data = cell(number_of_bad_files,1);
    for i = 1:1:number_of_bad_files
        input_data{i} = [num2str(i) '/' num2str(number_of_bad_files) ' file to discard'];
        bad_files{i} = '0';
    end
    input_title = 'Selecting bad files';
    resize = 'on'; dim = [1 80];
    valdef = bad_files;
    bad_files = inputdlg(input_data,input_title,dim,valdef,resize);

    for i = 1:1:number_of_bad_files
        [~,index_bad] = find(files_to_average == str2double(bad_files{i}));
        files_to_average(index_bad) = [];
    end
end

files_to_average = files_to_plot;
[index_bad, ~] = listdlg('PromptString', 'Select bad files(s)',...
    'SelectionMode', 'multiple', ...
    'ListString', file_name,...
    'InitialValue', []);
files_to_average(index_bad) = [];


%% AVERAGING DATA
% *************************************************************************
% make sure all measurements have the same time array
data_average = zeros(size(raw_data{i}));
for i = files_to_average
    data_average = data_average + raw_data{i}; 
end
data_average = data_average / size(files_to_average,2);

%% PLOTTING FIGURES
% *************************************************************************
close all
figure_picoscope = figure('Units','normalized','Position',[0.01 0.07 0.75 0.8],'tag','figure_picoscope');
subplot_channels = subplot(1,1,1);

menu_channels = 1;
% menu_channels = menu('Channels to plot', 'A & B', 'Just A (ref)', 'Just B (signal)');
if menu_channels == 1
    channels_to_plot = 2:3; % A and B
elseif menu_channels == 2
    channels_to_plot = 2; % A
elseif menu_channels == 3
    channels_to_plot = 3; % B
end

menu_subplots = 2;
if size(channels_to_plot,2) > 1
%     menu_subplots = menu('Plot each channel in a different subplot?', 'NO', 'YES (horizontal)', 'YES (vertical)');
end
if menu_subplots == 2 % horizontal
    layout = [1,2];
elseif menu_subplots == 3 % vertical
    layout = [2,1];
end

% plotting raw or smoothed channels
subplot(subplot_channels)
legend_A = {}; % initialising the legend cell
legend_B = {}; % initialising the legend cell
for i = files_to_average % files
    for j = channels_to_plot % channels
        if strfind(channel_name{i,j}, 'Channel A')
            if menu_subplots == 1
                subplot(subplot_channels)
                if size(channels_to_plot,2)>1
                    yyaxis left
                end
            elseif menu_subplots == 2
                subplot(layout(1),layout(2),1)
            end
            
            ylabel('Channel A (V)')
            grid on
            legend_A{end+1} = [file_name{i}(1:end-4) ...
                ' // ' channel_name{i,j} ...
                ];
%             pause(0.1)

        elseif strfind(channel_name{i,j}, 'Channel B')
            if menu_subplots == 1
                subplot(subplot_channels)
                if size(channels_to_plot,2)>1
                    yyaxis right
                end
            elseif menu_subplots == 2
                subplot(layout(1),layout(2),2)
            end
            
            ylabel('Channel B (V)')
            grid on
            legend_B{end+1} = [file_name{i}(1:end-4) ...
                ' // ' channel_name{i,j} ...
                ];
        end
        
        h(i,j) = plot(plot_data{i}(:,1), ...
            plot_data{i}(:,j),...
            'LineWidth', 1); hold all
        
%         pause(0.5)
    end
end

%% colour scheme
colour_type = {'DEFAULT', 'yellow' 'red', 'green', 'aqua', 'blue', 'purple', 'gray'};
for j = channels_to_plot
    menu_colour = 2;
%     menu_colour = menu([channel_name{1,j} ' colour scheme'], colour_type);
    for i = files_to_average   
        if menu_colour > 1 
            colour_RGB = colour_gradient(i, number_of_files, colour_type(menu_colour));
%             colour_RGB = colour_gradient(find(waveplate_angle_unique == waveplate_angle(i)), ...
%                 size(waveplate_angle_unique,2), ...
%                 colour_type(menu_colour));
            h(i,j).Color = colour_RGB;  
        end
        h(i,j).MarkerSize = 1;
        h(i,j).LineStyle = '-';
        h(i,j).LineWidth = 0.5;
    end
end

%% averaged
for j = channels_to_plot % channels
    if strfind(channel_name{i,j}, 'Channel A')
        if menu_subplots == 1
            subplot(subplot_channels)
            if size(channels_to_plot,2)>1
                yyaxis left
            end
        elseif menu_subplots == 2
            subplot(layout(1),layout(2),1)
        end

        ylabel('Channel A (V)')
        grid on
        legend_A{end+1} = 'average';
        
    elseif strfind(channel_name{i,j}, 'Channel B')
        if menu_subplots == 1
            subplot(subplot_channels)
            if size(channels_to_plot,2)>1
                yyaxis right
            end
        elseif menu_subplots == 2
            subplot(layout(1),layout(2),2)
        end

        ylabel('Channel B (V)')
        grid on
        legend_B{end+1} = 'average';
    end

    h_averaged(j) = plot(plot_data{1}(:,1), ...
        data_average(:,j),...
        'k', 'LineWidth', 2); hold all
    h_average.LineStyle = '-k';
    h_average.LineWidth = 2;
end


%% axis limits
% legend_location = 'NEO';
% x_limits = [-10,20];

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



%% SAVING DATA
% *************************************************************************
menu_save_average = 2;
% menu_save_average = menu('Save Averaged Data?', 'NO', 'YES');

if menu_save_average == 2
    file_name_save = file_name{1};
    file_name_save = strrep(file_name_save, '_01.txt', '_average.txt');
    file_name_save = strrep(file_name_save, '_1.txt', '_average.txt');
    [file_name_save,file_path_save,~] = uiputfile(['.' 'txt'],...
                'New File to Save the Averaged Data',[file_path file_name_save]); % choosing the file name
    file_identifier = fopen([file_path_save file_name_save], 'w');
    for j = 1:1:size(channel_name,2)
        fprintf(file_identifier, [channel_name{1,j} '\t']);
    end
%     fprintf(file_identifier, [channel_name{1,1} '\t']);
%     fprintf(file_identifier, [channel_name{1,2} '\t']);
%     fprintf(file_identifier, [channel_name{1,3} '\n']);
    fprintf(file_identifier, '\n');
    fprintf(file_identifier, '(ms)\t');
    for j = 2:1:size(channel_name,2)
        fprintf(file_identifier, '(V)\t');
    end
%     fprintf(file_identifier, '(V)\t');
%     fprintf(file_identifier, '(V)\n\n');
    fprintf(file_identifier, '\n\n');
    fclose(file_identifier);
    
    dlmwrite([file_path_save, file_name_save], data_average, '-append',...
        'delimiter', '\t', 'precision', '%.8f', 'newline', 'pc');
end


%% SAVING FIGURES
% *************************************************************************
menu_save_figures = 2;
menu_save_figures = menu('Save Figures?', 'NO', 'YES');
if menu_save_figures == 2    
    if findobj('tag','figure_picoscope') ~= 0
        figure_save = figure_picoscope;
%         file_name_save = file_name{1}(1:end-4);
%         file_name_save = strrep(file_name_save, '_01', '_average');
        file_name_save = file_name_save(1:end-4);
        
        figure(figure_save)
        pause(0.1)
%         file_path_save = file_path;
        [file_name_save,file_path_save,~] = uiputfile(['.' 'png'],...
            'File to Save the Figure',[file_path_save file_name_save]);
        hgexport(figure_save, [file_path_save file_name_save], hgexport('factorystyle'), 'Format', 'png')
        file_name_save = strrep(file_name_save, 'png', 'fig');    
        saveas(figure_save, [file_path_save file_name_save], 'fig');
        
    end
end

