% Created by Ana (aa938)
% Reads txt files exported with the Malvern Zetaseizer software

clc
clear
close all
figures = {};

%% choosing files
% *************************************************************************

% specify default path
folder_path = 'R:\aa938\NanoPhotonics\Laboratory\2017.06.13 - Au NR DLS\';
file_name{1} = '2017.06.13_AuNR_1_diffusions.txt';

file_path{1} = [folder_path file_name{1}];

number_of_folders = 1; % if 0 it will not ask for user input
% number_of_folders = menu('Where are the files located?', 'SINGLE folder', 'MULTIPLE folders');

if number_of_folders == 1
    % pop up window to choose the file(s) to read from a SINGLE FOLDER
    [file_name, folder_path, ~] = uigetfile('.txt',...
                                          'Files to Read (use CTRL to select multiple files)',...
                                          folder_path,...
                                          'MultiSelect','on');
    file_name = cellstr(file_name); % convert to cell array of strings
    file_path = cell(size(file_name));
    for i = 1:1:size(file_name,2)
        file_path{i} = [folder_path file_name{i}];
    end
elseif number_of_folders == 2
    % pop up window to choose the file(s) to read from MULTIPLE FOLDERS
    paths = uipickfiles('FilterSpec', [folder_path '*.txt']);
    file_path = {};
    for i = 1:1:size(paths,2)
        if strfind(paths{i}, 'txt')
            file_path{end+1} = paths{i};
        else
            directory_contents = dir(paths{i});
            for j = 3:1:size(directory_contents,1)
                file_path{end+1} = [paths{i} '\' directory_contents(j).name];
                if strfind(file_path{end}, 'txt')
                else % delete all non-csv files
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
    folder_name = folder_path(strfind(folder_path, 'Laboratory')+11:slash(slash_index(1))-1);
elseif strfind(folder_path, 'aa938')
    slash = strfind(folder_path, '\');
    slash_index = find(slash > strfind(folder_path, 'aa938')+6);
    folder_name = folder_path(strfind(folder_path, 'aa938')+6:slash(slash_index(1))-1);
else
    folder_name = folder_path;
end  

file_path_save = folder_path;

%% reading data
% *************************************************************************

% reading the data from the csv files
raw_data = cell(size(file_name));
intensities = [];
diffusions = [];
sizes = [];
T = [];
viscosity = [];
for i = 1:1:number_of_files
    disp(['Reading File ' num2str(i) '/' num2str(number_of_files)])
    raw_data{i} = readtable(file_path{i}, 'ReadVariableNames', 1, 'Delimiter', '\t');    
    
    index = strfind(raw_data{i}.Properties.VariableNames, 'Intensities');
    intensities = [intensities; table2array(raw_data{i}(:,not(cellfun('isempty', index))))];

    index = strfind(raw_data{i}.Properties.VariableNames, 'Diffusions');
    diffusions = [diffusions; table2array(raw_data{i}(:,not(cellfun('isempty', index))))];

    index = strfind(raw_data{i}.Properties.VariableNames, 'Sizes');
    sizes = [sizes; table2array(raw_data{i}(:,not(cellfun('isempty', index))))];
    
    index = strfind(raw_data{i}.Properties.VariableNames, 'Temperature');
    T = [T; table2array(raw_data{i}(:,not(cellfun('isempty', index))))];
    
    index = strfind(raw_data{i}.Properties.VariableNames, 'Viscosity');
    viscosity = [viscosity; table2array(raw_data{i}(:,not(cellfun('isempty', index))))];
end
disp('Finished reading all files!')
% disp(raw_data{1}.Properties.VariableNames)



%% plotting data
% *************************************************************************

figures{end+1} = figure('Units','normalized','Position',[0.1 0.09 0.8 0.8]);
p = cell(size(intensities,1),2);
p_legend = [];
for i = 1:1:size(intensities,1)
    subplot(2,1,1)
    p{i,1} = semilogx(diffusions(i,:), intensities(i,:), '.-'); hold all
    subplot(2,1,2)
    p{i,2} = semilogx(sizes(i,:), intensities(i,:), '.-'); hold all
    
    p_legend{end+1} = ['T = ' num2str(T(i)) '\circC'];
    p_legend{end} = [p_legend{end} ', \eta = ' num2str(viscosity(i)) ' cP'];
end

subplot(2,1,1)
title(file_name{1}, 'interpreter', 'none')
xlabel('Diffusion Coefficient (\mum^2/s)')
ylabel('Intensity (%)')
set(gca,'FontSize',16)
grid on
legend(p_legend, 'Location', 'EO')

subplot(2,1,2)
xlabel('Size (nm)')
ylabel('Intensity (%)')
set(gca,'FontSize',16)
grid on
legend(p_legend, 'Location', 'EO')

%% fitting data
% *************************************************************************


%% SAVING FIGURES
% *************************************************************************
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