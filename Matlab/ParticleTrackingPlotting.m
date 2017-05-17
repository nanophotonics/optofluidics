% Created by Ana (aa938)
% Reads csv files created with the Python particle tracking package

clc
clear
close all
figures = {};

%% CHOOSING FILES
% *************************************************************************

% specify default path
% folder_path = 'R:\aa938\NanoPhotonics\Laboratory\';
% folder_path = 'R:\3-Temporary\aa938\';
% folder_path = 'R:\3-Temporary\os354\';

folder_path = 'R:\aa938\NanoPhotonics\Laboratory\2017.05.12 - Omid particle tracking\';
file_name{1} = 'Substack_1_2999.csv';
file_path{1} = [folder_path file_name{1}];

number_of_folders = 0; % if 0 it will not ask for user input
% number_of_folders = menu('Where are the files located?', 'SINGLE folder', 'MULTIPLE folders');

if number_of_folders == 1
    % pop up window to choose the file(s) to read from a SINGLE FOLDER
    [file_name, folder_path, ~] = uigetfile('.csv',...
                                          'PicoScope Files to Read (use CTRL to select multiple files)',...
                                          folder_path,...
                                          'MultiSelect','on');
    file_name = cellstr(file_name); % convert to cell array of strings
    file_path = cell(size(file_name));
    for i = 1:1:size(file_name,2)
        file_path{i} = [folder_path file_name{i}];
    end
elseif number_of_folders == 2
    % pop up window to choose the file(s) to read from MULTIPLE FOLDERS
    paths = uipickfiles('FilterSpec', [folder_path '*.csv']);
    file_path = {};
    for i = 1:1:size(paths,2)
        if strfind(paths{i}, 'csv')
            file_path{end+1} = paths{i};
        else
            directory_contents = dir(paths{i});
            for j = 3:1:size(directory_contents,1)
                file_path{end+1} = [paths{i} '\' directory_contents(j).name];
                if strfind(file_path{end}, 'csv')
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

%% READING DATA
% *************************************************************************

% reading the data from the csv files
raw_data = cell(size(file_name));
for i = 1:1:number_of_files
    disp(['Reading File ' num2str(i) '/' num2str(number_of_files)])
    raw_data{i} = readtable(file_path{i});    
end
disp('Finished reading all files!')
% disp(raw_data{i}.Properties.VariableNames)
% size: units of pixels
% mass: brightness of the blog

%% SORT DATA BY PARTICLES
% *************************************************************************
clear particle
i = 1;
particle_numbers = unique(raw_data{i}.particle)';
% particle_numbers = 1:1:20;
particle = zeros(numel(particle_numbers),3);
for j = 1:1:numel(particle_numbers)
    indices = find(raw_data{i}.particle == particle_numbers(j));
    particle(j,1) = particle_numbers(j); % number
    particle(j,2) = mean(raw_data{i}.mass(indices)); % mass
    particle(j,3) = std(raw_data{i}.mass(indices)); % mass standard deviation
    particle(j,4) = mean(raw_data{i}.size(indices)); % size
    particle(j,5) = std(raw_data{i}.size(indices)); % size standard deviation
    particle(j,6) = numel(indices); % number of frames
end
particle = array2table(particle);
particle.Properties.VariableNames = {'number','mass','mass_std','size','size_std','frames'};

figures{end+1} = figure;
errorbar(particle.number, particle.mass, particle.mass_std, '.')
xlabel('particle number')
ylabel('mass')
title(file_name, 'interpreter', 'none')
set(gca, 'FontSize', 16)

figures{end+1} = figure;
histogram(raw_data{i}.mass), hold all
histogram(particle.mass,20), hold all
xlabel('mass')
ylabel('counts')
title(file_name, 'interpreter', 'none')
set(gca, 'FontSize', 16)
legend('all data points', 'one point per particle')

figures{end+1} = figure;
histogram(raw_data{i}.size), hold all
histogram(particle.size,20), hold all
xlabel('size')
ylabel('counts')
title(file_name, 'interpreter', 'none')
set(gca, 'FontSize', 16)
legend('all data points', 'one point per particle')

figures{end+1} = figure;
plot(particle.mass, particle.size, '.', 'MarkerSize', 8)
ylabel('size')
xlabel('mass')
title(file_name, 'interpreter', 'none')
set(gca, 'FontSize', 16)

figures{end+1} = figure;
histogram(particle.frames, 50)
xlabel('frames per particle')
ylabel('counts')
title(file_name, 'interpreter', 'none')
set(gca, 'FontSize', 16)

%% CALCULATE VELOCTIY
% *************************************************************************
% close all

i = 1;
x = [];
y = [];
vx = [];
vy = [];
frame_rate = 159.22; % fps
for j = 2:1:size(raw_data{i},1)
%     if raw_data{i}.mass(j) > 1.5e5 && raw_data{i}.size(j) > 15
    index = find(particle.number == raw_data{i}.particle(j));
    if particle.frames(index) > 200
        disp(raw_data{i}.particle(j))
        disp(index)
        disp(particle.frames(index))
        if raw_data{i}.particle(j) == raw_data{i}.particle(j-1)
            x(end+1) = raw_data{i}.x(j);
            y(end+1) = raw_data{i}.y(j);
            vx(end+1) = (raw_data{i}.x(j) - raw_data{i}.x(j-1)) / ...
                (raw_data{i}.frame(j) - raw_data{i}.frame(j-1)); % pixels / frame
            vy(end+1) = (raw_data{i}.y(j) - raw_data{i}.y(j-1)) / ...
                (raw_data{i}.frame(j) - raw_data{i}.frame(j-1)); % pixels / frame
        end
    end
    % the number of particles is equal to size(raw_data{i},1) - size(x,2)
end
% velocity = table(x,y,vx,vy);
% velocity.Properties.VariableNames = {'x','y','vx','vy'};

figures{end+1} = figure;
histogram(vx,30), hold all
histogram(vy,30), hold all
legend('vx','vy')
xlabel('velocity (pixels/frame)')
ylabel('counts')
title(file_name, 'interpreter', 'none')
set(gca, 'FontSize', 16)

figures{end+1} = figure;
subplot(2,1,1)
plot(x,vx, '.-', 'MarkerSize', 8), hold all
plot(x,vy, '.-', 'MarkerSize', 8), hold all
legend('vx vs. x', 'vy vs. x')
xlabel('x position (pixels)')
ylabel('velocity (pixels/frame)')
title(file_name, 'interpreter', 'none')
set(gca, 'FontSize', 16)

subplot(2,1,2)
plot(y,vy, '.-', 'MarkerSize', 8), hold all
plot(y,vx, '.-', 'MarkerSize', 8), hold all
legend('vy vs. y', 'vx vs. y')
xlabel('y position (pixels)')
ylabel('velocity (pixels/frame)')
set(gca, 'FontSize', 16)

figures{end+1} = figure;
plot(x,y, '.-', 'MarkerSize', 8), hold all
xlabel('x position (pixels)')
ylabel('y position (pixels)')
set(gca, 'FontSize', 16)

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

