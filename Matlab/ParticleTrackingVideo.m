% Created by Ana (aa938)
% Reads csv files and creates a video to test the particle tracking software

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

% folder_path = 'R:\aa938\NanoPhotonics\Matlab\Particle Tracking\';
% file_name{1} = 'test_1.csv';

file_names{1} = '';
folder_path = 'R:\aa938\NanoPhotonics\Laboratory\2017.05.21 - Omid particle tracking\';
% file_names{1} = '780_Substack_0_1999_linked.csv';
% file_names{2} = '780_Substack_1_12236.tif';

for i = 1:1:numel(file_names)
    file_path{i} = [folder_path file_names{i}];
end

number_of_folders = 1; % if 0 it will not ask for user input
% number_of_folders = menu('Where are the files located?', 'SINGLE folder', 'MULTIPLE folders');

if number_of_folders == 1
    % pop up window to choose the file(s) to read from a SINGLE FOLDER
    [file_names, folder_path, ~] = uigetfile('*.*',...
                                          'Files to Read (use CTRL to select multiple files)',...
                                          folder_path,...
                                          'MultiSelect','on');
    file_names = cellstr(file_names); % convert to cell array of strings
    file_path = cell(size(file_names));
    for i = 1:1:size(file_names,2)
        file_path{i} = [folder_path file_names{i}];
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
    file_names = cell(size(file_path));
    for i = 1:1:size(file_path,2)
        slash = strfind(folder_path, '\');
        file_names{i} = file_path{i}(slash(end)+1:end);
    end
end
number_of_files = size(file_names,2);

% sort files by type
csv_files = [];
tif_files = [];
for i = 1:1:number_of_files
    if strfind(file_names{i}, '.csv')
        csv_files(end+1) = i;
    elseif strfind(file_names{i}, '.tif')
        tif_files(end+1) = i;
    end
end


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
raw_data = cell(size(file_names));
for i = csv_files
    raw_data{i} = readtable(file_path{i});    
end
disp(raw_data{i}.Properties.VariableNames)
% size: units of pixels
% mass: brightness of the blog

black_frame = imread('R:\aa938\NanoPhotonics\Matlab\Particle Tracking\black.png');

%% filter data
% *************************************************************************
selected_filter = [1];
filter_options = {'Frame', 'Particle'};
% [selected_filter, ~] = listdlg('PromptString', 'Filter data by:',...
%                            'SelectionMode', 'single', ...
%                            'ListString', filter_options,...
%                            'InitialValue', selected_filter);

selected_data = raw_data{csv_files(1)};

if find(strcmp(filter_options(selected_filter), 'Frame'))
    desired_frames = 1:1:5999;
    desired_frames_indices = [];
    for j = desired_frames
        desired_frames_indices = [desired_frames_indices; ...
            find(raw_data{csv_files(1)}.frame == j)];
    end
else
    desired_frames_indices = 1:1:size(raw_data{csv_files(1)},1);
end

if find(strcmp(filter_options(selected_filter), 'Particle'))
    desired_particles = 1:100;
    desired_particles_indices = [];
    for j = desired_particles
        desired_particles_indices = [desired_particles_indices; ...
            find(raw_data{csv_files(1)}.particle == j)];
    end
else
    desired_particles_indices = 1:1:size(raw_data{csv_files(1)},1);
end

desired_indices = intersect(desired_frames_indices, desired_particles_indices);
selected_data = raw_data{csv_files(1)}(desired_indices,:);  

desired_frames = unique(selected_data.frame)';
desired_particles = unique(selected_data.particle);
desired_colours = prism(numel(desired_particles));

%% sort data
% *************************************************************************
frames = cell(size(desired_frames));
for j = desired_frames
    desired_frames_indices = find(selected_data.frame == j);
    particle = selected_data.particle(desired_frames_indices); % particle numbers
    colour = desired_colours(ismember(desired_particles,particle),:);
    x = selected_data.x(desired_frames_indices); % x positions
    y = selected_data.y(desired_frames_indices); % y positions
    size = selected_data.size(desired_frames_indices); % sizes
    frames{j} = table(particle,x,y,size,colour);
end


%% CREATE VIDEO FRAMES
% *************************************************************************

close all
% the figure size doesn't work with imshow
figure('Units','normalized','Position',[0.01 0.07 0.95 0.8]);
% imshow(black_frame), hold off
imshow(imread(file_path{tif_files(1)}, 'tif', 'Index', 2)), hold on
video_frames = getframe;
for j = desired_frames
    clc
    imshow(imread(file_path{tif_files(1)}, 'tif', 'Index', j+1)), hold on
    if find(desired_frames == j)
        for k = 1:1:height(frames{j})
            plot(frames{j}.x(k), frames{j}.y(k), 'o', ...
                'MarkerSize',frames{j}.size(k), ...
                'LineWidth', 1, 'Color', frames{j}.colour(k,:))
            text('Position',[frames{j}.x(k),frames{j}.y(k)-frames{j}.size(k)*3], ...
                 'FontSize', 7, 'VerticalAlignment', 'middle', ...
                 'HorizontalAlignment', 'center', ...
                 'Color', frames{j}.colour(k,:), ...
                 'String' , num2str(frames{j}.particle(k)))
        end
    end
    title([num2str(j/max(desired_frames)*100,'%02.0f') '%: frame ' ...
        num2str(j) ' of ' num2str(max(desired_frames))])
    hold off
    video_frames(end+1) = getframe;
end

%% WRITE VIDEO
% *************************************************************************
menu_save_video = 1;
menu_save_video = menu('Save Video?', 'NO', 'YES: avi');
% menu_save_video = menu('Save Video?', 'NO', 'YES: avi', 'YES: tif');
if menu_save_video == 2 
    video_path_save = file_path{2};
    video_path_save = strrep(video_path_save,'test','video');
    video_path_save = strrep(video_path_save,'csv','avi');
    video_path_save = strrep(video_path_save,'tif','avi');
    
    [video_name_save,folder_path_save,~] = uiputfile(['.' 'avi'],...
                'File to Save the Video',video_path_save);
    video_path_save = [folder_path_save video_name_save];
    
    v = VideoWriter(video_path_save, 'Uncompressed AVI');
    open(v)
    writeVideo(v,video_frames) 
    close(v)
elseif menu_save_video == 3 % DOESN'T WORK YET
    video_path_save = file_path{tif_files(1)};    
    [video_name_save,folder_path_save,~] = uiputfile(['.' 'tif'],...
                'File to Save the Video',video_path_save);
    video_path_save = [folder_path_save video_name_save];
    imwrite(video_frames(1).cdata, video_path_save)
    for k = 2:size(stack,3)
        imwrite(stack(:,:,k), video_path_save, 'writemode', 'append');
    end
end

disp('Finished saving video :)')