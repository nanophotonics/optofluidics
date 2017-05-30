% Created by Ana (aa938)
% Reads csv files created with the Python particle tracking package

clc
clear
close all
figures = {};

colour_type = {'DEFAULT', ...
               'parula', 'jet', 'hsv', 'cool', ...
               'spring', 'summer', 'autumn', 'autumn reversed', 'winter', ...
               'gray', 'copper',...
               'red', 'green', 'aqua', 'blue', 'purple',...
               };
selected_colour = 2;


%% choosing files
% *************************************************************************

% specify default path
% folder_path = 'R:\aa938\NanoPhotonics\Laboratory\';
% folder_path = 'R:\3-Temporary\aa938\';
% folder_path = 'R:\3-Temporary\os354\';

% folder_path = 'R:\aa938\NanoPhotonics\Laboratory\2017.05.12 - Omid particle tracking\';
% file_name{1} = 'Substack_1_2999.csv';

folder_path = 'R:\aa938\NanoPhotonics\Laboratory\2017.05.21 - Omid particle tracking\';
% file_name{1} = '780_Substack_0_1999_linked.csv';

% folder_path = 'R:\os354\Year 1\Lab Data\20170530 - Linked Trajectories for 780, 800, 850 900\';
file_name{1} = '';

file_path{1} = [folder_path file_name{1}];

number_of_folders = 1; % if 0 it will not ask for user input
% number_of_folders = menu('Where are the files located?', 'SINGLE folder', 'MULTIPLE folders');

if number_of_folders == 1
    % pop up window to choose the file(s) to read from a SINGLE FOLDER
    [file_name, folder_path, ~] = uigetfile('.csv',...
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

%% reading data
% *************************************************************************

% reading the data from the csv files
raw_data = cell(size(file_name));
for i = 1:1:number_of_files
    disp(['Reading File ' num2str(i) '/' num2str(number_of_files)])
    raw_data{i} = readtable(file_path{i});    
end
disp('Finished reading all files!')
disp(raw_data{i}.Properties.VariableNames)
% size: units of pixels
% mass: brightness of the blog

%% sort data by particles
% *************************************************************************
% close all

clear particles
i = 1;
particle_numbers = unique(raw_data{i}.particle)';
% particle_numbers = 20;
particles = zeros(numel(particle_numbers),3);
for j = 1:1:numel(particle_numbers)
    indices = find(raw_data{i}.particle == particle_numbers(j));
    particles(j,1) = particle_numbers(j); % number
    particles(j,2) = mean(raw_data{i}.mass(indices)); % mass
    particles(j,3) = std(raw_data{i}.mass(indices)); % mass standard deviation
    particles(j,4) = mean(raw_data{i}.size(indices)); % size
    particles(j,5) = std(raw_data{i}.size(indices)); % size standard deviation
    particles(j,6) = numel(indices); % number of frames
    particles(j,7) = max(raw_data{i}.x(indices))-min(raw_data{i}.x(indices));
end
particles = array2table(particles);
particles.Properties.VariableNames = {'number','mass','mass_std','size','size_std','frames','x_travelled'};

% figures{end+1} = figure;
% errorbar(particles.number, particles.mass, particles.mass_std, '.')
% xlabel('particle number')
% ylabel('mass')
% title([folder_name ': ' file_name{1}], 'interpreter', 'none')
% set(gca, 'FontSize', 16)

% figures{end+1} = figure;
% histogram(raw_data{i}.mass,'Normalization','probability'), hold all
% histogram(particles.mass,20,'Normalization','probability'), hold all
% xlabel('mass')
% ylabel('probability')
% title([folder_name ': ' file_name{1}], 'interpreter', 'none')
% set(gca, 'FontSize', 16)
% legend([num2str(size(raw_data{1},1)) ' points'], [num2str(size(particles,1)) ' particles'])

% figures{end+1} = figure;
% histogram(raw_data{i}.size,'Normalization','probability'), hold all
% histogram(particles.size,20,'Normalization','probability'), hold all
% xlabel('particle size (px)')
% ylabel('probability')
% title([folder_name ': ' file_name{1}], 'interpreter', 'none')
% set(gca, 'FontSize', 16)
% legend([num2str(size(raw_data{1},1)) ' points'], [num2str(size(particles,1)) ' particles'])

% figures{end+1} = figure;
% plot(raw_data{i}.mass, raw_data{i}.size, '.', 'MarkerSize', 8), hold all
% plot(particles.mass, particles.size, '.', 'MarkerSize', 8), hold all
% ylabel('particle size (px)')
% xlabel('mass')
% title([folder_name ': ' file_name{1}], 'interpreter', 'none')
% set(gca, 'FontSize', 16)
% legend([num2str(size(raw_data{1},1)) ' points'], [num2str(size(particles,1)) ' particles'])

figures{end+1} = figure;
suptitle([folder_name ': ' file_name{1}])
subplot(1,2,1)
plot(raw_data{i}.size, raw_data{i}.mass, '.', 'MarkerSize', 8), hold all
plot(particles.size, particles.mass, '.', 'MarkerSize', 8), hold all
xlabel('particle size (px)')
ylabel('mass')
% title([folder_name ': ' file_name{1}], 'interpreter', 'none')
set(gca, 'FontSize', 12)
legend([num2str(size(raw_data{1},1)) ' points'], [num2str(size(particles,1)) ' particles'], 'Location', 'NW')

% figures{end+1} = figure;
% histogram(raw_data{i}.mass./raw_data{i}.size,'Normalization','probability'), hold all
% histogram(particles.mass./particles.size,'Normalization','probability'), hold all
% xlabel('particle mass / particle size (px)')
% ylabel('probability')
% title([folder_name ': ' file_name{1}], 'interpreter', 'none')
% set(gca, 'FontSize', 16)
% legend([num2str(size(raw_data{1},1)) ' points'], [num2str(size(particles,1)) ' particles'])

% figures{end+1} = figure;
% histogram(raw_data{i}.size./raw_data{i}.mass,'Normalization','probability'), hold all
% histogram(particles.size./particles.mass,'Normalization','probability'), hold all
% xlabel('particle size (px) / particle mass')
% ylabel('probability')
% title([folder_name ': ' file_name{1}], 'interpreter', 'none')
% set(gca, 'FontSize', 16)
% legend([num2str(size(raw_data{1},1)) ' points'], [num2str(size(particles,1)) ' particles'])

% figures{end+1} = figure;
% subplot(2,2,4)
% histogram(particles.frames, 50,'Normalization','probability')
% xlabel('frames per particle')
% ylabel('probability')
% title([folder_name ': ' file_name{1}], 'interpreter', 'none')
% set(gca, 'FontSize', 16)

% figures{end+1} = figure;
% histogram(particles.x_travelled,40,'Normalization','probability'), hold all
% xlabel('x travelled (px)')
% ylabel('probability')
% title([folder_name ': ' file_name{1}], 'interpreter', 'none')
% set(gca, 'FontSize', 16)
% legend([num2str(size(particles,1)) ' particles'])

% figures{end+1} = figure;
% errorbar(particles.x_travelled, particles.mass, particles.mass_std, '.')
% plot(particles.x_travelled, particles.mass, '.')
% xlabel('x travelled (px)')
% ylabel('mass')
% title([folder_name ': ' file_name{1}], 'interpreter', 'none')
% set(gca, 'FontSize', 16)

% figures{end+1} = figure;
% errorbar(particles.x_travelled, particles.size, particles.size_std, '.')
% plot(particles.x_travelled, particles.size, '.')
% xlabel('x travelled (px)')
% ylabel('particle size (px)')
% title([folder_name ': ' file_name{1}], 'interpreter', 'none')
% set(gca, 'FontSize', 16)

% figures{end+1} = figure;
subplot(1,2,2)
plot(particles.frames, particles.x_travelled, '.')
ylabel('x travelled (px)')
xlabel('frames')
% title([folder_name ': ' file_name{1}], 'interpreter', 'none')
set(gca, 'FontSize', 12)
legend([num2str(size(particles,1)) ' particles'])

% figures{end+1} = figure;
% plot(raw_data{i}.y, raw_data{i}.mass, '.', 'MarkerSize', 8), hold all
% xlabel('y particle position (px)')
% ylabel('mass')
% title([folder_name ': ' file_name{1}], 'interpreter', 'none')
% set(gca, 'FontSize', 16)
% legend([num2str(size(raw_data{1},1)) ' points'])


%% units
% *************************************************************************

time_units = {'frame', 's'};
position_units = {'px', '\mum'};
t = 2;
p = 2;

% [t, p] = dialog_two_lists('Select plot options:', ...
%                                       'Units of time:', time_units, t,...
%                                       'Units of position:', position_units, p);

frame_rate = 1;
pixel_scale = 1;

if find(strcmp(time_units(t), 's'))
    frame_rate = 159.22; % frames / second
%     input_title = 'Video parameters';
%     input_data = {'Frame rate (fps):'};
%     default_values = {num2str(frame_rate)};
%     dlg_options.WindowStyle = 'normal'; dlg_options.Resize = 'on'; dim = [1 80];
%     answer = inputdlg(input_data,input_title,dim,default_values,dlg_options);
%     frame_rate = str2double(answer{1}); % fps
end

if find(strcmp(position_units(p), '\mum'))
    pixel_scale = 0.043; % microns / pixel
%     input_title = 'Video parameters';
%     input_data = {'Pixel scale (microns / pixel):'};
%     default_values = {num2str(pixel_scale)};
%     dlg_options.WindowStyle = 'normal'; dlg_options.Resize = 'on'; dim = [1 80];
%     answer = inputdlg(input_data,input_title,dim,default_values,dlg_options);
%     pixel_scale = str2double(answer{1}); % micron / pixel
end

%% CALCULATE velocity
% *************************************************************************

i = 1;
x = [];
y = [];
masses = [];
sizes = [];
vx = [];
vy = [];
for j = 2:1:size(raw_data{i},1)
% for j = 2:1:18400
    index = find(particles.number == raw_data{i}.particle(j));
%     if raw_data{i}.mass(j) > 1.5e5 && raw_data{i}.size(j) > 15
%     if particles.frames(index) > 200
%         disp(raw_data{i}.particle(j))
%         disp(index)
%         disp(particles.frames(index))
%     if particles.x_travelled(index) > 20
        if raw_data{i}.particle(j) == raw_data{i}.particle(j-1)
            x(end+1) = raw_data{i}.x(j) * pixel_scale;
            y(end+1) = raw_data{i}.y(j) * pixel_scale;
            masses(end+1) = raw_data{i}.mass(j);
            sizes(end+1) = raw_data{i}.size(j);
            vx(end+1) = (raw_data{i}.x(j) - raw_data{i}.x(j-1)) / ...
                (raw_data{i}.frame(j) - raw_data{i}.frame(j-1)) * ...
                pixel_scale * frame_rate;
            vy(end+1) = (raw_data{i}.y(j) - raw_data{i}.y(j-1)) / ...
                (raw_data{i}.frame(j) - raw_data{i}.frame(j-1)) *...
                pixel_scale * frame_rate;
%         end
    end
    % the number of particles is equal to size(raw_data{i},1) - size(x,2)
end
% velocity = table(x,y,vx,vy);
% velocity.Properties.VariableNames = {'x','y','vx','vy'};

%% CALCULATE velocity maps
% *************************************************************************
menu_velocity = 1;
% menu_velocity = menu('Velocity:', 'Average', 'Maximum', 'Minimum', 'Mode');

x_sections = 19;
y_sections = 11;

x_edges = linspace(min(x), max(x), x_sections + 2);
x_bin = discretize(x,x_edges);
y_edges = linspace(min(y), max(y), y_sections + 2);
y_bin = discretize(y,y_edges);

masses_map = zeros(numel(x_edges)-1,numel(y_edges)-1);
sizes_map = zeros(numel(x_edges)-1,numel(y_edges)-1);
vx_grid = cell(numel(x_edges)-1,numel(y_edges)-1);
vx_map = zeros(numel(x_edges)-1,numel(y_edges)-1);
vy_map = zeros(numel(x_edges)-1,numel(y_edges)-1);
for i = 1:1:numel(x_edges)-1
    x_indices = find(x_bin == i);
    for j = 1:1:numel(y_edges)-1
        y_indices = find(y_bin == j);
        xy_indices = intersect(x_indices, y_indices);
%         disp(numel(xy_indices))
        
        if menu_velocity == 1 % average
            vx_grid{i,j} = vx(xy_indices);
            vx_map(i,j) = mean(vx(xy_indices));
            vy_map(i,j) = mean(vy(xy_indices));
            masses_map(i,j) = mean(masses(xy_indices));
            sizes_map(i,j) = mean(sizes(xy_indices));
            v_label = 'average';
            
        elseif menu_velocity == 2 % maximum
            if isempty(vx(xy_indices))
                vx_map(i,j) = NaN;
            else           
                vx_map(i,j) = max(vx(xy_indices));
            end
            if isempty(vy(xy_indices))
                vy_map(i,j) = NaN;
            else           
                vy_map(i,j) = max(vy(xy_indices));
            end
            v_label = 'max';
            
        elseif menu_velocity == 3 % mimimum
            if isempty(vx(xy_indices))
                vx_map(i,j) = NaN;
            else           
                vx_map(i,j) = min(vx(xy_indices));
            end
            if isempty(vy(xy_indices))
                vy_map(i,j) = NaN;
            else           
                vy_map(i,j) = min(vy(xy_indices));
            end
            v_label = 'min';
            
        elseif menu_velocity == 4 % mode
            vx_map(i,j) = mode(round(vx(xy_indices),0));
            vy_map(i,j) = mode(round(vy(xy_indices),0));
            v_label = 'mode';
        end
        
    end
end

%% PLOT velocity histogram
% *************************************************************************
% close all

figures{end+1} = figure;
legend_hv = {};
text_hv = {};

number_of_bins = 100;
% number_of_bins = round(size(particles,1)/10);

hvx = histogram(vx,number_of_bins,'Normalization','probability'); hold all, legend_hv{end+1} = 'vx';
hvy = histogram(vy,number_of_bins,'Normalization','probability'); hold all, legend_hv{end+1} = 'vy';

text_hv{end+1} = 'Gaussian fit:';
text_hv{end+1} = 'y =  a1*exp(-((x-b1)/c1)^2)';
text_hv{end+1} = '';

hvx_x = hvx.BinEdges(1:end-1) + hvx.BinWidth/2;
hvx_y = hvx.Values;
fit_vx = fit(hvx_x', hvx_y', 'gauss1');
confint_vx = confint(fit_vx);
pvx = plot(fit_vx); hold all, legend_hv{end+1} = 'vx fit';
pvx.Color = 'b';
pvx.LineWidth = 2;
% text_hv{end+1} = ['vx mean = ' num2str(mean(vx)) ' ' position_units{p} '/' time_units{t}];
text_hv{end+1} = ['vx: a1 = ' num2str(fit_vx.a1, '%.4f') ' \pm ' ...
    num2str(abs(confint_vx(1,1)-confint_vx(2,1))/2, '%.4f') ...
    ' '];
text_hv{end+1} = ['vx: b1 = ' num2str(fit_vx.b1, '%.1f') ' \pm ' ...
    num2str(abs(confint_vx(1,2)-confint_vx(2,2))/2, '%.1f') ' ' ...
    position_units{p} '/' time_units{t}];
text_hv{end+1} = ['vx: c1 = ' num2str(fit_vx.c1, '%.1f') ' \pm ' ...
    num2str(abs(confint_vx(1,3)-confint_vx(2,3))/2, '%.1f')  ' ' ...
    position_units{p} '/' time_units{t}];
text_hv{end+1} = '';

hvy_x = hvy.BinEdges(1:end-1) + hvy.BinWidth/2;
hvy_y = hvy.Values;
fit_vy = fit(hvy_x', hvy_y', 'gauss1');
confint_vy = confint(fit_vy);
pvy = plot(fit_vy); hold all, legend_hv{end+1} = 'vy fit';
pvy.Color = 'r';
pvy.LineWidth = 2;
% text_hv{end+1} = ['vy mean = ' num2str(mean(vy)) ' ' position_units{p} '/' time_units{t}];
text_hv{end+1} = ['vy: a1 = ' num2str(fit_vy.a1, '%.4f') ' \pm ' ...
    num2str(abs(confint_vy(1,1)-confint_vy(2,1))/2, '%.4f') ...
    ' '];
text_hv{end+1} = ['vy: b1 = ' num2str(fit_vy.b1, '%.1f') ' \pm ' ...
    num2str(abs(confint_vy(1,2)-confint_vy(2,2))/2, '%.1f') ...
    position_units{p} '/' time_units{t}];
text_hv{end+1} = ['vy: c1 = ' num2str(fit_vy.c1, '%.1f') ' \pm ' ...
    num2str(abs(confint_vy(1,3)-confint_vy(2,3))/2, '%.1f') ...
    position_units{p} '/' time_units{t}];

legend(legend_hv)
xlabel(['velocity (' position_units{p} '/' time_units{t} ')'])
ylabel('probability')
title([folder_name ': ' file_name{1}], 'interpreter', 'none')
set(gca, 'FontSize', 16)


text('Units','normalized','Position',[0.08,0.95], ...
    'FontSize', 12, 'VerticalAlignment', 'top', 'String' , text_hv)


%% PLOT velocity histogram by sections
% *************************************************************************
% close all

if max(size(vx_map)) < 5
    figures{end+1} = figure;
    legend_vx_grid = {};
    text_hv = {};
    text_hv{end+1} = 'Gaussian fit:';
    text_hv{end+1} = 'y =  a1*exp(-((x-b1)/c1)^2)';
    text_hv{end+1} = '';
    
    grid_colours = parula(size(vx_grid,1)*size(vx_grid,2));
    k = 1;
    vx_map_centre = zeros(size(vx_grid));
    for i = 1:1:size(vx_grid,1)
        for j = 1:1:size(vx_grid,2)
            ph(i,j) = histogram(vx_grid{i,j},30); hold all
%             ph(i,j).Normalization = 'count';
            ph(i,j).Normalization = 'probability';
            ph(i,j).FaceColor = grid_colours(k,:);
            legend_vx_grid{end+1} = ['x = ' num2str(x_edges(i),'%02.0f') ' ' position_units{p} ...
                ', y = ' num2str(y_edges(j),'%02.0f') ' ' position_units{p}];
            
            hvx_x = ph(i,j).BinEdges(1:end-1) + ph(i,j).BinWidth/2;
            hvx_y = ph(i,j).Values;
            fit_vx_grid{i,j} = fit(hvx_x', hvx_y', 'gauss1');
            confint_vx_grid{i,j} = confint(fit_vx_grid{i,j});
            vx_map_centre(i,j) = fit_vx_grid{i,j}.b1;
            
            pvx_grid{i,j} = plot(fit_vx_grid{i,j}); hold all
            pvx_grid{i,j}.Color = ph(i,j).FaceColor;
            pvx_grid{i,j}.LineWidth = 2;
            
            legend_vx_grid{end+1} = [num2str(k) ': Gausian fit'];
            
            text_hv{end+1} = [num2str(k) ': a1 = ' num2str(fit_vx_grid{i,j}.a1, '%.2g') ' \pm ' ...
                num2str(abs(confint_vx_grid{i,j}(1,1)-confint_vx_grid{i,j}(2,1))/2, '%.2g') ...
                ' '];
            text_hv{end+1} = [num2str(k) ': b1 = ' num2str(fit_vx_grid{i,j}.b1, '%.1f') ' \pm ' ...
                num2str(abs(confint_vx_grid{i,j}(1,2)-confint_vx_grid{i,j}(2,2))/2, '%.1f') ...
                position_units{p} '/' time_units{t}];
            text_hv{end+1} = [num2str(k) ': c1 = ' num2str(fit_vx_grid{i,j}.c1, '%.1f') ' \pm ' ...
                num2str(abs(confint_vx_grid{i,j}(1,3)-confint_vx_grid{i,j}(2,3))/2, '%.1f') ...
                position_units{p} '/' time_units{t}];
            text_hv{end+1} = '';
            k = k + 1;
        end
    end
    xlabel(['x velocity (' position_units{p} '/' time_units{t} ')'])
    ylabel(ph(i,j).Normalization)
    legend(legend_vx_grid)
    set(gca, 'FontSize', 16)
    title([folder_name ': ' file_name{1}], 'interpreter', 'none')
    text('Units','normalized','Position',[0.08,0.95], ...
    'FontSize', 12, 'VerticalAlignment', 'top', 'String' , text_hv)

    % for i = 1:1:size(vx_map,2)
    %     h{1}(i) = plot(x_edges(1:end-1),vx_map(:,i), '.-', 'MarkerSize', 8); hold all
    %     legend_cell{1}{end+1} = ['y = ' num2str(y_edges(i),'%02.0f') ' ' position_units{p}];
    % end
    % xlabel(['x position (' position_units{p} ')'])
    % ylabel([v_label ' x velocity (' position_units{p} '/' time_units{t} ')'])
end

%% PLOT velocity scatter plot
% *************************************************************************
% figures{end+1} = figure;
% subplot(2,1,1)
% plot(x,vx, '.', 'MarkerSize', 8), hold all
% plot(x,vy, '.', 'MarkerSize', 8), hold all
% legend('vx vs. x', 'vy vs. x')
% xlabel(['x position (' position_units{p} ')'])
% ylabel(['velocity (' position_units{p} '/' time_units{t} ')'])
% title([folder_name ': ' file_name{1}], 'interpreter', 'none')
% set(gca, 'FontSize', 16)
% 
% subplot(2,1,2)
% plot(y,vy, '.', 'MarkerSize', 8), hold all
% plot(y,vx, '.', 'MarkerSize', 8), hold all
% legend('vy vs. y', 'vx vs. y')
% xlabel(['y position (' position_units{p} ')'])
% ylabel(['velocity (' position_units{p} '/' time_units{t} ')'])
% set(gca, 'FontSize', 16)

figures{end+1} = figure;
plot(x,y, '.', 'MarkerSize', 8), hold all
xlabel('x position (px)')
ylabel('y position (px)')
set(gca, 'FontSize', 16)


%% PLOT velocity color maps
% *************************************************************************
% close all

if min(size(vx_map)) > 1
    figures{end+1} = figure;
    subplot(2,1,1)
%     subplot(3,1,1)
    
    contour_levels = linspace(-35, 15, 100);
%     contour_levels = linspace(min(min(vx_map)), max(max(vx_map)), 100);
%     % contour_levels = linspace(min(min(vx_map)), 10, 100);
%     % contour_levels = linspace(min(min(min(vx_map)),min(min(vy_map))), max(max(max(vx_map)),max(max(vy_map))), 100);

%     contourf(x_edges(1:end-1), y_edges(1:end-1),vx_map',...
%                 'LineStyle', 'none',...
%                 'LevelListMode', 'manual', ...
%                 'LevelList', contour_levels)

    p_vx = pcolor(vx_map');
    p_vx.EdgeColor = 'none';
    % FIX X AND Y AXIS
    
    colormap(flipud(jet))
    colorbar
    caxis([min(contour_levels), max(contour_levels)])
    
    xlabel(['x position (' position_units{p} ')'])
    ylabel(['y position (' position_units{p} ')'])
    title([file_name{1} ' // ' v_label ' vx (' position_units{p} '/' time_units{t} ')'], 'interpreter', 'none')

%     subplot(3,1,2)
%     p_vx = pcolor(vx_map_centre');
%     p_vx.EdgeColor = 'none';
%     
%     colormap(flipud(jet))
%     colorbar
%     caxis([min(contour_levels), max(contour_levels)])
%     
%     xlabel(['x position (' position_units{p} ')'])
%     ylabel(['y position (' position_units{p} ')'])
%     title([file_name{1} ' // centre vx (' position_units{p} '/' time_units{t} ')'], 'interpreter', 'none')

    subplot(2,1,2)
%     subplot(3,1,3)
    % contour_levels = linspace(-10, 10, 100);
    % contour_levels = linspace(min(min(vy_map)), max(max(vy_map)), 100);
%     contourf(x_edges(1:end-1), y_edges(1:end-1),vy_map',...
%                 'LineStyle', 'none',...
%                 'LevelListMode', 'manual', ...
%                 'LevelList', contour_levels)
%     colormap(flipud(jet))
    
    p_vy = pcolor(vy_map');
    p_vy.EdgeColor = 'none';
    
    colormap(flipud(jet))
    colorbar
    caxis([min(contour_levels), max(contour_levels)])
    
    xlabel(['x position (' position_units{p} ')'])
    ylabel(['y position (' position_units{p} ')'])
    title([file_name{1} ' // ' v_label ' vy (' position_units{p} '/' time_units{t} ')'], 'interpreter', 'none')
end




%% PLOT mass and size color maps
% *************************************************************************
% close all
if min(size(vx_map)) > 30
    figures{end+1} = figure;
    subplot(2,1,1)
    % subplot(3,1,1)
    % contour_levels = linspace(-35, 5, 100);
    contour_levels = linspace(min(min(masses_map)), max(max(masses_map)), 100);
    % contour_levels = linspace(min(min(mass_map)), 10, 100);

%     contourf(x_edges(1:end-1), y_edges(1:end-1),masses_map',...
%                 'LineStyle', 'none',...
%                 'LevelListMode', 'manual', ...
%                 'LevelList', contour_levels)
    
    p_mass = pcolor(masses_map');
    p_mass.EdgeColor = 'none';
    
    colormap(jet)
    colorbar
    caxis([min(contour_levels), max(contour_levels)])
    
    xlabel(['x position (' position_units{p} ')'])
    ylabel(['y position (' position_units{p} ')'])
    title([file_name{1} ' //  mass'], 'interpreter', 'none')


    subplot(2,1,2)
    % subplot(3,1,1)
    % contour_levels = linspace(-35, 5, 100);
    contour_levels = linspace(min(min(sizes_map)), max(max(sizes_map)), 100);
    % contour_levels = linspace(min(min(size_map)), 10, 100);
    
%     contourf(x_edges(1:end-1), y_edges(1:end-1),sizes_map',...
%                 'LineStyle', 'none',...
%                 'LevelListMode', 'manual', ...
%                 'LevelList', contour_levels)
            
    p_size = pcolor(sizes_map');
    p_size.EdgeColor = 'none';
    
    colormap(jet)
    colorbar
    caxis([min(contour_levels), max(contour_levels)])
    
    xlabel(['x position (' position_units{p} ')'])
    ylabel(['y position (' position_units{p} ')'])
    title([file_name{1} ' //  size (px)'], 'interpreter', 'none')
end

%% PLOT velocity quiver map
% *************************************************************************
if min(size(vx_map)) > 3
    figures{end+1} = figure;
    % subplot(3,1,3)
    quiver(x_edges(1:end-1),y_edges(1:end-1),vx_map',vy_map','LineWidth',2)
    xlabel(['x position (' position_units{p} ')'])
    ylabel(['y position (' position_units{p} ')'])
    title([file_name{1} ' // ' v_label ' v (' position_units{p} '/' time_units{t} ')'], 'interpreter', 'none')
    set(gca,'FontSize',16)
    xlim([-5,35])
    ylim([-2,16])
end

%% PLOT velocity line plots
% *************************************************************************
% close all
if min(size(vx_map)) > 3
    figures{end+1} = figure('Units','normalized','Position',[0.01 0.07 0.95 0.8]);
    figure_velocity_line = figures{end};
    suptitle([folder_name ': ' file_name{1}])
    subplot(2,2,1)
    legend_cell{1} = {};
    for i = 1:1:size(vx_map,2)
        h{1}(i) = plot(x_edges(1:end-1),vx_map(:,i), '.-', 'MarkerSize', 8); hold all
        legend_cell{1}{end+1} = ['y = ' num2str(y_edges(i),'%02.0f') ' ' position_units{p}];
    end
    xlabel(['x position (' position_units{p} ')'])
    ylabel([v_label ' x velocity (' position_units{p} '/' time_units{t} ')'])

    subplot(2,2,2)
    legend_cell{2} = {};
    for i = 1:1:size(vx_map,1)
        h{2}(i) = plot(y_edges(1:end-1),vx_map(i,:), '.-', 'MarkerSize', 8); hold all
        legend_cell{2}{end+1} = ['x = ' num2str(x_edges(i),'%02.0f') ' ' position_units{p}];
    end
    xlabel(['y position (' position_units{p} ')'])
    ylabel([v_label ' x velocity (' position_units{p} '/' time_units{t} ')'])


    subplot(2,2,3)
    legend_cell{3} = {};
    for i = 1:1:size(vy_map,2)
        h{3}(i) = plot(x_edges(1:end-1),vy_map(:,i), '.-', 'MarkerSize', 8); hold all
        legend_cell{3}{end+1} = ['y = ' num2str(y_edges(i),'%02.0f') ' ' position_units{p}];
    end
    xlabel(['x position (' position_units{p} ')'])
    ylabel([v_label ' y velocity (' position_units{p} '/' time_units{t} ')'])

    subplot(2,2,4)
    legend_cell{4} = {};
    for i = 1:1:size(vy_map,1)
        h{4}(i) = plot(y_edges(1:end-1),vy_map(i,:), '.-', 'MarkerSize', 8); hold all
        legend_cell{4}{end+1} = ['x = ' num2str(x_edges(i),'%02.0f') ' ' position_units{p}];
    end
    xlabel(['y position (' position_units{p} ')'])
    ylabel([v_label ' y velocity (' position_units{p} '/' time_units{t} ')'])

    % format velocity line plots
    % *************************************************************************
    figure(figure_velocity_line)

    % [selected_colour, ~] = listdlg('PromptString', 'Colour scheme:',...
    %                            'SelectionMode', 'single', ...
    %                            'ListString', colour_type,...
    %                            'InitialValue', selected_colour);
    for i = 1:1:4
        subplot(2,2,i)
        legend(legend_cell{i}, 'Location', 'EO')
%         title([folder_name ': ' file_name{1}], 'interpreter', 'none')
        set(gca, 'FontSize', 12)
        for j = 1:1:numel(h{i})
            colour_RGB = colour_gradient(j, ...
                            numel(h{i}), ...
                            colour_type(selected_colour));
            h{i}(j).Color = colour_RGB;
        end
    end
%     suptitle([folder_name ': ' file_name{1}])
end

%% PLOT size and mass line plots
% *************************************************************************
% close all

if min(size(vx_map)) > 30
    figures{end+1} = figure('Units','normalized','Position',[0.01 0.07 0.95 0.8]);
    figure_velocity_line = figures{end};
    suptitle([folder_name ': ' file_name{1}])
    subplot(2,2,1)
    legend_cell{1} = {};
    for i = 1:1:size(masses_map,2)
        h{1}(i) = plot(x_edges(1:end-1),masses_map(:,i), '.-', 'MarkerSize', 8); hold all
        legend_cell{1}{end+1} = ['y = ' num2str(y_edges(i),'%02.0f') ' ' position_units{p}];
    end
    xlabel(['x position (' position_units{p} ')'])
    ylabel('mass')

    subplot(2,2,2)
    legend_cell{2} = {};
    for i = 1:1:size(masses_map,1)
        h{2}(i) = plot(y_edges(1:end-1),masses_map(i,:), '.-', 'MarkerSize', 8); hold all
        legend_cell{2}{end+1} = ['x = ' num2str(x_edges(i),'%02.0f') ' ' position_units{p}];
    end
    xlabel(['y position (' position_units{p} ')'])
    ylabel('mass')


    subplot(2,2,3)
    legend_cell{3} = {};
    for i = 1:1:size(sizes_map,2)
        h{3}(i) = plot(x_edges(1:end-1),sizes_map(:,i), '.-', 'MarkerSize', 8); hold all
        legend_cell{3}{end+1} = ['y = ' num2str(y_edges(i),'%02.0f') ' ' position_units{p}];
    end
    xlabel(['x position (' position_units{p} ')'])
    ylabel('size (px)')

    subplot(2,2,4)
    legend_cell{4} = {};
    for i = 1:1:size(sizes_map,1)
        h{4}(i) = plot(y_edges(1:end-1),sizes_map(i,:), '.-', 'MarkerSize', 8); hold all
        legend_cell{4}{end+1} = ['x = ' num2str(x_edges(i),'%02.0f') ' ' position_units{p}];
    end
    xlabel(['y position (' position_units{p} ')'])
    ylabel('size (px)')

    % format velocity line plots
    % *************************************************************************
    figure(figure_velocity_line)

    % [selected_colour, ~] = listdlg('PromptString', 'Colour scheme:',...
    %                            'SelectionMode', 'single', ...
    %                            'ListString', colour_type,...
    %                            'InitialValue', selected_colour);
    for i = 1:1:4
        subplot(2,2,i)
        legend(legend_cell{i}, 'Location', 'EO')
%         title([folder_name ': ' file_name{1}], 'interpreter', 'none')
        set(gca, 'FontSize', 12)
        for j = 1:1:numel(h{i})
            colour_RGB = colour_gradient(j, ...
                            numel(h{i}), ...
                            colour_type(selected_colour));
            h{i}(j).Color = colour_RGB;
        end
    end
end

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

