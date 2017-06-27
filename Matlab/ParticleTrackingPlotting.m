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

% invert the x axis so the particle travells from left to right
% raw_data{i}.x = max(raw_data{i}.x) - raw_data{i}.x;


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

%% PLOT

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
plot(raw_data{1}.size, raw_data{1}.mass, '.', 'MarkerSize', 8), hold all
plot(particles.size, particles.mass, '.', 'MarkerSize', 8), hold all
% plot(1:1:20, 3e3*(((1:1:20)-6).^2+0), '-k', 'LineWidth', 2), hold all
xlabel('particle size (px)')
ylabel('mass')
% title([folder_name ': ' file_name{1}], 'interpreter', 'none')
set(gca, 'FontSize', 12)
legend([num2str(size(raw_data{1},1)) ' points'], [num2str(size(particles,1)) ' particles'], 'Location', 'NW')

subplot(1,2,2)
plot(particles.frames, particles.x_travelled, '.')
ylabel('x travelled (px)')
xlabel('frames')
% title([folder_name ': ' file_name{1}], 'interpreter', 'none')
set(gca, 'FontSize', 12)
legend([num2str(size(particles,1)) ' particles'])

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
% plot(raw_data{i}.y, raw_data{i}.mass, '.', 'MarkerSize', 8), hold all
% xlabel('y particle position (px)')
% ylabel('mass')
% title([folder_name ': ' file_name{1}], 'interpreter', 'none')
% set(gca, 'FontSize', 16)
% legend([num2str(size(raw_data{1},1)) ' points'])

% figures{end+1} = figure;
% % subplot(2,2,4)
% hframes = histogram(particles.frames, 100);
% hframes.Normalization = 'count';
% xlabel('frames per particle')
% ylabel(hframes.Normalization)
% title([folder_name ': ' file_name{1}], 'interpreter', 'none')
% set(gca, 'FontSize', 16)



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
    frame_rate = 159.22; % frames / second --> 6.28 ms temporal resolution
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
        end
%     end
    % the number of particles is equal to size(raw_data{i},1) - size(x,2)
end
% velocity = table(x,y,vx,vy);
% velocity.Properties.VariableNames = {'x','y','vx','vy'};

%% PLOT velocity histogram
% *************************************************************************
% close all

figures{end+1} = figure;
% subplot(1,2,1)
legend_hv = {};
text_hv = {};

number_of_bins = 100;
% number_of_bins = round(size(particles,1)/10);

hvx = histogram(vx,number_of_bins); hold all, legend_hv{end+1} = 'vx';
hvy = histogram(vy,number_of_bins); hold all, legend_hv{end+1} = 'vy';

hvx.Normalization = 'count';
% hvx.Normalization = 'probability';

hvy.Normalization = hvx.Normalization;

text_hv{end+1} = 'Gaussian fit:';
text_hv{end+1} = 'y =  a1*exp(-((x-b1)/c1)^2)';
text_hv{end+1} = '';

hvx_x = hvx.BinEdges(1:end-1) + hvx.BinWidth/2;
hvx_y = hvx.Values;
fit_vx = fit(hvx_x', hvx_y', 'gauss1');
confint_vx = confint(fit_vx);
c1_error = abs((confint_vx(1,3)-confint_vx(2,3))/2);
pvx = plot(fit_vx); hold all, legend_hv{end+1} = 'vx fit';
pvx.Color = 'b';
pvx.LineWidth = 2;
text_hv{end+1} = ['mean(vx) = ' num2str(mean(vx)) ' ' position_units{p} '/' time_units{t}];
% text_hv{end+1} = ['mean(abs(vx)) = ' num2str(mean(abs(vx))) ' ' position_units{p} '/' time_units{t}];
text_hv{end+1} = ['vx: a1 = ' num2str(fit_vx.a1, '%.3g') ' \pm ' ...
    num2str(abs(confint_vx(1,1)-confint_vx(2,1))/2, '%.2g') ...
    ' '];
text_hv{end+1} = ['vx: b1 = ' num2str(fit_vx.b1, '%.1f') ' \pm ' ...
    num2str(abs(confint_vx(1,2)-confint_vx(2,2))/2, '%.1f') ' ' ...
    position_units{p} '/' time_units{t}];
text_hv{end+1} = ['vx: c1 = ' num2str(fit_vx.c1, '%.1f') ' \pm ' ...
    num2str(abs(confint_vx(1,3)-confint_vx(2,3))/2, '%.1f')  ' ' ...
    position_units{p} '/' time_units{t}];
text_hv{end+1} = ['vx: D = ' ...
    num2str(1/frame_rate/4*(fit_vx.c1)^2, '%.3g') ' \pm ' ...
    num2str(1/frame_rate*(fit_vx.c1)/2*c1_error, '%.2g') ' \mum^2/s'];
text_hv{end+1} = '';

hvy_x = hvy.BinEdges(1:end-1) + hvy.BinWidth/2;
hvy_y = hvy.Values;
fit_vy = fit(hvy_x', hvy_y', 'gauss1');
confint_vy = confint(fit_vy);
c1_error = abs((confint_vy(1,3)-confint_vy(2,3))/2);
pvy = plot(fit_vy); hold all, legend_hv{end+1} = 'vy fit';
pvy.Color = 'r';
pvy.LineWidth = 2;
text_hv{end+1} = ['mean(vy) = ' num2str(mean(vy)) ' ' position_units{p} '/' time_units{t}];
% text_hv{end+1} = ['mean(abs(vy)) = ' num2str(mean(abs(vy))) ' ' position_units{p} '/' time_units{t}];
text_hv{end+1} = ['vy: a1 = ' num2str(fit_vy.a1, '%.3g') ' \pm ' ...
    num2str(abs(confint_vy(1,1)-confint_vy(2,1))/2, '%.2g') ...
    ' '];
text_hv{end+1} = ['vy: b1 = ' num2str(fit_vy.b1, '%.1f') ' \pm ' ...
    num2str(abs(confint_vy(1,2)-confint_vy(2,2))/2, '%.1f') ...
    position_units{p} '/' time_units{t}];
text_hv{end+1} = ['vy: c1 = ' num2str(fit_vy.c1, '%.1f') ' \pm ' ...
    num2str(abs(confint_vy(1,3)-confint_vy(2,3))/2, '%.1f') ...
    position_units{p} '/' time_units{t}];
text_hv{end+1} = ['vy: D = ' ...
    num2str(1/frame_rate/4*(fit_vy.c1)^2, '%.3g') ' \pm ' ...
    num2str(1/frame_rate*(fit_vy.c1)/2*c1_error, '%.2g') ' \mum^2/s'];

legend(legend_hv)
xlabel(['velocity (' position_units{p} '/' time_units{t} ')'])
ylabel(hvx.Normalization)
title([folder_name ': ' file_name{1}], 'interpreter', 'none')
set(gca, 'FontSize', 16)


text('Units','normalized','Position',[0.08,0.95], ...
    'FontSize', 12, 'VerticalAlignment', 'top', 'String' , text_hv)

% subplot(1,2,2)
% semilogy(hvx_x,hvx_y, '.-'), hold all
% semilogy(hvy_x,hvy_y, '.-'), hold all
% xlabel(['velocity (' position_units{p} '/' time_units{t} ')'])
% ylabel(hvx.Normalization)
% set(gca, 'FontSize', 16)
% legend('vx','vy')
% grid on
% suptitle([folder_name ': ' file_name{1}])

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

% figures{end+1} = figure;
% hvx_x = hvx.BinEdges(1:end-1) + hvx.BinWidth/2;
% hvx_y = hvx.Values;
% vx_colours = jet(numel(hvx.Values));
% for i = numel(hvx.Values):-1:1
%     indices_larger = find(vx >= hvx.BinEdges(i));
%     indices_smaller = find(vx < hvx.BinEdges(i+1));
%     indices = intersect(indices_larger, indices_smaller);
%     plot(x(indices), y(indices), ...
%         '.', 'MarkerSize', 8, ...
%         'Color', vx_colours(i,:))
%     hold all
% end
% 
% xlabel('x position (px)')
% ylabel('y position (px)')
% set(gca, 'FontSize', 16)
% title([folder_name ': ' file_name{1}], 'interpreter', 'none')


%% CALCULATE velocity maps
% *************************************************************************
clear map_labels

menu_vmap = 5;
menu_vmap = menu('Velocity:', 'Average', 'Maximum', 'Minimum', 'Mode', 'Gaussian Fit');

x_sections = 19;
y_sections = 11;

% x_sections = 1;
% y_sections = 3;

x_edges = linspace(min(x), max(x), x_sections + 1);
x_centres = (x_edges(1:end-1)+x_edges(2:end))/2;
x_bin = discretize(x,x_edges);
y_edges = linspace(min(y), max(y), y_sections + 1);
y_centres = (y_edges(1:end-1)+y_edges(2:end))/2;
y_bin = discretize(y,y_edges);

masses_map = zeros(numel(x_centres),numel(y_centres));
sizes_map = zeros(numel(x_centres),numel(y_centres));
vx_grid = cell(numel(x_centres),numel(y_centres));

vx_map = zeros(numel(x_centres),numel(y_centres),1);
vy_map = zeros(numel(x_centres),numel(y_centres),1);

for i = 1:1:numel(x_centres)
    x_indices = find(x_bin == i);
    for j = 1:1:numel(y_centres)
        y_indices = find(y_bin == j);
        xy_indices = intersect(x_indices, y_indices);
%         disp(numel(xy_indices))

        if numel(xy_indices) < 30
            vx_map(i,j,:) = NaN;
            vy_map(i,j,:) = NaN;
            masses_map(i,j) = NaN;
            sizes_map(i,j) = NaN;
        else
            if menu_vmap == 1 % average
                vx_grid{i,j} = vx(xy_indices);
                vx_map(i,j,1) = mean(vx(xy_indices));
                vy_map(i,j,1) = mean(vy(xy_indices));
                masses_map(i,j) = mean(masses(xy_indices));
                sizes_map(i,j) = mean(sizes(xy_indices));
                map_labels{1} = ['average v (' position_units{p} '/' time_units{t} ')'];

            elseif menu_vmap == 2 % maximum
                if isempty(vx(xy_indices))
                    vx_map(i,j,1) = NaN;
                else           
                    vx_map(i,j,1) = max(vx(xy_indices));
                end
                if isempty(vy(xy_indices))
                    vy_map(i,j,1) = NaN;
                else           
                    vy_map(i,j,1) = max(vy(xy_indices));
                end
                map_labels{1} = ['max v (' position_units{p} '/' time_units{t} ')'];

            elseif menu_vmap == 3 % mimimum
                if isempty(vx(xy_indices))
                    vx_map(i,j,1) = NaN;
                else           
                    vx_map(i,j,1) = min(vx(xy_indices));
                end
                if isempty(vy(xy_indices))
                    vy_map(i,j,1) = NaN;
                else           
                    vy_map(i,j,1) = min(vy(xy_indices));
                end
                map_labels{1} = ['min v (' position_units{p} '/' time_units{t} ')'];

            elseif menu_vmap == 4 % mode
                vx_map(i,j,1) = mode(round(vx(xy_indices),0));
                vy_map(i,j,1) = mode(round(vy(xy_indices),0));
                map_labels{1} = ['mode v (' position_units{p} '/' time_units{t} ')'];
            
            elseif menu_vmap == 5 % gaussian fit
                clc
                disp(['Calculating map section: (' num2str(i) ',' num2str(j) ') of (' ...
                    num2str(numel(x_centres)) ',' num2str(numel(y_centres)) ')'])
                
                number_of_bins = min(round(numel(xy_indices)/10), 100);
                
                [vx_counts,vx_edges] = histcounts(vx(xy_indices),number_of_bins);
                vx_centres = (vx_edges(1:end-1) + vx_edges(2:end)) / 2;
                
                fit_parameters = fit(vx_centres', vx_counts', 'gauss1');
                vx_map(i,j,1) = fit_parameters.b1; % gaussian centre
                vx_map(i,j,2) = fit_parameters.c1; % gaussian width
                vx_map(i,j,3) = fit_parameters.a1; % gaussian height
                vx_map(i,j,4) = 1/frame_rate/4*fit_parameters.c1^2; % D (units of velocity^2/time)
                
                [vy_counts,vy_edges] = histcounts(vy(xy_indices),number_of_bins);
                vy_centres = (vy_edges(1:end-1) + vy_edges(2:end)) / 2;
                
                fit_parameters = fit(vy_centres', vy_counts', 'gauss1');
                vy_map(i,j,1) = fit_parameters.b1; % gaussian centre
                vy_map(i,j,2) = fit_parameters.c1; % gaussian width
                vy_map(i,j,3) = fit_parameters.a1; % gaussian height               
                vy_map(i,j,4) = 1/frame_rate/4*fit_parameters.c1^2; % D (units of velocity^2/time)

                vx_map(i,j,5) = numel(xy_indices); % number of data points
                vy_map(i,j,5) = numel(xy_indices); % number of data points
                
                vx_map(i,j,6) = number_of_bins; % number of bins
                vy_map(i,j,6) = number_of_bins; % number of bins
                
%                 vx_map(i,j,5) = numel(xy_indices); % number of data points
%                 vy_map(i,j,5) = number_of_bins; % number of bins
                
                map_labels{1} = ['Gauss. centre v (' position_units{p} '/' time_units{t} ')'];
                map_labels{2} = ['Gauss. width v (' position_units{p} '/' time_units{t} ')'];
                map_labels{3} = 'Gauss. height (counts)';
                map_labels{4} = ['D (' position_units{p} '^2/' time_units{t} ')'];
                map_labels{5} = 'No. of points';
                map_labels{6} = 'No. of bins';
                
%                 plot(vx_centres,vx_counts), hold all
%                 plot(vy_centres,vy_counts), hold all
            end
        end
    end
end
disp('Finished the map calculations.')



%% PLOT velocity color maps
% *************************************************************************
% close all

axis_labels = {'x','y'};
menu_vmap_axis = [1,2];
menu_vmap_plot = 1;

% [menu_vmap_plot, ~] = listdlg('PromptString', 'Velocity colour maps to plot:',...
%                            'SelectionMode', 'multiple', ...
%                            'ListString', map_labels,...
%                            'InitialValue', menu_vmap_plot);

% [menu_vmap_axis, ~] = listdlg('PromptString', 'Axis to plot:',...
%                            'SelectionMode', 'multiple', ...
%                            'ListString', axis_labels,...
%                            'InitialValue', menu_vmap_axis);
                       
[menu_vmap_plot, menu_vmap_axis] = dialog_two_lists('Select plot options:', ...
    'Variable:', map_labels, menu_vmap_plot,...
    'Axis:', axis_labels, menu_vmap_axis);
                       
map_variables = {vx_map,vy_map};

if min(size(vx_map(:,:,1))) > 5
    for i = menu_vmap_plot
        figures{end+1} = figure;
        k = 0;
        for j = menu_vmap_axis            
            k = k + 1;
            subplot(numel(menu_vmap_axis),1,k)

            % set NaN to blank
            imAlpha=ones(size(map_variables{j}(:,:,i)'));
            imAlpha(isnan(map_variables{j}(:,:,i)'))=0;

%             contour_levels = linspace(-36.2922, 22.1115, 100);
%             contour_levels = linspace(min(min(map_variables{j}(:,:,i))), max(max(map_variables{j}(:,:,i))), 100);
%             contour_levels = linspace(min(min(min(map_variables{1}(:,:,i))),min(min(map_variables{2}(:,:,i)))), ...
%                 max(max(max(map_variables{1}(:,:,i))),max(max(map_variables{2}(:,:,i)))), 100);
            contour_levels = linspace(min(min(min(map_variables{1}(2:end-1,2:end-1,i))), ...
                min(min(map_variables{2}(2:end-1,2:end-1,i)))), ...
                max(max(max(map_variables{1}(2:end-1,2:end-1,i))), ...
                max(max(map_variables{2}(2:end-1,2:end-1,i)))), 100);

            p_vx = imagesc(x_centres, y_centres, map_variables{j}(:,:,i)','AlphaData',imAlpha); 
%             set(gca,'YDir','normal')
    %         p_vx = pcolor(x_edges(1:end-1), y_edges(1:end-1), map_variables{j}(:,:,i)');
    %         p_vx.EdgeColor = 'none';

            colormap(flipud(jet))
%             colormap(jet)
            c = colorbar;
    %         c.Label.String = [map_labels{i} ' vx (' position_units{p} '/' time_units{t} ')'];
            c.Label.String = [axis_labels{j} ' ' map_labels{i}];
            c.FontSize = 16;

            caxis([min(contour_levels), max(contour_levels)])

            xlabel(['x position (' position_units{p} ')'])
            ylabel(['y position (' position_units{p} ')'])
            title([file_name{1} ' // average = ' ...
                num2str(nanmean(nanmean(map_variables{j}(:,:,i))), '%.3g')])
            set(gca, 'FontSize', 16)
        end
        
    end
end

%% PLOT velocity histogram by sections
% *************************************************************************
% close all

if max(size(vx_map(:,:,1))) < 5
    figures{end+1} = figure;
    legend_vx_grid = {};
    text_hv = {};
    text_hv{end+1} = 'Gaussian fit:';
    text_hv{end+1} = 'y =  a1*exp(-((x-b1)/c1)^2)';
    text_hv{end+1} = '';
    
    grid_colours = parula(size(vx_grid,1)*size(vx_grid,2));
    k = 1;
    for i = 1:1:size(vx_grid,1)
        for j = 1:1:size(vx_grid,2)
            ph(i,j) = histogram(vx_grid{i,j},30); hold all
            ph(i,j).Normalization = 'count';
%             ph(i,j).Normalization = 'probability';
            ph(i,j).FaceColor = grid_colours(k,:);
            legend_vx_grid{end+1} = ['x = ' num2str(x_edges(i),'%02.0f') ' ' position_units{p} ...
                ', y = ' num2str(y_edges(j),'%02.0f') ' ' position_units{p}];
            
            hvx_x = ph(i,j).BinEdges(1:end-1) + ph(i,j).BinWidth/2;
            hvx_y = ph(i,j).Values;
            fit_vx_grid{i,j} = fit(hvx_x', hvx_y', 'gauss1');
            confint_vx_grid{i,j} = confint(fit_vx_grid{i,j});
            
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
    %     h{1}(i) = plot(x_edges(1:end-1),vx_map(:,i,1), '.-', 'MarkerSize', 8); hold all
    %     legend_cell{1}{end+1} = ['y = ' num2str(y_edges(i),'%02.0f') ' ' position_units{p}];
    % end
    % xlabel(['x position (' position_units{p} ')'])
    % ylabel([map_labels{1} ' x velocity (' position_units{p} '/' time_units{t} ')'])
end


%% PLOT mass and size color maps
% *************************************************************************
% close all
if min(size(vx_map(:,:,1))) > 10 && 1
    figures{end+1} = figure;
    subplot(2,1,1)
    % subplot(3,1,1)
    % set NaN to blank
    imAlpha=ones(size(masses_map'));
    imAlpha(isnan(masses_map'))=0;
    % contour_levels = linspace(-35, 5, 100);
    contour_levels = linspace(min(min(masses_map)), max(max(masses_map)), 100);
    % contour_levels = linspace(min(min(mass_map)), 10, 100);

%     contourf(x_edges(1:end-1), y_edges(1:end-1),masses_map',...
%                 'LineStyle', 'none',...
%                 'LevelListMode', 'manual', ...
%                 'LevelList', contour_levels)
    
%     p_mass = pcolor(x_edges(1:end-1), y_edges(1:end-1), masses_map');
%     p_mass.EdgeColor = 'none';
    p_mass = imagesc(x_centres, y_centres, masses_map','AlphaData',imAlpha); 
    
    colormap(jet)
    c = colorbar;
    c.Label.String = 'mass';
    c.FontSize = 16;
    caxis([min(contour_levels), max(contour_levels)]);
    
    set(gca,'FontSize',16)
    xlabel(['x position (' position_units{p} ')'])
    ylabel(['y position (' position_units{p} ')'])
%     title([file_name{1} ' //  mass'], 'interpreter', 'none')
    title(file_name{1} , 'interpreter', 'none')
    

    subplot(2,1,2)
    % subplot(3,1,1)
    % contour_levels = linspace(-35, 5, 100);
    % set NaN to blank
    imAlpha=ones(size(sizes_map'));
    imAlpha(isnan(sizes_map'))=0;
    contour_levels = linspace(min(min(sizes_map)), max(max(sizes_map)), 100);
    % contour_levels = linspace(min(min(size_map)), 10, 100);
    
%     contourf(x_edges(1:end-1), y_edges(1:end-1),sizes_map',...
%                 'LineStyle', 'none',...
%                 'LevelListMode', 'manual', ...
%                 'LevelList', contour_levels)
            
%     p_size = pcolor(x_edges(1:end-1), y_edges(1:end-1), sizes_map');
%     p_size.EdgeColor = 'none';
    p_size = imagesc(x_centres, y_centres, sizes_map','AlphaData',imAlpha); 
    
    colormap(jet)
    c = colorbar;
    c.Label.String = 'size (px)';
    c.FontSize = 16;
    caxis([min(contour_levels), max(contour_levels)])
    caxis([min(contour_levels), max(contour_levels)])
    
    set(gca,'FontSize',16)
    xlabel(['x position (' position_units{p} ')'])
    ylabel(['y position (' position_units{p} ')'])
%     title([file_name{1} ' //  size (px)'], 'interpreter', 'none')
end

%% PLOT velocity quiver map
% *************************************************************************
if min(size(vx_map(:,:,1))) > 3 && 0
    figures{end+1} = figure;
    % subplot(3,1,3)
    quiver(x_edges(1:end-1),y_edges(1:end-1),vx_map(:,:,1)',vy_map(:,:,1)','LineWidth',2)
    xlabel(['x position (' position_units{p} ')'])
    ylabel(['y position (' position_units{p} ')'])
    title([file_name{1} ' // ' map_labels{1} ' v (' position_units{p} '/' time_units{t} ')'], 'interpreter', 'none')
    set(gca,'FontSize',16)
    xlim([-5,35])
    ylim([-2,16])
end

%% PLOT velocity line plots
% *************************************************************************
% close all
if min(size(vx_map(:,:,1))) > 3 && 0
    figures{end+1} = figure('Units','normalized','Position',[0.01 0.07 0.95 0.8]);
    figure_velocity_line = figures{end};
    suptitle([folder_name ': ' file_name{1}])
    subplot(2,2,1)
    legend_cell{1} = {};
    for i = 1:1:size(vx_map,2)
        h{1}(i) = plot(x_edges(1:end-1),vx_map(:,i,1), '.-', 'MarkerSize', 8); hold all
        legend_cell{1}{end+1} = ['y = ' num2str(y_edges(i),'%02.0f') ' ' position_units{p}];
    end
    xlabel(['x position (' position_units{p} ')'])
    ylabel([map_labels{1} ' x velocity (' position_units{p} '/' time_units{t} ')'])

    subplot(2,2,2)
    legend_cell{2} = {};
    for i = 1:1:size(vx_map,1)
        h{2}(i) = plot(y_edges(1:end-1),vx_map(i,:,1), '.-', 'MarkerSize', 8); hold all
        legend_cell{2}{end+1} = ['x = ' num2str(x_edges(i),'%02.0f') ' ' position_units{p}];
    end
    xlabel(['y position (' position_units{p} ')'])
    ylabel([map_labels{1} ' x velocity (' position_units{p} '/' time_units{t} ')'])


    subplot(2,2,3)
    legend_cell{3} = {};
    for i = 1:1:size(vy_map,2)
        h{3}(i) = plot(x_edges(1:end-1),vy_map(:,i,1), '.-', 'MarkerSize', 8); hold all
        legend_cell{3}{end+1} = ['y = ' num2str(y_edges(i),'%02.0f') ' ' position_units{p}];
    end
    xlabel(['x position (' position_units{p} ')'])
    ylabel([map_labels{1} ' y velocity (' position_units{p} '/' time_units{t} ')'])

    subplot(2,2,4)
    legend_cell{4} = {};
    for i = 1:1:size(vy_map,1)
        h{4}(i) = plot(y_edges(1:end-1),vy_map(i,:,1), '.-', 'MarkerSize', 8); hold all
        legend_cell{4}{end+1} = ['x = ' num2str(x_edges(i),'%02.0f') ' ' position_units{p}];
    end
    xlabel(['y position (' position_units{p} ')'])
    ylabel([map_labels{1} ' y velocity (' position_units{p} '/' time_units{t} ')'])

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

if min(size(vx_map(:,:,1))) > 30
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

