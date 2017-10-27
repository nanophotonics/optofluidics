clc
clear
close all

%% CHOOSING FILES
% *************************************************************************

% specify default path
folder_path = 'R:\aa938\NanoPhotonics\Output\Proposals\2017 Part III\PBG modes\';


% pop up window to choose the file(s) to read from a SINGLE FOLDER
[file_name, folder_path, ~] = uigetfile('.tif',...
                                      'Images to Read (use CTRL to select multiple files)',...
                                      folder_path,...
                                      'MultiSelect','on');
file_name = cellstr(file_name); % convert to cell array of strings
number_of_files = size(file_name,2);

%% READING DATA
% *************************************************************************


for ifn = 1:1:number_of_files
    file_path = [folder_path file_name{ifn}];
    disp(['Image ' num2str(ifn) '/' num2str(number_of_files)])

    [raw_data,map] = imread(file_path);

    disp(' ')
    im = double(raw_data(:,:,1));
    
    im = im-min(min(im));
    im = im/max(max(im));
    
    fig = figure('Units','normalized','Position',[0 0 1 1]);
    level_list = 0:0.01:1;
    step_1 = 10;
    X = 1:step_1:size(im,1);
    Y = 1:step_1:size(im,2);

    contourf(X, Y, im(X,Y)', ...
        'edgecolor','none', 'levellist', level_list)
    colormap(jet)
%     c = colorbar('location', 'southoutside');
%     caxis([0,1])
%     c.Label.String = 'normalised intensity';
    axis equal
    xlabel('x axis (pixels)')
    ylabel('y axis (pixels)')    
    title('Data')
    
    saveas(fig, strrep(file_path, '.tif', '_jet.png'), 'png'); % saving the .png file
end

close all