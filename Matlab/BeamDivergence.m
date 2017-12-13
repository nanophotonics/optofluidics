clc
clear
close all
figures = {};

%% CHOOSING FILES
% *************************************************************************

% specify default path
folder_path = 'R:\aa938\NanoPhotonics\Laboratory\';
% folder_path = 'R:\3-Temporary\aa938\';
% folder_path = 'R:\3-Temporary\sjd87\';
% folder_path = 'R:\3-Temporary\pc594\30817\';


% pop up window to choose the file(s) to read from a SINGLE FOLDER
[file_name, folder_path, ~] = uigetfile('.png',...
                                      'Images to Read (use CTRL to select multiple files)',...
                                      folder_path,...
                                      'MultiSelect','on');

file_name = cellstr(file_name); % convert to cell array of strings
file_path = cell(size(file_name));
for i = 1:1:size(file_name,2)
    file_path{i} = [folder_path file_name{i}];
end
number_of_files = size(file_name,2);


% folder_path_save = folder_path;
folder_path_save = [folder_path 'matlab fit\'];

if strfind(folder_path, 'Laboratory')
    slash = strfind(folder_path, '\');
    slash_index = find(slash > strfind(folder_path, 'Laboratory')+11);
%     folder_name = folder_path(strfind(folder_path, 'Laboratory')+11:end-1);
    folder_name = folder_path(strfind(folder_path, 'Laboratory')+11:slash(slash_index(1))-1);
elseif strfind(folder_path, 'aa938')
    slash = strfind(folder_path, '\');
    slash_index = find(slash > strfind(folder_path, 'aa938')+6);
%     folder_name = folder_path(strfind(folder_path, 'aa938')+6:end-1);
    folder_name = folder_path(strfind(folder_path, 'aa938')+6:slash(slash_index(1))-1);
elseif strfind(folder_path, 'sjd87')
    slash = strfind(folder_path, '\');
    slash_index = find(slash > strfind(folder_path, 'sjd87')+6);
%     folder_name = folder_path(strfind(folder_path, 'aa938')+6:end-1);
    folder_name = folder_path(strfind(folder_path, 'sjd87')+6:slash(slash_index(1))-1);
else
    folder_name = folder_path;
end    


%% READING DATA
% *************************************************************************

w1_mm = zeros(size(file_name));
w2_mm = zeros(size(file_name));
d_cm = zeros(size(file_name));
lambda_nm = zeros(size(file_name));
p_mW = zeros(size(file_name));
sample = cell(size(file_name));

for ifn = 1:1:number_of_files
    disp(['Image ' num2str(ifn) '/' num2str(number_of_files)])

    raw_data = imread(file_path{ifn}, 'png');

    disp(' ')
    im = double(raw_data);

    % z = reshape(im,[],1);
    % figure
    % h = histogram(z); hold all
    % background = h.BinEdges(h.Values==max(h.Values))*0.8;
    % background = 0;
    % im = im-background;
    im = im(:,:,1);
    im = im - mean(mean(im(1:100,1:100)))*0.9;
    im = im/max(max(im));

    %% 2D Gaussian fit

    fittype_Gauss_xy = fittype(@(A,a,b,c,x0,y0,z0,x,y) ...
            A*exp(-(a*(x-x0).^2 + 2*b*(x-x0).*(y-y0) + c*(y-y0).^2)) + z0, ...
            'independent', {'x', 'y'}, ...
            'dependent', 'z' );

    [X,Y] = ndgrid(1:size(im,1),1:1:size(im,2));

    x = reshape(X,[],1);
    y = reshape(Y,[],1);
    z = reshape(im,[],1);

    % A,a,b,c,x0,y0,z0
    start_point = [1,1e-6,1e-7,1e-6,450,758,0];
    lower_limit = [0,0,0,0,0,0,0];
    upper_limit = [2,1e-3,1e-3,1e-3,size(im,1),size(im,2),0];

    disp('Fitting 2D Gaussian...')
    fitobject_Gauss_xy = fit([x,y],z,fittype_Gauss_xy,...
        'StartPoint', start_point,...
        'Lower', lower_limit, ...
        'Upper', upper_limit);    
    disp('Done!')

    disp(fitobject_Gauss_xy)
    
    A = fitobject_Gauss_xy.A;

    a = fitobject_Gauss_xy.a;
    b = fitobject_Gauss_xy.b;
    c = fitobject_Gauss_xy.c;

    x0 = fitobject_Gauss_xy.x0;
    y0 = fitobject_Gauss_xy.y0;
    z0 = fitobject_Gauss_xy.z0;

    z_fit = A*exp(-(a*(x-x0).^2 + 2*b*(x-x0).*(y-y0) + c*(y-y0).^2)) + z0;
    im_fit = reshape(z_fit,size(im));
    
    %% Solving equations
    
    syms sx sy th
    eqns = [a == cos(th)^2/2/sx^2+sin(th)^2/2/sy^2, ...
        b == -sin(2*th)/4/sx^2+sin(2*th)/4/sy^2, ...
        c == sin(th)^2/2/sx^2+cos(th)^2/2/sy^2];
    solution = vpasolve(eqns, [sx sy th]);

    wx = 2*abs(double(solution.sx));
    wy = 2*abs(double(solution.sy));
    th = rem(double(solution.th),2*pi);

    pixel = 5.3e-6;
    
%     w1_mm(ifn) = max(wx,wy)*pixel*1e3;
%     w2_mm(ifn) = min(wx,wy)*pixel*1e3;
    
    w1_mm(ifn) = wx*pixel*1e3;
    w2_mm(ifn) = wy*pixel*1e3;
    
    d_cm(ifn) = str2double(file_name{ifn}(1:end-6));
    dash = strfind(file_name{ifn}, '-');
%     sample{ifn} = file_name{ifn}(dash(1)+1:dash(2)-1);
    sample{ifn} = file_name{ifn}(1:dash(1)-1);
    lambda_nm(ifn) = str2double(file_name{ifn}(dash(1)+1:dash(2)-3));
    p_mW(ifn) = str2double(file_name{ifn}(dash(2)+1:end-6));
    disp(' ')
    disp(['w1 = ' num2str(w1_mm(ifn)) ' mm'])
    disp(['w2 = ' num2str(w2_mm(ifn)) ' mm'])
    disp(['th = ' num2str(th/pi) 'pi rad'])
    disp(' ')

    %% Ploting data and fits
    plot_data = 1;    
    if plot_data       

        figure_handle = figure('Units','normalized','Position',[0.01 0.07 0.95 0.55]);
        level_list = 0:0.01:1;
        step_1 = 10;
        X = 1:step_1:size(im,1);
        Y = 1:step_1:size(im,2);

        subplot(1,4,1)
        contourf(X, Y, im(X,Y)', ...
            'edgecolor','none', 'levellist', level_list)
        c = colorbar('location', 'southoutside');
        caxis([0,1])
        axis equal
        xlabel('x axis (pixels)')
        ylabel('y axis (pixels)')
        c.Label.String = 'normalised intensity';
        title('Data')

        subplot(1,4,2)
        contourf(X, Y, im_fit(X,Y)', ...
            'edgecolor', 'none', 'levellist', level_list)
        c = colorbar('location', 'southoutside');
        caxis([0,1])
        axis equal
        xlabel('x axis (pixels)')
        ylabel('y axis (pixels)')
        c.Label.String = 'normalised intensity';
        title('Gaussian 2D fit')

        step_1 = 1;
        step_2 = 200;
        subplot(1,4,3)
        for i = 1:step_2:size(im,2)
            plot(1:step_1:size(im,1),im(1:step_1:end,round(i))); hold all
            plot(1:step_1:size(im_fit,1),im_fit(1:step_1:end,round(i)), ...
                'k-', 'LineWidth', 1); hold all
        end
        xlabel('x axis (pixels)')
        ylabel('normalised intensity')
        title('Line profiles')
        xlim([0,size(im,1)])
        ylim([0,1])

        subplot(1,4,4)
        for i = 1:step_2:size(im,1)
            plot(1:step_1:size(im,2),im(round(i),1:step_1:end)); hold all
            plot(1:step_1:size(im_fit,2),im_fit(round(i),1:step_1:end), 'k-', ...
                'LineWidth', 1); hold all
        end
        xlabel('y axis (pixels)')
        ylabel('normalised intensity')
        title('Line profiles')
        xlim([0,size(im,2)])
        ylim([0,1])

        suptitle([folder_name ' - ' file_name{ifn}(1:end-4) ...
            ': w1 = ' num2str(w1_mm(ifn), '%.2f') ' mm' ...
            ', w2 = ' num2str(w2_mm(ifn), '%.2f') ' mm' ])
        
        hgexport(figure_handle, [folder_path_save file_name{ifn}], hgexport('factorystyle'), 'Format', 'png')
        close all
    end
end

%% Plotting beam divergence

% folder_path = 'R:\3-Temporary\pc594\30817\';
% [file_name, folder_path, ~] = uigetfile('.xlsx',...
%     'Beam divergence data to read)',...
%     folder_path,...
%     'MultiSelect','off');
% data = readtable([folder_path file_name]);
% d_cm = data.d_cm_;
% w1_mm = data.w1_mm_;
% w2_mm = data.w2_mm_;
% number_of_files = 2;

if number_of_files > 1
                                  
    figure
    legend_w = {};
    
%     offset = 136; % cm
    offset = 155; % cm
%     d_cm = d_cm + offset;
    
%     plot(d_cm, w1_mm, '.', 'MarkerSize', 20, 'LineWidth', 1), hold all
%     plot(d_cm, w2_mm, '.', 'MarkerSize', 20, 'LineWidth', 1), hold all
    for s = unique(sample)
        indices_s = find(~cellfun(@isempty,strfind(sample,s)));
        for l = unique(lambda_nm)
%             indices_l = find(~cellfun(@isempty,strfind(lambda_nm,l)));
            indices_l = find(lambda_nm == l);
            
            indices = intersect(indices_s,indices_l);
            plot(p_mW(indices), w1_mm(indices), '.-', 'MarkerSize', 20, 'LineWidth', 1), hold all
            legend_w{end+1} = ['wx, \lambda = ' num2str(l(1)) 'nm, ' s{1}];
            plot(p_mW(indices), w2_mm(indices), 'o-', 'MarkerSize', 6, 'LineWidth', 1), hold all
            legend_w{end+1} = ['wy, \lambda = ' num2str(l(1)) 'nm, ' s{1}];
        end
    end
    xlabel('P (mW)')
%     xlabel('d (cm)')
    ylabel('w (mm)')
%     xlim([0,250])
%     ylim([0,6])
    grid on

%     lambda = 800e-9;

%     fittype_Rayleigh = fittype('w0*sqrt(1 + (x*lambda/pi/w0^2)^2)', ...
%         'problem', 'lambda');
% 
%     start_point = 0.0004;
% 
%     fitobject_w1 = fit(d_cm'*1e-2, w1_mm'*1e-3, ...
%         fittype_Rayleigh, 'problem', lambda, 'startpoint', start_point);
%     disp(fitobject_w1)
%     w01 = fitobject_w1.w0;
% 
%     fitobject_w2 = fit(d_cm'*1e-2, w2_mm'*1e-3, ...
%         fittype_Rayleigh, 'problem', lambda, 'startpoint', start_point);
%     disp(fitobject_w2)
%     w02 = fitobject_w2.w0;
% 
%     x = 0:0.005:2.5;
%     y1 = w01*sqrt(1 + (x*lambda/pi/w01^2).^2);
%     y2 = w02*sqrt(1 + (x*lambda/pi/w02^2).^2);
% 
%     plot(x*1e2,y1*1e3, 'LineWidth', 1)
%     plot(x*1e2,y2*1e3, 'LineWidth', 1)
% 
%     title([folder_name ': w01 = ' num2str(w01*1e3,'%.3f') ' mm, '...
%         'w02 = ' num2str(w02*1e3,'%.3f') ' mm'])

%     legend('w1 data', 'w2 data', 'w1 fit', 'w2 fit', 'Location', 'NW')
    legend(legend_w, 'Location', 'best')
end

%% Save data
data_table = [cell2table(sample','VariableNames',{'sample'}), ...
              array2table(lambda_nm','VariableNames',{'lambda_nm'}), ...
              array2table(p_mW','VariableNames',{'p_mW'}), ...
              array2table(w1_mm','VariableNames',{'wx_mm'}), ...
              array2table(w2_mm','VariableNames',{'wy_mm'}), ...
              ];

disp(data_table)
save_data = 1;
if save_data
    [file_name_save, folder_path_save,~] = uiputfile(['.' 'csv'],...
               'File to Save the Data',folder_path_save);
    writetable(data_table, [folder_path_save file_name_save]);
end

%% end
        
disp('FINISHED EVERYTHING!')
disp('')