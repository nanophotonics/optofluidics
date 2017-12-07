clear
% close all
clc

%% READING FILES
folder_path = 'R:\3-Temporary\aa938\';
[file_names, folder_path, ~] = uigetfile('.csv',...
                                      'Choose file',...
                                      folder_path,...
                                      'MultiSelect','on');
%% read data
data_table = table;
for i = 1:1:size(file_names,2)
    read_table = readtable([folder_path file_names{i}]);
    data_table = [data_table; read_table];
end

%% SORTING DATA
samples = unique(data_table.sample);
wavelengths = unique(data_table.lambda_nm);
legend_location = 'EO';
% legend_location = 'NW';
s = 1;
fit_data_table = table;
fit_column_names = {'wx_mm','wy_mm'};
figure
for i = 1:size(samples,1)
    indices_s = find(~cellfun(@isempty,strfind(data_table.sample,samples(i))));
    
%     colours = jet(size(indices,1));
%     scatter(data_table.lambda_nm(indices), data_table.p_mW(indices), data_table.wx_mm*15)    
    legend_cell = {};
    for j = 1:size(wavelengths,1)
        indices_l = find(data_table.lambda_nm == wavelengths(j));
        indices = intersect(indices_s, indices_l);
        legend_cell{end+1} = [num2str(wavelengths(j)) ' nm'];
%         legend_cell{end+1} = [num2str(wavelengths(j)) ' nm, wx'];
%         legend_cell{end+1} = [num2str(wavelengths(j)) ' nm, wy'];
        
        p = cell(2,1);
        for k = 1:size(fit_column_names,2)
            subplot(2,2,s+2*(k-1))
            plot(data_table.p_mW(indices), data_table{indices,fit_column_names(k)}, ...
                '.', 'MarkerSize', 18), hold all
            grid on
            xlabel('power (mW)')
            ylabel([fit_column_names{k}(1:2) ' (a.u.)'])
            title(samples(i))
            legend(legend_cell,'Location',legend_location)
            
            p{k} = polyfit(data_table.p_mW(indices),data_table{indices,fit_column_names(k)},1);
            plot(data_table.p_mW(indices),polyval(p{k},data_table.p_mW(indices)),...
                '-k', 'LineWidth', 1)
            
        end
        fit_table = table(samples(i), wavelengths(j), ...
            p{1}(1), p{1}(2), p{2}(1), p{2}(2));
        fit_data_table = [fit_data_table; fit_table];
        legend_cell{end+1} = 'linear fit';
        legend(legend_cell,'Location',legend_location)
    end
    s = s+1;
end
fit_data_table.Properties.VariableNames = {'sample','lambda_nm',...
                    'p1x','p2x','p1y','p2y'};
                
%% plot fit data
samples = unique(data_table.sample);
wavelengths = unique(data_table.lambda_nm);
legend_location = 'NW';

fit_column_names = {'p1x','p1y'}; % slope
% fit_column_names = {'p2x','p2y'}; % offset

figure
legend_cell = {};
for i = 1:size(samples,1)
    indices = find(~cellfun(@isempty,strfind(fit_data_table.sample,samples(i))));
    
    for k = 1:size(fit_column_names,2)
        legend_cell{end+1} = [fit_column_names{k} ', ' samples{i}];

        plot(fit_data_table.lambda_nm(indices), fit_data_table{indices,fit_column_names(k)}, ...
            '.-', 'MarkerSize', 18), hold all
        grid on
        xlabel('wavelength (nm)')
        ylabel('slope (mm/mW)')
%         ylabel('offset (mm)')
        legend(legend_cell,'Location',legend_location)
    end
end
