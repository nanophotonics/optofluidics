clc
clear
close all

pathRead = 'R:\aa938\NanoPhotonics\Matlab\Optical Forces\2017.06.26\';
[nameRead, pathRead, ~] = uigetfile('.fig',...
    'Select figure:',pathRead,'MultiSelect','off');

h = open([pathRead nameRead]);

% h = gcf;
axes = get(h,'Children');
dataObjs = get(axes,'Children');
objTypes = get(dataObjs{2}, 'Type');
data = cell(size(objTypes,1),2);
i = 1;
for i = 1:1:size(data,1)
    data{i,1} = get(dataObjs{2}(i), 'XData');
    data{i,2} = get(dataObjs{2}(i), 'YData');
    data{i,3} = get(dataObjs{2}(i), 'ZData');
end

menu_save = 1;
menu_save = menu('Save data as .txt file?', 'NO', 'YES');

if menu_save == 2  
    pathSave = pathRead;
    nameSave = strrep(nameRead, 'fig', 'txt');
    
    pause(0.1)
    [nameSave,pathSave,~] = uiputfile(['.' 'txt'],...
        'Choose file to save the data in:',[pathSave nameSave]); % choosing the file name
    fid1 = fopen([pathSave nameSave], 'wt');
    
    fprintf(fid1, 'Title:\t%s\n', axes(2).Title.String{1});
    fprintf(fid1, 'Title:\t%s\n', axes(2).Title.String{2});
    
    fprintf(fid1, 'Axis Labels (x,y):\t');
    fprintf(fid1, '%s\t%s\t', axes(2).XLabel.String, axes(2).YLabel.String);
    fprintf(fid1, '\n');
    
    fprintf(fid1, 'Legend:\t');
    for i = 1:1:size(data,1)
        fprintf(fid1, '%s\t', axes(1).String{i});
    end
    fprintf(fid1, '\n');
    
    fprintf(fid1, 'Data (x1,y1, x2,y2, ...):\n');
    % assuming al data sets are of the same langth
    for j = 1:1:size(data{1,1},2)       
        for i = 1:1:size(data,1)
            fprintf(fid1, '%f\t%f\t', data{i,1}(j), data{i,2}(j));
        end
        fprintf(fid1, '\n');
    end
    fclose(fid1);
end