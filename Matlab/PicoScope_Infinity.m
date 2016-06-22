clc
clear
close all

% specify default path
file_path = 'C:\Users\Ana Andres\Documents\NanoPhotonics\Laboratory\2016.06.21 - 60 nm Au PicoScope\';

% pop up window to choose the file(s) to read
[file_name, file_path, ~] = uigetfile('.txt',...
    'PicoScope Fil to Read (use CTRL to select multiple files)',file_path,'MultiSelect','on');

% if just a single file is selected, change the variable file_names
% from a string to a cell so the loops indexing will work
if isa(file_name,'char')
    temporary = file_name;
    clear file_name
    file_name{1} = temporary;
    clear temporary;
end

number_of_files = size(file_name,2);
% 269547
% reading the data from the txt files
% it takes a while because the files are large
raw_data = cell(size(file_name));
for i = 1:1:number_of_files
    disp(['Reading File ' num2str(i) '/' num2str(number_of_files)])
%     raw_data{i} = dlmread([file_path file_name{i}], '\t',2,0);

    T = readtable([file_path file_name{i}],...
        'ReadRowNames', false,...
        'Delimiter', '\t');
%     T.Properties.VariableUnits = T{1,:};
%     T(1,:) = [];
%     C = table2array(T);
    
%     for k = 1:1:size(T,2)
%         disp(k)
%         for j = 1:1:size(T,1)
% %         for j = 269540:1:269555
%             if strcmp(T{j,k}, 'Infinity')
%                 T{j,k} = {'Inf'};
%             elseif strcmp(T{j,k}, '-Infinity')
%                 T{j,k} = {'-Inf'};
%             end
%         end
%     end
    
    % just remove the 3rd channel instead of modifiying the infinities
    T(:,4) = [];
    
    writetable(T, [file_path file_name{i}(1:end-4) '-AB.txt'], 'Delimiter', '\t')

end
disp('Finished reading all files!')