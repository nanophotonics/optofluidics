function TimelapseInfo = h5_TimelapseInfo(FilePath)
% read info from h5 file
% FilePath = 'C:\Users\Ana Andres\Documents\NanoPhotonics\Laboratory\2016.04.27 - 60nm Au aggregates\2016-04-27-60nmAu.h5';
FileInfo = h5info(FilePath);

% menu_group = 1;
menu_group = menu('Choose group:', FileInfo.Groups.Name);
GroupInfo = h5info(FilePath, FileInfo.Groups(menu_group).Name);

% read info from selected group
TimelapseOptions = cell(size(GroupInfo.Groups));
for i = 1:1:size(GroupInfo.Groups,1)
    TimelapseOptions{i} = [strrep(GroupInfo.Groups(i).Name, [FileInfo.Groups(menu_group).Name '/'], '')...
        ' // ' ...
        h5readatt(FilePath,...
            [GroupInfo.Groups(i).Name '/' ...
            GroupInfo.Groups(i).Datasets(1).Name], 'information')];      
end
% menu_timelapse = 1;
% menu_timelapse = menu('Choose timelapse to plot:', groupInfo.Groups.Name);
menu_timelapse = menu('Choose timelapse:', TimelapseOptions);
% timelapseInfo = h5info(FilePath, GroupInfo.Groups(menu_timelapse).Name);
TimelapsePath = GroupInfo.Groups(menu_timelapse).Name;
TimelapseInfo = h5info(FilePath, TimelapsePath);

% % reading the data
% number_of_spectra = size(GroupInfo.Groups(menu_timelapse).Datasets,1);
% number_of_wavelengths = GroupInfo.Groups(menu_timelapse).Datasets(1).Dataspace(1).Size;
% data = zeros(number_of_spectra, number_of_wavelengths);
% for i = 1:1:number_of_spectra
%     data(i,:) = h5read(FilePath, ...
%         [timelapseInfo.Name '/' timelapseInfo.Datasets(i).Name]);
% end
end