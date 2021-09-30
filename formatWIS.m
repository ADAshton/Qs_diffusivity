function formatWIS(inputFolder_WISformat,outputFolder_WISData)

% code that formats WIS data


% % Get a list of all files in the folder with the desired file name pattern.
infilePattern = fullfile(inputFolder_WISformat, '*.nc'); % Change to whatever pattern you need.
input_files = dir(infilePattern);

for runs = 1 : length(input_files)
    baseFileName = input_files(runs).name;
    fullFileName = fullfile(inputFolder_WISformat, baseFileName);
    [pathstr,filename,ext] = fileparts(fullFileName);
    savename = [outputFolder_WISData filename '.mat'];
    
    if runs == 1
        info_1 = ncinfo(fullFileName);
        data = zeros(info_1.Dimensions.Length,length(input_files));
        station_ID = str2double(info_1.Attributes(51).Value);
    end
    
    data(:,1) = ncread(fullFileName,'waveHs');
    data(:,2) = ncread(fullFileName,'waveMeanDirection');
    data(:,3) = ncread(fullFileName,'waveTm');
end
    
    wis_station_info = readmatrix('wis_station_list.xlsx');
    allstations = wis_station_info(:,1);
    idx_station = find(allstations == station_ID);
    depth = wis_station_info(idx_station,5);
    save(savename,'data','depth')
end