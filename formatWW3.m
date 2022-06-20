function formatWW3(inputFolder_ww3format,outputFolder_WW3Data)

% code that formats WW3 data


% % Get a list of all files in the folder with the desired file name pattern.
infilePattern = fullfile(inputFolder_ww3format, '*.nc'); % Change to whatever pattern you need.
input_files = dir(infilePattern);

for runs = 1 : length(input_files)
    baseFileName = input_files(runs).name;
    fullFileName = fullfile(inputFolder_ww3format, baseFileName);
    [pathstr,filename,ext] = fileparts(fullFileName);
    savename = strcat(outputFolder_WW3Data,filename,'.mat');
    
    if runs == 1
        info_1 = ncinfo(fullFileName);
        data = zeros(info_1.Dimensions.Length,length(input_files));
    end
    data(:,1) = ncread(fullFileName,'hs')./1000; % mm to m
    data(:,2) = ncread(fullFileName,'dp')./100; %centi-degrees to degrees
    data(:,3) = ncread(fullFileName,'tp')./1000; % ms to s
    
end

depth = abs(info_1.Attributes(5).Value);


save(savename,'data','depth')
end