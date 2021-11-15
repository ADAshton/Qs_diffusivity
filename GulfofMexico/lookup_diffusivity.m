function [diffusivity] = lookup_diffusivity(sl_azimuth,stationID,waveclimatefolder)
if sl_azimuth>180
    sl_azimuth = sl_azimuth-180;
end
infilePattern_wc = fullfile(waveclimatefolder, '*.mat'); % Change to whatever pattern you need.
input_files_wc = dir(infilePattern_wc);

for runs = 1 : length(input_files_wc)
    baseFileName_wc = input_files_wc(runs).name;
    fullFileName_wc = fullfile(waveclimatefolder, baseFileName_wc);
    if contains(fullFileName_wc,num2str(stationID))
        load(fullFileName_wc)
    else
        continue;
    end
% find intersection between gammas and shoreline azimuth
[~,diffusivity] = polyxpoly(ShAngs,Gammas(:,1),[sl_azimuth sl_azimuth],[-1 1]);
% plot(ShAngs,Gammas(:,1),'g','linewidth',1.5)
% hold on
% xline(sl_azimuth)
end