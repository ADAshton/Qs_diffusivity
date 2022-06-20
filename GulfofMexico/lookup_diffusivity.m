function [diffusivity,WaveAngle] = lookup_diffusivity(sl_azimuth,stationID,waveclimatefolder)
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
% figure()
plot(ShAngs,TrueDiffs,'g','linewidth',1.5)
xlabel('hypothetical shoreline angle')
ylabel('diffusivity')
hold on
xline(sl_azimuth)


[~,diffusivity] = polyxpoly(ShAngs,TrueDiffs,[sl_azimuth sl_azimuth],[-1e10 1e10]);

% load('binarray_WIS_rose.mat')
[~,indmax] = max(EnergyFluxIn);
WaveAngle = binarray(indmax);

end