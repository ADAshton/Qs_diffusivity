% Gulf of Mexico barriers diffusivity

%% format and save WIS data
inputFolder_WIS = ['/Users/rosepalermo/Documents/GitHub/Qs_diffusivity/GulfofMexico/WIS_GOM/test/WIS_save/']; % where is your WW3 data?
% % formatWIS(inputFolder_WIS,inputFolder_WIS)
ShorelineAngle_init = 0;
saveWISdata(inputFolder_WIS,inputFolder_WIS,ShorelineAngle_init)
%% compute wave climate
inputFolder_Waves = ['/Users/rosepalermo/Documents/GitHub/Qs_diffusivity/GulfofMexico/WIS_GOM/test/WIS_save/'];
outputFolder_waveclimate = ['/Users/rosepalermo/Documents/GitHub/Qs_diffusivity/GulfofMexico/WIS_GOM/test/waveclimate/'];
plot_on = 0; % 1 to make plots, 0 to not
computewaveclimate_CERC(inputFolder_Waves,outputFolder_waveclimate,plot_on)

%% look up diffusivity for each barrier given sl angle and closest station

load('coord_az.mat')
load('wis_info.mat')
diffusivity = zeros(length(sl_az),1);
for i = 1:length(diffusivity)
    diffusivity(i) = lookup_diffusivity(sl_az(i),wis_station(i),outputFolder_waveclimate);
end
