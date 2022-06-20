% This code takes in a shoreline wave buoy data and outputs wave climate info, Qs, and diffusivity
% input is long and lat of shoreline, WaveWatch3, and WIS data
% only difference for using UTM is that "stretching" step is not needed.

% Last updated Sept. 8, 2021
% Rose Palermo

%% add main folder to path
% cd ..
% addpath(pwd)

%% Part 1: WAVES

%%%%% WIS and WW3 -- This only needs to be done once for each shoreline. Then saved data can be reused.
% here, we take WIS data and wave watch 3 data and save them as wave roses
% input is WIS and WW3 data as a table format saved in its own folder
% output is matlab files of computed wave climate info using CERC formula
% 
% WaveWatch3
inputFolder_WW3 = '/Users/rosepalermo/Documents/GitHub/Qs_diffusivity/barnegat_bay/WW3_out/'; % where is your WW3 data?
outputFolder_Waves = '/Users/rosepalermo/Documents/GitHub/Qs_diffusivity/barnegat_bay/Wave_Data/'; % where do you want the wave roses to go?
ShorelineAngle_init = 204; % where the shoreline angle loop starts
saveWW3data_barnegatbay(inputFolder_WW3,outputFolder_Waves,ShorelineAngle_init)

%%%%% compute wave climate -- this only needs to be done one time for a set of WIS & WW3 stations
outputFolder_waveclimate = "/Users/rosepalermo/Documents/GitHub/Qs_diffusivity/barnegat_bay/Waveclimate/";
plot_on = 1; % 1 to make plots, 0 to not
computewaveclimate_CERC(outputFolder_Waves,outputFolder_waveclimate,plot_on)

%% look up diffusivity for each barrier given sl angle and closest station

sl_az = 24.1669;
stationID = 842;
% Ang = (mod(sl_az,180));
% SLAng =  - atan(SLAng)/degtorad - RotationAngle;
[diffusivity,WaveAngle] = lookup_diffusivity_bb(sl_az,stationID,outputFolder_waveclimate)

% %% Part 2 shoreline + diffusivity and Qs calculations
% %%%%% load shoreline data
% % need long and lat points of the open ocean coast
% % load data as data_raw, an nx2 matrix of long and lat data
% load('CapeCodsl_example.mat')
% data_raw = SLdata;
% 
% %%%%% shorefixer
% % this function makes continuous coastline from fragmented segments by
% % joining expoints of segments less than TOL distance apart.
% tol = 0.01; % in degrees, change to m if UTM
% SLdata = join_cst(data_raw,tol);
% 
% 
% %%%%% open ocean only
% % here, create a variable of the indices of open ocean shoreline points (not in an estuary
% % or on the backbarrier). I haven't figured out a way to do this that isn't
% % by hand yet.
% ind_all = 1:length(SLdata);
% % loaded in CapeCodsl
% % ind_open = ind_all; % if they are the same;
% 
% %%%%% stretch
% % convert long and lat to x and y
% %the below conversion factors are calculated for 41.8 degrees latitude, a
% %central point along the cape, by a National Geospatial-Intelligence Agency
% %calculator:
% %http://msi.nga.mil/MSISiteContent/StaticFiles/Calculators/degree.html
% lat = 111.069; % cape cod example
% lon = 83.110; % cape cod example
% 
% SLdata(:,1) = (SLdata(:,1)-min(SLdata(:,1))) * lon;
% SLdata(:,2) = (SLdata(:,2)-min(SLdata(:,2))) * lat;
% 
% %%%%% detrend
% % rotate shoreline so that it is flat, land is down and ocean is up
% 
% % Confirm that the shoreline is mostly oriented up. If not, first rotate by 180
% % degrees using the rotate180 function.
% % [SLdata] = rotate180(SLdata); % uncomment if needed
% 
% [SLdata_detrend,RotationAngle] = detrendshoreline(SLdata,ind_open);
% 
% 
% %% Part 3: calculate sediment flux and diffusivity given shoreline and wave climate
% 
% % %Rose's example here
% % load('/Users/rosepalermo/Documents/Research/Luis Stuff/Files Rose Changed/caperot.mat')
% % SLdata_detrend = [Xrot Yrot];
% % ind_open = good;
% 
% plot_on = 1;
% outputFolder_waveclimate = 'FOLDER';
% outputFolder = 'FOLDER';
% 
% [Diffusivity_save] = calculateQs_Dif(outputFolder_waveclimate,outputFolder,SLdata_detrend,ind_open,RotationAngle,plot_on);
% 
