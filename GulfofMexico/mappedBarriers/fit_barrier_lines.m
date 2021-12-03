folder = '/Users/rosepalermo/Documents/GitHub/Qs_diffusivity/GulfofMexico/mappedBarriers/';
infilePattern = fullfile(folder, '*.kml'); % Change to whatever pattern you need.
input_files = dir(infilePattern);
sl_az = zeros(length(input_files),1);
Lat_all = cell(length(input_files),1);
Lon_all = cell(length(input_files),1);

for i = 1 : length(input_files)
    baseFileName = input_files(i).name;
    fullFileName = fullfile(folder, baseFileName);
    [Lat,Lon] = kml2latlong(baseFileName);
    pp = polyfit(Lon([1 length(Lon)]),Lat([1 length(Lon)]),1);
%     figure()
    hold on
    plot(Lon,Lat)
    plot(Lon,polyval(pp,Lon))
    Lat_all{i} = Lat;
    Lon_all{i} = Lon;
    sl_az(i) = rad2deg(atan(1/pp(1)));
end
sl_az(sl_az<0) = sl_az(sl_az<0)+180;
% save('barrier_coord_az','BarrierIsland','Lat_all','Lon_all','sl_az')