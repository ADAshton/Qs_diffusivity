folder = '/Users/rosepalermo/Documents/GitHub/Qs_diffusivity/GulfofMexico/mappedBarriers/';
infilePattern = fullfile(folder, '*.kml'); % Change to whatever pattern you need.
input_files = dir(infilePattern);
sl_az = zeros(length(input_files),1);
Lat_all = cell(length(input_files),1);
Lon_all = cell(length(input_files),1);
figure()
flip_to_east_west = [4 5 9];
for i = 1 : length(input_files)
    baseFileName = input_files(i).name;
    fullFileName = fullfile(folder, baseFileName);
    [Lat,Lon] = kml2latlong(baseFileName);
    if ismember(i,flip_to_east_west)
        Lat = flipud(Lat);
        Lon = flipud(Lon);
    end
    pp = polyfit(Lon,Lat,1);
    if i ==4 %Chandeleur is near vertical and doesn't fit well with the polyfit
        pp = polyfit(Lon([1 length(Lon)]),Lat([1 length(Lon)]),1);
    end
    subplot(5,4,i)
    hold on
    plot(Lon,Lat)
    scatter(Lon(1),Lat(1),'g')
    scatter(Lon(end),Lat(end),'r')
    plot(Lon,polyval(pp,Lon))
    axis equal
    Lat_all{i} = Lat;
    Lon_all{i} = Lon;
    sl_az(i) = rad2deg(atan(pp(1)));
    if sl_az(i)<0
        sl_az(i) = -sl_az(i)+270;
    else
        sl_az(i) = 90-sl_az(i);
    end
    if sl_az(i)>180
        sl_az(i) = sl_az(i)-180;
    end
    title(BarrierIsland_name{i})
%     title(sl_az(i))
end
% ind_neg = sl_az<0;
% ind_pos = sl_az>0;
% sl_az(ind_neg) = sl_az(ind_neg)-90;
% sl_az(ind_pos) = sl_az(ind_pos)+270;
% save('barrier_coord_az','BarrierIsland','Lat_all','Lon_all','sl_az')