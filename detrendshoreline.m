function [SLdata_detrend,RotationAngle] = detrendshoreline(SLdata)

[p,S,mu] = polyfit(SLdata(:,1),SLdata(:,2)); % calculate data slope
y_poly = polyval(p,SLdata(:,1),[],mu); % calculate y with slope
y_detrend = SLdata(:,2) - y_polyval;

% detrend and translate to origin
SLdata_detrend = [SLdata(:,1)-min(SLdata(:,1)) y_detrend-min(y_detrend)];

RotationAngle = rad2deg(atan(p(1)));
end