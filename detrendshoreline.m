function [SLdata_detrend,RotationAngle] = detrendshoreline(SLdata,ind_open)
[r_notnan,~] = find(~isnan(SLdata)); % get rid of nans
[p,S,mu] = polyfit(SLdata(ind_open,1),SLdata(ind_open,2),1); % calculate data slope
y_polyval = polyval(p,SLdata(:,1),[],mu); % calculate y with slope
y_detrend = SLdata(:,2) - y_polyval;

% detrend and translate to origin
SLdata_detrend = [SLdata(:,1)-min(SLdata(:,1)) y_detrend-min(y_detrend)];

RotationAngle = rad2deg(atan(p(1)));
end