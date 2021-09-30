function [Xrot,Yrot] = rotate180(SLdata)

RotationAngle = 180; % degrees.

degtorad = pi/180;
RotAng = RotationAngle * degtorad;


Xrot = zeros(length(SLdata),1);
Yrot = zeros(length(SLdata),1);


for i = 1:length(SLdata)
    Xrot(i) = SLdata(i,1) * cos(RotAng) + SLdata(i,2) * sin(RotAng);
    Yrot(i) = -SLdata(i,1) * sin(RotAng) + SLdata(i,2) * cos(RotAng);
     
end

Xrot = Xrot -min(Xrot); % rotated X
Yrot = Yrot - min(Yrot); % rotated Y

end