% WW3 barrier diffusivity

load('barrier_diffusivity_ww3_inputs_v2.mat')

k = 0.17;
dsf = 10;
ShAng = 5:5:180;
Ang = zeros(length(wave_angle0),1);
diffusivity_WW3 = zeros(length(wave_angle0),1);
sl_az(sl_az>180) = sl_az(sl_az>180)-180;


for j = 1:length(wave_angle0)
    if wave_angle0(j)>180 % to use the wave crest orientation -/+ 90 from the orientation reported by WIS...
        Ang(j,:) = (mod(wave_angle0(j)+90-sl_az(j)+180,360)-180);
    else
        Ang(j,:) = (mod(wave_angle0(j)-90-sl_az(j)+180,360)-180);
    end
end
asign = sign(Ang);
Angle = abs(Ang);
diffusivity_WW3 = -k./dsf.*T.^(1/5).*hs.^(12/5).*((cos(deg2rad(Angle).^(1/5))).*(((6/5)*sin(deg2rad(Angle)).^2)-cos(deg2rad(Angle)).^2));
diff_km_yr = diffusivity_WW3*365*24*3600/1e6;