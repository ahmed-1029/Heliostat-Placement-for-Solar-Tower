clear
latitude = 19.7633; 
longitude = 96.0785; 
time_zone_offset = 6.5; 
standard_meridian = time_zone_offset * 15; 

hours_per_year = 24 * 365; 
solar_time = zeros(1, hours_per_year); 


[day_of_year, hour_of_day] = ndgrid(1:365, 0:23);
day_of_year = day_of_year(:); 
hour_of_day = hour_of_day(:); 

B = 2 * pi * (day_of_year - 81) / 364; 
EoT = 9.87 * sin(2 * B) - 7.53 * cos(B) - 1.5 * sin(B); 

longitude_correction = 4 * (longitude - standard_meridian); 

for i = 1:hours_per_year
    LST = hour_of_day(i) + time_zone_offset; 
    
    EoT_hours = EoT(day_of_year(i)) / 60; 
    long_corr_hours = longitude_correction / 60; 

    solar_time(i) = LST + EoT_hours + long_corr_hours;
end


save('solar_time.mat', 'solar_time');
h_s = (solar_time - 12)*.15;
N = 1:365;
S = 23.45*sind((360/365)*(284 + N));
S_expanded = repelem(S, 24);
solar_altitude_angle = zeros(1,length(S_expanded));
for i = 1:length(h_s)
    solar_alt_angle = asind(cosd(h_s(i)) .* cosd(S_expanded(i)) .* cosd(latitude) ...
        + sind(S_expanded(i)) .* sind(latitude));
    solar_altitude_angle(i) = solar_alt_angle;
end

solar_zenith = zeros(1,length(S_expanded));
for i = 1:length(S_expanded)
    solar_zenith(i) = 90 - solar_altitude_angle(i);
end

solar_azimuth_factor = zeros(1, 8760);
for i = 1:8760
solar_azimuth_factor(i) = asind(cosd(S_expanded(i))*sind(h_s(i))/sind(solar_zenith(i)));
end

solar_azimuth_angle = zeros(1,8760);
for i = 1:8760
    if cosd(h_s(i)) >= tand(S_expanded(i))/tand(latitude)
        solar_azimuth_angle(i) = 180 - solar_azimuth_factor(i);
    else
        solar_azimuth_angle(i) = 180 + solar_azimuth_factor(i);
    end    

end 
surface_azimuth_angle = zeros(1,8760);
for i = 1:8760
    if solar_azimuth_angle(i)> solar_azimuth_factor(i)
        surface_azimuth_angle(i) = solar_azimuth_factor(i) + 90;
    else
        surface_azimuth_angle(i) = solar_azimuth_factor(i) - 90;
    end 
end


DM = sqrt(11^2 + 11^2) + 0.2; 
Rmin = DM * cosd(30);
R1 = 62.5267770487835;
R2 = 125.053554097567;
R3 = 250.107108195134;
rows1 = 4;
rows2 = 9;
rows3 = 6;


distance = 0; row = 0;
distanceS = zeros(1,4);
for zones = 1:3
    while row < 4 
        distance = R1 + Rmin*row;
        distanceS(row+1) = distance;
        row = row + 1;
    end
    z = 0;
    while row < 13
        distance = R2 + Rmin*z;
        z = z + 1;
        distanceS(row+1) = distance;
        row = row + 1;
    end
    z = 0;
    while row < 19
        distance = R3 + Rmin*z;
        z = z + 1;
        distanceS(row+1) = distance;
        row = row + 1;
    end   
end

number_of_helio = 1150;
height_h = 7;
height_t = 95;
solar_altitude_angle_tower = zeros(1,19);

for i = 1:19
    solar_altitude_angle_tower(i) = atand((height_t - height_h)/distanceS(i));
end    

distance_to_reciever = zeros(1,19);
for i = 1:length(distanceS)
    d = sqrt(distanceS(i)^2 + (height_t - height_h)^2);
    distance_to_reciever(i) = d;
end

repetitions = [25,25,25,25,50,50,50,50,50,50,50,50,50,100,100,100,100,100,100]; 

dist_2_reciever_final = repelem(distance_to_reciever, repetitions);

solar_alt_angle_tower_final = repelem(solar_altitude_angle_tower, repetitions);


d = 10;
l = 8;
L = 20;

acceptance_angle_r = asind((-2*d*l + (L)*sqrt(4*l^2 + L^2 - d^2))/(4*l^2 + L^2));
lamda = 90 - acceptance_angle_r;

cos_eff = zeros(1,8760);
for i=1:8760
    cos_eff(i) = (sqrt(2)/2)*(sind(solar_altitude_angle(i))*cosd(lamda) - cosd(surface_azimuth_angle(i) ...
        - solar_azimuth_angle(i))*cosd(solar_altitude_angle(i))*sind(lamda) + 1)^0.5;
end
cosine_mean = mean(cos_eff(:));


G_sc = 1367; 
G_o = zeros(1, 8760); 

for i = 1:8760
    N = ceil(i / 24); 
    G_o(i) = G_sc * (1 + 0.033 * cosd((360 * N) / 365)) * cosd(solar_zenith(i));
end

m_opt = zeros(1,8760);
for i = 1:8760
    m_opt(i) = 1/(sind(solar_altitude_angle(i)));
end

d_rm = zeros(1,8760);
for i = 1:8760
    d_rm(i) = (6.6296+1.7513*m_opt(i) - 0.12*m_opt(i)^2 + ...
        0.0065*m_opt(i)^3 - 0.00013*m_opt(i)^4)^-1;
end


Beta = zeros(1150,8760);
for i=1:1150
    for b = 1:8760
        Beta(i,b) = (solar_alt_angle_tower_final(i) + solar_altitude_angle(b))/2;
    end
end    

h_s2 = repmat(h_s, 1150, 1);
surface_azimuth_angle2 = repmat(surface_azimuth_angle, 1150, 1);
S_expanded2 = repmat(S_expanded, 1150, 1);

theta = zeros(1150,8760);
for i = 1:1150
    for b = 1:8760
    theta(i,b) = acosd( ...
        (sind(latitude) .* sind(S_expanded2(i,b)) .* cosd(Beta(i,b))) - ...
        (cosd(latitude) .* sind(S_expanded2(i,b)) .* sind(Beta(i,b)) .* cosd(surface_azimuth_angle2(i,b))) + ...
        (cosd(latitude) .* cosd(S_expanded2(i,b)) .* cosd(h_s2(i,b)) .* cosd(Beta(i,b))) + ...
        (sind(latitude) .* cosd(S_expanded2(i,b)) .* cosd(h_s2(i,b)) .* sind(Beta(i,b)) .* cosd(surface_azimuth_angle2(i,b))) + ...
        (cosd(S_expanded2(i,b)) .* sind(h_s2(i,b)) .* sind(Beta(i,b)) .* sind(surface_azimuth_angle2(i,b))))  ;
        
    end

end

Tlk = 2;
B_ic = zeros(1,8760);
for i = 1:8760
    B_ic(i) = ((G_o(i))*exp(-0.8662*Tlk*m_opt(i)*d_rm(i)))*sind(30);
end

n_opt = zeros(1150,8760);
n_ref = 0.95;
n_s = 0.944;

n_at = zeros(1150,1);
for i = 1:1150
    n_at(i) = 0.99321 - 0.0001176 * dist_2_reciever_final(i) + (1.97e-8) * dist_2_reciever_final(i)^2;
end
   
n_at2 = repmat(n_at,1,8760);
cos_eff2 = repmat(cos_eff,1150,1);
for i = 1:1150
    for b = 1:8760
        n_opt(i,b) = n_at2(i,b)*n_ref*n_s*cos_eff2(i,b);

    end    
end

B_ic2 = repmat(B_ic, 1150, 1);
A = 121;
Q = zeros(1150,8760);
for i = 1:1150
    for b = 1:8760
        Q(i,b) = n_opt(i,b)*A*B_ic2(i,b);

    end    

end
n_opt_mean = mean(n_opt(:));


