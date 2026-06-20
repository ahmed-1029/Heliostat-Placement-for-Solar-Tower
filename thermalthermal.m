clear
helio_total = 1150; 
desp = 0.2;
DM = sqrt(11^2 + 11^2) + desp; 
Rmin = DM * cosd(30); 
Nhel1 = 25;
R1 = Nhel1*(DM/(2*pi));
az1 = 2*asin(DM/(2*R1));
Nheli = 0;
Nrowi = 0;
i = 1;
R = zeros(1,4);
while true
    azi = az1/(2^(i-1));
    Ri = (2^(i-1))*(DM/az1);
    R(i) = Ri;
    Nheli =  round(2*pi/(azi));
    m = i + 1;
    Rm = (2^(m-1))*(DM/az1);
    
    R(m) = Rm;
    Nrowi = floor((Rm-Ri)/Rmin);
    if helio_total - Nheli*Nrowi < 0
        break
    end    
    helio_total = helio_total - Nheli*Nrowi;
    i = i + 1;
end

Nrowi = helio_total/Nheli;
helio_total = helio_total - Nheli*Nrowi;
