function[PAR_HR2000, NIR_HR2000]=calculate_PAR_NIR_HR2000(wvl,irradiance_data_HR2000)                   
index=find(wvl>=400 & wvl<=700);
par=single(nan(length(index),1));
irra=single(nan(length(index),1));
coeff=single(nan(length(index),1));
h=6.626*10^(-34);
c=3*10^(8);
for j=1:length(index)
    irra(j)=irradiance_data_HR2000(index(j))*((wvl(index(j)+1)-wvl(index(j)-1)))./2*10/1000;
    coeff(j)=1/(h*c./wvl(index(j))*10^9*6.02*10^23)*10^6;
    par(j)=irra(j)*coeff(j);
end

index1=find(wvl>=730 & wvl<=780);
nir=single(nan(length(index1),1));
for j=1:length(index1)
    nir(j)=irradiance_data_HR2000(index1(j))*((wvl(index1(j)+1)-wvl(index1(j)-1)))./2*10/1000;
end

PAR_HR2000=nansum(par);
NIR_HR2000=nansum(nir);
end