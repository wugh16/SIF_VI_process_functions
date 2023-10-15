function[NDVI, NIRv, EVI, CIrededge, CIgreen, PRI]=calculate_VI(wvl,irradiance_data_HR2000,radiance_data_HR2000)                   
spec=wvl;
irra_full=irradiance_data_HR2000;
rad_full=radiance_data_HR2000;

% NDVI
ndvi_temp_irra=radwv(irra_full,spec,[770,780;650,660]);
ndvi_temp_rad=radwv(rad_full,spec,[770,780;650,660]);
ndvi_temp=ndvi_temp_rad.*pi./ndvi_temp_irra;
ndvi_temp( ndvi_temp<0 |  ndvi_temp>1)=nan;
NDVI = (ndvi_temp(1)-ndvi_temp(2))/(ndvi_temp(1)+ndvi_temp(2));
       
%NIRv
nir=ndvi_temp(1);
NIRv=ndvi.*ndvi_temp(1);

% EVI
evi_temp_irra=radwv(irra_full,spec,[770,780;650,660;460,470]);
evi_temp_rad=radwv(rad_full,spec,[770,780;650,660;460,470]);
evi_temp=evi_temp_rad.*pi./evi_temp_irra;
evi_temp(evi_temp<0 | evi_temp>1)=nan;
EVI = 2.5*(evi_temp(1)-evi_temp(2))/(1+evi_temp(1)+6*evi_temp(2)-7.5*evi_temp(3));

% CIrededge
CIrededge_temp_irra=radwv(irra_full,spec,[770,780;720,730]);
CIrededge_temp_rad=radwv(rad_full,spec,[770,780;720,730]);                     
CIrededge_temp = CIrededge_temp_rad./(CIrededge_temp_irra./pi);
CIrededge_temp(CIrededge_temp<0 | CIrededge_temp>1)=nan;
CIrededge = CIrededge_temp(1)/CIrededge_temp(2)-1;
        
% CIgreen
CIgreen_temp_irra=radwv(irra_full,spec,[770,780;545,565]);
CIgreen_temp_rad=radwv(rad_full,spec,[770,780;545,565]);
CIgreen_temp =  CIgreen_temp_rad./(CIgreen_temp_irra./pi);
CIgreen_temp(CIgreen_temp<0 | CIgreen_temp>1)=nan;
CIgreen=CIgreen_temp(1)/CIgreen_temp(2)-1;  
  
%PRI
pri570_temp_irra=radwv(irra_full,spec,[530.5,531.5;569.5,570.5]);
pri570_temp_rad=radwv(rad_full,spec,[530.5,531.5;569.5,570.5]);
pri570_temp=pri570_temp_rad.*pi./pri570_temp_irra;
pri570_temp(pri570_temp<0 | pri570_temp>1)=nan;
PRI = (pri570_temp(1)-pri570_temp(2))/(pri570_temp(1)+pri570_temp(2));
end