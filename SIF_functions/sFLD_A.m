%Maria Pilar Cendrero Mateo
%sFLD
%%-------------------------------------------------------------------------

function[F,R]=sFLD_A(wvl,irradiance_data_QEPRO,radiance_data_QEPRO,do_plot,time,path_save_fig)
%WVL_IN - define inside absorption band
[~, wvl_758] = min(abs(wvl - 758));
[~, wvl_770] = min(abs(wvl - 770));
%wvl range absorption band
wvl_index_range_in = (wvl_758:wvl_770);
    
%WVL_OUT - two points outside absorption band
%out1
[~, wvl_745] = min(abs(wvl - 745));
wvl_index_out1 = (wvl_745:wvl_758);

%% Irradiance (E)
%Define  "E" PKS MAX outside absorption band
%wvl_index_out1
[value_E1,pks_E1] = findpeaks(irradiance_data_QEPRO(wvl_index_out1));
wvl_index_peaks_out1 =wvl_index_out1(pks_E1);

%Define "SAMPLE" MIN inside absorption band (only one value)
E_min = min(irradiance_data_QEPRO(wvl_index_range_in));
[~, wvl_E_index_in] = min(abs(irradiance_data_QEPRO - E_min));

%Value outside absorption band
%WR - retrieve fluorescence
E_out = irradiance_data_QEPRO(wvl_index_peaks_out1(end));

%sample - retrieve fluorescence
L_out = radiance_data_QEPRO(wvl_index_peaks_out1(end));

%Define 4 points inside abs band for wvl_E_index_in
pt1 = wvl_E_index_in;
if pt1>1
    pt2 = wvl_E_index_in-1;
    pt3 = wvl_E_index_in+1;
    pt4 = wvl_E_index_in+2;
    pt = [pt1 pt2 pt3 pt4];
else
    pt3 = wvl_E_index_in+1;
    pt4 = wvl_E_index_in+2;
    pt = [pt1 pt3 pt4];
end

%get wr_in and sample_in
E_in  = mean(irradiance_data_QEPRO(pt));
L_in = mean(radiance_data_QEPRO(pt));

%Compute Fluorescence and Reflectance
F=(E_out*L_in-L_out*E_in)./(E_out-E_in);
R=(L_out-L_in)./(E_out-E_in);

%% O2A - 760
if do_plot == 1 
    %Lin
    fig1 = figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(1,2,1)
    plot(wvl,irradiance_data_QEPRO,'-k')
    xlim ([740 780])
    hold on
    plot(wvl(wvl_index_peaks_out1(end)),irradiance_data_QEPRO(wvl_index_peaks_out1(end)),'ob')
    plot(wvl(pt),irradiance_data_QEPRO(pt),'or')
    title ('O2A - Interpolated Irradiance (Lin)');
    xlabel('Wavelengh (nm)');
    ylabel('Radiance (mW/m^2 sr nm)');
    %Lup
    subplot(1,2,2)
    plot(wvl,radiance_data_QEPRO,'-k')
    xlim ([740 780])
    hold on
    plot(wvl(wvl_index_peaks_out1(end)),radiance_data_QEPRO(wvl_index_peaks_out1(end)),'ob')
    plot(wvl(pt),radiance_data_QEPRO(pt),'or')   
    title ('O2A - Interpolated Radiance (Lup)');
    xlabel('Wavelengh (nm)');
    ylabel('Radiance (mW/m^2 sr nm)');
    
    %save name
    name_split      = strsplit(string(time),' ');    
    day             = strsplit(string(name_split{1}),'/');    
    day             = string([day{1} '_'  day{2} '_'  day{3}]);
    hour            = strsplit(string(name_split{2}),':');    
    hour            = string([hour{1} '_'  hour{2} '_'  hour{3} '_'  name_split{3}]);   
    name_save_1   = [day{1},'_',hour{1},'sFLD_Lin_Lup.png'];
    
    %save
    saveas(fig1,[path_save_fig name_save_1]);
    
    %close 
    close all;	% Close all figure windows except those created by imtool.
    imtool close all;	% Close all figure windows created by imtool.
end
end % end function