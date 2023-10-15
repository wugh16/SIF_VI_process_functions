function[F,R]=FLD3_A(wvl,E,L,do_plot,time,path_save_fig)   
% wvl,,,do_plot,,path_save_fig
% E=Lin(:,n_repli);
% L=Lup(:,n_repli);
% time=up_time{n_repli};
%WVL_IN - define inside absorption band
[~, wvl_758] = min(abs(wvl - 758));
[~, wvl_770] = min(abs(wvl - 770));
%wvl range absorption band
wvl_index_range_in = (wvl_758:wvl_770);
    
%WVL_OUT - two points outside absorption band
%out1
[~, wvl_745] = min(abs(wvl - 745));
wvl_index_out1 = (wvl_745:wvl_758);
%out2
[~, wvl_780] = min(abs(wvl - 780));
wvl_index_out2 = (wvl_770:wvl_780);   

% Irradiance (E)
%Define  "E" PKS MAX outside absorption band
%wvl_index_out1
[value_E1,pks_E1] = findpeaks(E(wvl_index_out1));
wvl_E_index_out1 = wvl_index_out1(pks_E1);
%wvl_index_out2
[value_E2,pks_E2] = findpeaks(E(wvl_index_out2));
wvl_E_index_out2 = wvl_index_out2([pks_E2]);
%Define "SAMPLE" MIN inside absorption band (only one value)
E_min = min(E(wvl_index_range_in));
[~, wvl_E_index_in] = min(abs(E - E_min));

%Interpolation inside the absorprion band
%WR - retrieve fluorescence
E_out = interp1(wvl([wvl_E_index_out1(end) wvl_E_index_out2(1)]),E([wvl_E_index_out1(end) wvl_E_index_out2(1)]),wvl(wvl_E_index_in),'linear');
%WR - plot it
E_out_line = interp1(wvl([wvl_E_index_out1(end) wvl_E_index_out2(1)]),E([wvl_E_index_out1(end) wvl_E_index_out2(1)]),wvl(wvl_E_index_out1(end):wvl_E_index_out2(1)),'linear');

%sample - retrieve fluorescence
%WR - retrieve fluorescence
L_out = interp1(wvl([wvl_E_index_out1(end) wvl_E_index_out2(1)]),L([wvl_E_index_out1(end) wvl_E_index_out2(1)]),wvl(wvl_E_index_in),'linear');
%WR - plot it
L_out_line = interp1(wvl([wvl_E_index_out1(end) wvl_E_index_out2(1)]),L([wvl_E_index_out1(end) wvl_E_index_out2(1)]),wvl(wvl_E_index_out1(end):wvl_E_index_out2(1)),'linear');

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


if ~isempty(pt)
    %get wr_in and sample_in
    E_in  = mean(E(pt));
    L_in = mean(L(pt));

    %Compute Fluorescence and Reflectance
    F=(E_out*L_in-L_out*E_in)./(E_out-E_in);
    R=(L_out-L_in)./(E_out-E_in);
end

%% O2A - 760
if do_plot == 1 
    %Lin
    fig1 = figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(1,2,1)
    plot(wvl,E,'-k')
    xlim ([740 780])
    hold on
    plot(wvl(wvl_E_index_out1),E(wvl_E_index_out1),'ob')
    plot(wvl(wvl_E_index_out2),E(wvl_E_index_out2),'ob')
    plot(wvl(wvl_E_index_out1(end):wvl_E_index_out2(1)),E_out_line,'-r')
    plot(wvl(wvl_E_index_in),E_out,'or')
    plot(wvl(wvl_E_index_in),E(wvl_E_index_in),'or')
    title ('O2A - Interpolated Irradiance (Lin)');
    xlabel('Wavelengh (nm)');
    ylabel('Radiance (mW/m^2 sr nm)');
    %Lup
    subplot(1,2,2)
    plot(wvl,L,'-k')
    xlim ([740 780])
    hold on
    plot(wvl(wvl_E_index_out1),L(wvl_E_index_out1),'ob')
    plot(wvl(wvl_E_index_out2),L(wvl_E_index_out2),'ob')
    plot(wvl(wvl_E_index_out1(end):wvl_E_index_out2(1)),L_out_line,'-r')
    plot(wvl(wvl_E_index_in),L_out,'or')
    plot(wvl(wvl_E_index_in),L(wvl_E_index_in),'or')    
    title ('O2A - Interpolated Radiance (Lup)');
    xlabel('Wavelengh (nm)');
    ylabel('Radiance (mW/m^2 sr nm)');

    %save name
    name_split      = strsplit(string(time),' ');    
    day             = strsplit(string(name_split{1}),'/');    
    day             = string([day{1} '_'  day{2} '_'  day{3}]);
    hour            = strsplit(string(name_split{2}),':');    
    hour            = string([hour{1} '_'  hour{2} '_'  hour{3} '_'  name_split{3}]);   
    name_save_1   = [day{1},'_',hour{1},'3FLD_Interpolated_Lin_Lup.png'];
    
    %save
    saveas(fig1,[path_save_fig name_save_1]);
    
    %close 
    close all;	% Close all figure windows except those created by imtool.
    imtool close all;	% Close all figure windows created by imtool.
end
end % end function