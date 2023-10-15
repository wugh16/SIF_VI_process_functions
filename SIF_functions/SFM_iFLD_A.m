function [F] = SFM_iFLD_A(wvl_all,Lin_all,Lup_all);
%input data
wvl = wvl_all;
irradiance_data   = Lin_all;
radiance_data   = Lup_all;

%% Step 1:O2-A - Define point outside absorption bands 
% Points for interpolation
[~, pos1A] = min(abs(wvl - 750));
[~, pos2A] = min(abs(wvl - 780));

%Input data O2A
wvl_inter = wvl(pos1A:pos2A);
E = irradiance_data(pos1A:pos2A);
L = radiance_data(pos1A:pos2A);

% EXCLUDE THE O2 ABSORPTION BAND
index_wvl_inside = find(wvl_inter>759&wvl_inter<770);
index_wvl_outside_1 = find(wvl_inter<759);
index_wvl_outside_2 = find(wvl_inter>770);
index_wvl_outside = [index_wvl_outside_1;index_wvl_outside_2];

%% Step 2: Irradiance (E) - define max point outside abs. band
%wvl_index_out
%Eliminate absorption band from ARHO data
E_NoABS = E; E_NoABS(index_wvl_inside) = [];
L_NoABS = L; L_NoABS(index_wvl_inside) = [];
%Eliminate absorption band from wvl data
dwvlNoABS = wvl_inter; dwvlNoABS(index_wvl_inside) =[];
%find peaks
[value_E1,pks_E1] = findpeaks(E_NoABS);

%% Step 3: IRRADIANCE (E) FIRST INTERPOLATION - Calculate and define the inputs needed to retrieve F using the peak height approach
%E_INTERPOLATED_1
[p,S,mu] = polyfit(dwvlNoABS(pks_E1),E_NoABS(pks_E1),2);
[E_INTERPOLATED_1,err_pol]  = polyval(p,wvl_inter,S,mu);

%% Step 4:  Find interpolated points that are lower than incomming irradiance
%Calculate difference between E_points_used_interpolation and E_points_result_Interpolation
E_INTERPOLATED_1_ABS = E_INTERPOLATED_1; E_INTERPOLATED_1_ABS(index_wvl_inside) = [];
diff_E_used_E_result = E_NoABS(pks_E1) - E_INTERPOLATED_1_ABS(pks_E1);

%Find where diff_E_used_E_result is positive - it is where E_points_used_interpolation > E_points_result_Interpolation
diff_E_used_E_result_positive = find (diff_E_used_E_result >= 0);

%define final wvl interpolation 
index_E_outside_final = pks_E1(diff_E_used_E_result_positive);

%% Step 5: SECOND INTERPOLATION - Calculate and define the inputs needed to retrieve F using the peak height approach
%E_INTERPOLATED_2
[p,S,mu] = polyfit(dwvlNoABS(index_E_outside_final),E_NoABS(index_E_outside_final),2);
[E_INTERPOLATED_2,err_pol]  = polyval(p,wvl_inter,S,mu);

%% Step 6: APPARENT REFLECTANCE (ARHO) FIRST INTERPOLATION - Calculate and define the inputs needed to retrieve F using the peak height approach
%Using E maximuns
%Calculate ARHO
ARHO = L./E;
%Define interpolation range/values
%Eliminate absorption band from ARHO data
RappNoABS = ARHO; RappNoABS(index_wvl_inside) = [];
RappNoABS_peak = RappNoABS(pks_E1); %index_E_outside_final
%Eliminate absorption band from wvl data
dwvlNoABS = wvl_inter; dwvlNoABS(index_wvl_inside) =[];
dwvlNoABS_peak = dwvlNoABS(pks_E1);%index_E_outside_final

%Spline interpolation ARHO_INTERPOLATED_1
% knots vector for 15 piecewise spline
knots = aptknt(dwvlNoABS_peak(round(linspace(1,numel(dwvlNoABS_peak),3-1+4))),5); % 4
sp = fastBSpline.lsqspline(knots,3,dwvlNoABS_peak,RappNoABS_peak);
p = sp.weights'; 
%Equation
ARHO_INTERPOLATED_1 = sp.evalAt(wvl_inter);

%% Step 6: Retrieve fluorescence in 3 different points inside the abs. band 
%F = (ARHO - ARHO_INTERPOLATED).*(E_INTERPOLATED.*E)./(pi.*(E_INTERPOLATED-E))

%wvl_index_in
E_min = min(E(index_wvl_inside));
[~, index_E_in] = min(abs(E - E_min));

%Define 4 points inside abs band for wvl_E_index_in
pt1 = index_E_in;
pt2 = index_E_in-1;
pt3 = index_E_in+1;
pt4 = index_E_in+2;
pt = [pt1 pt2 pt3 pt4];

%Compute F = (ARHO - ARHO_INTERPOLATED).*(E_INTERPOLATED.*E)./(pi.*(E_INTERPOLATED-E))
ARHO_mean                = mean(ARHO(pt));
ARHO_INTERPOLATED_1_mean = mean(ARHO_INTERPOLATED_1(pt));
E_INTERPOLATED_2_mean    = mean(E_INTERPOLATED_2(pt));
E_in_mean                = mean(E(pt));
L_in_mean                = mean(L(pt));
E_out_mean               = mean(E_NoABS(index_E_outside_final));
L_out_mean               = mean(L_NoABS(index_E_outside_final));

%Reflectance factor // Fr = (L(wvl_out1)/E(wvl_out1))/ARHO_INTERPOLATED_2(pt_inter)   
Fr = (L_out_mean/E_out_mean)/ARHO_INTERPOLATED_1_mean;
%Fluorescence factor // Ff = Fr * E[wl_out1]/E_INTERPOLATED_2(pt_inter)   
Ff = Fr * E_out_mean/E_INTERPOLATED_2_mean;

%Compute fluorescence // F = (Fr*E(wvl_out1)*L(ndx_in) - E(ndx_in)*L(wvl_out1)) / (Fr*E(wvl_out1) - Ff*E(ndx_in))
F = (Fr*E_out_mean*L_in_mean - E_in_mean*L_out_mean) / (Fr*E_out_mean - Ff*E_in_mean);

%Compute reflectance // R = (L[ndx_in] - F) / E[ndx_in]
R = (L_in_mean - F) / E_in_mean;
end