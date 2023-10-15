function [f_wvl_SFM,r_wvl_SFM]=SFM_A(wvl,irradiance_data_QEPRO,radiance_data_QEPRO,do_plot,time,path_save_fig,alg)
%% RESAMPLE TO REGULAR GRID
% wvl_all=wvl;
% Lin_all=Lin(:,n_repli);
% Lup_all=Lup(:,n_repli);
% wvl=wvl(pos1A:pos2A);
% Lin=Lin(pos1A:pos2A,n_repli);
% Lup=Lup(pos1A:pos2A,n_repli);

[~, pos1A] = min(abs(wvl - 750));
[~, pos2A] = min(abs(wvl - 780));
wvl1 = wvl(pos1A:pos2A);
ssi=diff(wvl1); ssiMin=min(ssi);%resampling
wvlRes=wvl1(1):ssiMin:wvl1(end); wvlRes=wvlRes';
Lup   = interp1(wvl1,Lup,wvlRes);
Lin   = interp1(wvl1,Lin,wvlRes);
%Define weight
w = abs((Lin-max(Lin(:))) ./ (max(Lin(:)-min(Lin(:)))));
%w = ones(size(Lin));

%% SCALE X-AXES
dwvl = wvl1-wvl1(1); %no lo utilizo

%% APPARENT REFLECTANCE  (poly2 - 3 degree freedom) 
Rapp=Lup./Lin;
% EXCLUDE THE O2 ABSORPTION BAND
index_wvl_interpolation_R = find(wvlRes>759&wvlRes<770);
%Define interpolation range/values
dwvlNoABS = wvlRes; dwvlNoABS(index_wvl_interpolation_R) =[];
RappNoABS = Rapp; RappNoABS(index_wvl_interpolation_R) = [];
%Spline interpolation
% knots vector for 15 piecewise spline
knots = aptknt(dwvlNoABS(round(linspace(1,numel(dwvlNoABS),3-1+4))),5); % 4
sp = fastBSpline.lsqspline(knots,3,dwvlNoABS,RappNoABS);
p = sp.weights'; 
%Equation
Ref_app = sp.evalAt(wvlRes);
% Plot        
% plot(wvlRes,Rapp)
% hold on
% plot(wvlRes,Ref_app) 

%% FIRST GUESS - FLUORESCENCE 
%First guess come from PEAK_HEIGHT retrieval method
[fluoFG]=SFM_iFLD_A(wvl,irradiance_data_QEPRO,radiance_data_QEPRO);
if (~isnan(fluoFG)) && (fluoFG>=0) && (~isnan(p(2)))
%% BOUNDARIES 
%first guess
FG = [fluoFG,24,p];
%lower boundary
%lb = [0,0,p];
lb = [0,0,-Inf(1,numel(p))];
%upper boundary
%ub = [15,Inf,p];
ub = [15,Inf,Inf(1,numel(p))];

%% OBJECTIVE FUNCTION + PROBLEM SCALING
%weigths
Fw = @(x,Lin)SFM_funct_A(x,Lin,wvlRes,w,sp); 

%weighting 
Lup_w = Lup .* w; % con pesi

switch alg
        case 'tr'
            options = optimset('Display','iter','MaxFunEvals',50,'TolFun',1e-12,'TolX',1e-15);
            [x,resnorm,residual,exitflag,output] = lsqcurvefit(Fw,FG,Lin,Lup_w,lb,ub,options);
            if exitflag==-2; resnorm=nan(1);end;
        case 'lm'
            options = optimset('Algorithm','levenberg-marquardt', 'Display','iter','MaxFunEvals',300, 'ScaleProblem','jacobian');
            [x,resnorm,residual,exitflag,output] = lsqcurvefit(Fw,FG,Lin,Lup_w,[],[],options);
            if exitflag==-2; resnorm=nan(1);end;
        case 'nlf'
            options = optimset('Algorithm','levenberg-marquardt', 'Display','final','TolFun',1e-18, 'TolX', 1e-18, 'ScaleProblem','jacobian');
            [x,resnorm,residual,exitflag,output] = lsqcurvefit(Fw,FG,Lin,Lup_w,[],[],options);
        otherwise
            disp('check the Optimization algorithm');
    end
%% Info optimset function
%https://www.mathworks.com/help/optim/ug/optimization-options-reference.html[TolFun and TolX]
%https://www.mathworks.com/help/optim/ug/tolerances-and-stopping-criteria.html[TolX]]
%% output R and F spectra
% 3 grados de libertad
%% Output R and F spectra
% Fs
f_wvl_SFM = x(1).*exp(-(((wvlRes-wvlRes(1))-(740-wvlRes(1))).^2/(2*x(2).^2)));
f_wvl_SFM = f_wvl_SFM';

% Re
sp.weights=x(3:end);
r_wvl_SFM = sp.evalAt(wvlRes);

%% Plots "measured" and "modeled" radiance and residuals
if do_plot == 1
    Lup_ret = SFM_funct_A(x,Lin,wvlRes,w,sp); 
    residual = (Lup.*w) - Lup_ret; resnorm = sum(residual.^2);
    %figure 1
    fig1 = figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(1,2,1)
    plot(wvlRes,Lup.*w,'k')
    hold on 
    plot(wvlRes,Lup_ret,'-.r')
    legend('Lup','Lup-ret')
    title('Lup retrived vs Lup measured')
    %figure 2
    subplot(1,2,2)
    plot(wvlRes,residual,'r')
    title('Lup retrived vs Lup measured residual')
   
    %save name
    name_split      = strsplit(string(time),' ');    
    day             = strsplit(string(name_split{1}),'/');    
    day             = string([day{1} '_'  day{2} '_'  day{3}]);
    hour            = strsplit(string(name_split{2}),':');    
    hour            = string([hour{1} '_'  hour{2} '_'  hour{3} '_'  name_split{3}]);   
    name_save_1   = [day{1},'_',hour{1},'_SFM_A_Lup_MvR_Residual.png'];
    
    %save
    saveas(fig1,[path_save_fig name_save_1]);
    
    %close 
    close all;	% Close all figure windows except those created by imtool.
    imtool close all;	% Close all figure windows created by imtool.
end
else
    f_wvl_SFM=nan(1,length(wvlRes));
    r_wvl_SFM=nan(length(wvlRes),1);
   
end