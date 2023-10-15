function[F, rsq]=SFM_linear(wvl,irradiance_data_QEPRO,radiance_data_QEPRO)                   
% prepare input
xdata = [wvl irradiance_data_QEPRO]';
ydata = radiance_data_QEPRO';

% SFM fitting with linear function
x0 = [0.1,0.0005,1,-0.001];
lb = [-inf,0,-inf,-inf];
ub = [inf,inf,inf,0];
options = optimset('TolX',1e-10,'TolFun',1e-10,'Display','off'); 

[x,SSresid] = lsqcurvefit(@radtrans,x0,xdata,ydata,lb,ub,options);
SStotal = (length(ydata)-1) * var(ydata);

F = 760.0*x(4)+x(3); 
rsq = 1 - SSresid/SStotal;
end
