function y = SFM_funct_A(x,Lin,wvl,w,sp)
%x=
% Fs - (gaussian) 
F = x(1).*exp(-(((wvl-wvl(1))-(740-wvl(1))).^2/(2*x(2).^2)));

% RHO (ploy4 - 5 degree freedom)
sp.weights=x(3:end);
R = sp.evalAt(wvl);

% Total canopy radiance
y = F + (R.*Lin);
y = y.*w;
end