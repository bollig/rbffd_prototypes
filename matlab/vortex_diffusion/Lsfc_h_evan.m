
function [val] = Lsfc_h_evan(phi, th, t, rho0, gamma)
%% DERIVED BY EVAN IN MATHEMATICA:
%
%rho0 = 3; 
%gamma = 5; 

rho_p = rho0*cos(th);

% General case
val = sech((rho0.*cos(th).*sin(phi - (3.*sqrt(3).*t.*sec(th).*sech(rho0.*cos(th)).^2.*tanh(rho0.*cos(th)))./(2.0.*rho0)))./gamma).^2.*tan(th).* ...
        (-((rho0.*sin(th).*sin(phi - (3.*sqrt(3).*t.*sec(th).*sech(rho0.*cos(th)).^2.*tanh(rho0.*cos(th)))./(2.0.*rho0)))./gamma) + ...
          (rho0.*cos(th).*cos(phi - (3.*sqrt(3).*t.*sec(th).*sech(rho0.*cos(th)).^2.*tanh(rho0.*cos(th)))./(2.0.*rho0)).* ...
             ((3.*sqrt(3).*t.*sech(rho0.*cos(th)).^4.*tan(th))./2.0 - (3.*sqrt(3).*t.*sec(th).*sech(rho0.*cos(th)).^2.*tan(th).*tanh(rho0.*cos(th)))./(2.0.*rho0) - ...
               3.*sqrt(3).*t.*sech(rho0.*cos(th)).^2.*tan(th).*tanh(rho0.*cos(th)).^2))./gamma) -  ...
       sech((rho0.*cos(th).*sin(phi - (3.*sqrt(3).*t.*sec(th).*sech(rho0.*cos(th)).^2.*tanh(rho0.*cos(th)))./(2.0.*rho0)))./gamma).^2.* ...
        (-((rho0.*cos(th).*sin(phi - (3.*sqrt(3).*t.*sec(th).*sech(rho0.*cos(th)).^2.*tanh(rho0.*cos(th)))./(2.0.*rho0)))./gamma) - ...
          (2.*rho0.*cos(phi - (3.*sqrt(3).*t.*sec(th).*sech(rho0.*cos(th)).^2.*tanh(rho0.*cos(th)))./(2.0.*rho0)).*sin(th).* ...
             ((3.*sqrt(3).*t.*sech(rho0.*cos(th)).^4.*tan(th))./2.0 - (3.*sqrt(3).*t.*sec(th).*sech(rho0.*cos(th)).^2.*tan(th).*tanh(rho0.*cos(th)))./(2.0.*rho0) - ...
               3.*sqrt(3).*t.*sech(rho0.*cos(th)).^2.*tan(th).*tanh(rho0.*cos(th)).^2))./gamma - ...
          (rho0.*cos(th).*sin(phi - (3.*sqrt(3).*t.*sec(th).*sech(rho0.*cos(th)).^2.*tanh(rho0.*cos(th)))./(2.0.*rho0)).* ...
             ((3.*sqrt(3).*t.*sech(rho0.*cos(th)).^4.*tan(th))./2.0 - (3.*sqrt(3).*t.*sec(th).*sech(rho0.*cos(th)).^2.*tan(th).*tanh(rho0.*cos(th)))./(2.0.*rho0) - ...
                3.*sqrt(3).*t.*sech(rho0.*cos(th)).^2.*tan(th).*tanh(rho0.*cos(th)).^2).^2)./gamma + ...
          (rho0.*cos(th).*cos(phi - (3.*sqrt(3).*t.*sec(th).*sech(rho0.*cos(th)).^2.*tanh(rho0.*cos(th)))./(2.0.*rho0)).* ...
             ((3.*sqrt(3).*t.*sec(th).^2.*sech(rho0.*cos(th)).^4)./2.0 + (3.*sqrt(3).*t.*sech(rho0.*cos(th)).^4.*tan(th).^2)./2.0 - ...
               (3.*sqrt(3).*t.*sec(th).^3.*sech(rho0.*cos(th)).^2.*tanh(rho0.*cos(th)))./(2.0.*rho0) + ...
              12.*sqrt(3).*rho0.*t.*sech(rho0.*cos(th)).^4.*sin(th).*tan(th).*tanh(rho0.*cos(th)) - ...
               (3.*sqrt(3).*t.*sec(th).*sech(rho0.*cos(th)).^2.*tan(th).^2.*tanh(rho0.*cos(th)))./(2.0.*rho0) - ...
               3.*sqrt(3).*t.*sec(th).^2.*sech(rho0.*cos(th)).^2.*tanh(rho0.*cos(th)).^2 - 3.*sqrt(3).*t.*sech(rho0.*cos(th)).^2.*tan(th).^2.*tanh(rho0.*cos(th)).^2 - ...
               6.*sqrt(3).*rho0.*t.*sech(rho0.*cos(th)).^2.*sin(th).*tan(th).*tanh(rho0.*cos(th)).^3))./gamma) + ...
       2.*sech((rho0.*cos(th).*sin(phi - (3.*sqrt(3).*t.*sec(th).*sech(rho0.*cos(th)).^2.*tanh(rho0.*cos(th)))./(2.0.*rho0)))./gamma).^2.* ...
        (-((rho0.*sin(th).*sin(phi - (3.*sqrt(3).*t.*sec(th).*sech(rho0.*cos(th)).^2.*tanh(rho0.*cos(th)))./(2.0.*rho0)))./gamma) + ...
           (rho0.*cos(th).*cos(phi - (3.*sqrt(3).*t.*sec(th).*sech(rho0.*cos(th)).^2.*tanh(rho0.*cos(th)))./(2.0.*rho0)).* ...
              ((3.*sqrt(3).*t.*sech(rho0.*cos(th)).^4.*tan(th))./2.0 - (3.*sqrt(3).*t.*sec(th).*sech(rho0.*cos(th)).^2.*tan(th).*tanh(rho0.*cos(th)))./(2.0.*rho0) - ...
                3.*sqrt(3).*t.*sech(rho0.*cos(th)).^2.*tan(th).*tanh(rho0.*cos(th)).^2))./gamma).^2.* ...
       tanh((rho0.*cos(th).*sin(phi - (3.*sqrt(3).*t.*sec(th).*sech(rho0.*cos(th)).^2.*tanh(rho0.*cos(th)))./(2.0.*rho0)))./gamma) + ...
       sec(th).^2.*((rho0.*cos(th).*sech((rho0.*cos(th).*sin(phi - (3.*sqrt(3).*t.*sec(th).*sech(rho0.*cos(th)).^2.*tanh(rho0.*cos(th)))./(2.0.*rho0)))./gamma).^2.* ...
             sin(phi - (3.*sqrt(3).*t.*sec(th).*sech(rho0.*cos(th)).^2.*tanh(rho0.*cos(th)))./(2.0.*rho0)))./gamma + ...
          (2.*rho0.^2.*cos(th).^2.*cos(phi - (3.*sqrt(3).*t.*sec(th).*sech(rho0.*cos(th)).^2.*tanh(rho0.*cos(th)))./(2.0.*rho0)).^2.* ...
             sech((rho0.*cos(th).*sin(phi - (3.*sqrt(3).*t.*sec(th).*sech(rho0.*cos(th)).^2.*tanh(rho0.*cos(th)))./(2.0.*rho0)))./gamma).^2.* ...
             tanh((rho0.*cos(th).*sin(phi - (3.*sqrt(3).*t.*sec(th).*sech(rho0.*cos(th)).^2.*tanh(rho0.*cos(th)))./(2.0.*rho0)))./gamma))./gamma.^2);

% handle Poles
% NOTE: when th == 0, rhop == 0 and h == 1; Lapl(h) == 0. 
val(abs(rho_p) < 4*eps) = 0;
