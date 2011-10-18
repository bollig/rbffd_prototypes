function [X,Y]=pr_mollweide(Long,Lat,R);
%--------------------------------------------------------
% pr_mollweide function    Project coordinates to the
%                     equal area Mollweide projection
%                     for set of longitude and latitude.
% Input  : - vector of Longitude, in radians.
%          - vector of Latitude, in radians.
%          - Scale radius, default is 1.
% Output : - vector of X position
%          - vector of Y position 
%    By : Eran O. Ofek             July 1999     
%--------------------------------------------------------
if (nargin==3),
   % no default
elseif (nargin==2),
   R = 1;
else
   error('Illigal number of argument');
end

% R is really R.*S (R-radius, S-scale factor)
X = 2.*Long.*sqrt(2).*R.*cos(Lat)./pi;
Y = sqrt(2).*R.*sin(Lat);