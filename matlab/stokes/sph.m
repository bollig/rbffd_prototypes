% Spherical harmonics.
function Y = sph(l,m,th,lam)

% Flatten and transpose th and lam so they work with the legendre function
sz = size(th); th = th(:)'; lam = lam(:)';

% Normalization
a = sqrt((2*l+1)/2/pi*factorial(l-abs(m))/factorial(l+abs(m)));
Y = legendre(l,sin(th));
% Get the right associated legendre function
Y = squeeze(Y(abs(m)+1,:,:));
% Determine if the cos or sin term should be added.
pos = abs(max(0,sign(m+1)));
% Compute the spherical harmonic
Y = (pos*cos(m*lam) + (1-pos)*sin(m*lam)).*(a*Y);
% Reshape so it is the same size as the th and lam that were passed in.
Y = reshape(Y,sz);