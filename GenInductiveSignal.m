clear
close all

load siginductiveref.mat
% pp=BSFK(x,y,3,15,[],struct('KnotRemoval','none','animation',1,'lambda',1e-6,'periodic',true,'pntcon',pntcon))

facq = 2.5e6;
Asig = 10000;
Avar = 0.1;
RPM = 3000;
nblades = 77;
nrounds = 10;
noisestd = 64;
filename = 'test.wav';

% facqorg = 25000;
% Asigorg = 2.333034178882874e+04;
% RPMorg = 1500;

dxi = (RPM*facqorg) / (RPMorg*facq);
dx = max(x)-min(x);
xlast = nblades*nrounds*dx;
xi = (0:dxi:xlast);
Env = 1-(Avar/2) * (1+sin(xi * (2*pi/(nblades*dx))));
xi = mod(xi,dx)+min(x);
yi = ppval(pp,xi) * (Asig/ Asigorg);
yi = Env.*yi + noisestd*randn(size(xi));
yi = int16(yi(:));

Curve_Preview(double(yi));

audiowrite(filename,yi(:),facq);
