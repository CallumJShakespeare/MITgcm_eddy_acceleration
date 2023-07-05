% This is a matlab script that generates the input data
prec='real*8';
ieee='b';

% Dimensions of grid
nx=1000;
ny=1000;
nz=200;
% Nominal depth of model (meters)
H=2000;
% Size of domain
Lx=500e3;
Ly=500e3;

x=((1:nx)-0.5)/(nx-1)*Lx; 
y=((1:ny)-0.5)/(ny-1)*Ly; 

z1=linspace(0,H,nz+1);
dz=diff(z1);
z=(z1(1:end-1)+z1(2:end))/2;

[X,Y,Z]=ndgrid(x,y,-z);

% physical params
tAlpha=2.43E-4;
rho0=1026.5;
N2=2e-5;
R=100e3;
h0=300;
f=1e-4;
g=9.81;
dT=2;
omega0=2*pi/12/3600;

% initial fields.
Tinit=18+N2*Z/(g*tAlpha)+dT*exp(-((X-Lx/2).^2+(Y-Ly/2).^2)/R^2).*exp(Z/h0);
vinit=1/f*h0*dT*g*tAlpha*-2/R^2*(X-Lx/2).*exp(-((X-Lx/2).^2+(Y-Ly/2).^2)/R^2).*exp(Z/h0);
uinit=-1/f*h0*dT*g*tAlpha*-2/R^2*(Y-Ly/2).*exp(-((X-Lx/2).^2+(Y-Ly/2).^2)/R^2).*exp(Z/h0);


fid=fopen('tinit.field','w',ieee); fwrite(fid,Tinit,prec); fclose(fid);
fid=fopen('uinit.field','w',ieee); fwrite(fid,uinit,prec); fclose(fid);
fid=fopen('vinit.field','w',ieee); fwrite(fid,vinit,prec); fclose(fid);


% tide: exactly 12 hours period; reference case is i=4
for i=[1,2,4,8]
F0=1e-6*i;
u0=omega0/(f^2-omega0^2)*F0;
fid=fopen(sprintf('uinit_tide_plus_eddy_%i.field',i),'w',ieee); fwrite(fid,uinit+u0,prec); fclose(fid);
end

% Goff topo
data=load('hdata.mat'); % generated from Goff (2010) spectra using parameters given in paper.
hGoff=data.h;
hGoff=hGoff-min(min(hGoff));
h=-H+hGoff;
fid=fopen('Goff_min5km_rms100.field','w',ieee); fwrite(fid,h,prec); fclose(fid);
hGoff=data.h;
hGoff=hGoff/100*50/sqrt(2);
hGoff=hGoff-min(min(hGoff));
h=-H+hGoff;
fid=fopen('Goff_min5km_rms35.field','w',ieee); fwrite(fid,h,prec); fclose(fid);



