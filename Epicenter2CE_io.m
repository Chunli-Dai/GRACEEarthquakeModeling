% Compute the ground surface of Japan region to the center of earth
% Epicenter by USGS is given in WGS84 reference frame, 
% check http://earthquake.usgs.gov/earthquakes/glossary.php
% Depth may be relative to mean sea-level or the average elevation of the seismic stations
% which provided arrival-time data for the earthquake location.

data=load('tmpin');
latgcd=data(2); %38.05; %geocentric coordinate in degree
%hothi=data(1)*1e3; %3.9433e3; % ocean thickness in meter
%elev=-hothi;
elev=data(1); % Elevation of surface or ocean bottom

Re=6371e3;
% lat=37.52*pi/180. ;% epicenter latitude, geocentric coordinate
latgc=latgcd*pi/180.;
%WGS84 parameter
a=6378137.0 ; %meter
f=1./298.257223563;
e2=2*f-f^2;
% WGS84 geodetic latitude:
lat=atan(tan(latgc)/(1.-e2 ));%lat=atan(tan(35.3*pi/180.)/(1.-e2 ));%lat=atan(tan(40.37190*pi/180.)/(1.-e2 ));
z=a*(1-e2)*sin(lat)/sqrt(1-e2*(sin(lat))^2);
x=a*cos(lat)/sqrt(1-e2*(sin(lat))^2);
Rellip=sqrt(x^2+z^2);
NGeoid=0.;
Rgeoid=Rellip+NGeoid;
% Rgeoid=6370e3;
Rob=Rgeoid+elev;
Remob=Re-Rob;
fprintf('Rgeoid elev Rob Remob %9.3f  %9.3f %9.3f %9.3f',Rgeoid,elev,Rob,Remob);
fid=fopen('tmp','w');
fprintf(fid,'%9.3f %9.3f',Remob,Rob);
fclose(fid);

% format short g
% disp([ Rgeoid Rob Remob])
% a*sqrt((1-2*e2*(sin(lat))^2+e2^2*(sin(lat))^2)/(1-e2*(sin(lat))^2))


