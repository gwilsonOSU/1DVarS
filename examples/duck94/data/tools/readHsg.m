function [data,meta]=readHsg(hsgfn)
%
% [H,meta]=readHsg(hsgfn)
%
% Reads Duck94 SPUV 'hsg' file, time-averaged significant wave height data
%
% OUTPUT: an array of structs:
%
% H(i).name - string with name of instrument, e.g. 'p02'
% H(i).x    - cross-shore location, m
% H(i).y    - along-shore location, m
% H(i).zobs - sensor z-position, m
% H(i).data - 10x 512s Hsig data points
% H(i).tstart_est - 10x start-times for the data points
%

fid=fopen(hsgfn);

% metadata line, contains four numbers:
%  - the first number tells how many sensors are in the file
%  - the second number tells how many sensors (irrelevant)
%  - the third number tells how many 512-s values are listed (always=10)
%  - the fourth number is our date/time, mmddhhhh
meta.nsensor=fscanf(fid,'%f',1);
fscanf(fid,'%f',1);  % ignore
meta.ndata  =fscanf(fid,'%f',1);
dstr=fscanf(fid,'%s',1);
meta.dnum_est=datenum(['1994-' dstr(1:2) '-' dstr(3:4) ' ' ...
                       dstr(5:6) ':' dstr(7:8)]);
fgetl(fid);  % newline

% now ignore some un-needed lines... comments show an example
fgetl(fid);  %  SENS 
fgetl(fid);  %  XLOC  YLOC  ZLOC
fgetl(fid);  %  TIM
fgetl(fid);  %     0.    0.    0.
fgetl(fid);  %   .160809E+06  .160826E+06  .160843E+06  .160860E+06  .160877E+06
fgetl(fid);  %   .160894E+06  .160911E+06  .160928E+06  .160945E+06  .160962E+06

% Define time stamps for the 10 data points per sensor.  From
% SteveE's readme file:
% 
% The first mean is the mean from top of hour + 16 seconds until 512 s later.
% Second number is the mean from 1024+16 s after top of hour until 512 s later.
%
for i=1:10
  tstart_est(i)=meta.dnum_est+((i-1)*1024+16)/24/60/60;
end
averagingTimeSec=512;

% now extract data for each sensor
i=0;
while(~feof(fid))
  i=i+1;
  this=struct;
  this.name=strtrim(fgetl(fid));
  this.x   =fscanf(fid,'%f',1);
  this.y   =fscanf(fid,'%f',1);
  this.zobs=fscanf(fid,'%f',1)/100;
  fgetl(fid);  % newline
  this.data=fscanf(fid,'%e',meta.ndata)/100;
  fgetl(fid);  % newline
  this.tstart_est=tstart_est;
  this.averagingTimeSec=averagingTimeSec;
  data(i)=this;
end
