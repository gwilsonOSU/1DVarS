function [out1,out2,out3,out4,out5,out6]=readMea(meafn)
%
% [data,meta]=readMea(meafn)
%
% OR...
%
% [s,p,u,v,t,meta]=readMea(meafn)
%
% Reads Duck94 SPUV 'mea' file, time-averaged (MEAn) instrument data
%
% There are 2 ways this code can be run (see above)... v1, and v2...
%
% OUTPUT, v1: array of structs 'data', with the following fields:
%
% data(i).name - string with name of instrument, e.g. 'p02'
% data(i).x    - cross-shore location, m
% data(i).y    - along-shore location, m
% data(i).zobs - sensor z-position, m
% data(i).data - 5x 512s averaged data points
%
% OUTPUT, v2: array of structs 's','p','u','v',...
%
% In this version, the code further parses the 'data' struct from v1 above,
% to isolate different sensor types.  The outputs are reformatted into
% struct arrays with same fields as the 'data' struct in v1
%

fid=fopen(meafn);

% metadata line:
%  - the first number tells how many sensors are in the file
%  - the second number tells how many 512-s values are listed (always=6)
%  - the third number is our date/time, mmddhhhh
meta.nsensor=fscanf(fid,'%f',1);
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


% OPTIONAL: user can select "v2" output, see header comments.  In this
% output, outputs are separated into different instrument types
if(nargout>2)
  insttype='spuvt';
  for i=1:length(insttype)
    clear thisdata id
    cnt=0;
    for j=1:length(data)
      if(data(j).name(1)==insttype(i))
        cnt=cnt+1;
        thisdata(cnt)=data(j);
        id(cnt)=str2num(data(j).name(2:end));
      end
    end
    [~,ind]=sort(id);
    thisdata=thisdata(ind);  % sort by instrument name
    eval([insttype(i) '=thisdata;']);
  end
  out1=s;
  out2=p;
  out3=u;
  out4=v;
  out5=t;
  out6=meta;
else
  out1=data;
  out2=meta;
end
