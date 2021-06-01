function [out1,out2,out3,out4,out5]=readDofs(doffn)
%
% instr=readDofs(doffn)
%
% OR...
%
% [s,p,u,v,t]=readDofs(doffn)
%
% Reads Duck94 SPUV 'dof' file, to get instrument locations
%
% OUTPUT: array of structs 'instr', with the following fields:
%
% instr(i).name - string with name of instrument, e.g. 'p02'
% instr(i).x    - cross-shore location, m
% instr(i).y    - along-shore location, m
% instr(i).zobs - sensor z-position, m
% instr(i).zbed - bed z-position, m
%
% If multiple outputs are requested, then 'instr' is split into different
% struct-arrays, one for each data type s,p,u,v,t
%

fid=fopen(doffn);
i=0;
while(~feof(fid))
  i=i+1;
  this=struct;
  this.name=strtrim(fgetl(fid));
  this.x   =fscanf(fid,'%f',1);
  this.y   =fscanf(fid,'%f',1);
  this.zobs=fscanf(fid,'%f',1)/100;
  this.zbed=fscanf(fid,'%f',1)/100;
  fgetl(fid);  % trailing newline
  instr(i)=this;
end

% sometimes last one is empty, suspect there are trailing newlines.  Roll
% back until I find the last non-empty entry
while(isempty(instr(i).zbed))
  i=i-1;
end
instr=instr(1:i);

% OPTIONAL: user can select "v2" output, see header comments.  In this
% output, outputs are separated into different instrument types
if(nargout>2)
  insttype='spuvt';
  for i=1:length(insttype)
    clear thisdata id
    cnt=0;
    for j=1:length(instr)
      if(instr(j).name(1)==insttype(i))
        cnt=cnt+1;
        thisdata(cnt)=instr(j);
        id(cnt)=str2num(instr(j).name(2:end));
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
else
  out1=instr;
end
