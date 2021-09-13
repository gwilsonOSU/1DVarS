function out=mergestruct(in1,in2)
%
% out=mergestruct(in1,in2)
%
% Merges two structs together.  All fields in either 'in1' or 'in2' will be
% copied into 'out'.  In cases where the field exists in both 'in1' and
% 'in2', the field from 'in2' is copied into 'out.
%

out=in2;
fld=fields(in1);
for i=1:length(fld)
  if(isfield(in2,fld{i}))
    ;  % already copied from in2
  else
    out=setfield(out,fld{i},getfield(in1,fld{i}));
  end
end
