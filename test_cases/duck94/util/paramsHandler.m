function [v1o,v2o,v3o,v4o,v5o,v6o]=paramsHandler(unpack,sedmodel,v1i,v2i,v3i,v4i,v5i,v6i)
%
% USAGE-STYLE-1:
%
%   [fv,ks,n,m,xi,alpha]=paramsHandler(1,'vanderA',params)
%
% USAGE-STYLE-2:
%
%   params=paramsHandler(1,'vanderA',fv,ks,n,m,xi,alpha)
%
% Helper function to unpack and re-pack params struct for van der A model.
% This was needed as a hack to avoid "cannot be classified" error when using
% parfor in full-run assimilation tests
%
% The variable names in params struct change depending on 'sedmodel'.  The
% example above is for sedmodel=='vanderA'.
%

if(unpack)
  params=v1i;
  v1o=params.fv   ;
  v2o=params.ks   ;
  if(strcmp(sedmodel,'vanderA'))
    v3o=params.n    ;
    v4o=params.m    ;
    v5o=params.xi   ;
    v6o=params.alpha;
  elseif(strcmp(sedmodel,'dubarbier'))
    v3o=params.Cw;
    v4o=params.Cc;
    v5o=params.Cf;
    v6o=params.Ka;
  end
else
  params.fv=v1i;
  params.ks=v2i;
  if(strcmp(sedmodel,'vanderA'))
    params.n    =v3i;
    params.m    =v4i;
    params.xi   =v5i;
    params.alpha=v6i;
  elseif(strcmp(sedmodel,'dubarbier'))
    params.Cw=v3i;
    params.Cc=v4i;
    params.Cf=v5i;
    params.Ka=v6i;
  end
  v1o=params;
end
