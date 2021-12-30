function [v1o,v2o,v3o,v4o,v5o,v6o,v7o,v8o]=paramsHandler(unpack,sedmodel,vi)
%
% USAGE-STYLE-1:
%
%   [fv,ks,n,m,xi,alpha,Cc,Cf]=paramsHandler(1,'vanderA',params)
%
% USAGE-STYLE-2:
%
%   params=paramsHandler(1,'vanderA',[fv ks n m xi alpha Cc Cf])
%
% Helper function to unpack and re-pack params struct for van der A model.
% This was needed as a hack to avoid "cannot be classified" error when using
% parfor in full-run assimilation tests
%
% The variable names in params struct change depending on 'sedmodel'.  The
% example above is for sedmodel=='vanderA'.
%

if(isvector(vi))
  v1i=vi(1);
  v2i=vi(2);
  v3i=vi(3);
  v4i=vi(4);
  v5i=vi(5);
  v6i=vi(6);
  if(length(vi)>6)  % vanderA has 8 params
    v7i=vi(7);
    v8i=vi(8);
  end

if(unpack)
  params=v1i;
  v1o=params.fv   ;
  v2o=params.ks   ;
  if(strcmp(sedmodel,'vanderA'))
    v3o=params.n    ;
    v4o=params.m    ;
    v5o=params.xi   ;
    v6o=params.alpha;
    v7o=params.Cc   ;
    v8o=params.Cf   ;
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
    params.Cc   =v7i;
    params.Cf   =v8i;
  elseif(strcmp(sedmodel,'dubarbier'))
    params.Cw=v3i;
    params.Cc=v4i;
    params.Cf=v5i;
    params.Ka=v6i;
  end
  v1o=params;
end
