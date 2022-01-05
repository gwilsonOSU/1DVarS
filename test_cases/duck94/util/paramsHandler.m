function [v1o,v2o,v3o,v4o,v5o,v6o,v7o,v8o,v9o]=paramsHandler(unpack,sedmodel,vi)
%
% USAGE-STYLE-1:
%
%   [fv,ks,lambda,n,m,xi,alpha,Cc,Cf]=paramsHandler(1,'vanderA',params)
%
% USAGE-STYLE-2:
%
%   params=paramsHandler(1,'vanderA',[fv ks lambda n m xi alpha Cc Cf])
%
% Helper function to unpack and re-pack params struct for van der A model.
% This was needed as a hack to avoid "cannot be classified" error when using
% parfor in full-run assimilation tests
%
% The variable names in params struct change depending on 'sedmodel'.  The
% example above is for sedmodel=='vanderA'.
%

if(length(vi)>1)
  v1i=vi(1);
  v2i=vi(2);
  v3i=vi(3);
  v4i=vi(4);
  v5i=vi(5);
  v6i=vi(6);
  v7i=vi(7);
  if(length(vi)>7)  % vanderA has 9 params
    v8i=vi(8);
    v9i=vi(9);
  end
end
if(unpack)
  params=vi;
  v1o=params.fv   ;
  v2o=params.ks   ;
  v3o=params.lambda;
  if(strcmp(sedmodel,'vanderA'))
    v4o=params.n    ;
    v5o=params.m    ;
    v6o=params.xi   ;
    v7o=params.alpha;
    v8o=params.Cc   ;
    v9o=params.Cf   ;
  elseif(strcmp(sedmodel,'dubarbier'))
    v4o=params.Cw;
    v5o=params.Cc;
    v6o=params.Cf;
    v7o=params.Ka;
  end
else
  params.fv=v1i;
  params.ks=v2i;
  params.lambda=v3i;
  if(strcmp(sedmodel,'vanderA'))
    params.n    =v4i;
    params.m    =v5i;
    params.xi   =v6i;
    params.alpha=v7i;
    params.Cc   =v8i;
    params.Cf   =v9i;
  elseif(strcmp(sedmodel,'dubarbier'))
    params.Cw=v4i;
    params.Cc=v5i;
    params.Cf=v6i;
    params.Ka=v7i;
  end
  v1o=params;
end
