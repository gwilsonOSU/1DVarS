function [v1o,v2o,v3o,v4o,v5o,v6o,v7o,v8o,v9o]=paramsHandler(unpack,sedmodel,vi)
%
% USAGE-STYLE-1:
%
%   [fv,ks,lambda,n,m,xi,alpha,Cc,Cf] = paramsHandler(1,'vanderA',params_struct)
%
% USAGE-STYLE-1b:
%
%   params_vec = paramsHandler(1,'vanderA',params_struct);  % vector has same order of params as above
%
% USAGE-STYLE-2:
%
%   params_struct = paramsHandler(0,'vanderA',[fv ks lambda n m xi alpha Cc Cf])
%
% 1st input 0 or 1 says whether to:
%    0) pack the inputs into a struct for output, or
%    1) unpack the input-struct and output scalar parameters
%
% Helper function to unpack and re-pack params struct for van der A model.
% This is needed to avoid "cannot be classified" errors when using parfor in
% full-run assimilation tests.  It also helps keep the order of parameters
% clearly defined when switching between struct and vector representations,
% important for assimlation and inversion codes.
%
% The variable names in the params struct (input or output) change depending
% on 'sedmodel'.  The examples above are for sedmodel=='vanderA'.
%

if(length(vi)>1)
  v1i=vi(1);
  v2i=vi(2);
  v3i=vi(3);
  v4i=vi(4);
  v5i=vi(5);
  v6i=vi(6);
  v7i=vi(7);
  if(length(vi)>7)  % vanderA has up to 9 params
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
    if(isfield(params,'Cc'))
      v8o=params.Cc   ;
      v9o=params.Cf   ;
    end
  elseif(strcmp(sedmodel,'dubarbier'))
    v4o=params.Cw;
    v5o=params.Cc;
    v6o=params.Cf;
    v7o=params.Ka;
  end
  if(nargout==1)  % vector output
    v1o=[v1o v2o v3o];
    if(strcmp(sedmodel,'vanderA'))
      v1o=[v1o v4o v5o v6o v7o];
      if(isfield(params,'Cc'))
        v1o=[v1o v8o v9o];
      end
    elseif(strcmp(sedmodel,'dubarbier'))
      v1o=[v1o v4o v5o v6o v7o];
    end
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
    if(exist('v8i'))
      params.Cc   =v8i;
      params.Cf   =v9i;
    end
  elseif(strcmp(sedmodel,'dubarbier'))
    params.Cw=v4i;
    params.Cc=v5i;
    params.Cf=v6i;
    params.Ka=v7i;
  end
  v1o=params;
end
