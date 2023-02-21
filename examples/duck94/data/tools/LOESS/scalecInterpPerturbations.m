function [zi, msei, nmsei, msri, HWAITBAR] = scalecInterpPerturbations(x, z, s, xi, lx, filtername, nmseitol,Ntotal,Ndone, HWAITBAR);
%
% [zi, msei, nmsei, msri] = scalecInterpPerturbations(x, z, s, xi, lx, filtername, nmseitol,Ntotal,Ndone);
% 
% DOES NOT REMOVE A POLYNOMIAL TREND, AS TREND IS ASSUMED ALREADY REMOVED FROM
% PERTURBATION DATA, DEFAULTS TO ZERO VALUE IF NO DATA
%
% Input
%   x, the nxm location of the data- repeated indices not handled well (ignored)
%   z, the observations MINUS a TREND surface
%   s, the the observation errors (i.e., standard deviations, rms, etc.)
%      s is used to weight the observations as 1/s
%      choose s=0 if no information is available about the weights 
%   xi, the interpolation locations
%   lx, the interpolation weighting length scales
%   filtername, the name of a filter to analyze:
%      'quadloess'
%      'linloess'
%      'hanning'
%      'boxcar'
%   nmseitol, a maximum error level, if exceeded causes doubling of smoothing scales
%       NOTE: if nmseitol=1 then this means we accept result with input scales
%             if nmseitol<1 then this means interpolation will successive doubling of scales to reach desired noise reduction
%   Ntotal, the total number of interpolated points being processed (e.g., larger than N if doing tiles)
%   Ndone, the total number of interpolated points already processed (e.g., number done in previous tiles)
% 
% Output
%   zi, the estimate
%   msei, the mean square interpolation error estimate (units of z)
%   nmsei, the normalized mean square error
%   msri, the mean square residuals

% SUGGESTED MODIFICATIONS
% 1. remove unnecessary scaling code
% 2. used consistentWeight.m to weight by the observation errors
% 3. use trick for returning observation errors if ndof is small (see approach in consistentWeight.m)

% 17 jun 2004: modified nmsei such that if one obs was encountered, nmsei is 1-eps, such that post processors know that grid point was influenced by data
% 26 oct 2004: modified by Todd for HWAITBAR issue with compiler

%  keep the waitbar open - now done by passing to routine since persistent doesn't work in compiler
% persistent HWAITBAR
% clear the waitbar if no inputs - done in other routine

dowaitbar=1;

% deal with input
[N,m] = size(x);
[Ni,mi] = size(xi);
if (nargin<8 | isempty(Ntotal)) %(~exist('Ntotal'))
    Ntotal = Ni;
end
if (nargin<9 | isempty(Ndone)) %(~exist('Ndone'))
    Ndone = 0;
end

% fix up s, if constant
[ns,ms]=size(s);
if(max([ns,ms])==1)
    s = repmat(s, N,1);
else
    s = s(:);
end

% deal with nans
id = find(isfinite(sum([x,z,s],2)));
x=x(id,:);
z = z(id);
s = s(id);
if(length(lx)==N & length(lx)~=Ni)
    % got to remove the corresponding scales
    lx = lx(id,:);
end
[N,m] = size(x);

% use weighted calculations on DATA
wtol = 0.1; % tolerance for iterative convergence of weights
s = s.^2; % need variance, not std
[wt, var_z] = consistentWeight(z, s, wtol);
% normalize weights
wt = (wt+eps)./(eps+max(wt));

% eliminate useless variables
std_x = std(x);
id =find(std_x==0);
if(length(id)>0) % catch variables with zero variance (e.g., a profile)
    std_x(id) = 1+0*id;
end

% get the convolution kernel at adequate resolution
RMAX = 1;
DR = 1e-2; % need to descretize the continous kernel function   
if (nargin<7 | isempty(nmseitol)) % (~exist('nmseitol'))
    nmseitol = inf; % never invoke tolerance
    fprintf('scalecInterpPerturbations: Setting maximum nmse tolerance to %.3f\n', nmseitol);
end    

% do the right thing for the right name
switch filtername
    case 'quadloess'
        % MUST DO THIS CORRECTLY FOR N-D
        [ri,ai] = loess_kernelND(m,2);
        FC = 0.7; % here is the half-power point
    case 'linloess'
        [ri,ai] = loess_kernelND(m,1);
        FC = 0.4; % here is the half-power point
    case 'hanning'
        ri = [0:DR:(RMAX+DR)]';
        ai = hanning_wt(ri);
        FC = 0.4; 
    case 'boxcar'
        ri = [0:DR:(RMAX+DR)]';
        ai = ones(length(ri),1);
        FC = 0.25;
    case 'si'
        ri = [0:DR:(RMAX+DR)]';
        ai = ones(length(ri),1);
        FC = 1;% ?
    otherwise
	disp(sprintf('Filter window %s not found',filtername));
        return;
end
% fprintf('interpolation method is %s\n', filtername)

% initialize output
% note that we don't know a consistent value for variance, so case of no data gets one
zi = zeros(Ni,1);
nmsei = ones(Ni,1);
%msei = var_z*ones(Ni,1);
msei = nmsei;
%msri = msei;
msri = nmsei;

% scale the data for interpolation
if(length(lx)==m)
    % fprintf('using constant smoothing scales\n')
    % scale obs. and interp. locations
    L = diag(1./lx);
    x = x*L;
    xi = xi*L;
else
    % fprintf('using variable smoothing scales\n')
    % scaling takes place below
end

% show progress bar (lives from one call to next
if(dowaitbar)
  if(Ndone==0)
    HWAITBAR = waitbar(0,'Interpolating...');
    waitpos = get(HWAITBAR,'position'); % put it off to the side
    %set(HWAITBAR, 'position', waitpos.*[0.05 (1-0.05) 1 1],'tag', 'scalecInterpWaitbar');
    set(HWAITBAR, 'position', waitpos.*[0.05 (1-0.05) 1 1]);
  end
end

tic
for i=1:Ni
  if(dowaitbar)
    waitbar((i+Ndone)/Ntotal, HWAITBAR, sprintf('Finished %3d percent',round((i+Ndone)/Ntotal*100)));
  end

    % center on point of interest
    y = x - ones(N,1)*xi(i,:);
    
    % now scale if using individual smoothing scales
    if(length(lx)==m)
        % do nothing
        if(i==1)
            %fprintf('no further scaling applied\n')
        end
    elseif(length(lx)==N)
        % scale the data
        if(i==1)
            %fprintf('scaling applied at each observation location-- scale varies on observational array\n')
        end
        y = y./lx;
    elseif(length(lx)==Ni)
        % scale at the interp. location
        if(i==1)
            %fprintf('scaling applied at each interpolation location--scale varies on interpolation array\n')
        end
        y = y*diag(1./lx(i,:));
    else
        % did nothing yet, an error
        error(sprintf('smoothing scales not interpreted: N=%d,Ni=%d,Lx=%dx%d',N,Ni,size(lx)));
    end
    
    % convert input to radial
    r = sqrt(sum(y.^2,2));
    
    % expand smoothing scales until tolerance met
    cnt = 0;
    while((nmsei(i)>nmseitol & cnt<10) | cnt==0) % do it once, at least
        p = (2^cnt); % scale r by factor of 2 each time to catch more data
        cnt = cnt + 1;
        
        % allow for computationally intensive methods here
        switch filtername
            case 'si'
                % there is not an abrupt cutoff with si, so feed it more data
                aid = find(r<(8*p));  
                na = length(aid);
                a=zeros(na,1); % default to first guess (also called norm)
                % got to compute full thing
                a = si_wt(y(aid), wt(aid)); % note, wt is %var that is signal
            otherwise
                aid = find(r<p); 
                na = length(aid);
                a=zeros(na,1); % default to first guess (also called norm)
                a = interp1(ri,ai, r(aid)/p, 'linear*');
                % apply a priori weights
                a = a.*wt(aid);
                
                % normalize if we have some weights
                suma = sum(a);
                if(abs(suma)>0)
                    a=(a/suma);
                end
        end
        
        % check weights
        if(sum(a)>0)
            % compute error, noise passed, fraction of target variance not explained
            %          (noise) + (   lost signal    )
            %nmsei(i) = (a'*a)  + 0*(1-mean(wt(aid).^2)); % no benfit found for this formulation since we normalize wts to range [0,1]
            %nmsei(i) = (a'*a); % want to distinguish points with at least one nonzero weight
            nmsei(i) = (a'*a)*(1-eps); % want to distinguish points with at least one nonzero weight
            % actual error is the noised passed (nmsei*s) + signal reduced (1-nmsei)*s
        else
            % weights all zero
            nmsei(i) = 1;
         end
    end
    
    % account for deviation from target scales with fraction of wavenumber band used!
    q = p^(-m);
    nmsei(i) = (1-q*(1-nmsei(i)));
    
    % convolve against data 
    zi(i) = (a)'*z(aid);
    
    % get error estimate, if sensible
    if(nmsei(i)<1 & na>0)
        % and weighted residuals are
        y = (z(aid)-zi(i)).*a;
        msri(i) = (y'*y)/nmsei(i); % weighted mean square residual
        % which is not accurate if small dof, add estimate that converges to observation error for (na-1)->0
        msri(i) = ((na-1).*msri(i) + a'*s(aid))./na ;
        % weighted mean square residual
        msei(i) = msri(i)*nmsei(i)/(1-nmsei(i)); % predicted mean square error
    else
        % set it to one!
        nmsei(i) = 1;
    end
end

% kill waitbar when done with all
if(dowaitbar)
  if((i+Ndone)>=Ntotal)
    close(HWAITBAR);
  end
end
etime = toc;
%fprintf('interpolated %d points/second\n', fix(Ni/toc));
