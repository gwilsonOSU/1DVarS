function [offset, offseterr] = computeOffset(data, lx, gridxy)
%try
% 
% Input
%   data, an array of structures (length k) with
%      data.x, the independent vars (e.g., [x,y,t], NxM)
%      data.z, the dependent var (e.g., [z], Nx1)
%      data.e, a priori rms error estimates on z (e.g., [z], Nx1)
%           (variable or constant)
%   lx, the filter width for finding intercomparison
%   gridxy, Nx2 array of grid coordinates.  computs the offset over the grid if it exists
%
% Output
%   offset, the array of k offset estimates for each data set
%   offseterr, the corresponding rms error of the estimates

% 14-dec-2003: change the approach for computing dz to matrix solution suggested by holman
% 16 jun 2004: use sparse matrices and made othe mods, including removal of plane offset
%              try to handle poor inversion properties, require interpolator to return nmsei slightly less than one if a data point was encountered

% what kind of interpolator to make grids?
filtername = 'hanning';
nmseitol = 1; % force one-pass interpolation, rather than expanding scales

% how many datasets?
K = length(data);
M = size(data(1).x,2);
if(M>2)
    error('this program cannot handle dim>2 fields')
end

% init output
offset = zeros(1,K);
offseterr = offset;
if (K<2)
    % nothing to do
    return
end

% see if we can use the grid
if nargin == 3 & ~isempty(gridxy) % used to be (exist('gridxy')==1)
    fprintf('using grid region to find offsets\n')
    maxX = max(gridxy);
    minX = min(gridxy);
    clear gridxy
    % toss data outside grid
    for k=1:K
        id = find(data(k).x(:,1)>(minX(1)-lx(1)) & data(k).x(:,1)<(maxX(1)+lx(1))); 
        if(M==2)
            id = id(find(data(k).x(id,2)>(minX(2)-lx(2)) & data(k).x(id,2)<(maxX(2)+lx(2))));
        end
        data(k).n = length(id);
        data(k).x = data(k).x(id,:);
        data(k).z = data(k).z(id);
        data(k).e = data(k).e(id);
    end
    % find data sets that survived
    kid = find([data.n]>0);
    data = data(kid);
    K = length(kid);
else
    % need to use same starting point, and ending point so get min/max of all datasets
    minX = min(cat(1,data.x));
    maxX = max(cat(1,data.x));
end

% reformat smoothing scale data into a useful measure, calculate dx
% dx will be coarse here, as we will average over cells anyhow: we are not after higher resolution than lx!
if (isstruct(lx))
    disp('reformatting lx into array variable');
    lx1 = cat(1,lx.l);
    [nlr,nlc] = size(lx1)
    for i=1:nlc
        dx(:,i) = mean(lx1(:,i));
    end
else
    dx = lx;    
end
dx = lx/2; % make coarser/finer?

% get integer range to ensure that all data enclosed within grid
nx = ceil((maxX-minX)./dx);
dx = (maxX-minX)./nx;
[X,Y] = meshgrid(minX(1):dx(1):maxX(1), minX(2):dx(2):maxX(2));
[ny,nx]=size(X);

% use interpolator to map each dataset to corresponding grid
% use sparse grids, as there are many zeros
z = sparse(1,1,0,ny*nx,K);
wt = sparse(1,1,0,ny*nx,K);
for k=1:K
    if (isstruct(lx))
        disp('performing downsampling of z, lx, ly, and lt in computeOffset');
        subdata = [data(k).z lx(k).l];
        [data(k).x,subdata,data(k).e] = subsampleData(data(k).x,subdata,data(k).e,min(lx(k).l)/2);
        data(k).z = subdata(:,1);
        lx(k).l = subdata(:,2:3);
        lx1 = lx(k).l;
    else
        lx1 = lx;
        [data(k).x,data(k).z, data(k).e] = subsampleData(data(k).x,data(k).z,data(k).e,lx1/2);
    end
    [z(:,k), msei, wt(:,k)] = scalecInterpTile(data(k).x, data(k).z, data(k).e, [X(:) Y(:)], lx1, filtername, nmseitol);
    wt(:,k) = 1-wt(:,k);
    % set worthless points back to zero to keep matrices sparse
    id = find(wt(:,k)==0);
    z(id,k) = 0;
    disp(sprintf('Finished interpolating %d of %d sets',k,K));
end

% form the matrix equation: DZ(k,k') = [Delta(k) - Delta(k')]dz(k)=>D = [R]*dz
% here are observed differences between data sets
Nr = fix(ny*nx*(K^2)/2);
D = sparse(1,1,0,Nr,1); % this array is big enough to hold all pairs 
R = sparse(1,1,0,Nr,K); % this holds the correlation matrix
W = sparse(1,1,0,Nr,1); % this holds the joint weight, sqrt(wt(k)*wt(k'))
cnt = 0;
% loop through all offset pairs
for k=1:K
    for kprime=(k+1):K
        % get overlap
        %id = find(wt(:,k)>0 & wt(:,kprime)>0); % faster as logical if don't do > test, assume all positive is OK
        id = find(wt(:,k) & wt(:,kprime)); % these are points that contribute to offset
        %for i=1:length(id) % New way is significantly faster
        if length(id)
            %cnt = cnt +1;
            cnt = cnt(end)+1:cnt(end)+length(id);
            %D(cnt,1) = (z(id(i),k)-z(id(i),kprime));
            D(cnt,1) = (z(id,k)-z(id,kprime));
            %R(cnt,k) = 1;
            R(cnt,k) = ones(size(cnt))';
            %R(cnt,kprime) = -1;
            R(cnt,kprime) = -R(cnt,k);
            %W(cnt,1) = sqrt(wt(id(i),k)*wt(id(i),kprime));
            W(cnt,1) = sqrt(wt(id,k).*wt(id,kprime));
        end
    end
    disp(sprintf('Finished weighting %d of %d sets',k,K))
end
cnt = cnt(end);

% convert to weighted space (priestly p.315)
D = D(1:cnt).*W(1:cnt);
% avoid making huge matrix with repmat
R = R(1:cnt,:);
%for k=1:K
    %error here, index exceeds dims    
    %R(1:cnt,k) = R(1:cnt,k).*(W(1:cnt)+eps);
%end
R = R.*repmat((W(1:cnt)+eps),1,K); % replaced for loop for speed

n = sum(W);

% model-model correlation (mult. by weights here)
RtR = R'*R/n;

% model-data correlation
RtD = D'*R/n;

% it is possible that there is no overlap for some data sets, remove
id = find(diag(RtR)>0);
K = length(id);
if(K==0)
    warning('no overlap found')
    return;
end
RtR = RtR(id,id);
RtD = RtD(id);
kid = kid(id);

% augment with Lagrange Mult., B
B = 0; % sum of offsets equals this
RtD = [RtD,B];
RtR = [RtR,ones(K,1)];
RtR = [RtR;ones(1,K+1)];
RtR(end) = 0;

% invert (carefully)
cnt = 0;
RtRinv = nan;
dzsolution = nan;
while(( any(isnan(dzsolution)) | any(isinf(dzsolution)) | any(isnan(RtRinv(:))) | any(isinf(RtRinv(:)))) & cnt<=10)
    p = 0;
    if(cnt>0)
        warning('attempting to stabilize RtR')
        p = eps*(2^cnt);
    end
    % try direct inversion after stabilizing (cnt>0)
    RtRtmp = RtR;
    RtRtmp(1:K,1:K) = RtRtmp(1:K,1:K) + p*eye(K);
    RtRinv = inv(RtRtmp);

    % solve
    dzsolution = RtRinv*RtD';
    if(any(isnan(dzsolution)) | any(isinf(dzsolution)))
        % try direct solution
        warning('matrix solution failed, solving directly')
        RtRinv = ones(size(RtRinv));
        dzsolution = RtRtmp\RtD';
    end
    cnt = cnt + 1;
end

% toss Lagrange Mult
offset(kid) = dzsolution(1:K);

% get error
% model residuals
msz = (D'*D)/n;
msr = msz - (dzsolution')*RtR*(dzsolution);
sk = 1-msr/msz;

% and perhaps we want all variance estimates
offseterr = offseterr + sqrt(msz); % default error
errsolution = (sqrt(diag(RtRinv)*msr/(n-M)));
if(any(imag(errsolution)))
    % failed to compute error
    errsolution = 0*errsolution + sqrt(msz);
end
offseterr(kid) = errsolution(1:K);
