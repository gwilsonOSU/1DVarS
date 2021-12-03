%
% AD symmetry test
%
clear

% define input variables
load ~wilsongr/work/unfunded_projects/sedimentTransport1D_TLAD/waveModel_jtech2018/example_inputOutput/assim_1dh_output_oct.mat
waves=posterior;
x=waves.x;
nx=length(x);
h=flipud(filtfilt(ones(5,1)/5,1,waves.h));

% background NL model run.  Need to run the full model b/c ddx_upwind.m
% requires q as input
dt=60*60;
nsubsteps=1;
H0=1;
theta0=0;
omega=2*pi/10;
ka_drag=0.01;
tau_wind=zeros(nx,2);
detady=zeros(nx,1);
dgamma=zeros(nx,1);
d50=ones(nx,1)*300e-6;
d90=ones(nx,1)*500e-6;
params.fv=0.1;
params.n=1.2;
params.m=11;
params.xi=1;  % ??? tuning parameter, O(1) according to Kranenburg (2013)
params.alpha=8.2;  % come in eqn 27-28, not the same as eqn 19
sedmodel='vanderA';
[Hrms,vbar,theta,kabs,q,hpout,workspc] = ...
    hydroSedModel(x,h,H0,theta0,omega,ka_drag,tau_wind,detady,dgamma,...
                  d50,d90,params,sedmodel,dt,nsubsteps);
dqdx=ddx_upwind(x,q,h);

% apply TL and ADJ models for n instances of random forcing/perturbations F
eps = 0.01;
n=5;
F = eps*rand(nx,n);  % 1st dim is number of tl input parameters
clear g
for i=1:n
  disp(['iter ' num2str(i) ' of ' num2str(n)])

  % TL model: TL*F
  tl_q=F([1:nx],i);
  tl_dqdx=tl_ddx_upwind(tl_q,x,q,h);

  % AD model: g=AD*(TL*F)
  ad_q=ad_ddx_upwind(tl_dqdx,x,q,h);

  % create output vector
  g([1:nx],i) = ad_q;

end

% test whether result makes sense
[nparam,nens]=size(g);
A = F(1:nparam,:)'*g(1:nparam,:);

(A-A')/sum(diag(A))  % should be zeros to within roundoff

%pcolor(A)
format short g
eig(A)  % should be +'ve
