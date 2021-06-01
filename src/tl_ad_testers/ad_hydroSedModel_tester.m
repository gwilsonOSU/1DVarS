%
% AD symmetry test
%
clear

load ../../waveModel_jtech2018/example_inputOutput/assim_1dh_output_oct.mat
waves=posterior;
physicalConstants;

% select dubarbier or vanderA qtrans model
% sedmodel='dubarbier';
sedmodel='vanderA';
% sedmodel='soulsbyVanRijn';

if(strcmp(sedmodel,'dubarbier'))
  params.Cw=0.00483 ;
  params.Cc=0.02002 ;
  params.Cf=0.01173 ;
  params.Ka=0.631e-4;
elseif(strcmp(sedmodel,'vanderA'))
  params.n=1.2;
  params.m=11;
  params.xi=1;  % ??? tuning parameter, O(1) according to Kranenburg (2013)
  params.alpha=8.2;  % come in eqn 27-28, not the same as eqn 19
elseif(strcmp(sedmodel,'soulsbyVanRijn'))
  params.alphab=1.6;
  params.facua =1.0;
else
  error(['invalid sedmodel=' sedmodel])
end

% define input variables
params.dt=60*60;  % one hour time step
params.fv=0.1;  % breaking-induced eddy viscosity calibration parameter
x=waves.x;
h=flipud(filtfilt(ones(5,1)/5,1,waves.h));
H0=waves.H0;
theta0=waves.theta0;
omega=waves.sigma;
ka_drag=waves.ka_drag;
windW=zeros(length(x),2);
detady=zeros(length(x),1);
d50=180e-6;
d90=400e-6;
nx=length(h);

% background NL model run
[Hrms,vbar,theta,kabs,Q,hp,bkgd]=hydroSedModel(x,h,H0,theta0,omega,ka_drag,...
                                               windW,detady,d50,d90,params,sedmodel);

% apply TL and ADJ models for n instances of random forcing/perturbations F
eps = 0.01;
n=5;
F = eps*rand(4*nx+11,n);  % 1st dim is number of tl input parameters
clear g
for i=1:n
  disp(['iter ' num2str(i) ' of ' num2str(n)])

  % TL model: TL*F
  tl_h      =F(0*nx+[1:nx],i);
  tl_H0     =F(1*nx+1,i);
  tl_theta0 =F(1*nx+2,i);
  tl_omega  =F(1*nx+3,i);
  tl_ka_drag=F(1*nx+4,i);
  tl_windW=reshape(F(1*nx+4+[1:(2*nx)],i),[nx 2]);
  tl_detady=reshape(F(3*nx+4+[1:nx],i),[nx 1]);
  tl_d50    =F(4*nx+5,i);
  tl_d90    =F(4*nx+6,i);
  tl_params.fv   =F(4*nx+7,i);
  if(strcmp(sedmodel,'dubarbier'))
    tl_params.Cw   =F(4*nx+8,i);
    tl_params.Cc   =F(4*nx+9,i);
    tl_params.Cf   =F(4*nx+10,i);
    tl_params.Ka   =F(4*nx+11,i);
  elseif(strcmp(sedmodel,'vanderA'))
    tl_params.n    =F(4*nx+8,i);
    tl_params.m    =F(4*nx+9,i);
    tl_params.xi   =F(4*nx+10,i);
    tl_params.alpha=F(4*nx+11,i);
  elseif(strcmp(sedmodel,'soulsbyVanRijn'))
    tl_params.alphab=F(4*nx+8,i);
    tl_params.facua =F(4*nx+9,i);
  end
  [tl_Hrms,tl_vbar,tl_theta,tl_kabs,tl_Q,tl_hp] = ...
      tl_hydroSedModel(tl_h,tl_H0,tl_theta0,tl_omega,tl_ka_drag,...
                       tl_windW,tl_detady,tl_d50,tl_d90,tl_params,bkgd);

  % AD model: g=AD*(TL*F)
  [ad_h,ad_H0,ad_theta0,ad_omega,ad_ka_drag,...
   ad_windW,ad_detady,ad_d50,ad_d90,ad_params] = ...
      ad_hydroSedModel(tl_Hrms,tl_vbar,tl_theta,tl_kabs,tl_Q,tl_hp,bkgd);

  % create output vector
  g(0*nx+[1:nx],i)      =ad_h(:)     ;
  g(1*nx+1,i)           =ad_H0       ;
  g(1*nx+2,i)           =ad_theta0   ;
  g(1*nx+3,i)           =ad_omega    ;
  g(1*nx+4,i)           =ad_ka_drag  ;
  g(1*nx+4+[1:(2*nx)],i)=ad_windW(:) ;
  g(3*nx+4+[1:(nx)],i)=ad_detady(:);
  g(4*nx+5,i)           =ad_d50      ;
  g(4*nx+6,i)           =ad_d90      ;
  g(4*nx+7,i)           =ad_params.fv;
  if(strcmp(sedmodel,'dubarbier'))
    g(4*nx+8,i) =ad_params.Cw;
    g(4*nx+9,i) =ad_params.Cc;
    g(4*nx+10,i)=ad_params.Cf;
    g(4*nx+11,i)=ad_params.Ka;
  elseif(strcmp(sedmodel,'vanderA'))
    g(4*nx+ 8,i)=ad_params.n    ;
    g(4*nx+ 9,i)=ad_params.m    ;
    g(4*nx+10,i)=ad_params.xi   ;
    g(4*nx+11,i)=ad_params.alpha;
  elseif(strcmp(sedmodel,'soulsbyVanRijn'))
    g(4*nx+8,i)=ad_params.alphab;
    g(4*nx+9,i)=ad_params.facua ;
  end

end

% test whether result makes sense
[nparam,nens]=size(g);
A = F(1:nparam,:)'*g(1:nparam,:);

(A-A')/sum(diag(A))  % should be zeros to within roundoff

%pcolor(A)
format short g
eig(A)  % should be +'ve
