% % % % %  Hyperbolic metamaterial based absorber with PSO optimization
% % % % %   m is filling fraction of the metal and PCM medium
% % % % %  By Mohammad Pourmand
clear ;
close all;
clc;
% Initializing the paramters
global phi
global c
global mu0
global epsilon0
global parameters
tic;
%% define the global parameters
T=300;                                  % The absolute temperature per Kelvin
e=1.60217662e-19;                       % The electron Charge per Columbs
h_bar=1.054571817e-34;                  % The Plank constant equlas to h/2pi 6.62607004e-34/2pi, per m^2 Kg S^-1(6.582119569e-19)
Kb=1.3806452e-23;                       % The Boltzmann's constant per m^2 Kg S^-2 K^-1
epsilon0=8.85187812813e-12;              % The vacuum permittivity constant per Fm^-1
mu0=pi*4e-7;                            % H/m Vacuum Permeability  
c=3e8;                                  % m.s^-2 speed of light in vacuum
%%
%  Define frequency range
um=1e-6;
nm=1e-9;
p=1e12;
l_ini=0.2*um;l_final=1.5*um; %frequency range
Step=10000;
lambda_margin=l_ini:(l_final-l_ini)/Step:l_final;
lambda_c=1.55*um;
j=sqrt(-1);
c=3e8;
epsilonAg_r=1;
omegaAg_p=1.38e16;                    % rad/s
tauAg=16.8e-15;                          % rad/s, tau=16.8 fs


muA=1;                                 % relative permeability of first layer
muB=1;                                 % relative permeability of second layer
muD=1;                                  % relative permeability of 3rd layer

MA=zeros(2);
MB=zeros(2);
MG=zeros(2);
T_T=zeros(2);
Td=zeros(2);
M=zeros(2);

R_TE=zeros(1,Step+1);R_TM=zeros(1,Step+1);R_Td=zeros(1,Step+1);T_Td=zeros(1,Step+1);
T_TM=zeros(1,Step+1);T_TE=zeros(1,Step+1);Ab_TE=zeros(1,Step+1);Ab_TM=zeros(1,Step+1);
g_ei=zeros(1,Step+1);g_er=zeros(1,Step+1);K=zeros(1,Step+1);Fr=zeros(1,Step+1);
Sg=zeros(1,Step+1);Tr=zeros(1,Step+1);R=zeros(1,Step+1);Ab=zeros(1,Step+1);Tetha=zeros(1,Step+1);

%%
%     % % % % % % % % % % % % % % % % SiO2 % % % % % % % % % % % % % % % %
%     % % % % %
%   
         epsilonS=(1+0.6961663./(1-(0.0684043./lambda_margin/1e6).^2)+0.4079426./(1-(0.1162414./lambda_margin/1e6).^2)+0.8974794./(1-(9.896161./lambda_margin/1e6).^2));
% %    2) C. Z. Tan. Determination of refractive index of silica glass for
% infrared wavelengths by IR spectroscopy, J. Non-Cryst. Solids 223,
% 158-163 (1998) %    Sellmeier formula is reported in Ref. 1 for the
% 0.21-3.71 ?m wavelength range. Ref. 2 verifies the validity of the
% formula up to 6.7 ?m.
%     % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%     % % % % % % % % % % % % % % %
%% Definig the Electrical properties of the mediums
load '<\dir\*>\MATLAB\6-TiNABsorbers\ellipsometry\GeSe.mat';
lambda=n_a(:,1);
  n_a_raw=n_a(:,2);
  n_a_r=interp1(lambda*nm, smooth(n_a_raw), lambda_margin,'makima','extrap');
  n_a_r(find(n_a_r==-inf))=-50;  % check if there are -infinity data points
lambda=k_a(:,1);
  k_a_raw=k_a(:,2);
  k_a_r=interp1(lambda*nm, smooth(k_a_raw), lambda_margin,'makima','extrap');
  k_a_r(find(k_a_r==-inf))=-50;  % check if there are -infinity data points
lambda=n_c(:,1);
  n_c_raw=n_c(:,2);
  n_c_r=interp1(lambda*nm, smooth(n_c_raw), lambda_margin,'makima','extrap');
  n_c_r(find(n_c_r==-inf))=-50;  % check if there are -infinity data points
lambda=k_c(:,1);
  k_c_raw=k_c(:,2);
  k_c_r=interp1(lambda*nm, smooth(k_c_raw), lambda_margin,'makima','extrap');
  k_c_r(find(k_c_r==-inf))=-50;  % check if there are -infinity data points
% % % %   Calculating the overall permittivity of GeSe
epsilon_PCM_a=(n_a_r+j*k_a_r).^2;
epsilon_PCM_c=(n_c_r+j*k_c_r).^2;
% %
% ==========================================================================================
fname='<\dir\*>\MATLAB\6-TiNABsorbers\ellipsometry\Al_nm_n_k.txt';
[n_Al,k_Al]=Extract_n_k(fname,nm,lambda_margin);

% % % %   Calculating the overall permittivity of Silver
    epsilon_m=(n_Al+j*k_Al).^2;
% ==========================================================================================

    Mt=eye(2,2);


    tm=[0.12 0.13 0.13 0 0 0 1 0.9 0.8];
    dz=[50 60 60 55 60 36 55 50 50]*nm;

% dz=50*nm*ones(size(tm));
%%% Initializing the TMM algorithim

    parameters.IncidenceAngle.theta=0;
    parameters.IncidenceAngle.phi=0;

    parameters.HMM.thicknesses=dz;
    parameters.HMM.fillingfactors=tm;

    parameters.capping.epsilon=epsilonS;
    parameters.capping.thickness=100*nm;

    parameters.metal.epsilon=epsilon_m;
    parameters.metal.thickness=50*nm;

    parameters.dielectric.epsilon=epsilon_PCM_a;
    parameters.dielectric.thickness=50*nm;

    parameters.substrate.epsilon=epsilonS;
    parameters.substrate.thickness=100*um;

    parameters.lambda=lambda_margin;
% Problem definition
    nVars=3;
    lb=[0 0 0.3];
    ub=[0.2 0.5 0.6];
   
    options = optimoptions(@particleswarm,...
                            'SwarmSize',20,'OutputFcn','psooutfcn');
options.Display='iter';
options.CreationFcn=@pswcreation;
options.UseVectorized=true;
options.SocialAdjustmentWeight=2;
options.SelfAdjustmentWeight=2;
options.UseParallel=true;
options.MaxIterations=50;
rng default;
   fun=@(x) Objectivefunction(x);

%     rho=Ratio(1,jj);  
%   epsilon_melt=ComputeEpsilonMelt(epsilon_SbS_a,epsilon_SbS_c,rho);
% parameters.dielectric.epsilon=epsilon_melt;

   [Gbest,fval,exitflag,output]=particleswarm(fun,nVars,lb,ub,options);


%%
% plotRTA(TM,param);
toc;
