% % % % %  Hyperbolic dispersion phase diagram for Au-MaPbI3  script
% % % % %   m is filling fraction of the metal and Perovskite medium
% % % % %  By Mohammad Pourmand
% clear ;
% close all;
% clc;
global phi
global c
global mu0
global epsilon0
tic;
%  Define frequency range
um=1e-6;
nm=1e-9;
p=1e12;
l_ini=1.2*um;l_final=1.8*um; %frequency range
Step=10000;
lambda_margin=l_ini:(l_final-l_ini)/Step:l_final;
lambda_c=1.55*um;
j=sqrt(-1);

%% Input characterisitcs of Layers
% epsilonA=2.04;                     % relative permittivity of first layer Teflon= 2.04
% epsilonB=11.9;                     % relative permittivity of second layer Si=11.9
% epsilonD=epsilonA;%2.2;             % relative permittivity of 3rd layer
% epsilonS=2.04;                     % relative permittivity of Substrate

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
%%
fname='<\dir\*>\MATLAB\6-TiNABsorbers\ellipsometry\Amorphous_Sb2S3.txt';%Sb2S3
[n_SbS_a,k_SbS_a]=Extract_n_k(fname,nm,lambda_margin);

fname='<\dir\*>\MATLAB\6-TiNABsorbers\ellipsometry\Crystalline_Sb2S3.txt';
[n_SbS_c,k_SbS_c]=Extract_n_k(fname,nm,lambda_margin);

% % % %   Calculating the overall permittivity of Sb2S3
        
  epsilon_SbS_a=(n_SbS_a+j*k_SbS_a).^2;

  epsilon_SbS_c=(n_SbS_c+j*k_SbS_c).^2;
% %
% ==========================================================================================
fname='<\dir\*>\MATLAB\6-TiNABsorbers\ellipsometry\Ag_Christy.txt';
[n_Ag,k_Ag]=Extract_n_k(fname,nm,lambda_margin);

% fname='<\dir\*>\MATLAB\Refractive Indices\Ag_20nm_k.txt';
% [k_Ag,~]=Extract_n_k(fname,lambda_margin);

% % % %   Calculating the overall permittivity of Silver
    epsilon_Ag=(n_Ag+j*k_Ag).^2;
% ==========================================================================================

nm=1e-9; 
% Ratio=linspace(0,1,11);
Ratio=[1];

Mt=eye(2,2);

param='\rho='+string(Ratio);
zci = @(v) find(diff(sign(v)));

tm=[0.12 0.13 0.13 0 0 0 1 0.9 0.8];
dz=[50 60 60 55 60 36 55 50 50]*nm;
% dz=50*nm*ones(size(tm));
%%% Initializing the TMM algorithim

parameters.IncidenceAngle.theta=0;
parameters.IncidenceAngle.phi=0;

parameters.HMM.thicknesses=dz;
parameters.HMM.fillingfactors=tm;

parameters.capping.epsilon=epsilonS;
parameters.capping.thickness=70*nm;

parameters.metal.epsilon=epsilon_Ag;
parameters.metal.thickness=50*nm;


parameters.dielectric.thickness=70*nm;

parameters.substrate.epsilon=epsilonS;
parameters.substrate.thickness=100*um;


for jj=1:length(Ratio)
    
    rho=Ratio(1,jj);  
  epsilon_melt=ComputeEpsilonMelt(epsilon_SbS_a,epsilon_SbS_c,rho);

parameters.dielectric.epsilon=epsilon_melt;

[TM,TE]= TMM_HMM(lambda_margin, parameters);
       Results.TM.R(jj,:)=TM.R;
            Results.TM.T(jj,:)=TM.T;
              Results.TM.Ab(jj,:)=TM.Ab;
        Results.TE.R(jj,:)=TE.R;
          Results.TE.T(jj,:)=TE.T;
             Results.TE.Ab(jj,:)=TE.Ab;
 Results.TM.Lambda=lambda_margin;
    % Loss(jj,1)=fitfunc(Ab_TM(jj,:),T_TM(jj,:),lambda_margin,1.5501e-6);

end


     

%%
% plotRTA(Results.TM,param);
toc;
