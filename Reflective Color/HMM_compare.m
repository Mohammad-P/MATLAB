% % % % %  Hyperbolic dispersion phase diagram for Au-MaPbI3  script
% % % % %   m is filling fraction of the metal and Perovskite medium
% % % % %  By Mohammad Pourmand
clear ;
close all;
global phi
global c
global mu0
global epsilon
%  Initialization

% Gamma=1/tau;                             % The scattering rate per 1/s
T=300;                                   % The absolute temperature per Kelvin
e=1.60217662e-19;                       % The electron Charge per Columbs
h_bar=1.054571817e-34;                  % The Plank constant equlas to h/2pi 6.62607004e-34/2pi, per m^2 Kg S^-1(6.582119569e-19)
Kb=1.3806452e-23;                       % The Boltzmann's constant per m^2 Kg S^-2 K^-1
epsilon=8.85187812813e-12;              % The vacuum permittivity constant per Fm^-1
mu0=pi*4e-7;                            % H/m Vacuum Permeability  
c=3e8;                                  % m.s^-2 speed of light in vacuum
nLayer=1;
% filename=['<\dir\*>\MATLAB\pi75.mat'];
%  Define frequency range
um=1e-6;
nm=1e-9;
p=1e12;
l_ini=0.2*um;l_final=1.6*um; %frequency range
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

 %%%%%%%%%%%%           width of two layers
% dA=lambda_c/4/sqrt(epsilonA);%35e-7;    % First Dielectric thickness
% dB=lambda_c/4/sqrt(epsilonB);
% dD=1*dA;%0.61*dA;%93.6355e-7
% dS=100e-9;%dS=lambda_c/4/sqrt(epsilonS);
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

% %%
% 
% fname='<\dir\*>\MATLAB\6-TiNABsorbers\ellipsometry\Amorphous_Sb2S3.txt';%Sb2S3
% [n_SbS_a,k_SbS_a]=Extract_n_k(fname,nm,lambda_margin);
% 
% fname='<\dir\*>\MATLAB\6-TiNABsorbers\ellipsometry\Crystalline_Sb2S3.txt';
% [n_SbS_c,k_SbS_c]=Extract_n_k(fname,nm,lambda_margin);
% 
% % % % %   Calculating the overall permittivity of Sb2S3
%         
%   epsilon_SbS_a=(n_SbS_a+j*k_SbS_a).^2;
% 
%   epsilon_SbS_c=(n_SbS_c+j*k_SbS_c).^2;
% % ==========================================================================================
fname='<\dir\*>\MATLAB\6-TiNABsorbers\ellipsometry\Ag_Christy.txt';
[n_Ag,k_Ag]=Extract_n_k(fname,nm,lambda_margin);

% fname='<\dir\*>\MATLAB\Refractive Indices\Ag_20nm_k.txt';
% [k_Ag,~]=Extract_n_k(fname,nm,lambda_margin);

% % % %   Calculating the overall permittivity of Silver
    epsilon_Ag=(n_Ag+j*k_Ag).^2;
% ==========================================================================================
% % ==========================================================================================
fname='<\dir\*>\MATLAB\6-TiNABsorbers\ellipsometry\Al_nm_n_k.txt';
[n_Al,k_Al]=Extract_n_k(fname,nm,lambda_margin);

% fname='<\dir\*>\MATLAB\Refractive Indices\Ag_20nm_k.txt';
% [k_Ag,~]=Extract_n_k(fname,nm,lambda_margin);

% % % %   Calculating the overall permittivity of Silver
    epsilon_Al=(n_Al+j*k_Al).^2;
% ==========================================================================================
% % =======================================================================================
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
epsilon_GeSe_a=(n_a_r+j*k_a_r).^2;
epsilon_GeSe_c=(n_c_r+j*k_c_r).^2;
% ==========================================================================================
nm=1e-9; 
% Ratio=linspace(0,1,3);
Ratio=[0];
pks=zeros(1,length(Ratio));
N=2;
dS=100*um;
dm=50*nm;
% d_PCM=25*nm;
theta=00;
phi=00;
dB=50e-9;
m=0.34;         % Filling factor
ds1=60* nm; % Capping SiO2 layer at top
ds2=10 * nm; % SiO2 layer at the bottom  
rho=1.00;    % MelingRatio
% t=2.136;
t=1.4;
Mt=eye(2,2);
% %      l_ini:(l_final-l_ini)/Step:l_final  
% % i=1;
param='\rho='+string(Ratio);
zci = @(v) find(diff(sign(v)));

% tm=[0.11 0 1 0.81 0.78 0.96 0.31 0.47 0.54 0.53 0.54 0];
tm=[0.34 0 1 0.81 0.78 0.96 0.31 0.47 0.54 0.53 0.54 0];
dz=50*nm;
for jj=1:length(Ratio)
    rho=Ratio(1,jj);  
    
    for i=1:length(lambda_margin)
        lambda= lambda_margin(i);
        freq=c./lambda;
        omega=2*pi*freq;
        %%%% incident, reflected and transmitted angle
        thetai=theta*pi/180;
            phi=phi*pi/180;
        thetar=thetai;
        thetat=thetai;
        ki=omega/c;
%         KbT=T*Kb;
% 
%% 
%     % % % % % % % % % % % % % % % % SiO2 % % % % % % % % % % % % % % % % % % % % %
%     
         epsilonS=(1+0.6961663./(1-(0.0684043./lambda/1e6).^2)+0.4079426./(1-(0.1162414./lambda/1e6).^2)+0.8974794./(1-(9.896161./lambda/1e6).^2));
% %    2) C. Z. Tan. Determination of refractive index of silica glass for infrared wavelengths by IR spectroscopy, J. Non-Cryst. Solids 223, 158-163 (1998)
% %    Sellmeier formula is reported in Ref. 1 for the 0.21-3.71 ?m wavelength range. Ref. 2 verifies the validity of the formula up to 6.7 ?m.
%     % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% 
%         %%
% %     % % % %     Ag layer refractive index calculation  %%%%%%%%%%%
%         epsilon_m=epsilonAg_r-1*omegaAg_p.^2/(omega^2+1i*omega/tauAg);
        epsilon_m=epsilon_Al(i);
%         dB=50e-9;

% %     % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

%     % % % % % % %  Phase Change Material Layer optical n , k     calculations
        d_PCM=(1/m-1)*dm;
%%
% % % % %         EMT approximation
epsilon_PCM_a=epsilon_GeSe_a;
epsilon_PCM_c=epsilon_GeSe_c;
        Beta=rho*(epsilon_PCM_c(i)-1)./(epsilon_PCM_c(i)+2)...
                -(rho-1).*(epsilon_PCM_a(i)-1)./(epsilon_PCM_a(i)+2);
        epsilon_melt=(2*Beta+1)./(1-Beta);
        epsilon_PCM=epsilon_melt;
%   epsilon_PCM=epsilonS;
[epsilon_t,epsilon_z]=EMT_e(epsilon_m,epsilon_PCM,dm,d_PCM);
%%
        for ind=1:size(tm,2)
            dmb=tm(ind)*dz;
                [epsilon_tm,epsilon_zm]=EMT_e(epsilon_m,epsilonS,dmb,dz-dmb); 
        epsilon_z_rm(ind,i)=real(epsilon_zm);
        epsilon_z_im(ind,i)=imag(epsilon_zm);
        epsilon_t_rm(ind,i)=real(epsilon_tm);
        epsilon_t_im(ind,i)=imag(epsilon_tm);
          end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
        Lambda(i)=lambda;
        epsilon_z_r(jj,i)=real(epsilon_z);
        epsilon_z_i(jj,i)=imag(epsilon_z);
        epsilon_t_r(jj,i)=real(epsilon_t);
        epsilon_t_i(jj,i)=imag(epsilon_t);


    end
end
    
 %%
figure;    
%      plot (Lambda/um, 10*log10(R_TM(:,:)),'LineWidth',3);
   plot(Lambda/um,epsilon_z_r,'-*',Lambda/um,epsilon_t_r,'-.','LineWidth',2);
hold on;
   plot(Lambda/um,epsilon_z_rm(1,:),'-k',Lambda/um,epsilon_t_rm(1,:),'-.','LineWidth',2);

set(gca,'fontsize',18,'FontWeight','normal','LineWidth',2,'XMinorTick','on','YMinorTick','on');

   ylabel('Effective permittivity','FontSize',18,'FontName','TimesNewRoman');
%     ylim([-25 35])
    set(gca,'fontsize',18,'FontWeight','normal','LineWidth',1,'XMinorTick','on','YMinorTick','on');
    xlabel('\lambda (\mum)','FontSize',18,'FontName','TimesNewRoman');

    grid on;
    box on;
            xlim([0.2 0.8]);

%             xlim([l_ini/um l_final/um]);

    grid on;
    hold off;
    box on;

%%
 %%
figure;    
%      plot (Lambda/um, 10*log10(R_TM(:,:)),'LineWidth',3);
   plot(Lambda/nm,k_c_r,'-',k_c(:,1),k_c(:,2),'-.','LineWidth',2);
hold on;
   plot(Lambda/nm,n_c_r,'-*',n_c(:,1),n_c(:,2),'-.','LineWidth',2);
   plot(Lambda/um,k_a_r,'-k',k_a(:,1),k_a(:,2),'-.','LineWidth',2);

set(gca,'fontsize',18,'FontWeight','normal','LineWidth',2,'XMinorTick','on','YMinorTick','on');

   ylabel('Refractive index n','FontSize',18,'FontName','TimesNewRoman');
%     ylim([-25 35])
    set(gca,'fontsize',18,'FontWeight','normal','LineWidth',1,'XMinorTick','on','YMinorTick','on');
    xlabel('\lambda (\mum)','FontSize',18,'FontName','TimesNewRoman');

    grid on;
    box on;
            xlim([0.2 0.8]);

%             xlim([l_ini/um l_final/um]);

    grid on;
    hold off;
    box on;
%%
function M=Mi_TE(kiz,di,omega)
mu0=pi*4e-7; 
M(1,1)=cos(kiz*di);
M(2,2)=M(1,1);
pim=kiz/omega/mu0;
M(1,2)=-1i/pim*sin(kiz*di);
M(2,1)=-1i*pim*sin(kiz*di);
end
function M=Mi_TM(kiz,di,omega,epsiloni)
mu0=pi*4e-7; 
epsilon=8.85187812813e-12;     
M(1,1)=cos(kiz*di);
M(2,2)=M(1,1);
pim=omega*epsilon*epsiloni/kiz;
M(1,2)=-1i/pim*sin(kiz*di);
M(2,1)=-1i*pim*sin(kiz*di);
end
