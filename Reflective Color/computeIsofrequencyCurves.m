% % % % %  Hyperbolic dispersion phase diagram for Au-MaPbI3  script
% % % % %   m is filling fraction of the metal and Perovskite medium
% % % % %  By Mohammad Pourmand
clear ;
close all;
%  Initialization
    
% Gamma=1/tau;                          % The scattering rate per 1/s
T=300;                                  % The absolute temperature per Kelvin
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

% 

fname='<\dir\*>\MATLAB\6-TiNABsorbers\ellipsometry\Amorphous_Sb2S3.txt';
[n_SbS_a,k_SbS_a]=Extract_n_k(fname,nm,lambda_margin);

fname='<\dir\*>\MATLAB\6-TiNABsorbers\ellipsometry\Crystalline_Sb2S3.txt';
[n_SbS_c,k_SbS_c]=Extract_n_k(fname,nm,lambda_margin);

% % % %   Calculating the overall permittivity if Sb2S3
        
  epsilon_SbS_a=(n_SbS_a+j*k_SbS_a).^2;

  epsilon_SbS_c=(n_SbS_c+j*k_SbS_c).^2;
% ==========================================================================================

nm=1e-9; 
% Ratio=cat(2,linspace(00,0.96,21),...
%             linspace(0.9625,0.976,11),...
%             linspace(0.97801,0.98201,19),...
%             linspace(0.9825,0.986,5),...
%             linspace(0.9865,1,5));
Ratio=linspace(-10,10,21);
pks=zeros(1,length(Ratio));
N=2;
dS=100*um;
dm=10*nm;
% d_PCM=25*nm;
theta=90;
dB=50e-9;
m=0.5;

rho=0.0;

% %      l_ini:(l_final-l_ini)/Step:l_final  
% % i=1;
lambda_margin=400*nm;
for jj=1:length(Ratio)
    rho=Ratio(1,jj);  
        i=1;
%     for i=1:length(lambda_margin)
        lambda= lambda_margin(i);
        freq=c./lambda;
        omega=2*pi*freq;
%         KbT=T*Kb;
% 
% 
% %     % % % %     Ag layer refractive index calculation  %%%%%%%%%%%
        epsilon_m=epsilonAg_r-1*omegaAg_p.^2/(omega^2+1i*omega/tauAg);
%         epsilonB=epsilonAg;
%         dB=50e-9;
% %     % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

%     % % % % % % %  Phase Change Material Layer optical n , k     calculations
        d_PCM=(1/m-1)*dm;
% % % % %         EMT approximation
        Beta=rho*(epsilon_SbS_c(i)-1)./(epsilon_SbS_c(i)+2)...
                -(rho-1).*(epsilon_SbS_a(i)-1)./(epsilon_SbS_a(i)+2);
        epsilon_melt=(2*Beta+1)./(1-Beta);
%         epsilon_PCM=epsilon_SbS_a(i);
        epsilon_PCM=epsilon_melt;
%         epsilon_m=epsilon_vn(i);
     fm=dm/(dm+d_PCM);
       fd=d_PCM/(dm+d_PCM);
        epsilon_t=fm*epsilon_m+fd*epsilon_PCM;
            epsilon_z=1/...
                (fm./epsilon_m+fd./epsilon_PCM);
  

 
%         %% 
        %%%% incident, reflected and transmitted angle
        thetai=theta*pi/180;
        thetar=thetai;
        thetat=thetai;
        ki=omega/c;
        kix=ki*sin(thetai);
        kiM=sqrt(ki*ki*epsilon_m-kix*kix);
 
%         kiP=sqrt(ki*ki*epsilon_PCM-kix*kix);  
%         kiPc=sqrt(ki*ki*epsilon_melt-kix*kix);
%         kiPa=sqrt(ki*ki*epsilon_SbS_a(i)-kix*kix);
%         kiVN=sqrt(ki*ki*epsilon_vn-kix*kix);
 % % % % % % % % %          Wavenumber's for EMT approximated Layer
            kiTE=sqrt(ki*ki*epsilon_t-kix*kix);
                kiTM=sqrt(ki*ki*epsilon_t-kix*kix*epsilon_t/epsilon_z);        
        


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
        epsilon_z_r(jj,i)=real(epsilon_z);
        epsilon_z_i(jj,i)=imag(epsilon_z);
        epsilon_t_r(jj,i)=real(epsilon_t);
        epsilon_t_i(jj,i)=imag(epsilon_t);

        k_z_TM(jj,i)=(kiTM./ki);
        k_x_TM(jj,i)=(ki);
        k_z_TE(jj,i)=(kiTE./ki);

%     end

end
% %%
%   ind=find(lambda_margin/um>0.5 & lambda_margin/um<.7);
%     R_p=R_TM(:,ind);
%     lambda=lambda_margin(ind)./um;
% for jj= 1:size(Ratio,2)
% 
% 
%         [~,x_value,~]=findpeaks(1-R_p(jj,:), lambda,'MinPeakHeight',0.5,...
%             'WidthReference','halfheight');
%         
%         
%         if  ~isempty(x_value)
%             %       [pks(jj),I]=max(pk);
%             ind=find(lambda>x_value-0.01 & lambda<x_value+0.01);
%             %             [pk,xx_value,width]=findpeaks(-R_TM(jj,ind), lambda_margin/um(ind),'Annotate','extents','WidthReference','halfheight');
%             [pk,xx_value,width]=findpeaks(-R_p(jj,ind), lambda(ind));
%             
%             if ~isempty(pk)
%                 [pks(jj),I]=max(pk);
%                 widths(jj)=width(I);
%                 x_values(jj)=xx_value(I);
%             end
%         end
% end
 %%
figure;
plot( k_x_TM,k_z_TM);
%       figure;
% hold on;
% plot(Ratio, x_values,':b');
%     for kk=1: 2:size(R_TM,1)
% plot(Ratio(1,kk), x_values(1,kk),'o','Color',c(kk,:),...
%           'MarkerEdgeColor',c(kk,:),'MarkerFaceColor',c(kk,:),'LineWidth',2);
%     end
% 
%       ylabel('Wavelength (\mum)','FontSize',12,'FontName','TimesNewRoman');
%     xlabel('Crystallinity ratio {\it m}','FontSize',18,'FontName','TimesNewRoman');
% %       xlim([0 90]);
%       yyaxis right;
% % 
%    for kk=1: 2:size(R_TM,1)
%       plot(Ratio(1,kk),10*log10(pks(1,kk)),'s-',...
%         'MarkerSize',9,'Color','k','LineWidth',1);
%    end
% plot(Ratio,10*log10(-pks),'-.k');
% hold off;
% grid on;
% box on;
%       ylabel('Reflection (dB)','FontSize',12,'FontName','TimesNewRoman');
%     set(gca,'fontsize',18,'FontWeight','normal','LineWidth',2,'XMinorTick','on','YMinorTick','on');
% % set(gcf, 'color', 'none');   
% % set(gca, 'color', 'none');
% % exportgraphics(gca,'fig8_inset_n.png','BackgroundColor','none');
%%

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
% 
% function PlotSigma(freq,p,realsg,imagsg)
% %  Plot Real and Imaginary parts of effective permittivity against
% %  wavelength
% h=figure(4);
% set(h,'Position',[100 180 560 420]);
% ax1=subplot(2,1,1);
% plot(freq/p,realsg,':b','LineWidth',3);
% % ylim([-10 10]);
% legend('\epsilon_r');xlabel('Frequency (THz)','FontSize',12,'FontName','TimesNewRoman');
% ylabel('Real \epsilon_r','FontSize',12,'FontName','TimesNewRoman');
% grid on;
% 
% % % % % % % % % % % Imainary Part
% ax2=subplot(2,1,2);
% plot(freq/p,imagsg,':r','LineWidth',3);
% legend('\epsilon_r');xlabel('Frequency (THz)','FontSize',12,'FontName','TimesNewRoman');
% % ylim([-10 10]);
% ylabel('Imaginary \epsilon_r','FontSize',12,'FontName','TimesNewRoman');
% grid on;
% % % %
% end
% 
% % function [Beta,K]=PlotBlochwavenumberTE(epsilon,d,freq)
% % sigma=sg_g(0.5,freq);
% % c=3e8;
% % k0=2*pi*freq/c;
% % kix=-150*k0:150*k0;
% % Beta=sqrt(k0*k0*epsilon-kix*kix);
% % 
% % K=(cos(Beta*d)-(2i*pi*sigma*Beta/k0/c/epsilon).*sin(Beta*d))/d;
% % Beta=Beta/k0;
% % K=K/k0;
% % plot(Beta,K);
% % end