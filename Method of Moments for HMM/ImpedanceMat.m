clear all;
%% Initialzation
j=sqrt(-1);
epsilon0=8.85187812813e-12;             % The vacuum permittivity constant per Fm^-1
mu0=pi*4e-7;                            % H/m Vacuum Permeability  
c=3e8;                                  % m.s^-2 speed of light in vacuum
thetai_degree=60;                       % Angle of incident in degree
freq=10e9;                              % Radiation frequency in Hz
N=10;                                    % Mesh numbers
E_inc=1;                                % Incidence electric field in V/m
%% Layer definition
mu_r0=1;
epsilon_r0=4;                           % Spatial dependence permitivitty epsilon=epsilon_r0*exp(kappa*z/t_d)
t_d=0.1;                               % layer thickness in m
kappa=1;
sigma_0=0;
%%
thetai=thetai_degree*pi/180;
omega=2*pi*freq;
k_x=omega/c*sin(thetai);
delta_z=t_d/N;
n=N;
m=N;
Anm=zeros(1,n);
Bnm=zeros(1,n);
Cnm=zeros(1,n);
Dnm=zeros(1,n);
Zmn=[];
Z_s=sqrt(mu0/epsilon0)/cos(thetai);         % For TE radiations
Z_l=Z_s;                                    % Transmission layer 
S_mat=repmat([Z_l;1],[m,1])./(Z_s+Z_l);
kron_del=@(n) double(n==0);                 % kroniker delta function
%% Integral calculations based on Linearly variations
Impd_z_TE=j*omega*mu0*mu_r0;
Admit_z_TE= @(z) sigma_0+j*omega*epsilon0*epsilon_r0*exp(kappa*z/t_d)+k_x*k_x/Impd_z_TE;
Admit_int_TE=@(lo,u,d_z) ((1-0.5*(u-lo)*Admit_z_TE(lo*d_z))+0.5*(u-lo)*Admit_z_TE((lo+1)*d_z))*(u-lo)*d_z;
Impd_int_TE=@(l,u,d_z) ((1-0.5*(u-l)*Impd_z_TE)+0.5*(u-l)*Impd_z_TE)*(u-l)*d_z;
%%
for ii=1:m
    for jj=1:n
        z_m=(ii-0.5)*delta_z;
        %% Impedance matrix elements
        Anm(1,jj)=kron_del(jj-ii)+Z_s*Z_l./(Z_s+Z_l)* Admit_int_TE(jj,jj-1,delta_z);
        Cnm(1,jj)=Admit_int_TE(jj-1+heaviside(ii-jj),jj-1,delta_z)...
                        -Z_l./(Z_s+Z_l)* Admit_int_TE(jj,jj-1,delta_z);
        Bnm(1,jj)=Impd_int_TE(jj-1+heaviside(ii-jj),jj-1,delta_z)-Z_s./(Z_s+Z_l)* Impd_int_TE(jj,jj-1,delta_z);
        Dnm(1,jj)=kron_del(jj-ii)+1./(Z_s+Z_l)* Impd_int_TE(jj,jj-1,delta_z);
    end
    Zmn=cat(1,Zmn,[Anm, Bnm;Cnm, Dnm]);
end
Fields=Zmn\S_mat;
E_field(1:m,:)=Fields(1:m,:);
H_field(1:m,:)=Fields(m+1:end,:);
figure;
plot(1:n,abs(E_field));
