function [Loss,F,Fu,Fc,Fl]= Costfunc(R_TM,T_TM,lambda_margin,lambda_center)
FWHM=50e-9;
sigma=FWHM/2/sqrt(2*log(2));
eta=0.9;
N=size(lambda_margin,2);
% lambda_center=1310e-9;
Loss=100/N*sqrt(...
    sum(abs(...
        R_TM-Gaussian(lambda_margin).^2)))+...
            sqrt(sum(...
                    abs(T_TM).^2));
% ind_low=(lambda_margin>lambda_center-10*FWHM & lambda_margin<lambda_center-FWHM/2);
ind_low=(lambda_margin<lambda_center-FWHM/2);
nom_low=sum(R_TM(ind_low).*Gaussian(lambda_margin(ind_low)));
denom_low=sum(Gaussian(lambda_margin(ind_low)));

Fl=nom_low./denom_low;

ind_center=(lambda_margin>lambda_center-FWHM/2 & lambda_margin<lambda_center+FWHM/2);
nom_center=sum(R_TM(ind_center).*Gaussian(lambda_margin(ind_center)));
denom_center=sum(Gaussian(lambda_margin(ind_center)));

Fc=nom_center./denom_center;

% ind_up=(lambda_margin>lambda_center+FWHM/2 & lambda_margin<lambda_center+10*FWHM);
ind_up=(lambda_margin>lambda_center+FWHM/2);
nom_up=R_TM(ind_up)*Gaussian(lambda_margin(ind_up))';
denom_up=sum(Gaussian(lambda_margin(ind_up)));

Fu=nom_up./denom_up;
F=1-(Fc-Fl-Fu);

function [y]=Gaussian(lambda)
    y=eta*...
         exp(-0.5*(...
                      (lambda-lambda_center)./sigma).^2);
end
end
