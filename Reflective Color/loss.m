function [loss]= loss(Ab_TM,lambda_margin,lambda_center)
FWHM=50e-9;
sigma=FWHM/2/sqrt(2*log(2));
eta=1;
% N=size(lambda_margin,2);
% % lambda_center=1310e-9;
% Loss=10/N*sqrt(...
%     sum(abs(...
%         Ab_TM-Gaussian(lambda_margin).^2)))+...
%             sqrt(sum(...
%                     abs(T_TM).^2));
ind_low=(lambda_margin>lambda_center-10*FWHM & lambda_margin<lambda_center-FWHM/2);
Fl=sum(Ab_TM(ind_low).*Gaussian(lambda_margin(ind_low)))./sum(Gaussian(lambda_margin(ind_low)));
    ind_center=(lambda_margin>lambda_center-FWHM/2 & lambda_margin<lambda_center+FWHM/2);
    Fc=sum(Ab_TM(ind_center).*Gaussian(lambda_margin(ind_center)))/sum(Gaussian(lambda_margin(ind_center)));
        ind_up=(lambda_margin>lambda_center+FWHM/2 & lambda_margin<lambda_center+10*FWHM);
        Fu=sum(Ab_TM(ind_up).*Gaussian(lambda_margin(ind_up)))./sum(Gaussian(lambda_margin(ind_up)));
 F=2*Fc-0.3*Fl-0.3*Fu;
loss=1-F/1.4;

        function [y]=Gaussian(lambda)
        y=eta*...
                            exp(-0.5*(...
                                        (lambda-lambda_center)./sigma).^2);
        end
end
