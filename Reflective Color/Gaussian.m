function [y]=Gaussian(lambda,lambda_center,FWHM,eta)
% FWHM=50e-9;
sigma=FWHM/2/sqrt(2*log(2));
% eta=1;
        y=eta*...
                  exp(-0.5*(...
                               (lambda-lambda_center)./sigma).^2);
end