function [n,k]= Extract_n_k(filename,lambda_margin)
fid=fopen(filename,'rb');
string = fread(fid,inf,'uint8=>char')';
M=str2num(string);
lambda=M(:,1)*1e-9;
switch size(M,2)
 case 1
fprintf(' Not enough data in text file')';
case 3
  n_raw=M(:,2); k_raw=M(:,3);
    n=interp1(lambda, smooth(n_raw), lambda_margin,'makima','extrap');
     n(find(n==-inf))=-50;  % check if there are -infinity data points
        k=interp1(lambda, smooth(k_raw), lambda_margin,'makima','extrap');
          k(find(k==-inf))=-50;  % check if there are -infinity data points
case 2
  n_raw=M(:,2);
    n=interp1(lambda, smooth(n_raw), lambda_margin,'makima','extrap');
     n(find(n==-inf))=-50;  % check if there are -infinity data points
k=NaN;

end

   
end

