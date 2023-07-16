%% This script is used to calculate color from spectrum
%It can load spectrum from an excel, the format should as follow:
%column 1: wavelength data(unit: nm/m)
%column 2-n:reflectance/transmittance spectrum
%The default light source is D65,one can change light source in
%spectrum2XYZ function
%Author:Tingbiao Guo, Zhejiang University
clear all;
figure(1);
%this script can calculate multiple spectra at once
%s means the first column need to be calculated
%e menas the last column need to be calculated
s=2;    %start column
e=3;   %end column
name='examples.xlsx';
plotChromaticity();
%plot  CIE 1931 diagram
hold on;
% coor=spectrum2XYZ(s,e,name);
x=coor(:,1);
y=coor(:,2);
sz=30;
scatter(x,y,sz,'r','filled');
%plot x y coordinate on CIE1931 diagram
%% plot color patch
z=1-x-y;
%change xyz to rgb value
%convert to uint8 type
rgb=xyz2rgb([x,y,z],'OutputType','uint8');
%convert to 0-1 type
rgb=double(rgb)/255;
%%plot color patch
figure(2);
% for i=1:(e-s+1)
rectangle('Position',[5*(i-1),0,5,5],'Curvature', [0 0], 'FaceColor',rgb(i,:));
% end
% hold on;
