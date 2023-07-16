% This function calulates an array of [r ,g ,b] from the reflection spectrum based on D65 illuminant.
% the input must be in form of m*n where the first column is wavelength/nm and the rest
% of columns include the reflection coefficient ( 0 to 1)
function [xyz,rgb] = getrgb(input,range)
load cie1931xyz1nm;
load 'illD65';
%CIE 1931 2-deg CMFs (CIE, 1932)
%CIE. (1926). Commission Internationale de l'Eclairage Proceedings, 1924. Cambridge: Cambridge University Press.
%CIE. (1932). Commission Internationale de l'Eclairage Proceedings, 1931. Cambridge: Cambridge University Press.
%Guild, J. (1931). The colorimetric properties of the spectrum. Philosophical Transactions of the Royal Society of London, A230, 149-187.

%interpolate based on input and save variable
cie1931xyz1nm = horzcat(range,interp1(cie1931xyz1nm (:,1),cie1931xyz1nm (:,[2 3 4]),range,'linear'));

%find starting point based on range input
startval = find(cie1931xyz1nm(:,1) == min(range));
endval = find(cie1931xyz1nm(:,1) == max(range));

inputstartval = find(input(:,1) == min(range));
inputendval = find(input(:,1) == max(range));

xyzbar=cie1931xyz1nm (:,[1 4 3 2]);

I_D=interp1(ill_D65(:,1),ill_D65(:,2),range,'linear');
N= sum(cie1931xyz1nm(startval:endval,3).*I_D);

XYZ = zeros(size(input,2)-1,3);
for j=1:size(input,2)-1
    for i=1:size(cie1931xyz1nm,2)-1
        XYZ(j,i) = 100./N.*sum(cie1931xyz1nm(startval:endval,i+1).*(input(inputstartval:inputendval,j+1).*I_D));
        i = i + 1;
    end
    j = j + 1;
end
% xyz=getxyz(XYZ);
xyzMag = sum(XYZ,2);
x = XYZ(:,1)./xyzMag;
y= XYZ(:,2)./xyzMag;
z= XYZ(:,3)./xyzMag;
xyz=[x,y,z];
rgb1=xyz2rgb(xyz,'OutputType','uint8');
%convert to 0-1 type
rgb=double(rgb1)/255;
