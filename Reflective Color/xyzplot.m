xyz_primaries = rgb2xyz([1 0 0; 0 1 0; 0 0 1]);

xyzMag = sum(xyz_primaries,2);
x_primary = xyz_primaries(:,1)./xyzMag;
y_primary = xyz_primaries(:,2)./xyzMag;



wp = whitepoint('D65');
wpMag = sum(wp,2);
x_whitepoint = wp(:,1)./wpMag;
y_whitepoint = wp(:,2)./wpMag;

figure;
plotChromaticity();
%plot  CIE 1931 diagram

hold on
scatter(x_whitepoint,y_whitepoint,36,'black');
plot([x_primary; x_primary],[y_primary; y_primary],'k');

load('r_xyz.mat');
xyzMag = sum(xyz,2);
x_in = xyz(:,1)./xyzMag;
y_in= xyz(:,2)./xyzMag;
% scatter(x_in(1),y_in(1),36,'s','filled','black')
scatter(x_in,y_in,128,'s','MarkerEdgeColor','w','LineWidth',2)
plot(x_in,y_in,':w','LineWidth',3)
% % % % % % % % % % % % % % % % % % % % % 
load('g_xyz.mat');
xyzMag = sum(xyz,2);
x_in = xyz(:,1)./xyzMag;
y_in= xyz(:,2)./xyzMag;
% scatter(x_in(1),y_in(1),36,'s','filled','black')
scatter(x_in,y_in,128,'^','MarkerEdgeColor','w','LineWidth',2)
plot(x_in,y_in,':w','LineWidth',3)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
load('b_xyz.mat');
xyzMag = sum(xyz,2);
x_in = xyz(:,1)./xyzMag;
y_in= xyz(:,2)./xyzMag;
% scatter(x_in(1),y_in(1),36,'s','filled','black')
scatter(x_in,y_in,128,'o','MarkerEdgeColor','w','LineWidth',2)
plot(x_in,y_in,':w','LineWidth',3)
set(gca,'fontsize',18,'FontWeight','normal','LineWidth',2,'XMinorTick','on','YMinorTick','on');

hold off