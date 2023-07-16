

% input=Data001;
input=interp1(Data002(:,1), smooth(Data002(:,2)),380:780 ,'makima',0.4345);
figure;plot(Data002(:,1),Data002(:,2),'o',380:780,input,':');
input1=[rang(1,:)' input(1,:)'];
TS1=gettristimulus2degn(input1,rang');    
rgb1=getxyz(TS1)*255;
figure;rectangle('Position',[0.4 0.8 0.05 0.175],'FaceColor',rgb1/255,'EdgeColor',rgb1/255,'LineWidth',3);