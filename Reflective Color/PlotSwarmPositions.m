% function plotSwarmmovement(X,max_iteration,noP)

noP=options.SwarmSize;
% max_iteration=2;
% noP=10;
max_iteration=size(X,1)/noP;
figure;
for ind= 0: max_iteration-1
  n_row=ind*noP+1;
    clf;
  plot(X(n_row:n_row+noP-1,2),X(n_row:n_row+noP-1,1),'o','LineWidth',3);
       grid on
set(gca,'fontsize',18,'FontWeight','normal','LineWidth',2,'XMinorTick','on','YMinorTick','on');
titl='Iteration no. '+string(ind);
title(titl,'FontSize',18,'FontName','TimesNewRoman');
    xlabel('X1','FontSize',18,'FontName','TimesNewRoman');
    ylabel('X2','FontSize',18,'FontName','TimesNewRoman');
%         ylim([0 1]);
%             xlim([0 1]);  
      pause(1); %Pause a little bit before plotting the next iteration

end