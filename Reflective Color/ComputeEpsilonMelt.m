function epsilon_melt=ComputeEpsilonMelt(epsilon_PCM_a,epsilon_PCM_c,rho)
% % % % %         EMT approximation
          Beta=rho*(epsilon_PCM_c-1)./(epsilon_PCM_c+2)...
                -(rho-1).*(epsilon_PCM_a-1)./(epsilon_PCM_a+2);
           epsilon_melt=(2*Beta+1)./(1-Beta);
end