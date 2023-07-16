% %  This simple code computes the permeabilty and permittivity tensors
% % for the titled medium by phi degree
function [t_epsilon,t_mu]=tensor_per(phi,epsilon_z,epsilon_t,mu_z,mu_t)
    t_epsilon=zeros(3);
        t_mu=zeros(3);
                t_epsilon(2,2)=epsilon_z;
                    t_mu(2,2)=mu_z;
                t_epsilon(1,1)=epsilon_z*cos(phi)^2+epsilon_t*sin(phi)^2;
            t_epsilon(3,3)=epsilon_z*sin(phi)^2+epsilon_t*cos(phi)^2;
            t_epsilon(1,3)=(epsilon_t-epsilon_z)*sin(phi)*cos(phi);
                t_epsilon(3,1)=t_epsilon(1,3);
                    t_mu(1,1)=mu_z*cos(phi)^2+mu_t*sin(phi)^2;
                t_mu(3,3)=mu_z*sin(phi)^2+mu_t*cos(phi)^2;
            t_mu(1,3)=(mu_t-mu_z)*sin(phi)*cos(phi);
    t_mu(3,1)=t_mu(1,3);
end

