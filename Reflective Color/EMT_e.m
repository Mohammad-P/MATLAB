function [epsilon_t,epsilon_z]=EMT_e(epsilon_m,epsilon_PCM,dm,d_PCM)
     fm=dm/(dm+d_PCM);
       fd=d_PCM/(dm+d_PCM);
        epsilon_t=fm*epsilon_m+fd*epsilon_PCM;
            epsilon_z=1/...
                (fm./epsilon_m+fd./epsilon_PCM);
  end