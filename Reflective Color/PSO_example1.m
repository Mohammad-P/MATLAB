% PSO Params
m=3;
n=5;
Wmax=0.9;
Wmin=0.4;
c1=2;
c2=2;
% Population Initialization
x0=zeros(n,m);
Lb=[0,0,0];
Ub=[10,10,10];
Maxt=50;
for i=1:n
    for j=1:m
        x0(i,j)=round(Lb(j)+rand()*(...
                            Ub(j)-Lb(j)));
    end
end