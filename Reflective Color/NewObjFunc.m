function [o]=NewObjFunc(x,TM)
[nPop, nRows]=size(x);
o=zeros(1,nPop);
    for i=1:nPop
                o(i)=Newfitfunc(TM,x(i,1),x(i,2),x(i,3));
    end
end