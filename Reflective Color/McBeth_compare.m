sum=0;
for i=1:length(mcBeth_xyz)
    for j=1:length(xyz)
            if sqrt((mcBeth_xyz(i,1)-xyz(j,1))^2+(mcBeth_xyz(i,2)-xyz(j,2))^2)<10e-3
                sum=sum+1;
                fprintf('%d , %d\n',i,j);
            end
    end
end