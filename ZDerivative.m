function F = ZDerivative(array,dz)
length=size(array);
F=zeros(length(1),length(2));
for i=1:length(1)
    for j=1:length(2)
        if j==1
            F(i,j)=(-array(i,j+2)+4*array(i,j+1)-3*array(i,j))/(2*dz);
        else
            if j==length(2)
                F(i,j)=(3*array(i,j)-4*array(i,j-1)+array(i,j-2))/(2*dz);
            else
                F(i,j)=(array(i,j+1)-array(i,j-1))/(2*dz);
            end
            
        end
    end
end

end