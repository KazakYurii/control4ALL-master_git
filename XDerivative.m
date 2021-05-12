function F = XDerivative(y,dx)
lengths=size(y);
F=zeros(lengths(1), lengths(2));
for I=1:lengths(1)
    for J=1:lengths(2)
        if I==1
            F(I,J)=(-y(I+2,J)+4*y(I+1,J)-3*y(I,J))/(2*dx);
        else
            if I==lengths(1)
                F(I,J)=(3*y(I,J)-4*y(I-1,J)+y(I-2,J))/(2*dx);
            else
                F(I,J)=(y(I+1,J)-y(I-1,J))/(2*dx);
            end
        end
    end
end
end
