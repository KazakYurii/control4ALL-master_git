function F = intFdxdz(f,dx1,dy1)
[nx, ny]=size(f);
F=0;
for i=2:2:nx-1
    for j=2:2:ny-1
        F=F+...
            (f(i-1,j-1)+f(i+1,j-1)+f(i-1,j+1)+f(i+1,j+1)+4*(f(i,j+1)+f(i,j-1)+f(i-1,j)+f(i+1,j))+16*f(i,j))*dx1*dy1/9;
    end
end

end