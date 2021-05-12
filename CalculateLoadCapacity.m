function [FX,FY]= CalculateLoadCapacity(mesh_parameters,geometry_parameters, operational_parameters,initial_parameters)

X=initial_parameters(1);
Xdot=initial_parameters(2);
Y=initial_parameters(3);
Ydot=initial_parameters(4);

L=geometry_parameters(1);
R=geometry_parameters(2);
h0=geometry_parameters(3);

[MeshX,MeshDX, MeshZ, MeshDZ]=initialisemesh(geometry_parameters,mesh_parameters);
[h, dhdx, dhdz] = initialiseradialgap(MeshX, MeshDX, MeshZ, MeshDX, geometry_parameters, initial_parameters);

omega = operational_parameters(1)*2*pi/60;
mu=operational_parameters(2);
PresCond=operational_parameters(3);

m=mesh_parameters(1);
n=mesh_parameters(2);

dx=MeshDX;
dz=MeshDZ;

U=repmat(omega*R+Xdot*cos(MeshX(:)/R)-Ydot*sin(MeshX(:)/R),1,numel(MeshZ));
V=repmat(Xdot*sin(MeshX(:)/R)+Ydot*cos(MeshX(:)/R),1,numel(MeshZ));

dU=XDerivative(U,dx);

A=zeros(m,n); D=zeros(m,n);
E=zeros(m,n); F=zeros(m,n);
G=zeros(m,n); H=zeros(m,n);
K=zeros(m,n);

for i=1:m
    for j=1:n
        A(i,j)=(h(i,j)^3)/mu;
        D(i,j)=6*(h(i,j)*dU(i,j)+U(i,j)*dhdx(i,j))-12*V(i,j);
    end
end % determination of matrix A

B=XDerivative(A,dx);
C=ZDerivative(A,dz);

for i=1:m
    for j=1:n
        E(i,j)=(B(i,j)/(2*dx))+(A(i,j)/dx^2);
        F(i,j)=(C(i,j)/(2*dz))+(A(i,j)/dz^2);
        G(i,j)=(-2*A(i,j)/(dx^2))+(-2*A(i,j)/dz^2);
        H(i,j)=(-B(i,j)/(2*dx))+(A(i,j)/dx^2);
        K(i,j)=(-C(i,j)/(2*dz))+(A(i,j)/dz^2);
    end
end
% determination of matrices E, F, G, H, K

M=n-2;
N=M*m;
pp=zeros(N);
FT=zeros(N,1);

for I=1:N
    i=ceil(I/M);
    j=I-(i-1)*M;
    
    FT(I,1)=D(i,j);
    if rem(I,M)==0
        FT(I,1)=D(i,j)-F(i,j)*PresCond;
        FT(I-M+1,1)=D(i,j)-K(i,j)*PresCond;
    end
    for J=1:N
        if J==I
            pp(I,J)=G(i,j);
            if (J+M)<=N
                pp(I,J+M)=E(i,j);
            else
                pp(I,J-N+2*M)=E(i,j);
            end
            if I<=M
                pp(I,N-2*M+J)=H(i,j);
            else
                pp(I,J-M)=H(i,j);
            end
            if rem(I,M)~=0 && I<=N
                pp(I+1,J)=K(i,j);
                pp(I,J+1)=F(i,j);
            end
        end
    end
end% filling in the pp matrix
pr=pp\FT;

p1=zeros(m,n-2);
for I=1:N
    ii=ceil(I/M);
    jj=I-(ii-1)*M;
    p1(ii,jj)=pr(I);
end% unwrapping the p1 matrix-solution
p000=ones(m,1)*PresCond;
p0=[p000 p1 p000];
for i=1:m
    for j=1:n
        if p0(i,j)<0, p0(i,j)=0; end % Guembel's condition for cavitation
    end
end

subRX=zeros(m,n); subRY=zeros(m,n);
for i=1:m
    for j=1:n
        subRX(i,j)=p0(i,j)*sin(MeshX(i)/R);
        subRY(i,j)=p0(i,j)*cos(MeshX(i)/R);
    end
end
FX=-intFdxdz(subRX,dx,dz);
FY=-intFdxdz(subRY,dx,dz);
R=sqrt(FX^2+FY^2);

end