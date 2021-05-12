function dX = EOM(t,x, FdX, FdY,F, mesh_parameters, geometry_parameters, operational_parameters)

dX=zeros(4,1);

X=x(1);
Xdot=x(2);
Y=x(3);
Ydot=x(4);

[FX,FY]= CalculateLoadCapacity(mesh_parameters,geometry_parameters, operational_parameters,[X Xdot Y Ydot]);

g=9.8;
m=F/g;
dX(1)=x(2);
dX(2)=(FX+FdX)/m;
dX(3)=x(4);
dX(4)=(-m*g+FY+FdY)/m;
end