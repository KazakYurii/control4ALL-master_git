function [F,dFdx,dFdz] = initialiseradialgap(MeshX, MeshDX, MeshZ, MeshDZ, geometry_parameters, initial_parameters)
L=geometry_parameters(1);
R=geometry_parameters(2);
h0=geometry_parameters(3);

X0=initial_parameters(1);
Y0=initial_parameters(3);


F = repmat(h0 - X0*sin(MeshX(:)/R) - Y0*cos(MeshX(:)/R), 1, numel(MeshZ));

dFdx = XDerivative(F,MeshDX);
dFdz = ZDerivative(F,MeshDZ);
end

