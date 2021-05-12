function [MeshX,MeshDX, MeshZ, MeshDZ] = initialisemesh(geometry_parameters,mesh_parameters)
            MeshX=linspace(0,2*pi*geometry_parameters(2),mesh_parameters(1));
            MeshDX=MeshX(2)-MeshX(1);
            MeshZ=linspace(0,geometry_parameters(1),mesh_parameters(2));
            MeshDZ=MeshZ(2)-MeshZ(1);
            
end
