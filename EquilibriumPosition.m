function [x_new, y_new] = EquilibriumPosition(F2,mesh_parameters,geometry_parameters, operational_parameters)
Force=F2;
%gamma1 - convergence criterion gamma1>=W-F/W
%gamma2 - convergence criterion gamma2>=Fx/Fy
%in the equilibrium position Fx = 0; Fy = F = W;

thetta=0;
gamma1=0.015; gamma2=0.015;
maxiterations=100; eh=0.9; el=0.1; e=(eh+el)/2;
message='Calculating equilibrium position...';
disp(message);
for i=1:maxiterations
    thetta_prev=thetta;
    for j=1:maxiterations
        e_prev=e;
        y_new=(cos(thetta)*e)*geometry_parameters(3);
        x_new=(sin(thetta)*e)*geometry_parameters(3);
        [FX,FY]= CalculateLoadCapacity(mesh_parameters,geometry_parameters, operational_parameters,[x_new 0 y_new 0]);
        ConvergenceCriterion1=abs((Force-FY)/Force)<=gamma1;
        if ConvergenceCriterion1==1, break;
        else if FY>Force, eh=e_prev; e=(eh+el)/2; else el=e_prev; e=(eh+el)/2;
            end
        end
    end
    ConvergenceCriterion2=abs((FX/FY))<=gamma2;
    if ConvergenceCriterion2==1, break;
    else thetta=thetta_prev-atan(FX/FY);
    end
end

end

