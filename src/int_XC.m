function [ Vxc, Exc, rhoInt ] = int_XC( basis, P, DGrid )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%Slater Exchange functional
rhoInt = 0;
Vxc = zeros(size(out.basis));
Exc = 0;

for iGrid = 1:numel(DGrid)
    r = DGrid(iGrid).xyz;
    w = DGrid(iGrid).weights;
    rho = 0;

    base = zeros(numel(out.basis));
    
    for mu = 1:numel(out.basis)
        
        musum = 0;
        xminA = r(1) - out.basis(mu).A(1);
        yminA = r(2) - out.basis(mu).A(2);
        zminA = r(3) - out.basis(mu).A(3);

        for k = 1:numel(out.basis(mu).d)
            musum = musum + out.basis(mu).d(k)*out.basis(mu).N(k)*...
                exp(-out.basis(mu).alpha(k)*sqrt(xminA^2 + yminA^2 +zminA^2));
        end

        base(mu) = xminA^(out.basis(mu).a(1)).*yminA^(out.basis(mu).a(2)).*...
            zminA^(out.basis(mu).a(3)).*musum;
    end
    
    for mu = 1:numel(out.basis)
        for nu = 1:numel(out.basis)
            rho = rho + out.P(mu,nu)*base(mu)*base(nu);
        end
    end
    
    rhoInt = rhoInt + w*rho;
    
    

end


end

%TODO: move the slater functions to their calling locations
function exs = xeslater(rho)
    exs = -3/4*(3/pi)^(1/3)*rho^(1/3);
end

function vxs = xvslater(rho)
    vxs = -(3/pi*rho)^1/3;
end

function [ecvwn, vcvwn] = vwn(rho, flag)
    if strcmp(flag, 'SVWN3')
        A = 0.310907;
        b = 13.0720;
        c = 42.7198;
        x_0 = -0.409286;
    elseif strcmp(flag, 'SVWN5')
        A = 0.310907;
        b = 3.72744;
        c = 12.9352;
        x_0 = -0.10498;
    else
        'Unsupported functional specified'
    end
    Q = sqrt(4*c-b^2);
    xi_0 = x_0^2 + b*x_0 + c;

    x = sqrt((3/(4*pi*rho))^(1/3));
    xi = x^2 + b*x + c;
    
    eta = arctan(Q/(2*x+b));
    
end


