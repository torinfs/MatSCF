function [ Vxc, Exc, rhoInt ] = int_XC( basis, P, DGrid, flag )
% INT_XC() Calulates the exchange-correlation functional for DFT using the
% LDA approximation.
%   INPUTS:
%       basis = basis set from structure array of mocalc output
%       P = MxM Density matrix
%       DGrid = molcular grid
%       flag = string specifying functional i.e 'SVWN3'
%
%   OUTPUTS:
%       Vxc = MxM exchange-correlation potential matrix
%       Exc = exchange-correlation energy
%       rhoInt = Integral of density over allspace. Equal to N electrons.
%
%

rhoInt = 0;
Vxc = zeros(size(basis));
Exc = 0;

% Calculate rho(r) from grid
for iGrid = 1:numel(DGrid)
    r = DGrid(iGrid).xyz;
    w = DGrid(iGrid).weights;
    rho = 0;

    base = zeros(numel(basis));
    
    for mu = 1:numel(basis)
        
        musum = 0;
        xminA = r(1) - basis(mu).A(1);
        yminA = r(2) - basis(mu).A(2);
        zminA = r(3) - basis(mu).A(3);

        for k = 1:numel(basis(mu).d)
            musum = musum + basis(mu).d(k)*basis(mu).N(k)*...
                exp(-basis(mu).alpha(k)*sqrt(xminA^2 + yminA^2 +zminA^2));
        end

        base(mu) = xminA^(basis(mu).a(1)).*yminA^(basis(mu).a(2)).*...
            zminA^(basis(mu).a(3)).*musum;
    end
    
    for mu = 1:numel(basis)
        for nu = 1:numel(basis)
            rho = rho + P(mu,nu)*base(mu)*base(nu);
        end
    end
    
    % Final rhoInt value
    rhoInt = rhoInt + w*rho;
end

% get correlation values
[Ecvwn, Vcvwn] = VWN(rho, flag);

% XC energy value
Exc = xeslater(rho) + Ecvwn;




end

%TODO: move the slater functions to their calling locations
function exs = xeslater(rho)
    exs = -3/4*(3/pi)^(1/3)*rho^(1/3);
end

function vxs = xvslater(rho)
    vxs = -(3/pi*rho)^1/3;
end

% Correlation functionals
function [Ecvwn, Vcvwn] = VWN(rho, flag)
    if strcmp(flag, 'SVWN3')
        A = 0.310907;
        b = 13.0720;
        c = 42.7198;
        x0 = -0.409286;
    elseif strcmp(flag, 'SVWN5')
        A = 0.310907;
        b = 3.72744;
        c = 12.9352;
        x0 = -0.10498;
    else
        disp('Unsupported functional specified')
        break
    end
    Q = sqrt(4*c-b^2);
    xi_0 = x0^2 + b*x0 + c;
    x = sqrt((3/(4*pi*rho))^(1/3));
    xi = x^2 + b*x + c;
    eta = arctan(Q/(2*x+b));
    
    Ecvwn = A*( log((x^2)/xi) + 2*b*eta/Q - (b*x0/xi_0) * ...
                    ( log(((x-x0)^2)/xi) + (2*(2*x0 + b)*eta)/Q));
                
    Vcvwm = -(1/3) * ( (A*(c*(x-x0) - b*x*x0)) / (xi*(x-x0)) ) + Ecvwn;
    
end

















