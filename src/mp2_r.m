function [ out ] = mp2_r( in, nel )
%Caclulates the correlation energy according to 2nd order moller plesset PT

cAO = in.C;
eriAO = in.Vee;
eps = in.epsilon;
dim = length(eriAO);

fprintf('Solving MP2 Energy...\n')
eriMO = zeros(dim,dim,dim,dim);
for p = 1:dim
    for q = 1:dim
        for r = 1:dim
            for s = 1:dim
                temp_eri = 0;
                for mu = 1:dim
                    for nu = 1:dim
                        for lam = 1:dim
                            for kap = 1:dim
                                temp_eri = temp_eri + (cAO(mu, p)*cAO(nu, q)*...
                                    cAO(lam, r)*cAO(kap, s)*eriAO(mu,nu,lam,kap));
                            end
                        end
                    end
                end
                eriMO(p,q,r,s) = temp_eri;
            end
        end
    end
end

emp2 = 0;
for i = 1:nel/2
    for j = 1:nel/2
        for a = nel/2+1:dim
            for b = nel/2+1:dim
                emp2 = emp2 + eriMO(i,a,j,b)*(2*eriMO(i,a,j,b)...
                    - eriMO(i,b,j,a))/(eps(i) + eps(j) - eps(a) - eps(b));
            end
        end
    end
end


out = emp2;

end


function y = fs( x )
    y = floor((x+1)/2);
end