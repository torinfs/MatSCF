function [ out ] = tdhf( in, nel )
%TDHF Calculates resonant absorbtion frequencies
%
%   Inputs:
%       in
%       nel --  number of electrons in system
%
%   Outputs:
%       out
%

cAO = in.C;
eriAO = in.Vee;
for i=1:numel(in.epsilon)
    eps(2*i-1:2*i) = in.epsilon(i);
end
dim = length(eriAO);

temp = zeros(dim,dim,dim,dim);
temp2 = zeros(dim,dim,dim,dim);
temp3 = zeros(dim,dim,dim,dim);
eriMO = zeros(dim,dim,dim,dim);

fprintf('Transforming ERIs...\n')

% for p = 1:dim
%    for mu = 1:dim
%       temp(p,:,:,:) = temp(p,:,:,:) + cAO(mu,p)*eriAO(mu,:,:,:);
%    end
%    
%    for q = 1:dim
%        for nu = 1:dim
%            temp2(p,q,:,:) = temp2(p,q,:,:) + cAO(nu,q)*temp(p,nu,:,:);
%        end
%        
%        for r = 1:dim
%             for lam = 1:dim
%                 temp3(p,q,r,:) = temp3(p,q,r,:) + cAO(lam,r)*temp2(p,q,lam,:);
%             end
%             
%             for s = 1:dim
%                 for kap = 1:dim
%                     eriMO(p,q,r,s) = eriMO(p,q,r,s) + cAO(kap,s)*...
%                                                         temp3(p,q,r,kap);
%                 end
%                 
%             end
%        end
%    end
% end

for p = 1:dim
    for q = 1:dim
        for r = 1:dim
            for s = 1:dim
                for mu = 1:dim
                    for nu = 1:dim
                        for lam = 1:dim
                            for kap = 1:dim
                                eriMO(p,q,r,s) = cAO(mu, p)*cAO(nu, q)*...
                                    cAO(lam, r)*cAO(kap, s)*eriAO(mu,nu,lam,kap);
                            end
                        end
                    end
                end
            end
        end
    end
end


sdim = dim*2;

seri = zeros(sdim,sdim,sdim,sdim);

for p = 1:sdim
    ps = fs(p);
    for q = 1:sdim
        qs = fs(q);
        for r = 1:sdim
            rs = fs(r);
            for s = 1:sdim
                j = eriMO(ps, rs, qs, fs(s))*(mod(p,2) == mod(r,2))*...
                        (mod(q,2) == mod(s,2));
                k = eriMO(ps, qs, fs(s), rs)*(mod(p,2) == mod(s,2))*...
                        (mod(q,2) == mod(r,2));
                %spin eri is double bar integral; i.e. <pq||rs>
                seri(p,q,r,s) = j - k;
            end
        end
    end
end

% Begin TDHF
A = zeros(sdim);
B = zeros(sdim);

ia = 0;

fprintf('Forming A and B...\n')
for i = 1:nel
    for a = nel+1:sdim
        ia = ia + 1;
        jb = 0;
        
        for j = 1:nel
            for b = nel+1:sdim
                jb = jb + 1;
                %A = (e_a - e_i) d_{ij} d_{ab} * <aj||ib>
                A(ia,jb) = (eps(a) - eps(i)) * (i == j) * (a == b) + ...
                    seri(a,j,i,b);
                %B = <ab||ij>
                B(ia,jb) = seri(a,b,i,j);
            end
        end 
    end
end

fprintf('Solving TDHF...\n')
M = (A+B)*(A-B);
M = [A, B; -B, -A];
[ctd, etd] = eig(M);

% fprintf('Solving CIS...\n')
% [ctd, etd] = eig(A);

% sort
[etd, eI] = sort(diag(etd));
ctd = real(ctd(:,eI));

out = {ctd, etd};

fprintf('TDHF Done\n')
end

function y = fs( x )
    y = floor((x+1)/2);
end
