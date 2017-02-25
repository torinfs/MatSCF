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

eriAO = in.Vee;
dim = length(in.C);
cAO = in.C;
eps = zeros(2*dim,1);
for i = 1:dim
    eps(2*i-1:2*i) = in.epsilon(i);
end

temp = zeros(dim,dim,dim,dim);
temp2 = zeros(dim,dim,dim,dim);
temp3 = zeros(dim,dim,dim,dim);
eriMO = zeros(dim,dim,dim,dim);

fprintf('Transforming ERIs...\n')

for p = 1:dim
   for mu = 1:dim
      temp(p,:,:,:) = temp(p,:,:,:) + cAO(p,mu)*eriAO(mu,:,:,:);
   end
   
   for q = 1:dim
       for nu = 1:dim
           temp2(p,q,:,:) = temp2(p,q,:,:) + cAO(q,nu)*temp(p,nu,:,:);
       end
       
       for r = 1:dim
            for lam = 1:dim
                temp3(p,q,r,:) = temp(p,q,r,:) + cAO(r,lam)*temp2(p,q,lam,:);
            end
            
            for s = 1:dim
                for kap = 1:dim
                    eriMO(p,q,r,s) = eriMO(p,q,r,s) + cAO(s,kap)*...
                                                        temp3(p,q,r,kap);
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
                j = eriMO(ps, qs, rs, fs(s))*(mod(p,2) == mod(r,2))*...
                        (mod(q,2) == mod(s,2));
                k = eriMO(ps, qs, fs(s), rs)*(mod(p,2) == mod(s,2))*...
                        (mod(q,2) == mod(r,2));
                seri(p,q,r,s) = j - k;
            end
        end
    end
end

A = zeros(sdim);
B = zeros(sdim);

I = 0;

fprintf('Forming A and B...\n')
for i = 1:nel
    for a = nel+1:sdim
        I = I + 1;
        J = 0;
        
        for j = 1:nel
            for b = nel+1:sdim
                J = J + 1;
                A(I,J) = (eps(a) - eps(i)) * (i == j) * (a == b) + ...
                    seri(a,j,i,b);
                B(I,J) = seri(a,b,i,j);
            end
        end
        
    end
end

fprintf('Solving TDHF...\n')
[ctd, etd] = eig([A,B;-B,-A]);

%sort
[etd, eI] = sort(real(diag(etd)));
ctd = real(ctd(:,eI));

out = {ctd, etd};
end

function y = fs( x )
    y = floor((x+1)/2);
end

