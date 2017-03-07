function [ out ] = mp2( in, nel )
%Caclulates the correlation energy according to 2nd order moller plesset PT

cAO = in.C;
eriAO = in.Vee;
for i=1:numel(in.epsilon)
    eps(2*i-1:2*i) = in.epsilon(i);
end
dim = length(eriAO);

fprintf('Solving MP2 Energy...\n')

temp = zeros(dim,dim,dim,dim);
temp2 = zeros(dim,dim,dim,dim);
temp3 = zeros(dim,dim,dim,dim);
eriMO = zeros(dim,dim,dim,dim);

for p = 1:dim
   for mu = 1:dim
      temp(p,:,:,:) = temp(p,:,:,:) + cAO(mu, p)*eriAO(mu,:,:,:);
   end
   
   for q = 1:dim
       for nu = 1:dim
           temp2(p,q,:,:) = temp2(p,q,:,:) + cAO(nu, q)*temp(p,nu,:,:);
       end
       
       for r = 1:dim
            for lam = 1:dim
                temp3(p,q,r,:) = temp3(p,q,r,:) + cAO(lam, r)*temp2(p,q,lam,:);
            end
            
            for s = 1:dim
                for kap = 1:dim
                    eriMO(p,q,r,s) = eriMO(p,q,r,s) + cAO(kap, s)*...
                                                        temp3(p,q,r,kap);
                end
                
            end
            
       end
   end
end


sdim = 2*dim;
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
                k = eriMO(ps, fs(s), qs, rs)*(mod(p,2) == mod(s,2))*...
                        (mod(q,2) == mod(r,2));
                %spin eri is double bar integral; i.e. <pq||rs>
                seri(p,q,r,s) = j - k;
            end
        end
    end
end

emp2 = 0;
for i = 1:nel
    for j = 1:nel
        for a = nel+1:sdim
            for b = nel+1:sdim
                emp2 = emp2 + 0.25*seri(i,j,a,b)*seri(i,j,a,b)/...
                    (eps(i) + eps(j) - eps(a) - eps(b));
            end
        end
    end
end


out = emp2;

fprintf('----------------------------------------\n')
form2 = 'Correlation Energy (MP2): %d.\n';
str = sprintf(form2,emp2);
fprintf(str)
fprintf('----------------------------------------\n')

end


function y = fs( x )
    y = floor((x+1)/2);
end