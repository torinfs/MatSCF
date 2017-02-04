function T = int_kinenergy(basis)
% INT_KINENERGY: Write a function that calculates the MxM matrix of kinetic
%   energy integrals.  The function should return a MxM matrix of KE
%   integrals which is symmetric.  Math is pulled from "Molecular Orbitals 
%   from Scratch" handout.

% Preallocate for speed
T = zeros(length(basis));
blength = length(basis);

% Summing over basis set mu
for u = 1:blength
    A = basis(u).A;
    a = basis(u).a;
    alpha = basis(u).alpha;
    
    % Summing over basis set nu
    for v = u:blength
        B = basis(v).A;
        b = basis(v).a;
        beta = basis(v).alpha;
        
        % Summing over the length of thse basis sets
        for k = 1:length(alpha)
            for l = 1:length(beta)
                
                % Since we don't use p, P, Kab outside of each element-wise
                % calculation, they live within these sets of loops.
                p = alpha(k) + beta(l);
                P = (alpha(k)*A + beta(l)*B)/p;
                
                % Why is this always 1 ?????
                Kab = exp(-((alpha(k) * beta(l)) / p)...
                    * abs((A - B)*(A - B)')); 
                
                % Here w = [x, y, z]
                for w = 1:3
                    for i = 0:((a(w)+b(w))/2)
                        f = f_w(a(w), b(w), P(w), A(w), B(w), 2*i);
                        
                        % While this is not the best name, this function
                        % calculates each 1D primative integral per row of
                        % a matrix.  Therefore, Prim1D(1,:) is the Ix for
                        % the overlap matrix, and so on.
                        IntS_1D(w,:) = (f * prod(((2*i) -1):-2:1)) / ...
                            ((2*p)^i);   

                    end
                    for j = 0:((a(w)+b(w)+2)/2)
                        f_pos = f_w(a(w), b(w) + 2, P(w), A(w), B(w), 2*j);
                        IntS_1D_pos(w,:) = (f_pos * prod(((2*j) -1):-2:1)) / ...
                                    ((2*p)^j);
                    
                    end
                    
                    % Check for negative b(w)
                    if b(w) < 2
                        IntS_1D_neg(w,:) = 0;
                    else
                        for q = 0:((a(w)+b(w)-2)/2)
                            f_neg = f_w(a(w), b(w) - 2, P(w), A(w), B(w), 2*q);
                            IntS_1D_neg(w,:) = (f_neg * prod(((2*q) -1):-2:1)) / ...
                                                ((2*p)^q);
                    end
                    
                        
                    end
                end
                % The equation for IntS is pulled from Eqn. 25
                T_prefactors = basis(u).d(k) * ...
                    basis(v).d(l) * basis(u).N(k) ...
                    * basis(v).N(l);
                
                % Bracket integrals with different b values i.e. [a|(b + 2)]
                kalb =  ((pi./p).^(3/2)) * Kab * IntS_1D(1,:)...
                                * IntS_1D(2,:) * IntS_1D(3,:);
                kalb_pos = ((pi./p).^(3/2)) * Kab *IntS_1D_pos(1,:)...
                                * IntS_1D_pos(2,:) * IntS_1D_pos(3,:);
                kalb_neg = ((pi./p).^(3/2)) * Kab * IntS_1D_neg(1,:)...
                                * IntS_1D_neg(2,:) * IntS_1D_neg(3,:);
                
                % Kinetic Energy integrals defined in equation (31)
                Ix = beta(l)*(2*b(1)+1)*kalb - 2*(beta(l).^2)*kalb_pos...
                                        - (1/2)*b(1)*(b(1)-1)*kalb_neg;
                Iy = beta(l)*(2*b(2)+1)*kalb - 2*(beta(l).^2)*kalb_pos...
                                        - (1/2)*b(2)*(b(2)-1)*kalb_neg;
                Iz = beta(l)*(2*b(3)+1)*kalb - 2*(beta(l).^2)*kalb_pos...
                                        - (1/2)*b(3)*(b(3)-1)*kalb_neg;
                T(u,v) = T_prefactors*(Ix + Iy + Iz);
            end
        end
    end
    %Nvec(u) = sqrt(1/T(u,u));
end

T = T + triu(T,1)';
end

function [y] = f_w(a_w, b_w, P_w, A_w, B_w, K)
    y = 0;
    for j = (max(0,(K-a_w)):min(K,b_w))
        prefactor = nchoosek(a_w, K-j)*nchoosek(b_w, j);
        other = (P_w-A_w)^(a_w-K+j)*(P_w-B_w)^(b_w-j);
        y_j = prefactor*other;
        y = y + y_j;
    end
end