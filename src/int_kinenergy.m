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

                Kab = exp(-((alpha(k) * beta(l)) / p)...
                    * norm(A - B)^2);
                
                % Here w = [x, y, z]
                IntS_1D = zeros(3,1);
                IntS_1D_pos = zeros(3,1);
                IntS_1D_neg = zeros(3,1);
                for w = 1:3
                    
                    % regular S
                    for i = 0:((a(w)+b(w))/2)
                        f = f_w(a(w), b(w), P(w), A(w), B(w), 2*i);
                        IntS_1D(w,:) = IntS_1D(w,:) + (f * fact2nd(((2*i) -1)) / ...
                            ((2*p)^i));   

                    end
                    
                    % b + 2
                    for j = 0:((a(w)+b(w)+2)/2)
                        f_pos = f_w(a(w), b(w) + 2, P(w), A(w), B(w), 2*j);
                        IntS_1D_pos(w,:) = IntS_1D_pos(w,:) +... 
                            (f_pos * fact2nd(((2*j) -1)) /((2*p)^j));
                    
                    end
                    
                    % b - 2
                    % Check for negative b(w) term
                    if b(w) <= 1
                        IntS_1D_neg(w,:) = 0;
                    else
                        for q = 0:((a(w)+b(w)-2)/2)
                            f_neg = f_w(a(w), b(w) - 2, P(w), A(w), B(w), 2*q);
                            IntS_1D_neg(w,:) = IntS_1D_neg(w,:) + (f_neg * ...
                                fact2nd(((2*q) -1)) / ((2*p)^q));
                        end
                        
                    end
                    
                end
                % The equation for IntS is pulled from Eqn. 25
                T_prefactors = basis(u).d(k) * ...
                    basis(v).d(l) * basis(u).N(k) ...
                    * basis(v).N(l);
                
                
                % Bracket integrals with different b values i.e. [a|(b + 2)]
                
                % original overlap
                kalb =  ((pi./p).^(3/2)) * Kab * IntS_1D(1,:)...
                                * IntS_1D(2,:) * IntS_1D(3,:);
                
                % b + 2 for each x,y,z exponent            
                kalb_posx = ((pi./p).^(3/2)) * Kab * IntS_1D_pos(1,:)...
                                * IntS_1D(2,:) * IntS_1D(3,:);
                kalb_posy = ((pi./p).^(3/2)) * Kab * IntS_1D(1,:)...
                                * IntS_1D_pos(2,:) * IntS_1D(3,:);
                kalb_posz = ((pi./p).^(3/2)) * Kab * IntS_1D(1,:)...
                                * IntS_1D(2,:) * IntS_1D_pos(3,:);
                
                % b - 2 for each x,y,z exponent            
                kalb_negx = ((pi./p).^(3/2)) * Kab * IntS_1D_neg(1,:)...
                                * IntS_1D(2,:) * IntS_1D(3,:);
                kalb_negy = ((pi./p).^(3/2)) * Kab * IntS_1D(1,:)...
                                * IntS_1D_neg(2,:) * IntS_1D(3,:);
                kalb_negz = ((pi./p).^(3/2)) * Kab * IntS_1D(1,:)...
                                * IntS_1D(2,:) * IntS_1D_neg(3,:);
                
                % Kinetic Energy integrals defined in equation (31)
                Ix = beta(l)*(2*b(1)+1)*kalb - 2*(beta(l).^2)*kalb_posx...
                                        - (1/2)*b(1)*(b(1)-1)*kalb_negx;
                Iy = beta(l)*(2*b(2)+1)*kalb - 2*(beta(l).^2)*kalb_posy...
                                        - (1/2)*b(2)*(b(2)-1)*kalb_negy;
                Iz = beta(l)*(2*b(3)+1)*kalb - 2*(beta(l).^2)*kalb_posz...
                                        - (1/2)*b(3)*(b(3)-1)*kalb_negz;
                                    
                T(u,v) = T(u,v) + T_prefactors*(Ix + Iy + Iz);
            end
        end
    end
end

T = (T + triu(T,1)');
end

function y = f_w(a_w, b_w, P_w, A_w, B_w, K)
    y = 0;
    for j = (max(0,(K-a_w)):min(K,b_w))
        prefactor = nchoosek(a_w, K-j)*nchoosek(b_w, j);
        other = ((P_w-A_w)^(a_w-K+j))*((P_w-B_w)^(b_w-j));
        y_j = prefactor*other;
        y = y + y_j;
    end
end

function xfac2 = fact2nd(x)
    xfac2 = prod(x:-2:1);
end