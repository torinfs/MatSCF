function [ S ] = int_overlap( basis )
% INT_OVERLAP: Write a function that calculates the MxM matrix of overlap
%   integrals.  The function should return a MxM matrix of overlap
%   integrals which has one on the diagonal and is symmetric.  Math is
%   pulled from "Molecular Orbitals from Scratch" handout.

% Preallocating space for the matrix.  No idea how big of a difference this
% actually makes, I just hate having the 'warnings' in the margins of
% Matlab.
f = zeros(1, 3);
S = zeros(length(basis), length(basis));
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
        for alphaLength = 1:length(alpha)
            for betaLength = 1:length(beta)
                
                % Since we don't use p, P, Kab outside of each element-wise
                % calculation, they live within these sets of loops.
                p = alpha(alphaLength) + beta(betaLength);
                P = (alpha(alphaLength)*A + beta(betaLength)*B)/p;
                Kab = exp(-((alpha(alphaLength) * beta(betaLength)) / p)...
                    * ((A - B)*(A - B)')); 
                
                % Here w = [x, y, z] as per Eqn. 27.  Consiquentially, this
                % returns a 
                Prim1D = zeros(3,1);
                for w = 1:3
                    for i = 0:((a(w)+b(w))/2)
                        f = f_w(a(w), b(w), P(w), A(w), B(w), 2*i);
                        
                        % While this is not the best name, this function
                        % calculates each 1D primative integral per row of
                        % a matrix.  Therefore, Prim1D(1,:) is the Ix for
                        % the overlap matrix, and so on.
                        Prim1D(w,:) = Prim1D(w,:) + (f * prod(((2*i) -1):-2:1)) / ...
                            ((2*p)^i);

                    end
                end
                % The equation for S is pulled from Eqn. 25
                S_prefactors = basis(u).d(alphaLength) * ...
                    basis(v).d(betaLength) * basis(u).N(alphaLength) ...
                    * basis(v).N(betaLength);
                
                S(u,v) = S(u,v) + S_prefactors * ((pi./p).^(3/2)) * Kab ...
                    * Prim1D(1,:) * Prim1D(2,:) * Prim1D(3,:);
            end
        end
    end
    Nvec(u) = sqrt(1/S(u,u));
end

% Symmetrizing the matrix
S = S + triu(S,1)';
for i = 1:blength
    for j = i:blength
        N(i,j) = Nvec(i)*Nvec(j);
    end
end
N = N + triu(N, 1)';
S = N.*S;
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