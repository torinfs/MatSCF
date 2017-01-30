%Takes a list of atoms, their cartesian coordinates, and the basis set
%structure as provided by basisread(). Returns a structure array of all the
%basis functions.
%
%   Inputs:
%       atoms-       1xK list of atomic numbers (e.g. [6 8] for CO)
%       xyz-         Kx3 array of nulcear coordinates in angstroms. Must
%                    match order of nuclei in atoms list.
%       basissetdef- basis set structure as read in by basisread().
%                    For optimal memory performance, basisread() should be
%                    called in the argument list of makebasis(), so that it
%                    is deallocated after makebasis is finished.
%
%   Outputs:
%    basis-          1xM structure array (where M is the number of basis
%                    functions. basis contains the following fields.
%     basis(p).atom  Atomic number of central nucleus (e.g. 6 for C).
%     basis(p).A     Vector of cartesian coordinates of cetral nucleus in
%                    Bohr.
%     basis(p).a     Vector of x,y,z angular ("Cartesian") exponents
%     basis(p).alpha Array of radial exponents of primitives making up
%                    the basis function.
%     basis(p).d     Array of contraction coefficents of the primitives.
%     basis(p).N     Array of normalization constants.
%
%Support notes:
%   1) Only basis functions with 'S', 'P', 'D', or 'SP' type subshells are
%       currently supported.


function basis = makebasis(atoms, xyz, basissetdef)
%
%Stupid conversion factor I'll fix later
bohr2ang = 0.52917721067;
%Because of the various number of basis functions defined by each shell
%type, it is simpler to use a basis function counter, m.
m = 1;

%Loop through each atom in the atoms list, and generate the basis functions
%for each atom

%Preallocation loop. This loop just gets the number of basis functions and
%stores that number in n.
elatms = numel(atoms);
n=0;
for i = 1:elatms
    basisA = basissetdef{atoms(i)};
    elbasis = numel(basisA);
    for j = 1:elbasis
        switch basisA(j).shelltype
            case 'S'
                n = n+1;
            case 'P'
                n = n+3;
            case 'SP'
                n = n+4;
            case 'D'
        end
    end
end

%Actually preallocating the structure array. I'm not sure how MATLAB
%   implements its memory storage, but preallocating the maximum size
%   supported (8 primitives in cc-pVDZ) will likely cause less
%   fragmentation (even though some internal fragmentation is guaranteed)
%   and reduce time searching for an appropriate memory block.

basis(n) = struct('atom', 0,...
                  'A', [-99,-99,-99],...
                  'a', [-1,-1,-1],...
                  'alpha', [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],...
                  'd', [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],...
                  'N', [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]);

%NOTE: Consider timing makebasis with and without the preallocation loop


for i = 1:elatms
    
    %Get the information for a specific atom
    atomid = atoms(i);
    xyzi = xyz(i)/bohr2ang;     %Convert the xyz coordinates to bohr
    basisA = basissetdef{atomid};
    
    %Loop through the basissetdef for a specific atom and define each basis
    %   function.
    elbasis = numel(basisA);
    for j = 1:elbasis
        
        %Filter basis functions by shell type
        %NOTE - only basis functions with 'S', 'P', 'D', or 'SP' character
        %   can be constructed.
        switch basisA(j).shelltype
            
            %S shells only consist of one basis function
            case 'S'
                basis(m).atom = atomid;
                basis(m).A = xyzi;
                basis(m).a = [0,0,0];
                basis(m).alpha = basisA(j).exponents;
                basis(m).d = basisA(j).coeffs;
                basis(m).N = normConst([0,0,0], basisA(j).exponents);
                m = m+1;
                
            %P shells consist of three functions, but the only difference
            %   between the three is the angular exponents.
            case 'P'
                basis(m).atom = atomid;
                basis(m).A = xyzi;
                basis(m).a = [1,0,0];
                basis(m).alpha = basisA(j).exponents;
                basis(m).d = basisA(j).coeffs;
                basis(m).N = normConst([1,0,0], basisA(j).exponents);
                m = m+1;
                for k = 2:3
                    basis(m) = basis(m-1);
                    basis(m).a(k-1) = 0;
                    basis(m).a(k) = 1;
                    m = m+1;
                end
                
            %SP shells consist of four functions - One S type and two P.
            %The easiest way to deal with these is to define S separate
            %   from P and reuse the code above, with minor changes due to
            %   a different structure in basissetdef.
            case 'SP'
                %S definition
                basis(m).atom = atomid;
                basis(m).A = xyzi;
                basis(m).a = [0,0,0];
                basis(m).alpha = basisA(j).exponents;
                basis(m).d = basisA(j).coeffs(1,:);
                basis(m).N = normConst([0,0,0], basisA(j).exponents);
                m = m+1;
                
                %P definition
                basis(m).atom = atomid;
                basis(m).A = xyzi;
                basis(m).a = [1,0,0];
                basis(m).alpha = basisA(j).exponents;
                basis(m).d = basisA(j).coeffs(2,:);
                basis(m).N = normConst([1,0,0], basisA(j).exponents);
                m = m+1;
                for k = 2:3
                    basis(m) = basis(m-1);
                    basis(m).a(k-1) = 0;
                    basis(m).a(k) = 1;
                    m = m+1;
                end
                
            %D shells consist of 6 basis functions, but they can be grouped
            %   into two sets of three. The two sets differ in
            %   normalization constant, and the within each set, the only
            %   difference is in the angular exponents
            case 'D'
                %D basis functions quadratic in one coordinate
                basis(m).atom = atomid;
                basis(m).A = xyzi;
                basis(m).a = [2,0,0];
                basis(m).alpha = basisA(j).exponents;
                basis(m).d = basisA(j).coeffs;
                basis(m).N = normConst([2,0,0], basisA(j).exponents);
                m = m+1;
                for k = 2:3
                    basis(m) = basis(m-1);
                    basis(m).a(k-1) = 0;
                    basis(m).a(k) = 2;
                    m = m+1;
                end
                
                %D basis functions linear in two coordinates
                basis(m) = basis(m-1);
                basis(m).a = [0,1,1];
                basis(m).N = normConst([0,1,1], basisA(j).exponents);
                m = m+1;
                for k = 2:3
                    basis(m) = basis(m-1);
                    basis(m).a(k-1) = 1;
                    basis(m).a(k) = 0;
                    m = m+1;
                end    
        end
    end
end

end


%Takes a list of angular (cartesian) exponents, carts, a list of radial
%   exponents, alphas, and returns a list of normalization constants
%   corresponding to the radial exponents.
function N = normConst(carts, alphas)
    numa = numel(alphas);
    N = zeros(1,numa);
    for i = 1:numa
        const = (2/pi)^(3/2);
        numer = 2^(sum(carts))*alphas(i)^((2*(sum(carts))+3)/4);
        denom = sqrt(fact2(2*carts(1)-1)*...
                     fact2(2*carts(2)-1)*fact2(2*carts(3)-1));
        N(i) = const*numer/denom;
    end
end

%Takes a scalar, x, and returns the semifactorial of x
function xfac2 = fact2(x)
    xfac2 = prod(x:-2:1);
end