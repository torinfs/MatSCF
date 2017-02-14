function R = atomicradius(iElements)

% List of covalent radii (Cordero, 2008; based on 37k crystal structures)
% DOI: http://dx.doi.org/10.1039/B801115J
R = [0.31 0.28 1.28 0.96 0.84 (0.76+0.73+0.69)/3 0.71 0.66 0.57 0.58]; % H-He, Li-Ne
% (Becke uses Bragg-Slater radii (Slater, 1964))

R = R(iElements);

bohr = 0.52917721067; % Bohr radius, in Angstrom
R = R / bohr; % Angstrom -> bohr
