%Test script for makebasis()
%%
%Section 1 - Definition of test molecules

%H2, CH4, CO bond distances and angles from CRC Handbook of chemistry and 
%   Physics - 96th edition.
H2BondDist = 0.74144;
H2 = struct('atoms', [1,1], 'xyz', [0,0,0; 0,0,H2BondDist]);

CH4Dist = 1.0870;
CH4Ang = 109.5;
CH4 = struct('atoms', [6,1,1,1,1], 'xyz',[0, 0, 0;...
      0,    0,      CH4Dist;...
      CH4Dist*sin(deg2rad(CH4Ang)), 0, CH4Dist*cos(deg2rad(CH4Ang));...
      CH4Dist*sin(deg2rad(CH4Ang))*cos(2*pi/3),...
      CH4Dist*sin(deg2rad(CH4Ang))*sin(2*pi/3),...
      CH4Dist*cos(deg2rad(CH4Ang));...
      CH4Dist*sin(deg2rad(CH4Ang))*cos(4*pi/3),...
      CH4Dist*sin(deg2rad(CH4Ang))*sin(4*pi/3),...
      CH4Dist*cos(deg2rad(CH4Ang))]);

Be = struct('atoms', [3], 'xyz', [0,0,0]);

CODist = 1.283;
CO = struct('atoms', [6,8], 'xyz', [0,0,0; 0,0,CODist]);

%HHe+ distance from http://link.springer.com/article/10.1007%2Fs00894-008-0371-3
%   Accessed 1/29/17. doi:10.1007/s00894-008-0371-3
HHeDist = 0.772;
HHe = struct('atoms', [2,1], 'xyz', [0,0,0;0,0,HHeDist]);

AtomList = {Be, H2, HHe, CO, CH4};

%%
%Section 2 - Loading of basis sets
b6_31G = basisread('6-31G');
b6_311G = basisread('6-311G');
bSTO_3G = basisread('STO-3G');
bcc_pVDZ = basisread('cc-pVDZ');
bases = {bSTO_3G, b6_31G, b6_311G, bcc_pVDZ};

%%
%Section 3 - Definition of comparison values
%Rows are basis sets, columns are molecules. The value is expected number
%   of basis sets defined.
truth = [5,  2,  2,  10,  9;...
         9,  4,  4,  18,  17;...
         13,  6,  6,  26,  25;...
         15, 10, 10, 30,  35];
     
%%
%Section 4 - testing of function
runtime = zeros(4,5);
numBases = zeros(4,5);
for i = 1:4
    for j = 1:5
        tic
        testBasis = makebasis(AtomList{j}.atoms, AtomList{j}.xyz, bases{i});
        runtime(i,j) = toc;
        numBases(i,j) = numel(testBasis);
    end
end

if numBases == truth
    'makebasis() succesful'
else
    'makebasis() failed'
end
        