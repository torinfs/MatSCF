%% Kinenergy test
clear; clc;    

% Test with water molecule
basissetdef = basisread('STO-3G');

load('testcases_v03.mat');
h2 = testcase(1);
water = testcase(6);
waterbasis = makebasis(water.Elements,water.xyz,basissetdef);
h2basis = makebasis(h2.Elements,h2.xyz,basissetdef);

% Kinetic Energy
h2.T
int_kinenergy(h2basis)
int_kinenergy(waterbasis)
water.T
%int_overlap(waterbasis) - water.S


% 
% function N = normConst(carts, alphas)
%     numa = numel(alphas);
%     N = zeros(1,numa);
%     for i = 1:numa
%         const = (2/pi)^(3/2);
%         numer = 2^(sum(carts))*alphas(i)^((2*(sum(carts))+3)/4);
%         denom = sqrt(fact2(2*carts(1)-1)*...
%                      fact2(2*carts(2)-1)*fact2(2*carts(3)-1));
%         N(i) = const*numer/denom;
%     end
% end