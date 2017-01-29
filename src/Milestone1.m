%% Milestone 1
clear; clc;

% Test with water molecule
basissetdef = basisread('6-31G');
water = [8, 1, 1];
coords = [ 0.000000,  0.000000,  0.000000; 
           0.758602,  0.000000,  0.504284; 
           0.758602,  0.000000,  -0.504284]; 
        
makebasis(water,coords,basissetdef)
