%% Kinenergy test
clear; clc;    

% Test with water molecule
basissetdef = basisread('STO-3G');

load('testcases_v03.mat');
water = testcase(6);
waterbasis = makebasis(water.Elements,water.xyz,basissetdef);

% Kinetic Energy
int_kinenergy(waterbasis)
water.T