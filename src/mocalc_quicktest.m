%% mocalc_quicktest

clear; clc;
load('testcases_v04.mat');
water = testcase(6);

options = struct(...
                'basisset',     water.Basis,...
                'tolEnergy',    0.01,...
                'tolDensity',   0.01);
                


out = mocalc(water.Elements,water.xyz,water.TotalCharge,options);