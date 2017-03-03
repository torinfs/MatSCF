%%
%Test script for TDHF
clear; clc;
load('testcases_v04');

%%
tnum= 3;

options = struct('basisset', testcase(tnum).Basis,...
                 'tolEnergy', 1e-8,...
                 'tolDensity', 1e-8,...
                 'Method', 'HF');
             
xyz = [0 0 0; 0 0 0.929];
             
out = mocalc(testcase(tnum).Elements, testcase(tnum).xyz,...
         testcase(tnum).TotalCharge, options);
     
out2 = mp2(out, sum(testcase(tnum).Elements)-testcase(tnum).TotalCharge)