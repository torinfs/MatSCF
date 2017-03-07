%%
%Test script for TDHF
clear; clc;
load('testcases_v04');

%%
clc;
tnum = 6;

options = struct('basisset', testcase(tnum).Basis,...
                 'tolEnergy', 1e-8,...
                 'tolDensity', 1e-8,...
                 'Method', 'HF');

out = mocalc(testcase(tnum).Elements, testcase(tnum).xyz,...
         testcase(tnum).TotalCharge, options);
     
outMP2 = mp2(out, sum(testcase(tnum).Elements)-testcase(tnum).TotalCharge);
outTDHF = tdhf(out, sum(testcase(tnum).Elements)-testcase(tnum).TotalCharge);