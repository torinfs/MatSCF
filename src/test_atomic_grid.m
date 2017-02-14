%%
%test script for atomic_grid()
clear; clc;
load('testcases_v04.mat');

%%
i = 5;

for j = 1:numel(testcase(i).Elements)
    agrid = atomic_grid(j, testcase(i).xyz,testcase(i).Elements,...
        3,38);
end