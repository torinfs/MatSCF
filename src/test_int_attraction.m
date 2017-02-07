%Test script for int_attraction and int_overlap
clear; clc;
load('testcases_v04.mat');
numtest = numel(testcase);
Vpassed = true;
Spassed = true;
Tpassed = true;
ErrList1 = {};
ErrList2 = {};
ErrList3 = {};
for molnum = 1:numtest
    basis = makebasis(testcase(molnum).Elements, testcase(molnum).xyz, basisread(testcase(molnum).Basis));
    Vne = int_attraction(basis);
    S = int_overlap(basis);
    T = int_kinenergy(basis);
    test1 = find(abs(testcase(molnum).Vne - Vne) > 1e-12);
    test2 = find(abs(testcase(molnum).S - S) > 1e-7);
    test3 = find(abs(testcase(molnum).T - T) > 1e-7);
    %testcase(molnum).MolName
    %isempty(test3)
    if ~isempty(test1)
        Vpassed = false;
        ErrList1{end+1} = struct('MolName', testcase(molnum).MolName,...
                            'trueVne', testcase(molnum).Vne,...
                            'calcVne', Vne);
    end
    if ~isempty(test2)
        Spassed = false;
        ErrList2{end+1} = struct('MolName', testcase(molnum).MolName,...
                            'trueS', testcase(molnum).S,...
                            'calcS', S);
    end
    if ~isempty(test3)
        Tpassed = false;
        ErrList3{end+1} = struct('MolName', testcase(molnum).MolName,...
                            'trueT', testcase(molnum).T,...
                            'calcT', T);
    end
end

if Vpassed
    'int_attraction tests passed!'
else
    'More work to do on attraction'
end
if Spassed
    'int_overlap tests passed!'
else
    'More work to do on overlap'
end
if Tpassed
    'int_kinenergy tests passed!'
else
    'More work to do on KE'
end