clear, clc

load testcases_v04
nTests = numel(testcase);

tic
for iTest = 1:nTests
  t = testcase(iTest);
  
  % Get input data for moplot() from testcase data
  atoms = t.Elements;
  xyz_A = t.xyz; % angstroms
  charge = t.TotalCharge;
  options.basisset = t.Basis;
  options.tolEnergy = 1e-8; % hartrees
  options.tolDensity = 1e-7; % hartrees
  doDFT = t.DFT.doDFT;
  if doDFT
    options.Method = 'KS';
    options.ExchFunctional = t.DFT.ExchKernel;
    options.CorrFunctional = t.DFT.CorrKernel;
    options.nRadialPoints = 50;
    options.nAngularPoints = 302;
  else
    options.Method = 'HF';
    options.ExchFunctional = '';
    options.CorrFunctional = '';
  end
  
  % Run calculation
  result = mocalc(atoms,xyz_A,charge,options);
  
  % Calculate error
  err = @(field) max(abs(result.(field)(:)-t.(field)(:)));
  Etoterr(iTest) = err('Etot');
    
  % Print results
  fprintf('Testcase %d:  %s, %s, %s\n',iTest,t.MolName,options.Method, t.Basis);
  fprintf('  Etot error: %0.6f eV\n',Etoterr(iTest));
end
totaltime = toc;
fprintf('TOTAL TIME:  %g s\n',totaltime);
fprintf('MAX ERROR:   %g eV\n',max(Etoterr));
