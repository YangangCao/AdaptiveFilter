eval(['save ' rp.rname ' E Wx Xi rp'])
% make rp structure to be used for passing run
% parameters to the run_lms_pred algorithm
rp.Nruns = 100;
rp.Ndata = 100;
2
rp.mult = 200;
rp.verbose = 0;
rp.mu = 0.05;
rp.a = 0.99;
rp.rname = 'run1';
rp.var_v = 0.02;
rp.decay = 0;
Answers 5.20