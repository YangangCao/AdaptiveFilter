% batch.m
make_rp;
rp.TNRdB = 10; rp.INRdB = 40; rp.Nsnaps= 200; rp.mu = 1e-10; rp.name = 'run1';
rp.sin_theta_2 = 0.05;run_qrd_rls_mvdr(rp);
rp.TNRdB = 10; rp.INRdB = 40; rp.Nsnaps= 200; rp.mu = 1e-10; rp.name = 'run2';
rp.sin_theta_2 = 0.15;run_qrd_rls_mvdr(rp);
rp.TNRdB = 10; rp.INRdB = 40; rp.Nsnaps= 200; rp.mu = 1e-10; rp.name = 'run3';
rp.sin_theta_2 = 0.2;run_qrd_rls_mvdr(rp);
plot_15_11;