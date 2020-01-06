function [Psinew] = func_KS(Psi)

                                                                                                                             c
% solution of household model
[gridx,gridsav,gridass,cfun] = func_hh(Psi);

% Simulation
[KK] = func_stochsim(gridx,gridsav,cfun,gridass);

% Regression
[Psinew] = func_regr(KK);


end