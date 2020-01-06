function OLG_KS
%%% This code draws on a number of auxiliary functions. The three main
%%% .m-files are
%%% 1) func_hh: solves the HH problem for a given Psi
%%% 2) func_stochsim: (defcetive) simulates aggregate shocks and aggregates
%%%    for each shock to get a sequence of aggregate capital
%%% 3) func_regr: performs OLS regression to get new Psi
clear; clc;
close all



opt_det=false;          % 1=deterministic model
opt_nosr=false;         % 1=no survival risk
opt_ny = 2;             % 1=Markov chain with number of states, ny=5,
% 2=Markov chain with ny=2 (Krüger-Ludwig
% calibration)
% calibration

tol = 1e-4;
maxit = 100;
omega = 0.1;

% SOLUTION
func_calibr;
%%

for it=1:maxit,
    func_hh;
    func_stochsim;
    func_reg;
    
    diff = abs(PsiNew-Psi)
    if (max(abs(diff)<tol)),
        disp('convergence');
        break;
    else
        Psi = omega*PsiNew+(1-omega)*Psi
        disp(['iteration #', num2str(it)]);
        disp(['Current coefficients: ', num2str(Psi)]);
    end;
end;
if (it>=maxit),
    warning('no convergence');
end;
disp(['Equilibrium coefficients: ', num2str(PsiNew)]);




end











  
    

    
    
        
        
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
