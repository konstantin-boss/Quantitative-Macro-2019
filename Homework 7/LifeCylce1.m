%% Problem Set 7 - Exercise Part d) %%
% Converting the infinite horizon model into a life-cycle model %

clear all
close all

global betta tetta r g gridx vpfun epsi probepsi ne nx min_cons

% -------------------------------------------------------------------------
% SETTINGS
maxit = 100; 
tol=1e-4;
nt=1100;    % periods for simulation
dt=100;     % burn in periods
min_cons=1.0e-08;

% parameters
r = 0.02;
rho = 0.03;
g = 0.01;
tetta = 1;
betta = 1/(1+rho);
sr = readfile([],'MR.txt',3); % Surv rate
sr1 = sr(21:end,2);


% grid
nx=50;              % # of grid-points
curv=3.0;           % curvature of grid
xmax = 30;          % scaling factor of saving grid
xmin = sqrt(eps);
gridx=makegrid(xmin,xmax,nx,curv);
gridx=gridx';

% income shocks
ne = 7;
varepsi = 0.01;
muepsi = -varepsi/2;
% [epsi,probepsi] = qnwnorm(ne,muepsi,varepsi);
% mat=[[1:ne]',epsi,probepsi];
% save rn.txt mat -ascii -double -tabs
mat=load('rn.txt');
epsi=mat(:,2);
probepsi=mat(:,3);
epsi=exp(epsi);
if (abs(sum(epsi.*probepsi)-1.0)>sqrt(eps)),
    error('random numbers fucked up');
end;
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% SOLUTION 
% The solution mechanism follows Problem Set 6 - A recursion on the value
% function in the HH model implementing also a borrowing constraint
% Life cycle
param.t=80;
% consumption last period 
cons=zeros(nx,param.t);
cons(:,param.t)=gridx;

% finite period simulation
for it=param.t-1:-1:1;
    vpfun =margutil(cons(:,it+1));
    for xc=1:nx;
        % check binding constraint:
        mincons=gridx(xc);
        mu = foc(mincons,gridx(xc),sr1(xc));
        if (mu>=0.0),
            cons(xc,it)=mincons;
        else,
            [cons(xc,it),fval] = fzero('foc',cons(xc,it),[],gridx(xc));
            cons(xc,it)=max(cons(xc,it),min_cons);
        end ;
    end;
   
    end
    
  
% -------------------------------------------------------------------------

figure;
plot(gridx,cons(:,1),'LineWidth',2);
xlabel('x');
ylabel('c');
title('consumption policy');

%[gridx,cons]





