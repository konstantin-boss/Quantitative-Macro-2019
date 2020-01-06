alpha = 0.33;
z_state = [1.03 0.97];
kstar = 7;
rho = 0.04;
betta = 1/(1+rho);
grdfac = 60;
curv = 0.3;
Psi = [1, alpha; 0.9 alpha];
nj = 80;
mr = readfile([],'MR.txt',3);
sr = 1.0-mr(21:21+nj-1,1);
jr =45;
delta = 0.05;
nk = 5;
nz = 2;
kgrid = linspace(0.5*kstar, 1.5*kstar,nk);
kgrid1 = zeros(nz,nk);
ret = zeros(nz,nk);
wage = zeros(nz,nk);
ret1 = zeros(nz,nk);
wage1 = zeros(nz,nk);
Psi =[1.3 0.3;1.2 0.3];
netw=1.0;
tetta = 2;
nx = 30;
pens=0.4;%0;%
replrate=0.6;%0;%
epsi=ones(nj,1);
if (jr<nj),
    epsi(jr+1:nj)=0.0;
end;
pop=zeros(nj,1);
pop(1)=100;
for jc=2:nj,
    pop(jc)=pop(jc-1)*sr(jc-1);
end;
totpop=sum(pop);

% normalize population to one:
pop=pop/totpop;
totpop=1.0;
frac=pop./totpop;
L=sum(pop(1:45));
R=1-L;

pi_z = [0.95 0.05; 0.05 0.95];
pi_state = [1.03, 0.97];
nz = 2;
opt_ny = 2;
    
    if (opt_ny==1)
        % number of income shocks
        ny = 5;
        % transition probability
        rhoeta=0.98;
        % variance of "permanent" shock
        % taken from Campbell, Viceira, Ch. 7
        vareta=0.01;
        % Markov chain:
        [pi,gridy] = markovappr(rhoeta,sqrt(vareta),2,ny);
        
        pi_ez = kron(pi_z, pi);
        % compute invariant distribution
        pini = 1/ny*ones(ny,1);
        for tc=1:100,
            pini = pi'*pini;
        end;
        
        % take exponent and rescale income shocks such that mean is one:
        gridy=exp(gridy)';
        gridy = gridy / sum(gridy.*pini);
        
    else
        
        % Alternative -- taken from Krüger and Ludwig (2007) using only two
        % states
        ny = 2;
        
        % transition probability and variance
        rhoeta=0.97;
        vary=0.08;    % taken from Storesletten, Telmer, Yaron
        
        % shock
        epsil=sqrt(vary/(4.0*rhoeta*(1.0-rhoeta)));
        
        % Markov chain
        [pini,pi,gridy]=mchain(rhoeta,epsil);
        pi_ez = kron(pi_z, pi);
    end;
