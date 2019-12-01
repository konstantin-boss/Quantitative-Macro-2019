% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function geneq

close all

global nj ny

tic

opt_det=false;          % 1=deterministic model
opt_nosr=false;         % 1=no survival risk
opt_ny = 2;             % 1=Markov chain with number of states, ny=5,
% 2=Markov chain with ny=2 (Krüger-Ludwig
% calibration)

% -------------------------------------------------------------------------
% SOLUTION

% calibration
func_calibr(opt_det,opt_nosr,opt_ny);

ret = 0.04;
tol = 1e-4;
maxit = 100;
df = 0.1;

for it=1:maxit,
    [fval, ass, Y, Phi, gridass, cfun, vfun, PhiAss] = func_olg(ret);
    if (abs(fval)<tol),
        disp('convergence');
        break;
    else,
        ret = ret - df * fval;
        disp(['iteration #', num2str(it)]);
        disp(['guess of rate of return: ', num2str(ret)]);
    end;
end;
if (it>=maxit),
    warning('no convergence');
end;
disp(['equilibrium interest rate: ', num2str(ret)]);
% disp(['equilibrium contribution rate: ', num2str(tau)]);
disp(['equilibrium capital output ratio: ', num2str(ass/Y)]);



% average life-cycle profiles
[labinclife,inclife,asslife,conslife,vallife] = lcprofile(Phi,gridass,cfun,vfun);
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% PLOTS
% graph of A(r) and K(r)
retold = ret;

% graph of A(r)
retvec = [-delta+0.01:0.001:0.1]';
n = length(retvec);
assvec = zeros(n,1);
for i = 1:n,
    ret = retvec(i);
    [fval, assvec(i)] = func_olg(ret);
end;

% graph of K(r):
kvec = (alpha./(retvec+delta)).^(1/(1-alpha)).*L;

figure; hold on;
plot(retvec,assvec,'b-','LineWidth',3);
plot(retvec,kvec,'r--','LineWidth',3);
legend('demand', 'supply');
title('supply and demand in OLG economy');
print -depsc equil.eps

end     % end function func_main
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function func_calibr(opt_det,opt_nosr,opt_ny)
close all
global betta tetta nj jr nx ny pi gridy netw pens sr epsi curv pini frac pop totpop grdfac delta
global alpha  replrate tau  
global retage maxage L R N_age  

rho = 0.04;  % 0.01
delta = 0.05;
alpha = 0.33; % 0.33; % 0.1;
ret = 0.04;


maxage = 80;    % corresponds to real age 100
retage = 45;    % corresponds to actual age 65
replrate = 0.0;
betta = 1/(1+rho);
tetta = 2;

nj=80;
jr=45;


nx=30;         % # of grid-points
curv=3.0;       % curvature of grid
grdfac=60;      % scaling factor of saving grid

% deterministic income component:
netw=1.0;
pens=0.4;
epsi=ones(nj,1);
if (jr<nj),
    epsi(jr+1:nj)=0.0;
end;

% survival rates
if opt_nosr,
    sr = ones(nj,1);
else
    mr = readfile([],'MR.txt',3);
    sr = 1.0-mr(21:21+nj-1,1);
end;


% population and fraction living in year...
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

% # of income states
if (opt_det==1),
    ny = 1;
    pini = 1.0;
    gridy = 1.0;
    pi = 1.0;
else
    
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
    end;
    
L=sum(pop(1:45));
R=1-L;
    
end;


end     % end function func_calibr
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++

% ----------------------------------------
function [fval, ass, Y, Phi, gridass, cfun, vfun, PhiAss] = func_olg(ret)

global retage maxage L R N_age replrate tau sr lambda delta gridx

mpk = ret + delta;

wage = func_firm(mpk);

tau = func_pens(L,R,replrate);

inc = func_inc(wage,tau,replrate,retage,maxage);

% solution of household model
[gridx,gridsav,gridass,cfun,vfun] = func_hh(ret, wage);

% aggregation
[Phi,PhiAss,ass] = func_aggr(gridx,gridsav,cfun,gridass, ret, wage);


[mpk,Y] = func_mpk(ass, L);

retnew = mpk - delta; 

fval = ret - retnew;

end
% ----------------------------------------
% -------------------------------------------
function wage=func_firm(mpk)

global alpha delta

k = (alpha/mpk)^(1/(1-alpha));
wage = (1-alpha)*k^alpha;

end
% -------------------------------------------
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function [gridx,gridsav,gridass,cfun,vfun] = func_hh(ret, wage)

global betta tetta r nj nx ny pi gridy netw pens sr epsi curv tau replrate grdfac

disp('solution of household model');

% grids and decisions rules:
gridx = zeros(nj,ny,nx);
gridsav = zeros(nx,1);
gridass = zeros(nj,ny,nx);
cfun = zeros(nj,ny,nx);
vfun = zeros(nj,ny,nx);
vpfun = zeros(nx,1);
vptrans = zeros(nj,ny,nx);

% savings grid: hold it constant:
maxsav=grdfac;
gridsav(2:nx)=makegrid(0.0,grdfac,nx-1,curv);
gridsav(1)=0.0;

% income states
for yc=1:ny
    % cash-on-hand grid at nj:
    inc = epsi(nj)*wage*gridy(yc)*(1-tau)+(1-epsi(nj))*wage*replrate*(1-tau);
    
    % in case of no pension system, assume some minimum cash on hand:
    minx=max(inc,sqrt(eps));
    maxx=gridsav(nx)*(1.0+r)+inc;
    gridx(nj,yc,:)=linspace(minx,maxx,nx);
    
    % Final period Consumption function, asset holdings, value function, including derivative
    cfun(nj,yc,:)=gridx(nj,yc,:);
    gridass(nj,yc,:)=(gridx(nj,yc,:)-inc)/(1+ret);
    vfun(nj,yc,:)=U(cfun(nj,yc,:));
    vpfun(:)=MUc(cfun(nj,yc,:));
    vptrans(nj,yc,:)=vpfun.^(-1.0/tetta);
end;

% Iterate Backwards
for jc=nj-1:-1:1,
    
    for yc=1:ny
        
        for xc=2:nx,
            vp=zeros(2,1);
            
            for ycc=1:ny,
                % income tomorrow:
                incp1=epsi(jc+1)*wage*gridy(ycc)*(1-tau)+(1-epsi(jc+1))*wage*replrate*(1-tau);
                
                % Maximum cash on hand tomorrow:
                % in case of zero savings and no pension system assume some
                % minimum cash on hand
                cah=max(sqrt(eps),incp1+(1.0+ret)*gridsav(xc));
                
                % Interpolate derivative of value function
                if ( cah<gridx(jc+1,ycc,1)),
                    disp('how can this be?')
                end;
                if ( cah>gridx(jc+1,ycc,nx) ),
                    % if out of bounds simply set it to decision at nx:
                    vptr = vptrans(jc+1,ycc,nx);
                else
                    vptr = interp1(squeeze(gridx(jc+1,ycc,:)),squeeze(vptrans(jc+1,ycc,:)),cah);
                end;
                vp(ycc)=vptr.^(-tetta);
            end;
            
            % Euler equation: RHS
            expvp=betta*sr(jc)*(1.0+ret)*sum(pi(yc,:)*vp(:));
            
            % consumption
            cfun(jc,yc,xc)=invut(expvp);
            
            % endogenous x-grid:
            gridx(jc,yc,xc)=gridsav(xc)+cfun(jc,yc,xc);
        end;
        
        % income (wages and pensions) in current period/age:
        inc=epsi(jc)*wage*gridy(yc)*(1-tau)+(1-epsi(jc))*wage*replrate*(1-tau);
        
        % decision at minx
        % notice: correction required for welfare calculation
        % the above is actually slightly inefficient because xmin
        % can be explicitly computed, then gridsav would be age and
        % state dependent.
        minx=max(inc,sqrt(eps));
        if (minx<gridx(jc,yc,2)),
            gridx(jc,yc,1)=minx;
        else    % set it to some arbitrary fracion of x(2)
            gridx(jc,yc,1)=0.9*gridx(jc,yc,2);
        end;
        
        % Compute optimal consumption and leisure for minx
        cfun(jc,yc,1)=gridx(jc,yc,1);
        
        % assets at all xc:
        gridass(jc,yc,:)=(gridx(jc,yc,:)-inc)/(1+ret);
        
        % Update vfun and vpfun
        vpfun(:)=MUc(cfun(jc,yc,:));
        vptrans(jc,yc,:)=vpfun(:).^(-1.0/tetta);
        
        % Calculate value function
        for xc=1:nx,
            
            v=zeros(2,1);
            for ycc=1:ny,
                % income tomorrow:
                incp1=epsi(jc+1)*wage*gridy(ycc)*(1-tau)+(1-epsi(jc+1))*wage*replrate*(1-tau);
                
                % cah tomorrow
                cah=max(sqrt(eps),incp1+(1.0+ret)*gridsav(xc));
                
                % this should never be the case:
                if ((cah+0.0001)<gridx(jc+1,ycc,1)),
                    warning('How can this be ?');
                end;
                % linear interpolation:
                v(ycc)=func_intp(squeeze(gridx(jc+1,ycc,:)),squeeze(vfun(jc+1,ycc,:)),cah);
            end;    % end for ycc
            
            % update value function
            expv=sum(pi(yc,:)*v(:));
            vfun(jc,yc,xc)=U(cfun(jc,yc,xc))+betta*sr(jc)*expv;
        end;    % end for xc
        
    end;    % end for yc
    
end;    % end for jc


% ---------------------------------------------------------------------
    function fv = func_intp(x,func,xp)
        
        
        n = length(x);
        if ( xp>x(n) ),
            % fv = func(n);
            fv=func_extrapol(x(n-1),x(n),func(n-1),func(n),xp);
        elseif (xp<x(1)),
            % fv = func(1);
            fv=func_extrapol(x(1),x(2),func(1),func(2),xp);
        else
            fv = interp1(x,func,xp);
        end;
        
    end
% ---------------------------------------------------------------------


% ---------------------------------------------------------------------
    function y=func_extrapol(x1,x2,y1,y2,x)
        
        % simple linear extrapolation
        
        m = (y2-y1)/(x2-x1);
        y = y1 + m*(x-x1);
        
    end
% ---------------------------------------------------------------------

end     % end function func_hh
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function [Phi,PhiAss,ass]=func_aggr(gridx,gridsav,cfun,gridass, ret, wage)

global r nj nx ny pi gridy netw pens sr epsi pini frac totpop replrate tau

disp('aggregation and cross-sectional measure');

% Compute Cross sectional distributions and aggregate variables
Phi = zeros(nj,ny,nx);          % distribution of assets conditional by age and shock
PhiAss = zeros(nx,1);             % distribution of assets

% Distribution of newborns over cash at hand
for yc=1:ny
    
    % income (wages and pensions) in current period/age:
    inc=epsi(1)*wage*gridy(yc)*(1-tau)+(1-epsi(1))*wage*replrate*(1-tau);
    
    % initial cash-on-hand:
    cahini=inc;
    
    [vals,inds]=basefun(gridx(1,yc,:),cahini,nx);
    Phi(1,yc,inds(1))=vals(1)*pini(yc)*frac(1);
    Phi(1,yc,inds(2))=vals(2)*pini(yc)*frac(1);
end;

for jc=2:nj
    TT = zeros(ny,nx,ny,nx);    % transfer function
    
    for xc=1:nx
        for yc=1:ny
            for ycc=1:ny
                
                % income (wages and pensions) in current period/age:
                inc=epsi(jc)*wage*gridy(ycc)*(1-tau)+(1-epsi(jc))*wage*replrate*(1-tau);
                
                % cash on hand: x=a*(1+r)+y = s(-1)*(1+r)+y;
                cah=inc+(1.0+ret)*gridsav(xc);
                
                [vals,inds]=basefun(gridx(jc,ycc,:),cah,nx);
             
                TT(ycc,inds(1),yc,xc)=vals(1)*pi(yc,ycc);
                TT(ycc,inds(2),yc,xc)=vals(2)*pi(yc,ycc);
                end;
             for xcc=1:nx
                for yccc = 1:ny
%                          if (Phi(jc-1,ycc,xcc) == 0) && (jc>2)   This
%                          somehow does not work
%                              break;
%                          end;
                   Phi(jc,yccc,xcc)=Phi(jc,yccc,xcc)+Phi(jc-1,yc,xc)*TT(yccc,xcc,yc,xc)*sr(jc-1);
    
            end; end;    
        end;    
    end;    

end;    % end for jc


% Check that for each country distribution sums to 1
sumprob=sum(sum(sum(Phi(:,:,:))));
if ( ( sumprob < 0.999 ) || ( sumprob > 1.001) ),
    beep; beep; beep;
    warning('distribution fucked up');
end;

% Check if Grid is Big enough
sumprob=sum(sum(Phi(:,:,nx)));
if (sumprob > 0.001 ),
    beep; beep; beep;
    warning('grid too small -- increase your grid');
    pause
end;

ass=0.0;
cons=0.0;
% lab=0.0;
% ret=0.0;

% aggregation
for jc=1:nj
    for yc=1:ny
        for xc=1:nx,
            PhiAss(xc)=PhiAss(xc)+Phi(jc,yc,xc);
            
            % asset holdings = capital stock in general equilibrium
            ass=ass+totpop*Phi(jc,yc,xc)*gridsav(xc);
            
            cons=cons+totpop*Phi(jc,yc,xc)*cfun(jc,yc,xc);
            
%             lab=lab+totpop*Phi(jc,yc,xc)*gridy(yc)*epsi(jc);
%             ret=ret+totpop*Phi(jc,yc,xc)*gridy(yc)*(1.0-epsi(jc));
        end;
    end;
end;


% ---------------------------------------------------------------------
    function [vals,inds]=basefun(grid_x,x,nx)
        % this subroutine returns the values and the indices of the two basis
        % functions that are positive on a given x in the grid_x
        
        % MF function to lookup the current position
        i=lookup1(grid_x,x,0);
        
        if ( (i+1)>nx),
            inds(1)=nx;
            inds(2)=nx;
            vals(2)=0.0;
            vals(1)=1.0;
        elseif (i==0),
            inds(1)=1;
            inds(2)=1;
            vals(1)=1.0;
            vals(2)=0.0;
        else
            inds(1)=i;
            inds(2)=i+1;
            dist = grid_x(i+1)-grid_x(i);
            vals(2)=( x-grid_x(i) )/dist;
            vals(1)=( grid_x(i+1)-x )/dist;
        end;
        
    end 	% end function basefun
% ---------------------------------------------------------------------




end     % end function func_aggr
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function [labinclife,inclife,asslife,conslife,vallife] = lcprofile(Phi,gridass,cfun,vfun)

global r nj nx ny gridy netw pens epsi frac

disp('life-cycle profiles')

asslife=zeros(nj,1);
inclife=zeros(nj,1);
labinclife=zeros(nj,1);
conslife=zeros(nj,1);
vallife=zeros(nj,1);

for jc=1:nj,
    for yc=1:ny
        for xc=1:nx,
            asslife(jc)=asslife(jc)+Phi(jc,yc,xc)*gridass(jc,yc,xc)/frac(jc);
            conslife(jc)=conslife(jc)+Phi(jc,yc,xc)*cfun(jc,yc,xc)/frac(jc);
            
            inc=epsi(jc)*netw.*gridy(yc)+(1-epsi(jc))*pens;
            labinclife(jc)=labinclife(jc)+Phi(jc,yc,xc)*inc/frac(jc);
            inclife(jc)=inclife(jc)+Phi(jc,yc,xc)*(r*gridass(jc,yc,xc)+inc)/frac(jc);
            
            vallife(jc)=vallife(jc)+Phi(jc,yc,xc)*vfun(jc,yc,xc)/frac(jc);
        end;
    end;
end;

end     % end function lcprofile
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++



% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function [pini,pi,gridy]=mchain(rhoeta,epsil)

% Transition Probabilities
pi=rhoeta*ones(2,2);
pi(1,2)=1.0-rhoeta;
pi(2,1)=1.0-rhoeta;

% Initial Distribution
pini=0.5*ones(2,1);

gridy=zeros(2,1);
gridy(1)=exp(1.0-epsil);
gridy(2)=exp(1.0+epsil);
gridy=2.0*gridy/(sum(gridy));

end  % end function mchain
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++



% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function u = U(c)
global tetta

% utility
if (abs(tetta-1-0)<sqrt(eps)),
    u = log(c);
else
    u = c.^(1.0-tetta)/(1.0-tetta);
end;
end     % end function U
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function muc=MUc(c)
global tetta

% maringal utility
muc = c.^(-tetta);
end     % end function MUc
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++

% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function invut=invut(marg)
global tetta

% invert utility for c
invut=marg.^(-1.0/tetta);
end     % end function invut
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++

% ---------------------------------
function [N,L,R]=func_pop(sr,maxage,retage);

N = zeros(maxage,1);
L = zeros(maxage,1);
R = zeros(maxage,1);
N(1) = 100;
for age=2:maxage,
    N(age)=N(age-1)*sr(age-1);
end;
L = N(1:retage);
R = N(retage+1:maxage);

end
% ---------------------------------


% ----------------------------------------------
function tau = func_pens(L,R,replrate)

tau = replrate*R ./ (L + replrate * R);

end
% ----------------------------------------------

% ---------------------------------------------
function inc=func_inc(wage,tau,replrate,retage,maxage);

inc=zeros(maxage,1);
inc(1:retage)=wage*(1-tau);
inc(retage+1:maxage)=replrate*wage*(1-tau);

end
% ---------------------------------------------


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function grd = makegrid(x1,x2,n,c)
% makes curved grid according to curvature parameter c
scale=x2-x1;
grd(1)=x1;
grd(n)=x2;
for i=2:n-1
    grd(i)=x1+scale*((i-1.0)/(n-1.0))^c;
end;
end     % end function makegrid
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++


% ----------------------------------------
function [mpk,Y] = func_mpk(ass, L)

global alpha

Y = ass.^alpha * L.^(1-alpha);
ky = ass./Y;
mpk = alpha * ky.^(-1);

end
% ----------------------------------------