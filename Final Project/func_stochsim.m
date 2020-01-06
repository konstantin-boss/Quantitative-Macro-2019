
 % Simulation of Z Markov Chain
rng(1);
T = 100; %Number of stoachastic simulations
transition_probabilities = [0.95 0.05;0.05 0.95]; starting_value = 1; chain_length = T;
    chain = zeros(1,chain_length);
    chain(1)=starting_value;
    for i=2:chain_length
        this_step_distribution = transition_probabilities(chain(i-1),:);
        cumulative_distribution = cumsum(this_step_distribution);
        r = rand();
        chain(i) = find(cumulative_distribution>r,1);
    end

z_sim = [];
for i=1:T
    z_sim(i)=z_state(chain(i));
end;


% Simulation of aggregate capital
%kstoch = 3;
KK = zeros(T+1,1);
KK(1) = 3;
KK1= zeros(T+1,1);

for j=1:T
    gridx_sim = squeeze(gridx(:,:,3,chain(j),:));
    %cfun_sim = squeeze(cfun(:,:,3,chain(j),:));

    % Compute Cross sectional distributions and aggregate variables
Phi = zeros(nj,ny,nx,nz,nk);          % distribution of assets conditional by age and shock
PhiAss = zeros(nx,nz);             % distribution of assets

% Distribution of newborns over cash at hand
for yc=1:ny
    for zz=1:nz
        for kk=1:nk
    
    % income (wages and pensions) in current period/age:
    inc=epsi(1)*(1-alpha)*KK(j)^(alpha)*z_sim(j)*gridy(yc)*(1-tau)+(1-epsi(1))*(1-alpha)*KK(j)^(alpha)*z_sim(j)*replrate*(1-tau);
    
    % initial cash-on-hand:
    cahini=inc;
    
     [vals,inds]=basefun(gridx(1,yc,:,zz,kk),cahini,nx);
     Phi(1,yc,inds(1),zz,kk)=vals(1)*pini(yc)*frac(1);
     Phi(1,yc,inds(2),zz,kk)=vals(2)*pini(yc)*frac(1);
end; end; end;


for jc=2:nj
    TT = zeros(ny,nx,ny,nx,nz,nk);    % transfer function
    
    for xc=1:nx
        for zzz=1:nz
            for kkk=1:nk
        for yc=1:ny
              for ycc=1:ny
              
                % income (wages and pensions) in current period/age:
                inc=epsi(jc)*(1-alpha)*KK(j)^(alpha)*z_sim(j)*gridy(ycc)*(1-tau)+(1-epsi(jc))*(1-alpha)*KK(j)^(alpha)*z_sim(j)*replrate*(1-tau);
                
                % cash on hand: x=a*(1+r)+y = s(-1)*(1+r)+y;
                cah=inc+(1.0+alpha*KK(j)^(alpha-1)*z_sim(j)-delta)*gridsav(xc,zzz);
                
                [vals,inds]=basefun(gridx(jc,ycc,:,zzz,kkk),cah,nx);
                
                
                TT(ycc,inds(1),yc,xc,zzz,kkk)=vals(1)*pi_z(zzz,ycc)*pi(yc,ycc);
                TT(ycc,inds(2),yc,xc,zzz,kkk)=vals(2)*pi_z(zzz,ycc)*pi(yc,ycc);
%                 
            end;
            
        end;    
    end;    end;end;
  
    for xc=1:nx
        for yc=1:ny
                for zzzz=1:nz
            for kkkk=1:nk
                for xcc=1:nx
                    for ycc=1:ny
                    % transfer distribution:
                    Phi(jc,ycc,xcc,zzzz,kkkk)=Phi(jc,ycc,xcc,zzzz,kkkk)+Phi(jc-1,yc,xc,zzzz,kkkk)*TT(ycc,xcc,yc,xc,zzzz,kkkk)*sr(jc-1);
                end;
            end;
        end;end; end;
    end;
    
end;    % end for jc


ass=0.0;
cons=0.0;


% aggregation
for jc=1:nj
    for yc=1:ny
        for xc=1:nx,
            for zzz=1:nz
                for kkk=1:nk
            PhiAss(xc,zzz)=PhiAss(xc,zzz)+Phi(jc,yc,xc,zzz,kkk);
            
            % asset holdings = capital stock in general equilibrium
            ass=ass+totpop*Phi(jc,yc,xc,zzz,kkk)*gridsav(xc,zzz);
            
            cons=cons+totpop*Phi(jc,yc,xc,zzz,kkk)*cfun(jc,yc,xc,zzz,kkk);
            
        
           
        end;
    end;
end;
    end; end;

 KK(j+1) = ass;
end;
 