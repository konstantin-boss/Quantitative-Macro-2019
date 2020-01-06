disp('solution of household model');

% grids and decisions rules:
gridx = zeros(nj,ny,nx,nz,nk);
gridsav = zeros(nx,nz);
gridass = zeros(nj,ny,nx,nz,nk);
cfun = zeros(nj,ny,nx,nz,nk);
vfun = zeros(nj,ny,nx,nz,nk);
vpfun = zeros(nx,nz);
vptrans = zeros(nj,ny,nx,nz,nk);

k1grid = zeros(nz,nk);
ret = zeros(nz,nk);
wage = zeros(nz,nk);
ret1 = zeros(nz,nk);
wage1 = zeros(nz,nk);

for i=1:nk
    for j=1:nz
        ret(j,i) = alpha*z_state(j)*kgrid(i)^(alpha-1)*L^(1-alpha)-delta;
        wage(j,i)= (1-alpha)*z_state(j)*kgrid(i)^(alpha)*L^(-alpha);
        kgrid1(j,i) = exp(Psi(j,1)+Psi(j,2)*log(kgrid(i)));
        ret1(j,i) = alpha*z_state(j)*kgrid1(j,i)^(alpha-1)*L^(1-alpha)-delta;    
        wage1(j,i)=(1-alpha)*z_state(j)*kgrid1(j,i)^(alpha)*L^(-alpha);
    end;
end;

tau = func_pens(L,R,replrate);

% savings grid: hold it constant:
maxsav=grdfac;
gridsav(1:nx,1)=makegrid(1,grdfac,nx,curv);
gridsav(1:nx,2)=makegrid(1,grdfac,nx,curv);


% income states
for yc=1:ny
    for zz=1:nz
        for kk=1:nk
    % cash-on-hand grid at nj:
        inc = epsi(nj)*wage(zz,kk)*gridy(yc)*(1-tau)+(1-epsi(nj))*wage(zz,kk)*replrate*(1-tau);
    
        % in case of no pension system, assume some minimum cash on hand:
        minx=max(inc,sqrt(eps));
        maxx=gridsav(nx,1)*(1.0+ret(zz,kk))+inc;
        gridx(nj,yc,:,zz,kk)=linspace(minx,maxx,nx);

        % Final period Consumption function, asset holdings, value function, including derivative
        cfun(nj,yc,:,zz,kk)=gridx(nj,yc,:,zz,kk);
        gridass(nj,yc,:,zz,kk)=(gridx(nj,yc,:,nz,kk)-inc)/(1+ret(zz,kk));
        vfun(nj,yc,:,zz,kk)=U(cfun(nj,yc,:,zz,kk));
        vpfun(:,zz)=MUc(cfun(nj,yc,:,zz,kk));
        vptrans(nj,yc,:,zz,kk)=vpfun(:,zz).^(-1.0/tetta);
    end;
  end;
end;

%%
% Iterate Backwards
for jc=nj-1:-1:1,
    
    for yc=1:ny
        for zzz=1:nz
              for kkk=1:nk
                   for xc=2:nx,
                                          
            vp=zeros(nz,1);
            for ycc=1:ny,
               
                % income tomorrow:
                incp1=epsi(jc+1)*wage1(zzz,kkk)*gridy(ycc)*(1-tau)+(1-epsi(jc+1))*wage1(zzz,kkk)*replrate*(1-tau);
                
                % Maximum cash on hand tomorrow:
                % in case of zero savings and no pension system assume some
                % minimum cash on hand
                cah=max(sqrt(eps),incp1+(1.0+ret1(zzz,kkk))*gridsav(xc,zzz));
                
                % Interpolate derivative of value function
                if ( cah<gridx(jc+1,ycc,1,zzz,kkk)),
                    disp('how can this be 11111?')
                    display(jc)
                end;
                if ( cah>gridx(jc+1,ycc,nx,zzz,kkk) ),
                    % if out of bounds simply set it to decision at nx:
                    vptr = vptrans(jc+1,ycc,nx,zzz,kkk);
                else
                    vptr = interp1(squeeze(gridx(jc+1,ycc,:,zzz,kkk)),squeeze(vptrans(jc+1,ycc,:,zzz,kkk)),cah);
                end;
                vp(ycc)=vptr.^(-tetta);
            end; % end for ycc
                          
                             
          
            % Euler equation: RHS
            expvp=betta*sr(jc)*(1.0+ret(zzz,kkk))*sum(sum(pi_z(zzz,:)*(pi(yc,:)*vp(:)))); %%%%%%%% Need yc to match 10x10!!!s
            
            % consumption
            cfun(jc,yc,xc,zzz,kkk)=invut(expvp);
            
            % endogenous x-grid:
            gridx(jc,yc,xc,zzz,kkk)=gridsav(xc,zzz)+cfun(jc,yc,xc,zzz,kkk);
        end; % end for xc = 2:
        
        
        % income (wages and pensions) in current period/age:
        inc=epsi(jc)*wage(zzz,kkk)*gridy(yc)*(1-tau)+(1-epsi(jc))*wage(zzz,kkk)*replrate*(1-tau);
        
        % decision at minx
        % notice: correction required for welfare calculation
        % the above is actually slightly inefficient because xmin
        % can be explicitly computed, then gridsav would be age and
        % state dependent.
        minx=max(inc,sqrt(eps));
        if (minx<gridx(jc,yc,2,zzz,kkk)),
            gridx(jc,yc,1,zzz,kkk)=minx;
        else    % set it to some arbitrary fracion of x(2)
            gridx(jc,yc,1,zzz,kkk)=0.9*gridx(jc,yc,2,zzz,kkk);
        end;
        
        % Compute optimal consumption and leisure for minx
        cfun(jc,yc,1,zzz,kkk)=gridx(jc,yc,1,zzz,kkk);
        
        % assets at all xc:
        gridass(jc,yc,:,zzz,kkk)=(gridx(jc,yc,:,zzz,kkk)-inc)/(1+ret(zzz,kkk));
        
        % Update vfun and vpfun
        vpfun(:,zzz)=MUc(cfun(jc,yc,:,zzz,kkk));
        vptrans(jc,yc,:,zzz,kkk)=vpfun(:,zzz).^(-1.0/tetta);
        
        % Calculate value function
        for xc=1:nx,
            
            v=zeros(nz*ny,1);
            for ycc=1:ny,
                % income tomorrow:
                incp1=epsi(jc+1)*wage(zzz,kkk)*gridy(ycc)*(1-tau)+(1-epsi(jc+1))*wage(zzz,kkk)*replrate*(1-tau);
                
                % cah tomorrow
                cah=max(sqrt(eps),incp1+(1.0+ret1(zzz,kkk))*gridsav(xc,zzz));
                
                % this should never be the case:
                if ((cah+0.0001)<gridx(jc+1,ycc,1,zzz,kkk)),
                    warning('How can this be 22222?');
                end;
                % linear interpolation:
                v(ycc*zzz)=func_intp(squeeze(gridx(jc+1,ycc,:,zzz,kkk)),squeeze(vfun(jc+1,ycc,:,zzz,kkk)),cah);
            end;    % end for ycc
            
            % update value function
            expv=sum(pi_ez(yc*zzz,:)*v(:));
            vfun(jc,yc,xc,zzz,kkk)=U(cfun(jc,yc,xc,zzz,kkk))+betta*sr(jc)*expv;
        end;    % end for xc
              end;
        end;
    end;    % end for yc
    
end;    % end for jc
