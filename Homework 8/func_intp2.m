function fv = func_intp2(x,func,xp)
        
       
        n = length(x);
        if ( xp>x(n) ),
            % fv = func(n);
            fv=func_extrapol2(x(n-1),x(n),func(n-1),func(n),xp);
        elseif (xp<x(1)),
            % fv = func(1);
            fv=func_extrapol2(x(1),x(2),func(1),func(2),xp);
        else
            fv = interp1(x,func,xp);
        end;
        
  end