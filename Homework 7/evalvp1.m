function vpp1 = evalvp1(cons,x)

global betta tetta r g gridx vpfun epsi probepsi ne nx 

vpp1 = zeros(ne,1);
for ec=1:ne,
    xp1 = ((x-cons)*(1+r)/(1+g)+epsi(ec))^(-1/tetta);
    vpp1(ec) = func_intp(gridx,vpfun,xp1);
end;
inv_vpp1 = vpp1.^(-tetta);
vpp1 = sum(inv_vpp1.*probepsi);

end


   
