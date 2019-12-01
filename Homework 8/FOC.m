function f = FOC(cons,x,t,yc)

global betta tetta r g gridx vpfun epsi probepsi ne nx sr

vpp1 = evalvp2(cons,x,t,yc);
margu = MUc(cons);
f = margu - betta*sr(t)*(1.0+r) * vpp1;

end