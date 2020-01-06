function wage=func_firm(mpk)

global alpha delta

k = (alpha./mpk).^(1/(1-alpha));
wage = (1-alpha)*k.^alpha;

end