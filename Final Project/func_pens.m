function tau = func_pens(L,R,replrate)

tau = replrate*R ./ (L + replrate * R);

end