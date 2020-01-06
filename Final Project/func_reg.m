        KKbad = [];
        KKbad1 = [];
        KKgood = [];
        KKgood1 = [];
        
        for jlo = 1:T
            if chain(jlo) == 1
                KKgood(jlo) = KK(jlo);
                               
            else
                KKbad(jlo) = KK(jlo);
                
            end;

        end;

        KKgood1 = KKgood(2:T)';
        KKgood = KKgood(1:T-1)';
        KKbad1 = KKbad(2:T-1)';
        KKbad = KKbad(1:T-2)';
%   

    XXg = [ones(length(KKgood),1), KKgood];
    yyg = KKgood1;
    XXb = [ones(length(KKbad),1), KKbad];
    yyb = KKbad1;
    
    PsiG = (inv(XXg'*XXg)*XXg'*yyg)';
    PsiB = (inv(XXb'*XXb)*XXb'*yyb)';
    PsiNew = [PsiG; PsiB];
