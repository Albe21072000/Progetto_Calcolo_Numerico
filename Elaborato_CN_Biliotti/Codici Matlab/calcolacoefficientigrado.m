function coef=calcolacoefficientigrado(n)
%
% coef=calcolacoefficientigrado(n)
%
% Input
%n: Grado (maggiore di 0) della formula di Newton-Cotes di cui vogliamo conoscere i pesi
%della quadratura.
%
%Output
%coef: pesi della quadratura della formula di grado desiderato
%
%Calcola i pesi della quadratura della formula di quadratura di
%Newton-Cotes di grado n.
if(n<=0), error("Valore del grado della formula di Newton-Cotes non valido"),end
coef=zeros(n+1,1);
if (mod(n,2) == 0)
    for i=0:n/2-1
        coef(i+1)=calcolacoefficienti(i,n);
    end
    coef(n/2+1)=n-sum(coef)*2;
    coef((n/2)+1:n+1)=coef((n/2)+1:-1:1);
else
    for i=0:round(n/2,0)-2
        coef(i+1)=calcolacoefficienti(i,n);
    end
    coef(round(n/2,0))=(n-sum(coef)*2)/2;
    coef(round(n/2,0)+1:n+1)=coef(round(n/2,0):-1:1);
end
return
end

function cin=calcolacoefficienti(i,n)
%calcola il peso della quadratura della formula di Newton-Cotes numero i di
%grado n
d=i-[0:i-1 i+1:n];
den=prod(d);
a=poly([0:i-1 i+1:n]);
a=[a./((n+1):-1:1) 0];
num=polyval(a,n);
cin=num/den;
end