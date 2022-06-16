function rad=radice(x)
%rad=radice(x)
%Input: x, numero positivo di cui si vuole conoscere la radice
%Output: rad, radice quadrata del numero in input
%Restituisce la radice quadrata del numero x passato in input.
if x==0
    rad=0;
    return;
end
if x<0, error("La radice quadrata di un numero negativo non è un numero reale!"), end
x0=x;
for i=0:10000
    rad=0.5*(x+x0/x);   %passo iterativo
    if(abs(rad-x)<eps*(1+abs(x))) %controllo di aver raggiunto la massima precisione possibile
        break;
    end
    x=rad;
end
return;