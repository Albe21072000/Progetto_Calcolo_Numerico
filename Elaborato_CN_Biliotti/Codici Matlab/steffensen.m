function [xfin,passi]=steffensen(x0, func, tol, max)
%[xfin,passi]=steffensen(x0, func, tol, max)
%Input: x0 punto d'innesco del metodo, func function handler della funzione di cui si vuole
%conoscere la radice, (opzionali) tol tolleranza desiderat, max massimo
%numero di iterazioni.
%Output: xfin radice trovata dal metodo, passi passi necessari per trovare
%la soluzione.
%Cerca la radice della funzione passata in input con il metodo di
%Steffensen.
passi=0;
if ~exist("max", "var"), max=100; end %gestico il caso di input mancanti
if ~exist("tol", "var"), tol=eps; end
for i=0:max
    x1=x0;
    f0=func(x0);
    f1=func(x0+f0)-f0;
    if f1==0  %evito di dividere per 0
        disp("Non si può dividere per 0, radice non ulteriormente approssimabile!"), break
    end
    f=(f0^2)/f1; %passo iterativo del metodo di Steffensen
    xfin=x0-f;
    passi=passi+1;
    if(abs(xfin-x0)<tol*(1+abs(x0))) %controllo se la tolleranza è rispettata
        break;
    end
    x0=xfin;
end
if abs(xfin-x1)>tol*(1+abs(x0)), disp("Il metodo non è convergente!"), end