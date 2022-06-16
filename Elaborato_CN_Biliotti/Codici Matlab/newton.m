function [xfin, passi]=newton(x0, func, der, tol, max)
%[xfin, passi]=newton(x0, func, der, tol, max)
%Output
%xfin: approsimazione della radice ottenuta, passi: passi effetutati dal metodo prima di raggiungere la tolleranza desiderata.
%Input
%x0: punto di partenza del metodo, func: function handler contenente una funzione(sufficientemente regolare) di cui si vuole conoscere la radice,
%der: derivata prima di tale funzione, (opzionali) tol: tolleranza desiderata(impostata di default sulla precisone di macchina) , max=
%numero massimo di passi effetuati dal metodo(di default 1000)
%
%Calcola un'approssimazione della radice della funzione passata in input con il metodo di Newton
%usando come punto d'innesco il valore passato in input
passi=0;
if ~exist("max", "var"), max=1000; end %gestico il caso di input mancanti
if ~exist("tol", "var"), tol=eps; end 
for i=1:max
    x1=x0;
    f0=func(x0);
    d0=der(x0);
    if(d0==0), disp("Derivata prima uguale a 0, usare un nuovo punto d'innesco o una funzione più regolare!"), xfin=x0; break,end
    passi=passi+1;
    xfin=x0-f0/d0;  %passo iterativo del metodo di Newton
    if(abs(xfin-x0)<tol*(1+abs(x0))) %controllo se la tolleranza è rispettata
        break;
    end
    x0=xfin;
end
if abs(xfin-x1)>tol*(1+abs(x0)), disp("il metodo non è convergente!"), end

