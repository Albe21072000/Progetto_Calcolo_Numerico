function [xfin,passi]=secanti(x0, x1, func, tol, max)
%[xfin, passi]=secanti(x0, x1, func, tol, max)
%Output
%xfin: radice ottenuta, passi: passi effetutati dal metodo prima di raggiungere la tolleranza desiderata.
%Input
%x0: punto d'innesco del metodo, func: funzione di cui si vuole conoscere la radice,
%x1: punto calcolato con un iterazione del metodo di Newton da x0, (opzionali) tol: tolleranza desiderata(impostata di default sulla precisone di macchina) , max=
%numero massimo di passi effetuati dal metodo(di default 1000)
%Calcola la radice della funzione passata in input con il metodo delle
%secanti usando come punto d'innesco il valore passato in input
passi=0;
if ~exist("max", "var"), max=1000; end %gestico il caso di input mancanti
if ~exist("tol", "var"), tol=eps; end 
f=func(x0);
xfin=x1;
for i=1:max
    passi=passi+1;
    if(abs(xfin-x0)<tol*(1+abs(x0))) %controllo se la tolleranza è rispettata
        break;
    end
    f0=f;
    f=func(xfin);
    if(f0==f), disp("Precisione massima raggiunta impossibile proseguire!"),break,end
    x1=(f*x0-f0*xfin)/(f-f0); %passo iterativo del metodo delle Secanti
    x0=xfin;
    xfin=x1;
end
if abs(xfin-x1)>tol*(1+abs(x0)), disp("il metodo non è convergente!"), end