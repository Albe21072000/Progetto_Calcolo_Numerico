function [x,nit] = newton2(fun,jacobian,x0,tol,maxit)
%[x,nit] = newton2(fun,jacobian,x0,tol,max)
%Input: fun vettore di handler di funzioni nonlineari di cui si vuole trovare le radici,
%jacobian: matrice di handler di funzioni contenenti il giacobiano di fun, x0 vettore d'innesco del metodo,
%(opzionali) tol tolleranza desiderata(di default precisione di macchina),
%max massimo numero di iterazioni del metodo(di default 1000)
%Output
%x: Vettore contente le radici del sistema non lineare, nit numero di
%iterazioni compiute dal metodo per trovare la radice
%Trova il vettore radice del sistema di equzioni non lineari passati in
%input col metodo di Newton partendo dal vettore d'innesco passato in input
nit=0;
if ~exist("maxit", "var"), maxit=1000; end %gestico il caso di input mancanti
if ~exist("tol", "var"), tol=eps; end 
x=x0;
for i=1:maxit
    f=-fun(x);
    j=jacobian(x);
    dx=mialu(j,f); %uso la function mialu vista in precedenza per ricavarmi il vettore della differenza dx
    x=x+dx;
    nit=nit+1;
    if(norm(dx)<=tol*(1+norm(x-dx)))   %controllo se la tolleranza è rispettata
        break;
    end
end
if (norm(dx)>tol*(1+norm(x-dx))), disp('il metodo non è convergente!'), end
