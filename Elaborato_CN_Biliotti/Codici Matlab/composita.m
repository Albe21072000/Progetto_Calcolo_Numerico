function [If,err,nfeval] = composita(fun,a,b,n,tol)
% [If,err,nfeval] = composita( fun, a, b, n, tol)
%
%Input
%fun: Function handler contenente la funzione di cui vogliamo calcolare l'approssimazione dell'integrale definito, 
%a,b: Estremi dell'intervallo d'integrazione, 
%n: Grado della formula di Newton-Cotes desiderato, 
%tol: Tolleranza richiesta
%Output
%If: approssimazione dell'integrale ottenuta, 
%err: stima dell'errore dell'approssimazione, 
%nfeval: numero di valutazioni funzionali effettuate,
%
%Calcola un'approssimazione dell'integrale definito della funzione fun
%con estremi a e b utilizzando la formula composita di Newton-Cotes
%di grado n con tollerenza tol
% 
if(n<=0), error('Grado della formula di Newton-Cotes non valido!'),end
if(tol<=0), error('Tolleranza non valida!'),end
x=linspace(a,b,n+1)';
coef=calcolacoefficientigrado(n);  %function vista nell'esercizio 24
y=fun(x);
nfeval=n+1;
u=1;
if (mod(n,2) == 0)
    u=2;
end

int=integrale(a,b,y,n,coef);
dn=2*(n+1);
int2=0;
for i=1:1000
        xi=linspace(a,b,dn-1)';
        y2=zeros(dn-1,1);
        y2(1:2:dn)=y;
        y2(2:2:dn-1)=(fun(xi((2:2:dn-1))));
        nfeval=nfeval+(dn/2-1); %aggiorno le valutazioni funzionali richieste
        inf=1;sup=n+1;
        for j=1:1:2^i
            int2=int2+integrale((j-1)*(b-a)/(2^i),(j)*(b-a)/(2^i),y2(inf:sup),n,coef); %sommo tutti i sottointervalli dell'integrale approssimato con le formule di Newton-Cotes
            inf=sup;
            sup=sup+n;
        end
        err=(int2-int)/(2^(n+u)-1); %calcolo una stima dell'errore
        if(abs(err)<tol), If=int2;break,end %se la stima rispetta la tolleranza interrompo il ciclo
        y=y2;
        dn=size(y,1)*2;
        int=int2;
        int2=0;
end
return
end

function int =integrale(a,b,y,n,coef)
%Calcola la formula di Newton-Cotes tra a e b di grado n con coefficienti
%coef e con y contenente i valori che la funzione assume nelle n+1 ascisse equistanti nell intervallo [a,b]
int=(b-a)/n*sum(y.*coef);       
return;
end