function [I2,nfeval] = adattivatrap( f, a, b, tol, fa, fb )
% [I2,nfeval] = adattivatrap( f, a, b, tol)
%
%Input
%f: Function handler contenente la funzione di cui vogliamo calcolare l'approssimazione dell'integrale definito, 
%a,b: Estremi dell'intervallo d'integrazione,  
%tol: Tolleranza richiesta 
%(fa,fb): parametri di lavoro da non specificare quando si vuole invocare
%la funzione
%
%Output
%If: approssimazione dell'integrale ottenuta,  
%nfeval: numero di valutazioni funzionali effettuate,
%
%Calcola un'approssimazione dell'integrale definito della funzione fun
%con estremi a e b utilizzando la formula adattiva di Newton-Cotes
%di grado 1 con tollerenza tol
% 
if(tol<=0), error('Tolleranza non valida!');end
nfeval=0;
if nargin<=4
fa=f(a); 
fb=f(b);
nfeval=2;
end 
nfeval=nfeval+1; %aggiorno il numero di valutazioni funzionali richieste
h=b-a; 
xm=(a+b)/2;  %calcolo il punto medio
fm=f(xm);
I1=(h/2)*(fa+fb); %Applico la formula dei trapezzi
I2=(I1+h*fm)/2;  %Applico la formula composita dei trapezzi su fa,fm e fb 
er=abs(I2-I1)/3; %mi ricavo una stima dell'errore
if er>tol %se la tolleranza non è rispettata riapplico la function sui due sottointervalli
    [lI,nfevall]=adattivatrap( f,a, xm, tol/2, fa, fm );
    [rI,nfevalr]=adattivatrap(f,xm, b, tol/2, fm, fb ); 
    I2=lI+rI;
    nfeval=nfeval+nfevalr+nfevall; %sommo le valutazioni funzionali richieste nei due sottointervalli
end
return
end


