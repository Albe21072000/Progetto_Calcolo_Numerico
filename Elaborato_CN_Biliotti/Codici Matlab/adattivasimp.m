function [I2,nfeval] = adattivasimp( f, a, b, tol, fa, fm,fb )
% [I2,nfeval] = adattivasimp( f, a, b, tol)
%
%Input
%f: Function handler contenente la funzione di cui vogliamo calcolare l'approssimazione dell'integrale definito, 
%a,b: Estremi dell'intervallo d'integrazione,  
%tol: Tolleranza richiesta 
%(fa,fm,fb): parametri di lavoro da non specificare quando si vuole invocare
%la funzione
%
%Output
%If: approssimazione dell'integrale ottenuta,  
%nfeval: numero di valutazioni funzionali effettuate,
%
%Calcola un'approssimazione dell'integrale definito della funzione f
%con estremi a e b utilizzando la formula adattiva di Newton-Cotes
%di grado 2 con tollerenza tol
% 
if(tol<=0), error('Tolleranza non valida!'),end
m = ( a + b )/2;  %calcolo il punto medio
nfeval=0;
if nargin<=4
fa = f(a); 
fb = f(b);   
fm = f(m);
nfeval=3;
end 
nfeval=nfeval+2; %aggiorno il numero di valutazioni funzionali richieste
ma=(a+m)/2; %calcolo il punto medio dell'intervallo [a,m]
mb=(b+m)/2;  %calcolo il punto medio dell'intervallo [m,b]
fma=f(ma);
fmb=f(mb);
h = b - a;
I1= ( h/6 )*( fa +4*fm+ fb); %Applico la formula di Simpson su fa,fm ed fb 
I2 = I1/2 + (h/6)*(2*fma + 2*fmb-fm); %Applico la formula composita di Simpson su fa,fma,fm,fmb e fb 
e = abs( I2 - I1 )/15;  %mi ricavo una stima dell'errore
if e>tol %se la tolleranza non è rispettata riapplico la function sui due sottointervalli
    [lI,nfevall]=adattivasimp(f,a, m, tol/2, fa, fma, fm );
    [rI,nfevalr]=adattivasimp(f,m, b, tol/2, fm, fmb, fb ); 
    I2=lI+rI;
    nfeval=nfeval+nfevalr+nfevall; %sommo le valutazioni funzionali richieste nei due sottointervalli
end
return
end