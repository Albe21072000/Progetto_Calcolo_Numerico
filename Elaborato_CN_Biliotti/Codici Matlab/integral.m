function int =integral(a,b,f,n)
%
% int =integral(a,b,f,n)
%
% Input
%
%f: Function handler contenente la funzione di cui vogliamo calcolare l'approssimazione dell'integrale definito, 
%a,b: Estremi dell'intervallo d'integrazione, 
%n: Grado della formula di Newton-Cotes desiderato, 
%Output
%int: approssimazione dell'integrale ottenuta
%
%Calcola un'approssimazione dell'integrale definito di f tra a e b con la
%formula di Newton-Cotes di grado n
%Newton-Cotes di grado n.
if(n<=0), error('Grado della formula di Newton-Cotes non valido!'),end
x=linspace(a,b,n+1)';
y=f(x);
int=(b-a)/n*sum(y.*calcolacoefficientigrado(n));
return;
end