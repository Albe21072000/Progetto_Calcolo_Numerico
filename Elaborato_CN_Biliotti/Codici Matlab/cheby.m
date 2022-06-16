function x =cheby(n,a,b)
%x =cheby(n,a,b)
%Input
%n: grado del polinomio interpolante, a,b: estremi dell'intervallo
%di interpolazione
%Output
%x: ascisse di Chebyshev ricavate
%
%Ricava le ascisse di Chebyshev nell'intervallo [a,b]
%per un polinomio interpolante di grado n 
if(n<0), error('Grado del polinomio interpolante non valido!'),end
if(a>=b), error('Intervallo definito in maniera non corretta!'),end
n=n+1;
x(n:-1:1)=(a+b)/2+((b-a)/2)*cos(((2*(1:n)-1)*pi)/(2*n));