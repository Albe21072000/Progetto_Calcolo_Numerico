function YQ = hermitte(X,Y, Y1, XQ)
% 
% YQ=hermitte(X,Y,Y1,XQ)
%
%Input
%X: Vettore colonna contenete le ascisse d'interpolazione distinte, 
%Y: Vettore colonna contenete i valori della funzione nelle ascisse d'interpolazion
%Y1: Vettore colonna contenente i valori che la derivata prima della funzione assume nelle ascisse
%d'interpolazione
%XQ: Vettore contenente le ascisse in cui vogliamo approssimare la funzione
%
%Output
%YQ: Valori approssimati della funzione con il polinomio interpolante di Hermitte in
%forma di Newton
%
%Calcola i valori approssimati della funzione(di cui conosciamo sia i valori Y che assume nelle ascisse X ed i valori Y1 la cui derivata prima assume nelle stesse ascisse) 
%calcolati attraverso il polinomio interpolante in
%forma di Lagrange nelle ascisse XQ.
% 
if(isempty(XQ)), error("Il vettore contenente le ascisse in cui interpolare la funzione è vuoto!"),end
if(length(X)~=length(Y)), error("Numero di ascisse d'interpolazione e di valori della funzione non uguale!"),end
if (length(X) ~= length(unique(X))), error("Ascisse d'interpolazione non distinte!"),end %uso la function unique che restuisce un vettore contenente i valori senza ripetizioni di X
if(length(Y1)~=length(Y)), error("Lunghezza dei dati in ingresso non compatibile!"),end
if(size(X,2)>1||size(Y,2)>1||size(Y1,2)>1),error("Inserire vettori colonna!"),end
n=(length(X));
fi(1:2:2*n-1)=Y;
fi(2:2:2*n)=Y1;
df=diffdivhermitte(X,fi');
n=length(df)-1;
YQ=df(n+1)*ones(size(XQ));    %algoritmo di horner per il calcolo dei valori di un polinomio
for i=n:-1:1
    YQ=YQ.*(XQ-X(ceil(i/2)))+df(i);
end
return
end

function f=diffdivhermitte(x,f)
%function per calcolare le differenze divise per il polinomio interpolante
%di hermitte
n=(length(f)/2)-1;
for i=2*n+1:-2:3
    f(i)=(f(i)-f(i-2))/(x(ceil(i/2))-x(ceil((i-1)/2)));
end
for j= 2:2*n+1
    for i=(2*n+2):-1:j+1
        f(i)=(f(i)-f(i-1))/(x(ceil(i/2))-x(ceil((i-j)/2)));
    end
end
end