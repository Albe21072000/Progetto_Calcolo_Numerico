function YQ = spline0(X, Y, XQ)
% YQ=spline=(X,Y,XQ)
%
%Input
%X: Vettore colonna contenente le ascisse d'interpolazione distinte, 
%Y: Vettore colonna contenente i valori che la funzione assume nelle ascisse d'interpolazione
%XQ: Vettore contente le ascisse in cui vogliamo approssimare la funzione
%Output
%YQ: Valori approssimati della funzione interpolata attraverso una funzione
%spline cubica naturale
%
%Calcola i valori approssimati della funzione(di cui conosciamo i valori Y
%che assume nelle ascisse X) calcolati attraverso una funzione spline
%cubica naturale interpolante la funzione.
% 
if(length(X)~=length(Y)), error("Lunghezza dei dati in ingresso non valida!"),end
if (length(X) ~= length(unique(X))), error("Ascisse d'interpolazione non distinte!"),end %uso la function unique che restuisce un vettore contenente i valori senza ripetizioni di X
if(isempty(XQ)), error("Il vettore contenente le ascisse in cui interpolare la funzione è vuoto!"),end
if(size(X,2)>1||size(Y,2)>1),error("Inserire vettori colonna!"),end
n=length(X)-1;
eps=zeros(n-1,1);
p=zeros(n-1,1);
for i = 1 : n - 1
    p(i)=(X(i+1)-X(i))/(X(i+2)-X(i));
    eps(i)=(X(i+2)-X(i+1))/(X(i+2)-X(i));
end
dif = diffdivspline(X, Y);
m = tridia(2*ones(n-1,1), p, eps,dif*6);
m=[0;m;0];
YQ = calcolapuntixq(X, Y, m, XQ);
end

function[ df ] = diffdivspline(x, f)
%calcola le differenze divise fermandosi pero' alla seconda iterazione
% per ottenere le differenze divise tra tre punti che ci serviranno per calcolare m 
n=size(x);
df=f;
n=n-1;
for j=1:2
    for i=n+1:-1:j+1
        df(i)=(df(i)-df(i-1))/(x(i)-x(i-j));
    end
end
df=df(3:n+1);
return;
end

function x=tridia(a,b,c,x)
%calcola il sistema tridiagonale in cui, nel nostro caso, il vettore a
%contiene solo 2 mentre i vettori b e c rispettivamente i valori p ed eps
%mentre x è il vettore delle differenze divise per ricavare il vettore
%contente i valori m
n=length(a);
for i=1:n-1                     % mi ricavo una fattorizzazione LU della matrice tridiagonale risolvendo contemporaneamente il fattore L 
   b(i)=b(i)/a(i);
   a(i+1)=a(i+1)-b(i)*c(i);
   x(i+1)=x(i+1)-b(i)*x(i);
end
x(n)=x(n)/a(n);
for i=n-1:-1:1
   x(i)=(x(i)-c(i)*x(i+1))/a(i);  %risolvo il fattore U rimanente della scomposizione LU
end
return
end

function s = calcolapuntixq( xi, fi, m, XQ)
%Calcola i valori della spline precedentemente ricavata grazie ai valori m
%nei punti XQ anche se essi si trovano al di fuori dell'intervallo d'interpolazione 
n = length(xi)-1;
z=length(XQ);
s=zeros(z,1);
[xi,p]=sort(xi);  %ordino i vettori delle ascisse in modo che formino una partizione vera e propria e di conseguenza ordino anche il vettore dei corrispondenti valori della funzione in tali ascisse
fi=fi(p);
for j=1:z
    for i=2:n+1 
        if((XQ(j)>=xi(i-1)&& XQ(j)<=xi(i))||XQ(j)<xi(1))
            hi=xi(i)-xi(i-1);
            ri=fi(i-1)-hi^2/6*m(i-1);
            qi=(fi(i)-fi(i-1))/hi-hi/6*(m(i)-m(i-1));
            s(j) =((XQ(j)-xi(i-1))^3*m(i)+(xi(i)-XQ(j))^3*m(i-1))/(6*hi)+qi*(XQ(j)-xi(i-1))+ri;
            break
        elseif(XQ(j)>xi(n+1))      %questa parte serve solo per fare in modo che, se anche cerco di approssimare un valore della funzione in un ascissa più grande dell'estremo superiore dell'intervallo d'interpolazione, 
            % la function restituisca comunque un valore che è calcolato nella funzione dell'ultimo tratto da cui è composta la spline (come avviene per la function spline) 
            hi=xi(n+1)-xi(n);
            ri=fi(n)-hi^2/6*m(n);
            qi=(fi(n+1)-fi(n))/hi-hi/6*(m(n+1)-m(n));
            s(j) =((XQ(j)-xi(n))^3*m(n+1)+(xi(n+1)-XQ(j))^3*m(n))/(6*hi)+qi*(XQ(j)-xi(n))+ri;
            break
        end
    end
end
end