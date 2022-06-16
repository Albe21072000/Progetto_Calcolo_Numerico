function x = mialu(A,b)
%x=mialu(A,b)
%Input
%A: Matrice quadrata nonsingolare(det(A)!=0) del sistema da risolvere
%b: Vettore colonna dei termini noti
%Output
%x: Soluzione del sistema lineare
%Risolve il sistema lineare Ax=b, passato in input, con la fattorizazione LU con pivoting parziale
if(size(A,1)~=size(A,2)), error("Matrice non quadrata!"), end
if(size(b,2)>1),error("Inserire un vettore colonna come vettore dei termini noti!"),end
if(size(A,1)~=size(b,1)), error("Dimensione del vettore dei termini noti non compatibile con la matrice!"), end

[LU,p]=LUpivoting(A); %Fattorizzo la matrice A con la fattorizzazione LU con pivoting parziale
x=b;
x=x(p);      %scambio le posizioni degli elementi del vettore dei termini noti in accordo con quelli effettuati per il pivoting
x=LUsolver(LU, x); %Risolvo il sistema lineare 
return
end

function [LU,p]= LUpivoting(A)
% [LU,p]=LUpivoting(A)
% Restiusce la fattorizzazione LU della matrice nonsingolare A con pivoting
% parziale di cui teniamo traccia nel vettore p
LU=A;
n=size(LU);
p=1:n; %genero un vettore contenente i numeri da 1 ad n
for i=1:n-1
    [ma,k]=max(abs(LU(i:n,i)));   %individuo il valore massimo del vettore che prendo in considerazione
    if(ma==0), error("Matrice singolare!"), end
    k=k+i-1;
    if k>i
        LU([i k],:)=LU([k i], :); %scambio le righe 
        p([i k])=p([k i]);  %tengo traccia dello scambio nel vettore delle permutazioni
    end
    LU(i+1:n,i)=LU(i+1:n,i)/LU(i,i);
    LU(i+1:n,i+1:n)= LU(i+1:n,i+1:n)- LU(i+1:n,i)*LU(i, i+1:n);
end
return;
end

function x=LUsolver(LU,b)
%x=LUsolver(LU,b) 
% Risolve il sistema lineare LUx=b gia decomposto con la
%fattorizazione LU
x=b;
n=size(LU);
for i=1:n-1
    x(i+1:n)=x(i+1:n)-LU(i+1:n,i)*x(i);  %risoluzione fattore L per colonne
end
for i=n:-1:1
    x(i)=x(i)/LU(i,i);
    x(1:i-1)=x(1:i-1)-LU(1:i-1,i)*x(i); %risoluzione fattore U per colonne
end
return;
end