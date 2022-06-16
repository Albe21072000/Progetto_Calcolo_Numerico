function [x,nr] = miaqr(A,b)
%[x,nr] = miaqr(A,b)
%Input: A matrice sovradimensionata a rango massimo(non quadrata, vedasi function mialu computazionalmente più efficiente), 
%b vettore colonna dei termini noti.
%Output: x soluzione del sistema lineare nel senso dei minimi quadrati, 
% nr: norma 2 del vettore residuo.
%Calcola la soluzione, nel senso dei minimi quadrati, del sistema lineare
%sovradeterminato passato in input e restiuisce anche la norma due del
%vettore dei residui (minimizzata grazie al metodo dei minimi quadrati).
[m,n] = size(A);
k=size(b,1);
if (n>=m),error("Sistema non sovradeterminato, conviene usare un'altro algoritmo più efficiente (vedasi function mialu)!"),end
if(m~=k), error("Le lunghezze della matrice e del vettore dei termini noti non coincidono!"), end
if(size(b,2)>1),error("Inserire un vettore colonna come vettore dei termini noti!"),end
QR=QRdecompose(A); %decompongo 
[x,nr]=QRsolver(QR,b);
return;
end


function QR=QRdecompose(A)
%QR=QRdecompose(A)
%Input: A matrice sovradimensionata a rango massimo
%Output: QR: matrice fattorizzata QR
%Decompone A matrice sovradimensionata nella sua fattorizazione QR dove Q è una matrice ortogonale
% e R una matrice composta nelle prime n righe da una matrice triangolare superiore di dimensione nXn ed il resto delle righe nullo
[m,n] = size(A);
norma2=norm(A);
QR=A;
for i = 1 : n
alpha = norm(QR(i : m, i));
if abs(alpha)<=eps*norma2
error('Matrice non a rango massimo!');
end
if QR(i,i)>=0
alpha=-alpha;
end
v1=QR(i,i)-alpha;
QR(i,i)=alpha;
QR(i+1:m,i)=QR(i+1:m,i)/v1;
beta =-v1/alpha;
QR(i:m,i+1:n)=QR(i: m,i+1:n)-(beta*[1;QR(i+1:m,i)])*([1 QR(i+1:m,i)']*QR(i:m,i+1:n));
end
end



function [x,nr] = QRsolver(QR,b)
%[x,nr] = QRsolver(QR,b)
%Input: QR: matrice precedentemente fattorizzata QR, b: vettore dei termini
%noti
%Output: x: soluzione del sistema lineare in input
%Risolve il sistema lineare con la scomposizione QR della matrice
%precedentemnte ottenuta.
x=b;
[m,n] = size(QR);
for i=1:n
v=[1;QR(i+1:m,i)];
x(i:m)=x(i:m)-((2*(v*v'))/(v'*v))*x(i:m); %Non necessito di "costruirmi" la matrice Q ma utilizzo solo i vettori v
end
x(1:n)=trsup(QR(1:n , 1 : n), x(1 : n));
nr=norm(x(n+1:m)); %mi ricavo la norma euclidea del vettore residuo
x=x(1:n);
end

function [b] = trsup(A,b)
%[b] = trsup(A,b)
%Input: A matrice quadrata triangolare superiore, b vettore dei termini
%noti
%Output: x: soluzione del sistema lineare in input
%Risolve il sistema lineare Ax=b passato in input dove A è una matrice triangolare superiore.
n=size(A,1);
for i=n:-1:1
    if(n>1)
    b(i)=b(i)-A(i,i+1:n)*b(i+1:n);
    end
    b(i)=b(i)/A(i,i);
end
end