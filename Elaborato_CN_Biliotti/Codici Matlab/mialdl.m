function x=mialdl(A,b)
%x=mialdlt(A,b)
%INPUT
%A: Matrice simmetrica definita positiva(sdp), b: Vettore dei termini noti
%OUTPUT
%x: Soluzione del sistema lineare Ax=b
%Risolve il sistema lineare passato in input grazie alla fattorizazione LDLt.  
if(size(A,1)~=size(A,2)), error("Matrice non quadrata!"), end
if(size(A,1)~=size(b)), error("Dimensione del vettore dei termini noti non compatibile con la matrice!"), end
LDLt=LDLtdecompose(A); %decompongo la matrice tramite fattorizzazione LDLt
x=trinfdiaguni(LDLt,b); %Risolvo Lx=b
x=diagonal(LDLt,x);  %Risolvo Dx=b dove b è la soluzione x ricavata in precedenza
x=trsupdiaguni(LDLt',x); %Risolvo L'x=b  dove b è la soluzione x ricavata in precedenza
return;
end

function LDLt = LDLtdecompose(A)
%LDLt = LDLtdecompose(A)
%decompone la matrice A simmetrica definita positiva tramite
%fattorizzazione LDLt dove L è una matrice triangolare inferiore a
%diagonale unitaria e D una matrice diagonale
n= size(A);
LDLt=A;
if(LDLt(1,1)<=0), error('Matrice non sdp!'),end
LDLt(2:n,1)=LDLt(2:n,1)/LDLt(1,1);
for i=2:n
v=(LDLt(i,1:i-1)').*diag(LDLt(1:i-1,1:i-1));
LDLt(i,i)=LDLt(i,i)-LDLt(i,1:(i-1))*v;
if(LDLt(i,i)<=0),error('Matrice non SDP!'),end
LDLt(i+1:n,i)=(LDLt(i+1:n,i)-LDLt(i+1:n,1:i-1)*v)/LDLt(i,i);
end
end

function b = trinfdiaguni(A,b)
%risolve il sistema lineare Ax=b con A matrice triangolare inferiore a diagonale unitaria 
n=size(A,1);
for i=1:n
b(i)=b(i)-A(i,1:i-1)*b(1:i-1);
end
end

function x = diagonal(D, b)
%x = diagonal(D, b)
%risolve il sistema lineare Dx=b con D matrice diagonale
n=size(D);
x=b;
for i=1:n
    x(i)=x(i)./D(i,i);
end
end

function x = trsupdiaguni(A,b)
%x = trsupdiaguni(A,b)
%risolve il sistema lineare Ax=b con A matrice triangolare superiore a diagonale unitaria 
n=size(A);
x=b;
for i=n:-1:1
x(i)=x(i)-A(i,i+1:n)*x(i+1:n);
end
end