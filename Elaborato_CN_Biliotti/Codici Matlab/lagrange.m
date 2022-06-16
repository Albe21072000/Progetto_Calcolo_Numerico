function YQ=lagrange(X,Y,XQ)
% YQ=lagrange(X,Y,XQ)
%
%Input
%X: Vettore colonna contenete le ascisse d'interpolazione che devono essere distinte l'una dall'altra, 
%Y: Vettore colonna contenete i valori della funzione nelle ascisse d'interpolazione
%XQ: Vettore colonna contenente le ascisse in cui vogliamo approssimare la funzione
%Output
%YQ: Valori approssimati della funzione con il polinomio interpolante in
%forma di Lagrange
%
%Calcola i valori approssimati della funzione(di cui conosciamo i valori Y che assume nelle ascisse X) calcolati attraverso il polinomio interpolante in
%forma di Lagrange nelle ascisse XQ.
% 
if(length(X)~=length(Y)), error("Numero di ascisse d'interpolazione e di valori della funzione non uguale!"),end
if (length(X) ~= length(unique(X))), error("Ascisse d'interpolazione non distinte!"),end %uso la function unique che restuisce un vettore contenente i valori senza ripetizioni di X
if(isempty(XQ)), error("Il vettore contenente le ascisse in cui interpolare la funzione è vuoto!"),end
if(size(X,2)>1||size(Y,2)>1||size(XQ,2)>1),error("Inserire vettori colonna!"),end
n=size(X,1);
L=ones(size(XQ,1),n);
   for i=1:n
      for j=1:n
         if (i~=j)
            L(:,i)=L(:,i).*((XQ-X(j))/(X(i)-X(j)));   %calcolo i polinomi di base di lagrange Lin(x)
         end
      end
   end
   YQ=zeros(size(XQ));
   for i=1:n
      YQ=YQ+Y(i).*L(:,i);    %calcolo la sommatoria dei prodotti fi*Lin(x) (con i=0,...,n))
   end
end