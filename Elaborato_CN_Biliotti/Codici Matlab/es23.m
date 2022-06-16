f=@(x) sin(pi.*x.^2); %funzione non perturbata
fper=@(x) f(x)+0.1*rand(size(x)); %funzione perturbata
n=10^4;
x=(1:n)./n;
yper=zeros(size(x));
for i=1:n
    yper(i)=fper(x(i));
end
van=fliplr(vander(x)); %creo la matrice di vandermonde 
coef=zeros(15,15);
errmax=zeros(15,1);
for m=1:15
    coef(1:m+1,m)=miaqr(van(:,1:m+1),yper'); %ottengo il vetttore dei coefficienti dei vari polinomi di approssimazione ai minimi quadrati
end
p=zeros(size(x,2),15);
error=zeros(size(p));
y=f(x);
for m=1:15
p(:,m)=horner(coef(:,m),x')'; %sfrutto l'algoritmo di Horner per calcolare i valori dei vari polinomi nei punti che mi interessano
error(1:n,m)=abs(p(1:n,m)-y');
errmax(m)=max(error(1:n,m)); %mi ricavo il massimo valore dell'errore per ciascuno dei polinomi di vario grado
end
semilogy(1:15,errmax); %disegno il grafico semilogy dell'errore massimo in base al grado m del polinomio