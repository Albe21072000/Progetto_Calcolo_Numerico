function p=horner(a,x)
%p=horner(a,x)
%
%calcola il polinomio con coefficienti a nei valori x
n=length(a);
p=ones(length(x),1).*a(n);
for k=n-1:-1:1
    p=p.*x+a(k);
end