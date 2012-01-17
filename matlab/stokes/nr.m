function nr(a) 

% plot the numerical range
% Author: Panos Psarrakos

[n,m]=size(a);
e1=1;r=2;p=[];q=[];i=sqrt(-1);normdif=norm(full(a),2);
while  e1>0.01
    r=2*r;
    d=(2*pi)/r;
    for  k=1:r
      B=(cos(k*d)+i*sin(k*d))*a;
      H=0.5*(B+B');
      [X,L]=eig(H);
      g=diag(real(L));
      for  c=1:n
           if  g(c)==max(g)
                  l(k)=g(c);
                  IV(:,k)=X(:,c);
           end
      end
    end
    l(k+1)=l(1);
    for k=1:r
      p(k)=IV(:,k)'*a*IV(:,k);
      f=(l(k)*cos(d)-l(k+1))/sin(d);
      q(k)=(cos(k*d)-i*sin(k*d))*(l(k)+(i*f));
    end
    p(r+1)=p(1);
    q(r+1)=q(1);
    s1=0;
    s2=0;
    for  k=1:r
       s1=s1+q(k)'*q(k+1);
       s2=s2+p(k)'*p(k+1);
    end
    e1=0.5*imag(s2-s1)/normdif;
    plot(p, 'b');
    fill(real(p),imag(p),'c');
end
   