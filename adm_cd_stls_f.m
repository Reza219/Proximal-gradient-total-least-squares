function [er2,er0a,er0b,x] = adm_cd_stls_f(A,b,M,N,K,lam,h,ni)
%alternating direction minimization, coordinate descent

x=zeros(N,1);                     %initialization of solution
E=zeros(M,N);                     %initialization of perturbation matrix
er2=zeros(ni,1);                  %error
er0a=zeros(ni,1);                 %missed detections
er0b=zeros(ni,1);                 %wrong detections
lam=lam/2;

for nn=1:ni                       %iterations loop
    A=A+E;
    
    for ii=1:N                    %coordinate-descent loop
        %calculate e(ii)
        e=b;
        for jj=setdiff(1:N,ii)
            e=e-A(:,jj)*x(jj);
        end
        
        %update x(ii)
        a=A(:,ii);
        if     e'*a>lam
            x(ii)=(e'*a-lam)/(a'*a);
        elseif e'*a<-lam
            x(ii)=(e'*a+lam)/(a'*a);
        else
            x(ii)=0;
        end
    end
    A=A-E;
    E=(b-A*x)*x'/(x'*x+1);
    er2(nn)=norm(x-h)^2;
    ll=length(intersect(find(h),find(x)));
    er0a(nn)=K-ll;
    er0b(nn)=length(find(x))-ll;
end