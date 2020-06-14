function [er2,er0a,er0b,x] = ass_pg_stls_f(A,b,N,K,lam,h,ni)
%adaptive-step-size proximal-gradient

AA=A'*A;
Ab=A'*b;
er2=zeros(ni,1);                    %error
er0a=zeros(ni,1);                   %missed detections
er0b=zeros(ni,1);                   %wrong detections
xo=zeros(N,1);                      %initialization of solution
g=-2*Ab;%-Ab/beta;
mu0=.2;%1;%.5;  %%%%%%%%%
x=wthresh(-mu0*g,'s',mu0*lam);
y=1/(x'*x+1);
c=y*norm(A*x-b)^2;
muo=mu0;
%mui=zeros(ni,1);

for nn=1:ni                         %iterations loop
    
    %calculate gradient
    go=g;
    co=c;
    g=2*y*(AA*x-Ab-co*x);
    
    %calculate step-size
    if (x-xo)'*(g-go)==0            % && (g-go)'*(g-go)~=0
        mu=muo;
    else
        mus=((x-xo)'*(x-xo))/((x-xo)'*(g-go));
        mum=((x-xo)'*(g-go))/((g-go)'*(g-go));
        if mum/mus>.5
            mu=mum;
        else
            mu=mus-mum/2;
        end
        if mu<=0
            mu=muo;
        end
    end
    
    %backtracking line-search
    while 1
        %proximal-gradient
        z=wthresh(x-mu*g,'s',mu*lam);
        y=1/(z'*z+1);
        c=y*norm(A*z-b)^2;
        if c<=co+(z-x)'*g+norm(z-x)^2/(2*mu)
            break
        end
        mu=mu/2;
    end
    %mui(nn)=c;    %%%%%%%%%%%
    muo=mu;
    xo=x;
    x=z;
    
    %calculate errors
    er2(nn)=norm(x-h)^2;
    ll=length(intersect(find(h),find(x)));
    er0a(nn)=K-ll;
    er0b(nn)=length(find(x))-ll;
end