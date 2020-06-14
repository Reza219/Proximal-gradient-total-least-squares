clear
load('A0h-gau.mat')

N=40;                                 %length of x
M=20;                                 %rows of A
K=5;                                  %support size
sA=.01/M;                             %A noise variance
sb=sA;                                %b noise variance
lam=.02;                              %regularization parameter
ni=150;                               %no. of iterations
nr=100;                               %no. of Monte Carlo runs

m21=zeros(ni,1);                      %mean-square error
m22=zeros(ni,1);                      %mean-square error

for ii=1:nr
    ii
    A=A0  +sqrt(sA)*randn(M,N);       %noisy A matrix
    b=A0*h+sqrt(sb)*randn(M,1);       %noisy b vector
    
    [e21,~,~,~]=adm_cd_stls_f(A,b,M,N,K,lam,h,ni);
    [e23,~,~,~]=ass_pg_stls_f(A,b,N,K,lam,h,ni);
    
    m21=m21+e21;
    m22=m22+e23;
end

m21=m21/nr;
m22=m22/nr;

figure
plot(10*log10(m21),'b','linewidth',3)
hold on
plot(10*log10(m22),'g','linewidth',3)
legend('AD-CD','proposed')