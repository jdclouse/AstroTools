close all
dtr = pi/180;

n = 500
x = [randn(n,1) randn(n,1)];

mx = 1.5*max(max(x));
plot(x(:,1),x(:,2),'g.'),xlim([-mx mx]),ylim([-mx mx])
axis('square'),grid
hold on


me=mean(x)
sig=std(x)
xe=x-ones(n,1)*me;

raderr=sqrt(sum(sum(xe.*xe))/n)
drawellipse(raderr,raderr,0,me(1),me(2),'r-.');

P=cov(xe)
[evec,ev]=eig(P);
semimaj=sqrt(ev(1,1)), evec(:,1)'
semimin=sqrt(ev(2,2)), evec(:,2)'
th=atan2(evec(2,1),evec(1,1));
thd=th/dtr
C=drawellipse(semimaj,semimin,th,me(1),me(2),'k');
