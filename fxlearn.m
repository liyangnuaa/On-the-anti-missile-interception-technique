function y=fxlearn(x)

global C1 C2 EXP0
N=length(EXP0(:,1));
V0=500;
s0=5000;

y=zeros(size(x));
y(1)=V0/s0*x(3)*cos(x(4));
y(2)=V0/s0*x(3)*sin(x(4));
for i=1:N
    y(3)=y(3)+C1(i)*x(2)^EXP0(i,1)*x(3)^EXP0(i,2)*x(4)^EXP0(i,3);
    y(4)=y(4)+C2(i)*x(2)^EXP0(i,1)*x(3)^EXP0(i,2)*x(4)^EXP0(i,3);
end

