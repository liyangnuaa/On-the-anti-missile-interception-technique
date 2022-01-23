function y=fx(x)

V0=500;
s0=5000;
S=0.0606987;
m=250;
L=2.8;
D=0.278;
Jz=78.45;
R=6371000;
Vs=20.046*sqrt(288.34-5.86*1e-3*x(2)*s0);     %% sound speed
g=9.806*(1-2*x(2)*s0/R);
rou=1.225*(1-2.0323*1e-5*x(2)*s0)^4.83;
alpha=x(6)-x(4);
ommegaz_=x(5)*D/(V0*x(3));

if x(3)*V0/Vs<0.6
    Cx0=0.1423;
    Cxa2=0.00178;
    Cya=0.06896;
    mzoz=-13.1058;
    mza=-0.1143;
elseif x(3)*V0/Vs<0.8
    k1=(0.8-x(3)*V0/Vs)/(0.8-0.6);
    k2=(x(3)*V0/Vs-0.6)/(0.8-0.6);
    Cx0=0.1423*k1+0.1556*k2;
    Cxa2=0.00178*k1+0.00189*k2;
    Cya=0.06896*k1+0.06932*k2;
    mzoz=-13.1058*k1-12.9125*k2;
    mza=-0.1143*k1-0.1129*k2;
elseif x(3)*V0/Vs<0.9
    k1=(0.9-x(3)*V0/Vs)/(0.9-0.8);
    k2=(x(3)*V0/Vs-0.8)/(0.9-0.8);
    Cx0=0.1556*k1+0.1687*k2;
    Cxa2=0.00189*k1+0.00191*k2;
    Cya=0.06932*k1+0.06959*k2;
    mzoz=-12.9125*k1-13.9052*k2;
    mza=-0.1129*k1-0.1095*k2;
elseif x(3)*V0/Vs<1.0
    k1=(1-x(3)*V0/Vs)/(1-0.9);
    k2=(x(3)*V0/Vs-0.9)/(1-0.9);
    Cx0=0.1687*k1+0.2771*k2;
    Cxa2=0.00191*k1+0.00198*k2;
    Cya=0.06959*k1+0.07099*k2;
    mzoz=-13.9052*k1-14.5589*k2;
    mza=-0.1095*k1-0.1219*k2;
elseif x(3)*V0/Vs<1.1
    k1=(1.1-x(3)*V0/Vs)/(1.1-1);
    k2=(x(3)*V0/Vs-1)/(1.1-1);
    Cx0=0.2771*k1+0.3365*k2;
    Cxa2=0.00198*k1+0.00207*k2;
    Cya=0.07099*k1+0.07395*k2;
    mzoz=-14.5589*k1-15.3402*k2;
    mza=-0.1219*k1-0.1297*k2;
else
    Cx0=0.3365;
    Cxa2=0.00207;
    Cya=0.07395;
    mzoz=-15.3402;
    mza=-0.1297;
end

cx=Cx0+Cxa2*alpha^2;
cy=Cya*alpha;
X=cx*0.5*rou*V0^2*x(3)^2*S;
Y=cy*0.5*rou*V0^2*x(3)^2*S;
Mz=(mza*alpha+mzoz*ommegaz_)*0.5*rou*V0^2*x(3)^2*S*L;

y=zeros(size(x));
y(1)=V0/s0*x(3)*cos(x(4));
y(2)=V0/s0*x(3)*sin(x(4));
y(3)=-(X+m*g*sin(x(4)))/(m*V0);
y(4)=(Y-m*g*cos(x(4)))/(m*V0*x(3));
y(5)=Mz/Jz;
y(6)=x(5);

