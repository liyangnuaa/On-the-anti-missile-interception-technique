clear;
clc;

z1min=0;
z1max=0;
z2min=0;
z2max=1.2;
z3min=0.1;
z3max=1.5;
z4min=-pi/2;
z4max=pi/2;
z5min=-1;
z5max=1;
z6min=-pi/2;
z6max=pi/2;
Nz=1e6;
z0x=zeros(6,Nz);
z0x(1,:)=rand(1,Nz)*(z1max-z1min)+z1min;
z0x(2,:)=rand(1,Nz)*(z2max-z2min)+z2min;
z0x(3,:)=rand(1,Nz)*(z3max-z3min)+z3min;
z0x(4,:)=rand(1,Nz)*(z4max-z4min)+z4min;
z0x(5,:)=rand(1,Nz)*(z5max-z5min)+z5min;
z0x(6,:)=rand(1,Nz)*(z6max-z6min)+z6min;

h=0.001;

%%% Generate data
% stochastic data
sigma=0.01;
Bh=sqrt(h)*randn(6,Nz);
zxf=zeros(6,Nz);
for i=1:Nz
    zxf(:,i)=z0x(:,i)+h*fx(z0x(:,i))+sigma*Bh(:,i);
end

% % deterministic data
% zxf=zeros(6,Nz);
% for i=1:Nz
%     zxf(:,i)=z0x(:,i)+h*fx(z0x(:,i));
% end

%%% Identify drift and diffusion terms using sparse learning
%%% EDMD algrithm
Ncoef=3;
p=(Ncoef+1)*(Ncoef+2)*(Ncoef+3)/6;
EXP=[0 0 0];
psiX = zeros(Nz,p);  %% observable for X
psiY = zeros(Nz,p);  %% observable for Y
psiX(:,1)=1;
psiY(:,1)=1;
for i=1:Ncoef
    for j=i:(-1):0
        for k=i-j:(-1):0
            psiX(:,i*(i+1)*(i+2)/6+(i-j)*(i-j+1)/2+i-j-k+1)=z0x(2,:)'.^j.*z0x(3,:)'.^k.*z0x(4,:)'.^(i-j-k);
            psiY(:,i*(i+1)*(i+2)/6+(i-j)*(i-j+1)/2+i-j-k+1)=zxf(2,:)'.^j.*zxf(3,:)'.^k.*zxf(4,:)'.^(i-j-k);
            EXP=[EXP;j k i-j-k];
        end
    end
end
G=psiX'*psiX/Nz;
A=psiX'*psiY/Nz;
K = pinv(G) * A;  %%Koopman operator
B3 = zeros(1,p);
B4 = zeros(1,p);
B3(3)=1;
B4(4)=1;

delta=1e-3;
L3 = ((K*B3')-B3')/h; %% Generator for v
I=abs(L3)<delta;
L3(I)=0;
L4 = ((K*B4')-B4')/h; %% Generator for theta
I=abs(L4)<delta;
L4(I)=0;

%%% test
T=35;
t=0:h:T;
nT=length(t);
theta=0.1;
x0=[0;0.96;0.6960;theta;0.001;0.00062];
y0=x0(1:4);
x=zeros(6,nT);
y=zeros(4,nT);
x(:,1)=x0;
y(:,1)=y0;

for i=1:nT-1
    x(:,i+1)=rk4(h,x(:,i));
end

global C1 C2 EXP0
C1=L3;
C2=L4;
EXP0=EXP;

for i=1:nT-1
    y(:,i+1)=rk4learn(h,y(:,i));
end

figure;
plot(t,x(1,:));

figure;
plot(t,x(2,:));

figure;
plot(t,x(3,:));

figure;
plot(t,x(4,:));

% figure;
% plot(t,x(5,:));
% 
% figure;
% plot(t,x(6,:));

figure;
plot(t,y(1,:));

figure;
plot(t,y(2,:));

figure;
plot(t,y(3,:));

figure;
plot(t,y(4,:));
