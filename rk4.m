function xout=rk4(h,x0)
k1=h*fx(x0);
k2=h*fx(x0+0.5*k1);
k3=h*fx(x0+0.5*k2);
k4=h*fx(x0+k3);
xout=x0+(k1+2*k2+2*k3+k4)/6;