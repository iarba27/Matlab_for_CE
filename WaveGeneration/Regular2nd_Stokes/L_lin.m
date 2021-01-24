function [L1]=L_lin(T,d)

i_max=1000;
i=0;
L0=9.81*T^2/(2*pi);
L1=L0;
L2=0;
while abs(L0*tanh(2*pi*d/L1)-L1)>0.000001 && i<i_max
    L2=L0*tanh(2*pi*d/L1);
    L1=L2;
    i=i+1;
end    
end