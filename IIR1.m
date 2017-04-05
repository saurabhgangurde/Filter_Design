
M=82;
q=floor(0.1*M);
r=M-10*q;

BL=2+0.7*q+2*r;
BH=BL+5;
B=BH-BL;
omega0_2=BL*BH;

syms s;
F=(s^2-omega0_2)/(B*s);

omega_p=double(subs(F,BH));

omega_s1=double(subs(F,BL-1));
omega_s2=double(subs(F,BH+1));

omega_s=min(abs(omega_s1),abs(omega_s2));

delta=0.1;

epsilon=sqrt(1/(1-delta)^2-1);
D2=1/delta^2-1;
N=ceil(acosh(sqrt(D2)/epsilon)/acosh(omega_s/omega_p));

coeff=zeros(1,N+1);
syms x;
cheby_coeffs=double(coeffs(chebyshevT(5,x)));
for i=1:1:N/2+1
    coeff(2*i-1)=cheby_coeffs(i);
end

denominator=epsilon^2*conv(coeff,coeff);
denominator(end)=denominator(end)+1;

H_s_2=tf([1],denominator);
