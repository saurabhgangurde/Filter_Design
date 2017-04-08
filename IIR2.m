
M=82;
q=floor(0.1*M);
r=M-10*q;

BL=2+0.6*q+1.5*r;
BH=BL+3;

omega_p2=BH+1;
omega_p1=BL-1;

B=omega_p2-omega_p1;
omega0_2=omega_p1*omega_p2;

syms s;
F=(B*s)/(omega0_2-s^2);

omega_p=abs(double(subs(F,omega_p2)));

omega_s1=double(subs(F,BL));
omega_s2=double(subs(F,BH));

omega_s=min(abs(omega_s1),abs(omega_s2));

delta1=0.1;
delta2=0.1;

D1=1/(1-delta1)^2-1;
D2=1/delta2^2-1;

N=ceil(0.5*log(D2/D1)/log(omega_s/omega_p));

omega_c=1/(2*N)*log(1/D1)*omega_p;

poles=[];

for j=0:N-1
    poles=[poles, 1i*omega_c*exp(1i*pi/(2*N)*(2*j+1))];
end


transformed_poles=[];


for j=1:N
    
    new_roots=roots([1,-1*B/poles(j),omega0_2]);
   transformed_poles=[transformed_poles,new_roots(1),new_roots(2)];  
end

numerator=conv([1,0,omega0_2],[1,0,omega0_2]);

for j=1:N-1
    numerator=conv(numerator,[1,0,omega0_2]);
end

Hbandpass=zpk([],transformed_poles,1)*tf(numerator,1);

bode(Hbandpass);

sysd=c2d(Hbandpass,1/90,'tustin');
[num,denom]=tfdata(sysd,'v');
freqz(num,denom);






