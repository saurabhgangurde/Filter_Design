A=-20*log10(0.1);

M=82;
q=floor(0.1*M);
r=M-10*q;

BL=2+0.6*q+1.5*r;
BH=BL+3;
Fs=90;

bl=BL*2*pi/Fs;
bh=BH*2*pi/Fs;

omegap=(bh-bl)/2;
omegas=omegap+2*pi/Fs;

omegac=(omegap+omegas)/2;
omegat=omegas-omegap;

N=(2+(A-8)/(2.285*omegat))/2;
N=ceil(N);

kaisercoff=kaiser(N,0);

n=0:N-1;
h1=sin(omegac*n)./(3.14*n);
h1(1)=1;

h=h1.*kaisercoff'.*cos((bl+bh)/2*n);
[trans_fun,w] = freqz(h,1);
plot(w/pi*45,abs(trans_fun));
freqz(h,1)