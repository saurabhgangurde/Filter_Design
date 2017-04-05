A=-20*log10(0.1);

BL=24.2;
BH=27.2;
Fs=90;

bl=BL*2*3.14/Fs;
bh=BH*2*3.14/Fs;

omegap=(bh-bl)/2;
omegas=omegap+2*3.14/Fs;

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