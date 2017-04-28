
global mu delta_x delta_y tm ic NX U_Disorder R_shift  eta pbc disorder_type E
pbc=1;
mu=0.2;  E=0.03;
% E and ky,kz switched global properties :April 10,2017 

R_shift=0;  kz=0.85;   ky=0.1;

% R_shift=0.8;kz=0.6;
%  R_shift=1.02;kz=0.2;

U_Disorder=0; 
disorder_type=1;
N_dirsorder=1;

ic=sqrt(-1);
tm=1;
eta=10^-4;

delta_x=0.4; delta_y= ic*delta_x ;


NX=80;
NY=10;
NZ=10;
NE=1000;

 

Kymin=0.35; Emin=0.3;
kz1=sqrt(1-mu+R_shift);
kz2=sqrt(1+mu-R_shift); 
kzmax=max(abs(kz1),abs(kz2));

kzrr=0.2;
ss=(kzmax+kzrr-(-kzmax-kzrr))/NZ;KZ_region=-kzmax-kzrr:ss:kzmax+kzrr-ss;
k1max=0.5;k2max=1.5; ss=(k2max-k1max)/NZ; KZ_region=k1max:ss:k2max-ss;

KY_region=-Kymin:2*Kymin/NY:Kymin-2*Kymin/NY;

E_region=-Emin:2*Emin/NE:Emin-2*Emin/NE;

 DOS=zeros(size(E_region));
LDOS_left=zeros(size(E_region));
LDOS_right=zeros(size(E_region));
