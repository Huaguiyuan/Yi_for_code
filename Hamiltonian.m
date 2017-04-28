function HH=Hamiltonian(ky,kz)
global mu delta_x delta_y tm ic NX R_shift pbc
Nx=NX;
s_x=[0,1;1,0];s_y=[0,-ic;ic,0];s_z=[1,0;0,-1];eye2=eye(2);
for n=1:Nx; Id(n,:)=(n-1)*2+(1:2);end

HK=0;HmK=0;HH=0;

% HK=sparse(zeros(2,2));      
% HmK=sparse(zeros(2,2));      
% HH=sparse(zeros(2,2));       

 HK(2*NX,2*NX)=1;  HK(2*NX,2*NX)=0;
 HmK(2*NX,2*NX)=1; HmK(2*NX,2*NX)=0;
 HH(4*NX,4*NX)=1; HH(4*NX,4*NX)=0;

H0=zeros(2,2);H1=zeros(2,2);
tp= -(ky^2+kz^2)-R_shift;vky=2*ky;
H0=tp*s_x + vky*s_z - mu*eye2;

H1=(tm*s_x - s_y*ic)  *1d0/2;
for n=1:Nx;   HK(Id(n,:),Id(n,:))=H0;     end;
for n=1:Nx-1; HK(Id(n,:),Id(n+1,:))=H1;    end; if pbc&&(NX>2);HK(Id(NX,:),Id(1,:))=H1;end;  
for n=2:Nx;   HK(Id(n,:),Id(n-1,:))=H1';   end; if pbc&&(NX>2);HK(Id(1,:),Id(NX,:))=H1';end;

H0=zeros(2,2);H1=zeros(2,2);

tp= -(ky^2+kz^2)-R_shift;vky=2*(-ky);
H0=tp*s_x + vky*s_z - mu*eye2;
H1=(tm*s_x - s_y*ic)  *1d0/2;

for n=1:Nx;   HmK(Id(n,:),Id(n,:))=H0;     end;
for n=1:Nx-1; HmK(Id(n,:),Id(n+1,:))=H1;    end; if pbc&&(NX>2);HmK(Id(NX,:),Id(1,:))=H1;end; 
for n=2:Nx;   HmK(Id(n,:),Id(n-1,:))=H1';   end; if pbc&&(NX>2);HmK(Id(1,:),Id(NX,:))=H1';end;


D=zeros(2*Nx,2*Nx);%paring between K and -K states
D0=zeros(2,2);D1=zeros(2,2);
D0=2*ic*delta_y*sin(ky) *s_x; 
%D0=2*delta_y*sin(ky) *s_y;
D1= 0.5* delta_x *eye2;  %if change 2 to 0.5, Yi's result is exactly reproduced.
for n=1:Nx;   D(Id(n,:),Id(n,:))=D0;     end;
for n=1:Nx-1; D(Id(n,:),Id(n+1,:))=D1;    end;  if pbc&&(NX>2); D(Id(NX,:),Id(1,:))=D1; end;
for n=2:Nx;   D(Id(n,:),Id(n-1,:))=-D1;   end;  if pbc&&(NX>2); D(Id(1,:),Id(NX,:))=-D1;end;


%H=[H0,D0;-conj(D0),-conj(H0)];  V=[H1,D1;-conj(D1),-conj(H1)];
 HId=zeros(2,2*Nx); for n=1:2; HId(n,:)=(n-1)*2*Nx+(1:2*Nx);end
HH(HId(1,:),HId(1,:))=HK;         HH(HId(2,:),HId(2,:))=-conj(HmK);
HH(HId(1,:),HId(2,:))=D;  HH(HId(2,:),HId(1,:))=D';
% 
% ky
% kz
% mu
% tp
% tm
% vky
% delta_x
% delta_y
% sss=sum(sum(abs(HH),1))
% stop
end

