function HH=Add_disorder_to_Hamiltonian(HH0)
global mu delta_x delta_y tm ic NX U_Disorder
% in this model, the hopping in x drection is independent of ky,kz
Nx=NX;

s_x=[0,1;1,0];s_y=[0,-ic;ic,0];s_z=[1,0;0,-1];
eye2=eye(2);
for n=1:Nx; Id(n,:)=(n-1)*2+(1:2);end
for n=1:Nx;   H0_disorder(Id(n,:),Id(n,:))= U_Disorder*diag(rand(1,2)-0.5);     end;

% D=zeros(2*Nx,2*Nx);%paring between K and -K states
% D0=zeros(2,2);D1=zeros(2,2);
% D0=2*ic*delta_y*sin(ky) *s_x; 
% D1= delta_x *eye2; %%D0=2*delta_y*sin(ky) *s_y;
% for n=1:Nx;   D(Id(n,:),Id(n,:))=D0;     end;
% for n=1:Nx-1; D(Id(n,:),Id(n+1,:))=D1;    end;
% for n=2:Nx;   D(Id(n,:),Id(n-1,:))=-D1;   end;

%H=[H0,D0;-conj(D0),-conj(H0)];  V=[H1,D1;-conj(D1),-conj(H1)];
%HH=zeros(Nx*4,Nx*4); 

HId=zeros(2,2*Nx); for n=1:2; HId(n,:)=(n-1)*2*Nx+(1:2*Nx);end
HH=HH0;
HH(HId(1,:),HId(1,:))=HH(HId(1,:),HId(1,:))+ H0_disorder;    
HH(HId(2,:),HId(2,:))=HH(HId(2,:),HId(2,:))-conj(H0_disorder);
%HH(HId(1,:),HId(2,:))=D; 
%HH(HId(2,:),HId(1,:))=D';

end

