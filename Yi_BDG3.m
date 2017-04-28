
clear; %smaller 2*Nx *2 BDG matrix:  Plot density of states for  surface states
       % used a criteria that can deside a eigenstates is a surface state
       % or not.
vars;

filename=strcat('TRY-','r= ',num2str(R_shift),',mu=',num2str(mu),',Disorder= ',...
    num2str(U_Disorder), ',NX=',num2str(NX),',NY=',num2str(NY),',NZ=',num2str(NZ),...
    ',DeltaX=',num2str(delta_x),',eta=',num2str(eta));
filename

tic;


for ns=1:N_dirsorder
     DOS=zeros(size(E_region));
% LDOS_left=zeros(size(E_region));
% LDOS_right=zeros(size(E_region));


%for ky=0.01:0.01:0.01   %KY_region;
nz=0;
for kz=  KZ_region;
    nz=nz+1;%ny=0;
nz
dos=zeros(NY,NE);

for ny=1:NY
    ky=KY_region(ny);
HH0=Hamiltonian(ky,kz);

HH=Add_disorder_to_Hamiltonian(HH0);

% cp,tp,tm,V,delt_x,delt_y
% sum(sum(abs(HH),1))
% stop

[v1,r]=eig(HH);r=diag(r); [r,Id]=sort(r); v1=v1(:,Id);

ind = (find(r>=-Emin&r<Emin))'; NumofStates(ny,nz)=size(ind,2);
for s=ind
%for s=1:4*NX
    nE=1;
   
%     a=abs(v1(:,s)').^2; b=a(1:2*NX)+a(2*NX+1:4*NX);
%     a=[sum(b(1:NX)),sum(b(NX+1:2*NX))];
%     if (abs(a(1)/a(2))>0.5)&(abs(a(1)/a(2))<2)
%       %  continue
%      end;
%     
%     surfden1=abs(v1(1,s))^2+abs(v1(2,s))^2+abs(v1(2*NX+1,s))^2+abs(v1(2*NX+2,s))^2;
%     surfden2=abs(v1(2*NX-1,s))^2+abs(v1(2*NX,s))^2+abs(v1(4*NX-1,s))^2+abs(v1(4*NX,s))^2;

%E=E_region;
%dos(ny,:)=dos(ny,:)+eta*((E_region-r(s)).^2+eta^2).^-1/pi;
DOS=DOS+eta*((E_region-r(s)).^2+eta^2).^-1/pi;
% for E=E_region;
   %  DOS(nE)=DOS(nE)+eta/((E-r(s))^2+eta^2)/pi;
% %    LDOS_left(nE)=LDOS_left(nE)+surfden1*eta/((E-r(s))^2+eta^2)/pi;
% %     LDOS_right(nE)=LDOS_right(nE)+surfden2*eta/((E-r(s))^2+eta^2)/pi;
%     nE=nE+1;
% end;
end;

% energy(ny,:)=real(r)'; Ky(ny)=ky;

                
% if abs(ky-ky00)<10^-4;
%     [v,rr]=eig(HH);
% end;
% 123;
% ny=ny+1;
end;
end;

% for i=NX*2-3:2*NX+3
% a=abs(v1(:,i))';b=a(1:2*NX)+a(2*NX+1:4*NX);
% figure; plot(b);title(num2str(i))
% end
 %plot(E_region,DOS,'-b',E_region,LDOS_left,'-r',E_region,LDOS_right,'-y');;
 

figure;  plot(E_region,DOS,'-b');; ylabel('DOS');xlabel('E');
title(filename);
savefig(strcat(filename,'DOS','.fig'));

figure;mesh(KZ_region,KY_region,NumofStates);xlabel('kz');ylabel('ky');
zlabel('Num of states within energy window');
movegui(gcf,'southwest');
title(filename);
savefig(strcat(filename,'number-of-States-Ky-Kz','.fig'));

S = struct('E_region',E_region,'DOS',DOS,'NumofStates',NumofStates);
filenamemat=strcat(filename,'.mat')
save(filenamemat,'S');
end;





toc