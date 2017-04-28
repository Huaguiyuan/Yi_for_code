
clear; %smaller 2*Nx *2 BDG matrix:  plot bands structure for pure system

vars;

ny=1;
for ky=KY_region;
HH=Hamiltonian(ky,kz);

[v,r]=eig(HH);r=diag(r); [r,Id]=sort(r); v=v(:,Id);


energy(ny,:)=sort(real(r))';

fny=1;
for s=1:4*NX
  a=abs(v(:,s))'.^2;b=a(1:2*NX)+a(2*NX+1:4*NX);
  a=[sum(b(1:NX)),sum(b(NX+1:2*NX))];
    
 fenergy(ny,fny)=real(r(s));
   
    if (abs(a(1)/a(2))<0.5)
     fny=fny-1;
    end;
    
 fny=fny+1;
end
if abs(ky-0.01)<10^-4;
    [v,rr]=eig(HH);
end;

ny=ny+1;
end;

% for i=NX*2-3:NX+3
% a=abs(v(:,i));
% figure; plot(a);title(num2str(i))
% end

figure; plot(KY_region,energy,'*'); ylim([-0.3,0.3]);xlim([-0.2,0.2]);movegui(gcf,'southwest');

figure;plot(KY_region,fenergy,'*'); ylim([-0.3,0.3]);xlim([-0.2,0.2]);
%subplot(2,1,2); plot(KY_region,energy,'*'); ylim([-0.3,0.3])