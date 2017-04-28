%Matlab????????
clear;
fclose('all');
filename='DOS.dat';
fd=fopen(filename);



i=0;

while ~feof(fd)
    i=i+1;
    aline=fgetl(fd)
  
    a=  str2num(aline);
    Eregion(i)=a(1);
    dos(i)=a(2);
    
end;
figure;plot( Eregion,dos); xlim([-0.3,0.3]); %ylim([-0.3,0.3]);xlim([-0.2,0.2]);

fclose(fd);