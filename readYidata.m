%Matlab????????
clear
fclose('all');
filename='My 0.85.dat';
fd=fopen(filename);



i=0;band1=zeros(101,160);band2=zeros(101,160);
while ~feof(fd)
    i=i+1;
    aline=fgetl(fd);
  
    a=  str2num(aline);
    kyregion(i)=a(1);
    band1(i,:)=a(2:161);band2(i,:)=a(162:321);
    
end;
figure;plot( kyregion,band1,kyregion,band2);ylim([-0.3,0.3]);xlim([-0.2,0.2]);

fclose(fd);