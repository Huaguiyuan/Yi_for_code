%Matlab????????
fclose('all');
filename='My 0.85eigenvec.dat';
fd=fopen(filename);



i=0;
cband1=zeros(101,2);
while ~feof(fd)
    i=i+1;
    aline=fgetl(fd);
  
    a=  str2num(aline);
    kyregion(i)=a(1);
    cband1(i,:)=a(2:3);
    
end;
figure;plot( kyregion,cband1);ylim([-0.3,0.3]);xlim([-0.2,0.2]);

fclose(fd);