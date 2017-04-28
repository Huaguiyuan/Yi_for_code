program main
use Hamilt
implicit none
integer:: iny


ALLOCATE(F(Dim,Dim),eval(Dim))
eval=0;

print*, 1




call DOSYi();return



! pbc=.false.

!--------------------------
!ALLOCATE(F(Dim,Dim),eval(Dim))
!eval=0;!print*,II,eval
!pbc=.false.;
!delx=0.2;
!dely= -delx*2;
!cp=0.2;

!Nfilemin=103;
!Nfilemax=105;
!fileincre=20;







write(filenamepara, '(A5,F5.2,A5,F5.2,A2,F5.2,A8)' )'pairx',delx,'pairy',dely,'cp',cp,'para.dat'
OPEN(unit=1000,file=filenamepara)
	Write(1000,*) "#periodic bc",pbc
 	Write(1000,*) "#L=",L_max,", Num_ky",Num_ky
	Write(1000,*) "#cp=",cp,", deltax=",delx,", deltay=",dely
	close(1000)



Do ifile=Nfilemin,Nfilemax,fileincre
kz=0.01d0*ifile
print*,kz
write (filename,'(A2,F5.2,A4)')'MYkz',kz,'.dat'

OPEN  (unit=ifile,file=filename)



!!vector start
jfile=ifile+2000
write (filenamevec,'(A2,F5.2,A12)') 'MYkz',kz,'eigenvec.dat'
open  (unit=jfile,file=filenamevec);
!!vector end


Do iny  = 0,Num_ky
    ky=kymin+iny*(kymax-kymin)/Num_ky
    tp=1d0
    tm=-(ky*ky+kz*kz)
    V=2d0*ky


call FillHam(F,cp,tp,tm,V,delx,dely,ky,tp_neg,tm_neg,V_neg,pbc)
call DiagMatrix(F,Dim,eval)
Write(ifile,100) ky,eval
!write(jfile,*) '============================='
!write(jfile,*) 'ky=',ky,', mideval1=',eval(L_max),', mideval2=',eval(L_max+1)
write(jfile,*) ky,eval(L_max),eval(L_max+1)
!write(jfile,*) '---------|mid evec1|, Dim=', Dim,'---------'
do jvec=1,Dim
!write(jfile,*) jvec,abs(F (jvec,L_max))
end do
!write(jfile,*) '---------|mid evec2|, Dim=', Dim,'---------'
do jvec=1,Dim
!write(jfile,*) jvec,abs(F (jvec,L_max+1))
end do
end do
!-----------------------
close(jfile)
close(ifile)
End do


!------------------------
100   Format(1000 (F9.6, 2X))


end program



