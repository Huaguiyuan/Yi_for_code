module Hamilt
implicit none
 REAL(kind=8),PARAMETER:: pi=3.141592653589793d0
      COMPLEX(kind=8):: II=(0d0,1d0)
      ! cp: chemical potential
      REAL(kind=8),parameter:: kymin=-0.25d0, kymax=0.25d0
      ! L_max: number of sites along x
      INTEGER ,PARAMETER:: L_max=160,Num_ky=100 !keep even
	! Ncell: number of cells, Dim: dimension of BdG Hamiltonian
      INTEGER ,PARAMETER:: Dim=2*L_max,Ncell=L_max/2


!----------------pairing-------------------------------
     	! bulk Weyl points (kx,ky,kz)=(0,0,+-1), finite size along x, good quantum #s: ky,kz.
	! tp,tm: hopping t+, t-
    	! V: on-site staggered potential
	! function_neg: functions of -ky,-kz
    	!REAL(kind=8):: kx,ky,kz,tp,tm,V,delx,dely,cp,tp_neg,tm_neg,V_neg
	REAL(kind=8), DIMENSION(:),allocatable :: eval
	COMPLEX (kind=8), DIMENSION(:,:),allocatable :: F	
!XSL
	character*8 string_pbc, string_delx, string_dely, string_cp, string_kstart, string_kend, string_kincrease

!----------output file control---------
	character (90) :: filename,filenamepara,filenamevec	!data and parameter files
	!ifile: fileloop variable; filemax: maximal # of sections of kz's; 
	integer:: ifile,jfile,jvec
	logical out_eigenvec


REAL(kind=8)::kx,ky,kz,tp,tm,V
REAL(kind=8),parameter::tp_neg=0,tm_neg=0,V_neg=0
real(kind=8),parameter::delx=0.2d0,dely=0.4d0,cp=0.2d0,U_disorder=0;
logical,parameter:: pbc=.true.
integer,parameter::Nfilemin=85,Nfilemax=105,fileincre=40


!    integer::lda,ipiv(dim),info,lwork
!    parameter(lwork=dim*10)
!    double complex work(NN8*10),work1(NN4*10)

! DOS calculation
  INTEGER ,PARAMETER:: NY=10;
  INTEGER ,PARAMETER:: NZ=10;
  INTEGER ,PARAMETER:: NE=1000;
  real(kind=8),parameter::R_shift=0
  real(kind=8),parameter::KKymin=0.35, Emin=0.3,eta=0.0001






!print*,II,eval



!Nfilemin=103;
!Nfilemax=105;
!fileincre=20;


!
!kz1=sqrt(1-cp+Rshift);
!kz2=sqrt(1+cp-Rshift);
!kzmax=max(abs(kz1),abs(kz2));
!
!kzrr=0.2;
!ss=(kzmax+kzrr-(-kzmax-kzrr))/NZ;
!KZ_region=-kzmax-kzrr:ss:kzmax+kzrr-ss;
!
!k1max=0.5;ak2max=1.5; ss=(k2max-k1max)/NZ; KZ_region=k1max:ss:k2max-ss;
!
!KY_region=-Kymin:2*Kymin/NY:Kymin-2*Kymin/NY;
!
!E_region=-Emin:2*Emin/NE:Emin-2*Emin/NE;



contains

!Inputting BdG Hamiltonian Matrix elements
subroutine FillHam(F,cp,tp,tm,V,delx,dely,ky,tp_neg,tm_neg,V_neg,pbc)
!use global
implicit none
complex (kind=8):: F(Dim,Dim)
REAL(kind=8):: cp,tp,tm,V,delx,dely,ky,tp_neg,tm_neg,V_neg
integer :: i,j
logical pbc
! ZERO
F=0d0
!------------------------------
! onsite potential and chemical potential
Do i=1,Ncell
F(2*i-1,2*i-1)			= V-cp	!-V-cp
F(2*i,2*i)				=-V-cp	!V-cp
F(L_max+2*i-1,L_max+2*i-1)	= V+cp	!V_neg+cp
F(L_max+2*i,L_max+2*i)	=-V+cp	!-V_neg+cp
Enddo
!------------------------------
!hopping along x
Do i=1,Ncell
F(	  2*i-1,		2*i)	 	= tm
F(	  2*i,		2*i-1)	= tm
F(L_max+2*i-1,	L_max+2*i)	 	=-tm
F(L_max+2*i,	L_max+2*i-1)	=-tm
!
If (i.lt.Ncell) then
F(	  2*i+1,		2*i)		= tp
F(	  2*i,		2*i+1)	= tp
F(L_max+2*i+1,  L_max+2*i)		=-tp
F(L_max+2*i,    L_max+2*i+1)	=-tp
Elseif (pbc) then
F(	  2*i,	      1)		= tp
F(	  1,		      2*i)		= tp
F(L_max+2*i,	L_max+1)		=-tp
F(L_max+1,		L_max+2*i)		=-tp
Endif
!
Enddo
!-----------------------------
!Pairing along y
Do i=1,Ncell
F(      2*i-1,	L_max+2*i  )	=-2d0*dely*dsin(ky)
F(      2*i,	L_max+2*i-1)	=-2d0*dely*dsin(ky)
F(L_max+2*i,          2*i-1)	=-2d0*dely*dsin(ky)
F(L_max+2*i-1,        2*i  )	=-2d0*dely*dsin(ky)
Enddo
!-----------------------------
!Pairing along x
Do i=1,Ncell-1
F(	  2*i-1,	L_max+2*i+1)	= delx
F(	  2*i+1,	L_max+2*i-1)	=-delx
F(	  2*i,	L_max+2*i+2)	= delx
F(	  2*i+2,	L_max+2*i  )	=-delx

F(L_max+2*i+1,	  2*i-1	)	= delx
F(L_max+2*i-1,	  2*i+1	)	=-delx
F(L_max+2*i+2,	  2*i   	)	= delx
F(L_max+2*i,	  2*i+2	)	=-delx
Enddo
If (pbc) then
F( 2*Ncell-1,	L_max+1  )	= delx
F( L_max+1,	2*Ncell-1)	= delx

F(	  1,	     2*Ncell-1+L_max)	=-delx
F(2*Ncell-1+L_max,	  1      )	=-delx

F( 2*Ncell  ,	L_max+2)	= delx
F( L_max+2  ,	2*Ncell)	= delx

F(	  2,		2*Ncell+L_max  )	=-delx
F(2*Ncell+L_max,		2  )		=-delx
Endif
!------------------------------

!-----------------------------------------------
end subroutine FillHam
!--------------------------------------------------

subroutine  add_disorder(F)
implicit none
complex (kind=8):: F(Dim,Dim)
real(kind=8)::r1,r2
integer::i
call random_seed()




Do i=1,Ncell
call random_number(r1);
call random_number(r2);
r1=(r1-0.5)*U_disorder;
r2=(r2-0.5)*U_disorder

F(	  2*i-1,		2*i-1)	 	= r1+ F(	  2*i-1,		2*i-1)
F(	  2*i,		2*i)	        = r2+ F(	  2*i,		2*i)
F(L_max+2*i,	L_max+2*i)	 	=-r1+ F(L_max+2*i,	L_max+2*i)
F(L_max+2*i-1,	L_max+2*i-1)	=-r2+ F(L_max+2*i-1,	L_max+2*i-1)
enddo
end subroutine


!---------------------------------------------
SUBROUTINE DiagMatrix(Matrix,Dimen,Eval)
implicit none
CHARACTER jobz,uplo
integer Dimen
integer enne,N,LDA,lwork,lrwork,liwork,info
COMPLEX (KIND=8)::Matrix(Dimen,Dimen)
DOUBLE PRECISION, dimension(Dimen) ::Eval(Dimen)
COMPLEX (kind=8), dimension(:), allocatable :: work
DOUBLE PRECISION, dimension(:), allocatable :: rwork
integer, dimension(:), allocatable :: iwork
enne=size(Matrix,1)
!     Set jobz='N' in order to compute eigenvalues only.
!     Set jobz='V' in order to compute eigenvalues and eigenvectors.
jobz='V'
!     Set uplo'=U' in order to give in input the upper triangle of the matrix
uplo='U'
N=enne
lda=enne
!     Dimensions of the workspaces
lwork=2*N+N*N
lrwork=1+5*N+2*N*N
liwork=3+5*N
ALLOCATE (work(lwork))
ALLOCATE (rwork(lrwork))
ALLOCATE (iwork(liwork))
eval=0d0
call zheevd(jobz,uplo,N,Matrix,LDA,Eval,work,lwork,rwork,lrwork,iwork,liwork,info)
!call ch(nm,n,ar,ai,w,matz,zr,zi,fv1,fv2,fm1,ierr)

!call zgetrf(Dimen, Dimen, Matrix, Dimen, ipiv, info)
!call zgetri(Dimen,Matrix,Dimen,ipiv,work,lwork,info )

IF (info.NE.0) then
STOP 'Problems with exact diagonalization'
end if
DEALLOCATE (work,rwork,iwork)
return
end subroutine


subroutine DOSYi()
implicit none
real(kind=8)::kz1,kz2,kzmax,kzrr,k1max,k2max,ss
integer::iky,ikz,iE,ir,ia,ib,count
real(kind=8)::KY_region(NY),E_region(NE),KZ_region(NZ),DOS(NE),sss,Fss(Dim,Dim)

!REAL(kind=8):: kx,ky,kz,tp,tm,V,tp_neg,tm_neg,V_neg

!KZ_region=-kzmax-kzrr:ss:kzmax+kzrr-ss;
kz1=sqrt(1-cp+R_shift);
kz2=sqrt(1+cp-R_shift);
kzmax=max(abs(kz1),abs(kz2));
kzrr=0.2;
k2max=kzmax+kzrr;k1max=-k2max;
k1max=0.5;k2max=1.5;  !for R_shift=0, with two fermi surfaces case
ss=(k2max-k1max)/NZ;

do ikz=1,NZ;
KZ_region(ikz)=k1max+(ikz-1)*ss;
enddo


!KY_region=-Kkymin:2*Kkymin/NY:Kkymin-2*Kkymin/NY;
ss=2*Kkymin/NY
do iky=1,NY;
KY_region(iky)=-Kkymin+(iky-1)*ss;
enddo

!E_region=-Emin:2*Emin/NE:Emin-2*Emin/NE;

ss=2*Emin/NE;
do iE=1,NE
E_region(iE)=-Emin+(iE-1)*ss;
enddo


DOS=0;count=0;

!LDOS_left=zeros(size(E_region));
!LDOS_right=zeros(size(E_region));

Do ikz=1,NZ
kz=KZ_region(ikz);

Do iky = 1,NY
ky=KY_region(iky);

count=count+1

tp=1d0
tm=-(ky*ky+kz*kz)
V=2d0*ky;


call FillHam(F,cp,tp,tm,V,delx,dely,ky,tp_neg,tm_neg,V_neg,pbc);
call add_disorder(F)

!fss=abs(F);
!sss=0;do ia=1,dim;do ib=1,dim;sss=sss+abs(f(ia,ib));enddo;enddo;
!print*,'cp=',cp,'tp=',tp,'tm=',tm,'V=',V,'delx',delx,'dely=',dely,'ky=',ky,'kz=',kz,'sss=', sss,'pbc=',pbc
!stop

call DiagMatrix(F,Dim,eval)



do ir=1,Dim
ss=eval(ir);
if ((ss>-Emin).and.(ss<Emin))then
do iE=1,NE;
DOS(iE)=DOS(iE)+eta*((E_region(iE)-ss)**2+eta**2)**(-1)/pi;
enddo
endif

enddo


!Write(ifile,100) ky,eval
!write(jfile,*) '============================='
!write(jfile,*) 'ky=',ky,', mideval1=',eval(L_max),', mideval2=',eval(L_max+1)
!write(jfile,*) ky,eval(L_max),eval(L_max+1)
!write(jfile,*) '---------|mid evec1|, Dim=', Dim,'---------'

!write(jfile,*) '---------|mid evec2|, Dim=', Dim,'---------'

end do
!-----------------------
End do

ifile=1001;
write (filename,'(A7)')'DOS.dat'

OPEN  (unit=ifile,file=filename)

do iE=1,NE
Write(ifile,*)E_region(iE) ,DOS(iE)
enddo

close(ifile)

print*,'count=',count
end subroutine





end module


