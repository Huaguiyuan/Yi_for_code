	MODULE global
      REAL(kind=8),PARAMETER:: pi=3.141592653589793d0
      COMPLEX(kind=8):: II=(0d0,1d0)
      ! cp: chemical potential
      REAL(kind=8),parameter:: kymin=-0.25d0, kymax=0.25d0
      ! L_max: number of sites along x
      INTEGER ,PARAMETER:: L_max=160,Num_ky=100 !keep even
	! Ncell: number of cells, Dim: dimension of BdG Hamiltonian
      INTEGER ,PARAMETER:: Dim=2*L_max,Ncell=L_max/2
	END MODULE global
    
	Program WeylBDG
    	USE global
    	IMPLICIT NONE
	INTEGER:: i,j,ny
!----------------pairing-------------------------------
     	! bulk Weyl points (kx,ky,kz)=(0,0,+-1), finite size along x, good quantum #s: ky,kz.
	! tp,tm: hopping t+, t-
    	! V: on-site staggered potential
	! function_neg: functions of -ky,-kz
    	REAL(kind=8):: kx,ky,kz,tp,tm,V,delx,dely,cp,tp_neg,tm_neg,V_neg,u,v,uhalf,vhalf,weight,charge
	REAL(kind=8), DIMENSION(:),allocatable :: eval
	COMPLEX (kind=8), DIMENSION(:,:),allocatable :: F	
!XSL
	character*8 string_pbc, string_delx, string_dely, string_cp, string_kstart, string_kend, string_kincrease
	logical pbc
!----------output file control---------
	character (90) :: filename,filenamepara,filenamevec	!data and parameter files
	!ifile: fileloop variable; filemax: maximal # of sections of kz's; 
	integer:: ifile,jfile,Nfilemin,Nfilemax,fileincre,jvec
	logical out_eigenvec
!---------------------------
	ALLOCATE(F(Dim,Dim),eval(Dim))

!XSL
        CALL GETARG(1,string_pbc)
        read(string_pbc,*)pbc

        CALL GETARG(2,string_delx)
        CALL GETARG(3,string_dely)
        read(string_delx,*)delx
        read(string_dely,*)dely

        CALL GETARG(4,string_cp)
        read(string_cp,*)cp

        CALL GETARG(5,string_kstart)
        CALL GETARG(6,string_kend)
        CALL GETARG(7,string_kincrease)
        read(string_kstart,*)Nfilemin
        read(string_kend,*)Nfilemax
        read(string_kincrease,*)fileincre

!	print*, 'input boundary condtion: F - open, T - periodic'
!   	read*, pbc
!	print*, pbc

!	print*, 'input deltax, deltay'
!   	read*,  delx,dely
!	print*, delx,dely
	
!	print*, 'input chemical potential'
!   	read*, cp
!	print*, cp	
	
!	print*, 'input integers only'
!	print*, 'input kz start value, kz end value, kz increment (all in unit 0.1)'
!   	read*,  Nfilemin, Nfilemax, fileincre
!	print*, Nfilemin, Nfilemax, fileincre

    	Write(*,*)  "---- L=",L_max,", Num_ky",Num_ky," ----"
	Write(*,*)  "---- chemical potential cp=",cp,", periodic bc",pbc," ----"
!-------------------------
	write(filenamepara, '(A5,F5.2,A5,F5.2,A2,F5.2,A8)' )'pairx',delx,'pairy',dely,'cp',cp,'para.dat'
 	OPEN(unit=1000,file=filenamepara)
	Write(1000,*) "#periodic bc",pbc
 	Write(1000,*) "#L=",L_max,", Num_ky",Num_ky
	Write(1000,*) "#cp=",cp,", deltax=",delx,", deltay=",dely
	close(1000)
!
!
    	Do ifile=Nfilemin,Nfilemax,fileincre
		kz=0.01d0*ifile
		write (filename,'(A2,F5.2,A4)')'kz',kz,'.dat' 
 		OPEN  (unit=ifile,file=filename)
!!vector start
		jfile=ifile+2000
		write (filenamevec,'(A2,F5.2,A12)') 'kz',kz,'eigenloc.dat'
		open  (unit=jfile,file=filenamevec)
!!vector end
		Do ny=0,Num_ky	!middle point,Num_ky/2
			ky=kymin+ny*(kymax-kymin)/Num_ky
			tp=1d0
			tm=-(ky*ky+kz*kz)
!			V=ky

			V=2d0*ky
			call FillHam(F,cp,tp,tm,V,delx,dely,ky,tp_neg,tm_neg,V_neg,pbc)
			call DiagMatrix(F,Dim,eval)
			Write(ifile,100) ky,eval
!!vector start
			u=0d0
			v=0d0
			uhalf=0d0
			vhalf=0d0
			do jvec=1,Ncell
				u=u+abs(F (jvec,L_max))*abs(F (jvec,L_max))
			end do
			uhalf=u
			do jvec=Ncell+1,2*Ncell
				u=u+abs(F (jvec,L_max))*abs(F (jvec,L_max))
			end do
			do jvec=2*Ncell+1,3*Ncell
				v=v+abs(F (jvec,L_max))*abs(F (jvec,L_max))
			end do
			vhalf=v
			do jvec=3*Ncell+1,4*Ncell
				v=v+abs(F (jvec,L_max))*abs(F (jvec,L_max))
			end do
			weight=uhalf+vhalf
			charge=u-v
			write(jfile,100) ky,eval(L_max),weight,charge,u,v,uhalf,vhalf

			u=0d0
			v=0d0
			uhalf=0d0
			vhalf=0d0
			do jvec=1,Ncell
				u=u+abs(F (jvec,L_max+1))*abs(F (jvec,L_max+1))
			end do
			uhalf=u
			do jvec=Ncell+1,2*Ncell
				u=u+abs(F (jvec,L_max+1))*abs(F (jvec,L_max+1))
			end do
			do jvec=2*Ncell+1,3*Ncell
				v=v+abs(F (jvec,L_max+1))*abs(F (jvec,L_max+1))
			end do
			vhalf=v
			do jvec=3*Ncell+1,4*Ncell
				v=v+abs(F (jvec,L_max+1))*abs(F (jvec,L_max+1))
			end do
			weight=uhalf+vhalf
			charge=u-v
			write(jfile,100) ky,eval(L_max+1),weight,charge,u,v,uhalf,vhalf
!-----------------------
		close(jfile)
!!vector end
		close(ifile)
	End do
!------------------------
 100   Format(1000 (F9.6, 2X))
!	1000: #of repeating patterns, F: float, 9:total words, 6:E-6, 2X: 2 spaces
	Deallocate(F,eval) 
	STOP
	END
!-----------------------------------------------
!-----------------------------------------------
! Inputting BdG Hamiltonian Matrix elements
    subroutine FillHam(F,cp,tp,tm,V,delx,dely,ky,tp_neg,tm_neg,V_neg,pbc)
    use global
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
	end subroutine FillHam
!--------------------------------------------------
!-----------------------------------------------
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
      IF (info.NE.0) then
	    STOP 'Problems with exact diagonalization'
	    end if
      DEALLOCATE (work,rwork,iwork)
      return 
      end      