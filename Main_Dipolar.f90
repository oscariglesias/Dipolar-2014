!_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
!	- Simulation of Dipolar fields generated in an assembly of NP placed in a regular latticecor random positions.
!	- Energy terms are: 
!		1) Uniaxial anisotropy. 
!		2) Zeeman energy, 
!		3) Dipolar interactions between all spins.
!		y = a2/(sqrt(2*pi)*a0)*exp(-(x-a1)^2/(2*a0^2))
!		y = a2/(sqrt(2*pi)*a0*x)*exp(-(ln(x)-a1)^2/(2*a0^2))

!_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
!	Compilation with g95:
!	->	g95 Main.f90 -o Main 
!	->	g95 Main.f90 -o Main -L/cygdrive/c/cygwin/release/lib/ -lfftw3 -lm  (for FFTW)
!	->   pgf95  -fast -O4 -o Main Main.f90
!	->	ifort -fast -O4 -o Main_Dipolar Main_Dipolar.f90
! 	-> Configuration files must be translated to Mathematica input format file *.nb 
!		with the program Tranconf_Dip2014.f90
!
! 	-> Files with dipolar fields must be translated to Mathematica input format file *.nb 
!		with the program Tranconf_Hdip2014.f90
!-----------------------------------------------------------------------
!	- To be completed:
!		o) 
!-----------------------------------------------------------------------
!	Author: Òscar Iglesias, UB, March 2014
!-----------------------------------------------------------------------
      include 'Common_Dipolar.h'
      include 'SUBROUTINES_Dipolar.f90'
!      include 'Sorting_2.f90'
      program DIPOLAR2014
      use COMMON_DIPOLAR
      use SUBROUTINES_DIPOLAR
!		use QSORT_C_MODULE
		
      implicit none
!-------------------- Internal variables -----------------------
      integer :: Ngmax
      integer :: i,k,itemp,ifield,Nh,Ntemp,imcs
      integer :: acc_i
      real :: acc_t
      real*8 :: Magn_av(DimeS),Ener_av(DimeE)
      real,dimension(Dime) :: sumx 
      character*80 :: filedat
      character(len=mediumchar) :: Name,Char_ifield			

!      call rmarin(176,317)
!      call rmarin(16,3177)
      call rmarin(11231,4522)

      filedat='Main_Dipolar.dat'
      print *,filedat
      call READDATA(filedat,Ngmax)
		Hext=Hini
		Hdip_fac=mu0*Msat      
      call INIT_VOLUMES
      call LATTICE
      print *,'Entering INIT_CONF'
      call INIT_CONF
      call INIT_ANIS
		if (IcompW.eq."NoDipolar") then
		elseif (trim(Net).eq.'SC_1D') then
			call WS_1D(Ngmax)
		elseif (trim(Net).eq.'SC_2D') then
			call WS_2D(Ngmax)
		else
			if (trim(ISimW).eq.'SimW') then
				print *,'Entering WS_TOP'
				call WS_TOP(Ngmax)
			elseif (trim(ISimW).eq.'NoSimW') then
				print *,'Entering WS'
				call WS(Ngmax)
			endif
		endif
		print *,'Entering HDIPOLAR'
		if (IcompW.ne."NoDipolar") then
			call HDIPOLAR
	      call HDIPSELEC
			call WRITEHDIP('_Hdipini',8)
			call WRITEHDIPSELEC('_HdipiniSel',11)			
			
			print '(1x,A,3(2x,e12.6))','Mean dipolar field modulus: ', sum(Hmod)/Nspin,Minval(Hmod),Maxval(Hmod)
			print '(1x,A,3(2x,e12.6))','Mean dipolar field modulus: ', sum(Hdip_fac*Vol*Hmod)/Nspin,Minval(Hdip_fac*Vol*Hmod),&
	  		Maxval(Hdip_fac*Vol*Hmod)
			print '(1x,A,3(2x,e12.6))','Mean dipolar field components: ', sum(Hdip,dim=1)/Nspin
			print '(1x,A,3(2x,e12.6))','Mean dipolar field selec modulus: ', sum(Hmod_Sel)/Nspin_Sel,Minval(Hmod_Sel),Maxval(Hmod_Sel)
			print '(1x,A,3(2x,e12.6))','Mean dipolar field selec components: ', sum(Hdip_Sel,dim=1)/Nspin_Sel
	      open(3,file=trim(FileHdip)//'_MeanHDip'//'.out', status='unknown')
	      open(33,file=trim(FileHdip)//'_MeanHDipSel'//'.out', status='unknown')
			write(3,22) Hini,sum(Hmod)/Nspin,Minval(Hmod),Maxval(Hmod)				
			write(33,22) Hini,sum(Hmod_Sel)/Nspin,Minval(Hmod_Sel),Maxval(Hmod_Sel)

			call WRITEHISTOINI(Hmod,Hdip,Nspin,'_Histoini.out',13)
			call WRITEHISTOINI(Hmod_Sel,Hdip_Sel,Nspin_Sel,'_HistoSelini.out',16)
			deallocate (Hdip_Sel,Hmod_Sel)
		endif

		print *,'Writing initial configuration '
		call WRITECONF('_Confini',8)
		print *,'Writing diameter histogram '
		if (trim(Ivolume).eq.'Logn') then
			allocate (Hist(Nbars),Bin(Nbars))
			call HISTO(Diam(:),'Diam',Nspin,Nbars,Bin,Hist)
			open(1,file=trim(FileHdip)//'_HistoDiam.out')
				do i=1,Nbars		
					write(1,*) Bin(i),Hist(i)
				enddo
			close(1)
			call HISTO(Ka(:),'Eba ',Nspin,Nbars,Bin,Hist)
			open(1,file=trim(FileHdip)//'_HistoEb0.out')
				do i=1,Nbars		
					write(1,*) Bin(i),Hist(i)
				enddo
			close(1)
			deallocate (Hist,Bin)
		endif

      call EBARRIER
		call WRITEHISTOEBINI(Eb,Nspin,'_HistoEbini.out',15)
		stop

		call TOTALENERGY
				
		Nh=abs(int((Hend-Hini))/abs(Dh))  ! Number of fields
		print *,'Nfileds=',Nh
		Ntemp=nint((Tempi-Tempf)/Tcoef)+1
      Temp=Tempi+Tcoef

      open(2,file=trim(FileHdip)//'_MH'//'.out', status='unknown')
		do itemp=1,Ntemp
			Temp=Temp-Tcoef
			print *,'Ntemp=',Ntemp,'   T= ',Temp,'   Thermal energy= ',kB*Temp
			if (Temp.gt.0.0) then
				if (Units.eq."RealUnits") then
					Beta=1./kB/Temp
				else
					Beta=1./Temp 
				endif
			else
				Beta=1.
			endif
			
			open(4,file=trim(FileHdip)//'_HistoMod'//'.out')
			open(5,file=trim(FileHdip)//'_HistoHz'//'.out')
			open(7,file=trim(FileHdip)//'_HistoHy'//'.out')
			open(8,file=trim(FileHdip)//'_HistoHx'//'.out')
			open(44,file=trim(FileHdip)//'_HistoSelMod'//'.out')
			open(55,file=trim(FileHdip)//'_HistoSelHz'//'.out')
			open(77,file=trim(FileHdip)//'_HistoSelHy'//'.out')
			open(88,file=trim(FileHdip)//'_HistoSelHx'//'.out')
			do ifield=0,Nh
	         Hext=Hini+Dh*ifield
	         print *,'Hext= ',Hext
				acc_t=0.
				print *,'Starting MONTE '
				Magn_av=0.
				Ener_av=0.
				Ener(3)=0.
				do i=1,Nspin
		         Ener(3)=Ener(3)-dot_product(H_n(:),S(i,:))*Vol(i) ! Zeeman energy
				enddo
				Ener(3)=Hext*Msat*LattConstant**3*Ener(3)
				Ener(4)=Ener(1)+Ener(2)+Ener(3)

				do imcs=1,Nterma
					call MONTE(imcs,acc_i)
!!					print *,'MCS, Magn, Ener= ',imcs,Magn(3)/Nspin,Ener(4)/Nspin
					Ener_av=Ener_av+Ener/Nspin
					Magn_av=Magn_av+Magn/Nspin
					acc_t=acc_t+acc_i
				enddo
				print *,'End of MONTE, Acc rate= ',acc_t/Nterma/Nspin*100,' %'
				print *,'Hext, Energy, Magn = ',Hext,Ener_av(4)/Nterma,Magn_av(3)/Nterma
				write(2,22) Hext,Magn_av/Nterma
				
!				write(unit=Char_ifield,fmt=*)ifield
				
				if (mod(ifield,Ihist).eq.0) call WRITEHISTO(ifield)
				write(3,22) Hext,sum(Hmod)/Nspin,Minval(Hmod),Maxval(Hmod)				
		      call HDIPSELEC
				if (mod(ifield,Ihist).eq.0) call WRITEHISTOSELEC(ifield)
				write(33,22) Hext,sum(Hmod_Sel)/Nspin_Sel,Minval(Hmod_Sel),Maxval(Hmod_Sel)				
				deallocate (Hdip_Sel,Hmod_Sel)
			enddo		! End of field loop
		enddo			! End of temperature loop
		close(2)
		close(3)
		close(33)
		close(4)
		close(5)
		close(7)
		close(8)
		close(44)
		close(55)
		close(77)
		close(88)

22   format(6(1xd12.6))

      end program DIPOLAR2014
