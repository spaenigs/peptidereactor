c     Predicting secondary structure
c	Eshel Faraggi -- Apr. 2007
c	21 window 
c	g77 -ffixed-line-length-none -Wall -W -Wsurprising  -O2 -o phipsi_ss1.e phipsi_ss1.f
c	ifort -132 -watch all -warn all -O2 -o phipsi_ss1.e phipsi_ss1.f
c	ifort -132 -O2 -o phipsi_ss1.e phipsi_ss1.f
c	ifort -132 -g -o phipsi_ss1.e phipsi_ss1.f
      program phipsi_ss1
      implicit real*8(a-h,o-z)
      PARAMETER(nrpx=8000,nrdx=2000000,nrx=21,nip=30,nw1=101,nw2=101,nrxo=1,numav=5,nfilt=21)

      common/weights/ wi(nrx+1,nip,nw1),wh(nw1+1,nw2),wo2(nw2+1,nrxo,3)
      common/filtwei/ fwi(nrx+1,3,nfilt), fwo(nfilt+1,nrxo,3)
      common/filin/ protin(nrdx,nip+6), iprotin(nrdx,3)
      character phia*1
      common/filout/ phia(nrdx)
      common/inout/ resv(nrx+1,nip)
      common/neurons/ hnu1(nw1+1),hnu2(nw2+1),onu2(nrdx,3),hnuf(nfilt+1), finp(nrx,3)
      common/edges/ ibegia(nrdx),iendia(nrdx),ibegoa(nrdx),iendoa(nrdx), nisf, nosf
      common/dcays/ dij1(nrx,nw1),dij2(nw2,nrxo),dfac,dijh(nw1,nw2),dijfi(nrx,nfilt),dijfo(nfilt,nrxo)
      common/misc/ alpha
      common/afilt/ onu2p(nrdx,3,numav)
      character seqres(nrdx)
      common seqres
      common/probss/ tss(nrdx,3)
      character rtrainf(nrpx)*150
      common/sinout/ ntf,nta,rtrainf

      character*150 spec, spxdir, workdir, listfl, flout

      call getarg(1,spec)
      call getarg(2,spxdir)
      call getarg(3,workdir)
      call getarg(4,listfl)
c      read(spec,'(I20)') inum

c* Definition of the general input/output files
      call getdefs(listfl)
      call netinit()
      call setio(spec,workdir)	! initialize input/output
C      "-------------------Run 1-------------------"
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      do 287 iav=1,numav
        call readwei(0,iav,spxdir)
        do ix=1,nta
            call readwin(ix,0)
            call calcnet(ix,0,iav)
        enddo

C      "-------------------Filter-------------------"
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        call readwei(10,iav,spxdir)
        do ix=1,nta
            call readwin(ix,10)
            call calcnet(ix,10,iav)
        enddo
287   continue

C      "-------------------Average-------------------"
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      call getangs()

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCC--- OUTPUT PREDICTION  ---CCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      flout = trim(workdir)//"/out_ss1"
      open(unit=99,file=flout,status='unknown')
      do i=1,nta
        write(99,499) i,seqres(i),phia(i),tss(i,1),tss(i,2),tss(i,3)
      enddo
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
499   format(I9,4x,A1,3x,A1,2x,F8.5,2x,F8.5,2x,F8.5)
      stop
      end
c_______________________________________________________________________
      
c------------- subroutine getdefs-------------------------------
c--
c* Get initial definitions
c--
      subroutine getdefs(listfl)
      implicit real*8(a-h,o-z)
      PARAMETER(nrpx=8000,nrdx=2000000,nrx=21,nip=30,nw1=101,nw2=101,nrxo=1,numav=5,nfilt=21)
      common/dcays/ dij1(nrx,nw1),dij2(nw2,nrxo),dfac,dijh(nw1,nw2),dijfi(nrx,nfilt),dijfo(nfilt,nrxo)
      character rtrainf(nrpx)*150
      common/sinout/ ntf,nta,rtrainf
      common/misc/ alpha

      character listfl*150

c* Definitions
c*********** alpha -- parameter for transfer function
      alpha = 0.2d0
c*********** dfac -- parameter for distance decay matrix
      dfac = 2.d0

      open(unit=3,file=listfl,status='old')
c* read names of files containing training proteins
      i = 1
      do 245 while (i.gt.0)
        read(3,'(A)',end=246)rtrainf(i)
	i = i + 1
245   continue
246   continue
      close(3)
      ntf = i-1

      end
c___________________________________________________________________

c------------- subroutine netinit -------------------------------
c--
c* Initialize the neural network (biases and guiding factors)
c--
      subroutine netinit()
      implicit real*8(a-h,o-z)
      PARAMETER(nrpx=8000,nrdx=2000000,nrx=21,nip=30,nw1=101,nw2=101,nrxo=1,numav=5,nfilt=21)
      common/inout/ resv(nrx+1,nip)
      common/neurons/ hnu1(nw1+1),hnu2(nw2+1),onu2(nrdx,3),hnuf(nfilt+1), finp(nrx,3)
      common/dcays/ dij1(nrx,nw1),dij2(nw2,nrxo),dfac,dijh(nw1,nw2),dijfi(nrx,nfilt),dijfo(nfilt,nrxo)
      
c***   setting up the biases
      resv(nrx+1,1) = 1.d0
      hnu2(nw2+1) = 1.d0
      hnu1(nw1+1) = 1.d0
      hnuf(nfilt+1) = 1.d0

        delhi = (nrx-1.d0)/(nw1-1.d0)
c***    initialize distance matrix
	do i=1,nrx
          do k=1,nw1
	    dij1(i,k) = dfac / dsqrt(1.d0+((k-1)*delhi-i+1)*
     *                          ((k-1)*delhi-i+1))
	  enddo
	enddo

        jc = (nw1+1)/2
        jc2 = (nw2+1)/2
        mc = (nrxo+1)/2
	do i=1,nw1
          do j=1,nw2
	    dijh(i,j) = dfac / dsqrt(1.d0+delhi*delhi*(i-jc-j+jc2)*
     *                          (i-jc-j+jc2))
          enddo
	enddo

	do i=1,nw2
          do j=1,nrxo
	    dij2(i,j) = dfac / dsqrt(1.d0+((i-jc)*delhi-j+mc)*
     *                          ((i-jc)*delhi-j+mc))
          enddo
	enddo

        delhi = (nrx-1.d0)/(nfilt-1.d0)
c***    initialize distance matrix
	do i=1,nrx
          do k=1,nfilt
	    dijfi(i,k) = dfac / dsqrt(1.d0+((k-1)*delhi-i+1)*
     *                          ((k-1)*delhi-i+1))
	  enddo
	enddo

        jc = (nfilt+1)/2
        mc = (nrxo+1)/2
	do i=1,nfilt
          do j=1,nrxo
	    dijfo(i,j) = dfac / dsqrt(1.d0+((i-jc)*delhi-j+mc)*
     *                          ((i-jc)*delhi-j+mc))
          enddo
	enddo

      end
c___________________________________________________________________

c------------- subroutine setio -------------------------------
c--
c* Setup the input/output arrays
c--
      subroutine setio(ftfl,workdir)
      implicit real*8(a-h,o-z)
      PARAMETER(nrpx=8000,nrdx=2000000,nrx=21,nip=30,nw1=101,nw2=101,nrxo=1,numav=5,nfilt=21)
      common/filin/ protin(nrdx,nip+6), iprotin(nrdx,3)
      common/edges/ ibegia(nrdx),iendia(nrdx),ibegoa(nrdx),iendoa(nrdx), nisf, nosf
      character seqres(nrdx)
      common seqres
      character rtrainf(nrpx)*150
      common/sinout/ ntf,nta,rtrainf
      
      character ftfl*150, inflnm*150, workdir*150, tmpchr*1

      double precision tprot(53)

c* midpoint for input/output windows
      nisf = (nrx+1)/2
      nosf = (nrxo+1)/2

      inflnm = trim(workdir)//"/out_ss0"
      open(unit=44,file=inflnm,status='old')
      inflnm = trim(workdir)//"/out_rsa0"
      open(unit=45,file=inflnm,status='old')
      inflnm = trim(workdir)//"/out_phipsi0"
      open(unit=46,file=inflnm,status='old')
      nta = 0
        i = 1
        do 265 while (i.gt.0)
          nta = nta + 1
          read(44,*,end=266) ntmp, tmpchr, tmpchr,atmp, atmp, atmp, (tprot(j),j=1,27),
     *      ibegia(nta),iendia(nta),ibegoa(nta),iendoa(nta)
          read(45,*,end=266) ntmp, tmpchr, atmp, protin(nta,nip)
          read(46,*,end=266) ntmp, tmpchr, atmp, atmp, protin(nta,nip+1), protin(nta,nip+2)
          seqres(nta) = tmpchr
          do j=1,27
            protin(nta,j) = tprot(j)
          enddo
265     continue
266     continue
      close(45)
      close(46)
      nta = nta - 1

      end
c___________________________________________________________________

c------------- subroutine readwei -------------------------------
c--
c* Read optimized weights
c--
      subroutine readwei(ishift,iav2,spxdir)
      implicit real*8(a-h,o-z)
      PARAMETER(nrpx=8000,nrdx=2000000,nrx=21,nip=30,nw1=101,nw2=101,nrxo=1,numav=5,nfilt=21)
      common/weights/ wi(nrx+1,nip,nw1),wh(nw1+1,nw2),wo2(nw2+1,nrxo,3)
      common/filtwei/ fwi(nrx+1,3,nfilt), fwo(nfilt+1,nrxo,3)
      
      character weifl*150,weishft*2,avo,weidir*150, spxdir*150

      weidir = trim(spxdir)//"/weights/ss1/"
      if (ishift.eq.0) then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCC--- INPUT WEIGHTS  ---CCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        itmp = 77
        write(weishft,'(I2)') itmp
        write(avo,'(I1)') iav2
        weifl = trim(weidir)//"/av"//avo//"/fort."//weishft
        open(unit=77,file=weifl,status='unknown')
        do i=1,nrx+1
          do j=1,nip
            do k=1,nw1
              read(77,*) wi(i,j,k)
            enddo
          enddo
        enddo
        close(77)
        itmp = itmp + 1
        write(weishft,'(I2)') itmp
        weifl = trim(weidir)//"/av"//avo//"/fort."//weishft
        open(unit=77,file=weifl,status='unknown')
        do i=1,nw1+1
          do j=1,nw2
            read(77,*) wh(i,j)
          enddo
        enddo
        close(77)
        itmp = itmp + 1
        write(weishft,'(I2)') itmp
        weifl = trim(weidir)//"/av"//avo//"/fort."//weishft
        open(unit=77,file=weifl,status='unknown')
        do i=1,nw2+1
          do j=1,nrxo
            read(77,*) (wo2(i,j,k),k=1,3)
          enddo
        enddo
        close(77)
      else
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCC--- FILTER WEIGHTS  ---CCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        itmp = 87
        write(weishft,'(I2)') itmp
        write(avo,'(I1)') iav2
        weifl = trim(weidir)//"/av"//avo//"/fort."//weishft
        open(unit=77,file=weifl,status='unknown')
        do i=1,nrx+1
          do j=1,3
            do k=1,nfilt
              read(77,*) fwi(i,j,k)
            enddo
          enddo
        enddo
        close(77)
        itmp = itmp + 1
        write(weishft,'(I2)') itmp
        weifl = trim(weidir)//"/av"//avo//"/fort."//weishft
        open(unit=77,file=weifl,status='unknown')
        do i=1,nfilt+1
          do j=1,3
            read(77,*) fwo(i,1,j)
          enddo
        enddo
        close(77)
      endif

      end
c___________________________________________________________________

c------------- subroutine readwin -------------------------------
c--
c* Read the current protein window
c--
      subroutine readwin(iwin,itp)
      implicit real*8(a-h,o-z)
      PARAMETER(nrpx=8000,nrdx=2000000,nrx=21,nip=30,nw1=101,nw2=101,nrxo=1,numav=5,nfilt=21)
      common/filin/ protin(nrdx,nip+6), iprotin(nrdx,3)
      common/inout/ resv(nrx+1,nip)
      common/neurons/ hnu1(nw1+1),hnu2(nw2+1),onu2(nrdx,3),hnuf(nfilt+1), finp(nrx,3)
      common/edges/ ibegia(nrdx),iendia(nrdx),ibegoa(nrdx),iendoa(nrdx), nisf, nosf
     
      integer iwin,itp

      if (itp.le.0) then
        do i=ibegia(iwin),iendia(iwin)
          irw = iwin-nisf+i
          resv(i,1) = protin(irw,nip+1)
          resv(i,2) = protin(irw,nip+2)
          do j=1,20
            resv(i,j+2) = protin(irw,j)
          enddo
          do j=1,7
            resv(i,j+22) = protin(irw,j+20)
          enddo
          resv(i,30) = protin(irw,nip)
        enddo
      else
        do i=ibegia(iwin),iendia(iwin)
          irw = iwin-nisf+i
          do j=1,3
            finp(i,j) = onu2(irw,j)
          enddo
        enddo
      endif
      
      end
c___________________________________________________________________

c------------- subroutine calcnet -------------------------------
c--
c* Calculate the neural network's output angle given input residue
c* sequence parameters in resv matrix
c--
      subroutine calcnet(iwin,itp,iav2)
      implicit real*8(a-h,o-z)
      PARAMETER(nrpx=8000,nrdx=2000000,nrx=21,nip=30,nw1=101,nw2=101,nrxo=1,numav=5,nfilt=21)

      common/weights/ wi(nrx+1,nip,nw1),wh(nw1+1,nw2),wo2(nw2+1,nrxo,3)
      common/filtwei/ fwi(nrx+1,3,nfilt), fwo(nfilt+1,nrxo,3)
      character phia*1
      common/filout/ phia(nrdx)
      common/inout/ resv(nrx+1,nip)
      common/neurons/ hnu1(nw1+1),hnu2(nw2+1),onu2(nrdx,3),hnuf(nfilt+1), finp(nrx,3)
      common/edges/ ibegia(nrdx),iendia(nrdx),ibegoa(nrdx),iendoa(nrdx), nisf, nosf
      common/dcays/ dij1(nrx,nw1),dij2(nw2,nrxo),dfac,dijh(nw1,nw2),dijfi(nrx,nfilt),dijfo(nfilt,nrxo)
      common/afilt/ onu2p(nrdx,3,numav)
      integer iwin,itp,iav2
      double precision tfunc !transfer function

      if (itp.le.0) then
c** Calculate nueron values of the first hidden layer  ---  NN 1
        do i=1,nw1
          ztemp=0.d0
          do j=ibegia(iwin),iendia(iwin)
            do k=1,nip
              ztemp = ztemp + resv(j,k) * wi(j,k,i) * dij1(j,i)
            enddo
          enddo
          ztemp = ztemp + resv(nrx+1,1) * wi(nrx+1,1,i)
          hnu1(i) = tfunc(ztemp)
        enddo

        do i=1,nw2
          ztemp=0.d0
          do j=1,nw1
            ztemp = ztemp + hnu1(j) * wh(j,i) * dijh(j,i)
          enddo
          ztemp = ztemp + hnu1(nw1+1) * wh(nw1+1,i)
          hnu2(i) = tfunc(ztemp)
        enddo

c** Calculate output nueron values
        do i=ibegoa(iwin),iendoa(iwin)
          do k=1,3
            ztemp=0.d0
            do j=1,nw2
              ztemp = ztemp + hnu2(j) * wo2(j,i,k) * dij2(j,i)
            enddo
            ztemp = ztemp + hnu2(nw2+1) * wo2(nw2+1,i,k)
            onu2(iwin,k) = tfunc(ztemp)
          enddo
        enddo
      else
c** Calculate nueron values of the first hidden layer     ---  FILTER NN
        do i=1,nfilt
          ztemp=0.d0
          do j=ibegia(iwin),iendia(iwin)
            do k=1,3
              ztemp = ztemp + finp(j,k) * fwi(j,k,i) * dijfi(j,i)
            enddo
          enddo
          ztemp = ztemp + fwi(nrx+1,1,i)
          hnuf(i) = tfunc(ztemp)
        enddo

c** Calculate output nueron values
        do i=ibegoa(iwin),iendoa(iwin)
          do k=1,3
            ztemp=0.d0
            do j=1,nfilt
              ztemp = ztemp + hnuf(j) * fwo(j,i,k) * dijfo(j,i)
            enddo
            ztemp = ztemp + hnuf(nfilt+1) * fwo(nfilt+1,i,k)
            onu2p(iwin,k,iav2) = tfunc(ztemp)
          enddo
        enddo
      endif

      end
c___________________________________________________________________

c------------- subroutine getangs -------------------------------
c--
c* Get the peak prediction
c--
      subroutine getangs()
      implicit real*8(a-h,o-z)
      PARAMETER(nrpx=8000,nrdx=2000000,nrx=21,nip=30,nw1=101,nw2=101,nrxo=1,numav=5,nfilt=21)
      common/filin/ protin(nrdx,nip+6), iprotin(nrdx,3)
      character phia*1
      common/filout/ phia(nrdx)
      common/afilt/ onu2p(nrdx,3,numav)
      common/probss/ tss(nrdx,3)
      character rtrainf(nrpx)*150
      common/sinout/ ntf,nta,rtrainf

c** Initialize averages
      do i=1,nta
        do k=1,3
          tss(i,k) = 0.d0
        enddo
      enddo

c** Calculate averages
      do j=1,numav
        do i=1,nta
          v3sm = 1.d0 / (3.d0 + onu2p(i,1,j) + onu2p(i,2,j) + onu2p(i,3,j))
          do k=1,3
            tss(i,k) = tss(i,k) + v3sm * (1.d0 + onu2p(i,k,j))
          enddo
        enddo
      enddo

      do i=1,nta
        do k=1,3
          tss(i,k) = tss(i,k) / numav
        enddo
        if (tss(i,3).ge.tss(i,2)) then
          if (tss(i,3).ge.tss(i,1)) then
            phia(i) = "H"
          else
            phia(i) = "E"
          endif
        else
          if (tss(i,1).ge.tss(i,2)) then
            phia(i) = "E"
          else 
            phia(i) = "C"
          endif
        endif
      enddo

      end
c___________________________________________________________________

c------------- function tfunc(xoutp) -------------------------------
c--
c* This is the transfer function f(x) = 1/(1+exp(alpha*x))
c* Note: the derivative f'(x) = alpha*f(x)*(1-f(x)). This fact is used
c* in the computation
c--
      real*8 function tfunc(xoutp)
      implicit real*8(a-h,o-z)
      common/misc/ alpha

      tfunc = tanh(alpha * xoutp)

      return
      end function tfunc
c___________________________________________________________________
