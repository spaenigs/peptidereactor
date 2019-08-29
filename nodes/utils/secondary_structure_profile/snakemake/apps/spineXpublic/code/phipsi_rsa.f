c     Predicting RSA using guided neural network
c	Eshel Faraggi -- Jan. 2008
c	21 window 
c	g77 -ffixed-line-length-none -Wall -W -Wsurprising  -O2 -o phipsi_rsa.e phipsi_rsa.f
c	ifort -132 -watch all -warn all -O2 -o phipsi_rsa.e phipsi_rsa.f
c	ifort -132 -O2 -o phipsi_rsa.e phipsi_rsa.f
c	ifort -132 -g -o phipsi_rsa.e phipsi_rsa.f
      program GUIDE_r
      implicit real*8(a-h,o-z)
      PARAMETER(nrpx=8000,nrdx=2000000,nrx=21,nip=30,nw1=101,nw2=101,nrxo=1,numav=5)

      common/weights/ wi(nrx+1,nip,nw1),wh(nw1+1,nw2),wo(nw2+1,nrxo)
      common/filin/ protin(nrdx,nip+3)
      common/filout/ phia(nrdx),icphi(nrdx),phiav(nrdx,numav), phiastd(nrdx)
      common/inout/ resv(nrx+1,nip)
      common/neurons/ hnu1(nw1+1),hnu2(nw2+1),onu(nrxo)
      common/edges/ ibegia(nrdx),iendia(nrdx),ibegoa(nrdx),iendoa(nrdx), nisf, nosf
      common/dcays/ dij1(nrx,nw1),dij2(nw2,nrxo),dfac,dijh(nw1,nw2)
      common/misc/ alpha, itsd
      character seqres(nrdx)
      common/seqr/ seqres
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

CCCCC  ASA:
C      "-------------------ANGLES-------------------"
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      do iav=1,numav
        call initphi(phia,icphi)
        call readwei("rsa",0,iav,spxdir)
        do ix=1,nta
            call readwin(ix,0)
            call calcnet(ix)
        enddo
        call normphiz(iav)
      enddo
      call getangs()

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCC--- OUTPUT ANGLE FILES  ---CCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      flout = trim(workdir)//"/out_rsa0"
      open(unit=99,file=flout,status='unknown')
      do i=1,nta
        tphia = (0.5d0 * phia(i) + 0.5d0) * asafact(seqres(i))
        write(99,499) i,seqres(i),tphia,phia(i),phiastd(i)
      enddo
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
499   format(I9,4x,A1,1x,F6.1,1x,F12.4,1x,F12.4)
      stop
      end
c_______________________________________________________________________
      
c------------- subroutine getdefs-------------------------------
c--
c* Get initial definitions
c--
      subroutine getdefs(listfl)
      implicit real*8(a-h,o-z)
      PARAMETER(nrpx=8000,nrdx=2000000,nrx=21,nip=30,nw1=101,nw2=101,nrxo=1,numav=5)
      common/dcays/ dij1(nrx,nw1),dij2(nw2,nrxo),dfac,dijh(nw1,nw2)
      character rtrainf(nrpx)*150
      common/sinout/ ntf,nta,rtrainf
      common/misc/ alpha, itsd

      character listfl*150

c* Definitions
c*********** alpha -- parameter for transfer function
      alpha = 0.25d0
c*********** dfac -- parameter for distance decay matrix
      dfac = 1.1d0

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
c* Initialize the neural network (the weights are randomized)
c--
      subroutine netinit()
      implicit real*8(a-h,o-z)
      PARAMETER(nrpx=8000,nrdx=2000000,nrx=21,nip=30,nw1=101,nw2=101,nrxo=1,numav=5)
      common/weights/ wi(nrx+1,nip,nw1),wh(nw1+1,nw2),wo(nw2+1,nrxo)
      common/inout/ resv(nrx+1,nip)
      common/neurons/ hnu1(nw1+1),hnu2(nw2+1),onu(nrxo)
      common/dcays/ dij1(nrx,nw1),dij2(nw2,nrxo),dfac,dijh(nw1,nw2)
      common/misc/ alpha, itsd
      
c***   setting up the biases
      resv(nrx+1,1) = 1.d0
      hnu2(nw2+1) = 1.d0
      hnu1(nw1+1) = 1.d0

        delhi = (nrx-1.d0)/(nw1-1.d0)
c***    initialize distance matrix
	do i=1,nrx
          do k=1,nw1
	    dij1(i,k) = dfac / (1.d0+dabs((k-51)*delhi-i+11))
	  enddo
	enddo

        jc = (nw1+1)/2
        jc2 = (nw2+1)/2
        mc = (nrxo+1)/2
	do i=1,nw1
          do j=1,nw2
	    dijh(i,j) = dfac / (1.d0+dabs(delhi*(i-jc-j+jc2)))
          enddo
	enddo

	do i=1,nw2
          do j=1,nrxo
	    dij2(i,j) = dfac / (1.d0+dabs((i-jc)*delhi-j+mc))
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
      PARAMETER(nrpx=8000,nrdx=2000000,nrx=21,nip=30,nw1=101,nw2=101,nrxo=1,numav=5)
      common/filin/ protin(nrdx,nip+3)
      common/edges/ ibegia(nrdx),iendia(nrdx),ibegoa(nrdx),iendoa(nrdx), nisf, nosf
      character seqres(nrdx)
      common/seqr/ seqres
      character rtrainf(nrpx)*150
      common/sinout/ ntf,nta,rtrainf
      
      character ftfl*150, inflnm*150, inflnm2*150, workdir*150, tmpchr*1, tmpchr2*1

      double precision tprot(49)

c* midpoint for input/output windows
      nisf = (nrx+1)/2
      nosf = (nrxo+1)/2

      inflnm2 = trim(workdir)//"/out_ss0"
      open(unit=44,file=inflnm2,status='old')
      nta = 0
        i = 1
        do 265 while (i.gt.0)
          nta = nta + 1
          read(44,*,end=266) ntmp, tmpchr, tmpchr2, protin(nta,nip+1), protin(nta,nip+2), protin(nta,nip+3), (tprot(j),j=1,27),
     *      ibegia(nta),iendia(nta),ibegoa(nta),iendoa(nta)
          seqres(nta) = tmpchr
          do j=1,7
            protin(nta,j) = tprot(j+20)
          enddo
          do j=1,20
            protin(nta,j+7) = tprot(j)
          enddo
265     continue
266     continue
      close(44)
      nta = nta - 1

      end
c___________________________________________________________________

c------------- subroutine initphi -------------------------------
c--
c* Initialize output angle array
c--
      subroutine initphi(phit,icphit)
      implicit real*8(a-h,o-z)
      PARAMETER(nrpx=8000,nrdx=2000000,nrx=21,nip=30,nw1=101,nw2=101,nrxo=1,numav=5)
      common/filin/ protin(nrdx,nip+3)
      character rtrainf(nrpx)*150
      common/sinout/ ntf,nta,rtrainf

      double precision phit(nrdx)
      integer icphit(nrdx)

      do i=1,nta
        phit(i) = 0.d0
        icphit(i) = 0
      enddo

      end
c___________________________________________________________________

c------------- subroutine readwei -------------------------------
c--
c* Read optimized weights
c--
      subroutine readwei(angid,ishift,iav2,spxdir)
      implicit real*8(a-h,o-z)
      PARAMETER(nrpx=8000,nrdx=2000000,nrx=21,nip=30,nw1=101,nw2=101,nrxo=1,numav=5)
      common/weights/ wi(nrx+1,nip,nw1),wh(nw1+1,nw2),wo(nw2+1,nrxo)

      character angid*3,weifl*150,weishft*2,avo,weidir*150, spxdir*150

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCC--- INPUT WEIGHTS  ---CCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      itmp = ishift + 77
      write(weishft,'(I2)') itmp
      write(avo,'(I1)') iav2
      
      weidir = trim(spxdir)//"/weights/"
      weifl = trim(weidir)//angid//"/av"//avo//"/fort."//weishft

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
      weifl = trim(weidir)//angid//"/av"//avo//"/fort."//weishft
      open(unit=77,file=weifl,status='unknown')
      do i=1,nw1+1
        do j=1,nw2
	  read(77,*) wh(i,j)
        enddo
      enddo
      close(77)
      itmp = itmp + 1
      write(weishft,'(I2)') itmp
      weifl = trim(weidir)//angid//"/av"//avo//"/fort."//weishft
      open(unit=77,file=weifl,status='unknown')
      do i=1,nw2+1
        do j=1,nrxo
	  read(77,*) wo(i,j)
        enddo
      enddo
      close(77)

      end
c___________________________________________________________________

c------------- subroutine readwin -------------------------------
c--
c* Read the current protein window
c--
      subroutine readwin(iwin,itp)
      implicit real*8(a-h,o-z)
      PARAMETER(nrpx=8000,nrdx=2000000,nrx=21,nip=30,nw1=101,nw2=101,nrxo=1,numav=5)
      common/filin/ protin(nrdx,nip+3)
      common/inout/ resv(nrx+1,nip)
      common/edges/ ibegia(nrdx),iendia(nrdx),ibegoa(nrdx),iendoa(nrdx), nisf, nosf

      integer iwin,itp

      do i=ibegia(iwin),iendia(iwin)
        irw = iwin-nisf+i
        resv(i,1) = protin(irw,nip+1)
        resv(i,2) = protin(irw,nip+2)
        resv(i,3) = protin(irw,nip+3)
        do j=4,30
          resv(i,j) = protin(irw,j-3)
        enddo
      enddo
      
      end
c___________________________________________________________________

c------------- subroutine calcnet -------------------------------
c--
c* Calculate the neural network's output angle given input residue
c* sequence parameters in resv matrix
c--
      subroutine calcnet(iwin)
      implicit real*8(a-h,o-z)
      PARAMETER(nrpx=8000,nrdx=2000000,nrx=21,nip=30,nw1=101,nw2=101,nrxo=1,numav=5)

      common/weights/ wi(nrx+1,nip,nw1),wh(nw1+1,nw2),wo(nw2+1,nrxo)
      common/filout/ phia(nrdx),icphi(nrdx),phiav(nrdx,numav), phiastd(nrdx)
      common/inout/ resv(nrx+1,nip)
      common/neurons/ hnu1(nw1+1),hnu2(nw2+1),onu(nrxo)
      common/edges/ ibegia(nrdx),iendia(nrdx),ibegoa(nrdx),iendoa(nrdx), nisf, nosf
      common/dcays/ dij1(nrx,nw1),dij2(nw2,nrxo),dfac,dijh(nw1,nw2)
      integer iwin
      double precision tfunc !transfer function

c** Calculate nueron values of the first hidden layer
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
        ztemp=0.d0
        do j=1,nw2
          ztemp = ztemp + hnu2(j) * wo(j,i) * dij2(j,i)
        enddo
        ztemp = ztemp + hnu2(nw2+1) * wo(nw2+1,i)
        onu(i) = tfunc(ztemp)
        phia(iwin+i-nosf) = phia(iwin+i-nosf) + onu(i)
        icphi(iwin+i-nosf) = icphi(iwin+i-nosf) + 1
      enddo

      end
c___________________________________________________________________

c------------- subroutine normphiz -------------------------------
c--
c* Normalize output angle (averaging)
c--
      subroutine normphiz(iav2)
      implicit real*8(a-h,o-z)
      PARAMETER(nrpx=8000,nrdx=2000000,nrx=21,nip=30,nw1=101,nw2=101,nrxo=1,numav=5)
      common/filin/ protin(nrdx,nip+3)
      common/filout/ phia(nrdx),icphi(nrdx),phiav(nrdx,numav), phiastd(nrdx)
      character rtrainf(nrpx)*150
      common/sinout/ ntf,nta,rtrainf

      do i=1,nta
        phiav(i,iav2) = phia(i) / icphi(i)
      enddo

      end
c___________________________________________________________________

c------------- subroutine getangs -------------------------------
c--
c* Get the peak prediction
c--
      subroutine getangs()
      implicit real*8(a-h,o-z)
      PARAMETER(nrpx=8000,nrdx=2000000,nrx=21,nip=30,nw1=101,nw2=101,nrxo=1,numav=5)
      common/filin/ protin(nrdx,nip+3)
      common/filout/ phia(nrdx),icphi(nrdx),phiav(nrdx,numav), phiastd(nrdx)
      character seqres(nrdx)
      common/seqr/ seqres
      character rtrainf(nrpx)*150
      common/sinout/ ntf,nta,rtrainf
      
      real*8 tmphia(nrdx)
      
      do i=1,nta
        phia(i) = phiav(i,1)
        do j=2,numav
          phia(i) = phia(i) + phiav(i,j)
        enddo
        phia(i) = phia(i) / numav
        tmphia(i) = (0.5d0 * phia(i) + 0.5d0) * asafact(seqres(i))
      enddo
      
      

      do i=1,nta
        tmp = (0.5d0 * phiav(i,1) + 0.5d0) * asafact(seqres(i))
        phiastd(i) = (tmp - tmphia(i)) * (tmp - tmphia(i))
        do j=2,numav
          tmp = (0.5d0 * phiav(i,j) + 0.5d0) * asafact(seqres(i))
          phiastd(i) = phiastd(i) + (tmp - tmphia(i)) * (tmp - tmphia(i))
        enddo
        phiastd(i) = phiastd(i) / numav
        phiastd(i) = sqrt(phiastd(i))
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
      common/misc/ alpha, itsd

      tfunc = tanh(alpha * xoutp)

      return
      end function tfunc
c___________________________________________________________________

c------------- function asafact(scr) -------------------------------
c--
c* Get ASA factor for AA
c--
      real*8 function asafact(scr)
      implicit real*8(a-h,o-z)

      character scr

      if (scr.eq."R") then
        asafact = 271.0
      elseif (scr.eq."K") then
        asafact = 257.0
      elseif (scr.eq."D") then
        asafact = 183.0
      elseif (scr.eq."E") then
        asafact = 286.0
      elseif (scr.eq."N") then
        asafact = 188.0
      elseif (scr.eq."Q") then
        asafact = 215.0
      elseif (scr.eq."H") then
        asafact = 238.0
      elseif (scr.eq."Y") then
        asafact = 250.0
      elseif (scr.eq."W") then
        asafact = 260.0
      elseif (scr.eq."S") then
        asafact = 181.0
      elseif (scr.eq."T") then
        asafact = 192.0
      elseif (scr.eq."G") then
        asafact = 136.0
      elseif (scr.eq."P") then
        asafact = 170.0
      elseif (scr.eq."A") then
        asafact = 169.0
      elseif (scr.eq."M") then
        asafact = 236.0
      elseif (scr.eq."C") then
        asafact = 139.0
      elseif (scr.eq."F") then
        asafact = 221.0
      elseif (scr.eq."L") then
        asafact = 221.0
      elseif (scr.eq."V") then
        asafact = 171.0
      elseif (scr.eq."I") then
        asafact = 210.0
      endif 

      return
      end function asafact
c___________________________________________________________________

