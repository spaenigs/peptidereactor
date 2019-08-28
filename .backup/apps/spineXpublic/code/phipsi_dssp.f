c     Retraining Protein Structure Prediction by Artificial Neural Network
c	21 window 
c	g77 -ffixed-line-length-none -Wall -W -Wsurprising  -O2 -o spiner.e spiner.f
c	ifort -132 -watch all -warn all -O2 -o spiner.e spiner.f
c	ifort -132 -O2 -o spiner.e spiner.f
c	this code originated from /data3/efaraggi/protein/spine/crossvald/new/t6.peak/t2.average/corrected
c	later from: /data3/efaraggi/protein/spine2.0/t7.sprt/predmach/step1/combine/tatctdteti/h150
c	08/31/07:   /data3/efaraggi/papers/2peak/data/weak/j2pk/weaker3/boosting/t5
c	11/07/07:   retraining on missed peak clasifications.
c	from:	    /data3/efaraggi/papers/2peak/data/weak/j2pk/weaker3/boosting/t6.new/linear/nodisca/h14/corrected
c	12/24/07    /data3/efaraggi/papers/2peak/data/CV10/t1.step1/cv01/psi/av01
c

      program SPINE_r
      implicit real*8(a-h,o-z)
      PARAMETER(nrdx=1300,nrpx=2640,naax=600000,nrx=21,nip=30,nw1=101,nw2=101,nrxo=1,nfilt=21)

      common/weights/ wi(nrx+1,nip,nw1),wh(nw1+1,nw2)
      common/pwei/ wo2(nw2+1,nrxo,3)
      common/filtwei/ fwi(nrx+1,3,nfilt), fwo(nfilt+1,nrxo,3)
      common/filin/ protin(naax,nip+2),phir2(naax,3),phip(naax),phir(naax),dphir(naax)
      common/filout/ phia(naax),icphi(naax)
      common/inout/ resv(nrx+1,nip),phi2(nrxo,3)
      common/neurons/ hnu1(nw1+1),hnu2(nw2+1),onu2(naax,3),hnuf(nfilt+1),onuf(nrxo,3), finp(nrx,3)
      common/learn/ alpha, ratel, rmoment
      common/edges/ ibegia(naax),iendia(naax),ibegoa(naax),iendoa(naax), nisf, nosf
      common/misc/ rmae, ccoef, q10scr, q10prh, cvfrac, itsd
      common/dcays/ dij1(nrx,nw1),dij2(nw2,nrxo),dfac,dijh(nw1,nw2),dijfi(nrx,nfilt),dijfo(nfilt,nrxo)
      common/sinout/ nta
      common/afilt/ onu2p(naax,3)

      character fylevld*120, flprof*150, flspx*150, cassi(naax)*1, aanm(naax), tmpfl*150
      double precision kponu2p(naax,3)
      
c* Get input/native files
      if (iargc().ge.2) then
        nrandm = 7
        call getarg(1,flprof)
        call getarg(2,flspx)
      else
        write(*,*) "use: phipsi_dssp.e profile_file spineXpred_file"
        write(*,*) "To get SS prediction according to DSSP classification"
        stop
      endif
      
c* Definition of the general input/output files
      call getdefs()
      call netinit(nrandm)    ! initialize weights
      call setio(flprof,flspx,aanm)	! initialize input/output
      do i=1,nta
        do j=1,3
          kponu2p(i,j) = 0.d0
        enddo
      enddo
      do iav=1,5
        call readwei(iav)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCC--- OUTPUT CROSS VALIDATION  ---CCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        call initphi(phia,icphi)
        do ix=1,nta
          call readwin(ix,0,0)
          call calcnet(ix,0)
        enddo
c        call normphi(phia)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      call readweif(iav)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCC--- OUTPUT CROSS VALIDATION  ---CCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        call initphi(phia,icphi)
        do ix=1,nta
          call readwinf(ix)
          call calcnetf(ix,0)
        enddo
c        call normphi(phia)
        do i=1,nta
          do j=1,3
            kponu2p(i,j) = kponu2p(i,j) + onu2p(i,j)
          enddo
        enddo
      enddo
      do i=1,nta
        do j=1,3
          onu2p(i,j) = kponu2p(i,j) / 5.d0
        enddo

        if (onu2p(i,3).ge.onu2p(i,2)) then
          if (onu2p(i,3).ge.onu2p(i,1)) then
            cassi(i) = "H"
          else
            cassi(i) = "E"
          endif
        else
          if (onu2p(i,1).ge.onu2p(i,2)) then
            cassi(i) = "E"
          else 
            cassi(i) = "C"
          endif
        endif
      enddo
      
      tmpfl = trim(flspx)//".dssp"
      open(unit=99,file=tmpfl,status='unknown')
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      do i=1,nta
        v3sm = 1.d0 / (3.d0 + onu2p(i,1) + onu2p(i,2) + onu2p(i,3))
        v1 = v3sm * (1.d0 + onu2p(i,1))
        v2 = v3sm * (1.d0 + onu2p(i,2))
        v3 = v3sm * (1.d0 + onu2p(i,3))
        write(99,499) i,aanm(i),cassi(i),v1,v2,v3
      enddo
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
399   format (I6,2x,D12.4,2x,D12.4,2x,D12.4)
499   format(I8,2x,A2,2x,A2,2x,D12.4,2x,D12.4,2x,D12.4)
      
      close(99)
      write(*,*) "DSSP prediction in ", trim(tmpfl)
      
      stop
      end
c_______________________________________________________________________
      
c------------- subroutine getdefs-------------------------------
c--
c* Get initial definitions
c--
      subroutine getdefs()
      implicit real*8(a-h,o-z)
      PARAMETER(nrdx=1300,nrpx=2640,naax=600000,nrx=21,nip=30,nw1=101,nw2=101,nrxo=1,nfilt=21)
      common/learn/ alpha, ratel, rmoment
      common/dcays/ dij1(nrx,nw1),dij2(nw2,nrxo),dfac,dijh(nw1,nw2),dijfi(nrx,nfilt),dijfo(nfilt,nrxo)
      common/misc/ rmae, ccoef, q10scr, q10prh, cvfrac, itsd
      common/sinout/ nta

c*********** alpha -- parameter for transfer function
      alpha = 0.2d0
c*********** dfac -- parameter for distance decay matrix
      dfac = 2.d0
c*********** icv -- set for 10 fold cross validation

      end
c___________________________________________________________________

c------------- subroutine netinit -------------------------------
c--
c* Initialize the neural network (the weights are randomized)
c--
      subroutine netinit(krndm)
      implicit real*8(a-h,o-z)
      PARAMETER(nrdx=1300,nrpx=2640,naax=600000,nrx=21,nip=30,nw1=101,nw2=101,nrxo=1,nfilt=21)

      common/weights/ wi(nrx+1,nip,nw1),wh(nw1+1,nw2)
      common/pwei/ wo2(nw2+1,nrxo,3)
      common/filtwei/ fwi(nrx+1,3,nfilt), fwo(nfilt+1,nrxo,3)
      common/inout/ resv(nrx+1,nip),phi2(nrxo,3)
      common/neurons/ hnu1(nw1+1),hnu2(nw2+1),onu2(naax,3),hnuf(nfilt+1),onuf(nrxo,3), finp(nrx,3)
      common/misc/ rmae, ccoef, q10scr, q10prh, cvfrac, itsd
      common/dcays/ dij1(nrx,nw1),dij2(nw2,nrxo),dfac,dijh(nw1,nw2),dijfi(nrx,nfilt),dijfo(nfilt,nrxo)
      
      integer*4 timeArray(3)    ! Holds the hour, minute, and second
      double precision randnm !random LCG

        call itime(timeArray)     ! Get the current time for seed
        ta1 = timeArray(1) / 1.74
        ta2 = timeArray(2) / 4.35
        ta3 = timeArray(3) / 2.8
        jta3 = timeArray(3)
        ita3 = exp(ta3)
        ta3 = (mod(ita3,61)/2.8)
        ta4 = exp(ta1) + exp(ta2) + exp(ta3) + krndm
        ita4log = log10(ta4)
        ita4 = ta4
        dta4 = 10**ita4log * (ta4 - ita4)
        itsd = abs(dta4 - ta4) + ta3
c      itsd = 44027
        write(*,*) "random seed = ",itsd

c***   setting up the biases
      resv(nrx+1,1) = 1.d0
      hnu2(nw2+1) = 1.d0
      hnu1(nw1+1) = 1.d0
      hnuf(nfilt+1) = 1.d0

c***    randomizing weights
	do i=1,nrx+1
	  do j=1,nip
	    do k=1,nw1
	      atmp = randnm() - 0.5
	      wi(i,j,k) = atmp
	    enddo
	  enddo
	enddo

	do i=1,nw1+1
          do j=1,nw2
	    atmp = randnm() - 0.5
	    wh(i,j) = atmp
          enddo
	enddo

	do i=1,nw2+1
          do j=1,nrxo
            do k=1,3
	      atmp = randnm() - 0.5
   	      wo2(i,j,k) = atmp
            enddo
          enddo
	enddo

	do i=1,nrx+1
	  do j=1,3
	    do k=1,nfilt
	      atmp = randnm() - 0.5
	      fwi(i,j,k) = atmp
	    enddo
	  enddo
	enddo
	do i=1,nfilt+1
          do j=1,3
	    atmp = randnm() - 0.5
	    fwo(i,1,j) = atmp
          enddo
	enddo

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

c------------- subroutine readwei -------------------------------
c--
c* Read optimized weights
c--
      subroutine readwei(iav)
      implicit real*8(a-h,o-z)
      PARAMETER(nrdx=1300,nrpx=2640,naax=600000,nrx=21,nip=30,nw1=101,nw2=101,nrxo=1,nfilt=21)
      common/weights/ wi(nrx+1,nip,nw1),wh(nw1+1,nw2)
      common/pwei/ wo2(nw2+1,nrxo,3)
      
      character avo*1,weidir*150,predfl*150

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCC--- INPUT WEIGHTS  ---CCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      write(avo,'(I1)') iav
      weidir = "/data3/efaraggi/server/spX/spineX/weights/dssp/av"
      predfl = trim(weidir)//avo//"/fort.77"
      open(unit=76,file=predfl,status='unknown')
      do i=1,nrx+1
	do j=1,nip
	  do k=1,nw1
	    read(76,*) wi(i,j,k)
	  enddo
	enddo
      enddo
      close(76)
      predfl = trim(weidir)//avo//"/fort.78"
      open(unit=76,file=predfl,status='unknown')
      do i=1,nw1+1
        do j=1,nw2
	  read(76,*) wh(i,j)
        enddo
      enddo
      close(76)
      predfl = trim(weidir)//avo//"/fort.79"
      open(unit=76,file=predfl,status='unknown')
      do i=1,nw2+1
        do j=1,nrxo
	  read(76,*) (wo2(i,j,k),k=1,3)
        enddo
      enddo
      close(76)

      end
c___________________________________________________________________

c------------- subroutine readweif -------------------------------
c--
c* Read optimized weights
c--
      subroutine readweif(iav)
      implicit real*8(a-h,o-z)
      PARAMETER(nrdx=1300,nrpx=2640,naax=600000,nrx=21,nip=30,nw1=101,nw2=101,nrxo=1,nfilt=21)
      common/filtwei/ fwi(nrx+1,3,nfilt), fwo(nfilt+1,nrxo,3)

      character avo*1,weidir*150,predfl*150

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCC--- INPUT WEIGHTS  ---CCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      write(avo,'(I1)') iav
      weidir = "/data3/efaraggi/server/spX/spineX/weights/dssp/av"
      predfl = trim(weidir)//avo//"/fort.87"
      open(unit=76,file=predfl,status='unknown')
      do i=1,nrx+1
	do j=1,3
	  do k=1,nfilt
	    read(76,*) fwi(i,j,k)
	  enddo
	enddo
      enddo
      close(76)
      predfl = trim(weidir)//avo//"/fort.88"
      open(unit=76,file=predfl,status='unknown')
      do i=1,nfilt+1
        do j=1,3
	  read(76,*) fwo(i,1,j)
        enddo
      enddo
      close(76)

      end
c___________________________________________________________________

c------------- subroutine setio -------------------------------
c--
c* Setup the input/output arrays
c--
      subroutine setio(flprof, flspx, aanm)
      implicit real*8(a-h,o-z)
      PARAMETER(nrdx=1300,nrpx=2640,naax=600000,nrx=21,nip=30,nw1=101,nw2=101,nrxo=1,nfilt=21)
      common/filin/ protin(naax,nip+2),phir2(naax,3),phip(naax),phir(naax),dphir(naax)
      common/edges/ ibegia(naax),iendia(naax),ibegoa(naax),iendoa(naax), nisf, nosf
      common/sinout/ nta

      character tempc, aa2, aa3, predfl*120, flprof*150, flspx*150, aanm(naax)
      double precision tprot(nrdx,60)

      character aapp(22)
      double precision physp(7,22), sevec(90,7)
      data aapp /"R","K","D","E","N","Q","H","Y","W","S","T","G","P","A","M","C","F","L","V","I","B","Z"/
      data physp /0.105, 0.373, 0.466,-0.900, 0.900, 0.528,-0.371,
     $ -0.088, 0.066, 0.163,-0.889, 0.727, 0.279,-0.265,
     $ -0.213,-0.417,-0.281,-0.767,-0.900,-0.155,-0.635,
     $ -0.230,-0.241,-0.058,-0.696,-0.868, 0.900,-0.582,
     $ -0.213,-0.329,-0.243,-0.674,-0.075,-0.403,-0.529,
     $ -0.230,-0.110,-0.020,-0.464,-0.276, 0.528,-0.371,
     $ 0.384, 0.110, 0.138,-0.271, 0.195,-0.031,-0.106,
     $ 0.363, 0.417, 0.541, 0.188,-0.274,-0.155, 0.476,
     $ 0.479, 0.900, 0.900, 0.900,-0.209, 0.279, 0.529,
     $ -0.337,-0.637,-0.544,-0.364,-0.265,-0.466,-0.212,
     $ 0.402,-0.417,-0.321,-0.199,-0.288,-0.403, 0.212,
     $ -0.900,-0.900,-0.900,-0.342,-0.179,-0.900,-0.900,
     $ 0.247,-0.900,-0.294, 0.055,-0.010,-0.900, 0.106,
     $ -0.350,-0.680,-0.677,-0.171,-0.170, 0.900,-0.476,
     $ 0.110, 0.066, 0.087, 0.337,-0.262, 0.652,-0.001,
     $ -0.140,-0.329,-0.359, 0.508,-0.114,-0.652, 0.476,
     $ 0.363, 0.373, 0.412, 0.646,-0.272, 0.155, 0.318,
     $ 0.213,-0.066,-0.009, 0.596,-0.186, 0.714,-0.053,
     $ 0.677,-0.285,-0.232, 0.331,-0.191,-0.031, 0.900,
     $ 0.900,-0.066,-0.009, 0.652,-0.186, 0.155, 0.688,
     $ -0.213,-0.417,-0.281,-0.767,-0.900,-0.155,-0.635, 
     $ -0.230,-0.241,-0.058,-0.696,-0.868, 0.900,-0.582/

c* midpoint for input/output windows
      nisf = (nrx+1)/2
      nosf = (nrxo+1)/2

c* set up aa parameter matrix:
      do i=1,22
        iaa = ichar(aapp(i))
        do j=1,7
          sevec(iaa,j) = physp(j,i)
        enddo
      enddo

      nta = 0
      ntf = 1
      do ix=1,ntf
        ntao = nta
        predfl = trim(flprof)
        open(unit=45,file=predfl,status='old')
        read(45,*) ; read(45,*) ; read(45,*)
        predfl = trim(flspx)
        open(unit=46,file=predfl,status='old')
253     continue
        read(46,'(A1)') tempc
        if (tempc.ne." ") goto 253
        backspace(46)
        i = 1
        do 265 while (i.gt.0)
          read(45,*,end=266,err=266) tempc,aa2,(tprot(1,j),j=1,20)
          read(46,*,end=266,err=266) itemp,aa3,tempc,(tprot(2,j),j=1,8)
          nta = nta + 1

          if (aa2.ne.aa3) then
            write(*,*) "aa23 mismatch, aborting"
            write(*,*) "i=",i,"    nta=", nta, "    ix=", ix
            write(*,*) "aa2: ",aa2, "     aa3: ", aa3
            write(*,*) "file: ", trim(flspx)
            stop
          endif
          
          iaa = ichar(aa2)
          aanm(nta) = aa2
          do j=1,7
            protin(nta,j) = sevec(iaa,j)/0.9d0
          enddo

          do j=1,20
            protin(nta,7+j) = tprot(1,j)/9.d0
          enddo

          protin(nta,28) = 2.d0 * (tprot(2,8)/asafact(aa2) - 0.5d0)

          protin(nta,29) = tprot(2,1) / 180.d0
          tphi = (180.d0 + tprot(2,2)) / 360.d0
          if (tphi.gt.0.1389d0) then
            tphi = (tphi-0.1389d0)
          else 
            tphi = (tphi+0.8611d0)
          endif
          tphi = 2.d0*(tphi - 0.5d0)
          protin(nta,30) = tphi

          tphi = 0.5d0
          phip(nta) = tphi

          if (dabs(tphi+0.6).le.0.1d0) then
            phir2(nta,1) = 1.d0
            phir2(nta,2) = -1.d0
            phir2(nta,3) = -1.d0
          elseif (dabs(tphi+0.2).le.0.1d0) then
            phir2(nta,1) = -1.d0
            phir2(nta,2) = 1.d0
            phir2(nta,3) = -1.d0
          elseif (dabs(tphi-0.5).le.0.1d0) then
            phir2(nta,1) = -1.d0
            phir2(nta,2) = -1.d0
            phir2(nta,3) = 1.d0
          endif

          i = i + 1
265     continue
266     continue
        close(4)
        numres = i-1  ! number of residues in current protein

c* Setup the edges for each amino acid
        do i=1,numres
          inta = ntao + i
          ib = 1 + nisf - i
          if(ib.lt.1) ib = 1
          ibegia(inta) = ib
          ib = numres + nisf - i
          if(ib.gt.nrx) ib = nrx
          iendia(inta) = ib
          ib = 1 + nosf - i
          if(ib.lt.1) ib = 1
          ibegoa(inta) = ib
          ib = numres + nosf - i
          if(ib.gt.nrxo) ib = nrxo
          iendoa(inta) = ib
        enddo
      enddo

      end
c___________________________________________________________________

c------------- subroutine initphi -------------------------------
c--
c* Initialize output angle array
c--
      subroutine initphi(phit,icphit)
      implicit real*8(a-h,o-z)
      PARAMETER(nrdx=1300,nrpx=2640,naax=600000,nrx=21,nip=30,nw1=101,nw2=101,nrxo=1,nfilt=21)
      common/sinout/ nta

      double precision phit(naax)
      integer icphit(naax)

      do i=1,nta
        phit(i) = 0.d0
        icphit(i) = 0
      enddo

      end
c___________________________________________________________________

c------------- subroutine readwin -------------------------------
c--
c* Read the current protein window
c--
      subroutine readwin(iwin,itp,itp2)
      implicit real*8(a-h,o-z)
      PARAMETER(nrdx=1300,nrpx=2640,naax=600000,nrx=21,nip=30,nw1=101,nw2=101,nrxo=1,nfilt=21)
      common/filin/ protin(naax,nip+2),phir2(naax,3),phip(naax),phir(naax),dphir(naax)
      common/inout/ resv(nrx+1,nip),phi2(nrxo,3)
      common/edges/ ibegia(naax),iendia(naax),ibegoa(naax),iendoa(naax), nisf, nosf
      
      integer iwin,itp,itp2

      do i=ibegia(iwin),iendia(iwin)
        irw = iwin-nisf+i
        do j=1,7
          resv(i,j) = protin(irw,j)	!physical parameters
        enddo
        do j=1,20
          resv(i,j+7) = protin(irw,j+7)	!sequence profile
        enddo
        resv(i,28) = protin(irw,28)	!ASA prediction
        resv(i,29) = protin(irw,29)	!phi prediction
        resv(i,30) = protin(irw,30)	!psi prediction
      enddo
      
      do i=ibegoa(iwin),iendoa(iwin)
        irw = iwin-nosf+i
        if (itp2.le.0) then		! PEAK (itp2=0)
          phi2(i,1) = phir2(irw,1)
          phi2(i,2) = phir2(irw,2)
          phi2(i,3) = phir2(irw,3)
        else				! dPHI (itp2=1)
          phi2(i,1) = dphir(irw)
          phi2(i,2) = dphir(irw)
        endif
      enddo

      end
c___________________________________________________________________

c------------- subroutine calcnet -------------------------------
c--
c* Calculate the neural network's output angle given input residue
c* sequence parameters in resv matrix
c--
      subroutine calcnet(iwin,itp)
      implicit real*8(a-h,o-z)
      PARAMETER(nrdx=1300,nrpx=2640,naax=600000,nrx=21,nip=30,nw1=101,nw2=101,nrxo=1,nfilt=21)
      common/weights/ wi(nrx+1,nip,nw1),wh(nw1+1,nw2)
      common/pwei/ wo2(nw2+1,nrxo,3)
      common/filout/ phia(naax),icphi(naax)
      common/inout/ resv(nrx+1,nip),phi2(nrxo,3)
      common/neurons/ hnu1(nw1+1),hnu2(nw2+1),onu2(naax,3),hnuf(nfilt+1),onuf(nrxo,3), finp(nrx,3)
      common/edges/ ibegia(naax),iendia(naax),ibegoa(naax),iendoa(naax), nisf, nosf
      common/dcays/ dij1(nrx,nw1),dij2(nw2,nrxo),dfac,dijh(nw1,nw2),dijfi(nrx,nfilt),dijfo(nfilt,nrxo)
      integer iwin,itp
      
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
        do k=1,3
          ztemp=0.d0
          do j=1,nw2
            ztemp = ztemp + hnu2(j) * wo2(j,i,k) * dij2(j,i)
          enddo
          ztemp = ztemp + hnu2(nw2+1) * wo2(nw2+1,i,k)
          onu2(iwin,k) = tfunc(ztemp)
        enddo

CC!CC: Need to modify next lines and normphi to run more than one ouput (i.e., if nrxo>1)
CCCCC! PEAK (itp=0)
        if (itp.le.0) then
          if (onu2(iwin,3).ge.onu2(iwin,2)) then
            if (onu2(iwin,3).ge.onu2(iwin,1)) then
              r2phin = 0.5d0
            else
              r2phin = -0.6d0
            endif
          else
            if (onu2(iwin,1).ge.onu2(iwin,2)) then
              r2phin = -0.6d0
            else 
              r2phin = -0.2d0
            endif
          endif
          phia(iwin+i-nosf) = phia(iwin+i-nosf) + r2phin
          icphi(iwin+i-nosf) = icphi(iwin+i-nosf) + 1
CCCCC! dPHI (itp=1)
        endif
      enddo

      end
c___________________________________________________________________

c------------- subroutine normphi -------------------------------
c--
c* Normalize peak assignement (averaging)
c--
      subroutine normphi(phit)
      implicit real*8(a-h,o-z)
      PARAMETER(nrdx=1300,nrpx=2640,naax=600000,nrx=21,nip=30,nw1=101,nw2=101,nrxo=1,nfilt=21)
      common/sinout/ nta

      double precision phit(naax)

      do i=1,nta
        if (phit(i).lt.0) then
          phit(i) = -1.d0
        else
          phit(i) = 1.d0
        endif
      enddo

      end
c___________________________________________________________________

c------------- subroutine readwinf -------------------------------
c--
c* Read the current protein window
c--
      subroutine readwinf(iwin)
      implicit real*8(a-h,o-z)
      PARAMETER(nrdx=1300,nrpx=2640,naax=600000,nrx=21,nip=30,nw1=101,nw2=101,nrxo=1,nfilt=21)
      common/filin/ protin(naax,nip+2),phir2(naax,3),phip(naax),phir(naax),dphir(naax)
      common/inout/ resv(nrx+1,nip),phi2(nrxo,3)
      common/neurons/ hnu1(nw1+1),hnu2(nw2+1),onu2(naax,3),hnuf(nfilt+1),onuf(nrxo,3), finp(nrx,3)
      common/edges/ ibegia(naax),iendia(naax),ibegoa(naax),iendoa(naax), nisf, nosf
      
      integer iwin

      do i=ibegia(iwin),iendia(iwin)
        irw = iwin-nisf+i
        do j=1,3
          finp(i,j) = onu2(irw,j)
        enddo
      enddo
      
      do i=ibegoa(iwin),iendoa(iwin)
        irw = iwin-nosf+i
        do j=1,3
          phi2(i,j) = phir2(irw,j)
        enddo
      enddo

      end
c___________________________________________________________________

c------------- subroutine calcnetf -------------------------------
c--
c* Calculate the neural network's output 
c--
      subroutine calcnetf(iwin,itp)
      implicit real*8(a-h,o-z)
      PARAMETER(nrdx=1300,nrpx=2640,naax=600000,nrx=21,nip=30,nw1=101,nw2=101,nrxo=1,nfilt=21)
      common/filtwei/ fwi(nrx+1,3,nfilt), fwo(nfilt+1,nrxo,3)
      common/filout/ phia(naax),icphi(naax)
      common/neurons/ hnu1(nw1+1),hnu2(nw2+1),onu2(naax,3),hnuf(nfilt+1),onuf(nrxo,3), finp(nrx,3)
      common/edges/ ibegia(naax),iendia(naax),ibegoa(naax),iendoa(naax), nisf, nosf
      common/dcays/ dij1(nrx,nw1),dij2(nw2,nrxo),dfac,dijh(nw1,nw2),dijfi(nrx,nfilt),dijfo(nfilt,nrxo)
      integer iwin,itp
      common/afilt/ onu2p(naax,3)
      
      double precision tfunc !transfer function

c** Calculate nueron values of the first hidden layer
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
          onuf(i,k) = tfunc(ztemp)
          onu2p(iwin,k) = onuf(i,k)
        enddo

CC!CC: Need to modify next lines and normphi to run more than one ouput (i.e., if nrxo>1)
CCCCC! PEAK (itp=0)
        if (itp.le.0) then
          if (onuf(i,3).ge.onuf(i,2)) then
            if (onuf(i,3).ge.onuf(i,1)) then
              r2phin = 0.5d0
            else
              r2phin = -0.6d0
            endif
          else
            if (onuf(i,1).ge.onuf(i,2)) then
              r2phin = -0.6d0
            else 
              r2phin = -0.2d0
            endif
          endif
          phia(iwin+i-nosf) = phia(iwin+i-nosf) + r2phin
          icphi(iwin+i-nosf) = icphi(iwin+i-nosf) + 1
CCCCC! dPHI (itp=1)
        endif
      enddo

      end
c___________________________________________________________________

c------------- subroutine normphiz -------------------------------
c--
c* Normalize output angle (averaging)
c--
      subroutine normphiz()
      implicit real*8(a-h,o-z)
      PARAMETER(nrdx=1300,nrpx=2640,naax=600000,nrx=21,nip=30,nw1=101,nw2=101,nrxo=1,nfilt=21)
      common/filout/ phia(naax),icphi(naax)
      common/sinout/ nta

      do i=1,nta
        phia(i) = phia(i) / icphi(i)
      enddo

      end
c___________________________________________________________________

c------------- subroutine calcerr -------------------------------
c--
c* Calculate the neural network's error
c--
      subroutine calcerr(phit,istart,ifinish)
      implicit real*8(a-h,o-z)
      PARAMETER(nrdx=1300,nrpx=2640,naax=600000,nrx=21,nip=30,nw1=101,nw2=101,nrxo=1,nfilt=21)
      common/misc/ rmae, ccoef, q10scr, q10prh, cvfrac, itsd
      common/filin/ protin(naax,nip+2),phir2(naax,3),phip(naax),phir(naax),dphir(naax)
      
      double precision phit(naax)

c** Calculate error
      q10prh = 0.d0
      do i=istart,ifinish
        if (dabs(phit(i)-phip(i)).le.0.1d0) then
          q10prh = q10prh + 1.d0
        endif
      enddo
      q10prh = q10prh / (ifinish - istart + 1)

      end
c___________________________________________________________________

c------------- subroutine calcor -------------------------------
c--
c* Calculate the neural network's correlation and Q10 errors
c--
      subroutine calcor(phit,istart,ifinish)
      implicit real*8(a-h,o-z)
      PARAMETER(nrdx=1300,nrpx=2640,naax=600000,nrx=21,nip=30,nw1=101,nw2=101,nrxo=1,nfilt=21)
      common/misc/ rmae, ccoef, q10scr, q10prh, cvfrac, itsd
      common/filin/ protin(naax,nip+2),phir2(naax,3),phip(naax),phir(naax),dphir(naax)
      
      double precision phit(naax)
      double precision tphi(naax),tphip(naax)

c** Calculate error, correlation and Q10
      rmae = 0.d0
      avphi = 0.d0
      avphip = 0.d0
      q10scr = 0.d0
      q10prh = 0.d0
      do i=istart,ifinish
        erro1 = dabs(phir(i) - phit(i))
        rmae = rmae + erro1
        avphi = avphi + phir(i)
	avphip = avphip + phit(i)
        iqrh = (1.d0+phir(i)) / 0.2d0
        iqph = (1.d0+phit(i)) / 0.2d0
        if (iqrh.eq.iqph) q10scr = q10scr + 1.d0
        if (erro1.le.0.2d0) q10prh = q10prh + 1
      enddo
      rmae = rmae / (ifinish - istart + 1)
      avphi = avphi / (ifinish - istart + 1)
      avphip = avphip / (ifinish - istart + 1)
      q10scr = q10scr / (ifinish - istart + 1)
      q10prh = q10prh / (ifinish - istart + 1)

      do i=istart,ifinish
        tphi(i) = phir(i) - avphi
        tphip(i) = phit(i) - avphip
      enddo

      crhn = 0.d0
      cnor1hn = 0.d0
      cnor2hn = 0.d0
      do i=istart,ifinish
	crhn = crhn + tphi(i)*tphip(i)
	cnor1hn = cnor1hn + tphi(i)*tphi(i)
	cnor2hn = cnor2hn + tphip(i)*tphip(i)
      enddo
      ccoef = crhn / (sqrt(cnor1hn)*sqrt(cnor2hn))

      end
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
      elseif (scr.eq."X") then
        asafact = 0.0
      else
        write(*,*) "can't find residue: ", scr
        stop
      endif 

      return
      end function asafact
c___________________________________________________________________

c------------- function tfunc(xoutp) -------------------------------
c--
c* This is the transfer function f(x) = 1/(1+exp(alpha*x))
c* Note: the derivative f'(x) = alpha*f(x)*(1-f(x)). This fact is used
c* in the computation
c--
      real*8 function tfunc(xoutp)
      implicit real*8(a-h,o-z)
      common/learn/ alpha, ratel, rmoment

      tfunc = tanh(alpha * xoutp)

      return
      end function tfunc
c___________________________________________________________________

c------------- function randnm() -------------------------------
c--
c* LCG psuedorandom 
c--
      real*8 function randnm()
      implicit real*8(a-h,o-z)
      common/misc/ rmae, ccoef, q10scr, q10prh, cvfrac, itsd
      
      dvs = 2147483647.d0
      dnm = 2147483711.d0
      dml = 16807.d0

      if(itsd .lt. 1) itsd = itsd + dvs
      ditsd = itsd * dml
      ditsd = mod(ditsd, dvs)
      randnm = itsd / dnm
      itsd = ditsd

      return
      end function randnm
c___________________________________________________________________
