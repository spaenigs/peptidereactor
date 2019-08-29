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
      common/oldwei/ oldwi(nrx+1,nip,nw1),oldwh(nw1+1,nw2)
      common/pwei/ wo2(nw2+1,nrxo,3),oldwo2(nw2+1,nrxo,3)
      common/filtwei/ fwi(nrx+1,3,nfilt), fwo(nfilt+1,nrxo,3)
      common/ofiltwei/ ofwi(nrx+1,3,nfilt), ofwo(nfilt+1,nrxo,3)
      common/filin/ protin(naax,nip+2),phir2(naax,3),phip(naax),phir(naax),dphir(naax)
      common/filout/ phia(naax),icphi(naax)
      common/inout/ resv(nrx+1,nip),phi2(nrxo,3)
      common/neurons/ hnu1(nw1+1),hnu2(nw2+1),onu2(naax,3),hnuf(nfilt+1),onuf(nrxo,3), finp(nrx,3)
      common/learn/ alpha, ratel, rmoment
      common/edges/ ibegia(naax),iendia(naax),ibegoa(naax),iendoa(naax), nisf, nosf
      common/misc/ rmae, ccoef, q10scr, q10prh, cvfrac, itsd
      common/dcays/ dij1(nrx,nw1),dij2(nw2,nrxo),dfac,dijh(nw1,nw2),dijfi(nrx,nfilt),dijfo(nfilt,nrxo)
      character rtrainf(nrpx)*120
      common/sinout/ numtrf,ntotf,ntf,numtra,ntota,nta,idxprt(nrpx),rtrainf
      common/peaks/ peks(naax)
      common/afilt/ onu2p(naax,3)

      character fylevld*120, crandm*120
      
c* Get random number initiator
      if (iargc().ge.1) then
        call getarg(1,crandm)
        read(crandm,*) atemp
	nrandm=nint(atemp)
      else
        nrandm = 7
      endif
      
c* Definition of the general input/output files
      call getdefs(fylevld,nte,npe,icv)
      open(unit=15,file=fylevld)
      call netinit(nrandm)    ! initialize weights
      call setio()	! initialize input/output
      cormax = -20.d0
      icntinc = 0
      write(*,*) "-------------------seconds-------------------"
      write(15,*) "-------------------seconds-------------------"
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCC---  TRAINING SECONDS  ---CCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C**   EPOCH:
      do ie=1,nte
        call initphi(phia,icphi)
        do ix=1,numtra
          call readwin(ix,0,0)
          call calcnet(ix,0)
          call backprop(ix)
        enddo
        do ix=numtra+1,ntota
          call readwin(ix,1,0)
          call calcnet(ix,0)
        enddo
        do ix=ntota+1,nta
          call readwin(ix,0,0)
          call calcnet(ix,0)
          call backprop(ix)
        enddo
c        call normphi(phia)
        call calcerr(phia,1,numtra)
        q10prho1 = q10prh
        call calcerr(phia,numtra+1,ntota)
        q10prho2 = q10prh
        call calcerr(phia,ntota+1,nta)
        write(15,399) ie, q10prho1, q10prho2, q10prh
C* END OF EXEMPLAR.
C* VALIDATION:
        if (q10prho2.gt.cormax) then
          call outwei(0)
          cormax = q10prho2
          icntinc = 0
        else
          icntinc = icntinc + 1
        endif
        if (icntinc.gt.100) then
          write(*,*) ie, npe, cormax
          goto 407
        endif
      enddo
C* END OF EPOCHS.
407   continue
      call readwei(0)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCC--- OUTPUT CROSS VALIDATION  ---CCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        call initphi(phia,icphi)
        do ix=1,numtra
          call readwin(ix,0,0)
          call calcnet(ix,0)
        enddo
        do ix=numtra+1,nta
          call readwin(ix,1,0)
          call calcnet(ix,0)
        enddo
c        call normphi(phia)
        call calcerr(phia,ntota+1,nta)
        write(*,*) "Q10P:", q10prh
	write(15,*) 
	write(15,*) "#Cross Validation"
	write(15,*) "#----------------"
	write(15,*) ie,q10prh
	write(15,*) 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      write(15,*) "#After filter"
      write(15,*) "#------------"
      write(*,*) "#After filter"
      write(*,*) "#------------"
      cormax = -20.d0
      icntinc = 0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCC---  TRAINING FILTER  ---CCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C**   EPOCH:
      do ie=1,nte
        call initphi(phia,icphi)
        do ix=1,numtra
          call readwinf(ix)
          call calcnetf(ix,0)
          call backpropf(ix)
        enddo
        do ix=numtra+1,ntota
          call readwinf(ix)
          call calcnetf(ix,0)
        enddo
        do ix=ntota+1,nta
          call readwinf(ix)
          call calcnetf(ix,0)
          call backpropf(ix)
        enddo
c        call normphi(phia)
        call calcerr(phia,1,numtra)
        q10prho1 = q10prh
        call calcerr(phia,numtra+1,ntota)
        q10prho2 = q10prh
        call calcerr(phia,ntota+1,nta)
        write(15,399) ie, q10prho1, q10prho2, q10prh
C* END OF EXEMPLAR.
C* VALIDATION:
        if (q10prho2.gt.cormax) then
          call outweif(0)
          cormax = q10prho2
          icntinc = 0
        else
          icntinc = icntinc + 1
        endif
        if (icntinc.gt.100) then
          write(*,*) ie, npe, cormax
          goto 411
        endif
      enddo
C* END OF EPOCHS.
411   continue
      call readweif(0)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCC--- OUTPUT CROSS VALIDATION  ---CCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        call initphi(phia,icphi)
        do ix=1,numtra
          call readwinf(ix)
          call calcnetf(ix,0)
        enddo
        do ix=numtra+1,nta
          call readwinf(ix)
          call calcnetf(ix,0)
        enddo
c        call normphi(phia)
        call calcerr(phia,ntota+1,nta)
        write(*,*) "Q10P:", q10prh
	write(15,*) 
	write(15,*) "#Cross Validation"
	write(15,*) "#----------------"
	write(15,*) ie,q10prh
	write(15,*) 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      close(15)
      do i=1,nta
        v3sm = 1.d0 / (3.d0 + onu2p(i,1) + onu2p(i,2) + onu2p(i,3))
        v1 = v3sm * (1.d0 + onu2p(i,1))
        v2 = v3sm * (1.d0 + onu2p(i,2))
        v3 = v3sm * (1.d0 + onu2p(i,3))
c        write(99,499) i,phip(i),phia(i),v1,v2,v3
      enddo
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
399   format (I6,2x,D12.4,2x,D12.4,2x,D12.4)
499   format(I8,2x,D16.4,2x,D16.4,2x,D12.4,2x,D12.4,2x,D12.4)
      stop
      end
c_______________________________________________________________________
      
c------------- subroutine getdefs-------------------------------
c--
c* Get initial definitions
c--
      subroutine getdefs(fylevld,nte,npe,icv)
      implicit real*8(a-h,o-z)
      PARAMETER(nrdx=1300,nrpx=2640,naax=600000,nrx=21,nip=30,nw1=101,nw2=101,nrxo=1,nfilt=21)
      common/learn/ alpha, ratel, rmoment
      common/dcays/ dij1(nrx,nw1),dij2(nw2,nrxo),dfac,dijh(nw1,nw2),dijfi(nrx,nfilt),dijfo(nfilt,nrxo)
      common/misc/ rmae, ccoef, q10scr, q10prh, cvfrac, itsd
      character rtrainf(nrpx)*120
      common/sinout/ numtrf,ntotf,ntf,numtra,ntota,nta,idxprt(nrpx),rtrainf

      character fyle*6, fylein*120, fylevld*120
      character rtrain*120, rtraing(nrpx)*120

c* Definition of the general input/output files
      fyle='spiner'
      fylein=fyle//'.in'
      fylevld=fyle//'.valdav'
      open(unit=2,file=fylein)
c****  The following parameters are read from the input file progname.in
c*********** nte -- number of training epochs
      read(2,*)nte
c*********** npe -- size of protein set
      read(2,*)npe
c*********** ipe -- iterations per exemplar
      read(2,*)ipe
c*********** ratel -- learning rate
      read(2,*)ratel
c*********** rmoment -- initial momentum
      read(2,*)rmoment
c*********** alpha -- parameter for transfer function
      read(2,*)alpha
c*********** dfac -- parameter for distance decay matrix
      read(2,*)dfac
c*********** icv -- set for 10 fold cross validation
      read(2,*)icv
c*********** cvfrac -- fraction for cross validation
      read(2,*)cvfrac
c*********** rtrain -- file containing names of training proteins
      read(2,*)rtrain, fyle
      close(2)

      open(unit=3,file=rtrain)
c* read names of files containing training proteins
      i = 1
      do 245 while (i.gt.0)
        read(3,'(A)',end=246)rtraing(i)
	i = i + 1
245   continue
246   continue
      close(3)
      ntotp = i-1
      if (npe.le.ntotp) then
        ntotf = npe
      else
        write(*,*) "Protein set too small for number of proteins selected, using full set."
        ntotf = ntotp
      endif
      numval2 = nint(cvfrac*ntotf)
      ntf = ntotf
      ntotf = ntf-numval2
      numval = nint(0.05d0 * ntotf)
      numtrf = ntotf - numval

      call plst(ntf,ntotf,numval2,icv,idxprt)

      do i=1,ntf
        rtrainf(i) = rtraing(idxprt(i))
      enddo

      end
c___________________________________________________________________

c------------- subroutine plst-------------------------------
c--
c* Initialize the protein lists
c--
      subroutine plst(itot,ipt,ipv,indcv,lstt)
      implicit real*8(a-h,o-z)
      PARAMETER(nrdx=1300,nrpx=2640,naax=600000,nrx=21,nip=30,nw1=101,nw2=101,nrxo=1,nfilt=21)
      common/misc/ rmae, ccoef, q10scr, q10prh, cvfrac, itsd

      integer itot, ipt, ipv, indcv
      integer lstt(nrpx)
      
      iptv = ipt + ipv
      if (iptv.eq.itot) then
        ivf = nint(cvfrac*indcv*iptv)
        ivs = ivf - ipv
        jt = 1
        jv = ipt+1
        do i=1,iptv
          if ((ivs.lt.i).and.(i.le.ivf)) then
            lstt(jv) = i
            jv = jv + 1
          else
            lstt(jt) = i
            jt = jt + 1
          endif
        enddo
      else
        write(*,*) "Problem generating protein list"
        stop
      endif

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
      common/oldwei/ oldwi(nrx+1,nip,nw1),oldwh(nw1+1,nw2)
      common/pwei/ wo2(nw2+1,nrxo,3),oldwo2(nw2+1,nrxo,3)
      common/filtwei/ fwi(nrx+1,3,nfilt), fwo(nfilt+1,nrxo,3)
      common/ofiltwei/ ofwi(nrx+1,3,nfilt), ofwo(nfilt+1,nrxo,3)
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
	      oldwi(i,j,k) = atmp / 10.d0
	    enddo
	  enddo
	enddo

	do i=1,nw1+1
          do j=1,nw2
	    atmp = randnm() - 0.5
	    wh(i,j) = atmp
	    oldwh(i,j) = atmp / 10.d0
          enddo
	enddo

	do i=1,nw2+1
          do j=1,nrxo
            do k=1,3
	      atmp = randnm() - 0.5
   	      wo2(i,j,k) = atmp
	      oldwo2(i,j,k) = atmp / 10.d0
            enddo
          enddo
	enddo

	do i=1,nrx+1
	  do j=1,3
	    do k=1,nfilt
	      atmp = randnm() - 0.5
	      fwi(i,j,k) = atmp
	      ofwi(i,j,k) = atmp / 10.d0
	    enddo
	  enddo
	enddo
	do i=1,nfilt+1
          do j=1,3
	    atmp = randnm() - 0.5
	    fwo(i,1,j) = atmp
	    ofwo(i,1,j) = atmp / 10.d0
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

c------------- subroutine outwei -------------------------------
c--
c* Output weights
c--
      subroutine outwei(ishift)
      implicit real*8(a-h,o-z)
      PARAMETER(nrdx=1300,nrpx=2640,naax=600000,nrx=21,nip=30,nw1=101,nw2=101,nrxo=1,nfilt=21)
      common/weights/ wi(nrx+1,nip,nw1),wh(nw1+1,nw2)
      common/pwei/ wo2(nw2+1,nrxo,3),oldwo2(nw2+1,nrxo,3)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCC--- OUTPUT WEIGHTS  ---CCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      do i=1,nrx+1
	do j=1,nip
	  do k=1,nw1
	    write(ishift+77,*) wi(i,j,k)
	  enddo
	enddo
      enddo
      close(ishift+77)
      do i=1,nw1+1
        do j=1,nw2
	  write(ishift+78,*) wh(i,j)
        enddo
      enddo
      close(ishift+78)
      do i=1,nw2+1
        do j=1,nrxo
	  write(ishift+79,*) (wo2(i,j,k),k=1,3)
        enddo
      enddo
      close(ishift+79)

      end
c___________________________________________________________________

c------------- subroutine readwei -------------------------------
c--
c* Read optimized weights
c--
      subroutine readwei(ishift)
      implicit real*8(a-h,o-z)
      PARAMETER(nrdx=1300,nrpx=2640,naax=600000,nrx=21,nip=30,nw1=101,nw2=101,nrxo=1,nfilt=21)
      common/weights/ wi(nrx+1,nip,nw1),wh(nw1+1,nw2)
      common/pwei/ wo2(nw2+1,nrxo,3),oldwo2(nw2+1,nrxo,3)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCC--- INPUT WEIGHTS  ---CCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      do i=1,nrx+1
	do j=1,nip
	  do k=1,nw1
	    read(ishift+77,*) wi(i,j,k)
	  enddo
	enddo
      enddo
      close(ishift+77)
      do i=1,nw1+1
        do j=1,nw2
	  read(ishift+78,*) wh(i,j)
        enddo
      enddo
      close(ishift+78)
      do i=1,nw2+1
        do j=1,nrxo
	  read(ishift+79,*) (wo2(i,j,k),k=1,3)
        enddo
      enddo
      close(ishift+79)

      end
c___________________________________________________________________

c------------- subroutine outweif -------------------------------
c--
c* Output filter weights
c--
      subroutine outweif(ishift)
      implicit real*8(a-h,o-z)
      PARAMETER(nrdx=1300,nrpx=2640,naax=600000,nrx=21,nip=30,nw1=101,nw2=101,nrxo=1,nfilt=21)
      common/filtwei/ fwi(nrx+1,3,nfilt), fwo(nfilt+1,nrxo,3)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCC--- OUTPUT WEIGHTS  ---CCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      do i=1,nrx+1
	do j=1,3
	  do k=1,nfilt
	    write(ishift+87,*) fwi(i,j,k)
	  enddo
	enddo
      enddo
      close(ishift+87)
      do i=1,nfilt+1
        do j=1,3
	  write(ishift+88,*) fwo(i,1,j)
        enddo
      enddo
      close(ishift+88)

      end
c___________________________________________________________________

c------------- subroutine readweif -------------------------------
c--
c* Read optimized weights
c--
      subroutine readweif(ishift)
      implicit real*8(a-h,o-z)
      PARAMETER(nrdx=1300,nrpx=2640,naax=600000,nrx=21,nip=30,nw1=101,nw2=101,nrxo=1,nfilt=21)
      common/filtwei/ fwi(nrx+1,3,nfilt), fwo(nfilt+1,nrxo,3)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCC--- INPUT WEIGHTS  ---CCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      do i=1,nrx+1
	do j=1,3
	  do k=1,nfilt
	    read(ishift+87,*) fwi(i,j,k)
	  enddo
	enddo
      enddo
      close(ishift+87)
      do i=1,nfilt+1
        do j=1,3
	  read(ishift+88,*) fwo(i,1,j)
        enddo
      enddo
      close(ishift+88)

      end
c___________________________________________________________________

c------------- subroutine setio -------------------------------
c--
c* Setup the input/output arrays
c--
      subroutine setio()
      implicit real*8(a-h,o-z)
      PARAMETER(nrdx=1300,nrpx=2640,naax=600000,nrx=21,nip=30,nw1=101,nw2=101,nrxo=1,nfilt=21)
      common/filin/ protin(naax,nip+2),phir2(naax,3),phip(naax),phir(naax),dphir(naax)
      common/edges/ ibegia(naax),iendia(naax),ibegoa(naax),iendoa(naax), nisf, nosf
      character rtrainf(nrpx)*120
      common/sinout/ numtrf,ntotf,ntf,numtra,ntota,nta,idxprt(nrpx),rtrainf

      character tempc, ssc, aa1, aa2, aa3, predfl*120
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
      do ix=1,ntf
        ntao = nta
        if (ix.eq.(numtrf+1)) numtra = nta 
        if (ix.eq.(ntotf+1)) ntota = nta 
        predfl = "/data3/efaraggi/papers/seconds/new_paper1/old.db25/dssp_data/"//trim(rtrainf(ix))
        open(unit=44,file=predfl)
        predfl = "/data3/efaraggi/papers/seconds/new_paper1/old.db25/new_prf/"//trim(rtrainf(ix))//".mat"
        open(unit=45,file=predfl)
        read(45,*) ; read(45,*) ; read(45,*)
        predfl = "/data3/efaraggi/papers/seconds/new_paper1/old.db25/spXout/"//trim(rtrainf(ix))//".spXout"
        open(unit=46,file=predfl)
        read(46,*)
        i = 1
        do 265 while (i.gt.0)
          read(44,*,end=266,err=266) tempc,aa1,tempc,ssc
          read(45,*,end=266,err=266) tempc,aa2,(tprot(1,j),j=1,20)
          read(46,*,end=266,err=266) itemp,aa3,tempc,(tprot(2,j),j=1,8)
          nta = nta + 1

          if (aa1.ne.aa2) then
            write(*,*) "aa12 mismatch, aborting"
            write(*,*) "i=",i,"    nta=", nta, "    ix=", ix
            write(*,*) "aa1: ",aa1, "     aa2: ", aa2
            write(*,*) "file: ", rtrainf(ix)
            stop
          endif
          if (aa2.ne.aa3) then
            write(*,*) "aa23 mismatch, aborting"
            write(*,*) "i=",i,"    nta=", nta, "    ix=", ix
            write(*,*) "aa2: ",aa2, "     aa3: ", aa3
            write(*,*) "file: ", rtrainf(ix)
            stop
          endif
          
          iaa = ichar(aa1)
          do j=1,7
            protin(nta,j) = sevec(iaa,j)/0.9d0
          enddo

          do j=1,20
            protin(nta,7+j) = tprot(1,j)/9.d0
          enddo

          protin(nta,28) = 2.d0 * (tprot(2,8)/asafact(aa1) - 0.5d0)

          protin(nta,29) = tprot(2,1) / 180.d0
          tphi = (180.d0 + tprot(2,2)) / 360.d0
          if (tphi.gt.0.1389d0) then
            tphi = (tphi-0.1389d0)
          else 
            tphi = (tphi+0.8611d0)
          endif
          tphi = 2.d0*(tphi - 0.5d0)
          protin(nta,30) = tphi

          if (ssc.eq."H") then
            tphi = 0.5d0
          elseif (ssc.eq."C") then
            tphi = -0.2d0
          elseif (ssc.eq."E") then
            tphi = -0.6d0
          else
            write(*,*) "SS type not found, aborting."
            stop
          endif
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
      character rtrainf(nrpx)*120
      common/sinout/ numtrf,ntotf,ntf,numtra,ntota,nta,idxprt(nrpx),rtrainf

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
      common/pwei/ wo2(nw2+1,nrxo,3),oldwo2(nw2+1,nrxo,3)
      common/filout/ phia(naax),icphi(naax)
      common/inout/ resv(nrx+1,nip),phi2(nrxo,3)
      common/neurons/ hnu1(nw1+1),hnu2(nw2+1),onu2(naax,3),hnuf(nfilt+1),onuf(nrxo,3), finp(nrx,3)
      common/edges/ ibegia(naax),iendia(naax),ibegoa(naax),iendoa(naax), nisf, nosf
      common/dcays/ dij1(nrx,nw1),dij2(nw2,nrxo),dfac,dijh(nw1,nw2),dijfi(nrx,nfilt),dijfo(nfilt,nrxo)
      common/peaks/ peks(naax)
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
      character rtrainf(nrpx)*120
      common/sinout/ numtrf,ntotf,ntf,numtra,ntota,nta,idxprt(nrpx),rtrainf

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

c------------- subroutine backprop -------------------------------
c--
c* Calculate the neural network's output error given training output
c* phi. Backpropagate this error to the hidden neurons
c* and update the weights.
c* Subroutine also caluculates the RMSE as rmse
c--
      subroutine backprop(iwin)
      implicit real*8(a-h,o-z)
      PARAMETER(nrdx=1300,nrpx=2640,naax=600000,nrx=21,nip=30,nw1=101,nw2=101,nrxo=1,nfilt=21)
      common/weights/ wi(nrx+1,nip,nw1),wh(nw1+1,nw2)
      common/oldwei/ oldwi(nrx+1,nip,nw1),oldwh(nw1+1,nw2)
      common/pwei/ wo2(nw2+1,nrxo,3),oldwo2(nw2+1,nrxo,3)
      common/inout/ resv(nrx+1,nip),phi2(nrxo,3)
      common/neurons/ hnu1(nw1+1),hnu2(nw2+1),onu2(naax,3),hnuf(nfilt+1),onuf(nrxo,3), finp(nrx,3)
      common/learn/ alpha, ratel, rmoment
      common/edges/ ibegia(naax),iendia(naax),ibegoa(naax),iendoa(naax), nisf, nosf
      common/dcays/ dij1(nrx,nw1),dij2(nw2,nrxo),dfac,dijh(nw1,nw2),dijfi(nrx,nfilt),dijfo(nfilt,nrxo)
      integer iwin
      
      real*8 errh1(nrx,nw1),birh1(nw1),errh2(nw1,nw2),birh2(nw2),
     *biro2(nrxo,3),erro2(nw2,nrxo,3)
      
      ibegi = ibegia(iwin)
      iendi = iendia(iwin)
      ibego = ibegoa(iwin)
      iendo = iendoa(iwin)
c** Calculate errors at the output layer and the RMSE
      do i=ibego,iendo
        errotmp = phi2(i,1) - onu2(iwin,1)
        rmse = rmse + errotmp*errotmp
        biro2(i,1) = errotmp*alpha*(1.d0-onu2(iwin,1)*onu2(iwin,1))

        errotmp = phi2(i,2) - onu2(iwin,2)
        rmse = rmse + errotmp*errotmp
        biro2(i,2) = errotmp*alpha*(1.d0-onu2(iwin,2)*onu2(iwin,2))

        errotmp = phi2(i,3) - onu2(iwin,3)
        rmse = rmse + errotmp*errotmp
        biro2(i,3) = errotmp*alpha*(1.d0-onu2(iwin,3)*onu2(iwin,3))

        do j=1,nw2
          erro2(j,i,1) = biro2(i,1) * dij2(j,i)
          erro2(j,i,2) = biro2(i,2) * dij2(j,i)
          erro2(j,i,3) = biro2(i,3) * dij2(j,i)
        enddo
      enddo
     
c** Calculate the errors at the second hidden layer
      do k=1,nw2
        ztemp=0.d0
        do j=ibego,iendo
	  ztemp = ztemp + wo2(k,j,1)*erro2(k,j,1) + 
     *                    wo2(k,j,2)*erro2(k,j,2) + wo2(k,j,3)*erro2(k,j,3)
        enddo
        birh2(k) = ztemp * alpha  * (1 - hnu2(k) * hnu2(k))
        do j=1,nw1
	  errh2(j,k) = birh2(k) * dijh(j,k)
        enddo
      enddo

c** Calculate the errors at the first hidden layer
      do k=1,nw1
        ztemp=0.d0
        do j=1,nw2
	  ztemp = ztemp + wh(k,j)*errh2(k,j)
        enddo
        birh1(k) = ztemp * alpha  * (1 - hnu1(k) * hnu1(k))
        do j=ibegi,iendi
	  errh1(j,k) = birh1(k) * dij1(j,k)
        enddo

c** Update weights:
c** Update input/hidden layer weights (wi)
        do i=ibegi,iendi
          do j=1,nip
            ztemp = ratel*errh1(i,k)*resv(i,j) + rmoment*oldwi(i,j,k)
            wi(i,j,k) = wi(i,j,k) + ztemp
            oldwi(i,j,k) = ztemp
          enddo
        enddo

c** Update input bias weights (wi(nrx+1,1,1-nw1))
        ztemp = ratel*birh1(k) + rmoment*oldwi(nrx+1,1,k)
        wi(nrx+1,1,k) = wi(nrx+1,1,k) + ztemp
        oldwi(nrx+1,1,k) = ztemp
      enddo

c** Update hidden layer weights (wh)
      do k=1,nw2
        do j=1,nw1
          ztemp = ratel*errh2(j,k)*hnu1(j) + rmoment*oldwh(j,k)
          wh(j,k) = wh(j,k) + ztemp
          oldwh(j,k) = ztemp
        enddo

c** Update hidden bias weights (wh(nw1+1,1-nw2))
        ztemp = ratel*birh2(k) + rmoment*oldwh(nw1+1,k)
        wh(nw1+1,k) = wh(nw1+1,k) + ztemp
        oldwh(nw1+1,k) = ztemp

c** Update hidden/output layer weights (wo)
        do i=ibego,iendo
	  ztemp = ratel*erro2(k,i,1)*hnu2(k) + rmoment*oldwo2(k,i,1)
	  wo2(k,i,1) = wo2(k,i,1) + ztemp
	  oldwo2(k,i,1) = ztemp

	  ztemp = ratel*erro2(k,i,2)*hnu2(k) + rmoment*oldwo2(k,i,2)
	  wo2(k,i,2) = wo2(k,i,2) + ztemp
	  oldwo2(k,i,2) = ztemp

	  ztemp = ratel*erro2(k,i,3)*hnu2(k) + rmoment*oldwo2(k,i,3)
	  wo2(k,i,3) = wo2(k,i,3) + ztemp
	  oldwo2(k,i,3) = ztemp
        enddo
      enddo

c** Update hidden/output layer bias weights (wo)
        k=nw2+1
        do i=ibego,iendo
	  ztemp = ratel*biro2(i,1) + rmoment*oldwo2(k,i,1)
	  wo2(k,i,1) = wo2(k,i,1) + ztemp
	  oldwo2(k,i,1) = ztemp

	  ztemp = ratel*biro2(i,2) + rmoment*oldwo2(k,i,2)
	  wo2(k,i,2) = wo2(k,i,2) + ztemp
	  oldwo2(k,i,2) = ztemp

	  ztemp = ratel*biro2(i,3) + rmoment*oldwo2(k,i,3)
	  wo2(k,i,3) = wo2(k,i,3) + ztemp
	  oldwo2(k,i,3) = ztemp
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

c------------- subroutine backpropf -------------------------------
c--
c* Calculate the neural network's output error given training output
c* phi. Backpropagate this error to the hidden neurons
c* and update the weights.
c--
      subroutine backpropf(iwin)
      implicit real*8(a-h,o-z)
      PARAMETER(nrdx=1300,nrpx=2640,naax=600000,nrx=21,nip=30,nw1=101,nw2=101,nrxo=1,nfilt=21)
      common/filtwei/ fwi(nrx+1,3,nfilt), fwo(nfilt+1,nrxo,3)
      common/ofiltwei/ ofwi(nrx+1,3,nfilt), ofwo(nfilt+1,nrxo,3)
      common/inout/ resv(nrx+1,nip),phi2(nrxo,3)
      common/neurons/ hnu1(nw1+1),hnu2(nw2+1),onu2(naax,3),hnuf(nfilt+1),onuf(nrxo,3), finp(nrx,3)
      common/learn/ alpha, ratel, rmoment
      common/edges/ ibegia(naax),iendia(naax),ibegoa(naax),iendoa(naax), nisf, nosf
      common/dcays/ dij1(nrx,nw1),dij2(nw2,nrxo),dfac,dijh(nw1,nw2),dijfi(nrx,nfilt),dijfo(nfilt,nrxo)
      integer iwin
      
      real*8 errh1(nrx,nfilt),birh1(nfilt),biro2(nrxo,3),erro2(nfilt,nrxo,3)
      
      ibegi = ibegia(iwin)
      iendi = iendia(iwin)
      ibego = ibegoa(iwin)
      iendo = iendoa(iwin)
c** Calculate errors at the output layer and the RMSE
      do i=ibego,iendo
        errotmp = phi2(i,1) - onuf(i,1)
        biro2(i,1) = errotmp*alpha*(1.d0-onuf(i,1)*onuf(i,1))

        errotmp = phi2(i,2) - onuf(i,2)
        biro2(i,2) = errotmp*alpha*(1.d0-onuf(i,2)*onuf(i,2))

        errotmp = phi2(i,3) - onuf(i,3)
        biro2(i,3) = errotmp*alpha*(1.d0-onuf(i,3)*onuf(i,3))

        do j=1,nfilt
          erro2(j,i,1) = biro2(i,1) * dijfo(j,i)
          erro2(j,i,2) = biro2(i,2) * dijfo(j,i)
          erro2(j,i,3) = biro2(i,3) * dijfo(j,i)
        enddo
      enddo
     
c** Calculate the errors at the first hidden layer
      do k=1,nfilt
        ztemp=0.d0
        do j=ibego,iendo
	  ztemp = ztemp + fwo(k,j,1)*erro2(k,j,1) + 
     *                    fwo(k,j,2)*erro2(k,j,2) + fwo(k,j,3)*erro2(k,j,3)
        enddo
        birh1(k) = ztemp * alpha  * (1 - hnuf(k) * hnuf(k))
        do j=ibegi,iendi
	  errh1(j,k) = birh1(k) * dijfi(j,k)
        enddo

c** Update weights:
c** Update input/hidden layer weights (fwi)
        do i=ibegi,iendi
          do j=1,3
            ztemp = ratel*errh1(i,k)*finp(i,j) + rmoment*ofwi(i,j,k)
            fwi(i,j,k) = fwi(i,j,k) + ztemp
            ofwi(i,j,k) = ztemp
          enddo
        enddo

c** Update input bias weights (wi(nrx+1,1,1-nw1))
        ztemp = ratel*birh1(k) + rmoment*ofwi(nrx+1,1,k)
        fwi(nrx+1,1,k) = fwi(nrx+1,1,k) + ztemp
        ofwi(nrx+1,1,k) = ztemp
      enddo

c** Update output layer weights (fwo)
      do k=1,nfilt
        do i=ibego,iendo
          do j=1,3
	    ztemp = ratel*erro2(k,i,j)*hnuf(k) + rmoment*ofwo(k,i,j)
	    fwo(k,i,j) = fwo(k,i,j) + ztemp
	    ofwo(k,i,j) = ztemp
          enddo
        enddo
      enddo

c** Update hidden/output layer bias weights (wo)
        k=nw2+1
        do i=ibego,iendo
          do j=1,3
	    ztemp = ratel*biro2(i,j) + rmoment*ofwo(k,i,j)
	    fwo(k,i,j) = fwo(k,i,j) + ztemp
	    ofwo(k,i,j) = ztemp
          enddo
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
      character rtrainf(nrpx)*120
      common/sinout/ numtrf,ntotf,ntf,numtra,ntota,nta,idxprt(nrpx),rtrainf

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

c------------- subroutine getpks -------------------------------
c--
c* Get the peak prediction
c--
      subroutine getpks()
      implicit real*8(a-h,o-z)
      PARAMETER(nrdx=1300,nrpx=2640,naax=600000,nrx=21,nip=30,nw1=101,nw2=101,nrxo=1,nfilt=21)
      common/filin/ protin(naax,nip+2),phir2(naax,3),phip(naax),phir(naax),dphir(naax)
      common/filout/ phia(naax),icphi(naax)
      character rtrainf(nrpx)*120
      common/sinout/ numtrf,ntotf,ntf,numtra,ntota,nta,idxprt(nrpx),rtrainf
      common/peaks/ peks(naax)
      
c** Calculate peaks
      do i=1,nta
        if (phia(i).le.0.d0) then
          peks(i) = -0.5d0
        else
          peks(i) = 0.5d0
        endif
        dphir(i) = 2.d0*(phir(i)-phip(i)/2.d0)
      enddo

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
