c       ***************************************************************
c       CABSDOCK
c       authors: A. Kolinski(kolinski@chem.uw.edu.pl)
c                M. Kurcinski(mkurc@chem.uw.edu.pl)
c
c       updated: 19/4/2017
c
c       This version supports sg-sg restraints.
c       This code is part of the CABSdock2.0 package.
c
c       ***************************************************************

        IMPLICIT INTEGER(I-Z)
        LOGICAL LOOK,GOODC(800,800)

        character*3 aa(-1:20),NAME

        PARAMETER(NDIM=500)
        PARAMETER(NREPS=10)
        PARAMETER(NMOLS=2)
        PARAMETER(MAXRES=500)

c       The recommended
c       number of replicas (T1=T2=1.0) is about 20 for small
c       proteins up to 30 for N>200

        COMMON /CHAINS/  ICA(NDIM),X(NDIM),Y(NDIM),Z(NDIM)
        COMMON /THREE/   ICONF(800,800)
        COMMON /LENGHTS/ PROD(800,800)
        COMMON /SEQE/    SEQ(nmols,ndim),SEC(nmols,ndim)
        COMMON /VECTORS/ vx(800),vy(800),vz(800)
        COMMON /BETAS/   CBX(800,800),CBY(800,800),CBZ(800,800)
        COMMON /BISEC/   CAX(800,800),CAY(800,800),CAZ(800,800)
        COMMON /HB/      HBX(800,800),HBY(800,800),HBZ(800,800)

        common /arandom/ aarand,abrand,acrand,adrand
        COMMON /ENERGY/  EHBOND,ESC,EREP,BAT,
     &  EHBIJ(nmols,nmols,ndim,ndim)

        COMMON /pair/   apa(nmols,nmols,ndim,ndim,2,2),
     &                  app(nmols,nmols,ndim,ndim,2,2),
     &                  apm(nmols,nmols,ndim,ndim,2,2)
        COMMON /size/   arla(2,0:19,0:19,2,2),arlm(2,0:19,0:19,2,2),
     &                  arlp(2,0:19,0:19,2,2)
        COMMON /sizea/  ala(2,0:19,0:19,2,2),alm(2,0:19,0:19,2,2),
     &                  alp(2,0:19,0:19,2,2)
        
        COMMON /KBB/    kbin(800,800)
        COMMON /short/  IBIN(-500:500),asr(nmols,ndim,-12:12)
        COMMON /short1/ JBIN(800),bsr(nmols,ndim,16)
        COMMON /short0/ SBIN(200),csr(nmols,ndim,8)

        COMMON /one/    acrit(nmols),eoinp(0:19,0:100)
        common /onecm/  icmx(nmols),icmy(nmols),icmz(nmols)
        COMMON /FR/     FRG(nmols,ndim,ndim)
        COMMON /RCN/    ERESTA,NRESTA
        common /RRCA/   MRESA(nmols,ndim),KRESA(nmols,ndim,maxres),
     &  arca(nmols,ndim,maxres),awrca(nmols,ndim,maxres)
        common /cross/  mresax(nmols,ndim),kresax(nmols,ndim,maxres),
     &  nresax(nmols,ndim,maxres),arcax(nmols,ndim,maxres),
     &  awrcax(nmols,ndim,maxres)

        common /RSGm/ nRSG,eRSG
        common /RSG/  RSGcount(nmols,ndim),RSGind(nmols,ndim,maxres),
     &  distRSG(nmols,ndim,maxres),forceRSG(nmols,ndim,maxres)
        common /RSGX/ RSGXcount(nmols,ndim),RSGXind(nmols,ndim,maxres),
     &  RSGXmol(nmols,ndim,maxres),distRSGX(nmols,ndim,maxres),
     &  forceRSGX(nmols,ndim,maxres)
        
        COMMON /SG/ GX(800,800,0:19),GY(800,800,0:19),GZ(800,800,0:19)
        common /vec/ vector(-7:7,-7:7,-7:7)
        DIMENSION apabla(2,0:19,0:19,2,2),apablp(2,0:19,0:19,2,2)
        DIMENSION apablm(2,0:19,0:19,2,2)
        DIMENSION asrr(0:19,0:19,-12:12),bsrr(0:19,0:19,16)
        DIMENSION asrh(0:19,0:19,-12:12),bsrh(0:19,0:19,16)
        DIMENSION asre(0:19,0:19,-12:12),bsre(0:19,0:19,16)
        DIMENSION ATREP(NREPS),EREPLICA(NMOLS,NREPS),EREPLICAS(NREPS)
        DIMENSION REPNUM(NREPS)
        DIMENSION EINTER(nmols,nreps)
        common /reps/   xreps(ndim,NREPS,nmols),
     &                  yreps(ndim,NREPS,nmols),
     &                  zreps(ndim,NREPS,nmols)
        DIMENSION SWAPS(NREPS)
        DIMENSION axalf(0:19),ayalf(0:19),azalf(0:19)
        DIMENSION axbet(0:19),aybet(0:19),azbet(0:19)
        dimension icao(ndim),icab(ndim)

        common /lenfs/ mols,imol,itemp,lenf,lenf1,lenf2,lenf3,
     &  al2,al4,al5,al19,length(nmols)

        NFRAME=0

c       x,y,z - Ca lattice coordinates
c       [vx,vy,vz] 800 possible lattice Ca-Ca vectors
c       cax, cay, caz - normalyzed coordinates of (vi-vi+1)
c       cbx, cby, cbz - vectors from Ca to Cb (vi, vi+1 dependent)
c       gx,gy,gz - coordinates of the center of mas of sidechains
c       hbx,hby,hbz   - "hydrogen bond" vectors (vi, vi+1 dependent)
c       SEQ - sequence (LENF-2 numbers, 0=Gly, 1=Ala, ... 19=Trp)
c       SEC - secondary structure 1=coil, 2=helix, 3=turn, 4=beta
c       (coil and turn treated in the same way)
c       apa,apm,app (pairwise, orientation and position dependent
c       potentials, apabla, apablp, apablm- sequence independent)
c       arla, arlm, arlp - cut-off distances for pairwise interactions
c       axalf, ayalf, azalf - coeficients of the side group coordinates
c       in local othonormal coordinate system - for alpha type
c       main chain conformation (planar angle)
c       axbet, aybet, azbet - as above for an open angle of the backbone

c       FRG - short range distances deduced from SEC - i.e. helix and
c       beta target geometry
c       acrit - expected radius of gyration for a molecule

c       INPUT
c       IJKM  10  KKK 20          (RANDOM, ICYCLE, PHOT,  N_OF_REPLICAS,)
c       1.5  1.0  4.0             (T1,T2, REPULSIVE)
c       0.75  1.75  1.5 -1.5 0.375(gen, PAIRS, CENTRO, H_bonds,  SHORTR)
c       iiii, 0.5                 (#of_contacts,  scale)
c       jjjj, 0.5                 (#of_distances, scale)
c
c       Strenghts rescaled in this version

        data aa  /'BCK','GLY','ALA','SER','CYS','VAL',
     &                  'THR','ILE','PRO','MET','ASP',
     &                  'ASN','LEU','LYS','GLU','GLN',
     &                  'ARG','HIS','PHE','TYR','TRP','CYX'/

        OPEN(UNIT=7, FILE='SEQ',     STATUS='OLD')
        OPEN(UNIT=6, FILE='OUT')
        OPEN(UNIT=5, FILE='INP',     STATUS='OLD')
        OPEN(UNIT=10,FILE='FCHAINS', STATUS='OLD')
        OPEN(UNIT=9, FILE='TRAF')
        OPEN(UNIT=26,FILE='QUASI3S', STATUS='OLD')
        OPEN(UNIT=11,FILE='R14',     STATUS='OLD')
        OPEN(UNIT=12,FILE='R15',     STATUS='OLD')
        OPEN(UNIT=13,FILE='R13',     STATUS='OLD')
        OPEN(UNIT=29,FILE='CENTRO',  STATUS='OLD')
        OPEN(UNIT=33,FILE='SIDECENT',STATUS='OLD')

c       SEQ sequence   (i  ALA  SECONDARY_STRUCT_code=1,2,4)
c       INP input file with scaling parameters and restraints
c       ACHAINS - set of chains for Replicas, built by NatCA.f
c       QUASI3 - statistical pairwise potential, orient. dependent
c       R13, R14, R15 - short ramge statistical potentials
c       R14E, R14H, R15E, R15H are for known (predicted) fragments
c               of secondary structure, R14, R15 are generic
c
c       CENTRO - one-body centrosymmetric
c       PROFILE3 - Statistical, multibody, orientation dependent
c       SIDECENT - Coordinates of the side chains in respect to the
c       backbone, diferent for Helices and Expanded
c       PAIR3 - "homology" based, protein dependent pairwise

        do k=0,19
        read(33,*)
        read(33,*) axalf(k),ayalf(k),azalf(k),axbet(k),aybet(k),azbet(k)
        enddo
        CLOSE(33)

c       Centrosymmetric one body potential, energy of aminoacid
c       i in the shell j of a spherical potential, 3-bins up to
c       Acrit contain about 55% of aminoacids, the rest in two
c       additional bins up to 5*acrit/3
        do i=0,19
        read(29,3009) name, (eoinp(i,j), j=1,4)
        eoinp(i,0)=eoinp(i,1)
        eoinp(i,5)=eoinp(i,4)
        do j=6,100
        eoinp(i,j)=eoinp(i,j-1)+0.25
c       for larger distances energy grows in alinear fashion
        enddo
        enddo
3009    FORMAT(A3,4F6.2)
        CLOSE(29)

        READ(5,*) RANDOM, NCYCLE, PHOT, REPLICAS, MOLS
        READ(5,*) ATEMP2, ATEMP1, EREP, ARLLO, DTEM
        READ(5,*) ESC_0, ARLO, ENONE, EHBOND, ASRO

        do ii=0,19
        do j=0,100
        eoinp(ii,j)=enone*eoinp(ii,j)
        enddo
        enddo

        adrand=121.3421721
        abrand=0.3709112
        acrand=0.8021453
        aarand=float(random)/(12.3211+1.5321*float(RANDOM))

c       write(6,*) arand(ise)

        REWIND(11)

        do i=0,19
        do j=0,19
        read(11,*)
c       reading of "average" r14 interactions, it is chiral
        read(11,*) (asrr(i,j,k),k=-12,-5)
        read(11,*) (asrr(i,j,k),k=-4,3)
        read(11,*) (asrr(i,j,k),k=5,12)
        do k=4,1,-1
        asrr(i,j,k)=asrr(i,j,k-1)
        enddo
        enddo
        enddo

        do i=0,19
        do j=0,19
        read(12,*)
c       reading of "average" r15 interactions, k-type of aminoacid
        read(12,*) (bsrr(i,j,k),k=1,16)
c       read(12,*) (bsrr(i,j,k),k=1,8)
c       read(12,*) (bsrr(i,j,k),k=9,16)
        enddo
        enddo

        CLOSE(11)
        CLOSE(12)

        OPEN(UNIT=11,FILE='R14H',   STATUS='OLD')
        OPEN(UNIT=12,FILE='R15H',   STATUS='OLD')
        REWIND(11)
        REWIND(12)

c       When secondary structure known (SEC=2, or SEC=4) use
c       statistical potential thatis derived for specific secondary
c       structure elements H, or E  (read in two following segments)
        do i=0,19
        do j=0,19
        read(11,*)
        read(11,*) (asrh(i,j,k),k=-12,-5)
        read(11,*) (asrh(i,j,k),k=-4,3)
        read(11,*) (asrh(i,j,k),k=5,12)
        do k=4,1,-1
        asrh(i,j,k)=asrh(i,j,k-1)
        enddo
        enddo
        enddo

        do i=0,19
        do j=0,19
        read(12,*)
        read(12,*) (bsrh(i,j,k),k=1,16)
c       read(12,*) (bsrh(i,j,k),k=1,8)
c       read(12,*) (bsrh(i,j,k),k=9,16)
        enddo
        enddo

        CLOSE(11)
        CLOSE(12)

        OPEN(UNIT=11,FILE='R14E',   STATUS='OLD')
        OPEN(UNIT=12,FILE='R15E',   STATUS='OLD')

        REWIND(11)
        REWIND(12)

        do i=0,19
        do j=0,19
        read(11,*)
        read(11,*) (asre(i,j,k),k=-12,-5)
        read(11,*) (asre(i,j,k),k=-4,3)
        read(11,*) (asre(i,j,k),k=5,12)
        do k=4,1,-1
        asre(i,j,k)=asre(i,j,k-1)
        enddo
        enddo
        enddo

        do i=0,19
        do j=0,19
        read(12,*)
        read(12,*) (bsre(i,j,k),k=1,16)
c       read(12,*) (bsre(i,j,k),k=1,8)
c       read(12,*) (bsre(i,j,k),k=9,16)
        enddo
        enddo

        CLOSE(11)
        CLOSE(12)

        do i=1,800
        kk=int((sqrt(float(i))*0.61))+1
        if(kk.gt.16) kk=16
        JBIN(I) = kk
        ENDDO

        do i=1,500
        kk=int((sqrt(float(i))*0.61))+1
        if(kk.gt.12) kk=12
        IBIN(I) = kk
        IBIN(-I)=-kk
        ENDDO
        IBIN(0)=IBIN(1)

c       Pairwise interactions apablp ... and cut-off parmeters
c       arlp, orientation dependent, pairwise specific, sequence
c       independent  (two types of local conformations
c
c       Compact(H)  r13<6A
c       EXPANDED        r13>6A
c       Potentials  HH, EE, HE, EH

        DO K=1,2

        read(26,*)
        read(26,*)
        do i=0,19
        read(26,725) NAME, (apablp(k,i,j,1,1),j=0,19)
        enddo
        read(26,*)
        read(26,*)
        do i=0,19
        read(26,725) NAME, (apablp(k,i,j,2,2),j=0,19)
        enddo
        read(26,*)
        read(26,*)
        do i=0,19
        read(26,725) NAME, (apablp(k,i,j,1,2),j=0,19)
        enddo
        read(26,*)
        read(26,*)
        do i=0,19
        read(26,725) NAME, (apablp(k,i,j,2,1),j=0,19)
        enddo

        read(26,*)
        read(26,*)
        do i=0,19
        read(26,725) NAME, (apablm(k,i,j,1,1),j=0,19)
        enddo
        read(26,*)
        read(26,*)
        do i=0,19
        read(26,725) NAME, (apablm(k,i,j,2,2),j=0,19)
        enddo
        read(26,*)
        read(26,*)
        do i=0,19
        read(26,725) NAME, (apablm(k,i,j,1,2),j=0,19)
        enddo
        read(26,*)
        read(26,*)
        do i=0,19
        read(26,725) NAME, (apablm(k,i,j,2,1),j=0,19)
        enddo

        read(26,*)
        read(26,*)
        do i=0,19
        read(26,725) NAME, (apabla(k,i,j,1,1),j=0,19)
        enddo
        read(26,*)
        read(26,*)
        do i=0,19
        read(26,725) NAME, (apabla(k,i,j,2,2),j=0,19)
        enddo
        read(26,*)
        read(26,*)
        do i=0,19
        read(26,725) NAME, (apabla(k,i,j,1,2),j=0,19)
        enddo
        read(26,*)
        read(26,*)
        do i=0,19
        read(26,725) NAME, (apabla(k,i,j,2,1),j=0,19)
        enddo

        read(26,*)
        read(26,*)
        do i=0,19
        read(26,725) NAME, (arlp(k,i,j,1,1),j=0,19)
        enddo
        read(26,*)
        read(26,*)
        do i=0,19
        read(26,725) NAME, (arlp(k,i,j,2,2),j=0,19)
        enddo
        read(26,*)
        read(26,*)
        do i=0,19
        read(26,725) NAME, (arlp(k,i,j,1,2),j=0,19)
        enddo
        read(26,*)
        read(26,*)
        do i=0,19
        read(26,725) NAME, (arlp(k,i,j,2,1),j=0,19)
        enddo

        read(26,*)
        read(26,*)
        do i=0,19
        read(26,725) NAME, (arlm(k,i,j,1,1),j=0,19)
        enddo
        read(26,*)
        read(26,*)
        do i=0,19
        read(26,725) NAME, (arlm(k,i,j,2,2),j=0,19)
        enddo
        read(26,*)
        read(26,*)
        do i=0,19
        read(26,725) NAME, (arlm(k,i,j,1,2),j=0,19)
        enddo
        read(26,*)
        read(26,*)
        do i=0,19
        read(26,725) NAME, (arlm(k,i,j,2,1),j=0,19)
        enddo

        read(26,*)
        read(26,*)
        do i=0,19
        read(26,725) NAME, (arla(k,i,j,1,1),j=0,19)
        enddo
        read(26,*)
        read(26,*)
        do i=0,19
        read(26,725) NAME, (arla(k,i,j,2,2),j=0,19)
        enddo
        read(26,*)
        read(26,*)
        do i=0,19
        read(26,725) NAME, (arla(k,i,j,1,2),j=0,19)
        enddo
        read(26,*)
        read(26,*)
        do i=0,19
        read(26,725) NAME, (arla(k,i,j,2,1),j=0,19)
        enddo

        enddo

725     format(a3,1x,20f5.1)
        close(26)

C       ***************************************************************
C        PREPARATION OF THE CHAIN UNITS AND THEIR CORRELATION

        NWMAX=0
        DO ix=-7,7
        DO iy=-7,7
        DO iz=-7,7
        vector(ix,iy,iz)=0
        ir=ix*ix+iy*iy+iz*iz
        if(ir.ge.29) then
        if(ir.le.49) then
c       definition of Ca-Ca lattice vectors
        NWMAX=NWMAX+1
        VX(NWMAX)=ix
        VY(NWMAX)=iy
        VZ(NWMAX)=iz
        VECTOR(ix,iy,iz)=NWMAX
c       VECTOR()= 1,2,3....800
c       ix,iy,iz= -7, -6, -4,-3,  ...3,4,6,7 with  28< ||vx,vy,vz||<50
c       this is an arbitrary choice of the lattice geometry
        endif
        endif
        ENDDO
        ENDDO
        ENDDO

c       write(*,*) NWMAX

c       *******************   CREATING REPLICAS  ********************

        DO IMOL=1,MOLS

        do k=1,REPLICAS
        SWAPS(k)=0
        repnum(k)=k
        
        if(k.eq.1) then
        read(10,*) length(imol)
        else
        read(10,*)
        endif

        do i=1,length(imol)
        read(10,*) xreps(i,k,imol),yreps(i,k,imol),zreps(i,k,imol)
        enddo

        enddo

c       SEQUENCE INPUT

        do 121 i=2,length(imol)-1
        read(7,707) k, NAME, SEC(imol,i)
        if(sec(imol,i).eq.3) sec(imol,i)=1
        do j=0,20
        if(NAME.eq.aa(j)) then
        SEQ(imol,i)=j
c       sequence SEQ()=0 means GLY,  =19 means TRP
        go to 121
        endif
        enddo
121     continue

        ENDDO
707     format(i5,3x,a3,4x,i1)

c       EHBIJ - set-up secondary structure dependent
c       strength of the hyrogen bond network - stroger for helices
c       and beta-sheets

        EHSEC=EHBOND*1.5

        do imol=1,mols
        do jmol=1,mols

        do i=2,length(imol)-1
        is=sec(imol,i)
        do j=2,length(jmol)-1
        js=sec(jmol,j)

        if(imol.eq.jmol) then

        if(iabs(i-j).gt.2) then
        EHBIJ(imol,jmol,i,j)=0.75*EHBOND

        if(is.ne.2.AND.is*js.eq.4) EHBIJ(imol,jmol,i,j)=EHBOND
        if(iabs(i-j).eq.3.and.is.eq.2.and.js.eq.2)
     &  EHBIJ(imol,jmol,i,j)=EHSEC

        if(is.eq.4.and.js.eq.4) EHBIJ(imol,jmol,i,j)=EHSEC
        if(is*js.EQ.8) EHBIJ(imol,jmol,i,j)=0.0
        else
        EHBIJ(imol,jmol,i,j)=0.0
        endif

        else
        EHBIJ(imol,jmol,i,j)=0.75*EHBOND
        if(is.ne.2.AND.is*js.eq.4) EHBIJ(imol,jmol,i,j)=EHBOND
        if(is.eq.4.and.js.eq.4) EHBIJ(imol,jmol,i,j)=EHSEC
        if(is*js.EQ.8) EHBIJ(imol,jmol,i,j)=0.0
        endif

        enddo
        enddo

        enddo
        enddo

c       Combine the generic (asrr) short range potential with
c       secondary structure dependent potentials (asrh and asre)

        do imol=1,mols

c       R14 chiral potential

        do i=1,length(imol)-3
        do k=-12,12
        asr(imol,i,k)=ASRO*asrr(seq(imol,i+1),seq(imol,i+2),k)
        if(sec(imol,i+1).eq.2.AND.sec(imol,i+2).eq.2.AND.
     &  sec(imol,i+3).eq.2)then
        if(sec(imol,i).eq.2) then
        asr(imol,i,k)=ASRO*asrh(seq(imol,i+1),seq(imol,i+2),k)
        endif
        endif
        if(sec(imol,i+1).eq.4.AND.sec(imol,i+2).eq.4.AND.
     &  sec(imol,i+3).eq.4)then
        if(sec(imol,i).eq.4) then
        asr(imol,i,k)=ASRO*asre(seq(imol,i+1),seq(imol,i+2),k)
        endif
        endif
        enddo
        enddo

c       R15 potential - strength (ASRO) reduced as above

        do i=1,length(imol)-4
        do k=1,16
        bsr(imol,i,k)=ASRO*bsrr(seq(imol,i+1),seq(imol,i+3),k)
        if(sec(imol,i+1).eq.2.AND.sec(imol,i+2).eq.2.AND.
     &  sec(imol,i+3).eq.2)then
        if(sec(imol,i).eq.2.AND.sec(imol,i+4).eq.2) then
        bsr(imol,i,k)=ASRO*bsrh(seq(imol,i+1),seq(imol,i+3),k)
        endif
        endif
        if(sec(imol,i+1).eq.4.AND.sec(imol,i+2).eq.4.AND.
     &  sec(imol,i+3).eq.4)then
        if(sec(imol,i).eq.4.AND.sec(imol,i+4).EQ.4) then
        bsr(imol,i,k)=ASRO*bsre(seq(imol,i+1),seq(imol,i+3),k)
        endif
        endif
        enddo
        enddo

        enddo

        REWIND(13)
        do i=1,200
        kk=int((sqrt(float(i))*0.61))+1
        if(kk.gt.8) kk=8
        SBIN(I) = kk
        ENDDO

        REWIND(13)
        do i=0,19
        do j=0,19
        read(13,*)
c       reading of "average" r13 interactions, k-type of aminoacid
        read(13,*) (bsrr(i,j,k),k=1,8)
        enddo
        enddo

        CLOSE(13)
        OPEN(UNIT=13,FILE='R13H',   STATUS='OLD')
        REWIND(13)

        do i=0,19
        do j=0,19
        read(13,*)
c       reading of "average" r13 interactions, k-type of aminoacid
        read(13,*) (bsrh(i,j,k),k=1,8)
        enddo
        enddo

        CLOSE(13)
        OPEN(UNIT=13,FILE='R13E',   STATUS='OLD')
        REWIND(13)

        do i=0,19
        do j=0,19
        read(13,*)
c       reading of "average" r13 interactions, k-type of aminoacid
        read(13,*) (bsre(i,j,k),k=1,8)
        enddo
        enddo

        CLOSE(13)

c       R13 potential - strength ASRO

        do imol=1,mols

        do i=2,length(imol)-3
        do k=1,8
        csr(imol,i,k)=ASRO*bsrr(seq(imol,i),seq(imol,i+2),k)
        if(sec(imol,i).eq.2.AND.sec(imol,i+1).eq.2.AND.
     &  sec(imol,i+2).eq.2)then
        csr(imol,i,k)=(csr(imol,i,k)+ASRO*bsrh(seq(imol,i),
     &  seq(imol,i+2),k))/2.0
        endif
        if(sec(imol,i).eq.4.AND.sec(imol,i+1).eq.4.AND.
     &  sec(imol,i+2).eq.4)then
        csr(imol,i,k)=(csr(imol,i,k)+ASRO*bsre(seq(imol,i),
     &  seq(imol,i+2),k))/2.0
        endif
        enddo
        enddo

        do k=1,8
        csr(imol,1,k)=0.0
        csr(imol,length(imol)-2,k)=0.0
        enddo

        enddo

        do imol=1,mols
        do jmol=1,mols

        do i=2,length(imol)-1
        ii=SEQ(imol,i)
        do j=2,length(jmol)-1
        jj=SEQ(jmol,j)

        dd=1.0

        if(imol.eq.jmol) then
        k=1
        aarlo=arlo
        if(iabs(i-j).eq.4.OR.iabs(i-j).EQ.6) dd=0.0
        else
        k=2
        aarlo=arlo*arllo
        endif

        do k1=1,2
        do k2=1,2
        apa(imol,jmol,i,j,k1,k2)=aarlo*(apabla(k,ii,jj,k1,k2)-dd)
        apm(imol,jmol,i,j,k1,k2)=aarlo*(apablm(k,ii,jj,k1,k2)-dd)
        app(imol,jmol,i,j,k1,k2)=aarlo*(apablp(k,ii,jj,k1,k2)-dd)
        enddo
        enddo

        enddo
        enddo

        enddo
        enddo

        d1=1.5
        d2=2.0

c       d1+d2 the width of the square well pairwise potential

        do k=1,2
        do k1=1,2
        do k2=1,2
        do i=0,19
        do j=0,19
        d11=d1
        d22=d2

        ala(k,i,j,k1,k2)=((arla(k,i,j,k1,k2)+d11)/0.61)**2
        alm(k,i,j,k1,k2)=((arlm(k,i,j,k1,k2)+d11)/0.61)**2
        alp(k,i,j,k1,k2)=((arlp(k,i,j,k1,k2)+d11)/0.61)**2
        if(apabla(k,i,j,k1,k2).gt.0.)
     &  ala(k,i,j,k1,k2)=((arla(k,i,j,k1,k2))/0.61)**2
        if(apablm(k,i,j,k1,k2).gt.0.)
     &  alm(k,i,j,k1,k2)=((arlm(k,i,j,k1,k2))/0.61)**2
        if(apablp(k,i,j,k1,k2).gt.0.)
     &  alp(k,i,j,k1,k2)=((arlp(k,i,j,k1,k2))/0.61)**2
c
        arla(k,i,j,k1,k2)=(max(2.,(arla(k,i,j,k1,k2)-d22))/0.61)**2
        arlm(k,i,j,k1,k2)=(max(2.,(arlm(k,i,j,k1,k2)-d22))/0.61)**2
        arlp(k,i,j,k1,k2)=(max(2.,(arlp(k,i,j,k1,k2)-d22))/0.61)**2

        enddo
        enddo
        enddo
        enddo
        enddo

c       Ca-Ca distance restraints

        do i=1,mols
        do j=1,length(i)
        mresa(i,j)=0
        mresax(i,j)=0
        enddo
        enddo

        read(5,*) NRESTA,ERESTA
        if(NRESTA.gt.0) then
        do k=1,NRESTA
        read(5,*) imol,irest,jmol,jrest,ar,aw
        ar=ar/0.61
        i=irest+1
        j=jrest+1
        if(imol.eq.jmol) then
        ii=mresa(imol,i)+1
        jj=mresa(jmol,j)+1
        mresa(imol,i)=ii
        mresa(jmol,j)=jj
        kresa(imol,i,ii)=j
        kresa(jmol,j,jj)=i
        arca(imol,i,ii)=ar
        arca(jmol,j,jj)=ar
        awrca(imol,i,ii)=aw
        awrca(jmol,j,jj)=aw
        else 
        ii=mresax(imol,i)+1
        jj=mresax(jmol,j)+1
        mresax(imol,i)=ii
        mresax(jmol,j)=jj
        kresax(imol,i,ii)=j
        kresax(jmol,j,jj)=i
        nresax(imol,i,ii)=jmol
        nresax(jmol,j,jj)=imol
        arcax(imol,i,ii)=ar
        arcax(jmol,j,jj)=ar
        awrcax(imol,i,ii)=aw
        awrcax(jmol,j,jj)=aw
        endif
        enddo
        endif

c       sg-sg restraints

        read(5,*) nRSG,eRSG
        if(nRSG.gt.0) then
        do k=1,nRSG
        read(5,*) imol,irest,jmol,jrest,ar,aw
        ar=ar/0.61
        i=irest+1
        j=jrest+1
        if(imol.eq.jmol) then
        ii=RSGcount(imol,i)+1
        jj=RSGcount(jmol,j)+1
        RSGcount(imol,i)=ii
        RSGcount(jmol,j)=jj
        RSGind(imol,i,ii)=j
        RSGind(jmol,j,jj)=i
        distRSG(imol,i,ii)=ar
        distRSG(jmol,j,jj)=ar
        forceRSG(imol,i,ii)=aw
        forceRSG(jmol,j,jj)=aw
        else
        ii=RSGXcount(imol,i)+1
        jj=RSGXcount(jmol,j)+1
        RSGXcount(imol,i)=ii
        RSGXcount(jmol,j)=jj
        RSGXind(imol,i,ii)=j
        RSGXind(jmol,j,jj)=i
        RSGXmol(imol,i,ii)=jmol
        RSGXmol(jmol,j,jj)=imol
        distRSGX(imol,i,ii)=ar
        distRSGX(jmol,j,jj)=ar
        forceRSGX(imol,i,ii)=aw
        forceRSGX(jmol,j,jj)=aw
        endif
        enddo
        endif

        do imol=1,mols
        sec(imol,1)=1
        sec(imol,length(imol))=1

c       Set a generalized shell model- Acrit is radius of gyration
c       the equation is semiempirical regresion against chain length AL2
c       0.61 is the scaling to the lattice

        Acrit(imol)=2.2*exp(0.38*alog(length(imol)-2.0001))/0.61

        enddo

c       800 vectors representing possible representations of the
c       alpha carbon virtual bonds.  1 lattice units = 0.61 A

c       ===============================================================

c       Prepare the set of good pairs of vectors - exclude too narrow
c       pairs and the colinear ones - to enable Cb definition

        DO I=1,800
        ix=vx(i)
        iy=vy(i)
        iz=vz(i)

        DO J=1,800
        goodc(i,j)=.TRUE.
        jx=vx(j)
        jy=vy(j)
        jz=vz(j)
        prod(i,j)= ix*jx+iy*jy+iz*jz
        ICONF(i,j)=(ix+jx)**2+(iy+jy)**2+(iz+jz)**2
        if(iconf(i,j).lt.45) GOODC(i,j)=.FALSE.
        if(iconf(i,j).gt.145) GOODC(i,j)=.FALSE.

        kx=iy*jz-iz*jy
        ky=jx*iz-ix*jz
        kz=ix*jy-iy*jx
        product=(kx*kx+ky*ky+kz*kz)
        IF(product.EQ.0) goodc(I,J)=.FALSE.
c
c       Prepare Cb positions - this is temporary (USE A FILE INSTEAD)
c       ideal tetrahedral conformation of Ca assumed -  should be OK

        IF(GOODC(i,j)) THEN

        if(iconf(i,j).lt.95) then
        kbin(i,j)=1
        else
        kbin(i,j)=2
        endif

        cx=ix+jx
        cy=iy+jy
        cz=iz+jz
        a=sqrt(cx*cx+cy*cy+cz*cz)
        cx=cx/a
        cy=cy/a
        cz=cz/a

        ax=ix-jx
        ay=iy-jy
        az=iz-jz

        a=sqrt(ax*ax+ay*ay+az*az)
        ax=ax/a
        ay=ay/a
        az=az/a

        CAX(i,j)=ax
        CAY(i,j)=ay
        CAZ(i,j)=az

c       Here define also the H-bond (Ca-based vectors)- 4.5A

        b=sqrt(float(product))
        bx=float(kx)/b
        by=float(ky)/b
        bz=float(kz)/b

        HBX(i,j)=8.25*bx
        HBY(i,j)=8.25*by
        HBZ(i,j)=8.25*bz

        DO k=0,19
        IF(iconf(i,j).LT.75) THEN
c       write down the helical Side Group positions
c       for a helical and turn-like conformations
        GX(i,j,k)=(axalf(k)*cx+ayalf(k)*bx+azalf(k)*ax)/0.61
        GY(i,j,k)=(axalf(k)*cy+ayalf(k)*by+azalf(k)*ay)/0.61
        GZ(i,j,k)=(axalf(k)*cz+ayalf(k)*bz+azalf(k)*az)/0.61
        ELSE
c       for expanded conformations
        IF(iconf(i,j).gt.110) THEN
        GX(i,j,k)=(axbet(k)*cx+aybet(k)*bx+azbet(k)*ax)/0.61
        GY(i,j,k)=(axbet(k)*cy+aybet(k)*by+azbet(k)*ay)/0.61
        GZ(i,j,k)=(axbet(k)*cz+aybet(k)*bz+azbet(k)*az)/0.61
        else
c       an approximation of the orientational averaging-consider fixing it 
c       (75-110)  is the intermediate range of iconf
        fb=float(iconf(i,j)-75)/35.0
        fa=1-fb
        axb=axbet(k)*fb+fa*axalf(k)
        ayb=aybet(k)*fb+fa*ayalf(k)
        azb=azbet(k)*fb+fa*azalf(k)
        GX(i,j,k)=(axb*cx+ayb*bx+azb*ax)/0.61
        GY(i,j,k)=(axb*cy+ayb*by+azb*ay)/0.61
        GZ(i,j,k)=(axb*cz+ayb*bz+azb*az)/0.61
        endif
        ENDIF

        ENDDO
        CBX(i,j)=GX(i,j,1)
        CBY(i,j)=GY(i,j,1)
        CBZ(i,j)=GZ(i,j,1)
        ENDIF
        ENDDO

        ENDDO

c       COMPUTE THE SECONDARY FRAGMENT BIASES

        do imol=1,mols
        do i=1,length(imol)
        do j=1,length(imol)
        FRG(imol,j,i)=0.0
        enddo
        enddo

        do i=2,length(imol)-7
        q=0
        do j=i+1,i+5
        if(sec(imol,j).eq.4) q=q+4
        enddo
c	length of six bond beta-type fragment
        if(q.eq.20) FRG(imol,i,i+6)=19.1/0.61
        enddo

        do i=2,length(imol)-8
        q=0
        do j=i,i+7
        if(sec(imol,j).eq.2) q=q+2
        enddo
c	length of eight bond Helix-type fragment
        if(q.eq.16) FRG(imol,i,i+7)=10.75/0.61
        enddo
        enddo

c       *****************************************************

c	THE MOST EXTERNAL LOOP FOR THE MC DYNAMICS
c	twenty temperatures for the annealing - setting 
c	atemp2=atemp1 reduces program to a more standard REMC
c	dtem1 -annealing increment
        
        dtem1=(atemp2-atemp1)/20.
c       dtem=5.8/float(lenf)
c       dtem=0.1
c	dtem distance between replicas
c   dtem read from INP
        aatemp2=atemp2

        DO IDDDUM=1,20
        
        aatemp1=aatemp2-dtem1
        atemp=aatemp1

        do i=1, REPLICAS
        ATREP(i)=atemp
C	atemp=atemp+dtem*((float(i-1+REPLICAS))/float(REPLICAS))**3
        atemp=atemp+dtem
c	atemp temperature of a replica, lowest replica temp=aatemp1
        EREPLICAS(i)=0.0
        do imol=1,mols
        EREPLICA(imol,i)=0.0
        enddo
        enddo

        aatemp2=aatemp1

        DO 7000 ICYCLE=1,NCYCLE

        IF(ICYCLE.GT.1) THEN
        If(REPLICAS.GT.1) THEN

        DO ITEMP=REPLICAS-1,1,-1
        DO ITEMP1=REPLICAS,ITEMP+1,-1

c	exchanging replicas

        betam=1.0/ATREP(itemp)
        betan=1.0/ATREP(itemp1)
        delta=(betan-betam)*(EREPLICAS(itemp)-EREPLICAS(itemp1))

        if(exp(-delta).gt.arand(ise)) then
        rnum=repnum(itemp)
        repnum(itemp)=repnum(itemp1)
        repnum(itemp1)=rnum
        DO imol=1,mols
        do i=1,length(imol)
        xt=xreps(i,itemp,imol)
        yt=yreps(i,itemp,imol)
        zt=zreps(i,itemp,imol)
        xreps(i,itemp,imol)=xreps(i,itemp1,imol)
        yreps(i,itemp,imol)=yreps(i,itemp1,imol)
        zreps(i,itemp,imol)=zreps(i,itemp1,imol)
        xreps(i,itemp1,imol)=xt
        yreps(i,itemp1,imol)=yt
        zreps(i,itemp1,imol)=zt
        enddo
        attt=EREPLICA(imol,itemp)
        EREPLICA(imol,itemp)=EREPLICA(imol,itemp1)
        EREPLICA(imol,itemp1)=attt
        ENDDO
        SWAPS(itemp)=SWAPS(itemp)+1
        SWAPS(itemp1)=SWAPS(itemp1)+1
        attt=EREPLICAS(itemp)
        EREPLICAS(itemp)=EREPLICAS(itemp1)
        EREPLICAS(itemp1)=attt
        endif
        
        ENDDO
        ENDDO

        ENDIF
        ENDIF

c	******************* ITERATE THE REPLICAS ********************

        DO ITEMP=REPLICAS,1,-1

        ATEMP=ATREP(itemp)

c	SCALING UP AT HIGH TEMP THE SHORT RANGE INTERACTIONS (GENERIC)
c	ESC=ESC_0 *sqrt(atemp/aatemp1)*(aatemp1/ATEMP1)
c       ESC=ESC_0 *sqrt(aatemp1/ATEMP1)
        ESC=ESC_0

c	a new cycle starts here (internal "time" loop)

        DO 7777 IPHOT=1, PHOT

c	CENTER THE SYSTEM at [0,0,0] point of the lattice

        sx=0
        sy=0
        sz=0
        lenfsum=0
        do imol=1,mols
        ssx=0
        ssy=0
        ssz=0
        do i=1,length(imol)
        ssx=ssx+xreps(i,itemp,imol)
        ssy=ssy+yreps(i,itemp,imol)
        ssz=ssz+zreps(i,itemp,imol)
        enddo
        sx=sx+ssx
        sy=sy+ssy
        sz=sz+ssz
        lenfsum=lenfsum+length(imol)
        icmx(imol)=nint(ssx/float(length(imol)))
        icmy(imol)=nint(ssy/float(length(imol)))
        icmz(imol)=nint(ssz/float(length(imol)))
        enddo
        sx=nint(sx/float(lenfsum))
        sy=nint(sy/float(lenfsum))
        sz=nint(sz/float(lenfsum))
        do imol=1,mols
        do i=1,length(imol)
        xreps(i,itemp,imol)=xreps(i,itemp,imol)-sx
        yreps(i,itemp,imol)=yreps(i,itemp,imol)-sy
        zreps(i,itemp,imol)=zreps(i,itemp,imol)-sz
        enddo
        icmx(imol)=icmx(imol)-sx
        icmy(imol)=icmy(imol)-sy
        icmz(imol)=icmz(imol)-sz
        enddo

        DO 9999 IMOL=1,MOLS

        sx=0
        sy=0
        sz=0
        lenfsum=0

        do i=1,mols
         if(i.ne.imol) then
          sx=sx+icmx(i)*length(i)
          sy=sy+icmy(i)*length(i)
          sz=sz+icmz(i)*length(i)
          lenfsum=lenfsum+length(i)
         endif        
        enddo

        sx=nint(sx/float(lenfsum))
        sy=nint(sy/float(lenfsum))
        sz=nint(sz/float(lenfsum))

        ssx=sx-icmx(imol)
        ssy=sy-icmy(imol)
        ssz=sz-icmz(imol)

        assr=sqrt(float(sx*ssx+ssy*ssy+ssz*ssz))

        if(assr.gt.82) then
         assr=(assr-82)/assr
         ssx=nint(assr*ssx)
         ssy=nint(assr*ssy)
         ssz=nint(assr*ssz)

         do i=1,length(imol)
          xreps(i,itemp,imol)=xreps(i,itemp,imol)+ssx
          yreps(i,itemp,imol)=yreps(i,itemp,imol)+ssy
          zreps(i,itemp,imol)=zreps(i,itemp,imol)+ssz
         enddo
        endif

        lenf=length(imol)
        lenf1=lenf-1
        lenf2=lenf-2
        lenf3=lenf-3
        al2=lenf-2.0001
        al4=lenf-4.0001
        al5=lenf-5.0001
        al19=lenf-19.0001
        if(al19.lt.2) al19=2

        DO i=1,lenf
        x(i)=xreps(i,itemp,imol)
        y(i)=yreps(i,itemp,imol)
        z(i)=zreps(i,itemp,imol)
        enddo

        DO I=1,LENF1
        J=I+1
        WX=X(J)-X(I)
        WY=Y(J)-Y(I)
        WZ=Z(J)-Z(I)
        ICA(I)=VECTOR(WX,WY,WZ)
        ENDDO
        ICA(LENF)=1

        energ=EHB(2,lenf1)+ESHORT(2,lenf1)+EHL(2,lenf1)
     &  +ERESTR(2,lenf1)+ECROSS(2,lenf1)

        DO 7775 IDUM=1,LENF2
c
c	TWO BOND KINKS- THE SMALLEST MOVES ATTEMPTED 10-times 
c	more frequently

        DO 7774 IONEBEAD=1,10

        j=INT(arand(ise)*AL4)+3
        i=j-1
        ii=ica(i)
        jj=ica(j)
        k=J+1

3307    ix=int(arand(ise)*6.9999)-3 +x(j)
        iy=int(arand(ise)*6.9999)-3 +y(j)
        iz=int(arand(ise)*6.9999)-3 +z(j)
        ir=(x(i)-ix)**2+(y(i)-iy)**2+(z(i)-iz)**2
        if(ir.gt.49) go to 3307
        if(ir.lt.29) go to 3307
        ir=(x(k)-ix)**2+(y(k)-iy)**2+(z(k)-iz)**2
        if(ir.gt.49) go to 3307
        if(ir.lt.29) go to 3307

        nv1=vector((-x(i)+ix),(-y(i)+iy),(-z(i)+iz))
        if(nv1.eq.ii) go to 7774
        nv2=vector((x(k)-ix),(y(k)-iy),(z(k)-iz))
        if(.NOT.GOODC(nv1,nv2)) GO to 3307

c	nv1, nv2 - possible new values of ICA(j-1),ICA(j)
c	j- is the moved Ca unit
        
        ica1=ica(i-1)
c	check planar angles with the rest of the chain
        if(.NOT.GOODC(ica1,nv1)) go to 3307
        ica3=ica(k)
        if(.NOT.GOODC(nv2,ica3)) go to 3307

        px=X(i)+vx(nv1)
        py=Y(i)+vy(nv1)
        pz=Z(i)+vz(nv1)

        kx=x(j)
        ky=y(j)
        kz=z(j)
        ICA(I)=nv1
        ICA(j)=nv2
        x(j)=px
        y(j)=py
        z(j)=pz

        IF(LOOK(i,k)) THEN
c	No steric overlaps of Ca and Cb detected for the new conformer
c	THE MOVE COULD BE SUCESSFULL. APPLY METROPOLIS CRITERION HERE

        ENEW=EHB(i,k)+ESHORT(i,k)+EHL(i,k)+ERESTR(i,k)+ECROSS(i,k)

        ICA(I)=ii
        ICA(j)=jj
        x(j)=kx
        y(j)=ky
        z(j)=kz

        EOLD=EHB(i,k)+ESHORT(i,k)+EHL(i,k)+ERESTR(i,k)+ECROSS(i,k)

        DE=ENEW-EOLD
        IF(DE.GT.0.0) THEN
        if(arand(ise).gt.EXP(-DE/ATEMP)) then
        go to 7774
        endif
        ENDIF

c	THE MOVE SUCESSFULL. WRITE-IN THE NEW CONFORMATION

        ENERG=ENERG+DE
        ICA(I)=nv1
        ICA(j)=nv2
        x(j)=px
        y(j)=py
        z(j)=pz

        go to 7774
        ELSE

c	Move rejected - restore the initial state of the chain		
        ICA(I)=ii
        ICA(j)=jj
        x(j)=kx
        y(j)=ky
        z(j)=kz
        ENDIF

7774    CONTINUE

c	THREE-BOND MOVE -permutation

        I=INT(arand(ise)*AL5)+2
        a=arand(ise)

        ii=ica(i)
        j=i+1
        jj=ica(j)
        k=i+2
        kk=ica(k)
        l=i+3

        if(arand(ise).gt.0.5) then
        nv1=kk
        nv2=jj
        nv3=ii
        if(nv1.eq.ii) go to 8885
        else
        if(arand(ise).gt.0.5) then
        nv1=kk
        nv2=ii
        nv3=jj
        if(nv1.eq.ii) go to 8885
        else
        nv1=jj
        nv2=kk
        nv3=ii
        if(nv1.eq.ii) go to 8885
        endif
        endif

        ica1=ica(i-1)
        if(.NOT.GOODC(ica1,nv1)) go to 8885
        if(.NOT.GOODC(nv1,nv2)) go to 8885
        if(.NOT.GOODC(nv2,nv3)) go to 8885
        if(.NOT.GOODC(nv3,ica(l))) go to 8885

c	Store the old conformation

        jx=x(j)
        jy=y(j)
        jz=z(j)

        kx=x(k)
        ky=y(k)
        kz=z(k)

        NX=x(i)+vx(nv1)
        NY=y(i)+vy(nv1)
        NZ=z(i)+vz(nv1)
        MX=NX+vx(nv2)
        MY=NY+vy(nv2)
        MZ=NZ+vz(nv2)
        X(j)=NX
        Y(j)=NY
        Z(J)=NZ
        X(k)=MX
        Y(k)=MY
        Z(k)=MZ
        ICA(i)=NV1
        ICA(j)=NV2
        ICA(k)=NV3

        IF(LOOK(i,l)) THEN
c	No overalps of Cas and Cbs
c	THE MOVE COULD BE SUCESSFULL. APPLY METROPOLIS CRITERION HERE

        ENEW=EHB(i,l)+ESHORT(i,l)+EHL(i,l)+ERESTR(i,l)+ECROSS(i,l)

        ICA(i)=ii
        ICA(j)=jj
        ICA(k)=kk
        x(j)=jx
        y(j)=jy
        z(j)=jz
        x(k)=kx
        y(k)=ky
        z(k)=kz

        EOLD=EHB(i,l)+ESHORT(i,l)+EHL(i,l)+ERESTR(i,l)+ECROSS(i,l)

        DE=ENEW-EOLD

        IF(DE.GT.0.0) THEN
        if(arand(ise).gt.EXP(-DE/ATEMP)) then
        go to 8885
        endif
        ENDIF

c	THE MOVE SUCESSFULL. WRITE-IN THE NEW CONFORMATION

        ENERG=ENERG+DE
        X(j)=NX
        Y(j)=NY
        Z(J)=NZ
        X(k)=MX
        Y(k)=MY
        Z(k)=MZ
        ICA(i)=NV1
        ICA(j)=NV2
        ICA(k)=NV3

        go to 8885

        ELSE

c	Move rejected - restore the initial state of the chain			
        ICA(i)=ii
        ICA(j)=jj
        ICA(k)=kk
        x(j)=jx
        y(j)=jy
        z(j)=jz
        x(k)=kx
        y(k)=ky
        z(k)=kz

        ENDIF

8885    CONTINUE

c	"REPTATION" -type moves, two bond progress (4 to 22 bonds)

        j=INT(arand(ise)*AL19)+3
        ppp=INT(arand(ise)*18.9999)+2
        i=j-1
        k=j+ppp
        L=k+1
        if(L.gt.lenf1) go to 4001

        do ii=i,k
        icao(ii)=ica(ii)
        enddo

        IF(arand(ise).gt.0.5) THEN
c	a move forward
        if(.not.GOODC(ica(l-3),ica(l))) go to 4001
        if(.not.GOODC(ica(i-1),ica(k-1))) go to 4001
        if(.not.GOODC(ica(k),ica(i))) go to 4001
        ii=ica(k-1)
        kk=ica(k)
        do pp=k,j+1,-1
        ica(pp)=ica(pp-2)
        enddo
        ica(i)=ii
        ica(j)=kk
        ELSE

c	a move backward
        if(.not.GOODC(ica(i-1),ica(j+1))) go to 4001
        if(.not.GOODC(ica(k),ica(i))) go to 4001
        if(.not.GOODC(ica(j),ica(l))) go to 4001
        ii=ica(i)
        kk=ica(j)
        do pp=i,k-2
        ica(pp)=ica(pp+2)
        enddo
        ica(k-1)=ii
        ica(k)=kk

        ENDIF

        do kkkk=i,k
        icab(kkkk)=ica(kkkk)
        enddo
        do kkkk=j,k
        pp=kkkk-1
        ii=ica(pp)
        x(kkkk)=x(pp)+vx(ii)
        y(kkkk)=y(pp)+vy(ii)
        z(kkkk)=z(pp)+vz(ii)
        enddo

        IF(LOOK(i,l)) THEN
c	No overalps of Cas and Cbs
c	THE MOVE COULD BE SUCESSFULL. APPLY METROPOLIS CRITERION HERE

        ENEW=EHB(i,l)+ESHORT(i,l)+EHL(i,l)+ERESTR(i,l)+ECROSS(i,l)

        do kkkk=i,k
        ica(kkkk)=icao(kkkk)
        enddo
        do kkkk=j,k
        pp=kkkk-1
        ii=ica(pp)
        x(kkkk)=x(pp)+vx(ii)
        y(kkkk)=y(pp)+vy(ii)
        z(kkkk)=z(pp)+vz(ii)
        enddo

        EOLD=EHB(i,l)+ESHORT(i,l)+EHL(i,l)+ERESTR(i,l)+ECROSS(i,l)

        DE=ENEW-EOLD

        IF(DE.GT.0.0) THEN
        if(arand(ise).gt.EXP(-DE/ATEMP)) then
        go to 4001
        endif
        ENDIF

c	THE MOVE SUCESSFULL. WRITE-IN THE NEW CONFORMATION

        ENERG=ENERG+DE

        do kkkk=i,k
        ica(kkkk)=icab(kkkk)
        enddo
        do kkkk=j,k
        pp=kkkk-1
        ii=ica(pp)
        x(kkkk)=x(pp)+vx(ii)
        y(kkkk)=y(pp)+vy(ii)
        z(kkkk)=z(pp)+vz(ii)
        enddo

        go to 4001

        ELSE

c	Move rejected - restore the initial state of the chain

        do kkkk=i,k
        ica(kkkk)=icao(kkkk)
        enddo
        do kkkk=j,k
        pp=kkkk-1
        ii=ica(pp)
        x(kkkk)=x(pp)+vx(ii)
        y(kkkk)=y(pp)+vy(ii)
        z(kkkk)=z(pp)+vz(ii)
        enddo

        ENDIF

c	FOUR-to-22 BOOND MOVES 

4001    j=INT(arand(ise)*AL19)+3
        ppp=INT(arand(ise)*18.9999)+2
        i=j-1
        k=j+ppp
        L=k+1
        if(L.gt.lenf1) go to 7775

8809    IF(arand(ise).gt.0.25) THEN
        mx=INT(arand(ise)*2.999999)-1
        my=INT(arand(ise)*2.999999)-1
        mz=INT(arand(ise)*2.999999)-1
        if((mx*mx+my*my+mz*mz).eq.0) go to 8809
        ELSE
        mx=INT(arand(ise)*4.999999)-2
        my=INT(arand(ise)*4.999999)-2
        mz=INT(arand(ise)*4.999999)-2
        if((mx*mx+my*my+mz*mz).eq.0) go to 8809
        ENDIF

        ii=ica(i)
        wx=vx(ii)+mx
        wy=vy(ii)+my
        wz=vz(ii)+mz
        ir=wx*wx+wy*wy+wz*wz
        if(ir.lt.29) go to 8809
        if(ir.gt.49) go to 8809
        iii=vector(wx,wy,wz)
        if(.NOT.GOODC(ica(i-1),iii)) go to 7775
        if(.NOT.GOODC(iii,ica(j))) go to 7775

        kk=ica(k)
        wx=vx(kk)-mx
        wy=vy(kk)-my
        wz=vz(kk)-mz
        ir=wx*wx+wy*wy+wz*wz
        if(ir.lt.29) go to 7775
        if(ir.gt.49) go to 7775
        kkk=vector(wx,wy,wz)
        if(.NOT.GOODC(ica(k-1),kkk)) go to 7775
        if(.NOT.GOODC(kkk,ica(l))) go to 7775

        ICA(i)=iii
        ICA(k)=kkk
        do kkkk=j,k
        x(kkkk)=x(kkkk)+mx
        y(kkkk)=y(kkkk)+my
        z(kkkk)=z(kkkk)+mz
        enddo

        IF(LOOK(i,l)) THEN
c	No overalps of Cas and Cbs
c	THE MOVE COULD BE SUCESSFULL. APPLY METROPOLIS CRITERION HERE

        ENEW=EHB(i,l)+ESHORT(i,l)+EHL(i,l)+ERESTR(i,l)+ECROSS(i,l)

        ICA(i)=ii
        ICA(k)=kk
        do kkkk=j,k
        x(kkkk)=x(kkkk)-mx
        y(kkkk)=y(kkkk)-my
        z(kkkk)=z(kkkk)-mz
        enddo

        EOLD=EHB(i,l)+ESHORT(i,l)+EHL(i,l)+ERESTR(i,l)+ECROSS(i,l)

        DE=ENEW-EOLD

        IF(DE.GT.0.0) THEN
        if(arand(ise).gt.EXP(-DE/ATEMP)) then
        go to 7775
        endif
        ENDIF

c	THE MOVE SUCESSFULL. WRITE-IN THE NEW CONFORMATION

        ENERG=ENERG+DE

        do kkkk=j,k
        x(kkkk)=x(kkkk)+mx
        y(kkkk)=y(kkkk)+my
        z(kkkk)=z(kkkk)+mz
        enddo
        ICA(i)=iii
        ICA(k)=kkk

        go to 7775

        ELSE

c	Move rejected - restore the initial state of the chain

        ICA(i)=ii
        ICA(k)=kk
        do kkkk=j,k
        x(kkkk)=x(kkkk)-mx
        y(kkkk)=y(kkkk)-my
        z(kkkk)=z(kkkk)-mz
        enddo

        ENDIF

7775    CONTINUE

c	N TERMINUS

        JV3=ICA(3)
60      NV2=INT(arand(ise)*799.99 )+1
        NV1=INT(arand(ise)*799.99 )+1

        if(.NOT.GOODC(nv1,nv2)) go to 60
        if(.NOT.GOODC(nv2,jv3)) go to 60

        kx=x(1)
        ky=y(1)
        kz=z(1)
        px=x(2)
        py=y(2)
        pz=z(2)
        ica1=ica(1)
        ica2=ica(2)

        X(2)=X(3)-vx(nv2)
        Y(2)=Y(3)-vy(nv2)
        Z(2)=Z(3)-vz(nv2)

        X(1)=X(2)-vx(nv1)
        Y(1)=Y(2)-vy(nv1)
        Z(1)=Z(2)-vz(nv1)
        ICA(1)=nv1
        ICA(2)=nv2

        if(LOOK(2,3)) THEN

c	THE MOVE COULD BE SUCESSFULL. APPLY METROPOLIS CRITERION HERE

        ENEW=EHB(2,3)+ESHORT(2,3)+EHL(2,3)+ERESTR(2,3)+ECROSS(2,3)
        x(1)=kx
        y(1)=ky
        z(1)=kz
        x(2)=px
        y(2)=py
        z(2)=pz
        ICA(1)=ica1
        ICA(2)=ica2

        EOLD=EHB(2,3)+ESHORT(2,3)+EHL(2,3)+ERESTR(2,3)+ECROSS(2,3)

        DE=ENEW-EOLD

        IF(DE.GT.0.0) THEN
        if(arand(ise).gt.EXP(-DE/ATEMP)) then
        go to 81
        endif
        ENDIF

c	THE MOVE SUCESSFULL. WRITE-IN THE NEW CONFORMATION

        ENERG=ENERG+DE
        X(2)=X(3)-vx(nv2)
        Y(2)=Y(3)-vy(nv2)
        Z(2)=Z(3)-vz(nv2)

        X(1)=X(2)-vx(nv1)
        Y(1)=Y(2)-vy(nv1)
        Z(1)=Z(2)-vz(nv1)
        ICA(1)=nv1
        ICA(2)=nv2

        go to 81
        ELSE

c	Move rejected - restore the initial state of the chain
        x(1)=kx
        y(1)=ky
        z(1)=kz
        x(2)=px
        y(2)=py
        z(2)=pz
        ICA(1)=ica1
        ICA(2)=ica2
        ENDIF

c	C TERMINUS

81      JV3=ICA(LENF3)
80      NV2=INT(arand(ise)*799.99 )+1
        NV1=INT(arand(ise)*799.99 )+1
        if(.NOT.GOODC(nv2,nv1)) go to 80
        if(.NOT.GOODC(jv3,nv2)) go to 80

        kx=x(lenf)
        ky=y(lenf)
        kz=z(lenf)
        px=x(lenf1)
        py=y(lenf1)
        pz=z(lenf1)
        ica1=ica(lenf1)
        ica2=ica(lenf2)

        X(lenf1)=X(lenf2)+vx(nv2)
        Y(lenf1)=Y(lenf2)+vy(nv2)
        Z(lenf1)=Z(lenf2)+vz(nv2)

        X(lenf)=X(lenf1)+vx(nv1)
        Y(lenf)=Y(lenf1)+vy(nv1)
        Z(lenf)=Z(lenf1)+vz(nv1)
        
        ICA(lenf1)=nv1
        ICA(lenf2)=nv2

        if(LOOK(lenf2,lenf1)) THEN

c	THE MOVE COULD BE SUCESSFULL. APPLY METROPOLIS CRITERION HERE

        ENEW=EHB(lenf2,lenf1)+ESHORT(lenf2,lenf1)+EHL(lenf2,lenf1)
     &  +ERESTR(lenf2,lenf1)+ECROSS(lenf2,lenf1)

        x(lenf)=kx
        y(lenf)=ky
        z(lenf)=kz
        x(lenf1)=px
        y(lenf1)=py
        z(lenf1)=pz
        ICA(lenf1)=ica1
        ICA(lenf2)=ica2

        EOLD=EHB(lenf2,lenf1)+ESHORT(lenf2,lenf1)+EHL(lenf2,lenf1)
     &  +ERESTR(lenf2,lenf1)+ECROSS(lenf2,lenf1)

        DE=ENEW-EOLD

        IF(DE.GT.0.0) THEN
        if(arand(ise).gt.EXP(-DE/ATEMP)) then
        go to 9998
        endif
        ENDIF

c	THE MOVE SUCESSFULL. WRITE-IN THE NEW CONFORMATION

        ENERG=ENERG+DE
        X(lenf1)=X(lenf2)+vx(nv2)
        Y(lenf1)=Y(lenf2)+vy(nv2)
        Z(lenf1)=Z(lenf2)+vz(nv2)

        X(lenf)=X(lenf1)+vx(nv1)
        Y(lenf)=Y(lenf1)+vy(nv1)
        Z(lenf)=Z(lenf1)+vz(nv1)

        ICA(lenf1)=nv1
        ICA(lenf2)=nv2

        go to 9998

        ELSE

c	Move rejected - restore the initial state of the chain		
        x(lenf)=kx
        y(lenf)=ky
        z(lenf)=kz
        x(lenf1)=px
        y(lenf1)=py
        z(lenf1)=pz
        ICA(lenf1)=ica1
        ICA(lenf2)=ica2

        ENDIF

9998    do i=1,lenf
        xreps(i,itemp,imol)=x(i)
        yreps(i,itemp,imol)=y(i)
        zreps(i,itemp,imol)=z(i)
        enddo

        if(iphot.eq.phot) then
        eitemp=EHL(2,lenf1)+ECROSS(2,lenf1)
        EREPLICA(imol,itemp)=energ-eitemp
c        write(*,*) EREPLICA(imol,itemp),EHB(2,lenf1)+ESHORT(2,lenf1)
c     &  +ERESTR(2,lenf1)
        endif

9999    CONTINUE
7777    CONTINUE
 
        EREPLICAS(itemp)=0.0

        do imol=1,mols

        lenf=length(imol)
        lenf1=lenf-1

        DO i=1,lenf
        x(i)=xreps(i,itemp,imol)
        y(i)=yreps(i,itemp,imol)
        z(i)=zreps(i,itemp,imol)
        enddo

        DO I=1,LENF1
        J=I+1
        WX=X(J)-X(I)
        WY=Y(J)-Y(I)
        WZ=Z(J)-Z(I)
        ICA(I)=VECTOR(WX,WY,WZ)
        ENDDO
        ICA(LENF)=1

        EINTER(imol,itemp)=EHL(2,lenf1)+ECROSS(2,lenf1)
        EREPLICAS(itemp)=EREPLICAS(itemp)+EREPLICA(imol,itemp)
     &  +0.5*EINTER(imol,itemp)

        enddo

c	end of ITEMPs replica	
 
        ENDDO

c	Search for the lowest energy replica
        energ=EREPLICAS(1)
        ITEMP=1
        do i=1,REPLICAS
        if(ereplicas(i).lt.energ) then
        energ=ereplicas(i)
        ITEMP=i
        endif
        enddo

        NFRAME=NFRAME+1

        itemp1=itemp

        do itemp=1,replicas
        do imol=1,mols
        WRITE(9,716) nframe,length(imol),ereplica(imol,itemp),
     &  einter(imol,itemp),ereplicas(itemp)
     &  ,atrep(itemp),repnum(itemp)
        write(9,713) (xreps(i,itemp,imol),yreps(i,itemp,imol),
     &  zreps(i,itemp,imol),i=1,length(imol))
        enddo
        enddo

        itemp=itemp1

716     format(2I6,3F10.2,F8.3,I6)
713     format(12I5)

7000    CONTINUE
 
        ENDDO

c	===============================================================	
 
c	CHECK THE GEOMETRY and WRITE THE FINAL CHAIN

        kcc=0
        kcb=0
        kbb=0

        OPEN(UNIT=19,FILE='ACHAINS_NEW')

        do imol=1,mols

        lenf=length(imol)
        lenf1=lenf-1
        lenf2=lenf-2
        lenf3=lenf-3

        DO i=1,lenf
        x(i)=xreps(i,itemp,imol)
        y(i)=yreps(i,itemp,imol)
        z(i)=zreps(i,itemp,imol)
        enddo

        DO I=1,LENF1
        J=I+1
        WX=X(J)-X(I)
        WY=Y(J)-Y(I)
        WZ=Z(J)-Z(I)
        ICA(I)=VECTOR(WX,WY,WZ)
        ENDDO
        ICA(LENF)=1

        DO J=1,LENF
        I=J-1
        if(j.ne.1.and.j.ne.lenf) THEN
        II=ICA(I)
        JJ=ICA(J)
        IF(.NOT.GOODC(II,JJ)) THEN
        WRITE(6,8112)I,J,vx(ii),vy(ii),vz(ii),vx(jj),vy(jj),vz(jj)
8112    FORMAT(5X,'WARNING2 -WRONG INPUT CHAIN - VECTORS ',8I4)
        STOP
        ENDIF
        ENDIF
        ENDDO

c	PUT ALSO TEST FOR OVERLAPS OF Ca with Cb - initial conformations
c	from program NatCa.f may sometimes have overlaps - they relax
c	rapidly

        do i=2,lenf3
        ix=x(i)
        iy=y(i)
        iz=z(i)
        bx=ix+CBX(ica(i-1),ica(i))
        by=iy+CBY(ica(i-1),ica(i))
        bz=iz+CBZ(ica(i-1),ica(i))
        do j=i+2,lenf1
        jx=x(j)
        jy=y(j)
        jz=z(j)
        cx=jx+CBX(ica(j-1),ica(j))
        cy=jy+CBY(ica(j-1),ica(j))
        cz=jz+CBZ(ica(j-1),ica(j))
        ir=(ix-jx)**2+(iy-jy)**2+(iz-jz)**2
        if(ir.lt.36) kcc=kcc+1
        if(SEQ(imol,i).NE.0.AND.SEQ(imol,j).NE.0) then
        ar=(bx-cx)**2+(by-cy)**2+(bz-cz)**2
        if(ar.lt.36.0) kbb=kbb+1
        endif
        if(SEQ(imol,j).NE.0) THEN
        ar=(ix-cx)**2+(iy-cy)**2+(iz-cz)**2
        if(ar.lt.43.0) kcb=kcb+1
        ENDIF
        if(SEQ(imol,i).NE.0) THEN
        ar=(bx-jx)**2+(by-jy)**2+(bz-jz)**2
        if(ar.lt.43.0) kcb=kcb+1
        endif
        enddo
        enddo

        do kmol=imol+1,mols
        do i=2,lenf1
        ix=x(i)
        iy=y(i)
        iz=z(i)
        bx=ix+CBX(ica(i-1),ica(i))
        by=iy+CBY(ica(i-1),ica(i))
        bz=iz+CBZ(ica(i-1),ica(i))
        do j=2,length(kmol)-1
        jx=xreps(j,itemp,kmol)
        jy=yreps(j,itemp,kmol)
        jz=zreps(j,itemp,kmol)
        wx=xreps(j+1,itemp,kmol)-jx
        wy=yreps(j+1,itemp,kmol)-jy
        wz=zreps(j+1,itemp,kmol)-jz
        ica1=vector(wx,wy,wz)
        wx=jx-xreps(j-1,itemp,kmol)
        wy=jy-yreps(j-1,itemp,kmol)
        wz=jz-zreps(j-1,itemp,kmol)
        ica2=vector(wx,wy,wz)
        cx=jx+CBX(ica2,ica1)
        cy=jy+CBY(ica2,ica1)
        cz=jz+CBZ(ica2,ica1)
        ir=(ix-jx)**2+(iy-jy)**2+(iz-jz)**2
        if(ir.lt.36) kcc=kcc+1
        if(SEQ(imol,i).NE.0  .AND. SEQ(kmol,j).NE.0) then
        ar=(bx-cx)**2+(by-cy)**2+(bz-cz)**2
        if(ar.lt.36.0) kbb=kbb+1
        endif
        if(SEQ(kmol,j).NE.0) THEN
        ar=(ix-cx)**2+(iy-cy)**2+(iz-cz)**2
        if(ar.lt.43.0) kcb=kcb+1
        ENDIF
        if(SEQ(imol,i).NE.0) THEN
        ar=(bx-jx)**2+(by-jy)**2+(bz-jz)**2
        if(ar.lt.43.0) kcb=kcb+1
        endif
        enddo
        enddo
        enddo
        write(6,*) "CHAIN: ",imol
        write(6,819) EREPLICA(imol,itemp),EINTER(imol,itemp)

c	Compute from scratch the final energy for the consistency
c	comparison with the cumulated value of energy - TEST

        energ=EHB(2,lenf1)+ESHORT(2,lenf1)+ERESTR(2,lenf1)
        energi=EHL(2,lenf1)+ECROSS(2,lenf1)

        write(6,819) energ,energi
819     FORMAT('FINAL ENERGY', 2f10.2)

        write(6,*)
        write(6,*) 'Ca-Ca overlaps  ', kcc
        write(6,*) 'Ca-Cb overlaps  ', kcb
        write(6,*) 'Cb-Cb overlaps  ', kbb
        write(6,*)

        do k=1,REPLICAS
        write(19,*) LENF
        DO I=1,LENF
        WRITE(19,8000) XREPS(I,k,imol),YREPS(I,k,imol),ZREPS(I,k,imol)
        ENDDO
        enddo

8000    format(3i5)

        enddo

        WRITE(6,*)
        WRITE(6,*) '   	FINAL REPLICA DISTRIBUTION'
        do i=1, REPLICAS
        write(6,5005) i,ATREP(i),EREPLICAS(i),SWAPS(i)
        enddo
5005    format('   replica #',I4,'   T =', f6.2,'  E =',f8.1, I10)

        STOP
        END

c	END OF THE MAIN PROGRAM

c	***************************************************************		

        FUNCTION LOOK(jj,kk)

        IMPLICIT INTEGER(I-Z)
        LOGICAL LOOK,LOOKL

        PARAMETER(NDIM=500)
        PARAMETER(NMOLS=2)

        COMMON /CHAINS/  ICA(NDIM),X(NDIM),Y(NDIM),Z(NDIM)
        COMMON /LENGHTS/ PROD(800,800)
        COMMON /SEQE/    SEQ(nmols,ndim),SEC(nmols,ndim)
        COMMON /BETAS/   CBX(800,800),CBY(800,800),CBZ(800,800)
        common /lenfs/ mols,imol,itemp,lenf,lenf1,lenf2,lenf3,
     &  al2,al4,al5,al19,length(nmols)

c	CHECKS Ca-Ca, Ca-Cb and Cb-Cb distances and overlaps

        DO k=jj,kk

        px=x(k)
        py=y(k)
        pz=z(k)
        if(SEQ(imol,k).NE.0) THEN
        nv1=ica(k-1)
        nv2=ica(k)
        bx=CBX(nv1,nv2)+px
        by=CBY(nv1,nv2)+py
        bz=CBZ(nv1,nv2)+pz
        ENDIF

        istart=k-2
        if(istart.gt.1) THEN
        do i=2,istart
        irx=(px-x(i))**2
        if(irx.lt.120) then
        iry=(py-y(i))**2
        if(iry.lt.120) then
        ir=irx+iry+(pz-z(i))**2
        if(ir.lt.120) then

c	detect detailed overlaps Ca-Ca, Cb-Cb, Ca-Cb, Cb-Ca

        if(ir.lt.36) then
        LOOK=.FALSE.
        RETURN
        endif

        if(SEQ(imol,k).NE.0) THEN
        ar=(bx-x(i))**2+(by-y(i))**2+(bz-z(i))**2
        if(ar.lt.43.0) then
        LOOK=.FALSE.
        RETURN
        endif
        ENDIF

        IF(SEQ(imol,i).NE.0) THEN
        icai=ICA(i-1)
        icaj=ICA(i)
        ax= CBX(icai,icaj)+x(i)
        ay= CBY(icai,icaj)+y(i)
        az= CBZ(icai,icaj)+z(i)

        IF(SEQ(imol,k).NE.0) THEN
        ar=(ax-bx)**2+(ay-by)**2+(az-bz)**2
        if(ar.lt.36.0) then
        LOOK=.FALSE.
        RETURN
        endif
        ENDIF

        ar=(ax-px)**2+(ay-py)**2+(az-pz)**2
        if(ar.lt.43.0) then
        LOOK=.FALSE.
        RETURN
        endif
        ENDIF

        endif
        endif
        endif
        enddo
        ENDIF

        iend=max(k+2,kk+1)
        if(iend.lt.lenf) THEN
        do i=iend,lenf1
        irx=(px-x(i))**2
        if(irx.lt.120) then
        iry=(py-y(i))**2
        if(iry.lt.120) then
        ir=irx+iry+(pz-z(i))**2
        if(ir.lt.120) then

c	detect detailed overlaps Ca-Ca, Cb-Cb, Ca-Cb, Cb-Ca

        if(ir.lt.36) then
        LOOK=.FALSE.
        RETURN
        endif

        if(SEQ(imol,k).NE.0) THEN
        ar=(bx-x(i))**2+(by-y(i))**2+(bz-z(i))**2
        if(ar.lt.43.0) then
        LOOK=.FALSE.
        RETURN
        endif
        ENDIF

        if(SEQ(imol,i).NE.0) THEN
        icai=ICA(i-1)
        icaj=ICA(i)
        ax= CBX(icai,icaj)+x(i)
        ay= CBY(icai,icaj)+y(i)
        az= CBZ(icai,icaj)+z(i)

        IF(SEQ(imol,k).NE.0) THEN
        ar=(ax-bx)**2+(ay-by)**2+(az-bz)**2
        if(ar.lt.36.0) then
        LOOK=.FALSE.
        RETURN
        endif
        ENDIF

        ar=(ax-px)**2+(ay-py)**2+(az-pz)**2
        if(ar.lt.43.0) then
        LOOK=.FALSE.
        RETURN
        endif

        endif
        endif
        ENDIF

        endif
        enddo
        ENDIF

        ENDDO

c	NO HARD CORE OVERLAPPS DETECTED

        LOOK=LOOKL(jj,kk)
        
        RETURN
        END
c	***************************************************************

        FUNCTION LOOKL(jj,kk)
        IMPLICIT INTEGER(I-Z)
        LOGICAL LOOKL

        PARAMETER(NDIM=500)
        PARAMETER(NREPS=10)
        PARAMETER(NMOLS=2)

        COMMON /CHAINS/  ICA(NDIM),X(NDIM),Y(NDIM),Z(NDIM)
        COMMON /LENGHTS/ PROD(800,800)
        COMMON /SEQE/    SEQ(nmols,ndim),SEC(nmols,ndim)
        COMMON /BETAS/   CBX(800,800),CBY(800,800),CBZ(800,800)
        common /reps/   xreps(ndim,NREPS,nmols),
     &                  yreps(ndim,NREPS,nmols),
     &                  zreps(ndim,NREPS,nmols)
        common /lenfs/ mols,imol,itemp,lenf,lenf1,lenf2,lenf3,
     &  al2,al4,al5,al19,length(nmols)
        common /vec/ vector(-7:7,-7:7,-7:7)

c	CHECKS Ca-Ca, Ca-Cb and Cb-Cb distances and overlaps

        DO k=jj,kk

        px=x(k)
        py=y(k)
        pz=z(k)
        if(SEQ(imol,k).NE.0) THEN
        nv1=ica(k-1)
        nv2=ica(k)
        bx=CBX(nv1,nv2)+px
        by=CBY(nv1,nv2)+py
        bz=CBZ(nv1,nv2)+pz
        ENDIF

        do 1983 kmol=1,mols

        if(kmol.eq.imol) goto 1983

        do i=2,length(kmol)-1
        irx=(px-xreps(i,itemp,kmol))**2
        if(irx.lt.120) then
        iry=(py-yreps(i,itemp,kmol))**2
        if(iry.lt.120) then
        ir=irx+iry+(pz-zreps(i,itemp,kmol))**2
        if(ir.lt.120) then
        
c	detect detailed overlaps Ca-Ca, Cb-Cb, Ca-Cb, Cb-Ca

        if(ir.lt.36) then
        LOOKL=.FALSE.
        RETURN
        endif

        if(SEQ(imol,k).NE.0) THEN
        ar=(bx-xreps(i,itemp,kmol))**2+(by-yreps(i,itemp,kmol))**2
     &  +(bz-zreps(i,itemp,kmol))**2
        if(ar.lt.43.0) then
        LOOKL=.FALSE.
        RETURN
        endif
        ENDIF

        IF(SEQ(kmol,i).NE.0) THEN
        wx=xreps(i,itemp,kmol)-xreps(i-1,itemp,kmol)
        wy=yreps(i,itemp,kmol)-yreps(i-1,itemp,kmol)
        wz=zreps(i,itemp,kmol)-zreps(i-1,itemp,kmol)
        icai=vector(wx,wy,wz)
        wx=xreps(i+1,itemp,kmol)-xreps(i,itemp,kmol)
        wy=yreps(i+1,itemp,kmol)-yreps(i,itemp,kmol)
        wz=zreps(i+1,itemp,kmol)-zreps(i,itemp,kmol)
        icaj=vector(wx,wy,wz)
        ax=CBX(icai,icaj)+xreps(i,itemp,kmol)
        ay=CBY(icai,icaj)+yreps(i,itemp,kmol)
        az=CBZ(icai,icaj)+zreps(i,itemp,kmol)

        IF(SEQ(imol,k).NE.0) THEN
        ar=(ax-bx)**2+(ay-by)**2+(az-bz)**2
        if(ar.lt.36.0) then
        LOOKL=.FALSE.
        RETURN
        endif
        ENDIF

        ar=(ax-px)**2+(ay-py)**2+(az-pz)**2
        if(ar.lt.43.0) then
        LOOKL=.FALSE.
        RETURN
        endif
        ENDIF

        endif
        endif
        endif
        enddo

1983    continue

        ENDDO

c	NO HARD CORE OVERLAPPS DETECTED

        LOOKL=.TRUE.

        RETURN
        END

c	***************************************************************	
c	A random number generator by AK - tested
        FUNCTION arand(ise)
        common/arandom/a,b,c,d
        g=(a+b+c)*d
        g=g-aint(g)
        a=b
        b=c
        c=g
        arand=g
        return
        end

C	***************************************************************	

        FUNCTION EHB(jjjj,kkkk)

        IMPLICIT INTEGER(I-Z)

        PARAMETER(NDIM=500)
        PARAMETER(NMOLS=2)

        COMMON /CHAINS/  ICA(NDIM),X(NDIM),Y(NDIM),Z(NDIM)
        COMMON /LENGHTS/ PROD(800,800)
        COMMON /SEQE/    SEQ(nmols,ndim),SEC(nmols,ndim)
        COMMON /BETAS/   CBX(800,800),CBY(800,800),CBZ(800,800)
        COMMON /BISEC/   CAX(800,800),CAY(800,800),CAZ(800,800)
        COMMON /HB/      HBX(800,800),HBY(800,800),HBZ(800,800)

        COMMON /ENERGY/  EHBOND,ESC,EREP,BAT,
     &  EHBIJ(nmols,nmols,ndim,ndim)
        COMMON /pair/   apa(nmols,nmols,ndim,ndim,2,2),
     &                  app(nmols,nmols,ndim,ndim,2,2),
     &                  apm(nmols,nmols,ndim,ndim,2,2)
        COMMON /size/   arla(2,0:19,0:19,2,2),arlm(2,0:19,0:19,2,2),
     &                  arlp(2,0:19,0:19,2,2)
        COMMON /sizea/  ala(2,0:19,0:19,2,2),alm(2,0:19,0:19,2,2),
     &                  alp(2,0:19,0:19,2,2)
        
        COMMON /KBB/    kbin(800,800)
        COMMON /one/    acrit(nmols),eoinp(0:19,0:100)
        COMMON /SG/ GX(800,800,0:19),GY(800,800,0:19),GZ(800,800,0:19)
        common /lenfs/ mols,imol,itemp,lenf,lenf1,lenf2,lenf3,
     &  al2,al4,al5,al19,length(nmols)

        EHB=0.0

        DO k=jjjj,kkkk
        im=seq(imol,k)
        px=x(k)
        py=y(k)
        pz=z(k)

        nv1=ica(k-1)
        nv2=ica(k)

        bx=HBX(nv1,nv2)
        by=HBY(nv1,nv2)
        bz=HBZ(nv1,nv2)
        hx=CAX(nv1,nv2)
        hy=CAY(nv1,nv2)
        hz=CAZ(nv1,nv2)
        agx=px+GX(nv1,nv2,im)
        agy=py+GY(nv1,nv2,im)
        agz=pz+GZ(nv1,nv2,im)

        kb= kbin(nv1,nv2)

        istart=k-2
        if(istart.gt.1) THEN
        do 1001 i=2,istart

        irx=(px-x(i))**2
        if(irx.lt.400) then
        iry=(py-y(i))**2
        if(iry.lt.400) then
        ir=irx+iry+(pz-z(i))**2
        if(ir.lt.400) then

        xi=x(i)
        yi=y(i)
        zi=z(i)

        ica1=ica(i-1)
        ica2=ica(i)
        ib=kbin(ica1,ica2)
        in=seq(imol,i)

        bgx=xi+GX(ica1,ica2,in)
        bgy=yi+GY(ica1,ica2,in)
        bgz=zi+GZ(ica1,ica2,in)

c	Repulsion between the peptide bond and Ca

        IF(iabs(i-k).gt.3) THEN
        k1=k+1
        i1=i+1
        dx=(px+x(k1))/2.0-xi
        dy=(py+y(k1))/2.0-yi
        dz=(pz+z(k1))/2.0-zi
        arrr=dx*dx+dy*dy+dz*dz
        if(arrr.lt.54.0) EHB=EHB+EREP

        dx=px-(xi+x(i1))/2.0
        dy=py-(yi+y(i1))/2.0
        dz=pz-(zi+z(i1))/2.0
        arrr=dx*dx+dy*dy+dz*dz
        if(arrr.lt.54.0) EHB=EHB+EREP

        k1=k-1
        i1=i-1
        dx=(px+x(k1))/2.0-xi
        dy=(py+y(k1))/2.0-yi
        dz=(pz+z(k1))/2.0-zi
        arrr=dx*dx+dy*dy+dz*dz
        if(arrr.lt.54.0) EHB=EHB+EREP

        dx=px-(xi+x(i1))/2.0
        dy=py-(yi+y(i1))/2.0
        dz=pz-(zi+z(i1))/2.0
        arrr=dx*dx+dy*dy+dz*dz
        if(arrr.lt.54.0) EHB=EHB+EREP

c	side group-Ca distance
        if(im.gt.1) then
        aks=(agx-xi)**2+(agy-yi)**2+(agz-zi)**2
        if(aks.lt.43.0) EHB=EHB+EREP
        endif

c	side group-Ca distance
        if(in.gt.1) then
        aks=(bgx-px)**2+(bgy-py)**2+(bgz-pz)**2
        if(aks.lt.43.0) EHB=EHB+EREP
        endif
        ENDIF

c	side group- side group distance		
        ar=(agx-bgx)**2+(agy-bgy)**2+(agz-bgz)**2

        fx=CAX(ica1,ica2)
        fy=CAY(ica1,ica2)
        fz=CAZ(ica1,ica2)

c	contact orientation (parallel, othogonal, or antiparallel)		
        aor=hx*fx+hy*fy+hz*fz

        IF(aor.gt.0.5) THEN
        if(ar.gt.arlp(1,in,im,ib,kb)) THEN
        if(ar.lt.alp(1,in,im,ib,kb)) then 
        EHB=EHB+app(imol,imol,i,k,ib,kb) 
        endif
        else
        EHB=EHB+EREP
c	repulsons from side group overlapps weaker than for Ca or Cb
        endif
        ELSE
        IF(aor.lt.-0.5) THEN
        if(ar.gt.arla(1,in,im,ib,kb)) THEN
        if(ar.lt.ala(1,in,im,ib,kb)) then
        EHB=EHB+apa(imol,imol,i,k,ib,kb) 
        endif
        else
        EHB=EHB+EREP
        endif
        ELSE
        if(ar.gt.arlm(1,in,im,ib,kb)) THEN
        if(ar.lt.alm(1,in,im,ib,kb)) then
        EHB=EHB+apm(imol,imol,i,k,ib,kb) 
        endif
        else
        EHB=EHB+EREP
        endif
        ENDIF
        ENDIF

        IF(aor.gt.0.0) THEN
        idist=iabs(i-k)
        if(idist.gt.2) THEN
        if(ir.lt.100) THEN

        IF(idist.ne.4) THEN
        IF(idist.gt.3) then
        if(sec(imol,k).eq.2.OR.sec(imol,i).eq.2) go to 1001
        endif

        cx=HBX(ica1,ica2)
        cy=HBY(ica1,ica2)
        cz=HBZ(ica1,ica2)

c	COOPERATIVITY OF THE H-bonds, plates of contacting fragments
c		must be almost parallel

        if(prod(nv1,ica1).GT.0.AND.prod(nv2,ica2).GT.0) THEN

c	chains must be parallel or (****)

        IF(idist.eq.3) then
        if(sec(imol,k-1).eq.4.OR.sec(imol,i+1).eq.4) go to 1001
        endif
        if(idist.gt.3.AND.idist.lt.20) go to 1001
c	distance dependent attraction of Ca in the H-bond type
c	geometry
        id=1
        ELSE
        if(prod(nv1,ica2).gt.0.OR.prod(nv2,ica1).gt.0) go to 1001
        id=-1
        ENDIF

        ax=(px-xi)
        ay=(py-yi)
        az=(pz-zi)

        delta=0.0
        k1=k+1
        i1=i+id
        dx=(px+x(k1)-xi-x(i1))/2.0
        dy=(py+y(k1)-yi-y(i1))/2.0
        dz=(pz+z(k1)-zi-z(i1))/2.0
        arrr=dx*dx+dy*dy+dz*dz
        if(arrr.lt.80.0) delta=delta+0.5

        k1=k-1
        i1=i-id
        dx=(px+x(k1)-xi-x(i1))/2.0
        dy=(py+y(k1)-yi-y(i1))/2.0
        dz=(pz+z(k1)-zi-z(i1))/2.0
        arrr=dx*dx+dy*dy+dz*dz
        if(arrr.lt.80.0) delta=delta+0.5

        ar=(bx-ax)**2+(by-ay)**2+(bz-az)**2
        if(ar.lt.7.0) then
        EHB=EHB+EHBIJ(imol,imol,i,k)*(0.5+delta)
        else
        ar=(bx+ax)**2+(by+ay)**2+(bz+az)**2
        if(ar.lt.7.0) then
        EHB=EHB+EHBIJ(imol,imol,i,k)*(0.5+delta)
        endif
        endif

        ar=(cx-ax)**2+(cy-ay)**2+(cz-az)**2
        if(ar.lt.7.0) then
        EHB=EHB+EHBIJ(imol,imol,i,k)*(0.5+delta)
        else
        ar=(cx+ax)**2+(cy+ay)**2+(cz+az)**2
        if(ar.lt.7.0) then
        EHB=EHB+EHBIJ(imol,imol,i,k)*(0.5+delta) 
        endif
        endif

        ENDIF
        ENDIF
        ENDIF
        ENDIF
        endif
        endif
        endif
1001    CONTINUE
        ENDIF

        iend=max(k+2,kkkk+1)
        if(iend.lt.lenf) THEN
        do 1002 i=iend,lenf1

        irx=(px-x(i))**2
        if(irx.lt.400) then
        iry=(py-y(i))**2
        if(iry.lt.400) then
        ir=irx+iry+(pz-z(i))**2
        if(ir.lt.400) then

        xi=x(i)
        yi=y(i)
        zi=z(i)
        ica1=ica(i-1)
        ica2=ica(i)
        ib=kbin(ica1,ica2)
        in=seq(imol,i)
        bgx=xi+GX(ica1,ica2,in)
        bgy=yi+GY(ica1,ica2,in)
        bgz=zi+GZ(ica1,ica2,in)

c	Repulsion between the peptide bond and Ca

        IF(iabs(i-k).gt.3) THEN
        k1=k+1
        i1=i+1
        dx=(px+x(k1))/2.0-xi
        dy=(py+y(k1))/2.0-yi
        dz=(pz+z(k1))/2.0-zi
        arrr=dx*dx+dy*dy+dz*dz
        if(arrr.lt.54.0) EHB=EHB+EREP

        dx=px-(xi+x(i1))/2.0
        dy=py-(yi+y(i1))/2.0
        dz=pz-(zi+z(i1))/2.0
        arrr=dx*dx+dy*dy+dz*dz
        if(arrr.lt.54.0) EHB=EHB+EREP

        k1=k-1
        i1=i-1
        dx=(px+x(k1))/2.0-xi
        dy=(py+y(k1))/2.0-yi
        dz=(pz+z(k1))/2.0-zi
        arrr=dx*dx+dy*dy+dz*dz
        if(arrr.lt.54.0) EHB=EHB+EREP

        dx=px-(xi+x(i1))/2.0
        dy=py-(yi+y(i1))/2.0
        dz=pz-(zi+z(i1))/2.0
        arrr=dx*dx+dy*dy+dz*dz
        if(arrr.lt.54.0) EHB=EHB+EREP

        if(im.gt.1) then
        aks=(agx-xi)**2+(agy-yi)**2+(agz-zi)**2
        if(aks.lt.43.0) EHB=EHB+EREP
        endif

        if(in.gt.1) then
        aks=(bgx-px)**2+(bgy-py)**2+(bgz-pz)**2
        if(aks.lt.43.0) EHB=EHB+EREP
        endif
        ENDIF

        ar=(agx-bgx)**2+(agy-bgy)**2+(agz-bgz)**2

        fx=CAX(ica1,ica2)
        fy=CAY(ica1,ica2)
        fz=CAZ(ica1,ica2)

        aor=hx*fx+hy*fy+hz*fz

        IF(aor.gt.0.5) THEN
        if(ar.gt.arlp(1,im,in,kb,ib)) THEN
        if(ar.lt.alp(1,im,in,kb,ib)) then
        EHB=EHB+app(imol,imol,k,i,kb,ib) 
        endif

        else
        EHB=EHB+EREP
        endif
        ELSE
        IF(aor.lt.-0.5) THEN
        if(ar.gt.arla(1,im,in,kb,ib)) THEN
        if(ar.lt.ala(1,im,in,kb,ib)) then
        EHB=EHB+apa(imol,imol,k,i,kb,ib) 
        endif
        else
        EHB=EHB+EREP
        endif
        ELSE
        if(ar.gt.arlm(1,im,in,kb,ib)) THEN
        if(ar.lt.alm(1,im,in,kb,ib)) then
        EHB=EHB+apm(imol,imol,k,i,kb,ib) 
        endif
        else
        EHB=EHB+EREP
        endif
        ENDIF
        ENDIF

        IF(aor.gt.0.0) THEN
        idist=iabs(i-k)
        if(idist.gt.2) THEN
        if(ir.lt.100) THEN

        IF(idist.ne.4) THEN
        IF(idist.gt.3) then
        if(sec(imol,k).eq.2.OR.sec(imol,i).eq.2) go to 1002
        endif

        cx=HBX(ica1,ica2)
        cy=HBY(ica1,ica2)
        cz=HBZ(ica1,ica2)

c		COOPERATIVITY

        if(prod(nv1,ica1).GT.0.AND.prod(nv2,ica2).GT.0) THEN
        IF(idist.eq.3) then
        if(sec(imol,k+1).eq.4.OR.sec(imol,i-1).eq.4) go to 1002
        endif
        if(idist.gt.3.AND.idist.lt.20) go to 1002
        id=1
        ELSE
        if(prod(nv1,ica2).gt.0.OR.prod(nv2,ica1).gt.0) go to 1002
        id=-1
        ENDIF

        ax=(px-xi)
        ay=(py-yi)
        az=(pz-zi)

        delta=0.0
        k1=k+1
        i1=i+id
        dx=(px+x(k1)-xi-x(i1))/2.0
        dy=(py+y(k1)-yi-y(i1))/2.0
        dz=(pz+z(k1)-zi-z(i1))/2.0
        arrr=dx*dx+dy*dy+dz*dz
        if(arrr.lt.80.0) delta=delta+0.5

        k1=k-1
        i1=i-id
        dx=(px+x(k1)-xi-x(i1))/2.0
        dy=(py+y(k1)-yi-y(i1))/2.0
        dz=(pz+z(k1)-zi-z(i1))/2.0
        arrr=dx*dx+dy*dy+dz*dz
        if(arrr.lt.80.0) delta=delta+0.5

        ar=(bx-ax)**2+(by-ay)**2+(bz-az)**2
        if(ar.lt.7.0) then
        EHB=EHB+EHBIJ(imol,imol,i,k)*(0.5+delta)
        else
        ar=(bx+ax)**2+(by+ay)**2+(bz+az)**2
        if(ar.lt.7.0) then
        EHB=EHB+EHBIJ(imol,imol,i,k)*(0.5+delta)
        endif
        endif

        ar=(cx-ax)**2+(cy-ay)**2+(cz-az)**2
        if(ar.lt.7.0) then
        EHB=EHB+EHBIJ(imol,imol,i,k)*(0.5+delta)
        else
        ar=(cx+ax)**2+(cy+ay)**2+(cz+az)**2
        if(ar.lt.7.0) then
        EHB=EHB+EHBIJ(imol,imol,i,k)*(0.5+delta) 
        endif
        endif

        ENDIF
        ENDIF
        ENDIF
        ENDIF

        endif
        endif
        endif
1002    CONTINUE
        ENDIF

        ENDDO

        RETURN
        END

C       ***************************************************************

        FUNCTION EHL(jjjj,kkkk)

        IMPLICIT INTEGER(I-Z)

        PARAMETER(NDIM=500)
        PARAMETER(NREPS=10)
        PARAMETER(NMOLS=2)

        COMMON /CHAINS/  ICA(NDIM),X(NDIM),Y(NDIM),Z(NDIM)
        COMMON /LENGHTS/ PROD(800,800)
        COMMON /SEQE/    SEQ(nmols,ndim),SEC(nmols,ndim)
        COMMON /BETAS/   CBX(800,800),CBY(800,800),CBZ(800,800)
        COMMON /BISEC/   CAX(800,800),CAY(800,800),CAZ(800,800)
        COMMON /HB/      HBX(800,800),HBY(800,800),HBZ(800,800)

        COMMON /ENERGY/  EHBOND,ESC,EREP,BAT,
     &  EHBIJ(nmols,nmols,ndim,ndim)
        COMMON /pair/   apa(nmols,nmols,ndim,ndim,2,2),
     &                  app(nmols,nmols,ndim,ndim,2,2),
     &                  apm(nmols,nmols,ndim,ndim,2,2)
        COMMON /size/   arla(2,0:19,0:19,2,2),arlm(2,0:19,0:19,2,2),
     &                  arlp(2,0:19,0:19,2,2)
        COMMON /sizea/  ala(2,0:19,0:19,2,2),alm(2,0:19,0:19,2,2),
     &                  alp(2,0:19,0:19,2,2)
        
        COMMON /KBB/    kbin(800,800)
        COMMON /one/    acrit(nmols),eoinp(0:19,0:100)
        common /vec/ vector(-7:7,-7:7,-7:7)
        COMMON /SG/ GX(800,800,0:19),GY(800,800,0:19),GZ(800,800,0:19)
        common /lenfs/ mols,imol,itemp,lenf,lenf1,lenf2,lenf3,
     &  al2,al4,al5,al19,length(nmols)
        common /reps/   xreps(ndim,NREPS,nmols),
     &                  yreps(ndim,NREPS,nmols),
     &                  zreps(ndim,NREPS,nmols)

        EHL=0.0

        DO k=jjjj,kkkk
        im=seq(imol,k)
        px=x(k)
        py=y(k)
        pz=z(k)

        nv1=ica(k-1)
        nv2=ica(k)

        bx=HBX(nv1,nv2)
        by=HBY(nv1,nv2)
        bz=HBZ(nv1,nv2)
        hx=CAX(nv1,nv2)
        hy=CAY(nv1,nv2)
        hz=CAZ(nv1,nv2)
        agx=px+GX(nv1,nv2,im)
        agy=py+GY(nv1,nv2,im)
        agz=pz+GZ(nv1,nv2,im)

        kb= kbin(nv1,nv2)

        do 1004 kmol=1,mols

        if(kmol.eq.imol) goto 1004

        do 1003 i=2,length(kmol)-1

        xi=xreps(i,itemp,kmol)
        yi=yreps(i,itemp,kmol)
        zi=zreps(i,itemp,kmol)
        xi1=xreps(i+1,itemp,kmol)
        yi1=yreps(i+1,itemp,kmol)
        zi1=zreps(i+1,itemp,kmol)
        xim1=xreps(i-1,itemp,kmol)
        yim1=yreps(i-1,itemp,kmol)
        zim1=zreps(i-1,itemp,kmol)

        irx=(px-xi)**2
        if(irx.lt.400) then
        iry=(py-yi)**2
        if(iry.lt.400) then
        ir=irx+iry+(pz-zi)**2
        if(ir.lt.400) then

        wx=xi-xim1
        wy=yi-yim1
        wz=zi-zim1
        ica1=vector(wx,wy,wz)
        wx=xi1-xi
        wy=yi1-yi
        wz=zi1-zi
        ica2=vector(wx,wy,wz)

        ib=kbin(ica1,ica2)
        in=seq(kmol,i)

        bgx=xi+GX(ica1,ica2,in)
        bgy=yi+GY(ica1,ica2,in)
        bgz=zi+GZ(ica1,ica2,in)

c	Repulsion between the peptide bond and Ca

        k1=k+1
        dx=(px+x(k1))/2.0-xi
        dy=(py+y(k1))/2.0-yi
        dz=(pz+z(k1))/2.0-zi
        arrr=dx*dx+dy*dy+dz*dz
        if(arrr.lt.54.0) EHL=EHL+EREP

        dx=px-(xi+xi1)/2.0
        dy=py-(yi+yi1)/2.0
        dz=pz-(zi+zi1)/2.0
        arrr=dx*dx+dy*dy+dz*dz
        if(arrr.lt.54.0) EHL=EHL+EREP

        k1=k-1
        dx=(px+x(k1))/2.0-xi
        dy=(py+y(k1))/2.0-yi
        dz=(pz+z(k1))/2.0-zi
        arrr=dx*dx+dy*dy+dz*dz
        if(arrr.lt.54.0) EHL=EHL+EREP

        dx=px-(xi+xim1)/2.0
        dy=py-(yi+yim1)/2.0
        dz=pz-(zi+zim1)/2.0
        arrr=dx*dx+dy*dy+dz*dz
        if(arrr.lt.54.0) EHL=EHL+EREP

c	side group-Ca distance
        if(im.gt.1) then
        aks=(agx-xi)**2+(agy-yi)**2+(agz-zi)**2
        if(aks.lt.43.0) EHL=EHL+EREP
        endif

c	side group-Ca distance
        if(in.gt.1) then
        aks=(bgx-px)**2+(bgy-py)**2+(bgz-pz)**2
        if(aks.lt.43.0) EHL=EHL+EREP
        endif

c	side group- side group distance		
        ar=(agx-bgx)**2+(agy-bgy)**2+(agz-bgz)**2

        fx=CAX(ica1,ica2)
        fy=CAY(ica1,ica2)
        fz=CAZ(ica1,ica2)

c	contact orientation (parallel, othogonal, or antiparallel)		
        aor=hx*fx+hy*fy+hz*fz

        IF(aor.gt.0.5) THEN
        if(ar.gt.arlp(2,in,im,ib,kb)) THEN
        if(ar.lt.alp(2,in,im,ib,kb)) then 
        EHL=EHL+app(kmol,imol,i,k,ib,kb)
        endif
        else
        EHL=EHL+EREP
c	repulsons from side group overlapps weaker than for Ca or Cb
        endif
        ELSE IF(aor.lt.-0.5) THEN
        if(ar.gt.arla(2,in,im,ib,kb)) THEN
        if(ar.lt.ala(2,in,im,ib,kb)) then
        EHL=EHL+apa(kmol,imol,i,k,ib,kb) 
        endif
        else
        EHL=EHL+EREP
        endif
        ELSE
        if(ar.gt.arlm(2,in,im,ib,kb)) THEN
        if(ar.lt.alm(2,in,im,ib,kb)) then
        EHL=EHL+apm(kmol,imol,i,k,ib,kb) 
        endif
        else
        EHL=EHL+EREP
        endif
        ENDIF

        IF(aor.gt.0.0) THEN

        cx=HBX(ica1,ica2)
        cy=HBY(ica1,ica2)
        cz=HBZ(ica1,ica2)

c       COOPERATIVITY OF THE H-bonds, plates of contacting fragments
c       must be almost parallel

        if(prod(nv1,ica1).GT.0.AND.prod(nv2,ica2).GT.0) THEN
c       chains must be parallel or (****)
c       distance dependent attraction of Ca in the H-bond type
c       geometry
        id=1
        ELSE
        if(prod(nv1,ica2).gt.0.OR.prod(nv2,ica1).gt.0) go to 1003
        id=-1
        ENDIF

        ax=(px-xi)
        ay=(py-yi)
        az=(pz-zi)

        delta=0.0
        k1=k+id
        dx=(px+x(k1)-xi-xi1)/2.0
        dy=(py+y(k1)-yi-yi1)/2.0
        dz=(pz+z(k1)-zi-zi1)/2.0
        arrr=dx*dx+dy*dy+dz*dz
        if(arrr.lt.80.0) delta=delta+0.5

        k1=k-id
        dx=(px+x(k1)-xi-xim1)/2.0
        dy=(py+y(k1)-yi-yim1)/2.0
        dz=(pz+z(k1)-zi-zim1)/2.0
        arrr=dx*dx+dy*dy+dz*dz
        if(arrr.lt.80.0) delta=delta+0.5

        ar=(bx-ax)**2+(by-ay)**2+(bz-az)**2
        if(ar.lt.7.0) then
        EHL=EHL+EHBIJ(kmol,imol,i,k)*(0.5+delta)
        else
        ar=(bx+ax)**2+(by+ay)**2+(bz+az)**2
        if(ar.lt.7.0) then
        EHL=EHL+EHBIJ(kmol,imol,i,k)*(0.5+delta)
        endif
        endif

        ar=(cx-ax)**2+(cy-ay)**2+(cz-az)**2
        if(ar.lt.7.0) then
        EHL=EHL+EHBIJ(kmol,imol,i,k)*(0.5+delta)
        else
        ar=(cx+ax)**2+(cy+ay)**2+(cz+az)**2
        if(ar.lt.7.0) then
        EHL=EHL+EHBIJ(kmol,imol,i,k)*(0.5+delta) 
        endif
        endif

        ENDIF

        endif
        endif
        endif
1003    CONTINUE

1004    continue
        ENDDO

        RETURN
        END

C       ***************************************************************

        FUNCTION ERESTR(iiii,jjjj)

        IMPLICIT INTEGER(I-Z)

        PARAMETER(NDIM=500)
        PARAMETER(NREPS=10)
        PARAMETER(NMOLS=2)
        PARAMETER(MAXRES=500)

        COMMON /CHAINS/ ICA(NDIM),X(NDIM),Y(NDIM),Z(NDIM)
        COMMON /SEQE/   SEQ(nmols,ndim),SEC(nmols,ndim)
        COMMON /RCN/    ERESTA,NRESTA
        common /RRCA/   MRESA(nmols,ndim),KRESA(nmols,ndim,maxres),
     &  arca(nmols,ndim,maxres),awrca(nmols,ndim,maxres)

        common /RSGm/ nRSG,eRSG
        common /RSG/  RSGcount(nmols,ndim),RSGind(nmols,ndim,maxres),
     &  distRSG(nmols,ndim,maxres),forceRSG(nmols,ndim,maxres)

        COMMON /SG/ GX(800,800,0:19),GY(800,800,0:19),GZ(800,800,0:19)
        common /reps/   xreps(ndim,NREPS,nmols),
     &                  yreps(ndim,NREPS,nmols),
     &                  zreps(ndim,NREPS,nmols)
        common /lenfs/ mols,imol,itemp,lenf,lenf1,lenf2,lenf3,
     &  al2,al4,al5,al19,length(nmols)

        ERESTR=0.0
        energyRSG=0.0

        DO i=iiii,jjjj

        MM=MRESA(imol,i)
        if(MM.GT.0) THEN
        ix=x(i)
        iy=y(i)
        iz=z(i)
        do kkk=1,MM
        j=kresa(imol,i,kkk)
        if(j.lt.i.OR.j.gt.jjjj) then
        ir=(ix-x(j))**2+(iy-y(j))**2+(iz-z(j))**2
        ar=sqrt(float(ir))
        arc=arca(imol,i,kkk)
        ar=abs(ar-arc)
        if(ar.gt.1.0) erestr=erestr+eresta*ar*awrca(imol,i,kkk)
        endif
        enddo
        endif

        ile=RSGcount(imol,i)
        if(ile.gt.0) then
        xi=x(i)
        yi=y(i)
        zi=z(i)
        nv1=ica(i-1)
        nv2=ica(i)
        im=SEQ(imol,i)
        agx=xi+GX(nv1,nv2,im)
        agy=yi+GY(nv1,nv2,im)
        agz=zi+GZ(nv1,nv2,im)
        do k=1,ile
        j=RSGind(imol,i,k)
        if(j.lt.i.OR.j.gt.jjjj) then
        xj=x(j)
        yj=y(j)
        zj=z(j)
        ica1=ica(j-1)
        ica2=ica(j)
        in=SEQ(imol,j)
        bgx=xj+GX(ica1,ica2,in)
        bgy=yj+GY(ica1,ica2,in)
        bgz=zj+GZ(ica1,ica2,in)
        ir=(agx-bgx)**2+(agy-bgy)**2+(agz-bgz)**2
        ar=sqrt(float(ir))
        arc=distRSG(imol,i,k)

c       ar=abs(ar-arc)
c       if(ar.gt.1.0) energyRSG=energyRSG+eRSG*ar*forceRSG(imol,i,k)

        if(ar.gt.arc) then
        energyRSG=energyRSG+eRSG*(ar-arc)*forceRSG(imol,i,k)
        endif

        endif
        enddo
        endif

        ENDDO

        ERESTR=ERESTR+energyRSG

        return
        end

C       ***************************************************************

        FUNCTION ECROSS(iiii,jjjj)

        IMPLICIT INTEGER(I-Z)

        PARAMETER(NDIM=500)
        PARAMETER(NREPS=10)
        PARAMETER(NMOLS=2)
        PARAMETER(MAXRES=500)

        common /vec/ vector(-7:7,-7:7,-7:7)
        COMMON /SEQE/ SEQ(nmols,ndim),SEC(nmols,ndim)
        COMMON /SG/ GX(800,800,0:19),GY(800,800,0:19),GZ(800,800,0:19)

        COMMON /CHAINS/ ICA(NDIM),X(NDIM),Y(NDIM),Z(NDIM)
        COMMON /RCN/    ERESTA,NRESTA
        common /cross/  mresax(nmols,ndim),kresax(nmols,ndim,maxres),
     &  nresax(nmols,ndim,maxres),arcax(nmols,ndim,maxres),
     &  awrcax(nmols,ndim,maxres)

        common /RSGm/ nRSG,eRSG
        common /RSGX/ RSGXcount(nmols,ndim),RSGXind(nmols,ndim,maxres),
     &  RSGXmol(nmols,ndim,maxres),distRSGX(nmols,ndim,maxres),
     &  forceRSGX(nmols,ndim,maxres)

        common /reps/   xreps(ndim,NREPS,nmols),
     &                  yreps(ndim,NREPS,nmols),
     &                  zreps(ndim,NREPS,nmols)
        common /lenfs/ mols,imol,itemp,lenf,lenf1,lenf2,lenf3,
     &  al2,al4,al5,al19,length(nmols)

        ECROSS=0.0
        energyRSG=0

        DO i=iiii,jjjj

        MM=MRESAX(imol,i)

        if(MM.GT.0) THEN
        ix=x(i)
        iy=y(i)
        iz=z(i)

        do kkk=1,MM
        j=kresax(imol,i,kkk)
        kmol=nresax(imol,i,kkk)
        ir=(ix-xreps(j,itemp,kmol))**2+(iy-yreps(j,itemp,kmol))**2
     &  +(iz-zreps(j,itemp,kmol))**2
        ar=sqrt(float(ir))
        arc=arcax(imol,i,kkk)
        ar=abs(ar-arc)
        if(ar.gt.1.0) ecross=ecross+eresta*ar*awrcax(imol,i,kkk)
        enddo
        endif

        ile=RSGXcount(imol,i)

        if(ile.gt.0) then

        xi=x(i)
        yi=y(i)
        zi=z(i)
        nv1=ica(i-1)
        nv2=ica(i)
        im=SEQ(imol,i)
        agx=xi+GX(nv1,nv2,im)
        agy=yi+GY(nv1,nv2,im)
        agz=zi+GZ(nv1,nv2,im)

        do k=1,ile
        j=RSGXind(imol,i,k)
        jmol=RSGXmol(imol,i,k)

        xj=xreps(j,itemp,jmol)
        yj=yreps(j,itemp,jmol)
        zj=zreps(j,itemp,jmol)

        wx1=xj-xreps(j-1,itemp,jmol)
        wy1=yj-yreps(j-1,itemp,jmol)
        wz1=zj-zreps(j-1,itemp,jmol)

        wx2=xreps(j+1,itemp,jmol)-xj
        wy2=yreps(j+1,itemp,jmol)-yj
        wz2=zreps(j+1,itemp,jmol)-zj

        ica1=vector(wx1,wy1,wz1)
        ica2=vector(wx2,wy2,wz2)

        in=SEQ(jmol,j)
        bgx=xj+GX(ica1,ica2,in)
        bgy=yj+GY(ica1,ica2,in)
        bgz=zj+GZ(ica1,ica2,in)
        ir=(agx-bgx)**2+(agy-bgy)**2+(agz-bgz)**2
        ar=sqrt(float(ir))
        arc=distRSGX(imol,i,k)

c       ar=abs(ar-arc)
c       if(ar.gt.1.0) energyRSG=energyRSG+eRSG*ar*forceRSG(imol,i,k)

        if(ar.gt.arc) then
        energyRSG=energyRSG+eRSG*(ar-arc)*forceRSGX(imol,i,k)
        endif

        enddo
        endif

        ENDDO

        ecross=ecross+energyRSG

        return
        end

C       ***************************************************************

        FUNCTION ESHORT(iiii,jjjj)

        IMPLICIT INTEGER(I-Z)

        PARAMETER(NDIM=500)
        PARAMETER(NREPS=10)
        PARAMETER(NMOLS=2)

        COMMON /CHAINS/  ICA(NDIM),X(NDIM),Y(NDIM),Z(NDIM)
        COMMON /LENGHTS/ PROD(800,800)
        COMMON /SEQE/    SEQ(nmols,ndim),SEC(nmols,ndim)
        COMMON /VECTORS/ vx(800),vy(800),vz(800)
        COMMON /BETAS/   CBX(800,800),CBY(800,800),CBZ(800,800)
        COMMON /BISEC/   CAX(800,800),CAY(800,800),CAZ(800,800)
        COMMON /HB/      HBX(800,800),HBY(800,800),HBZ(800,800)

        common /arandom/ aarand,abrand,acrand,adrand
        COMMON /ENERGY/  EHBOND,ESC,EREP,BAT,
     &  EHBIJ(nmols,nmols,ndim,ndim)
        COMMON /short/  IBIN(-500:500),asr(nmols,ndim,-12:12)
        COMMON /short1/ JBIN(800),bsr(nmols,ndim,16)
        COMMON /short0/ SBIN(200),csr(nmols,ndim,8)

        COMMON /one/    acrit(nmols),eoinp(0:19,0:100)
        common /onecm/  icmx(nmols),icmy(nmols),icmz(nmols)
        COMMON /FR/     FRG(nmols,ndim,ndim)
        
        COMMON /SG/ GX(800,800,0:19),GY(800,800,0:19),GZ(800,800,0:19)
        common /reps/   xreps(ndim,NREPS,nmols),
     &                  yreps(ndim,NREPS,nmols),
     &                  zreps(ndim,NREPS,nmols)
        common /lenfs/ mols,imol,itemp,lenf,lenf1,lenf2,lenf3,
     &  al2,al4,al5,al19,length(nmols)

        COMMON /THREE/   ICONF(800,800)
        COMMON /size/   arla(2,0:19,0:19,2,2),arlm(2,0:19,0:19,2,2),
     &                  arlp(2,0:19,0:19,2,2)

        ESHORT=0.0
        ESC2= ESC/2.0
        ESC1= ESC*4.0
        ESC4= ESC/8.0

        do i=iiii,jjjj
        k=seq(imol,i)
        ii=ICA(i-1)
        jj=ICA(i)
        
        aract2=float((x(i)-icmx(imol))**2+(y(i)-icmy(imol))**2
     &  +(z(i)-icmz(imol))**2)+0.01
        aract=sqrt(aract2)

        ff=aract/acrit(imol)
        ika=int(ff/0.333333)
c       one-body centrosymmetric
        ESHORT=ESHORT+eoinp(k,ika)

c       and the r13 potentials
        kkk=iconf(ii,jj)
        ESHORT=ESHORT+csr(imol,i-1,SBIN(kkk))

        enddo
               
        ii1=iiii-4
        if(ii1.lt.2) ii1=2
        ii2=jjjj+1
        if(ii2.gt.lenf1-4) ii2=lenf1-4
        do i=ii1,ii2
        i1=i+1
        i2=i+2
        i3=i+3
        j=i+4
        icam4=ica(i-1)
        icam3=ica(i)
        icam2=ica(i1)
        icam1=ica(i2)
        icam=ica(i3)
        icaj=ica(j)

c       ff=1.0
c       if(sec(i+1).eq.1) ff=ff-0.33
c       if(sec(i+2).eq.1) ff=ff-0.33 
c       if(sec(i+3).eq.1) ff=ff-0.33

        ax=float(x(i1)+x(i2)+x(i3)+x(j)-4*icmx(imol))/4.0
        ay=float(y(i1)+y(i2)+y(i3)+y(j)-4*icmy(imol))/4.0
        az=float(z(i1)+z(i2)+z(i3)+z(j)-4*icmz(imol))/4.0
        ff=min(1.0,((acrit(imol)*acrit(imol))/(ax*ax+ay*ay+az*az+0.01)))


c       Modiffied the generic stiffness algorithm (5/09/03)

      a=CAX(icam3,icam2)+CAX(icam2,icam1)+CAX(icam1,icam)+CAX(icam,icaj)
      b=CAY(icam3,icam2)+CAY(icam2,icam1)+CAY(icam1,icam)+CAY(icam,icaj)
      c=CAZ(icam3,icam2)+CAZ(icam2,icam1)+CAZ(icam1,icam)+CAZ(icam,icaj)
        aor4=max(sqrt(a*a+b*b+c*c), 0.5) -0.5
        if(aor4.gt.0.5) aor4=0.5
        ESHORT=ESHORT-ff*ESC*(0.5-aor4)

c       Additional generic (13/07/05)

        a=CAX(icam3,icam2)*CAX(icam2,icam1) 
        b=CAY(icam3,icam2)*CAY(icam2,icam1)
        c=CAZ(icam3,icam2)*CAZ(icam2,icam1)
        aor1=max(0.0, (a+b+c))

        a=CAX(icam2,icam1)*CAX(icam1,icam) 
        b=CAY(icam2,icam1)*CAY(icam1,icam)
        c=CAZ(icam2,icam1)*CAZ(icam1,icam)
        aor2=max(0.0, (a+b+c))

        ESHORT=ESHORT+ESC*(aor1+aor2)

        r14=(x(i3)-x(i))**2+(y(i3)-y(i))**2+(z(i3)-z(i))**2
        WX1=VX(icam3)
        WY1=VY(icam3)
        WZ1=VZ(icam3)
        WX2=VX(icam2)
        WY2=VY(icam2)
        WZ2=VZ(icam2)
        WX3=VX(icam1)
        WY3=VY(icam1)
        WZ3=VZ(icam1)
        px=wy1*wz2-wy2*wz1
        py=wz1*wx2-wz2*wx1
        pz=wx1*wy2-wx2*wy1
        ihand=px*wx3+py*wy3+pz*wz3
        if(ihand.lt.0) r14=-r14
        IBIN4=IBIN(r14)
c       contribution from chiral r14 potential
 
        ESHORT=ESHORT+asr(imol,i,IBIN4)
   
        ix=x(j)-x(i)
        iy=y(j)-y(i)
        iz=z(j)-z(i)

        iii=ix**2+iy**2+iz**2
c       contribution from r15 potential

        ESHORT=ESHORT+bsr(imol,i,JBIN(iii))

c       GENERIC BIAS TOWARDS (***) or (###)

        if(iii.gt.80.AND.iii.lt.160) THEN

        if(sec(imol,i+1).ne.4.AND.sec(imol,i+2).ne.4.AND.
     &  sec(imol,i+3).ne.4) then
        if(sec(imol,i).ne.4.AND.sec(imol,j).ne.4) then

        a=CAX(icam3,icam2)*CAX(icam1,icam) 
        b=CAY(icam3,icam2)*CAY(icam1,icam)
        c=CAZ(icam3,icam2)*CAZ(icam1,icam)
        aor2=a+b+c
        if(aor2.lt.0.0) then

c       a helix (***)

        if(IBIN4.gt.4.AND.IBIN4.lt.8) ESHORT=ESHORT-ESC-ff*ESC2
        endif
        endif
        endif

        ELSE
c       a beta (###)
        if(iii.gt.360.AND.iii.lt.530) then

        if(sec(imol,i+1).ne.2.AND.sec(imol,i+2).ne.2.AND.
     &  sec(imol,i+3).ne.2) then
        if(sec(imol,i).ne.2.AND.sec(imol,j).ne.2) then

        a=CAX(icam3,icam2)*CAX(icam1,icam) 
        b=CAY(icam3,icam2)*CAY(icam1,icam)
        c=CAZ(icam3,icam2)*CAZ(icam1,icam)
        aor2=a+b+c
        if(aor2.gt.0.5) then

        a=CAX(icam3,icam2)*CAX(icam2,icam1) 
        b=CAY(icam3,icam2)*CAY(icam2,icam1)
        c=CAZ(icam3,icam2)*CAZ(icam2,icam1)
        aor3=a+b+c
        if(aor3.lt.0.0) then
        a=CAX(icam2,icam1)*CAX(icam1,icam) 
        b=CAY(icam2,icam1)*CAY(icam1,icam)
        c=CAZ(icam2,icam1)*CAZ(icam1,icam)
        aor4=a+b+c
        if(aor4.lt.0.0) then
        if(iabs(IBIN4).gt.8) ESHORT=ESHORT-ESC-ff*ESC2
        endif
        endif
        endif
        endif
        endif
        endif

        ENDIF
        
        enddo

C       ***************************************************************
c       THE OLD STUFF, MODIFIED FOR THE NEW LATTICE

        I1=IIII-14
        if(i1.lt.2) i1=2
        i2=JJJJ
        if(i2.gt.lenf1-15) i2=lenf1-15

        do i=i1,i2
        j=i+5
        k=i+10
        wx1=x(j)-x(i)
        wy1=y(j)-y(i)
        wz1=z(j)-z(i)
        wx2=x(k)-x(j)
        wy2=y(k)-y(j)
        wz2=z(k)-z(j)
        if(wx1*wx2+wy1*wy2+wz1*wz2.lt.0) THEN
        kk=k+5
        wx3=x(kk)-x(k)
        wy3=y(kk)-y(k)
        wz3=z(kk)-z(k)
        if(wx1*wx3+wy1*wy3+wz1*wz3.gt.0) ESHORT=ESHORT+ESC1
c       penalty for "crumpling"
        ENDIF
        enddo

C       ***************************************************************
c       bias to predicted secondary structure geometry

        I1=IIII-7
        if(i1.lt.2) i1=2
        i2=JJJJ
        if(i2.gt.lenf1-7) i2=lenf1-7

        do i=i1,i2
        k=i+7
        m=k-6
        if(sec(imol,m+3).ne.2) then
        ag=FRG(imol,m,k)
        if(ag.gt.0.01) then
        ff=float((x(m)-x(k))**2+(y(m)-y(k))**2+(z(m)-z(k))**2)
        gg=sqrt(ff)
        df=abs(ag-gg)
        if(df.gt.2.0) ESHORT=ESHORT+ESC4*(df-1.0)+ESC2
        endif
        endif
        
        m=k-7
        if(sec(imol,m+3).ne.4) then
        ag=FRG(imol,m,k)
        if(ag.gt.0.01) then
        ff=float((x(m)-x(k))**2+(y(m)-y(k))**2+(z(m)-z(k))**2)
        gg=sqrt(ff)
        df=abs(ag-gg)
        if(df.gt.1.0) ESHORT=ESHORT+ESC4*df+ESC2
        endif 
        endif 

        enddo

        RETURN
        END
C       ***************************************************************
