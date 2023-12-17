                                                                        ACVI7838
      program rmegen                                                    ACVI7839
CB * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *ACVI7840
*                                                                      *ACVI7841
*               ***  Reduced Matrix Element Generator ***              *ACVI7842
*                             ** (rmegen) **                           *ACVI7843
*                                                                      *ACVI7844
CB * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *ACVI7845
C                                                                       ACVI7846
C Author:  Chairul Bahri                                                ACVI7847
C          Department of Physics and Astronomy                          ACVI7848
C          Louisiana State University                                   ACVI7849
C          Baton Rouge LA 70803 USA                                     ACVI7850
C                                                                       ACVI7851
C          e-mail address: phbahr @ lsuvm.sncc.lsu.edu (lsuvm.bitnet)   ACVI7852
C                           bahri @ rouge.phys.lsu.edu                  ACVI7853
C          phone: USA (504)-388-2261                                    ACVI7854
C                     (504)-388-6846                                    ACVI7855
C          fax:   USA (504)-388-5855                                    ACVI7856
C                                                                       ACVI7857
C ----------------------------------------------------------------------ACVI7858
C                                                                       ACVI7859
C Updates: 11/20/90: original ... IBM 3090/                             ACVI7860
C          12/13/90: test for a+ (fp)**8 ... SUCCESSFUL !               ACVI7861
C          01/13/91:          a+a (ds)**2 ... SUCCESSFUL !              ACVI7862
C          02/22/91: multi (lm mu) input                                ACVI7863
C          03/25/91: DONE !                                             ACVI7864
C          04/15/91: FINAL !! ID MUBARAK 1411 H                         ACVI7865
C          08/08/91:          a+a+aa (ds)**6 ... SUCCESSFUL !           ACVI7866
C          10/18/91: G. Rosensteel for U(6)>SU(3) t(2 2) in ds          ACVI7867
C          03/02/92: O. Castanos for double beta decay 238U             ACVI7868
C          11/28/92: IBM RS/6000                                        ACVI7869
C                                                                       ACVI7870
C ----------------------------------------------------------------------ACVI7871
C                                                                       ACVI7872
C General description: PROGRAM 3                                        ACVI7873
C                                                                       ACVI7874
C    Main program for calculating the SU(3) Reduced Coefficients of     ACVI7875
C    Fractional Parentage and/or the SU(3) Reduced Matrix Elements      ACVI7876
C    for a given set of quantum numbers of bra(L), ket(R), and operator ACVI7877
C    (T) using the scheme of SU(3) > SU(2) x U(1).                      ACVI7878
C                                                                       ACVI7879
C    The single barred matrix elements (overlaps) between L,T,R are     ACVI7880
C    evaluated by using bit manipulations in determining the phase.     ACVI7881
C                                                                       ACVI7882
C    The present calculation is for the triple barred                   ACVI7883
C        SU(3) x SU(2)S                                                 ACVI7884
C                                                                       ACVI7885
C ----------------------------------------------------------------------ACVI7886
C  References:                                                          ACVI7887
C    C. Bahri, research notes.                                          ACVI7888
C    J.P.Draayer and Y.Akiyama, J.Math.Phys. 14, 1904 (1973)            ACVI7889
C    Y.Akiyama and J.P.Draayer, Comp.Phys.Comm. 5, 405 (1973)           ACVI7890
C    D.Braunschweig, Comp.Phys.Comm. 14, 109 (1978)                     ACVI7891
C    M.F.O'Reilly, J.Math.Phys. 23, 2022 (1982)                         ACVI7892
C ----------------------------------------------------------------------ACVI7893
      implicit real*8(d)                                                ACVI7894
C                                                                       ACVI7895
      parameter ( IHWSFILE = 8,    ! unformatted hws file               ACVI7896
     1            IOPBFILE = 9,    ! unformatted opb file               ACVI7897
     2            IRMEFILE = 10 )  ! unformatted rme file               ACVI7898
      real title( 18 )             ! title for read only                ACVI7899
      character*20 cifile,         ! input filename                     ACVI7900
     1             chfile,         ! default hws filename               ACVI7901
     2             cofile,         ! default opb filename               ACVI7902
     3             crfile          ! default rme filename               ACVI7903
C                  infile          ! formatted input file unit for sets ACVI7904
C                                  !   of quantum numbers               ACVI7905
C                  lfile           ! log file for intermediate results  ACVI7906
C                  ieta            ! oscillator shell number            ACVI7907
      parameter ( MXKETS =  50000, ! max dimension of kets              ACVI7908
     1            MXSOLN = 100000, ! max dimension of dsolmat           ACVI7909
     2            MXTHWS = 140000 )! tree size                          ACVI7910
      common / HWSARR /                                                 ACVI7911
     1         dsolmat( MXSOLN ),  ! coefficients of hws for bit states ACVI7912
     2         istree(-10:MXTHWS ),! binary tree for quantum labels     ACVI7913
     3         kets( MXKETS )      ! bit states                         ACVI7914
      parameter ( MXWCOF =  50000, ! max dimension of dwsu3             ACVI7915
     1            MXTOPB = 450000, ! tree size for opb                  ACVI7916
     2            MXOPB  = 2000 )  ! max # of coupled operators         ACVI7917
      common / OPBARR /                                                 ACVI7918
     1         iqops( 2, MXOPB ),  ! quantum labels for operators       ACVI7919
     1         isp( MXOPB ),       !                                    ACVI7920
     2         iptr( MXOPB ),      ! pointers to iqops                  ACVI7921
     3         dwsu3( MXWCOF ),    ! coefficients of opb couplings      ACVI7922
     4         lwsu3(-2:MXWCOF ),  ! linked-lists of opb couplings      ACVI7923
     5         ibtree(-10:MXTOPB ),! tree for coupled tensors           ACVI7924
     6         ictree(-10:MXTOPB ) ! tree for uncoupled tensors         ACVI7925
      parameter ( MXRME  = 150000, ! max dimension of drme              ACVI7926
     1            MXTRME = 9*MXRME)! tree size for rme                  ACVI7927
      common / RMEARR /                                                 ACVI7928
     1         drme( MXRME ),      ! SU(3) rme                          ACVI7929
     2         irmtree(-10:MXTRME )! tree for rme                       ACVI7930
      dimension isolmat( 2*MXSOLN ), iwsu3( 2*MXWCOF ), idrme( 2*MXRME )ACVI7931
      equivalence ( isolmat, dsolmat )                                  ACVI7932
      equivalence ( iwsu3,   dwsu3   )                                  ACVI7933
      equivalence ( idrme,   drme    )                                  ACVI7934
C                                                                       ACVI7935
CW**********************************************************************ACVI7936
CW                                                                      ACVI7937
      write(6,'(/a,a/)') ' ENTERING *** REDUCED MATRIX ELEMENTS',       ACVI7938
     1                   ' (RMEGEN) ***'                                ACVI7939
CW                                                                      ACVI7940
CW**********************************************************************ACVI7941
C                                                                       ACVI7942
      call headfile( IHWSFILE, IOPBFILE, IRMEFILE, infile,              ACVI7943
     1               chfile,   cofile,   crfile  , cifile )             ACVI7944
      call heading ( infile, lfile, mltp, np, nh, iqops, iptr, MXOPB )  ACVI7945
      if(infile.ne.5) open( unit = infile, file = cifile ) ! old file   ACVI7946
      if(lfile.gt.0 .and. lfile.ne.6) open( unit = lfile,               ACVI7947
     1                                      file = 'rme.log' )          ACVI7948
      open( unit = IRMEFILE, file = crfile, form = 'unformatted' ) ! newACVI7949
1     format( 18a4 )                                                    ACVI7950
      if( infile .ne. 5 ) then                                          ACVI7951
        read(infile,1,end=999) title                                    ACVI7952
        read(infile,*,end=999) ieta                                     ACVI7953
      else                                                              ACVI7954
        write(6,*) ' Oscillator shell?'                                 ACVI7955
        read (5,*) ieta                                                 ACVI7956
      end if                                                            ACVI7957
      if( lfile .gt. 0 ) write(lfile,'(a/a,i6)')                        ACVI7958
     1  ' Reduced matrix elements for:', ' oscillator shell =', ieta    ACVI7959
      call blocks                      ! SU(3) package                  ACVI7960
C                                                                       ACVI7961
C     Setup trees                                                       ACVI7962
C       OPerator Bases                                                  ACVI7963
      open( unit = IOPBFILE, file = cofile, form = 'unformatted' ) ! oldACVI7964
      if( iptr( MXOPB-2 ) .le. 1 ) then                                 ACVI7965
        nkeys = 4 ! (lm1,mu1,lm2,mu2)(ro12,lm,mu,jt)(lvlp)(lvlh)        ACVI7966
        read(IOPBFILE) neta                             ! signature     ACVI7967
        if( neta.ne.ieta ) call error(' RMEGEN: ieta. Change FT09F001!')ACVI7968
        call wrdata(1, IOPBFILE, itdata, ibtree, -10)   ! labels        ACVI7969
        if( itdata.gt.MXTOPB ) call error(' RMEGEN: Increase MXTOPB!')  ACVI7970
        ibdata      = itdata                                            ACVI7971
        ibtree(-9 ) = MXTOPB / 9 - 1                                    ACVI7972
        ndat        = ibtree(-3 ) - ibtree(-4 )                         ACVI7973
        if( ibtree(-4).ne.nkeys .and. ndat.ne.2 ) call error(           ACVI7974
     1    ' RMEGEN: ibtree. Check FT09F001')                            ACVI7975
        call wrdata(1, IOPBFILE, itdata, ictree, -10)   ! labels        ACVI7976
        if( itdata.gt.MXTOPB ) call error(' RMEGEN: Increase MXTOPB!')  ACVI7977
        icdata      = itdata                                            ACVI7978
        ictree(-9 ) = MXTOPB / 6 - 1                                    ACVI7979
        ndat        = ictree(-3 ) - ictree(-4 )                         ACVI7980
        if( ictree(-4).ne.2 .and. ndat.ne.1 ) call error(               ACVI7981
     1    ' RMEGEN: ictree. Check FT09F001')                            ACVI7982
        call wrdata(1, IOPBFILE, itdata, iwsu3, 1)      ! SU(3) coeffs. ACVI7983
        if( itdata.gt.2*MXWCOF ) call error(' RMEGEN: Increase MXWCOF!')ACVI7984
        id = itdata / 2                                                 ACVI7985
      end if                                                            ACVI7986
      close( unit = IOPBFILE, status = 'keep' )                         ACVI7987
C       Highest Weight State                                            ACVI7988
      open( unit = IHWSFILE, file = chfile, form = 'unformatted' ) ! oldACVI7989
        nkeys = (mltp + 1)/4 + 1                                        ACVI7990
        read(IHWSFILE) neta                             ! signature     ACVI7991
        if( neta.ne.ieta ) call error(' RMEGEN: ieta. Change FT08F001!')ACVI7992
        call wrdata(1, IHWSFILE, isdata, istree, -10)   ! labels        ACVI7993
        if( isdata.gt.MXTHWS ) call error(' RMEGEN: Increase MXTHWS!')  ACVI7994
        istree(-9 ) = MXTHWS / (nkeys+5) - 1                            ACVI7995
        ndat        = istree(-3 ) - istree(-4 )                         ACVI7996
        if( ibtree(-4).ne.nkeys .and. ndat.ne.2 ) call error(           ACVI7997
     1    ' RMEGEN: istree. Check FT08F001')                            ACVI7998
        call wrdata(1, IHWSFILE, ikets, kets, 1)        ! bit states    ACVI7999
        if( ikets.gt.MXKETS ) call error(' RMEGEN: Increase MXKETS!')   ACVI8000
        call wrdata(1, IHWSFILE, isoln, isolmat, 1)     ! bit states    ACVI8001
        if( isoln.gt.2*MXSOLN ) call error(' RMEGEN: Increase MXSOLN!') ACVI8002
      close( unit = IHWSFILE, status = 'keep' )                         ACVI8003
C                                                                       ACVI8004
C     Statistics before rme                                             ACVI8005
      i = -10                         ! index for total elements        ACVI8006
      write(6,100)                                                      ACVI8007
      write(6,110) ' operators ibtree       :', ibtree(i), id           ACVI8008
      write(6,110) ' operators ictree       :', ictree(i), id           ACVI8009
      write(6,110) ' states    istree       :', isdata, ikets, isoln    ACVI8010
C                                                                       ACVI8011
C     Starting ...                                                      ACVI8012
      call rme( infile, lfile, np, nh, iqops, isp, iptr, MXOPB,         ACVI8013
     1  dsolmat,kets,istree, dwsu3,lwsu3,ibtree,ictree, MXWCOF,MXTOPB,  ACVI8014
     2  id, drme,irmtree, MXRME,MXTRME, irme,neta,mltp )                ACVI8015
C     Statistics after rme                                              ACVI8016
      write(6,100)                                                      ACVI8017
      write(6,110) ' operators ibtree       :', ibtree(i), id           ACVI8018
      write(6,110) ' operators ictree       :', ictree(i), id           ACVI8019
      write(6,110) ' rmes      irmtree      :', irmtree(i), irme        ACVI8020
100   format(/' Statistics from trees         :      Nodes       Data') ACVI8021
110   format(7x,a,i11,3x,4i8)                                           ACVI8022
C     ... dump the required data to external file/tree.                 ACVI8023
      iloc  = irmtree(-10)                                              ACVI8024
      iloc  = iloc * ( irmtree(-3 ) + 3 )                               ACVI8025
      irme  = 2 * irme                                                  ACVI8026
      write( IRMEFILE ) ieta                        ! signature         ACVI8027
      call wrdata(0, IRMEFILE, iloc, irmtree, -10)  ! labels            ACVI8028
      call wrdata(0, IRMEFILE, irme, idrme, 1)      ! bit states        ACVI8029
      close( unit = IRMEFILE, status = 'keep' )                         ACVI8030
      if(lfile.gt.0 .and. lfile.ne.6) close( unit = lfile,              ACVI8031
     1                                       status = 'keep' )          ACVI8032
      if(infile.ne.5) close( unit = infile, status = 'keep' )           ACVI8033
999   stop                                                              ACVI8034
C----*end of rmegen*----------------------------------------------------ACVI8035
      end                                                               ACVI8036
      subroutine headfile( ihwsfile, iopbfile, irmefile, infile,        ACVI8037
     1                     chfile,   cofile,   crfile,   cfile )        ACVI8038
C                                                                       ACVI8039
C     Read file structure of rme. (Can be omitted.)                     ACVI8040
C                                                                       ACVI8041
      implicit logical(t)                                               ACVI8042
C                   ihwsfile       ! unformatted input file unit for hwsACVI8043
C                   iopbfile       ! unformatted input file unit for opbACVI8044
C                   irmefile       ! unformatted output file unit fr rmeACVI8045
C                   infile         ! formatted input file unit for sets ACVI8046
C                                  !   of quantum numbers               ACVI8047
      character*(*) cfile,         ! input filename                     ACVI8048
     1              chfile,        ! default hws filename               ACVI8049
     2              cofile,        ! default opb filename               ACVI8050
     3              crfile         ! default rme filename               ACVI8051
      character*20  ctemp          ! temporary name                     ACVI8052
C                                                                       ACVI8053
1     format(a)                                                         ACVI8054
2     format(' =>',10i5)                                                ACVI8055
3     format(' =>',4x,a)                                                ACVI8056
      infile = 15                                                       ACVI8057
      cfile = ' '                                                       ACVI8058
C                                                                       ACVI8059
      write(6,1)         ' Enter input file unit: '                     ACVI8060
      read (5,*)         infile                                         ACVI8061
      write(6,2)         infile                                         ACVI8062
      write(6,*)                                                        ACVI8063
      infile = iabs( infile )                                           ACVI8064
      if( infile.eq.0 .or. infile.eq.6 .or.                             ACVI8065
     1    infile.eq.4 .or. infile.eq.ihwsfile .or.                      ACVI8066
     2    infile.eq.iopbfile .or. infile.eq.irmefile ) then             ACVI8067
        write(6,1) ' Sorry, file unit 15 is assigned. '                 ACVI8068
        infile = 15                                                     ACVI8069
      endif                                                             ACVI8070
      if( infile .ne. 5 ) then                                          ACVI8071
        write(6,1)       ' Enter input file name: '                     ACVI8072
        read (5,'(a)')   ctemp                                          ACVI8073
        write(6,3)       ctemp                                          ACVI8074
        call compchr( ctemp )                                           ACVI8075
        if( .not. tchkdt( ctemp ) ) then                                ACVI8076
          chfile = 'hwsw.' // ctemp(1:14)                               ACVI8077
          cofile = 'opbw.' // ctemp(1:14)                               ACVI8078
          crfile = 'rmew.' // ctemp(1:14)                               ACVI8079
          cfile  = 'hwsir.' // ctemp(1:14)                              ACVI8080
        else                                                            ACVI8081
          chfile = 'hwsw.out'                                           ACVI8082
          cofile = 'opbw.out'                                           ACVI8083
          crfile = 'rmew.out'                                           ACVI8084
          cfile  = ctemp                                                ACVI8085
        endif                                                           ACVI8086
      else                                                              ACVI8087
        chfile = 'hwsw.out'                                             ACVI8088
        cofile = 'opbw.out'                                             ACVI8089
        crfile = 'rmew.out'                                             ACVI8090
      endif                                                             ACVI8091
      write(6,3) chfile                                                 ACVI8092
      write(6,3) cofile                                                 ACVI8093
      write(6,3) crfile                                                 ACVI8094
      write(6,3) cfile                                                  ACVI8095
      write(6,*)                                                        ACVI8096
C                                                                       ACVI8097
      return                                                            ACVI8098
C----*end of headfile*--------------------------------------------------ACVI8099
      end                                                               ACVI8100
      subroutine heading( infile,lfile, mltp, np,nh, iqops,iptr,mxnumb )ACVI8101
C                                                                       ACVI8102
C     Read input parameters for rme.                                    ACVI8103
C                                                                       ACVI8104
      implicit logical(t)                                               ACVI8105
C               lfile              ! log file for intermediate results  ACVI8106
C               mltp               ! internal symmetry multiplicity     ACVI8107
C               np                 ! # of a+'s                          ACVI8108
C               nh                 ! # of a's                           ACVI8109
      dimension iqops( 2, * ),     ! temporary variables                ACVI8110
     1          iptr( * )          ! temporary variables                ACVI8111
C               mxnumb             ! temporary index                    ACVI8112
C                                                                       ACVI8113
      parameter ( IBIG = 1000 )                                         ACVI8114
1     format(a)                                                         ACVI8115
2     format(' =>',10i5)                                                ACVI8116
      write(6,1)         ' Enter log file unit number:'                 ACVI8117
      write(6,1)         '       0: none'                               ACVI8118
      write(6,1)         '     1-3: some intermediate results'          ACVI8119
      write(6,1)         '     >=7: all intermediate results'           ACVI8120
      read (5,*,end=999) lfile                                          ACVI8121
      write(6,2)         lfile                                          ACVI8122
      write(6,*)                                                        ACVI8123
      write(6,1)         ' Enter number of multiplicity for one level:' ACVI8124
      read (5,*,end=999) mltp                                           ACVI8125
      write(6,2)         mltp                                           ACVI8126
      write(6,1)         ' Enter operator structure:'                   ACVI8127
      write(6,1)         '       0: read everything from the file'      ACVI8128
      write(6,1)         '       1: read a+ only from the file'         ACVI8129
      write(6,1)         '     >=2: generate everything on the fly'     ACVI8130
      write(6,*)                                                        ACVI8131
      read (5,*,end=999) iop                                            ACVI8132
      write(6,2)         iop                                            ACVI8133
      write(6,1)         ' Enter # of a+_s, a_s:'                       ACVI8134
      read (5,*,end=999) np, nh                                         ACVI8135
      write(6,2)         np, nh                                         ACVI8136
      write(6,*)                                                        ACVI8137
      write(6,1)         ' Enter the tensor character of operators:'    ACVI8138
      write(6,1)         '       0: L=0 --> scalar'                     ACVI8139
      write(6,1)         '       1: L=1 --> vector'                     ACVI8140
      write(6,1)         '    else: all tensors'                        ACVI8141
      read (5,*,end=999) itensor                                        ACVI8142
      write(6,2)         itensor                                        ACVI8143
      write(6,*)                                                        ACVI8144
      write(6,1)         ' Do you want limited numbers of (lm mu)?'     ACVI8145
      write(6,1)         ' (0 for no, 1 for yes)'                       ACVI8146
      read (5,*,end=999) limitlm                                        ACVI8147
      write(6,2)         limitlm                                        ACVI8148
      write(6,*)                                                        ACVI8149
      write(6,1)         ' Enter minimum (lm mu):'                      ACVI8150
      read (5,*,end=999) lmopmn, muopmn                                 ACVI8151
      write(6,2)         lmopmn, muopmn                                 ACVI8152
      write(6,*)                                                        ACVI8153
      write(6,1)         ' Enter maximum (lm mu):'                      ACVI8154
      read (5,*,end=999) lmopmx, muopmx                                 ACVI8155
      write(6,2)         lmopmx, muopmx                                 ACVI8156
      write(6,*)                                                        ACVI8157
      write(6,1)         ' lm = mu ?'                                   ACVI8158
      read (5,*,end=999) iequal                                         ACVI8159
      write(6,2)         iequal                                         ACVI8160
      write(6,*)                                                        ACVI8161
      if( limitlm .le. 0 ) then                                         ACVI8162
        lmopmn = 0                                                      ACVI8163
        muopmn = 0                                                      ACVI8164
        lmopmx = IBIG                                                   ACVI8165
        muopmx = IBIG                                                   ACVI8166
        iequal = 0                                                      ACVI8167
      end if                                                            ACVI8168
C                                                                       ACVI8169
      if(lmopmn.lt.0 .or. muopmn.lt.0 .or. lmopmx.lt.0 .or. muopmx.lt.0)ACVI8170
     1  call error(' HEADING: positive (lm mu)s')                       ACVI8171
C                                                                       ACVI8172
      lfile = iabs( lfile )                                             ACVI8173
      if( (lfile.gt.3 .and. lfile.lt.6).and.(lfile.ge.8 .and. lfile.le. ACVI8174
     1     10) ) call error(' HEADING: Fileunit has been assigned')     ACVI8175
      iqops( 1, mxnumb-1 ) = lmopmn                                     ACVI8176
      iqops( 2, mxnumb-1 ) = muopmn                                     ACVI8177
      iqops( 1, mxnumb )   = lmopmx                                     ACVI8178
      iqops( 2, mxnumb )   = muopmx                                     ACVI8179
      iptr ( mxnumb-2 )    = iop                                        ACVI8180
      iptr ( mxnumb-1 )    = itensor                                    ACVI8181
      iptr ( mxnumb )      = iequal                                     ACVI8182
C                                                                       ACVI8183
      return                                                            ACVI8184
999   write(6,1) ' Please enter numbers, otherwise I quit.'             ACVI8185
C----*end of heading*-------------------------------------------------- ACVI8186
      end                                                               ACVI8187
      LOGICAL FUNCTION TCHKDT( CDOT )                                   ACVI8188
C                                                                       ACVI8189
C     Check whether there is a dot in a string.                         ACVI8190
C                                                                       ACVI8191
      CHARACTER*(*) CDOT                                                ACVI8192
      TCHKDT = .FALSE.                                                  ACVI8193
C     NLEN = 0                                                          ACVI8194
C     DO WHILE( .NOT. EOL )                                             ACVI8195
C       NLEN = NLEN + 1                                                 ACVI8196
      DO NLEN = 1, 20                                                   ACVI8197
        IF( CDOT( NLEN:NLEN ) .EQ. '.' ) THEN                           ACVI8198
          TCHKDT = .TRUE.                                               ACVI8199
          RETURN                                                        ACVI8200
        END IF                                                          ACVI8201
      END DO                                                            ACVI8202
      RETURN                                                            ACVI8203
      END                                                               ACVI8204
C                                                                       ACVI8205
      SUBROUTINE COMPCHR( CHARC )                                       ACVI8206
C                                                                       ACVI8207
C     Delete blanks before the first character.                         ACVI8208
C                                                                       ACVI8209
      CHARACTER*(*) CHARC                                               ACVI8210
      NLEN = 1                                                          ACVI8211
      DO WHILE( CHARC(NLEN:NLEN).EQ.' ' .AND. NLEN.LE.20 )              ACVI8212
        NLEN = NLEN + 1                                                 ACVI8213
      ENDDO                                                             ACVI8214
      CHARC = CHARC( NLEN: )                                            ACVI8215
      RETURN                                                            ACVI8216
      END                                                               ACVI8217
C                                                                       ACVI8218
      SUBROUTINE WRTOUT(LFILE,ICALL,MLTP,IQNUMS,DASH,LITER)             ACVI8219
C                                                                       ACVI8220
C     Write out the heading for each irreps.                            ACVI8221
C                                                                       ACVI8222
      DIMENSION IQNUMS(*)                                               ACVI8223
      CHARACTER DASH(78)                                                ACVI8224
      CHARACTER*(*) LITER                                               ACVI8225
      WRITE(LFILE,'(A)') LITER                                          ACVI8226
      WRITE(LFILE,'(A,I5)') ' STATE NUMBER ', ICALL                     ACVI8227
      WRITE(LFILE,*) DASH                                               ACVI8228
      WRITE(LFILE,1000) (' ',I,I=1,MLTP)                                ACVI8229
      WRITE(LFILE,1010) (IQNUMS(I),I=1,MLTP+2)                          ACVI8230
      WRITE(LFILE,*) DASH                                               ACVI8231
1000  FORMAT(' ',T6,'Lm',T12,'Mu',T15,10(A,2X,'f',I1,'~'))              ACVI8232
1010  FORMAT(' ',12I6)                                                  ACVI8233
      RETURN                                                            ACVI8234
      END                                                               ACVI8235
CB * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *ACVI8236
C                                                                      *ACVI8237
C                    *** REDUCED MATRIX ELEMENTS ***                   *ACVI8238
C                            *** ( RME ) ***                           *ACVI8239
C                                                                      *ACVI8240
CB * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *ACVI8241
C                                                                       ACVI8242
C Author:  Chairul Bahri                                                ACVI8243
C          Department of Physics and Astronomy                          ACVI8244
C          Louisiana State University                                   ACVI8245
C          Baton Rouge LA 70803 USA                                     ACVI8246
C                                                                       ACVI8247
C          E-MAIL: PHBAHR @ LSUVM.BITNET or LSUVM.SNCC.LSU.EDU          ACVI8248
C                  bahri  @ rouge.phys.lsu.edu                          ACVI8249
C          PHONE:  USA (504)-388-2261                                   ACVI8250
C                      (504)-388-6846                                   ACVI8251
C          FAX:    USA (504)-388-5855                                   ACVI8252
C                                                                       ACVI8253
C ----------------------------------------------------------------------ACVI8254
C                    WELCOME TO THE WORLD, MY SON: MUHIYYUDDIN          ACVI8255
C ----------------------------------------------------------------------ACVI8256
C                                                                       ACVI8257
C Updates: 11/20/90: original ... IBM 3090/                             ACVI8258
C          12/13/90: first test (CFP) ... SUCCESSFUL !                  ACVI8259
C          01/13/91: second test (one-body) ... SUCCESSFUL !            ACVI8260
C          02/22/91: multi (LM MU) input.                               ACVI8261
C          04/15/91: FINAL !! ID MUBARAK 1411                           ACVI8262
C          08/08/91: test for a+a+aa(ds)**6 ... SUCCESSFUL !            ACVI8263
C          12/16/91: SUN 4, subroutine form.                            ACVI8264
C                                                                       ACVI8265
C ----------------------------------------------------------------------ACVI8266
C                                                                       ACVI8267
C General description: SUB-PROGRAM 3                                    ACVI8268
C                                                                       ACVI8269
C    Main subprogram for calculating the SU(3) Reduced Coefficient of   ACVI8270
C    Fractional Parentage and/or the SU(3) Reduced Matrix Elements      ACVI8271
C    a given set of quantum numbers of BRA(L), KET(R), and OPERATOR(T)  ACVI8272
C    using the scheme SU(3) > SU(2) x U(1).                             ACVI8273
C                                                                       ACVI8274
C    The single barred matrix elements (overlaps) between L,R,T are     ACVI8275
C    evaluated by using bit manipulations (for determining the phase.)  ACVI8276
C                                                                       ACVI8277
C    The present calculation is for the triple barred                   ACVI8278
C        SU(3) x SU(2)                                                  ACVI8279
C                                                                       ACVI8280
C ----------------------------------------------------------------------ACVI8281
C References:                                                           ACVI8282
C    J.P.Draayer and Y.Akiyama, J.Math.Phys. 14, 1904 (1973)            ACVI8283
C    Y.Akiyama and J.P.Draayer, Comp.Phys.Comm. 5, 405 (1973)           ACVI8284
C    D.Braunschweig, Comp.Phys.Comm. 14, 109 (1978)                     ACVI8285
C    D.Braunschweig, Comp.Phys.Comm. 14, 109 (1978)                     ACVI8286
C    M.F.O'Reilly, J.Math.Phys. 23, 2022 (1982)                         ACVI8287
C    C. Bahri, research notes.                                          ACVI8288
C ----------------------------------------------------------------------ACVI8289
      subroutine rme( infile, logfile, np, nh, iqops, isp, iptr, mxopb, ACVI8290
     ]  dsolmat, kets, istree, dwsu3,lwsu3, ibtree,ictree, mxcop, mxtop,ACVI8291
     ]  id, drme, irmtree, mxrme, mxtrme, irme, neta, mltp )            ACVI8292
      implicit real*8(d), logical(t)                                    ACVI8293
C                                                                       ACVI8294
      character dash(68)           ! ---- ----                          ACVI8295
      logical   tfile              ! .true. need a log file             ACVI8296
C               infile             ! formatted input file unit for sets ACVI8297
C                                  !   of quantum numbers               ACVI8298
C               logfile            ! log file for intermediate results  ACVI8299
C               np, nh             ! # of a+'s and a's                  ACVI8300
      dimension iqops( 2, * ),     ! SU(3) coupled operators            ACVI8301
     1          isp( * ),          ! SU(2)spin coupled operators        ACVI8302
     2          iptr( * )          ! pointers to iqops(*,*)             ACVI8303
C               mxopb              ! max dimension of iqops(2,*),iptr(*)ACVI8304
      dimension dsolmat( * ),      ! coefficients of hws for bit states ACVI8305
     1          kets( * ),         ! bit states                         ACVI8306
     2          istree( -10:* )    ! binary tree for state labels       ACVI8307
      dimension dwsu3( * ),        ! SU(3) Wigner coefficients for op.  ACVI8308
     1          lwsu3( -2:* ),     ! linked-lists for operators         ACVI8309
     2          ibtree( -10:* ),   ! wst for SU(3) coupled tensors      ACVI8310
     3          ictree( -10:* )    ! wst for SU(3) uncoupled tensors    ACVI8311
C               mxcop              ! max dimension of dwsu3             ACVI8312
C               mxtop              ! max dimension of iotree            ACVI8313
      dimension drme( * ),         ! SU(3) rmes                         ACVI8314
     1          irmtree( -10:* )   ! binary tree for rmes               ACVI8315
C               mxrme              ! max dimension of drme              ACVI8316
C               mxtrme             ! max dimension of irmtree           ACVI8317
C               irme               ! # elements in drme                 ACVI8318
C               neta               ! oscillator shell number            ACVI8319
C               mltp               ! internal symmetry multiplicity     ACVI8320
      character*115 bin            ! character for binary reps.         ACVI8321
      common / LOGCON / tfile,     ! write to logfile (full)            ACVI8322
     1                  tsmfi,     ! write to logfile (small)           ACVI8323
     2                  tconj,     ! conjugate SU(3) label              ACVI8324
     3                  topbase    ! calculate SU(3) coupled opb        ACVI8325
      common / RMECON / lfile,     ! = logfile                          ACVI8326
     1                  ieta       ! = neta                             ACVI8327
C                                                                       ACVI8328
      real title( 18 )             ! title ( temporary )                ACVI8329
      dimension iqnums( 24 ),      ! quantum number array for read in   ACVI8330
     1          iqwrit( 65 ),      ! quantum number array for write out ACVI8331
     2          label( 4 )         ! temp array of packed q.n. for tree ACVI8332
C    2          mred( 10 )         ! subgroup chain of the int. group   ACVI8333
      dimension iket( 2 ),         ! kets indices for bra and ket       ACVI8334
     1          nbkets( 2 ),       ! kets-max for bra and ket           ACVI8335
     2          isol( 2 ),         ! dsolmat indices for bra and ket    ACVI8336
     3          nalpha( 2 )        ! alpha-max for bra and ket          ACVI8337
C                                                                       ACVI8338
      data nbit16 / 65535 /                                             ACVI8339
C     data nbit16 / zffff /        ! for IBM-3090 only                  ACVI8340
C     Unpacking function.                                               ACVI8341
      iupac( index, iover, ibits ) = iand( ishft( index, iover), ibits )ACVI8342
C                                                                       ACVI8343
C     File and tree initializations ...                                 ACVI8344
C                                                                       ACVI8345
      lfile = logfile                                                   ACVI8346
      tfile = lfile .gt. 0                                              ACVI8347
      iop   = iptr( mxopb-2 )                                           ACVI8348
      nop   = max0( np, nh )                                            ACVI8349
      topbase = iop .gt. 0                                              ACVI8350
C                                                                       ACVI8351
7     format(18a4)                                                      ACVI8352
      if( infile.ne.5 ) then                                            ACVI8353
        rewind(infile)                                                  ACVI8354
        read(infile,7,end=999) title         ! shell number             ACVI8355
        read(infile,*,end=999) ieta                                     ACVI8356
      else                                                              ACVI8357
        ieta = neta                                                     ACVI8358
      end if                                                            ACVI8359
      if( ieta .ne. neta ) call error(' RME: Wrong shell ieta !!')      ACVI8360
      if(tfile) write(lfile,'(a/a,i6)') ' Reduced matrix elements for:',ACVI8361
     1          ' oscillator shell =', ieta                             ACVI8362
C                                                                       ACVI8363
      n = mxtrme / 9 - 1                                                ACVI8364
      call tsetlf( irmtree, n, 5, 1 )                                   ACVI8365
C                              |  |                                     ACVI8366
C                              |  |-- data                              ACVI8367
C                              |----- key                               ACVI8368
C     initialize linked-lists for opb                                   ACVI8369
      lwsu3(-2 ) = mxcop                                                ACVI8370
      lwsu3(-1 ) = id + 1                                               ACVI8371
      lwsu3( 0 ) = mxcop - id                                           ACVI8372
      do i = 1, mxcop-1                                                 ACVI8373
        lwsu3( i ) = i + 1                                              ACVI8374
      end do                                                            ACVI8375
      lwsu3( mxcop ) = -1                                               ACVI8376
C                                                                       ACVI8377
C     Start ...                                                         ACVI8378
C                                                                       ACVI8379
      mltpp2  = mltp + 2                                                ACVI8380
      mltpmax = mltpp2 + mltpp2                                         ACVI8381
      n       = mod( mltpp2, 4 )                                        ACVI8382
      if( n .eq. 0 ) n = 4                                              ACVI8383
      nkeys   = (mltp + 1)/4 + 1                                        ACVI8384
      nkeym1  = nkeys - 1                                               ACVI8385
      ik1 = mltpp2 + 1                                                  ACVI8386
      do i = 1, 68                                                      ACVI8387
        dash( i ) = ' '                                                 ACVI8388
      end do                                                            ACVI8389
      idash = 4 * mltpp2 - 1                                            ACVI8390
      j     = idash + 2                                                 ACVI8391
      do i = 2, idash                                                   ACVI8392
        j = j + 1                                                       ACVI8393
        dash( i ) = '-'                                                 ACVI8394
        dash( j ) = '-'                                                 ACVI8395
      end do                                                            ACVI8396
      dash( 1 ) = '*'                                                   ACVI8397
      dash( j+1 ) = '*'                                                 ACVI8398
C     call setqns( mltp, mred, itype, jk1 )                             ACVI8399
C                                                                       ACVI8400
      irme = 0                                                          ACVI8401
      call su3u2fp( ieta, nop, iqops, isp, iptr, mxopb )                ACVI8402
      call su3tlm( np, -nh, iqops, isp, iptr, mxopb )                   ACVI8403
      call setrme( mltp, iqops, mxopb, mxjtck )                         ACVI8404
      call setbin( lfile, .true., bin, 11, -1, 35, 115 )                ACVI8405
C                                                                       ACVI8406
C     Read SU(3)xSU(M) quantum numbers.                                 ACVI8407
      if( infile.eq.5 ) then                                            ACVI8408
        write(6,*) ' Lm Mu f~''s'                                       ACVI8409
      else                                                              ACVI8410
        read(infile,7,end=999) title                                    ACVI8411
      endif                                                             ACVI8412
C                                                                       ACVI8413
      do while( irme.le.mxrme )   ! * * * * * * * * * * * * * * * * * * ACVI8414
C                                                                       ACVI8415
        if( infile.eq.5 ) then                                          ACVI8416
          read(infile,*) (iqnums(i),i=1,mltpmax)                        ACVI8417
        else                                                            ACVI8418
          read(infile,*,end=999) (iqnums(i),i=1,mltpmax)                ACVI8419
        end if                                                          ACVI8420
        call ritout( 6, dash, iqnums, mltpmax )                         ACVI8421
        if( tfile ) call ritout( lfile, dash, iqnums, mltpmax )         ACVI8422
C                                                                       ACVI8423
        npbra = 0                                                       ACVI8424
        npket = 0                                                       ACVI8425
        do i = 3, mltpp2                                                ACVI8426
          npbra = npbra + iqnums( i )                                   ACVI8427
          npket = npket + iqnums( mltpp2 + i )                          ACVI8428
        end do                                                          ACVI8429
        if( npbra .lt. np ) call attn(' RME: Trivial solution.',*120 )  ACVI8430
        if( npket .lt. nh ) call attn(' RME: Trivial solution.',*120 )  ACVI8431
        if( npbra-npket .ne. np-nh )                                    ACVI8432
     1      call attn(' RME: Number of particles not match.',*120 )     ACVI8433
C                                                                       ACVI8434
        iq = 0                                                          ACVI8435
C                                                                       ACVI8436
        do k = 1, 2                                                     ACVI8437
C         Pack labels: each 4.                                          ACVI8438
          do i = 1, nkeys                                               ACVI8439
            label( i ) = 0                                              ACVI8440
          end do                                                        ACVI8441
          do i = 1, nkeym1                                              ACVI8442
            do j = 1, 4                                                 ACVI8443
              iq = iq + 1                                               ACVI8444
              label( i ) = ior( ishft(label(i),8), iqnums(iq) )         ACVI8445
            end do                                                      ACVI8446
          end do                                                        ACVI8447
          do j = 1, n                                                   ACVI8448
            iq = iq + 1                                                 ACVI8449
            label(nkeys) = ior(ishft(label(nkeys),8),iqnums(iq))        ACVI8450
          end do                                                        ACVI8451
C         Retrieve ISTREE data                                          ACVI8452
          call tchk( label, istree, *110 )                              ACVI8453
            write(6,'(a,i3)') ' missing irrep', k                       ACVI8454
            go to 120                                                   ACVI8455
110       iloc = istree(-5 ) + istree(-4 ) + 1                          ACVI8456
          index = istree( iloc )                                        ACVI8457
          iket( k )   = iupac( index, -16, nbit16 ) + 1                 ACVI8458
          nbkets( k ) = iupac( index,   0, nbit16 )                     ACVI8459
          index = istree( iloc + 1 )                                    ACVI8460
          isol( k )   = iupac( index, -16, nbit16 ) + 1                 ACVI8461
          nalpha( k ) = iupac( index,   0, nbit16 )                     ACVI8462
        end do                                                          ACVI8463
C                                                                       ACVI8464
C       if( itype .eq. 1 ) then                                         ACVI8465
          call qnsu2( iqnums, iqwrit )        ! transfer qns            ACVI8466
C       else if( itype .eq. 2 ) then                                    ACVI8467
C         call qnsu422( ) then                                          ACVI8468
C       else                                                            ACVI8469
C       end if                                                          ACVI8470
C                                                                       ACVI8471
C       Calculate the rme of T(lm mu)                                   ACVI8472
        nal = nalpha( 1 )                                               ACVI8473
        nar = nalpha( 2 )                                               ACVI8474
        call wigmat( iqops, iqnums, index, mxjtck, mltp, *999 )         ACVI8475
        if( index .gt. 0 ) then                                         ACVI8476
          call me1bar( iqops, isp, index, nbkets(1), nbkets(2),         ACVI8477
     1                 nal, nar, kets(iket(1)), kets(iket(2)),          ACVI8478
     2                 dsolmat(isol(1)), dsolmat(isol(2)), mxjtck,      ACVI8479
     3                 mltp, ibtree, ictree, lwsu3, dwsu3, mxcop )      ACVI8480
          call rm3bar( iqops, isp, iqwrit, nal, nar,                    ACVI8481
     1                 drme, irmtree, mxrme, irme, mxjtck, *999 )       ACVI8482
        end if                                                          ACVI8483
C                                                                       ACVI8484
120     continue                                                        ACVI8485
        if( infile .eq. 5 ) write(6,*) ' Next ...'                      ACVI8486
C                                                                       ACVI8487
        id = max0( id, lwsu3(-1)-1 )                                    ACVI8488
      enddo                       ! * * * * * * * * * * * * * * * * * * ACVI8489
999   continue                                                          ACVI8490
      return                                                            ACVI8491
C ---*end of rme*-------------------------------------------------------ACVI8492
      end                                                               ACVI8493
      subroutine ritout( lfile, dash, iqnums, nmax )                    ACVI8494
      dimension iqnums( * )                                             ACVI8495
      character dash( 68 )                                              ACVI8496
      write(lfile,*)                                                    ACVI8497
      write(lfile,*) dash                                               ACVI8498
      write(lfile,'(24i4)') (iqnums( i ), i=1, nmax)                    ACVI8499
      write(lfile,*) dash                                               ACVI8500
      return                                                            ACVI8501
      end                                                               ACVI8502
C ----------------------------------------------------------------------ACVI8503
C                                                                       ACVI8504
      subroutine su3u2fp( ieta, nx, lmmua, ispta, iptr, mxnumb )        ACVI8505
C                                                                       ACVI8506
C ----------------------------------------------------------------------ACVI8507
C Author: Chairul Bahri (LSU 11/90 ... original)                        ACVI8508
C ----------------------------------------------------------------------ACVI8509
C     Sub-program to generate all SU(3)xSU(2) quantum numbers of coupledACVI8510
C     fermion a+'s.                                                     ACVI8511
C ----------------------------------------------------------------------ACVI8512
C                                                                       ACVI8513
C               ieta              ! shell number                        ACVI8514
C               nx                ! generation level                    ACVI8515
      dimension lmmua( 2, * ),    ! array of packed (lm,mu) of SU(3)    ACVI8516
     1          ispta( * ),       ! array of spin of SU(2)              ACVI8517
     2          iptr( * )         ! level pointers for generations      ACVI8518
C               mxnumb            ! max # (lm mu) pairs lmmua           ACVI8519
C     Packing functions                                                 ACVI8520
      data nbit8 / 255 /           ! zff                                ACVI8521
      ipack( i, j, k, l ) = ior(l,ishft(ior(k,ishft(ior(j,ishft(i,      ACVI8522
     1       8)),8)),8))                                                ACVI8523
      iupac( index, iover, ibits) = iand( ishft( index, iover), ibits ) ACVI8524
C                                                                       ACVI8525
C     Initialize lmmua                                                  ACVI8526
      do i = 1, mxnumb-2                                                ACVI8527
        lmmua( 1, i ) = 0                                               ACVI8528
        lmmua( 2, i ) = 0                                               ACVI8529
      end do                                                            ACVI8530
C                                                                       ACVI8531
C     First generation                                                  ACVI8532
      lmi = ieta                                                        ACVI8533
      mui = 0                                                           ACVI8534
                lmmua( 1, 1 ) = ipack( lmi, mui, 0, 0 )                 ACVI8535
                lmmua( 2, 1 ) = ipack( 1, lmi, mui, 0 )                 ACVI8536
                ispta( 1 )    = 1                                       ACVI8537
        iptr( 1 ) = 1                                                   ACVI8538
        l         = 1                                                   ACVI8539
        j         = l                                                   ACVI8540
C     write(6,110) 1                                                    ACVI8541
C     write(6,120) 1, 0,0, ieta,0, 1,ieta,0                             ACVI8542
C                                                                       ACVI8543
C     Second generation                                                 ACVI8544
      if( nx .lt. 2 ) go to 130                                         ACVI8545
      lmmuk  = ipack( lmi, mui, ieta, 0 )                               ACVI8546
      ispin0 = 0                                                        ACVI8547
      ispin1 = 2                                                        ACVI8548
      if( ieta .gt. mxnumb-2 ) call attn(' SU3U2FP: Too many ops.',*999)ACVI8549
              do mu = ieta, 0, -1                                       ACVI8550
                l  = l + 1                                              ACVI8551
                lm = 2 * (ieta - mu)                                    ACVI8552
                lmmua( 1, l ) = lmmuk                                   ACVI8553
                lmmua( 2, l ) = ipack( 1, lm, mu, 0 )                   ACVI8554
                if( btest(mu,0) ) then                                  ACVI8555
                  ispta( l ) = ispin1                                   ACVI8556
C                 ispta( l ) = ipack( 1,  0,  0, 2 )                    ACVI8557
                else                                                    ACVI8558
                  ispta( l ) = ispin0                                   ACVI8559
C                 ispta( l ) = ipack( 1,  0,  0, 0 )                    ACVI8560
                end if                                                  ACVI8561
              end do                                                    ACVI8562
        iptr( 2 ) = 2                                                   ACVI8563
        j         = l                                                   ACVI8564
C                                                                       ACVI8565
C     Next generations (waiting for ... UNU3)                           ACVI8566
      do i = 3, nx                                                      ACVI8567
CW      write(6,110) i                                                  ACVI8568
        do k = iptr( i-1 ), j                                           ACVI8569
          lmmuk = lmmua( 2, k )                                         ACVI8570
          lmk   = iupac( lmmuk, -16, nbit8 )                            ACVI8571
          muk   = iupac( lmmuk,  -8, nbit8 )                            ACVI8572
          lam   = lmi + mui + lmk + muk                                 ACVI8573
          do lm = 0, lam                                                ACVI8574
            do mu = 0, lam                                              ACVI8575
              call u3mult( lmi, mui, lmk, muk, lm, mu, kromax, *100 )   ACVI8576
                l = l + 1                                               ACVI8577
                if( l .gt. mxnumb-2 )                                   ACVI8578
     1            call attn(' SU3U2FP: Too many ops.',*999 )            ACVI8579
                lmmua( 1, l ) = ipack( lmi, mui, lmk, muk )             ACVI8580
                lmmua( 2, l ) = ipack( kromax, lm, mu, 0 )              ACVI8581
C**             ispta( l )    = ipack( k, 0, 0, ispt )                  ACVI8582
CW              write(6,120) l, lmi,mui, lmk,muk, kromax,lm,mu          ACVI8583
100           continue                                                  ACVI8584
            end do                                                      ACVI8585
          end do                                                        ACVI8586
        end do                                                          ACVI8587
        iptr( i ) = j + 1                                               ACVI8588
        j         = l                                                   ACVI8589
      end do                                                            ACVI8590
110   format(' Generation ',i3)                                         ACVI8591
120   format(8x,' # ',i3,' (',2i3,')x(',2i3,') .. ',i2,'(',2i3,')')     ACVI8592
130   iptr( nx + 1 ) = l + 1                                            ACVI8593
CW    write(6,'(a,i3,a,i4)') ' last pointer(',nx+1,') ',iptr(nx+1)      ACVI8594
999   return                                                            ACVI8595
C ---*end of su3fp*-----------------------------------------------------ACVI8596
      end                                                               ACVI8597
C ----------------------------------------------------------------------ACVI8598
C                                                                       ACVI8599
      subroutine su3tlm( n1, n2, lmmua, ispta, iptr, mxnumb )           ACVI8600
C                                                                       ACVI8601
C ----------------------------------------------------------------------ACVI8602
C Author: Chairul Bahri (LSU 11/90 ... original IBM 3090 ... SU3OP)     ACVI8603
C ----------------------------------------------------------------------ACVI8604
C     Sub-program to generate some SU(3) quantum numbers of coupled     ACVI8605
C     tensors (t1 x t2)T(lm mu).                                        ACVI8606
C ----------------------------------------------------------------------ACVI8607
C                                                                       ACVI8608
C               n1                ! number op. for t1                   ACVI8609
C               n2                ! number op. for t2                   ACVI8610
      logical tenchk                                                    ACVI8611
      dimension lmmua( 2, * ),    ! array of packed (lm,mu) of SU(3)    ACVI8612
     1          ispta( * ),       ! array of packed (2s1,2s2) of SU(2)  ACVI8613
     2          iptr( * )         ! level pointers for generations      ACVI8614
C               mxnumb            ! max # (lm mu) pairs lmmua           ACVI8615
C     Packing functions                                                 ACVI8616
      data nbit8 / 255 /           ! zff                                ACVI8617
      ipack( i, j, k, l ) = ior(l,ishft(ior(k,ishft(ior(j,ishft(i,      ACVI8618
     1       8)),8)),8))                                                ACVI8619
      iupac( index, iover, ibits) = iand( ishft( index, iover), ibits ) ACVI8620
C                                                                       ACVI8621
C     Transfer lmopmx and muopmx                                        ACVI8622
      lmopmn  = lmmua( 1, mxnumb-1 )                                    ACVI8623
      muopmn  = lmmua( 2, mxnumb-1 )                                    ACVI8624
      lmopmx  = lmmua( 1, mxnumb )                                      ACVI8625
      muopmx  = lmmua( 2, mxnumb )                                      ACVI8626
      itensor = iptr ( mxnumb-1 )                                       ACVI8627
      iequal  = iptr ( mxnumb )                                         ACVI8628
C                                                                       ACVI8629
      if( itensor.eq.0 .or. itensor.eq.1 ) then                         ACVI8630
        tenchk  = .false.                                               ACVI8631
      else                                                              ACVI8632
        tenchk  = .true.                                                ACVI8633
      end if                                                            ACVI8634
      n1a = iabs( n1 )                                                  ACVI8635
      n2a = iabs( n2 )                                                  ACVI8636
      nlo = min0( n1a, n2a )                                            ACVI8637
      nhi = max0( n1a, n2a )                                            ACVI8638
C                                                                       ACVI8639
      if( iequal .gt. 0 ) then                                          ACVI8640
        itemp  = max0( lmopmn, muopmn )                                 ACVI8641
C       lmopmn = itemp                                                  ACVI8642
        muopmn = itemp                                                  ACVI8643
        itemp  = min0( lmopmx, muopmx )                                 ACVI8644
        lmopmx = itemp                                                  ACVI8645
        muopmx = itemp                                                  ACVI8646
      end if                                                            ACVI8647
C                                                                       ACVI8648
C     Calculate some t1 x t2 = T(lm mu), where                          ACVI8649
C       t1 = a+a+..., t2 = a a ..., and t1 >= t2 if n2 < 0              ACVI8650
C       t1 = a+a+..., t2 = a+a+...               if n2 >= 0             ACVI8651
C     If t1 < t2 take the conjugate                                     ACVI8652
C                                                                       ACVI8653
      ibot = iptr( nhi )                                                ACVI8654
      itop = iptr( nhi+1 ) - 1                                          ACVI8655
      if( nlo .gt. 0 ) then                                             ACVI8656
        k = itop                                                        ACVI8657
        jbot = iptr( nlo )                                              ACVI8658
        jtop = iptr( nlo+1 ) - 1                                        ACVI8659
C                                                                       ACVI8660
        do i = ibot, itop                                               ACVI8661
          index = lmmua( 2, i )                                         ACVI8662
          lm1   = iupac( index, -16, nbit8 )                            ACVI8663
          mu1   = iupac( index,  -8, nbit8 )                            ACVI8664
          index = ispta( i )                                            ACVI8665
          ix1   = iupac( index, -24, nbit8 )                            ACVI8666
          ist1  = iupac( index,   0, nbit8 )                            ACVI8667
          do j = jbot, jtop                                             ACVI8668
            index = lmmua( 2, j )                                       ACVI8669
            if( n2 .lt. 0 ) then                                        ACVI8670
              mu2 = iupac( index, -16, nbit8 )                          ACVI8671
              lm2 = iupac( index,  -8, nbit8 )                          ACVI8672
            else                                                        ACVI8673
              lm2 = iupac( index, -16, nbit8 )                          ACVI8674
              mu2 = iupac( index,  -8, nbit8 )                          ACVI8675
            end if                                                      ACVI8676
            index = ispta( j )                                          ACVI8677
            ix2   = iupac( index, -24, nbit8 )                          ACVI8678
            ist2  = iupac( index,   0, nbit8 )                          ACVI8679
C                                                                       ACVI8680
            lam   = lm1 + mu1 + lm2 + mu2                               ACVI8681
            lmmax = min0( lmopmx, lam )                                 ACVI8682
            mumax = min0( muopmx, lam )                                 ACVI8683
            if( iequal .gt. 0 ) lmopmn = lmmax                          ACVI8684
            do ll = lmopmn, lmmax                                       ACVI8685
              do mu = muopmn, mumax                                     ACVI8686
                if( iequal .gt. 0 ) then                                ACVI8687
                  lm = mu                                               ACVI8688
                else                                                    ACVI8689
                  lm = ll                                               ACVI8690
                end if                                                  ACVI8691
                itest = iand( ior(lm,mu), 1 )                           ACVI8692
                if( itest.eq.itensor .or. tenchk ) then                 ACVI8693
                  call u3mult( lm1, mu1, lm2, mu2, lm, mu, kro, *100 )  ACVI8694
                  k = k + 1                                             ACVI8695
                  if(k.gt.mxnumb-1) call error(' SU3TLM: Too many ops.')ACVI8696
                  lmmua( 1, k ) = ipack( lm1, mu1, lm2, mu2 )           ACVI8697
                  lmmua( 2, k ) = ipack( kro,  lm,  mu,   0 )           ACVI8698
                  ispta( k )    = ipack( ix1,ist1, ix2,ist2 )           ACVI8699
100               continue                                              ACVI8700
                end if                                                  ACVI8701
              end do ! mu                                               ACVI8702
            end do ! ll                                                 ACVI8703
          end do ! j : a indices                                        ACVI8704
        end do ! i : a+ indices                                         ACVI8705
        ibot = itop + 1                                                 ACVI8706
        itop = k                                                        ACVI8707
      end if ! nlo check                                                ACVI8708
C                                                                       ACVI8709
C     lmmua( 1, mxnumb ) = ibot                                         ACVI8710
C     lmmua( 2, mxnumb ) = itop                                         ACVI8711
C                                                                       ACVI8712
C     Compress array lmmua                                              ACVI8713
      k = 0                                                             ACVI8714
      do i = ibot, itop                                                 ACVI8715
        k = k + 1                                                       ACVI8716
        lmmua( 1, k ) = lmmua( 1, i )                                   ACVI8717
        lmmua( 2, k ) = lmmua( 2, i )                                   ACVI8718
        ispta( k )    = ispta( i )                                      ACVI8719
      end do                                                            ACVI8720
      do i = k + 1, itop                                                ACVI8721
        lmmua( 1, i ) = 0                                               ACVI8722
        lmmua( 2, i ) = 0                                               ACVI8723
        ispta( i )    = 0                                               ACVI8724
      end do                                                            ACVI8725
C     Auxilliary indices                                                ACVI8726
      lmmua( 1, mxnumb-1 ) = ipack( 0, 0, n1a, n2a )                    ACVI8727
      lmmua( 2, mxnumb-1 ) = k                                          ACVI8728
      lmmua( 1, mxnumb ) = ipack( 0, 0, nlo, nhi )                      ACVI8729
      lmmua( 2, mxnumb ) = k                                            ACVI8730
999   return                                                            ACVI8731
C ---*end of su3tlm*----------------------------------------------------ACVI8732
      end                                                               ACVI8733
CB----------------------------------------------------------------------ACVI8734
C                                                                       ACVI8735
C                          **************                               ACVI8736
C                          *** SETRME ***                               ACVI8737
C                          **************                               ACVI8738
C                                                                       ACVI8739
C ----------------------------------------------------------------------ACVI8740
C Author:  Chairul Bahri                                                ACVI8741
C ----------------------------------------------------------------------ACVI8742
C Updates: 11/90 ==> Original                                           ACVI8743
C          12/13/90: test for a+ (fp)**8 & conjugate ... SUCCESSFUL !   ACVI8744
C          12/25/90: using uncoupled tensor T(LM LM)=a+ a               ACVI8745
C          01/13/91: test for a+a(ds)**2 ... SUCCESSFUL for RO<>1, too !ACVI8746
C                                (4 0) (0 2) (2 1)                      ACVI8747
C          02/25/91: multi irreps.                                      ACVI8748
C          01/20/92: Sun 4.                                             ACVI8749
C          01/07/93: RS/6000                                            ACVI8750
C ----------------------------------------------------------------------ACVI8751
C                                                                       ACVI8752
C    The subroutine SETRME setups the parameters for calculating the    ACVI8753
C    the Reduced Matrix Elements and/or the Coefficient of Fractional   ACVI8754
C    Parentages.                                                        ACVI8755
C                                                                       ACVI8756
C ----------------------------------------------------------------------ACVI8757
      subroutine setrme( nmax, iqops, mxopb, mxjtck )                   ACVI8758
      implicit real*8(d), logical(t)                                    ACVI8759
      parameter ( MXBITS = 32 )                                         ACVI8760
      dimension iqops( 2, * )      ! packed qn's for operators          ACVI8761
      parameter( MXJT = 1000 )     ! max # SU(3)>SU(2)xU(1) operators   ACVI8762
      common / LOGCON / tfile,     ! write to logfile (full)            ACVI8763
     1                  tsmfi,     ! write to logfile (small)           ACVI8764
     2                  tconj,     ! conjugate SU(3) label              ACVI8765
     3                  topbase    ! calculate SU(3) coupled opb        ACVI8766
      common / RMECON / lfile,     ! log file for intermediate results  ACVI8767
     1                  ieta       ! oscillator shell number            ACVI8768
      common / RMEBIT / nsize,     ! # of levels for each word          ACVI8769
     1                  nwords,    ! # of words for each bit state      ACVI8770
     2                  kact,      ! # of active words for bit man.     ACVI8771
     3                  nbits,     ! bit length for # levels in shell   ACVI8772
     4                  nbitz,     ! significant bits                   ACVI8773
     5                  nzxy( 66 ),! (nz,nx,ny) packed labels for 1-prt.ACVI8774
     6                  nmlvls( 3 )! pointers for next words            ACVI8775
      common / RMEIND / ilm,       ! # operators                        ACVI8776
     1                  i3ix,      ! last index for overlaps            ACVI8777
     2                  npmax,     ! # a+                               ACVI8778
     3                  nhmax,     ! # a                                ACVI8779
     4                  npop1      ! (npmax,nhmax)                      ACVI8780
      common / OPBSU2 / dwigsu2(9),! SU(2) coupled operator bases       ACVI8781
     1                  ismin,     ! 2 Smin                             ACVI8782
     2                  ismax      ! 2 Smax                             ACVI8783
      common / SU3U21 /                                                 ACVI8784
     1         dxwu3( 9*MXJT ),    ! Wigner coefficient matrices        ACVI8785
     2         ia( MXJT ),         ! indices for LU decomposition       ACVI8786
     3         lmjt( MXJT ),       ! packed ro(lm mu)jt labels          ACVI8787
     4         i3ptr( MXJT ),      ! pointers to lmjt(*)                ACVI8788
     5         ieop, mlmop         ! projection q.n.                    ACVI8789
      common / WORKAR / dwork(MXJT)! working arrays                     ACVI8790
      dimension mst( 9 )           ! SU(2) projected tensors            ACVI8791
C                                                                       ACVI8792
C     Packing and unpacking function.                                   ACVI8793
      data nbit8 / 255 /           ! zff                                ACVI8794
      ipack( i, j, k, l, m ) = ior(m, ishft( ior(l, ishft( ior(k,       ACVI8795
     1       ishft( ior(j, ishft( i,4 )),8 )),8 )),8 ))                 ACVI8796
      iupac( index, iover, ibits ) = iand( ishft( index,iover ), ibits )ACVI8797
C                                                                       ACVI8798
C     Initialize operator quantum numbers.                              ACVI8799
      npop   = iqops( 1, mxopb-1 )                                      ACVI8800
      npop1  = iqops( 1, mxopb )                                        ACVI8801
      ilm    = iqops( 2, mxopb )                                        ACVI8802
      nsize  = MXBITS                                                   ACVI8803
      if( ilm .gt. MXJT ) then                                          ACVI8804
        write(6,'(a,i5,a)') 'Sorry, we handle up to',MXJT,' operators.' ACVI8805
        ilm = MXJT                                                      ACVI8806
      end if                                                            ACVI8807
      mxjtck = MXJT                                                     ACVI8808
      tfile  = lfile .gt. 6                                             ACVI8809
      tsmfi  = lfile .gt. 0                                             ACVI8810
      i3ix   = 9 * ilm                                                  ACVI8811
      n      = (ieta + 1) * (ieta + 2) / 2                              ACVI8812
      nlevel = n                                                        ACVI8813
      nwords = (nlevel + nsize - 1) / nsize                             ACVI8814
C     Single particle quanta distributions                              ACVI8815
      i      = 0                                                        ACVI8816
      do nz = ieta, 0, -1                                               ACVI8817
        do nx = ieta-nz, 0, -1                                          ACVI8818
          ny        = ieta - (nz + nx)                                  ACVI8819
          i         = i + 1                                             ACVI8820
          nzxy( i ) = ipack( 0, 0, nz, nx, ny )                         ACVI8821
        end do                                                          ACVI8822
      end do                                                            ACVI8823
C     Pointers for applications of SU(M) actions.                       ACVI8824
      do i = 1, nwords                                                  ACVI8825
        nfinal      = min0( n, nsize )                                  ACVI8826
        nmlvls( i ) = nfinal - 1                                        ACVI8827
        n           = n - nfinal                                        ACVI8828
      end do                                                            ACVI8829
C                                                                       ACVI8830
      do i = 0, nsize-1                                                 ACVI8831
        if( btest( nlevel, i ) ) nbits = i + 1                          ACVI8832
      end do                                                            ACVI8833
      if((nmax*nbits).gt.nsize ) call error(' SETRME: too many conf.')  ACVI8834
      nbitz = ishft( 1, nbits ) - 1                                     ACVI8835
      nz    = mod( nlevel, nsize )                                      ACVI8836
      kact  = nmax * nwords                                             ACVI8837
      kactm = kact - nmax                                               ACVI8838
      do i = 1, ilm                                                     ACVI8839
        i3ptr( i ) = 0                                                  ACVI8840
      end do                                                            ACVI8841
C                                                                       ACVI8842
      index = iqops( 1, mxopb-1 )                                       ACVI8843
      npmax = iupac( index, -8, nbit8 )                                 ACVI8844
      nhmax = iupac( index,  0, nbit8 )                                 ACVI8845
      tconj = npop1.eq.npop .and. npmax.ne.nhmax                        ACVI8846
C     Transfer quantum numbers for operators.                           ACVI8847
C                                                                       ACVI8848
      if( tfile ) then                                                  ACVI8849
        write(lfile,'(/a/)') ' *** SU(3) operators ***'                 ACVI8850
        do i = 1, ilm                                                   ACVI8851
          index = iqops( 1, i )                                         ACVI8852
          lm1   = iupac( index, -24, nbit8 )                            ACVI8853
          mu1   = iupac( index, -16, nbit8 )                            ACVI8854
          lm2   = iupac( index,  -8, nbit8 )                            ACVI8855
          mu2   = iupac( index,   0, nbit8 )                            ACVI8856
          index = iqops( 2, i )                                         ACVI8857
          kro   = iupac( index, -24, nbit8 )                            ACVI8858
          lm    = iupac( index, -16, nbit8 )                            ACVI8859
          mu    = iupac( index,  -8, nbit8 )                            ACVI8860
          write(lfile,'(2(a,2i3),a,i3,a,2i3,a)')                        ACVI8861
     1      '      (',lm1,mu1,')x(',lm2,mu2,') =',kro,'(',lm,mu,')'     ACVI8862
          if( mod(i,5) .eq. 0 ) write(lfile,*)                          ACVI8863
        end do                                                          ACVI8864
        write(lfile,*)                                                  ACVI8865
      end if                                                            ACVI8866
C                                                                       ACVI8867
C     SU(2) coupled operator bases                                      ACVI8868
      is1min = iand( npmax, 1 )                                         ACVI8869
      is1max = npmax                                                    ACVI8870
      is2min = iand( nhmax, 1 )                                         ACVI8871
      is2max = nhmax                                                    ACVI8872
C                                                                       ACVI8873
      nhi = max0( npmax, nhmax )                                        ACVI8874
      if( nhi .gt. 4 ) call error(' SETRME: Np>4 not available.')       ACVI8875
      if( nhi .eq. 1 ) then                                             ACVI8876
        mst( 1 ) = -1                                                   ACVI8877
        mst( 3 ) = 1                                                    ACVI8878
        return                                                          ACVI8879
      end if                                                            ACVI8880
C                                                                       ACVI8881
      do i = 1, 9                                                       ACVI8882
        dwork( i ) = 0.d0                                               ACVI8883
      end do                                                            ACVI8884
C                                                                       ACVI8885
C     SU(2) Wigner coefficients for operator bases                      ACVI8886
      dwigsu2( 1 ) = dwr3( 1,1,2, -1,-1,-2 )   ! 00 2                   ACVI8887
      dwigsu2( 3 ) = dwr3( 1,1,2, -1, 1, 0 )   ! 01 2                   ACVI8888
      dwigsu2( 5 ) = dwr3( 1,1,2,  1, 1, 2 )   ! 11 2                   ACVI8889
      dwigsu2( 2 ) = dwr3( 1,1,0, -1, 1, 0 )   ! 01 0                   ACVI8890
      mst( 1 ) = -2                                                     ACVI8891
      mst( 3 ) =  0                                                     ACVI8892
      mst( 5 ) =  2                                                     ACVI8893
      mst( 2 ) =  0                                                     ACVI8894
      index = 5                                                         ACVI8895
      do i = 1, index                                                   ACVI8896
        dwork( i ) = dwigsu2( i )                                       ACVI8897
      end do                                                            ACVI8898
      ibits = 1                                                         ACVI8899
      do i = 3, nhi                                                     ACVI8900
        if( i .eq. 4 ) ibits = 2                                        ACVI8901
        do j = 1, index                                                 ACVI8902
          dwigsu2( j ) = 0.d0                                           ACVI8903
        end do                                                          ACVI8904
        is3min = iand( i, 1 )                                           ACVI8905
        do is3t = is3min, i, 2                                          ACVI8906
          k  = is3t - 1                                                 ACVI8907
          k2 = k / 2                                                    ACVI8908
          l  = is3t / 2                                                 ACVI8909
          do j = 0, i-1                                                 ACVI8910
            ix2 = ior( ishft(j,1), k )                                  ACVI8911
            if( ix2 .gt. 0 ) then                                       ACVI8912
              m2t = mst( k2 )                                           ACVI8913
              ix3 = ior( ishft(j,   ibits), l )                         ACVI8914
              ix4 = ior( ishft(j+1, ibits), l )                         ACVI8915
              mst( ix3 )     = m2t - 1                                  ACVI8916
              dwigsu2( ix3 ) = dwigsu2( ix3 ) + dwork( ix2 ) *          ACVI8917
     1                         dwr3( 1,k,is3t, -1,m2t,mst(ix3) )        ACVI8918
              mst( ix4 )     = m2t + 1                                  ACVI8919
              dwigsu2( ix4 ) = dwigsu2( ix4 ) + dwork( ix2 ) *          ACVI8920
     1                         dwr3( 1,k,is3t,  1,m2t,mst(ix3) )        ACVI8921
            end if                                                      ACVI8922
          end do                                                        ACVI8923
          index = ix4                                                   ACVI8924
          do j = 1, index                                               ACVI8925
            dwork( j ) = dwigsu2( j )                                   ACVI8926
          end do                                                        ACVI8927
        end do                                                          ACVI8928
      end do                                                            ACVI8929
      if( tfile ) then                                                  ACVI8930
        write(lfile,'(/a/)') ' SU(2)-op Wigner (CG) coefficients.'      ACVI8931
        j    = iand( nhi, 1 )                                           ACVI8932
        do i = 1, index - 2                                             ACVI8933
          ist = iand( i, ibits ) * 2 + j                                ACVI8934
          write(lfile,'(f7.3,i5,i3)') dwigsu2(i), ist, mst(i)           ACVI8935
        end do                                                          ACVI8936
          ist = iand( index, ibits ) * 2 + j                            ACVI8937
          write(lfile,'(f7.3,i5,i3)') dwigsu2(index), ist, mst(index)   ACVI8938
      end if                                                            ACVI8939
      return                                                            ACVI8940
C ---*end of setrme*----------------------------------------------------ACVI8941
      end                                                               ACVI8942
      subroutine qnsu2( iqnums, iqwrit , KS1,KS2)                       ACVI8943
      implicit real*8(d)                                                ACVI8944
      dimension iqnums( * ),                                            ACVI8945
     1          iqwrit( * )                                             ACVI8946
      common / OPBSU2 / dwigsu2(9),! SU(2) coupled operator bases       ACVI8947
     1                  ismin,     ! 2 Smin                             ACVI8948
     2                  ismax      ! 2 Smax                             ACVI8949
C     bras                                                              ACVI8950
      iqwrit( 1 ) = iqnums( 3 ) + iqnums( 4 )                           ACVI8951
      iqwrit( 3 ) = iqnums( 1 )                                         ACVI8952
      iqwrit( 4 ) = iqnums( 2 )                                         ACVI8953
      iqwrit( 5 ) = iqnums( 3 ) - iqnums( 4 )                           ACVI8954
      islt  = iqwrit( 5 )                                               ACVI8955
C     kets                                                              ACVI8956
      iqwrit( 18 ) = iqnums( 7 ) + iqnums( 8 )                          ACVI8957
      iqwrit( 20 ) = iqnums( 5 )                                        ACVI8958
      iqwrit( 21 ) = iqnums( 6 )                                        ACVI8959
      iqwrit( 22 ) = iqnums( 7 ) - iqnums( 8 )                          ACVI8960
      isrt  = iqwrit( 22 )                                              ACVI8961
      ismin = iabs( islt - isrt )                                       ACVI8962
      ismax = islt + isrt                                               ACVI8963
C                                                                       ACVI8964
      return                                                            ACVI8965
      end                                                               ACVI8966
C---------------------------------------------------------------------- ACVI8967
C                                                                       ACVI8968
C                           **************                              ACVI8969
C                           *** OPTLM2 ***                              ACVI8970
C                           **************                              ACVI8971
C                                                                       ACVI8972
C ----------------------------------------------------------------------ACVI8973
C Author:  Chairul Bahri (LSU)                                          ACVI8974
C References:  -*-                                                      ACVI8975
C ----------------------------------------------------------------------ACVI8976
C Updates: 11/90 ==> original (IBM 3090)                                ACVI8977
C          11/28/90: implemented in RMECODE.                            ACVI8978
C          04/09/90: phase included.                                    ACVI8979
C          11/17/92: IBM RS/6000-560.                                   ACVI8980
C ----------------------------------------------------------------------ACVI8981
C                                                                       ACVI8982
C    The subroutine OPTLM generates the matrix elements of coupled SU(3)ACVI8983
C    tensor operators t1xt2 (lm mu) between various SU(3) basis, recur- ACVI8984
C    sively.                                                            ACVI8985
C                                                                       ACVI8986
C ----------------------------------------------------------------------ACVI8987
C                                                                       ACVI8988
      subroutine optlm2( label,lctree, dwigsu3,dwsu3,mxwcof, n1,n2a, * )ACVI8989
      implicit real*8(d), logical(t)                                    ACVI8990
CD    character*4 cflag, cphase                                         ACVI8991
      parameter ( DZERO = 1.d-12 )                                      ACVI8992
C                                                                       ACVI8993
      dimension label( * ),        ! labels for trees of operators      ACVI8994
     1          lctree(-10:* ),    ! binary tree for uncoupled tensors  ACVI8995
     2          dwigsu3( * ),      ! coefficients for coupled SU(3)     ACVI8996
     3          dwsu3( * )         ! temporary coeff. for coupled SU(3) ACVI8997
C               mxwcof             ! max # operator basis               ACVI8998
C               n1                 ! # a+ in t1                         ACVI8999
C               n2a                ! # a+ (or a) in t2                  ACVI9000
      common / LOGCON / tfile,     ! yes, write to log file (full)      ACVI9001
     1                  tsmfi,     ! write to log file (small)          ACVI9002
     2                  tconj,     ! conjugate SU(3) label              ACVI9003
     3                  topbase    ! calculate SU(3) coupled opb        ACVI9004
      common / RMECON / lfile,     ! logfile                            ACVI9005
     1                  ieta       ! oscillator shell number            ACVI9006
      common / RMEBIT / nsize,     ! # of levels for each word          ACVI9007
     1                  nwords,    ! # words in a bit state             ACVI9008
     2                  kact,      ! # active words for bit manipulat.  ACVI9009
     3                  nbits,     ! bit length for # levels in shell   ACVI9010
     4                  nbitz,     ! significant bits                   ACVI9011
     5                  nzxy( 66 ),! (nz,nx,ny) packed labels for 1-prt.ACVI9012
     6                  nmlvls( 3 )! pointers for words                 ACVI9013
CD                      cflag,     ! flag for new SU(3) irreps          ACVI9014
CD                      cphase,    ! SU(3) phase                        ACVI9015
      character*4 cann             ! annihilation operator              ACVI9016
C                                                                       ACVI9017
C     SU(3) package                                                     ACVI9018
C                                                                       ACVI9019
      parameter ( NCE1 = 9, NCE2 = 13244, NCEX = 42 )                   ACVI9020
      parameter ( NCW1 = 9, NCW2 = 42, NCW3 = 9030 )                    ACVI9021
      parameter ( KIMAX1 = 3*NCE2, KIMAX2 = 3*NCW3 )                    ACVI9022
      parameter ( NCW22 = NCW2*NCW2 )                                   ACVI9023
C     common / SU3EXT / tj2ta( NCE2 )                                   ACVI9024
      common / SU3EXI / j1ta( NCE2 ), j2ta( NCE2 ), iea( NCE2 )         ACVI9025
      common / SU3I   / j1smax( NCW22 ), j1tmax( NCW22 ),               ACVI9026
     1                  j2smax( NCW2 ), j2tmax( NCW2 ),  indmat( NCW22 )ACVI9027
      common / SU3EXD / dewu3( KIMAX1 )                                 ACVI9028
      common / SU3D   / dwu3( KIMAX2 )                                  ACVI9029
C                                                                       ACVI9030
C     Packing functions                                                 ACVI9031
      data nbit4, nbit8 / 16, 255 /! zf, zff                            ACVI9032
      ipack( i, j, k, l ) = ior( l, ishft( ior( k, ishft( ior( j,       ACVI9033
     1        ishft( i, 8 )), 8 )), 8 ))                                ACVI9034
      ipack5( i, j, k, l, m ) = ior( m, ishft( ior( l, ishft( ior( k,   ACVI9035
     1        ishft( ior( j, ishft( i, 4 )), 8 )), 8 )), 8 ))           ACVI9036
      iupac(index,iover,ibits) = iand( ishft( index, iover), ibits )    ACVI9037
C                                                                       ACVI9038
C t(lm1 mu1) x t(lm2 mu2)                                               ACVI9039
C n1           n2                                                       ACVI9040
C              n2a<0 for annihilation operators.                        ACVI9041
C     Transfer labels                                                   ACVI9042
      n2     = iabs( n2a )                                              ACVI9043
      ntot   = n1 + n2                                                  ACVI9044
      n1eta  = n1 * ieta                                                ACVI9045
      n2eta  = n2 * ieta                                                ACVI9046
      index  = label( 1 )                                               ACVI9047
        label1 = index                                                  ACVI9048
        lm1    = iupac( index, -24, nbit8 )                             ACVI9049
        mu1    = iupac( index, -16, nbit8 )                             ACVI9050
        lm2    = iupac( index,  -8, nbit8 )                             ACVI9051
        mu2    = iupac( index,   0, nbit8 )                             ACVI9052
      index  = label( 2 )                                               ACVI9053
        label2 = index                                                  ACVI9054
        kro12m = iupac( index, -24, nbit8 )                             ACVI9055
        lm     = iupac( index, -16, nbit8 )                             ACVI9056
        mu     = iupac( index,  -8, nbit8 )                             ACVI9057
        jt     = iupac( index,   0, nbit8 )                             ACVI9058
      lvls   = label( 3 )                                               ACVI9059
        nz1  = 0                                                        ACVI9060
        nx1  = 0                                                        ACVI9061
        ny1  = 0                                                        ACVI9062
        do i = 1, n1                                                    ACVI9063
          level = iand( lvls, nbitz )                                   ACVI9064
          index = nzxy( level )                                         ACVI9065
          nz1   = nz1 + iupac( index, -16, nbit8 )                      ACVI9066
          nx1   = nx1 + iupac( index,  -8, nbit8 )                      ACVI9067
          ny1   = ny1 + iupac( index,   0, nbit8 )                      ACVI9068
          lvls  = ishft( lvls, -nbits )                                 ACVI9069
        end do                                                          ACVI9070
      lvls   = label( 4 )                                               ACVI9071
        nz2  = 0                                                        ACVI9072
        nx2  = 0                                                        ACVI9073
        ny2  = 0                                                        ACVI9074
        do i = 1, n2                                                    ACVI9075
          level = iand( lvls, nbitz )                                   ACVI9076
          index = nzxy( level )                                         ACVI9077
          nz2   = nz2 + iupac( index, -16, nbit8 )                      ACVI9078
          nx2   = nx2 + iupac( index,  -8, nbit8 )                      ACVI9079
          ny2   = ny2 + iupac( index,   0, nbit8 )                      ACVI9080
          lvls  = ishft( lvls, -nbits )                                 ACVI9081
        end do                                                          ACVI9082
      ie1 = 3 * nz1 - n1eta                                             ACVI9083
      ie2 = 3 * nz2 - n2eta                                             ACVI9084
      m1t = nx1 - ny1                                                   ACVI9085
      m2t = nx2 - ny2                                                   ACVI9086
      m1tabs = iabs( m1t )                                              ACVI9087
      m2tabs = iabs( m2t )                                              ACVI9088
      iphase = 0                                                        ACVI9089
      if( n2a .lt. 0 ) then                                             ACVI9090
        ie2    = -ie2                                                   ACVI9091
        m2t    = -m2t                                                   ACVI9092
        iphase = ( 2 * (mu2 - lm2) + ie2 - 3 * m2t ) / 6                ACVI9093
      end if                                                            ACVI9094
      ie  = ie1 + ie2                                                   ACVI9095
      mt  = m1t + m2t                                                   ACVI9096
C                                                                       ACVI9097
C     Check allowed chains                                              ACVI9098
      if( ie1.lt.-(lm1+2*mu1) .or. ie1.gt.(2*lm1+mu1) ) return 1        ACVI9099
      if( ie2.lt.-(lm2+2*mu2) .or. ie2.gt.(2*lm2+mu2) ) return 1        ACVI9100
      if( m1tabs.gt.lm1+mu1 ) return 1                                  ACVI9101
      if( m2tabs.gt.lm2+mu2 ) return 1                                  ACVI9102
C     Reset dwsu3                                                       ACVI9103
      do kro = 1, kro12m                                                ACVI9104
        dwsu3( kro ) = 0.d0                                             ACVI9105
      end do                                                            ACVI9106
CW----------------------------------------------------------------------ACVI9107
CW    if( tfile ) write(lfile,'(a)') ' Entering OPTLM2:'                ACVI9108
CW----------------------------------------------------------------------ACVI9109
C                                                                       ACVI9110
C         Operator COUPLINGS                                            ACVI9111
C                   <   t1   ;   t2   |   T  > v for all isoscalars     ACVI9112
          call xewu3( lm1,mu1, lm2,mu2, lm,mu, 1, NEC,                  ACVI9113
     1      kromax,indmax, dewu3,j1ta,j2ta,iea, NCE1,NCE2, KIMAX1 )     ACVI9114
          call xwu3( lm1,mu1, lm2,mu2, lm,mu, ie,jt, NEC, dewu3,        ACVI9115
     1      kromax,indmax, dwu3,j1smax,j1tmax,j2smax,j2tmax,            ACVI9116
     2      iesmax,ie2max, indmat, NCW1,NCW2,NCW3, KIMAX2 )             ACVI9117
          do ies = 1, iesmax                                            ACVI9118
            if( ie2 .eq. ie2max-3*(iesmax-ies) ) then                   ACVI9119
              do 150 j2s = 1, j2smax( ies )                             ACVI9120
                j2t = j2tmax(ies) - 2*(j2s-1)                           ACVI9121
                if( j2t .ge. m2tabs ) then                              ACVI9122
                  iesj2s = ies + NCW2*(j2s-1)                           ACVI9123
                  if( n2a.ge.0 ) then                                   ACVI9124
                    label( 1 ) = ipack( 1, lm2, mu2, j2t )              ACVI9125
                  else                                                  ACVI9126
                    label( 1 ) = ipack( 1, mu2, lm2, j2t )              ACVI9127
                  end if                                                ACVI9128
                  label( 2 ) = label( 4 )                               ACVI9129
C                 Load t2 from lctree                                   ACVI9130
                  call tchk( label, lctree, *120 )                      ACVI9131
CW                  write(lfile,*) 't2 not found'                       ACVI9132
                    go to 150                                           ACVI9133
120               iloc  = lctree(-5 ) + lctree(-4 ) + 1                 ACVI9134
                  dwig2 = dwigsu3( lctree(iloc) )                       ACVI9135
                  if( btest(iphase,0) ) dwig2 = -dwig2                  ACVI9136
                  do 140 j1s = 1, j1smax( iesj2s )                      ACVI9137
                    j1t = j1tmax(iesj2s) - 2*(j1s-1)                    ACVI9138
                    if( j1t .ge. m1tabs ) then                          ACVI9139
                      ind = (indmat(iesj2s) - j1t) / 2                  ACVI9140
                      label( 1 ) = ipack( 1, lm1, mu1, j1t )            ACVI9141
                      label( 2 ) = label( 3 )                           ACVI9142
C                     Load t1 from lctree                               ACVI9143
                      call tchk( label, lctree, *130 )                  ACVI9144
CW                      write(lfile,*) 't1 not found'                   ACVI9145
                        go to 140                                       ACVI9146
130                   iloc  = lctree(-5 ) + lctree(-4 ) + 1             ACVI9147
                      dwig1 = dwigsu3( lctree(iloc) )                   ACVI9148
                      dwig1 = dwig1*dwig2*dwr3( j1t,j2t,jt, m1t,m2t,mt )ACVI9149
C                                                                       ACVI9150
C                     Couple the tensors                                ACVI9151
                      indro = kromax * (ind - 1)                        ACVI9152
                      do kro = 1, kromax                                ACVI9153
                        indro = indro + 1                               ACVI9154
                        dwsu3( kro ) = dwsu3( kro ) + dwig1*dwu3(indro) ACVI9155
                      end do ! kro                                      ACVI9156
                    end if                                              ACVI9157
140               continue ! j1s                                        ACVI9158
                end if                                                  ACVI9159
150           continue ! j2s                                            ACVI9160
            end if                                                      ACVI9161
          end do ! ies                                                  ACVI9162
C     Reset label                                                       ACVI9163
      label( 1 ) = label1                                               ACVI9164
      label( 2 ) = label2                                               ACVI9165
      return                                                            ACVI9166
CS999 write(6,*) ' OPTLM2: Tree overloaded.'                            ACVI9167
CS    stop                                                              ACVI9168
C ---*end of optlm2*----------------------------------------------------ACVI9169
      end                                                               ACVI9170
CB----------------------------------------------------------------------ACVI9171
C                                                                       ACVI9172
C                          **************                               ACVI9173
C                          *** WIGMAT ***                               ACVI9174
C                          **************                               ACVI9175
C                                                                       ACVI9176
C ----------------------------------------------------------------------ACVI9177
C Author:  Chairul Bahri                                                ACVI9178
C ----------------------------------------------------------------------ACVI9179
C Updates: 11/90 ==> Original ... (in SETRME)                           ACVI9180
C          02/25/91: multi irreps.                                      ACVI9181
C          06/01/91: separated from SETRME.                             ACVI9182
C          06/18/91: a+ and a modified.                                 ACVI9183
C          08/17/91: WST and OPBASIS implemented                        ACVI9184
C ----------------------------------------------------------------------ACVI9185
C                                                                       ACVI9186
C    The subroutine WIGMAT calculates SU(3) Wigner coefficient matrices ACVI9187
C    being used to extract the SU(3) reduced matrix elements from       ACVI9188
C    the matrix elementsthe matrix elements.                            ACVI9189
C                                                                       ACVI9190
C    The present calculation is for the triple barred                   ACVI9191
C       SU(3) x SU(2) only.                                             ACVI9192
C --------------------------------------------------------------------- ACVI9193
      subroutine wigmat( iqops, iqnums, i3max, mxjtck, nmax, * )        ACVI9194
      implicit real*8(d), logical(t)                                    ACVI9195
      dimension iqops( 2, * ),     ! packed qn's for operators          ACVI9196
     1          iqnums( * )        ! quantum numbers for states         ACVI9197
C               i3max              ! # matrix elements                  ACVI9198
C               mxjtck             ! check MXJT                         ACVI9199
      parameter( MXJT = 1000 )     ! max # SU(3)>SU(2)xU(1) operators   ACVI9200
      common / SU3U21 /                                                 ACVI9201
     1         dxwu3( 9*MXJT ),    ! Wigner coefficient matrices        ACVI9202
     2         ia( MXJT ),         ! indices for LU decomposition       ACVI9203
     3         lmjt( MXJT ),       ! packed ro(lm mu)jt labels for ops. ACVI9204
     4         i3ptr( MXJT ),      ! pointers to lmjt(*)                ACVI9205
     5         ieop, mlmop         ! SU(3)?                             ACVI9206
      common / LOGCON / tfile,     ! write to logfile (full)            ACVI9207
     1                  tsmfi,     ! write to logfile (small)           ACVI9208
     2                  tconj,     ! conjugate SU(3) label              ACVI9209
     3                  topbase    ! calculate SU(3) coupled opb        ACVI9210
      common / RMECON / lfile,     ! log file for intermediate results  ACVI9211
     1                  ieta       ! oscillator shell number            ACVI9212
      common / RMEIND / ilm,       ! # operators                        ACVI9213
     1                  i3ix,      ! last index for overlaps            ACVI9214
     2                  npmax,     ! # a+                               ACVI9215
     3                  nhmax,     ! # a                                ACVI9216
     4                  npop1      ! (npmax,nhmax)                      ACVI9217
C                                                                       ACVI9218
C     SU(3) package                                                     ACVI9219
C                                                                       ACVI9220
      parameter ( NCE1 = 9, NCE2 = 13244, NCEX = 42 )                   ACVI9221
      parameter ( KIMAX1 = 3*NCE2 )                                     ACVI9222
      common / SU3EXT / tj2ta( NCE2 )                                   ACVI9223
      common / SU3EXI / j1ta( NCE2 ), j2ta( NCE2 ), iea( NCE2 )         ACVI9224
      common / SU3EXD / dewu3( KIMAX1 )                                 ACVI9225
C                                                                       ACVI9226
C     Unpacking function.                                               ACVI9227
      data nbit8 / 255 /           ! zff                                ACVI9228
      iupac( index, iover, ibits ) = iand( ishft( index,iover ), ibits )ACVI9229
C                                                                       ACVI9230
      if( mxjtck .ne. MXJT ) call error(' WIGMAT: Check MXJT!')         ACVI9231
C                                                                       ACVI9232
C     Transfer SU(3) quantum numbers for operators.                     ACVI9233
C     Transfer SU(3) quantum numbers for states.                        ACVI9234
      lml = iqnums( 1 )                                                 ACVI9235
      mul = iqnums( 2 )                                                 ACVI9236
      lmr = iqnums( nmax+3 )                                            ACVI9237
      mur = iqnums( nmax+4 )                                            ACVI9238
C                                                                       ACVI9239
C     Initialize operator quantum numbers.                              ACVI9240
C                                                                       ACVI9241
C     <|...(2Jop)|> = dxwu3(2Jop,kro) * <||...||>kro                    ACVI9242
C                                                                       ACVI9243
C     Calculate SU(3) Wigner matrix coefficients dxwu3.                 ACVI9244
      if( tfile ) write(lfile,'(/a)')' SU(3) Wigner matrix coefficients'ACVI9245
      ieop   = (2*lml + mul) - (2*lmr+mur)                              ACVI9246
      mlmop  = mul - mur                                                ACVI9247
      jtmin  = iabs( mlmop )                                            ACVI9248
      jtmax  = mul + mur                                                ACVI9249
      ix     = 0                                                        ACVI9250
C     iy     = 0                                                        ACVI9251
      i3max  = 0                                                        ACVI9252
      ixromn = 0                                                        ACVI9253
      do i = 1, ilm                                                     ACVI9254
        index = iqops( 2, i )                                           ACVI9255
C       write(6,'(/2z10)') iqops(1,i), iqops(2,i)                       ACVI9256
        kro12 = iupac( index, -24, nbit8 )                              ACVI9257
        if( tconj ) then                                                ACVI9258
          lm  = iupac( index,  -8, nbit8 )                              ACVI9259
          mu  = iupac( index, -16, nbit8 )                              ACVI9260
        else                                                            ACVI9261
          lm  = iupac( index, -16, nbit8 )                              ACVI9262
          mu  = iupac( index,  -8, nbit8 )                              ACVI9263
        end if                                                          ACVI9264
        call u3mult( lmr, mur, lm, mu, lml, mul, nro, *100 )            ACVI9265
C         iy = iy + 1                                                   ACVI9266
CW----------------------------------------------------------------------ACVI9267
          if( tfile ) write(lfile,*)                                    ACVI9268
          if( tsmfi ) write(lfile,'(a,2i3,a)') ' (',lm,mu,')'           ACVI9269
CW----------------------------------------------------------------------ACVI9270
C         Calculate the extremal Wigner coefficients.                   ACVI9271
C                   <   ket  ;   T  |   bra  > v for nz,nx,ny order     ACVI9272
          call xewu3( lmr,mur, lm,mu, lml,mul, 0, NEC,                  ACVI9273
     1      kromax, indmax, dewu3, j1ta,j2ta,iea, NCE1,NCE2,KIMAX1 )    ACVI9274
          ixmin = ix                                                    ACVI9275
          do ind = 1, indmax                                            ACVI9276
            tj2ta( ind ) = iea( ind ).eq.ieop .and. j1ta( ind ).eq.mur  ACVI9277
     1        .and. jtmin.le.j2ta( ind ) .and. j2ta( ind ).le.jtmax     ACVI9278
            if( tj2ta(ind) ) ix = ix + 1                                ACVI9279
          end do                                                        ACVI9280
          ixmax = ix - ixmin                                            ACVI9281
          if( ixmax .lt. kromax )                                       ACVI9282
     1      call attn(' WIGMAT: System equation under-determined',*999) ACVI9283
          indro = 0                                                     ACVI9284
          ix = ixmin                                                    ACVI9285
          do ind = 1, indmax                                            ACVI9286
            if( tj2ta(ind) ) then                                       ACVI9287
              j2t  = j2ta( ind )                                        ACVI9288
              dwu2 = dwr3( mur,j2t,mul, mur,mlmop,mul )                 ACVI9289
              ix   = ix + 1                                             ACVI9290
              if( ix .gt. MXJT ) then                                   ACVI9291
                i3max = 0                                               ACVI9292
                call attn(' WIGMAT: Increase MXJT!',*999)               ACVI9293
              end if                                                    ACVI9294
              ixromn     = ixromn + 1                                   ACVI9295
              lmjt( ix ) = ior( index, j2t )                            ACVI9296
              ixro = ixromn                                             ACVI9297
              do kro = 1, kromax                                        ACVI9298
                indro = indro + 1                                       ACVI9299
                dxwu3( ixro ) = dewu3( indro ) * dwu2                   ACVI9300
                ixro  = ixro + ixmax                                    ACVI9301
              end do                                                    ACVI9302
CW----------------------------------------------------------------------ACVI9303
CW    Write Wigner coefficients                                         ACVI9304
      if( tfile ) write(lfile,1000)(dxwu3(j), j=ixromn,ixro-ixmax,ixmax)ACVI9305
 1000 format(10f8.4)                                                    ACVI9306
CW----------------------------------------------------------------------ACVI9307
            end if                                                      ACVI9308
            indro = ind * kromax                                        ACVI9309
          end do ! ind                                                  ACVI9310
          ixromn = ixro - ixmax                                         ACVI9311
          iqops( 2, i ) = ior( index, nro )                             ACVI9312
C         i3ptr( iy )   = ixmax                                         ACVI9313
          i3ptr( i )    = ixmax                                         ACVI9314
          i3max  = i3max + ixmax                                        ACVI9315
100     continue                                                        ACVI9316
      end do                                                            ACVI9317
C     if( tsmfi ) then                                                  ACVI9318
C       write(lfile,*) 'Pointers of SU(3) operators'                    ACVI9319
C       write(lfile,1010) (i3ptr(i),iqops(1,i),iqops(2,i), i=1,ilm )    ACVI9320
C1010   format(4(i5,2z10))                                              ACVI9321
C     end if                                                            ACVI9322
C                                                                       ACVI9323
C     LU decomposition of Wigner Matrix Coefficients                    ACVI9324
      if( i3max .gt. 0 ) then                                           ACVI9325
        ix = 0                                                          ACVI9326
        ir = 1                                                          ACVI9327
        i3first = 1                                                     ACVI9328
        do i = 1, ilm                                                   ACVI9329
          kromax = iupac( iqops(2,i), 0, nbit8 )                        ACVI9330
          ixmax = i3ptr( i )                                            ACVI9331
          if( ixmax .gt. 0 ) then                                       ACVI9332
            do j = 1, ixmax                                             ACVI9333
              ix = ix + 1                                               ACVI9334
              ia( ix ) = j                                              ACVI9335
            end do                                                      ACVI9336
            call dlut( ixmax, kromax, ia(i3first), dxwu3(ir), ixmax )   ACVI9337
            i3first = i3first + ixmax                                   ACVI9338
            ir = ir + ixmax*kromax                                      ACVI9339
          end if                                                        ACVI9340
        end do                                                          ACVI9341
      end if                                                            ACVI9342
      return                                                            ACVI9343
999   return 1                                                          ACVI9344
C ---*end of wigmat*----------------------------------------------------ACVI9345
      end                                                               ACVI9346
CB----------------------------------------------------------------------ACVI9347
C                                                                       ACVI9348
C                          **************                               ACVI9349
C                          *** ME1BAR ***                               ACVI9350
C                          **************                               ACVI9351
C                                                                       ACVI9352
C ----------------------------------------------------------------------ACVI9353
C Author:  Chairul Bahri                                                ACVI9354
C ----------------------------------------------------------------------ACVI9355
C Updates: 11/90 ==> Original ... (in SETRME)                           ACVI9356
C          02/25/91: multi irreps.                                      ACVI9357
C          06/01/91: separated from SETRME.                             ACVI9358
C          06/18/91: a+ and a modified.                                 ACVI9359
C          08/17/91: WST and OPBASIS implemented                        ACVI9360
C ----------------------------------------------------------------------ACVI9361
C                                                                       ACVI9362
C    The subroutine ME1BAR calculates matrix elements for a given set ofACVI9363
C    quantum numbers for BRA(L), KET(R), and OPERATOR(T) using the      ACVI9364
C    SU(3) > SU(2) x U(1) scheme.                                       ACVI9365
C                                                                       ACVI9366
C    The single barred matrix elements (overlaps) between L,R,T are     ACVI9367
C    evaluated by using bit manipulations (for determining the phase.)  ACVI9368
C                                                                       ACVI9369
C    The present calculation is for the triple barred                   ACVI9370
C       SU(3) x SU(2) only.                                             ACVI9371
C --------------------------------------------------------------------- ACVI9372
      subroutine me1bar( iqops, ispta, i3max, nbras, nkets, nal, nar,   ACVI9373
     ]                   kvecl, kvecr, dsolnl, dsolnr, mxjtck, nmax,    ACVI9374
     ]                   ibtree,ictree, lwigsu3,dwigsu3, mxwcof )       ACVI9375
      implicit real*8(d), logical(t)                                    ACVI9376
      parameter ( DZERO = 0.d0 )                                        ACVI9377
      character*4 cphase                                                ACVI9378
      dimension iqops( 2, * ),     ! packed qn's for operators          ACVI9379
     1          ispta( * ),        ! SU(2)spin q.n.                     ACVI9380
     2          kvecl( nmax, * ),  !                                    ACVI9381
     3          kvecr( nmax, * ),  ! bit states                         ACVI9382
     4          dsolnl( nbras, * ),!                                    ACVI9383
     5          dsolnr( nkets, * ),! coefficients of bit states         ACVI9384
     6          ibtree(-10:* ),    ! SU(3)-coupled tensors              ACVI9385
     7          ictree(-10:* ),    ! SU(3)-uncoupled tensors            ACVI9386
     8          lwigsu3(-2:* ),    ! linked-lists for SU(3) tensor coef.ACVI9387
     9          dwigsu3( * )       ! SU(3) tensor coefficients          ACVI9388
C               nbras, nkets       ! # bit states in hws                ACVI9389
C               nal, nar           ! U(N)>SU(3) alpha multiplicity      ACVI9390
C               i3max              ! # matrix elements                  ACVI9391
C               mxjtck             ! check MXJT                         ACVI9392
C               nmax               ! multiplicity                       ACVI9393
      parameter( MXJT = 1000 )     ! max # SU(3)>SU(2)xU(1) operators   ACVI9394
      common / SU3U21 /                                                 ACVI9395
     1         dxwu3( 9*MXJT ),    ! Wigner coefficient matrices        ACVI9396
     2         ia( MXJT ),         ! indices for LU decomposition       ACVI9397
     3         lmjt( MXJT ),       ! packed ro(lm mu)jt labels for ops. ACVI9398
     4         i3ptr( MXJT ),      ! pointers to lmjt(*)                ACVI9399
     5         ieop, mlmop         ! SU(3)>SU(2)xU(1) operators (I/O)   ACVI9400
      common / OVRLAP /                                                 ACVI9401
     1         dovrlp( 100*MXJT )  ! overlap = < | | >                  ACVI9402
      common / LOGCON / tfile,     ! write to logfile (full)            ACVI9403
     1                  tsmfi,     ! write to logfile (small)           ACVI9404
     2                  tconj,     ! conjugate SU(3) label              ACVI9405
     3                  topbase    ! calculate SU(3) coupled opb        ACVI9406
      common / RMECON / lfile,     ! log file for intermediate results  ACVI9407
     1                  ieta       ! oscillator shell number            ACVI9408
      common / RMEBIT / nsize,     ! # of levels for each word          ACVI9409
     1                  nwords,    ! # of words in each bit state       ACVI9410
     2                  kact,      ! # active words for bit manipulat.  ACVI9411
     3                  nbits,     ! bit length for # levels in shell   ACVI9412
     4                  nbitz,     ! significant bits                   ACVI9413
     5                  nzxy( 66 ),! (nz,nx,ny) packed labels for 1-prt.ACVI9414
     6                  nmlvls( 3 )! pointers for next words            ACVI9415
      common / RMEIND / ilm,       ! # operators                        ACVI9416
     1                  i3ix,      ! last index for overlaps            ACVI9417
     2                  npmax,     ! # a+                               ACVI9418
     3                  nhmax,     ! # a                                ACVI9419
     4                  npop1      ! (npmax,nhmax)                      ACVI9420
      common / OPBSU2 / dwigsu2(9),! SU(2) coupled operator bases       ACVI9421
     1                  ismin,     ! 2 Smin                             ACVI9422
     2                  ismax      ! 2 Smax                             ACVI9423
C                                                                       ACVI9424
C     Bit manipulation arrays                                           ACVI9425
      dimension label( 4 ), mst( 2 )                                    ACVI9426
      dimension kvexl( 30 ), kvexr( 30 ), nstay( 30 ), lvls( 66 )       ACVI9427
      dimension idist( MXJT ), mdist( MXJT ), ntemp( MXJT ),            ACVI9428
     1          iapls( MXJT ), iamin( MXJT )                            ACVI9429
C                                                                       ACVI9430
C     SU(3) Wigner arrays                                               ACVI9431
      dimension dwsu3( 9 ), idat( 2 )                                   ACVI9432
C                                                                       ACVI9433
C     Unpacking function.                                               ACVI9434
      data nbit8 / 255 /           ! zff                                ACVI9435
      iupac( index, iover, ibits ) = iand( ishft( index,iover ), ibits )ACVI9436
C                                                                       ACVI9437
      if( mxjtck .ne. MXJT ) call error(' ME1BAR: Check MXJT!')         ACVI9438
      ibits  = 1                                                        ACVI9439
      maskl5 = nsize - 1                                                ACVI9440
      nlevel = (ieta + 1)*(ieta + 2)/2                                  ACVI9441
      nz     = mod( nlevel, nsize )                                     ACVI9442
      if( max0(npmax,nhmax) .eq. 4 ) ibits = 2                          ACVI9443
      kactm  = kact - nmax                                              ACVI9444
C                                                                       ACVI9445
C     Reset dovrlp( * )                                                 ACVI9446
      do i = 1, i3ix                                                    ACVI9447
        dovrlp( i ) = 0.d0                                              ACVI9448
      end do                                                            ACVI9449
C                                                                       ACVI9450
C     Reset phase                                                       ACVI9451
      nphas1 = npmax * ieta                                             ACVI9452
      nphas2 = nhmax * ieta                                             ACVI9453
      npop   = npmax - nhmax                                            ACVI9454
      npmin  = max0( 0, npop )                                          ACVI9455
      nhmin  = max0( 0,-npop )                                          ACVI9456
C                                                                       ACVI9457
      if( tfile ) then                                                  ACVI9458
        write(lfile,1000)                                               ACVI9459
        write(lfile,1010) i3max                                         ACVI9460
C       write(lfile,'(8z10)') ( lmjt(i), i=1,i3max )                    ACVI9461
 1000   format(//' Overlaps'/)                                          ACVI9462
 1010   format(  ' The size of SU(3) operators',i5)                     ACVI9463
      end if                                                            ACVI9464
C                                                                       ACVI9465
C     Estimate the size of dovrlp                                       ACVI9466
      i3max = 0                                                         ACVI9467
      do i = 1, ilm                                                     ACVI9468
        index   = iqops( 2, i )                                         ACVI9469
        kro12mx = iupac( index, -24, nbit8 )                            ACVI9470
        index   = ispta( i )                                            ACVI9471
        is1t    = iupac( index, -16, nbit8 )                            ACVI9472
        is2t    = iupac( index,   0, nbit8 )                            ACVI9473
        istmin  = max0( ismin, iabs(is1t-is2t) )                        ACVI9474
        istmax  = min0( ismax, is1t+is2t )                              ACVI9475
        nst     = max0( 0, (istmax-istmin)/2 + 1 )                      ACVI9476
        i3max   = i3max + nst*kro12mx*i3ptr( i )                        ACVI9477
      end do                                                            ACVI9478
      if( tfile ) write(lfile,1020) i3max                               ACVI9479
 1020 format(/' The size of matrix elements',i6)                        ACVI9480
C                                                                       ACVI9481
C     Starting ...                                                      ACVI9482
      il = 0                                                            ACVI9483
      do ibra = 1, nbras                                                ACVI9484
        k = 0                                                           ACVI9485
        do j = 1, nwords                                                ACVI9486
          il = il + 1                                                   ACVI9487
          do iun = 1, nmax                                              ACVI9488
            k = k + 1                                                   ACVI9489
            kvexl( k ) = kvecl( iun, il )                               ACVI9490
          end do                                                        ACVI9491
        end do                                                          ACVI9492
        ir = 0                                                          ACVI9493
        do iket = 1, nkets                                              ACVI9494
          k = 0                                                         ACVI9495
          do j = 1, nwords                                              ACVI9496
            ir = ir + 1                                                 ACVI9497
            do iun = 1, nmax                                            ACVI9498
              k = k + 1                                                 ACVI9499
              kvexr( k ) = kvecr( iun, ir )                             ACVI9500
            end do                                                      ACVI9501
          end do                                                        ACVI9502
C         Generate all possible actions of creation and annihilation    ACVI9503
C         operators.                                                    ACVI9504
CW        if( tfile ) write(lfile,'(/a/)') ' Second Quantized operators'ACVI9505
          do k = 1, kact                                                ACVI9506
            nmove      = ieor( kvexl( k ), kvexr( k ) )                 ACVI9507
            iamin( k ) = iand( nmove     , kvexr( k ) )                 ACVI9508
            iapls( k ) = iand( nmove     , kvexl( k ) )                 ACVI9509
            nstay( k ) = iand( kvexl( k ), kvexr( k ) )                 ACVI9510
          end do                                                        ACVI9511
C         Generate pointers for the actions of s.q.op.                  ACVI9512
          nh = 0                                                        ACVI9513
          np = 0                                                        ACVI9514
          ns = 0                                                        ACVI9515
          do iun = 1, nmax                                              ACVI9516
            k = iun                                                     ACVI9517
            l = 0                                                       ACVI9518
            do j = 1, nwords                                            ACVI9519
              iam = iamin( k )                                          ACVI9520
              iap = iapls( k )                                          ACVI9521
              nst = nstay( k )                                          ACVI9522
              do i = 0, nmlvls( j )                                     ACVI9523
                l = l + 1                                               ACVI9524
                if( btest( iam, i ) ) nh = nh + 1                       ACVI9525
                if( btest( iap, i ) ) np = np + 1                       ACVI9526
                if( btest( nst, i ) ) then                              ACVI9527
                  ns = ns + 1                                           ACVI9528
                  lvls( ns ) = ior( ishft(k,8), l )                     ACVI9529
                end if                                                  ACVI9530
              end do                                                    ACVI9531
              k = k + nmax                                              ACVI9532
            end do                                                      ACVI9533
          end do ! iun                                                  ACVI9534
          tpart = npmin.le.np .and. np.le.npmax .and.                   ACVI9535
     +            nhmin.le.nh .and. nh.le.nhmax                         ACVI9536
C                                                                       ACVI9537
C         Take the actions if within the possible range.                ACVI9538
          if( tpart ) then                                              ACVI9539
CW----------------------------------------------------------------------ACVI9540
CW    Write the overlap                                                 ACVI9541
      if( tfile ) write(lfile,1030) ibra, npmax, nhmax, iket            ACVI9542
 1030 format(/' < #',i4,'  | a+_',i1,' a_',i1,' | #',i4,'  >')          ACVI9543
CW----------------------------------------------------------------------ACVI9544
            nst  = npmax - np                                           ACVI9545
            irun = 0                                                    ACVI9546
C           Distribute particles among all possibilities.               ACVI9547
            do j = 1, nst                                               ACVI9548
              mdist( j ) = j                                            ACVI9549
            end do                                                      ACVI9550
C           Generate the first operator.                                ACVI9551
C           Label of level is                                           ACVI9552
C           xxxx xxxx xxxx xxxx xxxx xxxx xxxx xxxx                     ACVI9553
C                                    \__/    \____/                     ACVI9554
C                                      k        l                       ACVI9555
            do k = 1, kact                                              ACVI9556
              ntemp( k ) = 0                                            ACVI9557
            end do                                                      ACVI9558
            do i = 1, nst                                               ACVI9559
              level = lvls( i )                                         ACVI9560
              k     = ishft( level, -8 )                                ACVI9561
              l     = iand ( level, maskl5 ) - 1                        ACVI9562
              ntemp( k ) = ibset( ntemp( k ), l )                       ACVI9563
            end do                                                      ACVI9564
            do k = 1, kact                                              ACVI9565
              idist( k ) = ntemp( k )                                   ACVI9566
            end do                                                      ACVI9567
            irun = kact                                                 ACVI9568
100         do j = nst, 1, -1                                           ACVI9569
              nstmj = nst - j                                           ACVI9570
              if( mdist(j) .lt. (ns-nstmj) ) then                       ACVI9571
                mdist( j ) = mdist( j ) + 1                             ACVI9572
                do jj = 1, nstmj                                        ACVI9573
                  mdist( j + jj ) = mdist( j ) + jj                     ACVI9574
                end do                                                  ACVI9575
C               Generate all possible operators.                        ACVI9576
                do k = 1, kact                                          ACVI9577
                  ntemp( k ) = 0                                        ACVI9578
                end do                                                  ACVI9579
                do i = 1, nst                                           ACVI9580
                  level = lvls( mdist(i) )                              ACVI9581
                  k     = ishft( level, -8 )                            ACVI9582
                  l     = iand ( level, maskl5 ) - 1                    ACVI9583
                  ntemp( k ) = ibset( ntemp( k ), l )                   ACVI9584
                end do                                                  ACVI9585
                do k = 1, kact                                          ACVI9586
                  irun = irun + 1                                       ACVI9587
                  if( irun .gt. MXJT )                                  ACVI9588
     1              call error(' ME1BAR: # distributions > MXJT.')      ACVI9589
                  idist( irun ) = ntemp( k )                            ACVI9590
                end do                                                  ACVI9591
                go to 100                                               ACVI9592
              end if                                                    ACVI9593
            end do                                                      ACVI9594
C                                                                       ACVI9595
C           Generate all possible creation and annihilation operators.  ACVI9596
            do i = irun, 1, -1                                          ACVI9597
              k = mod( i, kact )                                        ACVI9598
              if( k .eq. 0 ) k = kact                                   ACVI9599
              iapls( i ) = ior( iapls( k ), idist( i ) )                ACVI9600
              iamin( i ) = ior( iamin( k ), idist( i ) )                ACVI9601
            end do                                                      ACVI9602
            irun = irun / kact                                          ACVI9603
            if( tfile ) write(lfile,1040) irun                          ACVI9604
 1040       format(' Number of possible operators',i10/)                ACVI9605
C           Generate T from the couplings of a+s and a's. Keep track    ACVI9606
C           of generations.                                             ACVI9607
            irn = 0                                                     ACVI9608
            do n = 1, irun                                              ACVI9609
              iphase = 0          ! Phase in normal ordering (Wick's)   ACVI9610
              lvlh   = 0                                                ACVI9611
              lvlp   = 0                                                ACVI9612
              nbtwnh = 0                                                ACVI9613
              nbtwnp = npmax / 2                                        ACVI9614
              nh   = 0                                                  ACVI9615
              np   = 0                                                  ACVI9616
              msh  = 0                                                  ACVI9617
              msp  = 0                                                  ACVI9618
              isun = 1                                                  ACVI9619
              do iun = 1, nmax                                          ACVI9620
                k  = iun                                                ACVI9621
                kr = irn + iun                                          ACVI9622
                l  = 0                                                  ACVI9623
                do j = 1, nwords                                        ACVI9624
                  ket    = kvexr( k )                                   ACVI9625
                  iam    = iamin( kr )                                  ACVI9626
                  iap    = iapls( kr )                                  ACVI9627
                  newket = ieor( ket, iam )                             ACVI9628
                  do i = 0, nmlvls(j)                                   ACVI9629
                    l = l + 1                                           ACVI9630
                    if( btest( iam, i ) ) then    ! a(i)                ACVI9631
                      nh     = nh + 1                                   ACVI9632
                      iphase = iphase + nbtwnh                          ACVI9633
                      lvlh   = ior( ishft(lvlh, nbits), l )             ACVI9634
                      msh    = msh + isun                               ACVI9635
                      nbtwnh = nbtwnh - 1                               ACVI9636
                    end if                                              ACVI9637
                    if( btest( iap, i ) ) then    ! a+(i)               ACVI9638
                      np     = np + 1                                   ACVI9639
                      iphase = iphase + nbtwnp                          ACVI9640
                      lvlp   = ior( ishft(lvlp, nbits), l )             ACVI9641
                      msp    = msp + isun                               ACVI9642
                    end if                                              ACVI9643
                    if( np.ge.npmax .and. nh.ge.nhmax ) go to 110       ACVI9644
                    if( btest(ket,i) )    nbtwnh = nbtwnh + 1           ACVI9645
                    if( btest(newket,i) ) nbtwnp = nbtwnp + 1           ACVI9646
                  end do                                                ACVI9647
                  k  = k + nmax                                         ACVI9648
                  kr = kr + nmax                                        ACVI9649
                end do ! j                                              ACVI9650
                isun = ishft( isun, nbits )                             ACVI9651
              end do ! iun                                              ACVI9652
110           continue                                                  ACVI9653
C             Space order preferred to save the storage for op.basis.   ACVI9654
              mphase = iphase                                           ACVI9655
              tphas1 = .false.                                          ACVI9656
              tphas2 = .false.                                          ACVI9657
              if( nhmax .gt. 1 )                                        ACVI9658
     1          call piksrt( nbits, nbitz, nhmax, lvlh, tphas2 )        ACVI9659
              if( npmax .gt. 1 )                                        ACVI9660
     1          call piksrt( nbits, nbitz, npmax, lvlp, tphas1 )        ACVI9661
C             Spin degrees of freedom.                                  ACVI9662
              index = msh                                               ACVI9663
              do k = 1, 2                                               ACVI9664
                mst( k ) = iand( index, nbitz )                         ACVI9665
                index    = ishft( index, -nbits )                       ACVI9666
              end do                                                    ACVI9667
              mshu = mst( 1 )                                           ACVI9668
              msh  = mst( 2 ) - mshu                                    ACVI9669
              index = msp                                               ACVI9670
              do k = 1, 2                                               ACVI9671
                mst( k ) = iand( index, nbitz )                         ACVI9672
                index    = ishft( index, -nbits )                       ACVI9673
              end do                                                    ACVI9674
              mspu = mst( 1 )                                           ACVI9675
              msp  = mspu - mst( 2 )                                    ACVI9676
              mso  = msp + msh                                          ACVI9677
C                                                                       ACVI9678
CW----------------------------------------------------------------------ACVI9679
CW    Write the overlaps                                                ACVI9680
      if( tfile ) then                                                  ACVI9681
      if( nlevel .gt. nsize ) then                                      ACVI9682
        write(lfile,'(39x,1h*)')                                        ACVI9683
        do i = 1, kactm                                                 ACVI9684
          call dbovrl ( kvexl(i),iapls(irn+i), iamin(irn+i), kvexr(i) ) ACVI9685
        end do                                                          ACVI9686
        do i = kactm + 1, kact                                          ACVI9687
          call dbovrp( nz,kvexl(i),iapls(irn+i),iamin(irn+i),kvexr(i) ) ACVI9688
        end do                                                          ACVI9689
      else                                                              ACVI9690
        do i = 1, kact                                                  ACVI9691
          call dbovrc( nz,kvexl(i),iapls(irn+i),iamin(irn+i),kvexr(i) ) ACVI9692
        end do                                                          ACVI9693
      end if                                                            ACVI9694
      write(lfile,1050)                                                 ACVI9695
      end if                                                            ACVI9696
 1050 format(' LxM 2J', 5x,                                             ACVI9697
     1 '< al  # | ( lm mu)2s( lm mu)2s ro( lm mu) ep 2J 2M;2S | al # >')ACVI9698
CW----------------------------------------------------------------------ACVI9699
C                                                                       ACVI9700
              if( tconj ) then                                          ACVI9701
                label( 3 ) = lvlh                                       ACVI9702
                label( 4 ) = lvlp                                       ACVI9703
              else                                                      ACVI9704
                label( 3 ) = lvlp                                       ACVI9705
                label( 4 ) = lvlh                                       ACVI9706
              end if                                                    ACVI9707
C             Run over SU(3) qns.                                       ACVI9708
              i3first = 1                                               ACVI9709
              i3last  = i3ptr( 1 )                                      ACVI9710
              i3x     = 0                                               ACVI9711
C             i3x     = 1                                               ACVI9712
              do i = 1, ilm                                             ACVI9713
                index      = iqops( 1, i )                              ACVI9714
                label( 1 ) = index                                      ACVI9715
                lm1        = iupac( index, -24, nbit8 )                 ACVI9716
                mu1        = iupac( index, -16, nbit8 )                 ACVI9717
                lm2        = iupac( index,  -8, nbit8 )                 ACVI9718
                mu2        = iupac( index,   0, nbit8 )                 ACVI9719
                index      = iqops( 2, i )                              ACVI9720
                label( 2 ) = index                                      ACVI9721
                kro12mx    = iupac( index, -24, nbit8 )                 ACVI9722
                lm         = iupac( index, -16, nbit8 )                 ACVI9723
                mu         = iupac( index,  -8, nbit8 )                 ACVI9724
                kromax     = iupac( index,   0, nbit8 )                 ACVI9725
                index      = ispta( i )                                 ACVI9726
                is1t       = iupac( index, -16, nbit8 )                 ACVI9727
                is2t       = iupac( index,   0, nbit8 )                 ACVI9728
C               write(lfile,'(2h #,2z10)') iqops(1,i), iqops(2,i)       ACVI9729
C               write(6,'(2h #,3z10)') iqops(1,i), iqops(2,i),ispta(i)  ACVI9730
                if( tconj ) then                                        ACVI9731
                  index = lm                                            ACVI9732
                  lm    = mu                                            ACVI9733
                  mu    = index                                         ACVI9734
                end if                                                  ACVI9735
C               SU(3) multi-particle phase                              ACVI9736
                iphase   = mphase                                       ACVI9737
                if( tphas1 ) iphase = iphase + nphas1 + lm1 + mu1       ACVI9738
                if( tphas2 ) iphase = iphase + nphas2 + lm2 + mu2       ACVI9739
                if( nhmax .ge. 1 ) iphase = iphase + (is2t+msh)/2       ACVI9740
                tphase = btest( iphase, 0 )                             ACVI9741
C               SU(2) spin coupling                                     ACVI9742
                istmin = max0( ismin, iabs(is1t-is2t) )                 ACVI9743
                istmax = min0( ismax, is1t+is2t )                       ACVI9744
                nst    = max0( 0, (istmax-istmin)/2 + 1 )               ACVI9745
                if( nst .le. 0 ) then                                   ACVI9746
                  nst     = 1                                           ACVI9747
                  kro12mx = 1                                           ACVI9748
                  go to 140                                             ACVI9749
                end if                                                  ACVI9750
C                                                                       ACVI9751
C               Retrieving Wigner coefficients                          ACVI9752
                dwu2 = 1.d0                                             ACVI9753
                if( npmax .gt. 1 ) then                                 ACVI9754
                  k = ior( ishft(mspu,ibits), is1t/2 )                  ACVI9755
                  dwu2 = dwu2 * dwigsu2( k )                            ACVI9756
                end if                                                  ACVI9757
                if( nhmax .gt. 1 ) then                                 ACVI9758
                  k = ior( ishft(mshu,ibits),is2t/2)                    ACVI9759
                  dwu2 = dwu2 * dwigsu2( k )                            ACVI9760
                end if                                                  ACVI9761
                i3size = i3ptr( i )                                     ACVI9762
                do i3 = i3first, i3last                                 ACVI9763
                  i3x = i3x + 1                                         ACVI9764
                  dcwu3 = 0.d0                                          ACVI9765
                  if( npop1 .eq. 1 ) then                               ACVI9766
                    dcwu3 = 1.d0                                        ACVI9767
                  else                                                  ACVI9768
C                   Retrieving SU(3) Wigner from tree.                  ACVI9769
                    index      = lmjt( i3 )                             ACVI9770
                    label( 2 ) = index                                  ACVI9771
                    jt         = iupac( index, 0, nbit8 )               ACVI9772
                    call tchk( label, ibtree, *120 )                    ACVI9773
                      if( topbase ) then                                ACVI9774
                        call optlm2( label, ictree, dwigsu3, dwsu3,     ACVI9775
     1                               mxwcof, npmax, -nhmax, *130 )      ACVI9776
                        call tadd(   label,idat, dwsu3, kro12mx,kro12mx,ACVI9777
     1                               ibtree, dwigsu3, lwigsu3 )         ACVI9778
                        dcwu3 = 1.d0                                    ACVI9779
                          go to 130                                     ACVI9780
                      else                                              ACVI9781
                          go to 130                                     ACVI9782
                      end if                                            ACVI9783
C                                                                       ACVI9784
120                   iloc  = ibtree(-5 ) + ibtree(-4 ) + 1             ACVI9785
                      ix    = ibtree( iloc )                            ACVI9786
                      ixmax = ibtree( iloc+1 )                          ACVI9787
                      if( ixmax .ne. kro12mx ) call error(' ME1BAR:?$#')ACVI9788
                      do kro12 = 1, kro12mx                             ACVI9789
                        dwsu3( kro12 ) = dwigsu3( ix )                  ACVI9790
                        ix = lwigsu3( ix )                              ACVI9791
                      end do                                            ACVI9792
                      dcwu3 = 1.d0                                      ACVI9793
130                 continue                                            ACVI9794
                  end if                                                ACVI9795
C                                                                       ACVI9796
                  i3alr = i3x - i3max                                   ACVI9797
                  dwig  = dwu2 * dcwu3                                  ACVI9798
C                                                                       ACVI9799
                  do ial = 1, nal                                       ACVI9800
                    dcl = dwig * dsolnl( ibra, ial )                    ACVI9801
                    do iar = 1, nar                                     ACVI9802
                      i3alr = i3alr + i3max                             ACVI9803
                      dclr  = dcl * dsolnr( iket, iar )                 ACVI9804
                      i3ix  = i3alr - i3size                            ACVI9805
C                                                                       ACVI9806
C                     Run over spin q.n.                                ACVI9807
                      do ist = istmin, istmax, 2                        ACVI9808
                        dcwu2 = dclr                                    ACVI9809
C                       Retrieve SU(2)S Wigner ...                      ACVI9810
                        if( npop1 .ge. 256 ) dcwu2 = dcwu2 *            ACVI9811
     1                      dwr3( is1t,is2t,ist, msp,msh,mso )          ACVI9812
                        do kro12 = 1, kro12mx                           ACVI9813
                          i3ix = i3ix + i3size                          ACVI9814
                          if( i3ix.gt.100*MXJT )                        ACVI9815
     1                      call error(' ME1BAR: dim of DOVRLP!')       ACVI9816
                          dweight = dcwu2                               ACVI9817
                          if( npop1.gt.1 ) dweight = dcwu2*dwsu3(kro12) ACVI9818
                          if( tphase ) then                             ACVI9819
                            dovrlp( i3ix ) = dovrlp(i3ix) - dweight     ACVI9820
                          else                                          ACVI9821
                            dovrlp( i3ix ) = dovrlp(i3ix) + dweight     ACVI9822
                          end if                                        ACVI9823
CW----------------------------------------------------------------------ACVI9824
CW    Write overlap                                                     ACVI9825
      if( tfile .and. dabs(dweight).gt.DZERO ) then                     ACVI9826
        if( tphase ) then                                               ACVI9827
          cphase = ' = -'                                               ACVI9828
        else                                                            ACVI9829
          cphase = ' = +'                                               ACVI9830
        end if                                                          ACVI9831
        if( kro12 .eq. 1 ) then                                         ACVI9832
          write(lfile,1060) i,i3,i3ix, ial,ibra,lm1,mu1,is1t,           ACVI9833
     1      lm2,mu2,is2t, kro12,lm,mu,ieop,jt,mlmop,ist,iar,iket,       ACVI9834
     2      cphase,dweight,ibtree(iloc)                                 ACVI9835
        else                                                            ACVI9836
          write(lfile,1070) i3ix, kro12, cphase,dweight                 ACVI9837
        end if                                                          ACVI9838
      end if                                                            ACVI9839
 1060 format(i3,i4,i5,'<',2i3,' | ', 2('(',2i3,')',i2), i3,'(',2i3,')', ACVI9840
     1       3i3,';', i2,' |',2i3,' >',a,f9.6, ' (#',i6,')')            ACVI9841
 1070 format(7x,i5,30x,i3,30x,a,f9.6)                                   ACVI9842
CW----------------------------------------------------------------------ACVI9843
                        end do ! kro12                                  ACVI9844
                      end do ! ist                                      ACVI9845
                    end do ! iar                                        ACVI9846
                  end do ! ial                                          ACVI9847
                end do ! i3: lm,mu,jt                                   ACVI9848
C                                                                       ACVI9849
140             continue                                                ACVI9850
                i3first = i3last + 1                                    ACVI9851
                i3last  = i3last + i3ptr( i+1 )                         ACVI9852
C               i3x = i3x + (nst-1) * (kro12mx-1) * i3size              ACVI9853
                i3x = i3x + (nst*kro12mx-1) * i3size                    ACVI9854
C               i3x = i3ix                                              ACVI9855
              end do ! i: lm,mu                                         ACVI9856
C                                                                       ACVI9857
              irn = irn + kact                                          ACVI9858
            end do ! n                                                  ACVI9859
C                                                                       ACVI9860
          end if                                                        ACVI9861
        end do                                                          ACVI9862
      end do                                                            ACVI9863
      if( tfile ) then                                                  ACVI9864
        write(lfile,'(/i5,a)') i3ix,' elements:'                        ACVI9865
        write(lfile,'(/10(10f8.3/))') ( dovrlp(i), i=1, i3ix )          ACVI9866
      end if                                                            ACVI9867
C                                                                       ACVI9868
      return                                                            ACVI9869
C ---*end of me1bar*----------------------------------------------------ACVI9870
      end                                                               ACVI9871
CB----------------------------------------------------------------------ACVI9872
C                                                                       ACVI9873
C                          **************                               ACVI9874
C                          *** RM3BAR ***                               ACVI9875
C                          **************                               ACVI9876
C                                                                       ACVI9877
C ----------------------------------------------------------------------ACVI9878
C Author:  Chairul Bahri                                                ACVI9879
C ----------------------------------------------------------------------ACVI9880
C Updates: 11/90 ==> Original                                           ACVI9881
C          12/13/90: test for a+ (fp)**8 & conjugate ... SUCCESSFUL !   ACVI9882
C          12/25/90: using uncoupled tensor T(LM LM)=a+ a               ACVI9883
C          01/13/91: test for a+a(ds)**2 ... SUCCESSFUL for RO<>1, too !ACVI9884
C                                (4 0) (0 2) (2 1)                      ACVI9885
C          02/25/91: multi irreps.                                      ACVI9886
C          01/20/92: Sun 4.                                             ACVI9887
C          01/07/93: RS/6000                                            ACVI9888
C ----------------------------------------------------------------------ACVI9889
C                                                                       ACVI9890
C    The subroutine RM3BAR calculates the SU(3) the Reduced Matrix      ACVI9891
C    Elements and/or the Coefficient of Fractional Parentages           ACVI9892
C    using LU decomposition as                                          ACVI9893
C                                                                       ACVI9894
C       <|...(J2T)|> = DXWU3(J2T,KRO) * <||...||>(KRO).                 ACVI9895
C                                                                       ACVI9896
C ----------------------------------------------------------------------ACVI9897
      subroutine rm3bar( iqops, ispta, iqwrit, nal, nar,                ACVI9898
     ]                   drme, irmtree, mxrme, irme, mxjtck, * )        ACVI9899
      implicit real*8(d), logical(t)                                    ACVI9900
      parameter ( DZERO = 1.d-12 )                                      ACVI9901
      dimension iqops( 2, * ),     ! packed qns for operators           ACVI9902
     1          ispta( * ),        !                                    ACVI9903
     2          iqwrit( * ),       ! quantum number array for write out ACVI9904
     3          drme( * ),         ! SU(3) rmes                         ACVI9905
     4          irmtree( -10:* )   ! binary tree for rmes               ACVI9906
C               nal, nar           ! U(N)>SU(3) alpha multiplicity      ACVI9907
C               mxrme              ! max dimension of drme              ACVI9908
C               irme               ! index for rme                      ACVI9909
C               mxjtck             ! check MXJT                         ACVI9910
      parameter( MXJT = 1000 )     ! max # SU(3)>SU(2)xU(1) operators   ACVI9911
      common / SU3U21 /                                                 ACVI9912
     1         dxwu3( 9*MXJT ),    ! Wigner coefficient matrices        ACVI9913
     2         ia( MXJT ),         ! indices for LU decomposition       ACVI9914
     3         lmjt( MXJT ),       ! packed ro(lm mu)jt labels for ops. ACVI9915
     4         i3ptr( MXJT ),      ! pointers to lmjt(*)                ACVI9916
     5         ieop, mlmop         ! SU(3)>SU(2)xU(1) operators         ACVI9917
      common / OVRLAP /                                                 ACVI9918
     1         dovrlp( 100*MXJT )  ! overlap = < | | >                  ACVI9919
      common / LOGCON / tfile,     ! write to logfile (full)            ACVI9920
     1                  tsmfi,     ! write to logfile (small)           ACVI9921
     2                  tconj,     ! conjugate SU(3) label              ACVI9922
     3                  topbase    ! calculate SU(3) coupled opb        ACVI9923
      common / RMECON / lfile,     ! = logfile                          ACVI9924
     1                  ieta       ! = neta                             ACVI9925
      common / RMEIND / ilm,       ! # operators                        ACVI9926
     1                  i3ix,      ! last index for overlaps            ACVI9927
     2                  npmax,     ! # a+                               ACVI9928
     3                  nhmax,     ! # a                                ACVI9929
     4                  npop1      ! (npmax,nhmax)                      ACVI9930
      common / OPBSU2 / dwigsu2(9),! SU(2) coupled operator bases       ACVI9931
     1                  ismin,     ! 2 Smin                             ACVI9932
     2                  ismax      ! 2 Smax                             ACVI9933
      common / WORKAR / dwork(MXJT)! working arrays                     ACVI9934
      dimension df( 0:7 )          ! antisymmetric factors              ACVI9935
      dimension label( 5 )         ! temp label for rme tree            ACVI9936
      data df / 1.d0, 1.d0, 2.d0, 6.d0, 24.d0, 120.d0, 720.d0, 5040.d0 /ACVI9937
C                                                                       ACVI9938
C     Packing and unpacking function.                                   ACVI9939
      data nbit8 / 255 /           ! zff                                ACVI9940
      ipack( i, j, k, l, m ) = ior(m, ishft( ior(l, ishft( ior(k,       ACVI9941
     1       ishft( ior(j, ishft( i,4 )),8 )),8 )),8 ))                 ACVI9942
      iupac( index, iover, ibits ) = iand( ishft( index,iover ), ibits )ACVI9943
C                                                                       ACVI9944
      if( mxjtck .ne. MXJT ) call error(' RM3BAR: Check MXJT!')         ACVI9945
      dmult = df( npmax ) * df( nhmax )                                 ACV19946
C                                                                       ACV19947
      jk = 18                                                           ACV19948
      jk2 = 19                                                          ACV19949
      japl = 6                                                          ACV19950
      jamn = 10                                                         ACV19951
      jtop = 14                                                         ACV19952
      label( 1 ) = iqwrit( 1 )                                          ACV19953
      label( 2 ) = iqwrit( jk )   ! jk = 18                             ACV19954
      islt       = iqwrit( 5 )                                          ACV19955
      isrt       = iqwrit( jk+4 )                                       ACV19956
      msopt      = islt - isrt                                          ACV19957
      iqwrit( japl ) = npmax      ! japl = 6                            ACV19958
      iqwrit( jamn ) = nhmax      ! jamn = 10                           ACV19959
C                                                                       ACV19960
      iq = 2                                                            ACV19961
      jq = jk2                                                          ACV19962
      do j = 1, 3                                                       ACV19963
        iq = iq + 1                                                     ACV19964
        jq = jq + 1                                                     ACV19965
        label( 1 ) = ior( ishft( label(1),8 ), iqwrit(iq) )             ACV19966
        label( 2 ) = ior( ishft( label(2),8 ), iqwrit(jq) )             ACV19967
      end do                                                            ACV19968
C                                                                       ACV19969
      ix = 0                                                            ACV19970
      do ial = 1, nal                                                   ACV19971
        iqwrit( 2 ) = ial                                               ACV19972
        do iar = 1, nar                                                 ACV19973
          iqwrit( jk2 ) = iar                                           ACV19974
C                                                                       ACV19975
CW----------------------------------------------------------------------ACV19976
CW    < bra ||| * ||| ket >                                             ACV19977
      write(6,1000)                                                     ACV19978
      write(6,1010) ( iqwrit(i), i=1,5 ), ( iqwrit(i), i=18,22 )        ACV19979
      if( tfile ) then                                                  ACV19980
        write(lfile,1000)                                               ACV19981
        write(lfile,1010) ( iqwrit(i), i=1,5 ), ( iqwrit(i), i=18,22 )  ACV19982
        write(lfile,1020)                                               ACV19983
      else                                                              ACV19984
        write(6,1020)                                                   ACV19985
      end if                                                            ACV19986
1000  format(/'  < NL : AL( LM MU) ST  |||  T  ||| NR : AL( LM MU) ST >'ACV19987
     2       )                                                          ACV19988
1010  format( 2(i6,i5,i4,i3,i4,10x) )                                   ACV19989
1020  format(/'    < bra ||| N1 (LM MU)ST(x)N2 (LM MU)ST;RO( LM MU)ST ' ACV19990
     2       ,'||| ket >RO')                                            ACV19991
CW----------------------------------------------------------------------ACV19992
          i3first = 1                                                   ACV19993
          ir = 1                                                        ACV19994
          inxold = -1                                                   ACV19995
C                                                                       ACV19996
          do i = 1, ilm                                                 ACV19997
C                                                                       ACV19998
            index  = iqops( 1, i )                                      ACV19999
            indu3c = index                                              ACV10000
C                                                                       ACV10001
            lm1    = iupac( index,-24, nbit8 )                          ACV10002
            mu1    = iupac( index,-16, nbit8 )                          ACV10003
            lm2    = iupac( index, -8, nbit8 )                          ACV10004
            mu2    = iupac( index,  0, nbit8 )                          ACV10005
C                                                                       ACV10006
            index  = ispta( i )                                         ACV10007
            is1t   = iupac( index, -16, nbit8 )                         ACV10008
            is2t   = iupac( index,   0, nbit8 )                         ACV10009
            istmin = max0( ismin, iabs(is1t-is2t) )                     ACV10010
            istmax = min0( ismax, is1t+is2t )                           ACV10011
                  iqwrit( japl+1 ) = lm1                                ACV10012
                  iqwrit( japl+2 ) = mu1                                ACV10013
                  iqwrit( japl+3 ) = is1t                               ACV10014
                  iqwrit( jamn+1 ) = lm2                                ACV10015
                  iqwrit( jamn+2 ) = mu2                                ACV10016
                  iqwrit( jamn+3 ) = is2t                               ACV10017
C                                                                       ACV10018
            index = iqops( 2, i )                                       ACV10019
            kro12mx = iupac( index,-24, nbit8 )                         ACV10020
            kromax  = iupac( index,  0, nbit8 )                         ACV10021
            if( tconj ) then                                            ACV10022
              lm = iupac( index, -8, nbit8 )                            ACV10023
              mu = iupac( index,-16, nbit8 )                            ACV10024
            else                                                        ACV10025
              lm = iupac( index,-16, nbit8 )                            ACV10026
              mu = iupac( index, -8, nbit8 )                            ACV10027
            end if                                                      ACV10028
                  iqwrit( jtop+1 ) = lm                                 ACV10029
                  iqwrit( jtop+2 ) = mu                                 ACV10030
            i3size  = i3ptr( i )                                        ACV10031
            if( i3size .gt. 0 ) then                                    ACV10032
C                                                                       ACV10033
              do ist = istmin, istmax, 2                                ACV10034
                dwu2 = dwr3( isrt,ist,islt, isrt,msopt,islt )           ACV10035
                  iqwrit( jtop + 3 ) = ist   ! jtop = 14                ACV10036
                do kro12 = 1, kro12mx                                   ACV10037
                  iqwrit( jtop )   = kro12                              ACV10038
                  call dbsr( kromax, dxwu3(ir), dovrlp(ix+1), dwork,    ACV10039
     1                     ia(i3first), i3size )                        ACV10040
CW----------------------------------------------------------------------ACV10041
CW    < *** ||| T ||| *** >                                             ACV10042
      iw = ix                                                           ACV10043
      do kro = 1, kromax                                                ACV10044
        iw = iw + 1                                                     ACV10045
        dme = dovrlp( iw )                                              ACV10046
        if( dabs( dme ) .gt. DZERO ) then                               ACV10047
          tflag = inxold .ne. indu3c                                    ACV10048
          dme = dme * dmult / dwu2                                      ACV10049
          if( tfile ) then                                              ACV10050
            if( tflag ) then                                            ACV10051
              write(lfile,1030)( iqwrit(j),j=japl,jtop+3 ), kro, dme    ACV10052
            else                                                        ACV10053
              write(lfile,1040)( iqwrit(j),j=jtop,jtop+3 ), kro, dme    ACV10054
            end if                                                      ACV10055
          else                                                          ACV10056
            if( tflag ) then                                            ACV10057
              write(6,1030)( iqwrit(j),j=japl,jtop+3 ), kro, dme        ACV10058
            else                                                        ACV10059
              write(6,1040)( iqwrit(j),j=jtop,jtop+3 ), kro, dme        ACV10060
            end if                                                      ACV10061
          end if                                                        ACV10062
C                                                                       ACV10063
          irme = irme + 1                                               ACV10064
          if(irme.gt.mxrme) call attn(' RM3BAR: Increase MXRME!',*999)  ACV10065
C                           (     4,    4,   8,   8,   8 )              ACV10066
          label( 3 ) = ipack( npmax, is1t, lm1, mu1, ial )              ACV10067
          label( 4 ) = ipack( nhmax, is2t, lm2, mu2, iar )              ACV10068
          label( 5 ) = ipack( ist,  kro12,  lm,  mu, kro )              ACV10069
C                                                                       ACV10070
          call tchk( label, irmtree, *100 )                             ACV10071
            iloc = irmtree(-5 )                                         ACV10072
            irmtree( iloc + irmtree(-2) )     = 1                       ACV10073
            irmtree( iloc + irmtree(-4) + 1 ) = irme                    ACV10074
            drme( irme ) = dme                                          ACV10075
            call tins( label, irmtree )                                 ACV10076
100       continue                                                      ACV10077
          inxold = indu3c                                               ACV10078
        end if                                                          ACV10079
      end do                                                            ACV10080
1030  format(14x, i2,' (',i2,i3,')',i2, '(x)', i2,' (',i2,i3,')',i2,    ACV10081
     1       ';', i2, '(',2i3,')',i2, 10x, i2, f18.10)                  ACV10082
1040  format(42x, i2, '(',2i3,')',i2, 10x, i2, f18.10)                  ACV10083
CW----------------------------------------------------------------------ACV10084
                  ix = ix + i3size                                      ACV10085
                end do ! kro12                                          ACV10086
              end do ! ist                                              ACV10087
              i3first = i3first + i3size                                ACV10088
              ir = ir + i3size * kromax                                 ACV10089
            end if                                                      ACV10090
          end do ! i                                                    ACV10091
        end do ! iar                                                    ACV10092
      end do ! ial                                                      ACV10093
      write(6,*)                                                        ACV10094
C     Reset iqops( * )                                                  ACV10095
      do i = 1, ilm                                                     ACV10096
        iqops( 2, i ) = ishft( ishft( iqops(2,i),-8 ), 8 )              ACV10097
        i3ptr( i )    = 0                                               ACV10098
      end do                                                            ACV10099
      return                                                            ACV10100
999   return 1                                                          ACV10101
C ---*end of rm3bar*----------------------------------------------------ACV10102
      end                                                               ACV10103
      SUBROUTINE DLUT(MA,NA,IA,DA,MD)                                   ACV10104
C     ------------------------------------------------------------------ACV10105
C     DECOMPOSITION OF A REAL MATRIX INTO THE PRODUCT OF A LOWER AND    ACV10106
C     AN UPPER TRIANGULAR MATRIX                                        ACV10107
C     ------------------------------------------------------------------ACV10108
C     REFERENCES--COMPUTER SOLUTION OF LINEAR ALGEBRAIC SYSTEMS,        ACV10109
C                 G.FORSYTHE AND C.MOLER, PRENTICE HALL                 ACV10110
C     PARAMETERS--MA  ACTUAL NUMBER OF ROWS                             ACV10111
C                 NA  ACTUAL NUMBER OF COLUMNS                          ACV10112
C                 MD  DIMENSIONED NUMBER OF ROWS                        ACV10113
C                 ND  DIMENSIONED NUMBER OF COLUMNS                     ACV10114
C                 DA  MATRIX TO BE DECOMPOSED                           ACV10115
C                 IA  INPUT AS IA(I)=I RETURNED IN SORTED ORDER         ACV10116
C     DIMENSIONS--DA(MD*ND),IA(MD)                                      ACV10117
C     ------------------------------------------------------------------ACV10118
      IMPLICIT REAL*8(D)                                                ACV10119
      DIMENSION DA(1),IA(1)                                             ACV10120
      LA=MIN0(MA,NA)                                                    ACV10121
      LQ=-MD                                                            ACV10122
      DO 35 L=1,LA                                                      ACV10123
      IF(L.LT.LA)GO TO 10                                               ACV10124
      IF(MA.LE.NA)RETURN                                                ACV10125
   10 LQ=LQ+MD                                                          ACV10126
      DBIG=0.D0                                                         ACV10127
      DO 20 I=L,MA                                                      ACV10128
      IL=I+LQ                                                           ACV10129
      IF(DABS(DA(IL))-DABS(DBIG))20,20,15                               ACV10130
   15 IBIG=I                                                            ACV10131
      DBIG=DA(IL)                                                       ACV10132
   20 CONTINUE                                                          ACV10133
      ISAVE=IA(L)                                                       ACV10134
      IA(L)=IA(IBIG)                                                    ACV10135
      IA(IBIG)=ISAVE                                                    ACV10136
      JQ=-MD                                                            ACV10137
      DO 25 J=1,NA                                                      ACV10138
      JQ=JQ+MD                                                          ACV10139
      LJ=L+JQ                                                           ACV10140
      IBIGJ=IBIG+JQ                                                     ACV10141
      DSAVE=DA(LJ)                                                      ACV10142
      DA(LJ)=DA(IBIGJ)                                                  ACV10143
   25 DA(IBIGJ)=DSAVE                                                   ACV10144
      K=L+1                                                             ACV10145
      DO 30 I=K,MA                                                      ACV10146
      IL=I+LQ                                                           ACV10147
   30 DA(IL)=DA(IL)/DBIG                                                ACV10148
      IF(L.EQ.LA)RETURN                                                 ACV10149
      DO 35 I=K,MA                                                      ACV10150
      IL=I+LQ                                                           ACV10151
      JQ=LQ                                                             ACV10152
      DO 35 J=K,NA                                                      ACV10153
      JQ=JQ+MD                                                          ACV10154
      IJ=I+JQ                                                           ACV10155
      LJ=L+JQ                                                           ACV10156
   35 DA(IJ)=DA(IJ)-DA(IL)*DA(LJ)                                       ACV10157
      RETURN                                                            ACV10158
      END                                                               ACV10159
      SUBROUTINE DBSR(MA,DA,DB,DC,IA,MD)                                ACV10160
C                                                                       ACV10161
C Modified by C.Bahri on January 13, 1991: including IA in the subrouti-ACV10162
C ne, rather than doing it externally.                                  ACV10163
C                                                                       ACV10164
C Modified by EJR on September 1, 1985.  The special case for MA = 1 hasACV10165
C been included in this subroutine, rather than done externally as in   ACV10166
C the original version of CRU3(...).                                    ACV10167
C                                                                       ACV10168
C     ------------------------------------------------------------------ACV10169
C     DOUBLE BACK SUBSTITUTION FOR SOLVING SIMULTANEOUS EQUATIONS       ACV10170
C     ------------------------------------------------------------------ACV10171
C     REFERENCES--COMPUTER SOLUTION OF LINEAR ALGEBRAIC SYSTEMS,        ACV10172
C                 G.FORSYTHE AND C.MOLER, PRENTICE HALL                 ACV10173
C     PARAMETERS--MA  ACTUAL NUMBER OF ROWS                             ACV10174
C                 MD  DIMENSIONED NUMBER OF ROWS                        ACV10175
C                 ND  DIMENSIONED NUMBER OF COLUMNS                     ACV10176
C                 DA  MATRIX OF COEFFICIENTS PREPARED IN DLUT           ACV10177
C                 DB  INPUT CONSTANTS REPLACED BY SOLUTIONS             ACV10178
C                 DC  WORK VECTOR                                       ACV10179
C     DIMENSIONS--DA(MD*ND),DB(ND),DC(ND)                               ACV10180
C     ------------------------------------------------------------------ACV10181
      IMPLICIT REAL*8(D)                                                ACV10182
      DIMENSION DA(1),DB(1),DC(1),IA(1)                                 ACV10183
      IF( MA.GT.1 ) GOTO 5                                              ACV10184
      DB(1) = DB(IA(1))/DA(1)                                           ACV10185
      GOTO 30                                                           ACV10186
5     DC(1)=DB(IA(1))                                                   ACV10187
      DO 15 I=2,MA                                                      ACV10188
      NA=I-1                                                            ACV10189
      DSUM=0.D0                                                         ACV10190
      JQ=-MD                                                            ACV10191
      DO 10 J=1,NA                                                      ACV10192
      JQ=JQ+MD                                                          ACV10193
      IJ=I+JQ                                                           ACV10194
   10 DSUM=DSUM+DA(IJ)*DC(J)                                            ACV10195
   15 DC(I)=DB(IA(I))-DSUM                                              ACV10196
      MAMA=IJ+MD                                                        ACV10197
      DB(MA)=DC(MA)/DA(MAMA)                                            ACV10198
      IQ=NA*MD                                                          ACV10199
      DO 25 IP=1,NA                                                     ACV10200
      IQ=IQ-MD                                                          ACV10201
      I=MA-IP                                                           ACV10202
      II=I+IQ                                                           ACV10203
      DSUM=0.D0                                                         ACV10204
      JQ=MA*MD                                                          ACV10205
      DO 20 JP=1,IP                                                     ACV10206
      JQ=JQ-MD                                                          ACV10207
      J=MA-JP+1                                                         ACV10208
      IJ=I+JQ                                                           ACV10209
   20 DSUM=DSUM+DA(IJ)*DB(J)                                            ACV10210
   25 DB(I)=(DC(I)-DSUM)/DA(II)                                         ACV10211
30    RETURN                                                            ACV10212
      END                                                               ACV10213
      subroutine piksrt( nbits, nbitx, npart, level, todd )             ACV10214
C                                                                       ACV10215
C     Sort particles                                                    ACV10216
C                                                                       ACV10217
      logical todd                                                      ACV10218
      dimension lvla( 8 )                                               ACV10219
C                                                                       ACV10220
      todd = .false.                                                    ACV10221
      do i = npart, 1, -1                                               ACV10222
        lvl       = iand( level, nbitx )                                ACV10223
        lvla( i ) = lvl                                                 ACV10224
        level     = ishft( level, -nbits )                              ACV10225
      end do                                                            ACV10226
C     Straight insertion ( Num. Recipes, p.227 )                        ACV10227
      do j = 2, npart           ! Pick out each element in turn.        ACV10228
        lvl = lvla( j )                                                 ACV10229
        do i = j-1, 1, -1       ! Look for the place to insert.         ACV10230
          if( lvla( i ) .le. lvl ) go to 10                             ACV10231
          todd = .not. todd                                             ACV10232
          lvla( i+1 ) = lvla( i )                                       ACV10233
        end do                                                          ACV10234
        i = 0                                                           ACV10235
10      lvla( i+1 ) = lvl                                               ACV10236
      end do                                                            ACV10237
      level = lvla( 1 )                                                 ACV10238
      do i = 2, npart                                                   ACV10239
        level = ior( ishft( level, nbits ), lvla(i) )                   ACV10240
      end do                                                            ACV10241
      return                                                            ACV10242
C----*end of piksrt*--------------------------------------------------- ACV10243
      end                                                               ACV10244
CB----------------------------------------------------------------------ACV10245
C                          ***************                              ACV10246
C                          **  UTILITY ***                              ACV10247
C                          ***************                              ACV10248
C ----------------------------------------------------------------------ACV10249
      SUBROUTINE SETBIN(LFILX,TWRITX,BINX,NEMPTX,NFIRSTX,NMAXX,NTOTALX) ACV10250
C     Generating the binary I/O.                                        ACV10251
C                                                                       ACV10252
      COMMON /SAVARG/ LFILE, TWRITE, NEMPTY, NFIRST, NMAX, NTOTAL       ACV10253
      COMMON /SAVCHR/ BIN                                               ACV10254
      COMMON /SAVINT/ IPM, NMAXP1, NMAX2, NMAXP2, NMAX3                 ACV10255
      LOGICAL TWRITE, TWRITX                                            ACV10256
      CHARACTER*(*) BINX                                                ACV10257
      CHARACTER BIN*115                                                 ACV10258
      LFILE = LFILX                                                     ACV10259
      TWRITE = TWRITX                                                   ACV10260
      NEMPTY = NEMPTX                                                   ACV10261
      NFIRST = NFIRSTX                                                  ACV10262
      NMAX = NMAXX                                                      ACV10263
      NTOTAL = NTOTALX                                                  ACV10264
      NMAXP1=NMAX+6                                                     ACV10265
      NMAX2=NMAX+NMAX+5                                                 ACV10266
      NMAXP2=NMAX2+6                                                    ACV10267
      NMAX3=3*NMAX+10                                                   ACV10268
      IF(NFIRST.LE.0)THEN                                               ACV10269
        IPM=1                                                           ACV10270
      ELSE                                                              ACV10271
        IPM=-1                                                          ACV10272
      ENDIF                                                             ACV10273
      DO J=1,115                                                        ACV10274
        BIN(J:J)=' '                                                    ACV10275
      ENDDO                                                             ACV10276
      RETURN                                                            ACV10277
C                                                                       ACV10278
      ENTRY DBOVRL(ILEFT,IL,IR,IRGHT)                                   ACV10279
C     Converting a number from decimal to binary. (1a)                  ACV10280
C-----------------------------------------------------------------------ACV10281
C     Reading bits in number.                                           ACV10282
      I=IPM                                                             ACV10283
      DO J=1,NMAX                                ! LEFT                 ACV10284
        IF(MOD(J,NEMPTY).NE.0)THEN                                      ACV10285
           IF(BTEST(ILEFT,NFIRST+I))THEN                                ACV10286
              BIN(J:J)='1'                                              ACV10287
           ELSE                                                         ACV10288
              BIN(J:J)='0'                                              ACV10289
           ENDIF                                                        ACV10290
           I=I+IPM                                                      ACV10291
        ENDIF                                                           ACV10292
      ENDDO                                                             ACV10293
                  IF(NMAX.EQ.NTOTAL) THEN                               ACV10294
C     Writing number in binary form.                                    ACV10295
      IF(TWRITE) THEN                                                   ACV10296
C       WRITE(LFILE,100) ILEFT,BIN                                      ACV10297
        write(lfile,100) ileft, (bin(i:i), i=1, j)                      ACV10298
      ELSE                                                              ACV10299
        DO J = 1, NTOTAL                                                ACV10300
          BINX(J:J) = BIN(J:J)                                          ACV10301
        END DO                                                          ACV10302
      END IF                                                            ACV10303
                  ELSE                                                  ACV10304
      I=IPM                                                             ACV10305
      NMIN=MOD(NMAXP1-1,NEMPTY)                                         ACV10306
      DO J=NMAXP1,NMAX2                          ! MIDDLE               ACV10307
        IF(MOD(J,NEMPTY).NE.NMIN)THEN                                   ACV10308
           NF=NFIRST+I                                                  ACV10309
           IF(BTEST(IL,NF))THEN                                         ACV10310
              IF(BTEST(IR,NF))THEN                                      ACV10311
                 BIN(J:J)='I'                                           ACV10312
              ELSE                                                      ACV10313
                 BIN(J:J)='V'                                           ACV10314
              ENDIF                                                     ACV10315
           ELSE                                                         ACV10316
              IF(BTEST(IR,NF))THEN                                      ACV10317
                 BIN(J:J)='A'                                           ACV10318
              ELSE                                                      ACV10319
                 BIN(J:J)='-'                                           ACV10320
              ENDIF                                                     ACV10321
           ENDIF                                                        ACV10322
           I=I+IPM                                                      ACV10323
        ENDIF                                                           ACV10324
      ENDDO                                                             ACV10325
      I=IPM                                                             ACV10326
      NMIN=MOD(NMAXP2-1,NEMPTY)                                         ACV10327
      DO J=NMAXP2,NMAX3                          ! RIGHT                ACV10328
        IF(MOD(J,NEMPTY).NE.NMIN)THEN                                   ACV10329
           IF(BTEST(IRGHT,NFIRST+I))THEN                                ACV10330
              BIN(J:J)='1'                                              ACV10331
           ELSE                                                         ACV10332
              BIN(J:J)='0'                                              ACV10333
           ENDIF                                                        ACV10334
           I=I+IPM                                                      ACV10335
        ENDIF                                                           ACV10336
      ENDDO                                                             ACV10337
C     Writing number in binary form.                                    ACV10338
C     IF(TWRITE) WRITE(LFILE,'(T2,A)') BIN                              ACV10339
      if( twrite ) write(lfile,'(t2,115a)') (bin(i:i), i=1, j)          ACV10340
                  ENDIF                                                 ACV10341
      RETURN                                                            ACV10342
C                                                                       ACV10343
      ENTRY DBOVRP(NBITS,ILEFT,IL,IR,IRGHT)                             ACV10344
C     Converting a number from decimal to binary. (1b)                  ACV10345
C-----------------------------------------------------------------------ACV10346
C     Reading bits in number.                                           ACV10347
      I=IPM                                                             ACV10348
      J=0                                                               ACV10349
      DO WHILE(J.LE.NMAX .AND. IABS(I).LE.NBITS)   ! LEFT               ACV10350
        J=J+1                                                           ACV10351
        IF(MOD(J,NEMPTY).NE.0) THEN                                     ACV10352
           IF(BTEST(ILEFT,NFIRST+I))THEN                                ACV10353
              BIN(J:J)='1'                                              ACV10354
           ELSE                                                         ACV10355
              BIN(J:J)='0'                                              ACV10356
           ENDIF                                                        ACV10357
           I=I+IPM                                                      ACV10358
        ENDIF                                                           ACV10359
      ENDDO                                                             ACV10360
      I=IPM                                                             ACV10361
      J=NMAXP1-1                                                        ACV10362
      NMIN=MOD(J,NEMPTY)                                                ACV10363
      DO WHILE(J.LE.NMAX2 .AND. IABS(I).LE.NBITS)   ! MIDDLE            ACV10364
        J=J+1                                                           ACV10365
        IF(MOD(J,NEMPTY).NE.NMIN) THEN                                  ACV10366
           NF=NFIRST+I                                                  ACV10367
           IF(BTEST(IL,NF))THEN                                         ACV10368
              IF(BTEST(IR,NF))THEN                                      ACV10369
                 BIN(J:J)='I'                                           ACV10370
              ELSE                                                      ACV10371
                 BIN(J:J)='A'                                           ACV10372
              ENDIF                                                     ACV10373
           ELSE                                                         ACV10374
              IF(BTEST(IR,NF))THEN                                      ACV10375
                 BIN(J:J)='V'                                           ACV10376
              ELSE                                                      ACV10377
                 BIN(J:J)='-'                                           ACV10378
              ENDIF                                                     ACV10379
           ENDIF                                                        ACV10380
           I=I+IPM                                                      ACV10381
        ENDIF                                                           ACV10382
      ENDDO                                                             ACV10383
      I=IPM                                                             ACV10384
      J=NMAXP2-1                                                        ACV10385
      NMIN=MOD(J,NEMPTY)                                                ACV10386
      DO WHILE(J.LE.NMAX3 .AND. IABS(I).LE.NBITS)   ! RIGHT             ACV10387
        J=J+1                                                           ACV10388
        IF(MOD(J,NEMPTY).NE.NMIN) THEN                                  ACV10389
           IF(BTEST(IRGHT,NFIRST+I))THEN                                ACV10390
              BIN(J:J)='1'                                              ACV10391
           ELSE                                                         ACV10392
              BIN(J:J)='0'                                              ACV10393
           ENDIF                                                        ACV10394
           I=I+IPM                                                      ACV10395
        ENDIF                                                           ACV10396
      ENDDO                                                             ACV10397
C     Writing number in binary form.                                    ACV10398
      IF(TWRITE) THEN                                                   ACV10399
C       WRITE(LFILE,'(T2,A)') BIN                                       ACV10400
        write(lfile,'(t2,115a)') (bin(i:i), i=1, j)                     ACV10401
      ELSE                                                              ACV10402
        DO J = 1, NTOTAL                                                ACV10403
          BINX(J:J) = BIN(J:J)                                          ACV10404
        END DO                                                          ACV10405
      END IF                                                            ACV10406
      RETURN                                                            ACV10407
C                                                                       ACV10408
      ENTRY DBOVRC(NBITS,ILEFT,IL,IR,IRGHT)                             ACV10409
C     Converting a number from decimal to binary. (2)                   ACV10410
C-----------------------------------------------------------------------ACV10411
C     Reading bits in number.                                           ACV10412
      I=IPM                                                             ACV10413
      J=0                                                               ACV10414
      DO WHILE(J.LE.NMAX .AND. IABS(I).LE.NBITS)   ! LEFT               ACV10415
        J=J+1                                                           ACV10416
        IF(MOD(J,NEMPTY).NE.0) THEN                                     ACV10417
           IF(BTEST(ILEFT,NFIRST+I))THEN                                ACV10418
              BIN(J:J)='1'                                              ACV10419
           ELSE                                                         ACV10420
              BIN(J:J)='0'                                              ACV10421
           ENDIF                                                        ACV10422
           I=I+IPM                                                      ACV10423
        ENDIF                                                           ACV10424
      ENDDO                                                             ACV10425
                  IF(NMAX.EQ.NTOTAL) THEN                               ACV10426
C     Writing number in binary form.                                    ACV10427
      IF(TWRITE) THEN                                                   ACV10428
C       WRITE(LFILE,100) ILEFT,BIN                                      ACV10429
        write(lfile,100) ileft, (bin(i:i), i=1, j)                      ACV10430
      ELSE                                                              ACV10431
        DO J = 1, NTOTAL                                                ACV10432
          BINX(J:J) = BIN(J:J)                                          ACV10433
        END DO                                                          ACV10434
      END IF                                                            ACV10435
                  ELSE                                                  ACV10436
      NMAXP1=NBITS+6                                                    ACV10437
      NMAX2=NBITS+NBITS+5                                               ACV10438
      NMAXP2=NMAX2+6                                                    ACV10439
      NMAX3=3*NBITS+10                                                  ACV10440
      I=IPM                                                             ACV10441
      J=NMAXP1-1                                                        ACV10442
      NMIN=MOD(J,NEMPTY)                                                ACV10443
      DO WHILE(J.LE.NMAX2 .AND. IABS(I).LE.NBITS)   ! MIDDLE            ACV10444
        J=J+1                                                           ACV10445
        IF(MOD(J,NEMPTY).NE.NMIN) THEN                                  ACV10446
           NF=NFIRST+I                                                  ACV10447
           IF(BTEST(IL,NF))THEN                                         ACV10448
              IF(BTEST(IR,NF))THEN                                      ACV10449
                 BIN(J:J)='I'                                           ACV10450
              ELSE                                                      ACV10451
                 BIN(J:J)='A'                                           ACV10452
              ENDIF                                                     ACV10453
           ELSE                                                         ACV10454
              IF(BTEST(IR,NF))THEN                                      ACV10455
                 BIN(J:J)='V'                                           ACV10456
              ELSE                                                      ACV10457
                 BIN(J:J)='-'                                           ACV10458
              ENDIF                                                     ACV10459
           ENDIF                                                        ACV10460
           I=I+IPM                                                      ACV10461
        ENDIF                                                           ACV10462
      ENDDO                                                             ACV10463
      I=IPM                                                             ACV10464
      J=NMAXP2-1                                                        ACV10465
      NMIN=MOD(J,NEMPTY)                                                ACV10466
      DO WHILE(J.LE.NMAX3 .AND. IABS(I).LE.NBITS)   ! RIGHT             ACV10467
        J=J+1                                                           ACV10468
        IF(MOD(J,NEMPTY).NE.NMIN) THEN                                  ACV10469
           IF(BTEST(IRGHT,NFIRST+I))THEN                                ACV10470
              BIN(J:J)='1'                                              ACV10471
           ELSE                                                         ACV10472
              BIN(J:J)='0'                                              ACV10473
           ENDIF                                                        ACV10474
           I=I+IPM                                                      ACV10475
        ENDIF                                                           ACV10476
      ENDDO                                                             ACV10477
C     Writing number in binary form.                                    ACV10478
C     IF(TWRITE) WRITE(LFILE,'(T2,A)') BIN                              ACV10479
      if( twrite ) write(lfile,'(t2,115a)') (bin(i:i), i=1, j)          ACV10480
                  ENDIF                                                 ACV10481
C100  FORMAT(' ',I23,2X,A)                                              ACV10482
100   format(' ',i23,2x,115a)                                           ACV10483
      RETURN                                                            ACV10484
C ---*end of SETBIN*----------------------------------------------------ACV10485
      END                                                               ACV10486
C                                                                       ACV10487
      SUBROUTINE ATTN(LITER,*)                                          ACV10488
C     Warning because of gross error.                                   ACV10489
      CHARACTER*(*) LITER                                               ACV10490
      WRITE(6,'(//A,A)') ' ***** ATTENTION ',LITER                      ACV10491
      RETURN 1                                                          ACV10492
C     RETURN                                                            ACV10493
C                                                                       ACV10494
      ENTRY ERROR(LITER)                                                ACV10495
C     Exit program because of gross error.                              ACV10496
      WRITE(6,'(//A,A)') ' ***** FATAL ERROR ',LITER                    ACV10497
      STOP                                                              ACV10498
C ---*end of ATTN*------------------------------------------------------ACV10499
      END                                                               ACV10500
C                                                                       ACV10501
      SUBROUTINE WRDATA(IOTYPE,IOFILE,NUMI,IRAY,NF)                     ACV10502
C                                                                       ACV10503
C     ORIGINAL: JPD ('PHDRYR.RESKE2.FORT(IODATA)')                      ACV10504
C     MODIFIED: CB (2/15/91) FOR GENERAL PURPOSE                        ACV10505
C                                                                       ACV10506
CB    IMPLICIT REAL*8(A-H,O-Z)                                          ACV10507
      DIMENSION IRAY(NF:1)                                              ACV10508
      IF (IOTYPE.EQ.0) THEN                                             ACV10509
         WRITE(IOFILE)NUMI                                              ACV10510
         WRITE(IOFILE)(IRAY(I),I=NF,0)                                  ACV10511
      ELSE                                                              ACV10512
         READ(IOFILE)NUMI                                               ACV10513
         READ(IOFILE)(IRAY(I),I=NF,0)                                   ACV10514
      ENDIF                                                             ACV10515
C     SPLIT: LRECL=1724 & BLKSIZE=32760 --> 430 4 BYTE SETS PER RECORD  ACV10516
C            (BUFFER SIZE IS 32760 SO SPLIT INTO 19 RECORDS WITH        ACV10517
C             4 BYTES PER RECORD PLUS 4 BYTES PER BLOCK OVERHEAD:       ACV10518
C             19*(1720+4)+4=32760 FOR OPTIMUM BUFFER UTILIZATION)       ACV10519
      NRUN=NUMI                                                         ACV10520
      NREC=NRUN/430                                                     ACV10521
      NBEG=1                                                            ACV10522
      NEND=430                                                          ACV10523
C     I/O FULL BLOCKS, INTEGER ARRAY                                    ACV10524
      DO 10 NR=1,NREC                                                   ACV10525
         IF (IOTYPE.EQ.0) THEN                                          ACV10526
            WRITE(IOFILE)(IRAY(N),N=NBEG,NEND)                          ACV10527
         ELSE                                                           ACV10528
            READ(IOFILE)(IRAY(N),N=NBEG,NEND)                           ACV10529
         ENDIF                                                          ACV10530
         NBEG=NEND+1                                                    ACV10531
10       NEND=NEND+430                                                  ACV10532
C     I/O RESIDUAL                                                      ACV10533
      IF (NBEG.LE.NRUN) THEN                                            ACV10534
         IF (IOTYPE.EQ.0) THEN                                          ACV10535
            WRITE(IOFILE)(IRAY(N),N=NBEG,NRUN)                          ACV10536
         ELSE                                                           ACV10537
            READ(IOFILE)(IRAY(N),N=NBEG,NRUN)                           ACV10538
         ENDIF                                                          ACV10539
      ENDIF                                                             ACV10540
      RETURN                                                            ACV10541
      END                                                               ACV10542
      SUBROUTINE RATIONAL(VALUE,INUMER,IDENOM,IFIND)                    ACV10543
C                                                                       ACV10544
C        This subroutine determines the rational form of a real number  ACV10545
C     (e.g. 11/17) provided the lesser of the numerator or denominator  ACV10546
C     is less than or equal to NTIMES. The method used is a simple      ACV10547
C     brute force multiplication search for an integer numerator or     ACV10548
C     denominator.                                                      ACV10549
C                                                                       ACV10550
      IMPLICIT REAL*8(A-H,O-Z)                                          ACV10551
      IMPLICIT INTEGER*4(I-N)                                           ACV10552
      PARAMETER(NTIMES=5000,TOLER=1.0D-5)                               ACV10553
C                                                                       ACV10554
      IF(DABS(VALUE) .LT. TOLER) THEN                                   ACV10555
         INUMER=0                                                       ACV10556
         IDENOM=0                                                       ACV10557
         RETURN                                                         ACV10558
      END IF                                                            ACV10559
      IFIND=0                                                           ACV10560
      IF(VALUE .GE. 1.D0) THEN                                          ACV10561
         I=0                                                            ACV10562
         DO WHILE(I .LE. NTIMES)                                        ACV10563
            I=I+1                                                       ACV10564
            VALNEW=VALUE*DFLOAT(I)                                      ACV10565
            VALUP=IDNINT(VALNEW)                                        ACV10566
            IF(DABS(VALUP-VALNEW) .LT. TOLER) VALNEW=VALUP              ACV10567
            IF((VALNEW-IDINT(VALNEW)) .LT. TOLER) THEN                  ACV10568
              INUMER=IDINT(VALNEW)                                      ACV10569
              IDENOM=I                                                  ACV10570
              IFIND=1                                                   ACV10571
              RETURN                                                    ACV10572
            END IF                                                      ACV10573
         END DO                                                         ACV10574
      ELSE                                                              ACV10575
         I=0                                                            ACV10576
         DO WHILE(I .LE. NTIMES)                                        ACV10577
            I=I+1                                                       ACV10578
            VALNEW=DFLOAT(I)/VALUE                                      ACV10579
            VALUP=IDNINT(VALNEW)                                        ACV10580
            IF(DABS(VALUP-VALNEW) .LT. TOLER) VALNEW=VALUP              ACV10581
            IF(DABS(VALNEW-IDINT(VALNEW)) .LT. TOLER) THEN              ACV10582
              INUMER=I                                                  ACV10583
              IDENOM=IDINT(VALNEW)                                      ACV10584
              IFIND=1                                                   ACV10585
              RETURN                                                    ACV10586
            END IF                                                      ACV10587
         END DO                                                         ACV10588
      END IF                                                            ACV10589
      RETURN                                                            ACV10590
      END                                                               ACV10591
C ----------------------------------------------------------------------ACV10592
C                                                                       ACV10593
C                       ***********************                         ACV10594
C                       ***   WST PACKAGE   ***                         ACV10595
C                       ***********************                         ACV10596
C                                                                       ACV10597
C                -------------------------------------                  ACV10598
C                                                                       ACV10599
C             *******************************************               ACV10600
C             ***    WEIGHTED SEARCH TREE ROUTINES    ***               ACV10601
C             ***                for                  ***               ACV10602
C             ***  SCIENTIFIC (FORTRAN) APPLICATIONS  ***               ACV10603
C             *******************************************               ACV10604
C                                                                       ACV10605
C ----------------------------------------------------------------------ACV10606
C                                                                       ACV10607
C Authors: Soon Park, C. Bahri, J. P. Draayer and S.-Q. Zheng           ACV10608
C          Department of Physics (Computer Science)                     ACV10609
C          Louisiana State University                                   ACV10610
C          Baton Rouge LA                                               ACV10611
C          USA 70803-4001                                               ACV10612
C                                                                       ACV10613
C                 BITNET:  PHDRYR @ LSUMVS or LSUVM                     ACV10614
C                 TELEX:   559184                                       ACV10615
C                 PHONE:   USA-504-388-2261                             ACV10616
C                 FAX:     USA-504-388-5855                             ACV10617
C                                                                       ACV10618
C Version: 1.1    LSU (07/01/91)                                        ACV10619
C Version: 2.1    LSU (04/01/94)                                        ACV10620
C                                                                       ACV10621
C ----------------------------------------------------------------------ACV10622
C                                                                       ACV10623
C Updates: 07/90 Original from a FORTRAN code written by Soon Park      ACV10624
C          04/94 Inclusion of TMRG subroutine by C. Bahri               ACV10625
C                                                                       ACV10626
C ----------------------------------------------------------------------ACV10627
C                                                                       ACV10628
C General comments on the package:                                      ACV10629
C                                                                       ACV10630
C   This package is written in FORTRAN since most scientific programs   ACV10631
C   require FORTRAN compatibility and for it the existing scientific    ACV10632
C   subroutine libraries are the most extensive and efficient.          ACV10633
C                                                                       ACV10634
C   The tree is a linear array ID(-10:*) consisting of eleven integers, ACV10635
C   ID(-10:0), that specifies the structure of the tree, see TSET for   ACV10636
C   documentation on this, and node information starting with ID(1:1).  ACV10637
C   The latter includes the key(s), data, balance factor and priority,  ACV10638
C   as well as the left and right child pointers. See the documentation ACV10639
C   on each subroutine for further details.                             ACV10640
C                                                                       ACV10641
C ----------------------------------------------------------------------ACV10642
C                                                                       ACV10643
C This numerical database package consists of 7 different subroutines:  ACV10644
C                                                                       ACV10645
C   1. TSET    -->  initializes a storage area for use as a binary tree ACV10646
C       TSETLL -->  ... entry in TSET for a tree with a linked-list     ACV10647
C       TSETLF -->  ... entry in TSET for a tree without a linked-list  ACV10648
C   2. TCHK    -->  performs the search operation on the data structure ACV10649
C   3. TADD    -->  add a new element to an existing WST data structure ACV10650
C   4. TINS    -->  called by TADD to insert new node into existing WST ACV10651
C   5. TDEL    -->  called by TADD to delete low priority node in a WST ACV10652
C   6. TOUT    -->  generates output information on specific tree nodes ACV10653
C   7. TMRG    -->  merges two trees that have the same structure       ACV10654
C                                                                       ACV10655
C The four routines TSET, TCHK, TADD and TOUT subject to user control   ACV10656
C while TINS and TDEL are only used by TADD and not otherwise needed.   ACV10657
C The type of information stored in the buffer may vary from applicaton ACV10658
C to application. This can be achieved by editing TADD. For example,    ACV10659
C if variable length integer data is to be stored, then BUFFER must be  ACV10660
C replaced by, INTGER, an integer array throughout. Likewise if BUFFER  ACV10661
C holds double precision or complex data, the statement that defines    ACV10662
C BUFFER in TADD must be changed accordingly. The program can also be   ACV10663
C used in other ways, for example, each node can refer to more than a   ACV10664
C single buffer. If this is done, both TSET and TADD must be modified   ACV10665
C in a rather obvious way to add the required multiple link-list and    ACV10666
C buffer arrays.                                                        ACV10667
C                                                                       ACV10668
C ----------------------------------------------------------------------ACV10669
C                                                                       ACV10670
C                           **************                              ACV10671
C                           ***  TSET  ***                              ACV10672
C                           **************                              ACV10673
C                                                                       ACV10674
C The subroutine TSET must be called before inserting the first item    ACV10675
C into the tree. This call fixes the first 11 values of the array ID:   ACV10676
C                                                                       ACV10677
C      ID(-10): Number of nodes currently in the tree    ==> ID(-10)    ACV10678
C      ID(-9) : Maximum nodes in in tree                 ==> MXNODE     ACV10679
C      ID(-8) : Parent node of ID(-7), see next entry    ==> NF         ACV10680
C      ID(-7) : Node to be balanced                      ==> NA         ACV10681
C      ID(-6) : Parent node pointer                      ==> NQ         ACV10682
C      ID(-5) : Current node pointer                     ==> Determined ACV10683
C      ID(-4) : Number of integers assigned the key      ==> NKEY       ACV10684
C      ID(-3) : Number of integers for key and data      ==> NSUM       ACV10685
C      ID(-2) : Position of the priority in a node       ==> NPR        ACV10686
C      ID(-1) : Position of the next available node      ==> Determined ACV10687
C      ID( 0) : Root pointer (-1 for empty tree)         ==> Determined ACV10688
C                                                                       ACV10689
C The subroutine TSET also initializes all entries of the array LLBUFF  ACV10690
C where the link-list information for BUFFER is stored. In particular,  ACV10691
C the first three which refer to the free space are assigned values:    ACV10692
C                                                                       ACV10693
C      LLBUFF(-2) : Tail of free space                                  ACV10694
C      LLBUFF(-1) : Head of free space                                  ACV10695
C      LLBUFF( 0) : Size of free space                                  ACV10696
C                                                                       ACV10697
C ----------------------------------------------------------------------ACV10698
C                                                                       ACV10699
      SUBROUTINE TSET                                                   ACV10700
C                                                                       ACV10701
      INTEGER ID(-10:*),LLBUFF(-2:*)                                    ACV10702
C                                                                       ACV10703
      ENTRY TSETLL(ID,MXNODE,NKEY,NDAT,LLBUFF,MXBUFF)                   ACV10704
C                                                                       ACV10705
C Initialize link-list array values                                     ACV10706
C                                                                       ACV10707
      DO 10 J=1,MXBUFF-1                                                ACV10708
10       LLBUFF(J)=J+1                                                  ACV10709
      LLBUFF(MXBUFF)=-1                                                 ACV10710
C                                                                       ACV10711
C Initialize free space pointer/counter                                 ACV10712
C                                                                       ACV10713
      LLBUFF(-2)=MXBUFF                                                 ACV10714
      LLBUFF(-1)=1                                                      ACV10715
      LLBUFF(0)=MXBUFF                                                  ACV10716
C                                                                       ACV10717
      ENTRY TSETLF(ID,MXNODE,NKEY,NDAT)                                 ACV10718
C                                                                       ACV10719
C Initialize WST parameters                                             ACV10720
C                                                                       ACV10721
      NSUM=NKEY+NDAT                                                    ACV10722
      IMAX=MXNODE*(NSUM+3)                                              ACV10723
      ID(-10)=0                                                         ACV10724
      ID(-9)=MXNODE                                                     ACV10725
      ID(-8)=-1                                                         ACV10726
      ID(-7)=-1                                                         ACV10727
      ID(-6)=-1                                                         ACV10728
      ID(-5)=0                                                          ACV10729
      ID(-4)=NKEY                                                       ACV10730
      ID(-3)=NSUM                                                       ACV10731
      ID(-2)=NSUM+1                                                     ACV10732
      ID(-1)=0                                                          ACV10733
      ID(0)=-1                                                          ACV10734
      DO 30 I=1,IMAX                                                    ACV10735
30       ID(I)=0                                                        ACV10736
      RETURN                                                            ACV10737
      END                                                               ACV10738
C ----------------------------------------------------------------------ACV10739
C                                                                       ACV10740
C                            ************                               ACV10741
C                            *** TCHK ***                               ACV10742
C                            ************                               ACV10743
C                                                                       ACV10744
C The subroutine TCHK constructs a weighted search tree or locates a    ACV10745
C specific node in a weighted search tree. The subroutine uses two      ACV10746
C arrays, NEWKEY and ID. A call to TCHK(,,) generates a search of ID forACV10747
C the key stored in NEWKEY. If the search is successful the priority of ACV10748
C the node will be increased. If not successful, a normal return is     ACV10749
C generated and the position, ID(-5) where a new node can be inserted   ACV10750
C is given. After this subroutine is called, the position(s) where the  ACV10751
C key(s) of the new node are to be assigned are ID(-5)+1, ID(-5)+2,     ACV10752
C ID(-5)+3, ... and the positions of the location and size of the new   ACV10753
C incoming data are ID(IPOS)+1 and ID(IPOS)+2, respectively, where      ACV10754
C IPOS = ID(-5)+ID(-4) and ID(-4) is the number of integer words set    ACV10755
C aside for the key.                                                    ACV10756
C                                                                       ACV10757
C   RETURN 1 --> Successful search (key already in the tree,            ACV10758
C                ID(-5) is set to point to the located node).           ACV10759
C                The priority of the located node is updated.           ACV10760
C                                                                       ACV10761
C ----------------------------------------------------------------------ACV10762
C                                                                       ACV10763
      SUBROUTINE TCHK(NEWKEY,ID,*)                                      ACV10764
C                                                                       ACV10765
      INTEGER NEWKEY(*),ID(-10:*)                                       ACV10766
      LOGICAL FLAG                                                      ACV10767
C                                                                       ACV10768
C Integer data for a DEC system                                         ACV10769
      DATA NF0,N2,N3/1073741823,-2147483648,-1073741824/                ACV10770
C Hexadecimal data for a IBM system                                     ACV10771
C                                                                       ACV10772
C      DATA NF0,N2,N3/Z3FFFFFFF,Z80000000,ZC0000000/                     ACV10773
      NKEY=ID(-4)                                                       ACV10774
      NPR=ID(-2)                                                        ACV10775
      NLC=NPR+1                                                         ACV10776
      NRC=NPR+2                                                         ACV10777
C                                                                       ACV10778
C Special case (empty tree)                                             ACV10779
C                                                                       ACV10780
      IF (ID(0).EQ.-1) THEN                                             ACV10781
         ID(-5)=0                                                       ACV10782
         RETURN                                                         ACV10783
      ENDIF                                                             ACV10784
C                                                                       ACV10785
C Normal case (non-empty tree)                                          ACV10786
C                                                                       ACV10787
      NF=-1                                                             ACV10788
      NA=ID(0)                                                          ACV10789
      NP=ID(0)                                                          ACV10790
      NQ=-1                                                             ACV10791
100   IF (NP.NE.-1) THEN                                                ACV10792
C                                                                       ACV10793
C Check if the balance factor of NP is 0 or not                         ACV10794
C                                                                       ACV10795
         IF (ID(NP+NPR).LT.0) THEN                                      ACV10796
            NA=NP                                                       ACV10797
            NF=NQ                                                       ACV10798
         ENDIF                                                          ACV10799
         DO 101 I=1,NKEY                                                ACV10800
            IF (NEWKEY(I).LT.ID(NP+I)) THEN                             ACV10801
               NQ=NP                                                    ACV10802
               NO=ID(NP+NLC)                                            ACV10803
               IF (NO.EQ.-1) THEN                                       ACV10804
                  NP=-1                                                 ACV10805
               ELSEIF (ID(NO+NRC).EQ.NP) THEN                           ACV10806
                  IF (ISHFT(ID(NP+NPR),-30).EQ.2) THEN                  ACV10807
                     NP=NO                                              ACV10808
                  ELSE                                                  ACV10809
                     NP=-1                                              ACV10810
                  ENDIF                                                 ACV10811
               ELSE                                                     ACV10812
                  NP=NO                                                 ACV10813
               ENDIF                                                    ACV10814
               GOTO 100                                                 ACV10815
            ELSEIF (NEWKEY(I).GT.ID(NP+I)) THEN                         ACV10816
               NQ=NP                                                    ACV10817
               NO=ID(NP+NLC)                                            ACV10818
               IF (NO.EQ.-1) THEN                                       ACV10819
                  NP=-1                                                 ACV10820
               ELSEIF (ID(NO+NRC).EQ.NP) THEN                           ACV10821
                  IF (ISHFT(ID(NP+NPR),-30).EQ.2) THEN                  ACV10822
                     NP=-1                                              ACV10823
                  ELSE                                                  ACV10824
                     NP=NO                                              ACV10825
                  ENDIF                                                 ACV10826
               ELSE                                                     ACV10827
                  NP=ID(NO+NRC)                                         ACV10828
               ENDIF                                                    ACV10829
               GOTO 100                                                 ACV10830
            ENDIF                                                       ACV10831
101      CONTINUE                                                       ACV10832
C                                                                       ACV10833
C If a match is found, increase the node priority, reconstruct the      ACV10834
C linear array, set the pointer to the node location, and RETURN 1      ACV10835
C                                                                       ACV10836
         NA=NP                                                          ACV10837
         FLAG=.TRUE.                                                    ACV10838
         ITEM=IAND(ID(NP+NPR),NF0)                                      ACV10839
         ITEMP=ITEM+ISHFT(ISHFT(ITEM,24),-16)                           ACV10840
C                                                                       ACV10841
C Saturation condition for the priority                                 ACV10842
C                                                                       ACV10843
         IF (ITEMP.LE.NF0) THEN                                         ACV10844
            ID(NP+NPR)=ITEMP+IAND(ID(NP+NPR),N3)                        ACV10845
            ITEMPR=ITEMP                                                ACV10846
         ELSE                                                           ACV10847
            ITEMPR=ITEM                                                 ACV10848
         ENDIF                                                          ACV10849
         NL=ID(-1)                                                      ACV10850
102      NB=NA+NA+NRC                                                   ACV10851
         IF (NB.LT.NL) THEN                                             ACV10852
            NC=NB+NRC                                                   ACV10853
            NBPR=IAND(ID(NB+NPR),NF0)                                   ACV10854
            NCPR=IAND(ID(NC+NPR),NF0)                                   ACV10855
            IF ((NBPR.LE.NCPR).OR.(NC.GE.NL)) THEN                      ACV10856
               IF (NBPR.LT.ITEMPR) THEN                                 ACV10857
                  IF (FLAG) THEN                                        ACV10858
                     NARC=ID(NA+NRC)                                    ACV10859
                     NALC=ID(NA+NLC)                                    ACV10860
                     IF (NARC.NE.-1) THEN                               ACV10861
                        IF (ID(NARC+NLC).EQ.NA) THEN                    ACV10862
                           ID(NARC+NLC)=NL                              ACV10863
                        ELSE                                            ACV10864
                           IF (ID(ID(NARC+NRC)+NLC).EQ.NA) THEN         ACV10865
                              ID(ID(NARC+NRC)+NLC)=NL                   ACV10866
                           ELSE                                         ACV10867
                              ID(ID(NARC+NLC)+NRC)=NL                   ACV10868
                           ENDIF                                        ACV10869
                        ENDIF                                           ACV10870
                     ENDIF                                              ACV10871
                     IF (NALC.NE.-1) THEN                               ACV10872
                        IF (ID(NALC+NRC).EQ.NA) THEN                    ACV10873
                           ID(NALC+NRC)=NL                              ACV10874
                        ELSE                                            ACV10875
                           ID(ID(NALC+NRC)+NRC)=NL                      ACV10876
                        ENDIF                                           ACV10877
                     ENDIF                                              ACV10878
                     DO 105 I=1,NRC                                     ACV10879
105                     ID(NL+I)=ID(NA+I)                               ACV10880
                     FLAG=.FALSE.                                       ACV10881
                  ENDIF                                                 ACV10882
                  NBRC=ID(NB+NRC)                                       ACV10883
                  NBLC=ID(NB+NLC)                                       ACV10884
                  IF (NBRC.EQ.-1) THEN                                  ACV10885
                     ID(0)=NA                                           ACV10886
                  ELSE                                                  ACV10887
                     IF (ID(NBRC+NLC).EQ.NB) THEN                       ACV10888
                        ID(NBRC+NLC)=NA                                 ACV10889
                     ELSE                                               ACV10890
                        IF (ID(ID(NBRC+NRC)+NLC).EQ.NB) THEN            ACV10891
                           ID(ID(NBRC+NRC)+NLC)=NA                      ACV10892
                        ELSE                                            ACV10893
                           ID(ID(NBRC+NLC)+NRC)=NA                      ACV10894
                        ENDIF                                           ACV10895
                     ENDIF                                              ACV10896
                  ENDIF                                                 ACV10897
                  IF (NBLC.NE.-1) THEN                                  ACV10898
                     IF (ID(NBLC+NRC).EQ.NB) THEN                       ACV10899
                        ID(NBLC+NRC)=NA                                 ACV10900
                     ELSE                                               ACV10901
                        ID(ID(NBLC+NRC)+NRC)=NA                         ACV10902
                     ENDIF                                              ACV10903
                  ENDIF                                                 ACV10904
                  DO 109 I=1,NRC                                        ACV10905
                     ID(NA+I)=ID(NB+I)                                  ACV10906
109               CONTINUE                                              ACV10907
                  NA=NB                                                 ACV10908
                  GOTO 102                                              ACV10909
               ENDIF                                                    ACV10910
            ELSE                                                        ACV10911
               IF (NCPR.LT.ITEMPR) THEN                                 ACV10912
                  IF (FLAG) THEN                                        ACV10913
                     NARC=ID(NA+NRC)                                    ACV10914
                     NALC=ID(NA+NLC)                                    ACV10915
                     IF (NARC.NE.-1) THEN                               ACV10916
                        IF (ID(NARC+NLC).EQ.NA) THEN                    ACV10917
                           ID(NARC+NLC)=NL                              ACV10918
                        ELSE                                            ACV10919
                           IF (ID(ID(NARC+NRC)+NLC).EQ.NA) THEN         ACV10920
                              ID(ID(NARC+NRC)+NLC)=NL                   ACV10921
                           ELSE                                         ACV10922
                              ID(ID(NARC+NLC)+NRC)=NL                   ACV10923
                           ENDIF                                        ACV10924
                        ENDIF                                           ACV10925
                     ENDIF                                              ACV10926
                     IF (NALC.NE.-1) THEN                               ACV10927
                        IF (ID(NALC+NRC).EQ.NA) THEN                    ACV10928
                           ID(NALC+NRC)=NL                              ACV10929
                        ELSE                                            ACV10930
                           ID(ID(NALC+NRC)+NRC)=NL                      ACV10931
                        ENDIF                                           ACV10932
                     ENDIF                                              ACV10933
                     DO 112 I=1,NRC                                     ACV10934
112                     ID(NL+I)=ID(NA+I)                               ACV10935
                     FLAG=.FALSE.                                       ACV10936
                  ENDIF                                                 ACV10937
                  NCRC=ID(NC+NRC)                                       ACV10938
                  NCLC=ID(NC+NLC)                                       ACV10939
                  IF (NCRC.EQ.-1) THEN                                  ACV10940
                     ID(0)=NA                                           ACV10941
                  ELSE                                                  ACV10942
                     IF (ID(NCRC+NLC).EQ.NC) THEN                       ACV10943
                        ID(NCRC+NLC)=NA                                 ACV10944
                     ELSE                                               ACV10945
                        IF (ID(ID(NCRC+NRC)+NLC).EQ.NC) THEN            ACV10946
                           ID(ID(NCRC+NRC)+NLC)=NA                      ACV10947
                        ELSE                                            ACV10948
                           ID(ID(NCRC+NLC)+NRC)=NA                      ACV10949
                        ENDIF                                           ACV10950
                     ENDIF                                              ACV10951
                  ENDIF                                                 ACV10952
                  IF (NCLC.NE.-1) THEN                                  ACV10953
                     IF (ID(NCLC+NRC).EQ.NC) THEN                       ACV10954
                        ID(NCLC+NRC)=NA                                 ACV10955
                     ELSE                                               ACV10956
                        ID(ID(NCLC+NRC)+NRC)=NA                         ACV10957
                     ENDIF                                              ACV10958
                  ENDIF                                                 ACV10959
                  DO 116 I=1,NRC                                        ACV10960
                     ID(NA+I)=ID(NC+I)                                  ACV10961
116               CONTINUE                                              ACV10962
                  NA=NC                                                 ACV10963
                  GOTO 102                                              ACV10964
               ENDIF                                                    ACV10965
            ENDIF                                                       ACV10966
         ENDIF                                                          ACV10967
118      IF (FLAG) GOTO 124                                             ACV10968
         NLRC=ID(NL+NRC)                                                ACV10969
         NLLC=ID(NL+NLC)                                                ACV10970
         IF (NLRC.EQ.-1) THEN                                           ACV10971
            ID(0)=NA                                                    ACV10972
         ELSE                                                           ACV10973
            IF (ID(NLRC+NLC).EQ.NL) THEN                                ACV10974
               ID(NLRC+NLC)=NA                                          ACV10975
            ELSE                                                        ACV10976
               IF (ID(ID(NLRC+NRC)+NLC).EQ.NL) THEN                     ACV10977
                  ID(ID(NLRC+NRC)+NLC)=NA                               ACV10978
               ELSE                                                     ACV10979
                  ID(ID(NLRC+NLC)+NRC)=NA                               ACV10980
               ENDIF                                                    ACV10981
            ENDIF                                                       ACV10982
         ENDIF                                                          ACV10983
         IF (NLLC.NE.-1) THEN                                           ACV10984
            IF (ID(NLLC+NRC).EQ.NL) THEN                                ACV10985
               ID(NLLC+NRC)=NA                                          ACV10986
            ELSE                                                        ACV10987
               ID(ID(NLLC+NRC)+NRC)=NA                                  ACV10988
            ENDIF                                                       ACV10989
         ENDIF                                                          ACV10990
         DO 120 I=1,NRC                                                 ACV10991
            ID(NA+I)=ID(NL+I)                                           ACV10992
            ID(NL+I)=0                                                  ACV10993
120      CONTINUE                                                       ACV10994
124      ID(-5)=NA                                                      ACV10995
         RETURN 1                                                       ACV10996
      ENDIF                                                             ACV10997
      ID(-8)=NF                                                         ACV10998
      ID(-7)=NA                                                         ACV10999
      ID(-6)=NQ                                                         ACV11000
      ID(-5)=ID(-1)+NRC                                                 ACV11001
      RETURN                                                            ACV11002
      END                                                               ACV11003
C ----------------------------------------------------------------------ACV11004
C                                                                       ACV11005
C                            ************                               ACV11006
C                            *** TADD ***                               ACV11007
C                            ************                               ACV11008
C                                                                       ACV11009
C The subroutine TADD allows the user to add elements to an existing    ACV11010
C database provided space is available. If it is not and the priority   ACV11011
C of the incoming element is greater than the lowest priority element   ACV11012
C currently in the database that element is eliminated to make room for ACV11013
C the new one. The last scenario is repeated as necessary to accommodateACV11014
C the new element. The logic of the program is the following:           ACV11015
C                                                                       ACV11016
C Step 1) check if either the tree or buffer is full                    ACV11017
C         if yes, check if the priority of incoming item is higher      ACV11018
C                 than that of the lowest priority node in the tree     ACV11019
C                 if yes, delete lowest priority node and goto Step 1)  ACV11020
C                 else RETURN                                           ACV11021
C         else    goto Step 2)                                          ACV11022
C Step 2) insert the incoming item into the database                    ACV11023
C                                                                       ACV11024
C ----------------------------------------------------------------------ACV11025
C                                                                       ACV11026
      SUBROUTINE TADD(NEWKEY,NEWDAT,BULOAD,NOSIZE,NPBASE,ID,BUFFER,     ACV11027
     &                LLBUFF)                                           ACV11028
C                                                                       ACV11029
      INTEGER NEWKEY(*),NEWDAT(*),ID(-10:*),LLBUFF(-2:*)                ACV11030
      REAL    BULOAD(*),BUFFER(*)                     ! *** storage typeACV11031
C     REAL*8  BULOAD(*),BUFFER(*)                     ! *** storage typeACV11032
      DATA    NF0/1073741823/                                            ACV11033
C                                                                       ACV11034
C ...define the number of keys in a node                                ACV11035
C                                                                       ACV11036
      NKEY=ID(-4)                                                       ACV11037
C                                                                       ACV11038
C ...calculate current priority                                         ACV11039
C                                                                       ACV11040
      NPCURR=NPBASE+ISHFT(NPBASE,8)                                     ACV11041
C                                                                       ACV11042
C ...check to see if the tree and buffer can receive the new item       ACV11043
C                                                                       ACV11044
      IF ((ID(-10).LE.ID(-9)).AND.(LLBUFF(0).GE.NOSIZE)) THEN           ACV11045
         GOTO 300                                                       ACV11046
C                                                                       ACV11047
C ...check if priority of new item is greater than lowest priority      ACV11048
C                                                                       ACV11049
      ELSE                                                              ACV11050
100      IF (NPCURR.GT.IAND(ID(ID(-2)),NF0)) THEN                       ACV11051
C                                                                       ACV11052
C Rearrange LLBUFF array which holds free space list information        ACV11053
C                                                                       ACV11054
            IX=ID(NKEY+1)                                               ACV11055
            IF (LLBUFF(0).EQ.0) THEN                                    ACV11056
               LLBUFF(0)=ID(NKEY+2)                                     ACV11057
               LLBUFF(-1)=IX                                            ACV11058
            ELSE                                                        ACV11059
               LLBUFF(0)=LLBUFF(0)+ID(NKEY+2)                           ACV11060
               LLBUFF(LLBUFF(-2))=IX                                    ACV11061
            ENDIF                                                       ACV11062
            LLBUFF(-2)=IX                                               ACV11063
            DO 200 I=1,ID(NKEY+2)-1                                     ACV11064
               IX=LLBUFF(IX)                                            ACV11065
200            LLBUFF(-2)=IX                                            ACV11066
C                                                                       ACV11067
            CALL TDEL(ID)                                               ACV11068
C                                                                       ACV11069
            IF (LLBUFF(0).LT.NOSIZE) GOTO 100                           ACV11070
            CALL TCHK(NEWKEY,ID)                                        ACV11071
            GOTO 300                                                    ACV11072
C                                                                       ACV11073
C ...doing nothing when new priority is less than lowest priority       ACV11074
C                                                                       ACV11075
         ELSE                                                           ACV11076
            RETURN                                                      ACV11077
         ENDIF                                                          ACV11078
      ENDIF                                                             ACV11079
300   CONTINUE                                                          ACV11080
C                                                                       ACV11081
C Load the key, data and its priority into the tree                     ACV11082
C                                                                       ACV11083
C ...last key position of the current node                              ACV11084
C                                                                       ACV11085
      IFIND=ID(-5)+NKEY                                                 ACV11086
C                                                                       ACV11087
C ...load integer data into tree after the last key                     ACV11088
C                                                                       ACV11089
      NEWDAT(1)=LLBUFF(-1)                                              ACV11090
      NEWDAT(2)=NOSIZE                                                  ACV11091
      NDAT=ID(-3)-ID(-4)                                                ACV11092
      DO 400 IKK=1,NDAT                                                 ACV11093
400   ID(IFIND+IKK)=NEWDAT(IKK)                                         ACV11094
C                                                                       ACV11095
C ...load real data into the BUFFER and update LLBUFF                   ACV11096
C                                                                       ACV11097
      DO 500 IKK=1,NOSIZE                                               ACV11098
         BUFFER(LLBUFF(-1))=BULOAD(IKK)                                 ACV11099
500      LLBUFF(-1)=LLBUFF(LLBUFF(-1))                                  ACV11100
C                                                                       ACV11101
C ...reduce size of free space                                          ACV11102
C                                                                       ACV11103
      LLBUFF(0)=LLBUFF(0)-NOSIZE                                        ACV11104
C                                                                       ACV11105
C ...load the priority                                                  ACV11106
C                                                                       ACV11107
      IFIND=ID(-5)+ID(-2)                                               ACV11108
      ID(IFIND)=NPCURR                                                  ACV11109
C                                                                       ACV11110
C ...call TINS to complete (balance tree and restructure heap)          ACV11111
C                                                                       ACV11112
      CALL TINS(NEWKEY,ID)                                              ACV11113
      RETURN                                                            ACV11114
      END                                                               ACV11115
C ----------------------------------------------------------------------ACV11116
C                                                                       ACV11117
C                            ************                               ACV11118
C                            *** TINS ***                               ACV11119
C                            ************                               ACV11120
C                                                                       ACV11121
C The subroutine TINS is used to insert a new item into a WST. A call   ACV11122
C to TINS should follow a call to TCHK since TCHK tell where to insert  ACV11123
C the new item. After the item is placed in the WST, TINS rebalances theACV11124
C the tree and reconstruct the priority heap.                           ACV11125
C                                                                       ACV11126
C ----------------------------------------------------------------------ACV11127
C                                                                       ACV11128
      SUBROUTINE TINS(NEWKEY,ID)                                        ACV11129
C                                                                       ACV11130
      INTEGER NEWKEY(*),ID(-10:*)                                       ACV11131
C                                                                       ACV11132
C Integer data for a DEC system                                         ACV11133
      DATA NF0,N2,N3/1073741823,-2147483648,-1073741824/                ACV11134
C Hexadecimal data for a IBM  system                                    ACV11135
C                                                                       ACV11136
C      DATA NF0,N2,N3/Z3FFFFFFF,Z80000000,ZC0000000/                     ACV11137
      ID(-10)=ID(-10)+1                                                 ACV11138
      NKEY=ID(-4)                                                       ACV11139
      NPR=ID(-2)                                                        ACV11140
      NLC=NPR+1                                                         ACV11141
      NRC=NPR+2                                                         ACV11142
C                                                                       ACV11143
C Special case (empty tree)                                             ACV11144
C                                                                       ACV11145
      IF (ID(0).EQ.-1) THEN                                             ACV11146
         ID(-1)=NRC                                                     ACV11147
         ID(0)=0                                                        ACV11148
C                                                                       ACV11149
C The first node position                                               ACV11150
C                                                                       ACV11151
         ID(NRC)=-1                                                     ACV11152
         ID(NLC)=-1                                                     ACV11153
         DO 10 I=1,NKEY                                                 ACV11154
10          ID(I)=NEWKEY(I)                                             ACV11155
         RETURN                                                         ACV11156
      ENDIF                                                             ACV11157
C                                                                       ACV11158
C Normal case (non-empty tree)                                          ACV11159
C                                                                       ACV11160
      NY=ID(-5)                                                         ACV11161
      NQ=ID(-6)                                                         ACV11162
      NA=ID(-7)                                                         ACV11163
      NF=ID(-8)                                                         ACV11164
C                                                                       ACV11165
C Set the pointer to indicate the position of the new node and load key ACV11166
C                                                                       ACV11167
      DO 130 I=1,NKEY                                                   ACV11168
130      ID(NY+I)=NEWKEY(I)                                             ACV11169
      ID(NY+NLC)=-1                                                     ACV11170
      IF (ID(NQ+NLC).EQ.-1) THEN                                        ACV11171
         ID(NY+NRC)=NQ                                                  ACV11172
         ID(NQ+NLC)=NY                                                  ACV11173
      ELSE                                                              ACV11174
         DO 131 I=1,NKEY                                                ACV11175
            IF (NEWKEY(I).LT.ID(NQ+I)) THEN                             ACV11176
               ID(NY+NRC)=ID(NQ+NLC)                                    ACV11177
               ID(NQ+NLC)=NY                                            ACV11178
               ID(NQ+NPR)=IAND(ID(NQ+NPR),NF0)                          ACV11179
               GOTO 999                                                 ACV11180
            ELSEIF (NEWKEY(I).GT.ID(NQ+I)) THEN                         ACV11181
               ID(NY+NRC)=NQ                                            ACV11182
               ID(ID(NQ+NLC)+NRC)=NY                                    ACV11183
               ID(NQ+NPR)=IAND(ID(NQ+NPR),NF0)                          ACV11184
               GOTO 999                                                 ACV11185
            ENDIF                                                       ACV11186
131      CONTINUE                                                       ACV11187
      ENDIF                                                             ACV11188
      DO 161 I=1,NKEY                                                   ACV11189
         IF (NEWKEY(I).LT.ID(NA+I)) THEN                                ACV11190
            NP=ID(NA+NLC)                                               ACV11191
            NB=NP                                                       ACV11192
            IND=N2                                                      ACV11193
            GOTO 200                                                    ACV11194
         ELSEIF (NEWKEY(I).GT.ID(NA+I)) THEN                            ACV11195
            NO=ID(NA+NLC)                                               ACV11196
            IF (ID(NO+NRC).EQ.NA) THEN                                  ACV11197
               NP=NO                                                    ACV11198
            ELSE                                                        ACV11199
               NP=ID(NO+NRC)                                            ACV11200
            ENDIF                                                       ACV11201
            NB=NP                                                       ACV11202
            IND=N3                                                      ACV11203
            GOTO 200                                                    ACV11204
         ENDIF                                                          ACV11205
161   CONTINUE                                                          ACV11206
200   IF (NP.NE.NY) THEN                                                ACV11207
         DO 211 I=1,NKEY                                                ACV11208
            IF (NEWKEY(I).LT.ID(NP+I)) THEN                             ACV11209
               ID(NP+NPR)=IAND(ID(NP+NPR),NF0)+N2                       ACV11210
               NP=ID(NP+NLC)                                            ACV11211
               GOTO 200                                                 ACV11212
            ELSEIF (NEWKEY(I).GT.ID(NP+I)) THEN                         ACV11213
               ID(NP+NPR)=IAND(ID(NP+NPR),NF0)+N3                       ACV11214
               NO=ID(NP+NLC)                                            ACV11215
               IF (ID(NO+NRC).EQ.NP) THEN                               ACV11216
                  NP=NO                                                 ACV11217
               ELSE                                                     ACV11218
                  NP=ID(NO+NRC)                                         ACV11219
               ENDIF                                                    ACV11220
               GOTO 200                                                 ACV11221
            ENDIF                                                       ACV11222
211      CONTINUE                                                       ACV11223
      ENDIF                                                             ACV11224
      IF (ID(NA+NPR).GE.0) THEN                                         ACV11225
         ID(NA+NPR)=IAND(ID(NA+NPR),NF0)+IND                            ACV11226
         GOTO 999                                                       ACV11227
      ELSEIF (IAND(ID(NA+NPR),N3).NE.IND) THEN                          ACV11228
         ID(NA+NPR)=IAND(ID(NA+NPR),NF0)                                ACV11229
         GOTO 999                                                       ACV11230
      ENDIF                                                             ACV11231
      IF (IND.EQ.N2) THEN                                               ACV11232
         IF (ISHFT(ID(NB+NPR),-30).EQ.2) THEN                           ACV11233
            IF (ID(ID(NA+NLC)+NRC).EQ.NA) THEN                          ACV11234
               ID(NB+NRC)=ID(NA+NRC)                                    ACV11235
               ID(ID(NB+NLC)+NRC)=NA                                    ACV11236
               ID(NA+NLC)=-1                                            ACV11237
               ID(NA+NRC)=NB                                            ACV11238
            ELSE                                                        ACV11239
               NC=ID(ID(NB+NLC)+NRC)                                    ACV11240
               ID(NC+NRC)=ID(NB+NRC)                                    ACV11241
               ID(NB+NRC)=ID(NA+NRC)                                    ACV11242
               ID(ID(NB+NLC)+NRC)=NA                                    ACV11243
               ID(NA+NLC)=NC                                            ACV11244
               ID(NA+NRC)=NB                                            ACV11245
            ENDIF                                                       ACV11246
            ID(NA+NPR)=IAND(ID(NA+NPR),NF0)                             ACV11247
            ID(NB+NPR)=IAND(ID(NB+NPR),NF0)                             ACV11248
         ELSE                                                           ACV11249
            IF (ID(ID(NB+NLC)+NRC).EQ.NB) THEN                          ACV11250
               NC=ID(NB+NLC)                                            ACV11251
            ELSE                                                        ACV11252
               NC=ID(ID(NB+NLC)+NRC)                                    ACV11253
            ENDIF                                                       ACV11254
            IF (ID(NC+NPR).GE.0) THEN                                   ACV11255
               ID(NC+NRC)=ID(NA+NRC)                                    ACV11256
               ID(NA+NLC)=-1                                            ACV11257
               ID(NA+NRC)=NC                                            ACV11258
               ID(NB+NLC)=-1                                            ACV11259
               ID(NC+NLC)=NB                                            ACV11260
               ID(NA+NPR)=IAND(ID(NA+NPR),NF0)                          ACV11261
               ID(NB+NPR)=IAND(ID(NB+NPR),NF0)                          ACV11262
            ELSE                                                        ACV11263
               IF (ID(ID(NC+NLC)+NRC).EQ.NC) THEN                       ACV11264
                  IF (ISHFT(ID(NC+NPR),-30).EQ.2) THEN                  ACV11265
                     ID(NA+NLC)=ID(NB+NRC)                              ACV11266
                     ID(ID(NB+NLC)+NRC)=ID(NC+NLC)                      ACV11267
                     ID(ID(NC+NLC)+NRC)=NB                              ACV11268
                  ELSE                                                  ACV11269
                     ID(NA+NLC)=ID(NC+NLC)                              ACV11270
                     ID(ID(NA+NLC)+NRC)=ID(NB+NRC)                      ACV11271
                     ID(ID(NB+NLC)+NRC)=NB                              ACV11272
                  ENDIF                                                 ACV11273
               ELSE                                                     ACV11274
                  ID(NA+NLC)=ID(ID(NC+NLC)+NRC)                         ACV11275
                  ID(ID(NA+NLC)+NRC)=ID(NB+NRC)                         ACV11276
                  ID(ID(NB+NLC)+NRC)=ID(NC+NLC)                         ACV11277
                  ID(ID(NC+NLC)+NRC)=NB                                 ACV11278
               ENDIF                                                    ACV11279
               ID(NC+NLC)=NB                                            ACV11280
               ID(NC+NRC)=ID(NA+NRC)                                    ACV11281
               ID(NA+NRC)=NC                                            ACV11282
               ID(NB+NRC)=NA                                            ACV11283
               IF (ISHFT(ID(NC+NPR),-30).EQ.2) THEN                     ACV11284
                  ID(NA+NPR)=IAND(ID(NA+NPR),NF0)+N3                    ACV11285
                  ID(NB+NPR)=IAND(ID(NB+NPR),NF0)                       ACV11286
               ELSE                                                     ACV11287
                  ID(NA+NPR)=IAND(ID(NA+NPR),NF0)                       ACV11288
                  ID(NB+NPR)=IAND(ID(NB+NPR),NF0)+N2                    ACV11289
               ENDIF                                                    ACV11290
               ID(NC+NPR)=IAND(ID(NC+NPR),NF0)                          ACV11291
            ENDIF                                                       ACV11292
            NB=NC                                                       ACV11293
         ENDIF                                                          ACV11294
      ELSE                                                              ACV11295
         IF (ISHFT(ID(NB+NPR),-30).EQ.3) THEN                           ACV11296
            IF (ID(ID(NA+NLC)+NRC).EQ.NA) THEN                          ACV11297
               ID(NA+NLC)=-1                                            ACV11298
               ID(NB+NRC)=ID(NA+NRC)                                    ACV11299
               ID(NA+NRC)=ID(NB+NLC)                                    ACV11300
               ID(NB+NLC)=NA                                            ACV11301
            ELSE                                                        ACV11302
               NC=ID(NB+NLC)                                            ACV11303
               ID(NB+NLC)=NA                                            ACV11304
               ID(NB+NRC)=ID(NA+NRC)                                    ACV11305
               ID(ID(NA+NLC)+NRC)=NC                                    ACV11306
               ID(NA+NRC)=ID(NC+NRC)                                    ACV11307
               ID(NC+NRC)=NA                                            ACV11308
            ENDIF                                                       ACV11309
            ID(NA+NPR)=IAND(ID(NA+NPR),NF0)                             ACV11310
            ID(NB+NPR)=IAND(ID(NB+NPR),NF0)                             ACV11311
         ELSE                                                           ACV11312
            NC=ID(NB+NLC)                                               ACV11313
            IF (ID(NC+NPR).GE.0) THEN                                   ACV11314
               ID(NC+NRC)=ID(NA+NRC)                                    ACV11315
               ID(NC+NLC)=NA                                            ACV11316
               ID(NA+NLC)=-1                                            ACV11317
               ID(NA+NRC)=NB                                            ACV11318
               ID(NB+NLC)=-1                                            ACV11319
               ID(NB+NRC)=NC                                            ACV11320
               ID(NA+NPR)=IAND(ID(NA+NPR),NF0)                          ACV11321
               ID(NB+NPR)=IAND(ID(NB+NPR),NF0)                          ACV11322
            ELSE                                                        ACV11323
               IF (ID(ID(NC+NLC)+NRC).EQ.NC) THEN                       ACV11324
                  IF (ISHFT(ID(NC+NPR),-30).EQ.3) THEN                  ACV11325
                     ID(ID(NA+NLC)+NRC)=NA                              ACV11326
                     ID(ID(NC+NLC)+NRC)=ID(NC+NRC)                      ACV11327
                     ID(NB+NLC)=ID(NC+NLC)                              ACV11328
                  ELSE                                                  ACV11329
                     ID(ID(NA+NLC)+NRC)=ID(NC+NLC)                      ACV11330
                     ID(ID(NC+NLC)+NRC)=NA                              ACV11331
                     ID(NB+NLC)=ID(NC+NRC)                              ACV11332
                  ENDIF                                                 ACV11333
               ELSE                                                     ACV11334
                  ID(ID(NA+NLC)+NRC)=ID(NC+NLC)                         ACV11335
                  ID(NB+NLC)=ID(ID(NC+NLC)+NRC)                         ACV11336
                  ID(ID(NB+NLC)+NRC)=ID(NC+NRC)                         ACV11337
                  ID(ID(NC+NLC)+NRC)=NA                                 ACV11338
               ENDIF                                                    ACV11339
               ID(NC+NLC)=NA                                            ACV11340
               ID(NC+NRC)=ID(NA+NRC)                                    ACV11341
               ID(NA+NRC)=NB                                            ACV11342
               ID(NB+NRC)=NC                                            ACV11343
               IF (ISHFT(ID(NC+NPR),-30).EQ.3) THEN                     ACV11344
                  ID(NA+NPR)=IAND(ID(NA+NPR),NF0)+N2                    ACV11345
                  ID(NB+NPR)=IAND(ID(NB+NPR),NF0)                       ACV11346
               ELSE                                                     ACV11347
                  ID(NA+NPR)=IAND(ID(NA+NPR),NF0)                       ACV11348
                  ID(NB+NPR)=IAND(ID(NB+NPR),NF0)+N3                    ACV11349
               ENDIF                                                    ACV11350
               ID(NC+NPR)=IAND(ID(NC+NPR),NF0)                          ACV11351
            ENDIF                                                       ACV11352
            NB=NC                                                       ACV11353
         ENDIF                                                          ACV11354
      ENDIF                                                             ACV11355
      IF (NF.EQ.-1) THEN                                                ACV11356
         ID(0)=NB                                                       ACV11357
      ELSEIF (NA.EQ.ID(NF+NLC)) THEN                                    ACV11358
         ID(NF+NLC)=NB                                                  ACV11359
      ELSE                                                              ACV11360
         ID(ID(NF+NLC)+NRC)=NB                                          ACV11361
      ENDIF                                                             ACV11362
999   CONTINUE                                                          ACV11363
C                                                                       ACV11364
C Reconstruct the priority array                                        ACV11365
C                                                                       ACV11366
      NL=ID(-1)                                                         ACV11367
      NRC2=NRC+NRC                                                      ACV11368
1000  IF (MOD(NL,NRC2).EQ.0) THEN                                       ACV11369
         NA=NL/2-NRC                                                    ACV11370
      ELSE                                                              ACV11371
         NA=(NL-NRC)/2                                                  ACV11372
      ENDIF                                                             ACV11373
      IF (IAND(ID(NA+NPR),NF0).GT.IAND(ID(NY+NPR),NF0)) THEN            ACV11374
         NARC=ID(NA+NRC)                                                ACV11375
         NALC=ID(NA+NLC)                                                ACV11376
         IF (NARC.EQ.-1) THEN                                           ACV11377
            ID(0)=NL                                                    ACV11378
         ELSE                                                           ACV11379
            IF (ID(NARC+NLC).EQ.NA) THEN                                ACV11380
               ID(NARC+NLC)=NL                                          ACV11381
            ELSE                                                        ACV11382
               IF (ID(ID(NARC+NRC)+NLC).EQ.NA) THEN                     ACV11383
                  ID(ID(NARC+NRC)+NLC)=NL                               ACV11384
               ELSE                                                     ACV11385
                  ID(ID(NARC+NLC)+NRC)=NL                               ACV11386
               ENDIF                                                    ACV11387
            ENDIF                                                       ACV11388
         ENDIF                                                          ACV11389
         IF (NALC.NE.-1) THEN                                           ACV11390
            IF (ID(NALC+NRC).EQ.NA) THEN                                ACV11391
               ID(NALC+NRC)=NL                                          ACV11392
            ELSE                                                        ACV11393
               ID(ID(NALC+NRC)+NRC)=NL                                  ACV11394
            ENDIF                                                       ACV11395
         ENDIF                                                          ACV11396
         DO 1400 I=1,NRC                                                ACV11397
            ID(NL+I)=ID(NA+I)                                           ACV11398
1400     CONTINUE                                                       ACV11399
         NL=NA                                                          ACV11400
         IF (NL.GT.0) GOTO 1000                                         ACV11401
      ENDIF                                                             ACV11402
      NYRC=ID(NY+NRC)                                                   ACV11403
      NYLC=ID(NY+NLC)                                                   ACV11404
      IF (NYRC.EQ.-1) THEN                                              ACV11405
         ID(0)=NL                                                       ACV11406
      ELSE                                                              ACV11407
         IF (ID(NYRC+NLC).EQ.NY) THEN                                   ACV11408
            ID(NYRC+NLC)=NL                                             ACV11409
         ELSE                                                           ACV11410
            IF (ID(ID(NYRC+NRC)+NLC).EQ.NY) THEN                        ACV11411
               ID(ID(NYRC+NRC)+NLC)=NL                                  ACV11412
            ELSE                                                        ACV11413
               ID(ID(NYRC+NLC)+NRC)=NL                                  ACV11414
            ENDIF                                                       ACV11415
         ENDIF                                                          ACV11416
      ENDIF                                                             ACV11417
      IF (NYLC.NE.-1) THEN                                              ACV11418
         IF (ID(NYLC+NRC).EQ.NY) THEN                                   ACV11419
            ID(NYLC+NRC)=NL                                             ACV11420
         ELSE                                                           ACV11421
            ID(ID(NYLC+NRC)+NRC)=NL                                     ACV11422
         ENDIF                                                          ACV11423
      ENDIF                                                             ACV11424
      DO 7000 I=1,NRC                                                   ACV11425
         ID(NL+I)=ID(NY+I)                                              ACV11426
7000  CONTINUE                                                          ACV11427
      ID(-1)=ID(-1)+NRC                                                 ACV11428
      ID(-5)=NL                                                         ACV11429
      RETURN                                                            ACV11430
      END                                                               ACV11431
C ----------------------------------------------------------------------ACV11432
C                                                                       ACV11433
C                            ************                               ACV11434
C                            *** TDEL ***                               ACV11435
C                            ************                               ACV11436
C                                                                       ACV11437
C The subroutine TDEL deletes a node from a weighted search tree,       ACV11438
C keeping the binary tree balanced. The node deleted from the tree      ACV11439
C is the first element of the linear array ID. Since the linear array   ACV11440
C also has a heap structure, the root of a heap, which is the first     ACV11441
C element of an array, is the lowest priority element (minheap). The    ACV11442
C subroutine, TDEL, also manages the free space in the array BUFFER.    ACV11443
C The pointer array LLBUFF is used for this purpose.                    ACV11444
C                                                                       ACV11445
C ----------------------------------------------------------------------ACV11446
C                                                                       ACV11447
      SUBROUTINE TDEL(ID)                                               ACV11448
C                                                                       ACV11449
      INTEGER ID(-10:*)                                                 ACV11450
      LOGICAL FLAG                                                      ACV11451
C                                                                       ACV11452
C Integer data for a DEC system                                         ACV11453
      DATA NF0,N2,N3/1073741823,-2147483648,-1073741824/                ACV11454
C Hexadecimal data for a IBM system                                     ACV11455
C                                                                       ACV11456
C      DATA NF0,N2,N3/Z3FFFFFFF,Z80000000,ZC0000000/                     ACV11457
C                                                                       ACV11458
C Initialize some integer constants                                     ACV11459
C                                                                       ACV11460
      NKEY=ID(-4)                                                       ACV11461
      NPR=ID(-2)                                                        ACV11462
      NLC=NPR+1                                                         ACV11463
      NRC=NPR+2                                                         ACV11464
C                                                                       ACV11465
C NA keeps track of most recent node with BF(0)                         ACV11466
C NF is the parent of NA                                                ACV11467
C NQ follows NP through the tree                                        ACV11468
C                                                                       ACV11469
      NF=-1                                                             ACV11470
      NA=ID(0)                                                          ACV11471
      NP=ID(0)                                                          ACV11472
      NQ=-1                                                             ACV11473
      NR=-1                                                             ACV11474
100   IF (NP.NE.-1) THEN                                                ACV11475
C                                                                       ACV11476
C Looking for the last node with BF(0) that is not a leaf               ACV11477
C                                                                       ACV11478
         IF ((ID(NP+NPR).GE.0).AND.(ID(NP+NLC).NE.-1)) THEN             ACV11479
            NA=NP                                                       ACV11480
            NF=NQ                                                       ACV11481
         ELSEIF (NQ.NE.-1) THEN                                         ACV11482
            IF (FLAG) THEN                                              ACV11483
               NO=ID(NQ+NLC)                                            ACV11484
               IF (ID(NO+NRC).EQ.NQ) THEN                               ACV11485
                  NRCHILD=NO                                            ACV11486
               ELSE                                                     ACV11487
                  NRCHILD=ID(NO+NRC)                                    ACV11488
               ENDIF                                                    ACV11489
               IF ((ISHFT(ID(NQ+NPR),-30).EQ.3)                         ACV11490
     *             .AND.(ID(NRCHILD+NPR).GE.0)) THEN                    ACV11491
                  NA=NQ                                                 ACV11492
                  NF=NR                                                 ACV11493
               ENDIF                                                    ACV11494
            ELSE                                                        ACV11495
               IF ((ISHFT(ID(NQ+NPR),-30).EQ.2)                         ACV11496
     *             .AND.(ID(ID(NQ+NLC)+NPR).GE.0)) THEN                 ACV11497
                  NA=NQ                                                 ACV11498
                  NF=NR                                                 ACV11499
               ENDIF                                                    ACV11500
            ENDIF                                                       ACV11501
         ENDIF                                                          ACV11502
         DO 101 I=1,NKEY                                                ACV11503
            IF (ID(I).LT.ID(NP+I)) THEN                                 ACV11504
               NR=NQ                                                    ACV11505
               NQ=NP                                                    ACV11506
               NP=ID(NP+NLC)                                            ACV11507
               FLAG=.TRUE.                                              ACV11508
               GOTO 100                                                 ACV11509
            ELSEIF (ID(I).GT.ID(NP+I)) THEN                             ACV11510
               NR=NQ                                                    ACV11511
               NQ=NP                                                    ACV11512
               NO=ID(NP+NLC)                                            ACV11513
               IF (ID(NO+NRC).EQ.NP) THEN                               ACV11514
                  NP=NO                                                 ACV11515
               ELSE                                                     ACV11516
                  NP=ID(NO+NRC)                                         ACV11517
               ENDIF                                                    ACV11518
               FLAG=.FALSE.                                             ACV11519
               GOTO 100                                                 ACV11520
            ENDIF                                                       ACV11521
101      CONTINUE                                                       ACV11522
C                                                                       ACV11523
C A match is found                                                      ACV11524
C                                                                       ACV11525
         GOTO 130                                                       ACV11526
      ENDIF                                                             ACV11527
C                                                                       ACV11528
C A match is not found                                                  ACV11529
C                                                                       ACV11530
      WRITE(6,*) 'TREE IS EMPTY.'                                       ACV11531
      GOTO 9999                                                         ACV11532
130   CONTINUE                                                          ACV11533
      ID(-10)=ID(-10)-1                                                 ACV11534
C                                                                       ACV11535
C Matched node does not have a child                                    ACV11536
C                                                                       ACV11537
      IF (ID(NP+NLC).EQ.-1) THEN                                        ACV11538
         IF (NQ.EQ.-1) THEN                                             ACV11539
            ID(0)=-1                                                    ACV11540
            ID(-1)=0                                                    ACV11541
            GOTO 9999                                                   ACV11542
         ELSEIF (ID(NP+NRC).EQ.NQ) THEN                                 ACV11543
            IF (ID(NQ+NLC).EQ.NP) THEN                                  ACV11544
               ID(NQ+NLC)=-1                                            ACV11545
            ELSE                                                        ACV11546
               ID(ID(NQ+NLC)+NRC)=NQ                                    ACV11547
            ENDIF                                                       ACV11548
         ELSE                                                           ACV11549
            ID(NQ+NLC)=ID(NP+NRC)                                       ACV11550
         ENDIF                                                          ACV11551
C                                                                       ACV11552
C Matched node has only one child                                       ACV11553
C                                                                       ACV11554
      ELSEIF (ID(ID(NP+NLC)+NRC).EQ.NP) THEN                            ACV11555
         NO=ID(NP+NLC)                                                  ACV11556
         IF (ID(NP+NRC).EQ.-1) THEN                                     ACV11557
            ID(0)=ID(NO)                                                ACV11558
            ID(NO+NRC)=-1                                               ACV11559
            GOTO 999                                                    ACV11560
         ELSE                                                           ACV11561
            IF (ID(NQ+NLC).EQ.NP) THEN                                  ACV11562
               ID(NQ+NLC)=NO                                            ACV11563
            ELSE                                                        ACV11564
               ID(ID(NQ+NLC)+NRC)=NO                                    ACV11565
            ENDIF                                                       ACV11566
            ID(NO+NRC)=ID(NP+NRC)                                       ACV11567
         ENDIF                                                          ACV11568
C                                                                       ACV11569
C Match node has both children                                          ACV11570
C                                                                       ACV11571
      ELSE                                                              ACV11572
         NI=NQ                                                          ACV11573
         NH=NP                                                          ACV11574
         NG=ID(NP+NLC)                                                  ACV11575
         IF (ID(NP+NPR).GE.0) THEN                                      ACV11576
            NA=NH                                                       ACV11577
            NF=NI                                                       ACV11578
         ELSEIF ((ISHFT(ID(NP+NPR),-30).EQ.3).AND.                      ACV11579
     *           (ID(ID(NG+NRC)+NPR).GE.0)) THEN                        ACV11580
            NA=NH                                                       ACV11581
            NF=NI                                                       ACV11582
         ENDIF                                                          ACV11583
         IF (ID(NG+NLC).EQ.-1) THEN                                     ACV11584
            ID(ID(NG+NRC)+NRC)=NG                                       ACV11585
            ID(NG+NLC)=ID(NG+NRC)                                       ACV11586
            ID(NG+NRC)=ID(NP+NRC)                                       ACV11587
         ELSEIF (ID(ID(NG+NLC)+NRC).EQ.NG) THEN                         ACV11588
           IF (ISHFT(ID(NG+NPR),-30).EQ.2) THEN                         ACV11589
               ID(ID(NG+NLC)+NRC)=ID(NG+NRC)                            ACV11590
               ID(ID(NG+NRC)+NRC)=NG                                    ACV11591
               ID(NG+NRC)=ID(NP+NRC)                                    ACV11592
            ELSE                                                        ACV11593
               NH=NG                                                    ACV11594
               NG=ID(NG+NLC)                                            ACV11595
               ID(NH+NLC)=-1                                            ACV11596
               ID(ID(NH+NRC)+NRC)=NG                                    ACV11597
               ID(NG+NLC)=NH                                            ACV11598
               ID(NG+NRC)=ID(NP+NRC)                                    ACV11599
            ENDIF                                                       ACV11600
         ELSE                                                           ACV11601
150         NI=NH                                                       ACV11602
            NH=NG                                                       ACV11603
            NG=ID(ID(NG+NLC)+NRC)                                       ACV11604
            IF(ID(NH+NPR).GE.0) THEN                                    ACV11605
               NA=NH                                                    ACV11606
               NF=NI                                                    ACV11607
            ELSEIF ((ISHFT(ID(NH+NPR),-30).EQ.2)                        ACV11608
     *          .AND.(ID(ID(NH+NLC)+NPR).GE.0)) THEN                    ACV11609
               NA=NH                                                    ACV11610
               NF=NI                                                    ACV11611
            ENDIF                                                       ACV11612
            IF (ID(NG+NLC).EQ.-1) THEN                                  ACV11613
               ID(ID(NH+NLC)+NRC)=NH                                    ACV11614
            ELSEIF (ID(ID(NG+NLC)+NRC).EQ.NG) THEN                      ACV11615
               IF (ISHFT(ID(NG+NPR),-30).EQ.2) THEN                     ACV11616
                  ID(ID(NH+NLC)+NRC)=ID(NG+NLC)                         ACV11617
                  ID(ID(NG+NLC)+NRC)=NH                                 ACV11618
               ELSE                                                     ACV11619
                  NH=NG                                                 ACV11620
                  NG=ID(NG+NLC)                                         ACV11621
                  ID(NH+NLC)=-1                                         ACV11622
               ENDIF                                                    ACV11623
            ELSE                                                        ACV11624
               GOTO 150                                                 ACV11625
            ENDIF                                                       ACV11626
            ID(ID(ID(NP+NLC)+NRC)+NRC)=NG                               ACV11627
            ID(NG+NLC)=ID(NP+NLC)                                       ACV11628
            ID(NG+NRC)=ID(NP+NRC)                                       ACV11629
         ENDIF                                                          ACV11630
         IF (NQ.EQ.-1) THEN                                             ACV11631
            ID(0)=NG                                                    ACV11632
         ELSEIF (ID(NQ+NLC).EQ.NP) THEN                                 ACV11633
            ID(NQ+NLC)=NG                                               ACV11634
         ELSE                                                           ACV11635
            ID(ID(NQ+NLC)+NRC)=NG                                       ACV11636
         ENDIF                                                          ACV11637
         ID(NG+NPR)=IAND(ID(NG+NPR),NF0)+IAND(ID(NP+NPR),N3)            ACV11638
         IF (NA.EQ.NP) NA=NG                                            ACV11639
         IF (NF.EQ.NP) NF=NG                                            ACV11640
         DO 160 I=1,NKEY                                                ACV11641
160         ID(I)=ID(NG+I)                                              ACV11642
      ENDIF                                                             ACV11643
C                                                                       ACV11644
C Balance the tree                                                      ACV11645
C                                                                       ACV11646
500   CONTINUE                                                          ACV11647
      IF (ID(NA+NLC).EQ.-1) THEN                                        ACV11648
         ID(NA+NPR)=IAND(ID(NA+NPR),NF0)                                ACV11649
         GOTO 999                                                       ACV11650
      ENDIF                                                             ACV11651
      DO 180 I=1,NKEY                                                   ACV11652
         IF (ID(I).LT.ID(NA+I)) THEN                                    ACV11653
            GOTO 190                                                    ACV11654
         ELSEIF (ID(I).GT.ID(NA+I)) THEN                                ACV11655
            NO=ID(NA+NLC)                                               ACV11656
            IF (NO.EQ.-1) THEN                                          ACV11657
               NP=-1                                                    ACV11658
            ELSEIF (ID(NO+NRC).EQ.NA) THEN                              ACV11659
               IF (ISHFT(ID(NA+NPR),-30).EQ.2) THEN                     ACV11660
                  NP=-1                                                 ACV11661
               ELSE                                                     ACV11662
                  NP=NO                                                 ACV11663
               ENDIF                                                    ACV11664
            ELSE                                                        ACV11665
               NP=ID(NO+NRC)                                            ACV11666
            ENDIF                                                       ACV11667
            NB=NO                                                       ACV11668
            IND=N2                                                      ACV11669
            GOTO 300                                                    ACV11670
         ENDIF                                                          ACV11671
180   CONTINUE                                                          ACV11672
190   NO=ID(NA+NLC)                                                     ACV11673
      IF (NO.EQ.-1) THEN                                                ACV11674
         NP=-1                                                          ACV11675
      ELSEIF (ID(NO+NRC).EQ.NA) THEN                                    ACV11676
         IF (ISHFT(ID(NA+NPR),-30).EQ.2) THEN                           ACV11677
            NP=NO                                                       ACV11678
         ELSE                                                           ACV11679
            NP=-1                                                       ACV11680
         ENDIF                                                          ACV11681
      ELSE                                                              ACV11682
         NP=NO                                                          ACV11683
      ENDIF                                                             ACV11684
      IF (ID(NO+NRC).EQ.NA) THEN                                        ACV11685
         NB=NO                                                          ACV11686
      ELSE                                                              ACV11687
         NB=ID(NO+NRC)                                                  ACV11688
      ENDIF                                                             ACV11689
      IND=N3                                                            ACV11690
300   IF (ID(NA+NPR).GE.0) THEN                                         ACV11691
         ID(NA+NPR)=IAND(ID(NA+NPR),NF0)+IND                            ACV11692
         GOTO 600                                                       ACV11693
      ELSEIF (IAND(ID(NA+NPR),N3).NE.IND) THEN                          ACV11694
         ID(NA+NPR)=IAND(ID(NA+NPR),NF0)                                ACV11695
         GOTO 600                                                       ACV11696
      ENDIF                                                             ACV11697
C                                                                       ACV11698
C Rotate                                                                ACV11699
C                                                                       ACV11700
      IF (IND.EQ.N2) THEN                                               ACV11701
         IF (ISHFT(ID(NB+NPR),-30).EQ.3) THEN                           ACV11702
            IF (ID(ID(NB+NLC)+NRC).EQ.NB) THEN                          ACV11703
               NC=ID(NB+NLC)                                            ACV11704
            ELSE                                                        ACV11705
               NC=ID(ID(NB+NLC)+NRC)                                    ACV11706
            ENDIF                                                       ACV11707
            IF (ID(NC+NPR).GE.0) THEN                                   ACV11708
               IF (ID(NB+NRC).EQ.NA) THEN                               ACV11709
                  ID(NC+NRC)=ID(NA+NRC)                                 ACV11710
                  ID(NA+NLC)=-1                                         ACV11711
                  ID(NA+NRC)=NC                                         ACV11712
                  ID(NB+NLC)=-1                                         ACV11713
                  ID(NC+NLC)=NB                                         ACV11714
                  ID(NA+NPR)=IAND(ID(NA+NPR),NF0)                       ACV11715
                  ID(NB+NPR)=IAND(ID(NB+NPR),NF0)                       ACV11716
               ELSE                                                     ACV11717
                  IDCL=ID(NC+NLC)                                       ACV11718
                  ID(NA+NLC)=ID(IDCL+NRC)                               ACV11719
                  ID(ID(NA+NLC)+NRC)=ID(NB+NRC)                         ACV11720
                  ID(ID(NB+NLC)+NRC)=IDCL                               ACV11721
                  ID(IDCL+NRC)=NB                                       ACV11722
                  ID(NC+NLC)=NB                                         ACV11723
                  ID(NC+NRC)=ID(NA+NRC)                                 ACV11724
                  ID(NA+NRC)=NC                                         ACV11725
                  ID(NB+NRC)=NA                                         ACV11726
                  ID(NA+NPR)=IAND(ID(NA+NPR),NF0)                       ACV11727
                  ID(NB+NPR)=IAND(ID(NB+NPR),NF0)                       ACV11728
               ENDIF                                                    ACV11729
            ELSE                                                        ACV11730
               IDCL=ID(NC+NLC)                                          ACV11731
               IF (ID(IDCL+NRC).EQ.NC) THEN                             ACV11732
                  IF (ISHFT(ID(NC+NPR),-30).EQ.2) THEN                  ACV11733
                     ID(NA+NLC)=ID(NB+NRC)                              ACV11734
                     ID(ID(NB+NLC)+NRC)=IDCL                            ACV11735
                     ID(IDCL+NRC)=NB                                    ACV11736
                  ELSE                                                  ACV11737
                     ID(NA+NLC)=IDCL                                    ACV11738
                     ID(ID(NA+NLC)+NRC)=ID(NB+NRC)                      ACV11739
                     ID(ID(NB+NLC)+NRC)=NB                              ACV11740
                  ENDIF                                                 ACV11741
               ELSE                                                     ACV11742
                  ID(NA+NLC)=ID(IDCL+NRC)                               ACV11743
                  ID(ID(NA+NLC)+NRC)=ID(NB+NRC)                         ACV11744
                  ID(ID(NB+NLC)+NRC)=IDCL                               ACV11745
                  ID(IDCL+NRC)=NB                                       ACV11746
               ENDIF                                                    ACV11747
               ID(NC+NLC)=NB                                            ACV11748
               ID(NC+NRC)=ID(NA+NRC)                                    ACV11749
               ID(NA+NRC)=NC                                            ACV11750
               ID(NB+NRC)=NA                                            ACV11751
               IF (ISHFT(ID(NC+NPR),-30).EQ.2) THEN                     ACV11752
                  ID(NA+NPR)=IAND(ID(NA+NPR),NF0)+N3                    ACV11753
                  ID(NB+NPR)=IAND(ID(NB+NPR),NF0)                       ACV11754
               ELSE                                                     ACV11755
                  ID(NA+NPR)=IAND(ID(NA+NPR),NF0)                       ACV11756
                  ID(NB+NPR)=IAND(ID(NB+NPR),NF0)+N2                    ACV11757
               ENDIF                                                    ACV11758
               ID(NC+NPR)=IAND(ID(NC+NPR),NF0)                          ACV11759
            ENDIF                                                       ACV11760
            NB=NC                                                       ACV11761
         ELSE                                                           ACV11762
            IF (ID(NB+NRC).EQ.NA) THEN                                  ACV11763
               IF (ID(NB+NPR).GE.0) THEN                                ACV11764
                  NC=ID(ID(NB+NLC)+NRC)                                 ACV11765
                  ID(NA+NLC)=NC                                         ACV11766
                  ID(NC+NRC)=NA                                         ACV11767
               ELSE                                                     ACV11768
                  ID(NA+NLC)=-1                                         ACV11769
               ENDIF                                                    ACV11770
               ID(NB+NRC)=ID(NA+NRC)                                    ACV11771
               ID(ID(NB+NLC)+NRC)=NA                                    ACV11772
               ID(NA+NRC)=NB                                            ACV11773
            ELSE                                                        ACV11774
               NC=ID(ID(NB+NLC)+NRC)                                    ACV11775
               ID(NC+NRC)=ID(NB+NRC)                                    ACV11776
               ID(NB+NRC)=ID(NA+NRC)                                    ACV11777
               ID(ID(NB+NLC)+NRC)=NA                                    ACV11778
               ID(NA+NLC)=NC                                            ACV11779
               ID(NA+NRC)=NB                                            ACV11780
            ENDIF                                                       ACV11781
            IF (ID(NB+NPR).GE.0) THEN                                   ACV11782
               ID(NA+NPR)=IAND(ID(NA+NPR),NF0)+N2                       ACV11783
               ID(NB+NPR)=IAND(ID(NB+NPR),NF0)+N3                       ACV11784
            ELSE                                                        ACV11785
               ID(NA+NPR)=IAND(ID(NA+NPR),NF0)                          ACV11786
               ID(NB+NPR)=IAND(ID(NB+NPR),NF0)                          ACV11787
            ENDIF                                                       ACV11788
         ENDIF                                                          ACV11789
      ELSE                                                              ACV11790
         IF (ISHFT(ID(NB+NPR),-30).EQ.2) THEN                           ACV11791
            NC=ID(NB+NLC)                                               ACV11792
            IF (ID(NC+NPR).GE.0) THEN                                   ACV11793
               IF (ID(NA+NLC).EQ.NB) THEN                               ACV11794
                  ID(NC+NRC)=ID(NA+NRC)                                 ACV11795
                  ID(NC+NLC)=NA                                         ACV11796
                  ID(NA+NLC)=-1                                         ACV11797
                  ID(NA+NRC)=NB                                         ACV11798
                  ID(NB+NLC)=-1                                         ACV11799
                  ID(NB+NRC)=NC                                         ACV11800
                  ID(NA+NPR)=IAND(ID(NA+NPR),NF0)                       ACV11801
                  ID(NB+NPR)=IAND(ID(NB+NPR),NF0)                       ACV11802
               ELSE                                                     ACV11803
                  IDCL=ID(NC+NLC)                                       ACV11804
                  ID(ID(NA+NLC)+NRC)=IDCL                               ACV11805
                  ID(NB+NLC)=ID(IDCL+NRC)                               ACV11806
                  ID(ID(NB+NLC)+NRC)=ID(NC+NRC)                         ACV11807
                  ID(IDCL+NRC)=NA                                       ACV11808
                  ID(NC+NLC)=NA                                         ACV11809
                  ID(NC+NRC)=ID(NA+NRC)                                 ACV11810
                  ID(NA+NRC)=NB                                         ACV11811
                  ID(NB+NRC)=NC                                         ACV11812
                  ID(NA+NPR)=IAND(ID(NA+NPR),NF0)                       ACV11813
                  ID(NB+NPR)=IAND(ID(NB+NPR),NF0)                       ACV11814
               ENDIF                                                    ACV11815
            ELSE                                                        ACV11816
               IDCL=ID(NC+NLC)                                          ACV11817
               IF (ID(IDCL+NRC).EQ.NC) THEN                             ACV11818
                  IF (ISHFT(ID(NC+NPR),-30).EQ.3) THEN                  ACV11819
                     ID(ID(NA+NLC)+NRC)=NA                              ACV11820
                     ID(IDCL+NRC)=ID(NC+NRC)                            ACV11821
                     ID(NB+NLC)=IDCL                                    ACV11822
                  ELSE                                                  ACV11823
                     ID(ID(NA+NLC)+NRC)=IDCL                            ACV11824
                     ID(IDCL+NRC)=NA                                    ACV11825
                     ID(NB+NLC)=ID(NC+NRC)                              ACV11826
                  ENDIF                                                 ACV11827
               ELSE                                                     ACV11828
                  ID(ID(NA+NLC)+NRC)=IDCL                               ACV11829
                  ID(NB+NLC)=ID(IDCL+NRC)                               ACV11830
                  ID(ID(NB+NLC)+NRC)=ID(NC+NRC)                         ACV11831
                  ID(IDCL+NRC)=NA                                       ACV11832
               ENDIF                                                    ACV11833
               ID(NC+NLC)=NA                                            ACV11834
               ID(NC+NRC)=ID(NA+NRC)                                    ACV11835
               ID(NA+NRC)=NB                                            ACV11836
               ID(NB+NRC)=NC                                            ACV11837
               IF (ISHFT(ID(NC+NPR),-30).EQ.3) THEN                     ACV11838
                  ID(NA+NPR)=IAND(ID(NA+NPR),NF0)+N2                    ACV11839
                  ID(NB+NPR)=IAND(ID(NB+NPR),NF0)                       ACV11840
               ELSE                                                     ACV11841
                  ID(NA+NPR)=IAND(ID(NA+NPR),NF0)                       ACV11842
                  ID(NB+NPR)=IAND(ID(NB+NPR),NF0)+N3                    ACV11843
               ENDIF                                                    ACV11844
               ID(NC+NPR)=IAND(ID(NC+NPR),NF0)                          ACV11845
            ENDIF                                                       ACV11846
            NB=NC                                                       ACV11847
         ELSE                                                           ACV11848
            IF (ID(NA+NLC).EQ.NB) THEN                                  ACV11849
               IF (ID(NB+NPR).GE.0) THEN                                ACV11850
                  NC=ID(NB+NLC)                                         ACV11851
                  ID(NB+NRC)=ID(NA+NRC)                                 ACV11852
                  ID(NA+NLC)=NC                                         ACV11853
                  ID(NA+NRC)=ID(NC+NRC)                                 ACV11854
                  ID(NC+NRC)=NA                                         ACV11855
               ELSE                                                     ACV11856
                  ID(NA+NLC)=-1                                         ACV11857
                  ID(NB+NRC)=ID(NA+NRC)                                 ACV11858
                  ID(NA+NRC)=ID(NB+NLC)                                 ACV11859
               ENDIF                                                    ACV11860
               ID(NB+NLC)=NA                                            ACV11861
            ELSE                                                        ACV11862
               NC=ID(NB+NLC)                                            ACV11863
               ID(NB+NLC)=NA                                            ACV11864
               ID(NB+NRC)=ID(NA+NRC)                                    ACV11865
               ID(ID(NA+NLC)+NRC)=NC                                    ACV11866
               ID(NA+NRC)=ID(NC+NRC)                                    ACV11867
               ID(NC+NRC)=NA                                            ACV11868
            ENDIF                                                       ACV11869
            IF (ID(NB+NPR).GE.0) THEN                                   ACV11870
               ID(NA+NPR)=IAND(ID(NA+NPR),NF0)+N3                       ACV11871
               ID(NB+NPR)=IAND(ID(NB+NPR),NF0)+N2                       ACV11872
            ELSE                                                        ACV11873
               ID(NA+NPR)=IAND(ID(NA+NPR),NF0)                          ACV11874
               ID(NB+NPR)=IAND(ID(NB+NPR),NF0)                          ACV11875
            ENDIF                                                       ACV11876
         ENDIF                                                          ACV11877
      ENDIF                                                             ACV11878
      IF (NF.EQ.-1) THEN                                                ACV11879
         ID(0)=NB                                                       ACV11880
      ELSEIF (NA.EQ.ID(NF+NLC)) THEN                                    ACV11881
         ID(NF+NLC)=NB                                                  ACV11882
      ELSE                                                              ACV11883
         ID(ID(NF+NLC)+NRC)=NB                                          ACV11884
      ENDIF                                                             ACV11885
600   NF=NA                                                             ACV11886
      NA=NP                                                             ACV11887
      IF (NP.NE.-1) GOTO 500                                            ACV11888
999   CONTINUE                                                          ACV11889
C                                                                       ACV11890
C Reconstruct the priority array                                        ACV11891
C                                                                       ACV11892
      NA=0                                                              ACV11893
      NL=ID(-1)-NRC                                                     ACV11894
      NLPR=IAND(ID(NL+NPR),NF0)                                         ACV11895
1000  NB=NA+NA+NRC                                                      ACV11896
      IF (NB.LT.NL)THEN                                                 ACV11897
         NC=NB+NRC                                                      ACV11898
         NBPR=IAND(ID(NB+NPR),NF0)                                      ACV11899
         NCPR=IAND(ID(NC+NPR),NF0)                                      ACV11900
         IF ((NBPR.LE.NCPR).OR.(NC.GT.NL)) THEN                         ACV11901
            IF (NBPR.LT.NLPR) THEN                                      ACV11902
               NBRC=ID(NB+NRC)                                          ACV11903
               NBLC=ID(NB+NLC)                                          ACV11904
               IF (NBRC.EQ.-1) THEN                                     ACV11905
                  ID(0)=NA                                              ACV11906
               ELSE                                                     ACV11907
                  IF (ID(NBRC+NLC).EQ.NB) THEN                          ACV11908
                     ID(NBRC+NLC)=NA                                    ACV11909
                  ELSE                                                  ACV11910
                     IF (ID(ID(NBRC+NRC)+NLC).EQ.NB) THEN               ACV11911
                        ID(ID(NBRC+NRC)+NLC)=NA                         ACV11912
                     ELSE                                               ACV11913
                        ID(ID(NBRC+NLC)+NRC)=NA                         ACV11914
                     ENDIF                                              ACV11915
                  ENDIF                                                 ACV11916
               ENDIF                                                    ACV11917
               IF (NBLC.NE.-1) THEN                                     ACV11918
                  IF (ID(NBLC+NRC).EQ.NB) THEN                          ACV11919
                     ID(NBLC+NRC)=NA                                    ACV11920
                  ELSE                                                  ACV11921
                     ID(ID(NBLC+NRC)+NRC)=NA                            ACV11922
                  ENDIF                                                 ACV11923
               ENDIF                                                    ACV11924
               DO 1400 I=1,NRC                                          ACV11925
                  ID(NA+I)=ID(NB+I)                                     ACV11926
1400           CONTINUE                                                 ACV11927
               NA=NB                                                    ACV11928
               GOTO 1000                                                ACV11929
            ENDIF                                                       ACV11930
         ELSE                                                           ACV11931
            IF (NCPR.LT.NLPR) THEN                                      ACV11932
               NCRC=ID(NC+NRC)                                          ACV11933
               NCLC=ID(NC+NLC)                                          ACV11934
               IF (NCRC.EQ.-1) THEN                                     ACV11935
                  ID(0)=NA                                              ACV11936
               ELSE                                                     ACV11937
                  IF (ID(NCRC+NLC).EQ.NC) THEN                          ACV11938
                     ID(NCRC+NLC)=NA                                    ACV11939
                  ELSE                                                  ACV11940
                     IF (ID(ID(NCRC+NRC)+NLC).EQ.NC) THEN               ACV11941
                        ID(ID(NCRC+NRC)+NLC)=NA                         ACV11942
                     ELSE                                               ACV11943
                        ID(ID(NCRC+NLC)+NRC)=NA                         ACV11944
                     ENDIF                                              ACV11945
                  ENDIF                                                 ACV11946
               ENDIF                                                    ACV11947
               IF (NCLC.NE.-1) THEN                                     ACV11948
                  IF (ID(NCLC+NRC).EQ.NC) THEN                          ACV11949
                     ID(NCLC+NRC)=NA                                    ACV11950
                  ELSE                                                  ACV11951
                     ID(ID(NCLC+NRC)+NRC)=NA                            ACV11952
                  ENDIF                                                 ACV11953
               ENDIF                                                    ACV11954
               DO 2400 I=1,NRC                                          ACV11955
                  ID(NA+I)=ID(NC+I)                                     ACV11956
2400           CONTINUE                                                 ACV11957
               NA=NC                                                    ACV11958
               GOTO 1000                                                ACV11959
            ENDIF                                                       ACV11960
         ENDIF                                                          ACV11961
      ENDIF                                                             ACV11962
      NLRC=ID(NL+NRC)                                                   ACV11963
      NLLC=ID(NL+NLC)                                                   ACV11964
      IF (NLRC.EQ.-1) THEN                                              ACV11965
         ID(0)=NA                                                       ACV11966
      ELSE                                                              ACV11967
         IF (ID(NLRC+NLC).EQ.NL) THEN                                   ACV11968
            ID(NLRC+NLC)=NA                                             ACV11969
         ELSE                                                           ACV11970
            IF (ID(ID(NLRC+NRC)+NLC).EQ.NL) THEN                        ACV11971
               ID(ID(NLRC+NRC)+NLC)=NA                                  ACV11972
            ELSE                                                        ACV11973
               ID(ID(NLRC+NLC)+NRC)=NA                                  ACV11974
            ENDIF                                                       ACV11975
         ENDIF                                                          ACV11976
      ENDIF                                                             ACV11977
      IF (NLLC.NE.-1) THEN                                              ACV11978
         IF (ID(NLLC+NRC).EQ.NL) THEN                                   ACV11979
            ID(NLLC+NRC)=NA                                             ACV11980
         ELSE                                                           ACV11981
            ID(ID(NLLC+NRC)+NRC)=NA                                     ACV11982
         ENDIF                                                          ACV11983
      ENDIF                                                             ACV11984
      DO 7000 I=1,NRC                                                   ACV11985
         ID(NA+I)=ID(NL+I)                                              ACV11986
         ID(NL+I)=0                                                     ACV11987
7000  CONTINUE                                                          ACV11988
      ID(-1)=NL                                                         ACV11989
9999  RETURN                                                            ACV11990
      END                                                               ACV11991
C ----------------------------------------------------------------------ACV11992
C                                                                       ACV11993
C                           **************                              ACV11994
C                           ***  TOUT  ***                              ACV11995
C                           **************                              ACV11996
C                                                                       ACV11997
C The subroutine TOUT traverses the weighted search tree in order to    ACV11998
C specific nodes and sends the outputs to an output devise designated   ACV11999
C by NFILE. The traversal is done in ascending order if NAD=1 and in    ACV12000
C descending order otherwise. A description of the arguments follows:   ACV12001
C                                                                       ACV12002
C   NFILE = File number assigned to output device where the results are ACV12003
C           to be written. In the special case NFILE = 0, TOUT does not ACV12004
C           output anything, however after upon returning ID(-5) points ACV12005
C           to the current node position.                               ACV12006
C     NAD = 1 if the traversal is to be in ascending order, otherwise   ACV12007
C           it is done in descending order.                             ACV12008
C     MIN = Number of the first node that is to be retrieved.           ACV12009
C     MAX = Number of the last node that is to be retrieved.            ACV12010
C   NSTEP = Step size; nodes that are integer multiples of NSTEP beyond ACV12011
C           MIN up to MAX will be retrieved.                            ACV12012
C      ID = Name of the tree array that is to be traversed.             ACV12013
C                                                                       ACV12014
C An internal array called STACK is used to store information about the ACV12015
C path traversed. It is dimensioned 32 since this is the maximum height ACV12016
C the tree can have.  This limit is set by the fact that 2**32-1 is the ACV12017
C largest integer pointer that can be used on a 32 bit machine.         ACV12018
C                                                                       ACV12019
C ----------------------------------------------------------------------ACV12020
C                                                                       ACV12021
      SUBROUTINE TOUT(NFILE,NAD,MIN,MAX,NSTEP,ID)                       ACV12022
C                                                                       ACV12023
      INTEGER ID(-10:*),STACK(32)                                       ACV12024
      NSUM=ID(-3)                                                       ACV12025
      NLC=NSUM+2                                                        ACV12026
      NRC=NLC+1                                                         ACV12027
      I=0                                                               ACV12028
      M=0                                                               ACV12029
      NUM=MIN                                                           ACV12030
      NODE=ID(0)                                                        ACV12031
      IF (NAD.EQ.1) THEN                                                ACV12032
200      IF (NODE.NE.-1) THEN                                           ACV12033
            I=I+1                                                       ACV12034
            STACK(I)=NODE                                               ACV12035
            NO=ID(NODE+NLC)                                             ACV12036
            IF (NO.EQ.-1) THEN                                          ACV12037
               NODE=-1                                                  ACV12038
            ELSEIF (ID(NO+NRC).EQ.NODE) THEN                            ACV12039
               IF (ISHFT(ID(NODE+NPR),-30).EQ.2) THEN                   ACV12040
                  NODE=NO                                               ACV12041
               ELSE                                                     ACV12042
                  NODE=-1                                               ACV12043
               ENDIF                                                    ACV12044
            ELSE                                                        ACV12045
               NODE=NO                                                  ACV12046
            ENDIF                                                       ACV12047
            GOTO 200                                                    ACV12048
300         M=M+1                                                       ACV12049
            ID(-5)=NODE                                                 ACV12050
            IF ((M.EQ.NUM).OR.(M.EQ.MAX)) THEN                          ACV12051
               IF (NFILE.GT.0) WRITE(NFILE,1000)                        ACV12052
     *           M,(ID(NODE+II),II=1,NSUM)                              ACV12053
               IF (M.EQ.MAX) GOTO 999                                   ACV12054
               NUM=NUM+NSTEP                                            ACV12055
            ENDIF                                                       ACV12056
            NO=ID(NODE+NLC)                                             ACV12057
            IF (NO.EQ.-1) THEN                                          ACV12058
               NODE=-1                                                  ACV12059
            ELSEIF (ID(NO+NRC).EQ.NODE) THEN                            ACV12060
               IF (ISHFT(ID(NODE+NPR),-30).EQ.2) THEN                   ACV12061
                  NODE=-1                                               ACV12062
               ELSE                                                     ACV12063
                  NODE=NO                                               ACV12064
               ENDIF                                                    ACV12065
            ELSE                                                        ACV12066
               NODE=ID(NO+NRC)                                          ACV12067
            ENDIF                                                       ACV12068
            GOTO 200                                                    ACV12069
         ENDIF                                                          ACV12070
         IF (I.NE.0) THEN                                               ACV12071
            NODE=STACK(I)                                               ACV12072
            I=I-1                                                       ACV12073
            GOTO 300                                                    ACV12074
         ENDIF                                                          ACV12075
      ELSE                                                              ACV12076
400      IF (NODE.NE.-1) THEN                                           ACV12077
            I=I+1                                                       ACV12078
            STACK(I)=NODE                                               ACV12079
            NO=ID(NODE+NLC)                                             ACV12080
            IF (NO.EQ.-1) THEN                                          ACV12081
               NODE=-1                                                  ACV12082
            ELSEIF (ID(NO+NRC).EQ.NODE) THEN                            ACV12083
               IF (ISHFT(ID(NODE+NPR),-30).EQ.2) THEN                   ACV12084
                  NODE=-1                                               ACV12085
               ELSE                                                     ACV12086
                  NODE=NO                                               ACV12087
               ENDIF                                                    ACV12088
            ELSE                                                        ACV12089
               NODE=ID(NO+NRC)                                          ACV12090
            ENDIF                                                       ACV12091
            GOTO 400                                                    ACV12092
500         M=M+1                                                       ACV12093
            ID(-5)=NODE                                                 ACV12094
            IF ((M.EQ.NUM).OR.(M.EQ.MAX)) THEN                          ACV12095
               IF (NFILE.GT.0) WRITE(NFILE,1000)                        ACV12096
     *           M,(ID(NODE+II),II=1,NSUM)                              ACV12097
               IF (M.EQ.MAX) GOTO 999                                   ACV12098
               NUM=NUM+NSTEP                                            ACV12099
            ENDIF                                                       ACV12100
            NO=ID(NODE+NLC)                                             ACV12101
            IF (NO.EQ.-1) THEN                                          ACV12102
               NODE=-1                                                  ACV12103
            ELSEIF (ID(NO+NRC).EQ.NODE) THEN                            ACV12104
               IF (ISHFT(ID(NODE+NPR),-30).EQ.2) THEN                   ACV12105
                  NODE=NO                                               ACV12106
               ELSE                                                     ACV12107
                  NODE=-1                                               ACV12108
               ENDIF                                                    ACV12109
            ELSE                                                        ACV12110
               NODE=NO                                                  ACV12111
            ENDIF                                                       ACV12112
            GOTO 400                                                    ACV12113
         ENDIF                                                          ACV12114
         IF (I.NE.0) THEN                                               ACV12115
            NODE=STACK(I)                                               ACV12116
            I=I-1                                                       ACV12117
            GOTO 500                                                    ACV12118
         ENDIF                                                          ACV12119
      ENDIF                                                             ACV12120
      WRITE(NFILE,*)                                                    ACV12121
      WRITE(NFILE,*)'  WARNING:  TOUT SEARCH HAS GONE OUT OF BOUNDS!'   ACV12122
999   RETURN                                                            ACV12123
1000  FORMAT(3X,I6,':',2X,14I5/(12X,14I5))                              ACV12124
      END                                                               ACV12125
C ----------------------------------------------------------------------ACV12126
C                                                                       ACV12127
C                            ************                               ACV12128
C                            *** TMRG ***                               ACV12129
C                            ************                               ACV12130
C                                                                       ACV12131
C The subroutine TMRG allows the user to merge two trees that have the  ACV12132
C same structure. Elements of the first tree are stored in the second.  ACV12133
C If the merge results in an overflow condition, notification is given  ACV12134
C and the lowest priority item is eliminated in favor of the incoming   ACV12135
C element with higher priority. If the trees have a different structure,ACV12136
C the merge is automatically aborted (alternate return) after giving an ACV12137
C appropriate warning message. If the element already exists in the     ACV12138
C second tree, the priority is set to the maximum of the two.           ACV12139
C                                                                       ACV12140
C   RETURN 1 --> Unsuccessful match (the trees have different           ACV12141
C                structures) so merge is aborted.                       ACV12142
C                                                                       ACV12143
C ----------------------------------------------------------------------ACV12144
C                                                                       ACV12145
      SUBROUTINE TMRG(ID1,BUFF1,LLBF1,ID2,BUFF2,LLBF2,*)                ACV12146
C                                                                       ACV12147
      INTEGER ID1(-10:*),LLBF1(-2:*),ID2(-10:*),LLBF2(-2:*)             ACV12148
      REAL    BUFF1(*),BUFF2(*)                       ! *** storage typeACV12149
C     REAL*8  BUFF1(*),BUFF2(*)                       ! *** storage typeACV12150
      PARAMETER (MAXDAT=1000)                                           ACV12151
      INTEGER NDAT(2)                                                   ACV12152
      REAL    BULOAD(MAXDAT)                          ! *** temp storageACV12153
C     REAL*8  BULOAD(MAXDAT)                          ! *** temp storageACV12154
C     DATA    NF0/Z3FFFFFFF/                                            ACV12155
      DATA    NF0/1073741823/                                           ACV12156
C                                                                       ACV12157
C Check the structure of the trees for compatibility                    ACV12158
C                                                                       ACV12159
      IF (ID1(-4).NE.ID2(-4).OR.ID1(-3).NE.ID2(-3)) THEN                ACV12160
         WRITE(6,*) ' TREE 1 AND 2 HAVE DIFFERENT STRUCTURES'           ACV12161
         RETURN 1                                                       ACV12162
      ENDIF                                                             ACV12163
      IF (ID1(-10)+ID2(-10).GE.ID2(-9).OR.LLBF1(-1).GE.LLBF2(0)) THEN   ACV12164
         WRITE(6,*) ' SIZE OF THE INCOMING NODE IS LARGER THAN MAXIMUM,'ACV12165
         WRITE(6,*) ' ONE OR MORE LOW PRIORITY NODES WILL BE DELETED.'  ACV12166
      ENDIF                                                             ACV12167
C                                                                       ACV12168
C List out information from tree 1                                      ACV12169
C                                                                       ACV12170
      DO INODE1=1,ID1(-10)                                              ACV12171
         CALL TOUT(0,0,INODE1,INODE1,1,ID1)                             ACV12172
C                                                                       ACV12173
C ...find the key and data indices in first tree                        ACV12174
C                                                                       ACV12175
         IPNODE=ID1(-5)                                                 ACV12176
         IPKEY=IPNODE+1                                                 ACV12177
         IPRIOR=ID1(IPNODE+ID1(-2))                                     ACV12178
C                                                                       ACV12179
C Check whether or not the incoming node is new                         ACV12180
C                                                                       ACV12181
         CALL TCHK(ID1(IPKEY),ID2,*100)                                 ACV12182
C                                                                       ACV12183
C ...find the data information from first tree                          ACV12184
C                                                                       ACV12185
         IPDATA=IPNODE+ID1(-4)+1                                        ACV12186
         IFIND=ID1(IPDATA)                                              ACV12187
         NOSIZE=ID1(IPDATA+1)                                           ACV12188
C                                                                       ACV12189
C ...abort if the size is too big                                       ACV12190
C                                                                       ACV12191
         IF (NOSIZE.GT.MAXDAT) THEN                                     ACV12192
            WRITE(6,*) ' SIZE OF TEMPORARY BUFFER BULOAD IS TOO SMALL.' ACV12193
            WRITE(6,*) ' INCREASE MAXDAT IN SUBROUTINE TMRG AND'        ACV12194
            WRITE(6,*) ' CONTINUE.'                                     ACV12195
            STOP                                                        ACV12196
         ENDIF                                                          ACV12197
C                                                                       ACV12198
         DO IKK=1,NOSIZE                                                ACV12199
            BULOAD(IKK)=BUFF1(IFIND)                                    ACV12200
            IFIND=LLBF1(IFIND)                                          ACV12201
         ENDDO                                                          ACV12202
C                                                                       ACV12203
C Store key, data, and information in second tree                       ACV12204
C                                                                       ACV12205
         CALL TADD(ID1(IPKEY),NDAT,BULOAD,NOSIZE,IPRIOR,ID2,BUFF2,LLBF2)ACV12206
         GOTO 200                                                       ACV12207
C                                                                       ACV12208
C Update the priority if the incoming node is old                       ACV12209
C                                                                       ACV12210
100      IFIND=ID2(-5)+ID2(-2)                                          ACV12211
         ID2(IFIND)=MAX0(ID2(IFIND),IPRIOR)                             ACV12212
C                                                                       ACV12213
200      CONTINUE                                                       ACV12214
      ENDDO                                                             ACV12215
      RETURN                                                            ACV12216
      END                                                               ACV12217
C ----------------------------------------------------------------------ACV12218
C                                                                       ACV12219
C                       ***********************                         ACV12220
C                       ***   SU3 PACKAGE   ***                         ACV12221
C                       ***********************                         ACV12222
C                                                                       ACV12223
C ----------------------------------------------------------------------ACV12224
      SUBROUTINE BLOCKS                                                 ACV12225
C     ------------------------------------------------------------------ACV12226
C     BINOMIAL COEFFICIENTS AND FACTORIALS ***** SEE COMMENT BELOW *****ACV12227
C     ------------------------------------------------------------------ACV12228
C     UPDATE/MOD: (MTS,06-76)  H.SATO             EXPANDED RANGE        ACV12229
C                 (LSU,11-78)  J.P.DRAAYER        LOG FACTORIALS        ACV12230
C                 (LSU,05-80)  J.P.DRAAYER        LOG BINOMIALS         ACV12231
C                 (LSU,08-81)  J.P.DRAAYER        EXPANDED RANGE        ACV12232
C                 (LSU,01-83)  J.P.DRAAYER        EXTENDED PRECISION    ACV12233
C                 (LSU,03-83)  J.P.DRAAYER        D,Q MIX & LOG CUTS    ACV12234
C                 (LSU,11-84)  J.P.DRAAYER        MODIFIED FACTORIALS   ACV12235
C                 (LSU,01-88)  J.P.DRAAYER        DLOGF RANGE/INDEX     ACV12236
C                 (LSU,10-89)  J.P.DRAAYER        BINOMIAL INVERSES     ACV12237
C                 (LSU,11-89)  J.P.DRAAYER        POWERS OF TWO ARRAY   ACV12238
C                                                                       ACV12239
C     BKDB--BINOMIAL (BINO) COEFFICIENTS (EXPANDED 6-76,8-81,1-83)      ACV12240
C           SCALE: BINO(I,J)=DBINO(I*(I+1)/2+J+1)                       ACV12241
C           RANGE: BINO(0,0)=DBINO(1) TO BINO(128,128)=DBINO(8385)      ACV12242
C           ADDED: 2**I = DTWOS(I) WHERE I=-128,128                     ACV12243
C     BKQB--BINOMIAL (BINO) COEFFICIENTS (EXPANDED 6-76,8-81,1-83)      ACV12244
C           SCALE: BINO(I,J)=QBINO(I*(I+1)/2+J+1)                       ACV12245
C           RANGE: BINO(0,0)=QBINO(1) TO BINO(192,192)=QBINO(18721)     ACV12246
C           ADDED: 2**I = QTWOS(I) WHERE I=-192,192                     ACV12247
C     BKDF--LOG FACTORIALS (FACT) (INSERTED 11-78, MODIFIED 01-88)      ACV12248
C           SCALE: LNFACT(I)=DLOGF(2*I)                                 ACV12249
C           RANGE: LNFACT(0)=DLOGF(0) TO LNFACT(1000)=DLOGF(2000)       ACV12250
C                                                                       ACV12251
C          ********************************************************     ACV12252
C          **  BLOCKS INPUT MUST BE PREGENERATED USING SU3GENBK  **     ACV12253
C          ********************************************************     ACV12254
C                                                                       ACV12255
C     ------------------------------------------------------------------ACV12256
      IMPLICIT REAL*8(D),REAL*16(Q)                                     ACV12257
      COMMON/BKDB/DBINO(8385),DBINV(8385),DTWOS(-128:128)               ACV12258
      COMMON/BKQB/QBINO(18721),QBINV(18721),QTWOS(-192:192)             ACV12259
      COMMON/BKDF/DLOGF(0:2000)                                         ACV12260
      READ(4)DLOGF                                                      ACV12261
      READ(4)DBINO                                                      ACV12262
      READ(4)DBINV                                                      ACV12263
      READ(4)QBINO                                                      ACV12264
      READ(4)QBINV                                                      ACV12265
      READ(4)DTWOS                                                      ACV12266
      READ(4)QTWOS                                                      ACV12267
      RETURN                                                            ACV12268
      END                                                               ACV12269
      FUNCTION MULTU3(LX1,MX1,LX2,MX2,LX3,MX3)                          ACV12270
C     ------------------------------------------------------------------ACV12271
C     MULTIPLICITY IN U3 COUPLING (SEE MULTTEST FOR VARIOUS VERSIONS)   ACV12272
C     ... FASTER THAN THE DRAAYER ORIGINAL AND THE MILLENER KAS FUNCTIONACV12273
C     ------------------------------------------------------------------ACV12274
C     UPDATE/MOD: (LSU,11-78)  J.P.DRAAYER        ORIGINAL VERSION      ACV12275
C                 (BNL,06-87)  J.MILLENER         MILLENER VERSION      ACV12276
C                 (LSU,10-89)  J.P.DRAAYER        PROFILER OPTIMIZED    ACV12277
C     ------------------------------------------------------------------ACV12278
      MULTU3=0                                                          ACV12279
      NX=LX1+LX2-LX3-MX1-MX2+MX3                                        ACV12280
      MX=NX/3                                                           ACV12281
      IF(3*MX.NE.NX)RETURN                                              ACV12282
      IF(MX.GE.0)THEN                                                   ACV12283
         L1=LX1                                                         ACV12284
         L2=LX2                                                         ACV12285
         L3=LX3                                                         ACV12286
         M1=MX1                                                         ACV12287
         M2=MX2                                                         ACV12288
         M3=MX3                                                         ACV12289
      ELSE                                                              ACV12290
         L1=MX1                                                         ACV12291
         L2=MX2                                                         ACV12292
         L3=MX3                                                         ACV12293
         M1=LX1                                                         ACV12294
         M2=LX2                                                         ACV12295
         M3=LX3                                                         ACV12296
         MX=-MX                                                         ACV12297
      ENDIF                                                             ACV12298
      NX=MX+M1+M2-M3                                                    ACV12299
      MU=MIN0(L1-MX,M2)                                                 ACV12300
      IF(MU.LT.0)RETURN                                                 ACV12301
      NU=MIN0(L2-MX,M1)                                                 ACV12302
      IF(NU.LT.0)RETURN                                                 ACV12303
      MULTU3=MAX0(MIN0(NX,NU)-MAX0(NX-MU,0)+1,0)                        ACV12304
      RETURN                                                            ACV12305
      END                                                               ACV12306
      SUBROUTINE U3MULT(LX1,MX1,LX2,MX2,LX3,MX3,MULTU3,*)               ACV12307
C     ------------------------------------------------------------------ACV12308
C     MULTIPLICITY IN U3 COUPLING (SEE MULTTEST FOR VARIOUS VERSIONS)   ACV12309
C     ... MULTU3 IS FASTER THAN THE DRAAYER ORIGINAL AND MILLENER KAS   ACV12310
C     ... SUBROUTINE FORM EXECUTES FASTEST ... RETURN 1 BRANCH OPTION   ACV12311
C     ------------------------------------------------------------------ACV12312
C     UPDATE/MOD: (LSU,11-78)  J.P.DRAAYER     ORIGINAL FUNCTION FORM   ACV12313
C                 (BNL,06-87)  J.MILLENER      MILLENER FUNCTION ... KASACV12314
C                 (LSU,10-89)  J.P.DRAAYER     PROFILER OPTIMIZED FORM  ACV12315
C                 (LSU,11-89)  J.P.DRAAYER     SUBROUTINE FORM FOR SPEEDACV12316
C     ------------------------------------------------------------------ACV12317
      MULTU3=0                                                          ACV12318
      NX=LX1+LX2-LX3-MX1-MX2+MX3                                        ACV12319
      MX=NX/3                                                           ACV12320
      IF(3*MX.NE.NX)RETURN 1                                            ACV12321
      IF(MX.GE.0)THEN                                                   ACV12322
         L1=LX1                                                         ACV12323
         L2=LX2                                                         ACV12324
         L3=LX3                                                         ACV12325
         M1=MX1                                                         ACV12326
         M2=MX2                                                         ACV12327
         M3=MX3                                                         ACV12328
      ELSE                                                              ACV12329
         L1=MX1                                                         ACV12330
         L2=MX2                                                         ACV12331
         L3=MX3                                                         ACV12332
         M1=LX1                                                         ACV12333
         M2=LX2                                                         ACV12334
         M3=LX3                                                         ACV12335
         MX=-MX                                                         ACV12336
      ENDIF                                                             ACV12337
      NX=MX+M1+M2-M3                                                    ACV12338
      MU=MIN0(L1-MX,M2)                                                 ACV12339
      IF(MU.LT.0)RETURN 1                                               ACV12340
      NU=MIN0(L2-MX,M1)                                                 ACV12341
      IF(NU.LT.0)RETURN 1                                               ACV12342
      MULTU3=MAX0(MIN0(NX,NU)-MAX0(NX-MU,0)+1,0)                        ACV12343
      IF(MULTU3.NE.0)RETURN                                             ACV12344
      RETURN 1                                                          ACV12345
      END                                                               ACV12346
      SUBROUTINE XEWU3S(INC,LAM1,MU1,LAM2,MU2,NEC,NNC,KR0A,KR0B,DEWU3P, ACV12347
     1J1TA,IAB,ICD,INDMAX,DEWU3,KR0MAX)                                 ACV12348
C     ------------------------------------------------------------------ACV12349
C     SUPPORT ROUTINE FOR XEWU3S (X PREFIX FOR 6-81 VERSION)            ACV12350
C     ------------------------------------------------------------------ACV12351
C     UPDATE/MOD: (LSU,06-81)  J.P.DRAAYER        INDEXING DEWU3        ACV12352
C                                                                       ACV12353
C     PARAMETERS--                                                      ACV12354
C       INC=0:<(LAM1,MU1)????;(LAM2,MU2)????::KR0(LAM3,MU3)HIGH>        ACV12355
C             FROM <(LAM1,MU1)????;(LAM2,MU2)HIGH::KR0(LAM3,MU3)HIGH>   ACV12356
C       INC=1:<(LAM1,MU1)????;(LAM2,MU2)????::KR0(LAM3,MU3)HIGH>        ACV12357
C             FROM <(LAM1,MU1)HIGH;(LAM2,MU2)????::KR0(LAM3,MU3)HIGH>   ACV12358
C     ------------------------------------------------------------------ACV12359
      IMPLICIT REAL*8(D)                                                ACV12360
      DIMENSION DEWU3(1),J1TA(1),DEWU3P(1),IAB(1),ICD(1)                ACV12361
      INDEX(J1TD,LAM1,J1T,J2TD,LAM2,J2T)=1+J2TD*(J2TD+1)*(3*J1TD+J2TD+5)ACV12362
     1/6+(J1TD+1)*(LAM2+J2TD-J2T)/2+(LAM1+J1TD-J1T)/2                   ACV12363
      INDQ=-KR0MAX                                                      ACV12364
      IF(INC.EQ.1)INDQ=INDQ+(INDMAX-NNC)*KR0MAX                         ACV12365
      INDPQ=-KR0MAX                                                     ACV12366
      DO 10 IIQ=1,NNC                                                   ACV12367
      INDQ=INDQ+KR0MAX                                                  ACV12368
      INDPQ=INDPQ+KR0MAX                                                ACV12369
      DO 10 KR0=KR0A,KR0B                                               ACV12370
      KI=KR0+INDQ                                                       ACV12371
      KIP=KR0+INDPQ                                                     ACV12372
   10 DEWU3(KI)=DEWU3P(KIP)                                             ACV12373
      IF(NEC.EQ.0)RETURN                                                ACV12374
      IF(INC.EQ.1)GO TO 15                                              ACV12375
      L1=LAM1                                                           ACV12376
      M1=MU1                                                            ACV12377
      L2=LAM2                                                           ACV12378
      M2=MU2                                                            ACV12379
      GO TO 20                                                          ACV12380
   15 L1=LAM2                                                           ACV12381
      M1=MU2                                                            ACV12382
      L2=LAM1                                                           ACV12383
      M2=MU1                                                            ACV12384
   20 LL1=L1+1                                                          ACV12385
      MM1=M1+1                                                          ACV12386
      LL2=L2+1                                                          ACV12387
      MM2=M2+1                                                          ACV12388
      LM1=LL1+MM1                                                       ACV12389
      LM2=LL2+MM2                                                       ACV12390
      DO 70 J2TD=1,NEC                                                  ACV12391
      J1TD=NEC-J2TD                                                     ACV12392
      J2D=J2TD-1                                                        ACV12393
      J1D=J1TD+1                                                        ACV12394
      IIQ2A=J2TD-M2                                                     ACV12395
      IF(IIQ2A.LT.0)IIQ2A=0                                             ACV12396
      IIQ2A=IIQ2A+1                                                     ACV12397
      IIQ2B=J2TD                                                        ACV12398
      IF(L2.LT.IIQ2B)IIQ2B=L2                                           ACV12399
      IIQ2B=IIQ2B+1                                                     ACV12400
      IIQ1A=J1TD-M1                                                     ACV12401
      IF(IIQ1A.LT.0)IIQ1A=0                                             ACV12402
      IIQ1A=IIQ1A+1                                                     ACV12403
      IIQ1B=J1TD                                                        ACV12404
      IF(L1.LT.IIQ1B)IIQ1B=L1                                           ACV12405
      IIQ1B=IIQ1B+1                                                     ACV12406
      DO 70 IIQ2=IIQ2A,IIQ2B                                            ACV12407
      IQ2=IIQ2-1                                                        ACV12408
      IP2=J2TD-IQ2                                                      ACV12409
      J2T=L2+IP2-IQ2                                                    ACV12410
      JJ2T=J2T+1                                                        ACV12411
      IQ=-1                                                             ACV12412
      IP=-1                                                             ACV12413
      IF(IP2.EQ.0)GO TO 25                                              ACV12414
      IQ2P=IQ2                                                          ACV12415
      IP2P=IP2-1                                                        ACV12416
      IQ2D=0                                                            ACV12417
      IP2D=1                                                            ACV12418
      IF(INC.EQ.1)IP=1                                                  ACV12419
      NM=IP2*(M2-IP2P)*(LL2+IP2)                                        ACV12420
      DN=DFLOAT(J2T)                                                    ACV12421
      GO TO 30                                                          ACV12422
   25 IQ2P=IQ2-1                                                        ACV12423
      IP2P=IP2                                                          ACV12424
      IQ2D=1                                                            ACV12425
      IP2D=0                                                            ACV12426
      IF(INC.EQ.0)IQ=1                                                  ACV12427
      NM=IQ2*(L2-IQ2P)*(LM2-IQ2)                                        ACV12428
      DN=DFLOAT(J2T+2)                                                  ACV12429
   30 J2TP=L2+IP2P-IQ2P                                                 ACV12430
      JTA=J2TD-IQ2P                                                     ACV12431
      JTB=NNC-JTA                                                       ACV12432
      NQD=NEC-IQ2P                                                      ACV12433
      DO 70 IIQ1=IIQ1A,IIQ1B                                            ACV12434
      IQ1=IIQ1-1                                                        ACV12435
      IP1=J1TD-IQ1                                                      ACV12436
      J1T=L1+IP1-IQ1                                                    ACV12437
      IF(INC.EQ.0)IND=INDEX(J1TD,L1,J1T,J2TD,L2,J2T)                    ACV12438
      IF(INC.EQ.1)IND=INDEX(J2TD,L2,J2T,J1TD,L1,J1T)                    ACV12439
      IF(J1TA(IND).LT.0)GO TO 70                                        ACV12440
      IF(IP1.EQ.M1)GO TO 50                                             ACV12441
      IQ1P=IQ1                                                          ACV12442
      IP1P=IP1+1                                                        ACV12443
      J1TP=J1T+1                                                        ACV12444
      IF(INC.EQ.0)GO TO 35                                              ACV12445
      INDP=INDEX(J2D,L2,J2TP,J1D,L1,J1TP)                               ACV12446
      IF(J1TA(INDP).LT.0)GO TO 50                                       ACV12447
      J123=JTB-IQ1                                                      ACV12448
      GO TO 40                                                          ACV12449
   35 INDP=INDEX(J1D,L1,J1TP,J2D,L2,J2TP)                               ACV12450
      IF(J1TA(INDP).LT.0)GO TO 50                                       ACV12451
      J123=JTA+IQ1                                                      ACV12452
   40 IF(IP2D.EQ.1)I=IAB(J123)                                          ACV12453
      IF(IQ2D.EQ.1)I=ICD(NQD-IQ1)                                       ACV12454
      I=JJ2T*IP1P*(MM1-IP1P)*(LL1+IP1P)*I                               ACV12455
      DC=DSQRT(DFLOAT(I)/(DFLOAT((J1T+2)*J1TP*NM)*DN))                  ACV12456
      IF(IQ.LT.0)DC=-DC                                                 ACV12457
      INDQ=(IND-1)*KR0MAX                                               ACV12458
      INDPQ=(INDP-1)*KR0MAX                                             ACV12459
      DO 45 KR0=KR0A,KR0B                                               ACV12460
      KI=KR0+INDQ                                                       ACV12461
      KIP=KR0+INDPQ                                                     ACV12462
   45 DEWU3(KI)=DC*DEWU3(KIP)                                           ACV12463
   50 IF(IQ1.EQ.L1)GO TO 70                                             ACV12464
      IQ1P=IQ1+1                                                        ACV12465
      IP1P=IP1                                                          ACV12466
      J1TP=J1T-1                                                        ACV12467
      IF(INC.EQ.0)GO TO 55                                              ACV12468
      INDP=INDEX(J2D,L2,J2TP,J1D,L1,J1TP)                               ACV12469
      IF(J1TA(INDP).LT.0)GO TO 70                                       ACV12470
      J123=JTB-IQ1                                                      ACV12471
      GO TO 60                                                          ACV12472
   55 INDP=INDEX(J1D,L1,J1TP,J2D,L2,J2TP)                               ACV12473
      IF(J1TA(INDP).LT.0)GO TO 70                                       ACV12474
      J123=JTA+IQ1                                                      ACV12475
   60 IF(IP2D.EQ.1)I=ICD(NQD-IQ1)                                       ACV12476
      IF(IQ2D.EQ.1)I=IAB(J123)                                          ACV12477
      I=JJ2T*IQ1P*(LL1-IQ1P)*(LM1-IQ1P)*I                               ACV12478
      DC=DSQRT(DFLOAT(I)/(DFLOAT((J1TP+2)*J1T*NM)*DN))                  ACV12479
      IF(IP.LT.0)DC=-DC                                                 ACV12480
      INDQ=(IND-1)*KR0MAX                                               ACV12481
      INDPQ=(INDP-1)*KR0MAX                                             ACV12482
      DO 65 KR0=KR0A,KR0B                                               ACV12483
      KI=KR0+INDQ                                                       ACV12484
      KIP=KR0+INDPQ                                                     ACV12485
   65 DEWU3(KI)=DEWU3(KI)+DC*DEWU3(KIP)                                 ACV12486
   70 CONTINUE                                                          ACV12487
      RETURN                                                            ACV12488
      END                                                               ACV12489
      SUBROUTINE XEWU3(LAM1X,MU1X,LAM2X,MU2X,LAM3X,MU3X,I3,NEC,KR0MAX,  ACV12490
     1INDMAX,DEWU3,J1TA,J2TA,IEA,N1,N2,KIMAX1)                          ACV12491
C     ------------------------------------------------------------------ACV12492
C     EXTREMAL WIGNER COEFFICIENTS FOR U3 (X PREFIX FOR 6-81 VERSION)   ACV12493
C     ------------------------------------------------------------------ACV12494
C     UPDATE/MOD: (LSU,05-80)  J.P.DRAAYER        LOG BINOMIALS         ACV12495
C                 (LSU,06-81)  J.P.DRAAYER        INDEXING DEWU3        ACV12496
C                 (LSU,03-83)  J.P.DRAAYER        SPACE SAVING MEASURE  ACV12497
C                 (LSU,02-87)  J.P.DRAAYER        OVERFLOW CORRECTION   ACV12498
C                 (LSU,10-89)  J.P.DRAAYER        ZERO OUT RELOCATED    ACV12499
C                                                                       ACV12500
C     REFERENCES--J.P.DRAAYER AND Y.AKIYAMA, J.MATH.PHYS.14(1973)1904   ACV12501
C                 K.T.HECHT, NUCL.PHYS.62(1965)1                        ACV12502
C     PARAMETERS--(I3) : (1)=GHW, (0)=GLW                               ACV12503
C       EXTERNAL--N1=MAX(KR0MAX)                                        ACV12504
C                 N2=MAX(INDMAX)=NX*(NX+1)*(NX+2)/6, NX=MAX(LAM2+MU2+1) ACV12505
C*                KIMAX1=MAX(KR0MAX*INDMAX)                             ACV12506
C       INTERNAL--X1=ABS(N1*NX)                                         ACV12507
C                 X2=ABS(NX)                                            ACV12508
C     EXTENSIONS--CHANGE EXTERNAL PARAMETERS IN CALL                    ACV12509
C                 ADJUST INTERNAL PARAMETERS BELOW                      ACV12510
C*    DIMENSIONS--DEWU3(N1*N2->KIMAX1),J1TA(N2),J2TA(N2),IEA(N2),       ACV12511
C                 DEWU3P(X1),DZ(X2),J1TAP(X2),IAB(X2),ICD(X2)           ACV12512
C       COMMENTS--ASSUME MAX N1=9,NX=42,N2=13244                        ACV12513
C                        SET X1=378,X2=42                               ACV12514
C                 DZ ARRAY ADDED FOR CORRECTING THE OVERFLOW PROBLEM    ACV12515
C     ------------------------------------------------------------------ACV12516
      IMPLICIT REAL*8(D)                                                ACV12517
      COMMON/BKDB/DBINO(8385),DBINV(8385),DTWOS(-128:128)               ACV12518
      COMMON/BKDF/DLOGF(0:2000)                                         ACV12519
      DIMENSION DEWU3(1),J1TA(1),J2TA(1),IEA(1),                        ACV12520
     1          DEWU3P(378),DZ(42),J1TAP(42),IAB(42),ICD(42)            ACV12521
      INDEX(J1TD,LAM1,J1T,J2TD,LAM2,J2T)=1+J2TD*(J2TD+1)*(3*J1TD+J2TD+5)ACV12522
     1/6+(J1TD+1)*(LAM2+J2TD-J2T)/2+(LAM1+J1TD-J1T)/2                   ACV12523
C     DIMENSION CHECKS (LSU,6-81)-START                                 ACV12524
      IF(N1.GT.9)GO TO 200                                              ACV12525
      NX=XLAM2+XMU2+1                                                   ACV12526
      IF(NX.GT.42)GO TO 210                                             ACV12527
C     DIMENSION CHECKS (LSU,6-81)-START                                 ACV12528
      KR0MAX=MULTU3(LAM1X,MU1X,LAM2X,MU2X,LAM3X,MU3X)                   ACV12529
      IF(KR0MAX.EQ.0)RETURN                                             ACV12530
      IF(I3.EQ.1)GO TO 10                                               ACV12531
      LAM1=MU1X                                                         ACV12532
      LAM2=MU2X                                                         ACV12533
      LAM3=MU3X                                                         ACV12534
      MU1=LAM1X                                                         ACV12535
      MU2=LAM2X                                                         ACV12536
      MU3=LAM3X                                                         ACV12537
      GO TO 15                                                          ACV12538
   10 LAM1=LAM1X                                                        ACV12539
      LAM2=LAM2X                                                        ACV12540
      LAM3=LAM3X                                                        ACV12541
      MU1=MU1X                                                          ACV12542
      MU2=MU2X                                                          ACV12543
      MU3=MU3X                                                          ACV12544
   15 NEC=(LAM1+LAM2-LAM3+2*(MU1+MU2-MU3))/3                            ACV12545
      IAH=(LAM2+LAM3-LAM1-NEC)/2                                        ACV12546
      IBH=(LAM3+LAM1-LAM2+NEC+2)/2                                      ACV12547
      ICH=(LAM1+LAM2-LAM3-NEC)/2                                        ACV12548
      IDH=(LAM1+LAM2+LAM3-NEC+2)/2                                      ACV12549
      DO 20 I=1,NEC                                                     ACV12550
      IAB(I)=(IAH+I)*(IBH-I)                                            ACV12551
   20 ICD(I)=(ICH+I)*(IDH+I)                                            ACV12552
      NCDMAX=MULTHY(LAM1,MU1,LAM2,MU2,LAM3,MU3)                         ACV12553
      NEC=NEC-NCDMAX                                                    ACV12554
      LAM2=LAM2-NCDMAX                                                  ACV12555
      MU2=MU2-NCDMAX                                                    ACV12556
      NCDMIN=1                                                          ACV12557
   25 IF(NCDMIN.EQ.NCDMAX)GO TO 30                                      ACV12558
      IF(MULTU3(LAM1,MU1,LAM2+1,MU2+1,LAM3,MU3).GT.0)GO TO 30           ACV12559
      NEC=NEC+1                                                         ACV12560
      LAM2=LAM2+1                                                       ACV12561
      MU2=MU2+1                                                         ACV12562
      NCDMIN=NCDMIN+1                                                   ACV12563
      GO TO 25                                                          ACV12564
C     DIMENSION MODIFICATION (LSU,6-81)-START                           ACV12565
   30 NNCMAX=NEC+NCDMAX-NCDMIN+2                                        ACV12566
      KITEST=KR0MAX*(NNCMAX)*(NNCMAX+1)*(NNCMAX+2)/6                    ACV12567
      IF(KITEST.GT.KIMAX1)GO TO 220                                     ACV12568
C     DIMENSION MODIFICATION (LSU,6-81)--STOP                           ACV12569
      DO I=1,KITEST                                                     ACV12570
         DEWU3(I)=0.D0                                                  ACV12571
      ENDDO                                                             ACV12572
      LL1=LAM1+1                                                        ACV12573
      MM1=MU1+1                                                         ACV12574
      LL2=LAM2+1                                                        ACV12575
      MM2=MU2+1                                                         ACV12576
      IA1=2*LAM1+4*MU1                                                  ACV12577
      IB1=4*LAM1+2*MU1                                                  ACV12578
      IC1=IB1-IA1                                                       ACV12579
      IA2=2*LAM2+4*MU2                                                  ACV12580
      IB2=4*LAM2+2*MU2                                                  ACV12581
      IC2=IB2-IA2                                                       ACV12582
      IS1=LL1+MM1                                                       ACV12583
      IS2=LL2+MM2                                                       ACV12584
      ISS=MM1+LAM2+MU2-NEC                                              ACV12585
      IE3=-(LAM3+2*MU3)                                                 ACV12586
      IEH=-(LAM2+2*MU2+3)                                               ACV12587
      KR0CNT=0                                                          ACV12588
      DO 135 NCD=NCDMIN,NCDMAX                                          ACV12589
      NEC=NEC+1                                                         ACV12590
      LAM2=LAM2+1                                                       ACV12591
      MU2=MU2+1                                                         ACV12592
      NNC=NEC+1                                                         ACV12593
      INDMAX=NNC*(NNC+1)*(NNC+2)/6                                      ACV12594
      IA2=IA2+6                                                         ACV12595
      IB2=IB2+6                                                         ACV12596
      IS2=IS2+2                                                         ACV12597
      ISS=ISS+1                                                         ACV12598
      IEH=IEH-3                                                         ACV12599
      LL2=LAM2+1                                                        ACV12600
      MM2=MU2+1                                                         ACV12601
      LN1=LAM1+NEC                                                      ACV12602
      LN2=LAM2+NEC                                                      ACV12603
      INN=NEC*NNC/2                                                     ACV12604
      IF(NCD.EQ.NCDMIN)GO TO 40                                         ACV12605
      DO 35 I=1,KITEST                                                  ACV12606
   35 DEWU3(I)=0.D0                                                     ACV12607
   40 DO 45 IND=1,INDMAX                                                ACV12608
      IEA(IND)=-1000                                                    ACV12609
      J2TA(IND)=-1000                                                   ACV12610
   45 J1TA(IND)=-1000                                                   ACV12611
      IE2=IEH                                                           ACV12612
      I=1000                                                            ACV12613
      DO 55 IIE=1,NNC                                                   ACV12614
      IE2=IE2+3                                                         ACV12615
      IE1=IE3-IE2                                                       ACV12616
      J2TD=IIE-1                                                        ACV12617
      J1TD=NNC-IIE                                                      ACV12618
      JJ2TA=IA2-IE2                                                     ACV12619
      JJ2TB=IB2+IE2                                                     ACV12620
      IF(JJ2TB.LT.JJ2TA)JJ2TA=JJ2TB                                     ACV12621
      JJ2TA=JJ2TA/3+1                                                   ACV12622
      JJ2TB=JJ2TA-IABS(IC2-IE2)/3                                       ACV12623
      JJ1TA=IA1-IE1                                                     ACV12624
      JJ1TB=IB1+IE1                                                     ACV12625
      IF(JJ1TB.LT.JJ1TA)JJ1TA=JJ1TB                                     ACV12626
      JJ1TA=JJ1TA/3+1                                                   ACV12627
      JJ1TB=JJ1TA-IABS(IC1-IE1)/3                                       ACV12628
      J=0                                                               ACV12629
      DO 50 JJ2T=1,JJ2TB,2                                              ACV12630
      J2T=JJ2TA-JJ2T                                                    ACV12631
      L=IABS(J2T-LAM3)                                                  ACV12632
      M=J2T+LAM3                                                        ACV12633
      DO 50 JJ1T=1,JJ1TB,2                                              ACV12634
      J1T=JJ1TA-JJ1T                                                    ACV12635
      IF(J1T.LT.L)GO TO 50                                              ACV12636
      IF(J1T.GT.M)GO TO 50                                              ACV12637
      IND=INDEX(J1TD,LAM1,J1T,J2TD,LAM2,J2T)                            ACV12638
      IEA(IND)=IE2                                                      ACV12639
      J2TA(IND)=J2T                                                     ACV12640
      J1TA(IND)=J1T                                                     ACV12641
      J=J+1                                                             ACV12642
   50 CONTINUE                                                          ACV12643
      IF(J.LT.I)I=J                                                     ACV12644
   55 CONTINUE                                                          ACV12645
      IF(I.EQ.0)GO TO 135                                               ACV12646
      IF(KR0CNT.EQ.0)GO TO 80                                           ACV12647
C     GENERATE <(LAM1,MU1)????;(LAM2,MU2)HIGH::KR0(LAM3,MU3)HIGH>       ACV12648
C     FROM <(LAM1,MU1)????;(LAM2-1,MU2-1)HIGH::KR0(LAM3,MU3)HIGH>       ACV12649
      INDQ=-KR0MAX                                                      ACV12650
      DO 75 IND=1,NNC                                                   ACV12651
      INDQ=INDQ+KR0MAX                                                  ACV12652
      J1T=J1TA(IND)                                                     ACV12653
      IF(J1T.LT.0)GO TO 75                                              ACV12654
      IQ1=(LN1-J1T)/2                                                   ACV12655
      IF(IQ1.LE.0)GO TO 65                                              ACV12656
      J1TP=J1T+1                                                        ACV12657
      INDP=(LN1-J1TP-1)/2+1                                             ACV12658
      IF(J1TAP(INDP).LT.0)GO TO 65                                      ACV12659
      I=IAB(IQ1)*IQ1*(LL1-IQ1)*(IS1-IQ1)                                ACV12660
      DC=-DSQRT(DFLOAT(I)/DFLOAT((J1T+2)*J1TP))                         ACV12661
      INDPQ=(INDP-1)*KR0MAX                                             ACV12662
      DO 60 KR0=1,KR0CNT                                                ACV12663
      KI=KR0+INDQ                                                       ACV12664
      KIP=KR0+INDPQ                                                     ACV12665
   60 DEWU3(KI)=DC*DEWU3P(KIP)                                          ACV12666
   65 IP1=NEC-IQ1                                                       ACV12667
      IF(IP1.LE.0)GO TO 75                                              ACV12668
      J1TP=J1T-1                                                        ACV12669
      INDP=(LN1-J1TP-1)/2+1                                             ACV12670
      IF(J1TAP(INDP).LT.0)GO TO 75                                      ACV12671
      I=ICD(NNC-IND)*IP1*(MM1-IP1)*(LL1+IP1)                            ACV12672
      DC=DSQRT(DFLOAT(I)/DFLOAT((J1TP+2)*J1T))                          ACV12673
      INDPQ=(INDP-1)*KR0MAX                                             ACV12674
      DO 70 KR0=1,KR0CNT                                                ACV12675
      KI=KR0+INDQ                                                       ACV12676
      KIP=KR0+INDPQ                                                     ACV12677
   70 DEWU3(KI)=DEWU3(KI)+DC*DEWU3P(KIP)                                ACV12678
   75 CONTINUE                                                          ACV12679
      INC=0                                                             ACV12680
      IF(KR0CNT.EQ.KR0MAX)GO TO 125                                     ACV12681
C     EVALUATE <(LAM1,MU1)HIGH;(LAM2,MU2)????::KR0(LAM3,MU3)HIGH>       ACV12682
C     WITH (LAM2,MU2) A MINIMUM FOR KR0=KR0CNT                          ACV12683
   80 KR0CNT=KR0CNT+1                                                   ACV12684
      I=0                                                               ACV12685
      IND=INDEX(0,LAM1,LAM1,NEC,LAM2,LN2)-1                             ACV12686
      INDQ=-KR0MAX                                                      ACV12687
      DO 85 IIQ2=1,NNC                                                  ACV12688
      INDQ=INDQ+KR0MAX                                                  ACV12689
      IND=IND+1                                                         ACV12690
      KI=KR0CNT+INDQ                                                    ACV12691
      DEWU3P(KI)=0.D0                                                   ACV12692
      IF(J1TA(IND).LT.0)GO TO 85                                        ACV12693
      I=I+1                                                             ACV12694
      IIQ2B=IIQ2                                                        ACV12695
   85 CONTINUE                                                          ACV12696
C                                                                       ACV12697
C     *****MODIFIED TO AVOID OVERFLOW (LSU,5-80)-START                  ACV12698
C      ****DLOGB CHANGED TO SAVE SPACE (LSU,3-83)****                   ACV12699
C      ****FURTHER OVERFLOW CORRECTION (LSU,2-87)****                   ACV12700
C                                                                       ACV12701
      IIQ2A=IIQ2B-I+1                                                   ACV12702
      IQ2B=IIQ2B-1                                                      ACV12703
      INDQ=(IIQ2A-2)*KR0MAX                                             ACV12704
      IZ=0                                                              ACV12705
      DO 115 IIQ2=IIQ2A,IIQ2B                                           ACV12706
      IZ=IZ+1                                                           ACV12707
      DZ(IZ)=1.D0                                                       ACV12708
      INDQ=INDQ+KR0MAX                                                  ACV12709
      L=LL2-IIQ2                                                        ACV12710
C                                                                       ACV12711
C --> START NUMERATOR PRODUCT LOOP                                      ACV12712
C                                                                       ACV12713
      IX=L-IIQ2+NNC+1                                                   ACV12714
      IF(IX.EQ.0)GO TO 120                                              ACV12715
      IY=IABS(IX)                                                       ACV12716
      IN=IX/IY                                                          ACV12717
      DN=DLOG(DFLOAT(IY))                                               ACV12718
      IF(IIQ2A.EQ.IIQ2B)GO TO 95                                        ACV12719
      DO 90 I=IIQ2A,IQ2B                                                ACV12720
      J=NNC-I                                                           ACV12721
      IF(I.LT.IIQ2)THEN                                                 ACV12722
      K=IAB(J)*ICD(J)*(IS2-I)                                           ACV12723
      ELSE                                                              ACV12724
      K=MM2-J                                                           ACV12725
      ENDIF                                                             ACV12726
      IF(K.EQ.0)GO TO 120                                               ACV12727
      IF(K.LT.0)IN=-IN                                                  ACV12728
   90 DN=DLOG(DFLOAT(IABS(K)))+DN                                       ACV12729
   95 DN=DN+DLOG(DBINO(INN+IIQ2))                                       ACV12730
C                                                                       ACV12731
C --> END NUMERATOR PRODUCT LOOP & START DENOMINATOR PRODUCT LOOP       ACV12732
C                                                                       ACV12733
      ID=1                                                              ACV12734
      DD=0.D0                                                           ACV12735
      DO 100 I=1,NNC                                                    ACV12736
      IX=I+L                                                            ACV12737
      IF(IX.LT.0)ID=-ID                                                 ACV12738
  100 DD=DLOG(DFLOAT(I+L))+DD                                           ACV12739
C                                                                       ACV12740
C --> END DENOMINATOR PRODUCT LOOP & START INNER PRODUCT/SUM LOOP       ACV12741
C                                                                       ACV12742
      IP2=NNC-IIQ2                                                      ACV12743
C                                                                       ACV12744
C     MULTIPLY BY SMALL NUMBER --> DEXP(-172) LIMIT FOR IBM SYSTEMS     ACV12745
C                                                                       ACV12746
      DZ(IZ)=DEXP(-DMIN1(DLOGF(2*IP2),172.D0))                          ACV12747
      IIP2=IP2+1                                                        ACV12748
      M=IP2*IIP2/2                                                      ACV12749
      DS=0.D0                                                           ACV12750
      DO 110 I=1,IIP2                                                   ACV12751
      DC=DZ(IZ)*DBINO(I+M)                                              ACV12752
      IF(IIP2.EQ.1)GO TO 110                                            ACV12753
      DO 105 J=1,IP2                                                    ACV12754
      IF(J.LT.I)THEN                                                    ACV12755
         K=(J+L)*(ISS+J)                                                ACV12756
      ELSE                                                              ACV12757
         K=IAB(J)                                                       ACV12758
      ENDIF                                                             ACV12759
  105 DC=DFLOAT(K)*DC                                                   ACV12760
  110 DS=DS+DC                                                          ACV12761
C                                                                       ACV12762
C --> END INNER PRODUCT/SUM LOOP & ASSIGN unnormalized DEWU3P VALUE     ACV12763
C                                                                       ACV12764
      IF(2*(IP2/2).NE.IP2)DS=-DS                                        ACV12765
      KI=KR0CNT+INDQ                                                    ACV12766
  115 DEWU3P(KI)=DFLOAT(IN*ID)*DS*DEXP((DN-DD)/2.D0)                    ACV12767
C                                                                       ACV12768
C --> START renormalization PROCEDURE                                   ACV12769
C                                                                       ACV12770
      DMIN=1.D0                                                         ACV12771
      IZ=0                                                              ACV12772
      DO 117 IIQ2=IIQ2A,IIQ2B                                           ACV12773
      IZ=IZ+1                                                           ACV12774
  117 DMIN=DMIN1(DMIN,DZ(IZ))                                           ACV12775
      INDQ=(IIQ2A-2)*KR0MAX                                             ACV12776
      IZ=0                                                              ACV12777
      DO 118 IIQ2=IIQ2A,IIQ2B                                           ACV12778
      IZ=IZ+1                                                           ACV12779
      INDQ=INDQ+KR0MAX                                                  ACV12780
      KI=KR0CNT+INDQ                                                    ACV12781
  118 DEWU3P(KI)=(DMIN/DZ(IZ))*DEWU3P(KI)                               ACV12782
C                                                                       ACV12783
C     *****MODIFIED TO AVOID OVERFLOW (LSU,5-80)--STOP                  ACV12784
C                                                                       ACV12785
      KR0A=KR0CNT                                                       ACV12786
      KR0B=KR0CNT                                                       ACV12787
C     GENERATE <(LAM1,MU1)????;(LAM2,MU2)????::KR0(LAM3,MU3)HIGH>       ACV12788
C     FROM <(LAM1,MU1)HIGH;(LAM2,MU2)????::KR0(LAM3,MU3)HIGH>           ACV12789
      CALL XEWU3S(1,LAM1,MU1,LAM2,MU2,NEC,NNC,KR0A,KR0B,DEWU3P,J1TA,    ACV12790
     1IAB,ICD,INDMAX,DEWU3,KR0MAX)                                      ACV12791
      INC=1                                                             ACV12792
      GO TO 125                                                         ACV12793
  120 KR0CNT=KR0CNT-1                                                   ACV12794
      WRITE(1,195)XLAM1,XMU1,XLAM2,XMU2,XLAM3,XMU3,KR0MAX,KR0CNT        ACV12795
  125 IF(KR0CNT.EQ.0)GO TO 135                                          ACV12796
      INDQ=-KR0MAX                                                      ACV12797
      DO 130 IND=1,NNC                                                  ACV12798
      INDQ=INDQ+KR0MAX                                                  ACV12799
      J1TAP(IND)=J1TA(IND)                                              ACV12800
      DO 130 KR0=1,KR0CNT                                               ACV12801
      KI=KR0+INDQ                                                       ACV12802
  130 DEWU3P(KI)=DEWU3(KI)                                              ACV12803
  135 CONTINUE                                                          ACV12804
      IF(KR0CNT.EQ.0)RETURN                                             ACV12805
      KR0A=1                                                            ACV12806
      KR0B=KR0CNT-INC                                                   ACV12807
      IF(KR0B.EQ.0)GO TO 140                                            ACV12808
C     GENERATE <(LAM1,MU1)????;(LAM2,MU2)????::KR0(LAM3,MU3)HIGH>       ACV12809
C     FROM <(LAM1,MU1)????;(LAM2,MU2)HIGH::KR0(LAM3,MU3)HIGH>           ACV12810
      CALL XEWU3S(0,LAM1,MU1,LAM2,MU2,NEC,NNC,KR0A,KR0B,DEWU3P,J1TA,    ACV12811
     1IAB,ICD,INDMAX,DEWU3,KR0MAX)                                      ACV12812
  140 CONTINUE                                                          ACV12813
C     RENORMALIZE VIA LARGEST ELEMENT TO AVOID OVERFLOW (LSU,5-80)-STARTACV12814
      DO 142 KR0=1,KR0CNT                                               ACV12815
      DC=1.D0                                                           ACV12816
      INDQ=-KR0MAX                                                      ACV12817
      DO 141 IND=1,INDMAX                                               ACV12818
      INDQ=INDQ+KR0MAX                                                  ACV12819
      KI=KR0+INDQ                                                       ACV12820
  141 DC=DMAX1(DC,DABS(DEWU3(KI)))                                      ACV12821
      INDQ=-KR0MAX                                                      ACV12822
      DO 142 IND=1,INDMAX                                               ACV12823
      INDQ=INDQ+KR0MAX                                                  ACV12824
      KI=KR0+INDQ                                                       ACV12825
  142 DEWU3(KI)=DEWU3(KI)/DC                                            ACV12826
C     RENORMALIZE VIA LARGEST ELEMENT TO AVOID OVERFLOW (LSU,5-80)--STOPACV12827
C     ORTHONORMALIZATION OF SOLUTIONS                                   ACV12828
      DO 165 KR0=1,KR0CNT                                               ACV12829
      KR0PA=1                                                           ACV12830
      IF(INC.EQ.1.AND.KR0.EQ.KR0CNT)KR0PA=KR0CNT                        ACV12831
      DO 165 KR0P=KR0PA,KR0                                             ACV12832
      DN=0.D0                                                           ACV12833
      INDQ=-KR0MAX                                                      ACV12834
      DO 145 IND=1,INDMAX                                               ACV12835
      INDQ=INDQ+KR0MAX                                                  ACV12836
      KI=KR0+INDQ                                                       ACV12837
      KIP=KR0P+INDQ                                                     ACV12838
  145 DN=DN+DEWU3(KI)*DEWU3(KIP)                                        ACV12839
      IF(KR0P.EQ.KR0)GO TO 155                                          ACV12840
      IF(DABS(DN).LT.1.D-12)GO TO 165                                   ACV12841
      INDQ=-KR0MAX                                                      ACV12842
      DO 150 IND=1,INDMAX                                               ACV12843
      INDQ=INDQ+KR0MAX                                                  ACV12844
      KI=KR0+INDQ                                                       ACV12845
      KIP=KR0P+INDQ                                                     ACV12846
  150 DEWU3(KI)=DEWU3(KI)-DN*DEWU3(KIP)                                 ACV12847
      GO TO 165                                                         ACV12848
  155 DN=1.D0/DSQRT(DN)                                                 ACV12849
      INDQ=-KR0MAX                                                      ACV12850
      DO 160 IND=1,INDMAX                                               ACV12851
      INDQ=INDQ+KR0MAX                                                  ACV12852
      KI=KR0+INDQ                                                       ACV12853
  160 DEWU3(KI)=DN*DEWU3(KI)                                            ACV12854
  165 CONTINUE                                                          ACV12855
C     SET PHASE CONVENTION (K.T.HECHT, NUCL.PHYS.62(1965)1)             ACV12856
      IE2=IE3+(LAM1+2*MU1)                                              ACV12857
      IPH=2*(LAM1+LAM2-LAM3+MU1+MU2-MU3+KR0MAX)                         ACV12858
      INDQ=-KR0MAX                                                      ACV12859
      DO 180 IND=1,INDMAX                                               ACV12860
      INDQ=INDQ+KR0MAX                                                  ACV12861
      IF(IEA(IND).NE.IE2)GO TO 180                                      ACV12862
      I=IPH+J1TA(IND)+J2TA(IND)-LAM3                                    ACV12863
      DO 175 KR0=1,KR0MAX                                               ACV12864
      I=I-2                                                             ACV12865
      KI=KR0+INDQ                                                       ACV12866
      J=I                                                               ACV12867
      IF(DEWU3(KI).LT.0.D0)J=J-2                                        ACV12868
      IF(4*(J/4).EQ.J)GO TO 175                                         ACV12869
      INDPQ=-KR0MAX                                                     ACV12870
      DO 170 INDP=1,INDMAX                                              ACV12871
      INDPQ=INDPQ+KR0MAX                                                ACV12872
      KIP=KR0+INDPQ                                                     ACV12873
  170 DEWU3(KIP)=-DEWU3(KIP)                                            ACV12874
  175 CONTINUE                                                          ACV12875
      GO TO 185                                                         ACV12876
  180 CONTINUE                                                          ACV12877
  185 IF(I3.EQ.1)RETURN                                                 ACV12878
      INDQ=-KR0MAX                                                      ACV12879
      DO 190 IND=1,INDMAX                                               ACV12880
      INDQ=INDQ+KR0MAX                                                  ACV12881
      IEA(IND)=-IEA(IND)                                                ACV12882
      I=IPH+J1TA(IND)+J2TA(IND)-LAM3                                    ACV12883
      DO 190 KR0=1,KR0MAX                                               ACV12884
      I=I-2                                                             ACV12885
      KI=KR0+INDQ                                                       ACV12886
      IF(4*(I/4).NE.I)DEWU3(KI)=-DEWU3(KI)                              ACV12887
  190 CONTINUE                                                          ACV12888
  195 FORMAT(28H *****U3 COUPLING ERROR*****,3X,3(4X,2I3),3X,           ACV12889
     115HRO(ABSOLUTE) = ,I2,3X,15HRO(RELATIVE) = ,I2,4X,                ACV12890
     226H*****REPORT TO AUTHOR*****)                                    ACV12891
      RETURN                                                            ACV12892
C     DIMENSION CHECKS (LSU,6-81)-START                                 ACV12893
  200 WRITE(6,205)N1                                                    ACV12894
  205 FORMAT(36H ***** XEWU3 DIMENSION OVERFLOW: N1=,I10)               ACV12895
      GO TO 230                                                         ACV12896
  210 WRITE(6,215)NX                                                    ACV12897
  215 FORMAT(36H ***** XEWU3 DIMENSION OVERFLOW: NX=,I10)               ACV12898
      GO TO 230                                                         ACV12899
  220 WRITE(6,225)KITEST,KIMAX1                                         ACV12900
  225 FORMAT(40H ***** XEWU3 DIMENSION OVERFLOW: KITEST=,I10,5X,        ACV12901
     17HKIMAX1=,I10)                                                    ACV12902
  230 STOP                                                              ACV12903
C     DIMENSION CHECKS (LSU,6-81)--STOP                                 ACV12904
      END                                                               ACV12905
      SUBROUTINE XWU3(LAM1,MU1,LAM2,MU2,LAM3,MU3,IE,JT,NEC,DEWU3,       ACV12906
     1KR0MAX,INDMAX,DWU3,J1SMAX,J1TMAX,J2SMAX,J2TMAX,IESMAX,IE2MAX,     ACV12907
     2INDMAT,N1,N2,N3,KIMAX2)                                           ACV12908
C     ------------------------------------------------------------------ACV12909
C     WIGNER COEFFICIENTS FOR U3 (X PREFIX FOR 6-81 VERSION)            ACV12910
C     ------------------------------------------------------------------ACV12911
C     UPDATE/MOD: (LSU,06-81)  J.P.DRAAYER        INDEXING OF DEWU3     ACV12912
C                 (LSU,11-89)  J.P.DRAAYER        DWU3 ZERO-OUT RANGE   ACV12913
C                                                                       ACV12914
C     REFERENCES--J.P.DRAAYER AND Y.AKIYAMA, J.MATH.PHYS.14(1973)1904   ACV12915
C                 K.T.HECHT, NUCL.PHYS.62(1965)1                        ACV12916
C     PARAMETERS--SEE ALSO XEWU3                                        ACV12917
C       EXTERNAL--N1=MAX(KR0MAX)                                        ACV12918
C                 N2=MAX(IESMAX) SAFE TO SET N2=MAX(LAM2+MU2+1)         ACV12919
C                 N3=MAX(DIM(LAM2,MU2))                                 ACV12920
C                 NA=MAX(INDMAX)=NX*(NX+1)*(NX+2)/6, NX=MAX(LAM2+MU2+1) ACV12921
C                 NB=MAX(J2SMAX) SAFE TO SET NB=MAX(LAM2+MU2+1)         ACV12922
C*                KIMAX1=MAX(KR0MAX*INDMAX) (SEE XEWU3)                 ACV12923
C*                KIMAX2=MAX(KR0MAX*DIM(LAM2,MU2))                      ACV12924
C       INTERNAL--X1=ABS(N1*N3)                                         ACV12925
C                 X2=ABS(N2)                                            ACV12926
C                 X3=ABS(N2*NB)                                         ACV12927
C     EXTENSIONS--CHANGE EXTERNAL PARAMETERS IN CALL STATEMENT          ACV12928
C                 ADJUST INTERNAL PARAMETERS BELOW                      ACV12929
C*    DIMENSIONS--DEWU3(N1*NA),DWU3(N1*N3),J1SMAX(N2*NB),               ACV12930
C                 J1TMAX(N2*NB),J2SMAX(N2),J2TMAX(N2),INDMAT(N2*NB),    ACV12931
C                 DWU3P(X1),J2TMAP(X2),INDMAP(X3)                       ACV12932
C*      COMMENTS--USE N1*NA->KIMAX1,N1*N3->KIMAX2                       ACV12933
C                 ASSUME MAX N1=9,N2=42,N3=9030                         ACV12934
C                        SET X1=27090,X2=42,X3=1764 (X1=3*N3,FIXED)     ACV12935
C     ------------------------------------------------------------------ACV12936
      IMPLICIT REAL*8(D),INTEGER(X)                                     ACV12937
C     IMPLICIT INTEGER(X)                                               ACV12938
      DIMENSION DEWU3(1),DWU3(1),J1SMAX(1),                             ACV12939
     1          J1TMAX(1),J2SMAX(1),J2TMAX(1),INDMAT(1),                ACV12940
     2          DWU3P(27090),J2TMAP(42),INDMAP(1764)                    ACV12941
      INDEX(J1TD,LAM1,J1T,J2TD,LAM2,J2T)=1+J2TD*(J2TD+1)*(3*J1TD+J2TD+5)ACV12942
     1/6+(J1TD+1)*(LAM2+J2TD-J2T)/2+(LAM1+J1TD-J1T)/2                   ACV12943
      IDM(LAM,MU)=(LAM+1)*(MU+1)*(LAM+MU+2)/2                           ACV12944
C     DIMENSION CHECKS (LSU,6-81)-START                                 ACV12945
      IF(N1.GT.9)GO TO 130                                              ACV12946
      IF(N2.GT.42)GO TO 140                                             ACV12947
      IF(N3.GT.9030)GO TO 150                                           ACV12948
      IDTEST=KR0MAX*IDM(LAM2,MU2)                                       ACV12949
      IF(IDTEST.GT.KIMAX2.OR.IDTEST.GT.27090)GO TO 160                  ACV12950
C     DIMENSION CHECKS (LSU,6-81)--STOP                                 ACV12951
      LL1=LAM1+1                                                        ACV12952
      MM1=MU1+1                                                         ACV12953
      LL2=LAM2+1                                                        ACV12954
      MM2=MU2+1                                                         ACV12955
      LL3=LAM3+1                                                        ACV12956
      MM3=MU3+1                                                         ACV12957
      LM1=LAM1+MU1                                                      ACV12958
      LM2=LAM2+MU2                                                      ACV12959
      LLMM1=LL1+MM1                                                     ACV12960
      LLMM2=LL2+MM2                                                     ACV12961
      LLMM3=LL3+MM3                                                     ACV12962
      JJTD=(IE+LAM3+2*MU3)/3+1                                          ACV12963
      IP=(JJTD+JT-LL3)/2                                                ACV12964
      NCC=NEC-1                                                         ACV12965
      INC=1                                                             ACV12966
      IQ3=0                                                             ACV12967
      IP3=-1                                                            ACV12968
      J3T=LAM3+IP3-IQ3                                                  ACV12969
      DO 125 JJ3TD=1,JJTD                                               ACV12970
C     DO 10 N=1,KIMAX2                                                  ACV12971
      DO 10 N=1,IDTEST                                                  ACV12972
   10 DWU3(N)=0.D0                                                      ACV12973
      NCC=NCC+1                                                         ACV12974
      IF(IP3.EQ.IP)INC=0                                                ACV12975
      IF(INC.EQ.1)GO TO 15                                              ACV12976
      IQ3=IQ3+1                                                         ACV12977
      J3T=J3T-1                                                         ACV12978
      NM=(LL3-IQ3)*IQ3*(LLMM3-IQ3)                                      ACV12979
      GO TO 20                                                          ACV12980
   15 IP3=IP3+1                                                         ACV12981
      J3T=J3T+1                                                         ACV12982
      NM=(MM3-IP3)*IP3*(LL3+IP3)                                        ACV12983
   20 JJ2TDA=NCC-LM1                                                    ACV12984
      IF(JJ2TDA.LT.0)JJ2TDA=0                                           ACV12985
      JJ2TDA=JJ2TDA+1                                                   ACV12986
      JJ2TDB=LM2                                                        ACV12987
      IF(NCC.LT.JJ2TDB)JJ2TDB=NCC                                       ACV12988
      JJ2TDB=JJ2TDB+1                                                   ACV12989
      JJ2TDC=JJ2TDA                                                     ACV12990
      IND=0                                                             ACV12991
      IES=0                                                             ACV12992
      DO 115 JJ2TD=JJ2TDA,JJ2TDB                                        ACV12993
      J2TD=JJ2TD-1                                                      ACV12994
      J1TD=NCC-J2TD                                                     ACV12995
      IES=IES+1                                                         ACV12996
      IIQ2A=J2TD-MU2                                                    ACV12997
      IF(IIQ2A.LT.0)IIQ2A=0                                             ACV12998
      IIQ2A=IIQ2A+1                                                     ACV12999
      IIQ2B=J2TD                                                        ACV13000
      IF(LAM2.LT.IIQ2B)IIQ2B=LAM2                                       ACV13001
      IIQ2B=IIQ2B+1                                                     ACV13002
      IIQ1A=J1TD-MU1                                                    ACV13003
      IF(IIQ1A.LT.0)IIQ1A=0                                             ACV13004
      IIQ1A=IIQ1A+1                                                     ACV13005
      IIQ1B=J1TD                                                        ACV13006
      IF(LAM1.LT.IIQ1B)IIQ1B=LAM1                                       ACV13007
      IIQ1B=IIQ1B+1                                                     ACV13008
      J2S=0                                                             ACV13009
      DO 105 IIQ2=IIQ2A,IIQ2B                                           ACV13010
      IQ2=IIQ2-1                                                        ACV13011
      IP2=J2TD-IQ2                                                      ACV13012
      J2T=LAM2+IP2-IQ2                                                  ACV13013
      J23S=J2T+J3T                                                      ACV13014
      J23D=J3T-J2T                                                      ACV13015
      J23H=IABS(J23D)                                                   ACV13016
      J1S=0                                                             ACV13017
      DO 100 IIQ1=IIQ1A,IIQ1B                                           ACV13018
      IQ1=IIQ1-1                                                        ACV13019
      IP1=J1TD-IQ1                                                      ACV13020
      J1T=LAM1+IP1-IQ1                                                  ACV13021
      IF(J1T.LT.J23H)GO TO 100                                          ACV13022
      IF(J1T.GT.J23S)GO TO 100                                          ACV13023
      J1TS=J1T                                                          ACV13024
      J2TS=J2T                                                          ACV13025
      INDQ=IND*KR0MAX                                                   ACV13026
      IND=IND+1                                                         ACV13027
      J1S=J1S+1                                                         ACV13028
      IF(JJ3TD.EQ.1)GO TO 90                                            ACV13029
      JA=(J23S-J1T)/2                                                   ACV13030
      JJA=JA+1                                                          ACV13031
      JB=(J23D+J1T)/2                                                   ACV13032
      JJB=JB+1                                                          ACV13033
      JC=(J1T+J23S)/2+1                                                 ACV13034
      JJC=JC+1                                                          ACV13035
      JD=(J1T-J23D)/2+1                                                 ACV13036
      JJD=JD-1                                                          ACV13037
      IESP=J2TD-JJ2TDP                                                  ACV13038
      DO 85 I=1,4                                                       ACV13039
      IF(I.EQ.1)IESP=IESP+1                                             ACV13040
      IF(I.EQ.3)IESP=IESP+1                                             ACV13041
      IF(IESP.LT.1)GO TO 85                                             ACV13042
      IF(IESP.GT.IESMAX)GO TO 85                                        ACV13043
      GO TO (25,35,45,55),I                                             ACV13044
   25 J2TP=J2T+1                                                        ACV13045
      J1TP=J1T                                                          ACV13046
      IF(J1TP.LT.IABS(J2TP-J3TP))GO TO 85                               ACV13047
      IF(J1TP.GT.J2TP+J3TP)GO TO 85                                     ACV13048
      M=IQ2                                                             ACV13049
      IF(M.EQ.0)GO TO 85                                                ACV13050
      N=LL2-M                                                           ACV13051
      N=(LLMM2-M)*N                                                     ACV13052
      J12TP=J2T+1                                                       ACV13053
      IF(INC.EQ.1)GO TO 30                                              ACV13054
      IAB=JJA                                                           ACV13055
      ICD=JJC                                                           ACV13056
      IPH=1                                                             ACV13057
      GO TO 65                                                          ACV13058
   30 IAB=JB                                                            ACV13059
      ICD=JD                                                            ACV13060
      IPH=-1                                                            ACV13061
      GO TO 65                                                          ACV13062
   35 J2TP=J2T-1                                                        ACV13063
      J1TP=J1T                                                          ACV13064
      IF(J1TP.LT.IABS(J2TP-J3TP))GO TO 85                               ACV13065
      IF(J1TP.GT.J2TP+J3TP)GO TO 85                                     ACV13066
      M=IP2                                                             ACV13067
      IF(M.EQ.0)GO TO 85                                                ACV13068
      N=MM2-M                                                           ACV13069
      N=(LLMM2-N)*N                                                     ACV13070
      J12TP=J2T                                                         ACV13071
      IF(INC.EQ.1)GO TO 40                                              ACV13072
      IAB=JJB                                                           ACV13073
      ICD=JJD                                                           ACV13074
      IPH=1                                                             ACV13075
      GO TO 65                                                          ACV13076
   40 IAB=JA                                                            ACV13077
      ICD=JC                                                            ACV13078
      IPH=1                                                             ACV13079
      GO TO 65                                                          ACV13080
   45 J2TP=J2T                                                          ACV13081
      J1TP=J1T+1                                                        ACV13082
      IF(J1TP.LT.IABS(J2TP-J3TP))GO TO 85                               ACV13083
      IF(J1TP.GT.J2TP+J3TP)GO TO 85                                     ACV13084
      M=IQ1                                                             ACV13085
      IF(M.EQ.0)GO TO 85                                                ACV13086
      N=LL1-M                                                           ACV13087
      N=(LLMM1-M)*N                                                     ACV13088
      J12TP=J1T+1                                                       ACV13089
      IF(INC.EQ.1)GO TO 50                                              ACV13090
      IAB=JJB                                                           ACV13091
      ICD=JJC                                                           ACV13092
      IPH=1                                                             ACV13093
      GO TO 65                                                          ACV13094
   50 IAB=JA                                                            ACV13095
      ICD=JD                                                            ACV13096
      IPH=1                                                             ACV13097
      GO TO 65                                                          ACV13098
   55 J2TP=J2T                                                          ACV13099
      J1TP=J1T-1                                                        ACV13100
      IF(J1TP.LT.IABS(J2TP-J3TP))GO TO 85                               ACV13101
      IF(J1TP.GT.J2TP+J3TP)GO TO 85                                     ACV13102
      M=IP1                                                             ACV13103
      IF(M.EQ.0)GO TO 85                                                ACV13104
      N=MM1-M                                                           ACV13105
      N=(LLMM1-N)*N                                                     ACV13106
      J12TP=J1T                                                         ACV13107
      IF(INC.EQ.1)GO TO 60                                              ACV13108
      IAB=JJA                                                           ACV13109
      ICD=JJD                                                           ACV13110
      IPH=-1                                                            ACV13111
      GO TO 65                                                          ACV13112
   60 IAB=JB                                                            ACV13113
      ICD=JC                                                            ACV13114
      IPH=1                                                             ACV13115
   65 IF(J12TP.GT.0)GO TO 70                                            ACV13116
      IF(INC.EQ.1)IAB=1                                                 ACV13117
      IF(INC.EQ.0)ICD=1                                                 ACV13118
      DC=1.D0                                                           ACV13119
      GO TO 75                                                          ACV13120
   70 DC=DFLOAT(J12TP*(J12TP+1))                                        ACV13121
   75 J2SP=(J2TMAP(IESP)-J2TP+2)/2                                      ACV13122
      INDP=(INDMAP(IESP+(J2SP-1)*N2)-J1TP)/2                            ACV13123
      DC=DSQRT(DFLOAT(IAB*ICD*M*N)/(DFLOAT(NM)*DC))                     ACV13124
      IF(IPH.LT.0)DC=-DC                                                ACV13125
      INDPQ=(INDP-1)*KR0MAX                                             ACV13126
      DO 80 KR0=1,KR0MAX                                                ACV13127
      KI=KR0+INDQ                                                       ACV13128
      KIP=KR0+INDPQ                                                     ACV13129
   80 DWU3(KI)=DWU3(KI)+DC*DWU3P(KIP)                                   ACV13130
   85 CONTINUE                                                          ACV13131
      GO TO 100                                                         ACV13132
   90 INDPQ=(INDEX(J1TD,LAM1,J1T,J2TD,LAM2,J2T)-1)*KR0MAX               ACV13133
      DO 95 KR0=1,KR0MAX                                                ACV13134
      KI=KR0+INDQ                                                       ACV13135
      KIP=KR0+INDPQ                                                     ACV13136
   95 DWU3(KI)=DEWU3(KIP)                                               ACV13137
  100 CONTINUE                                                          ACV13138
      IF(J1S.EQ.0)GO TO 105                                             ACV13139
      IESJ2S=IES+J2S*N2                                                 ACV13140
      J2S=J2S+1                                                         ACV13141
      J1SMAX(IESJ2S)=J1S                                                ACV13142
      J1TMAX(IESJ2S)=J1TS+2*(J1S-1)                                     ACV13143
      INDMAT(IESJ2S)=2*IND+J1TS                                         ACV13144
  105 CONTINUE                                                          ACV13145
      IF(J2S.NE.0)GO TO 110                                             ACV13146
      IES=IES-1                                                         ACV13147
      IF(IES.EQ.0)JJ2TDC=JJ2TDC+1                                       ACV13148
      GO TO 115                                                         ACV13149
  110 J2SMAX(IES)=J2S                                                   ACV13150
      J2TMAX(IES)=J2TS+2*(J2S-1)                                        ACV13151
  115 CONTINUE                                                          ACV13152
      IESMAX=IES                                                        ACV13153
      IF(JJ3TD.EQ.JJTD)GO TO 125                                        ACV13154
      J3TP=J3T                                                          ACV13155
      JJ2TDP=JJ2TDC                                                     ACV13156
      IND=0                                                             ACV13157
      DO 120 IES=1,IESMAX                                               ACV13158
      J2TMAP(IES)=J2TMAX(IES)                                           ACV13159
      J2SB=J2SMAX(IES)                                                  ACV13160
      J2SQ=-N2                                                          ACV13161
      DO 120 J2S=1,J2SB                                                 ACV13162
      J2SQ=J2SQ+N2                                                      ACV13163
      IESJ2S=IES+J2SQ                                                   ACV13164
      INDMAP(IESJ2S)=INDMAT(IESJ2S)                                     ACV13165
      J1SB=J1SMAX(IESJ2S)                                               ACV13166
      DO 120 J1S=1,J1SB                                                 ACV13167
      INDQ=IND*KR0MAX                                                   ACV13168
      IND=IND+1                                                         ACV13169
      DO 120 KR0=1,KR0MAX                                               ACV13170
      KI=KR0+INDQ                                                       ACV13171
  120 DWU3P(KI)=DWU3(KI)                                                ACV13172
  125 CONTINUE                                                          ACV13173
      INDMAX=IND                                                        ACV13174
      IE2MAX=-(LAM2+2*MU2)+3*(JJ2TDC-1)+3*(IES-1)                       ACV13175
      RETURN                                                            ACV13176
C     DIMENSION CHECKS (LSU,6-81)-START                                 ACV13177
  130 WRITE(6,135)N1                                                    ACV13178
  135 FORMAT(35H *****XWU3 DIMENSION OVERFLOW:  N1=,I10)                ACV13179
      GO TO 170                                                         ACV13180
  140 WRITE(6,145)N2                                                    ACV13181
  145 FORMAT(35H *****XWU3 DIMENSION OVERFLOW:  N2=,I10)                ACV13182
      GO TO 170                                                         ACV13183
  150 WRITE(6,155)N3                                                    ACV13184
  155 FORMAT(35H *****XWU3 DIMENSION OVERFLOW:  N3=,I10)                ACV13185
      GO TO 170                                                         ACV13186
  160 WRITE(6,165)IDTEST                                                ACV13187
  165 FORMAT(39H *****XWU3 DIMENSION OVERFLOW:  IDTEST=,I10)            ACV13188
  170 STOP                                                              ACV13189
C     DIMENSION CHECKS (LSU,6-81)--STOP                                 ACV13190
      END                                                               ACV13191
      FUNCTION MULTHY(L1,M1,L2,M2,L3,M3)                                ACV13192
C     ------------------------------------------------------------------ACV13193
C     MULTIPLICITY (THEORY) IN U3 COUPLING                              ACV13194
C     ------------------------------------------------------------------ACV13195
      DIMENSION IX(6)                                                   ACV13196
      MULTHY=0                                                          ACV13197
      IX(1)=L1+L2-L3+2*(M1+M2-M3)                                       ACV13198
      IX(2)=M1+M2-M3+2*(L1+L2-L3)                                       ACV13199
      IX(3)=2*L2+M2+M1-L1-M3+L3                                         ACV13200
      IX(4)=2*M2+L2+L1-M1-L3+M3                                         ACV13201
      IX(5)=L3+M2-L1+2*(M3+L2-M1)                                       ACV13202
      IX(6)=M3+L2-M1+2*(L3+M2-L1)                                       ACV13203
      IXMIN=1000                                                        ACV13204
      DO 10 I=1,6                                                       ACV13205
      IXDB3=IX(I)/3                                                     ACV13206
      IF(3*IXDB3.LT.IX(I))RETURN                                        ACV13207
      IF(IXDB3.LT.IXMIN)IXMIN=IXDB3                                     ACV13208
   10 CONTINUE                                                          ACV13209
      IF(IXMIN.LT.0)RETURN                                              ACV13210
      MULTHY=MIN0(IXMIN,L2,M2)+1                                        ACV13211
      RETURN                                                            ACV13212
      END                                                               ACV13213
C     ------------------------------------------------------------------ACV13214
C                                                                       ACV13215
C     AUTHOR: ORIGINAL CODE BY J. P. DRAAYER (U. OF MICHIGAN, 1970-1974)ACV13216
C                                                                       ACV13217
C     UPDATE: 02/25/79 --> MODIFIED MTS R3PACK (LOG FACTORIALS INSERTED)ACV13218
C             01/10/88 --> MODIFIED FOR ENHANCED APPLICATIONS (V/P WORK)ACV13219
C                          1) DLOGF REDEFINED & RANGE EXTENDED IN BLOCKSACV13220
C                          2) DLOGF IMPLEMENTED TO AVOID XFLOWS IN DELTAACV13221
C                          3) BTEST INSERTED TO SIMPLIFY & IMPROVE CODESACV13222
C                                                                       ACV13223
C     NOTICE: 1) BLKNEW MUST BE CALLED (ONCE) BEFORE USING THE PROGRAMS.ACV13224
C             2) THE RANGE OF DLOGF CAN BE EXTENDED FOR HIGHER J VALUES.ACV13225
C             3) A SIMILAR PACKAGE FOR VECTOR APPLICATIONS IS AVAILABLE.ACV13226
C                                                                       ACV13227
C     ------------------------------------------------------------------ACV13228
      FUNCTION DWR3(J1T,J2T,J3T,M1T,M2T,M3T)                            ACV13229
C     ------------------------------------------------------------------ACV13230
C     WIGNER COEFFICIENTS FOR R3--TRIANGLE RELATIONS CHECKED IN DELTA   ACV13231
C     REFERENCES--ELEMENTARY THEORY OF ANGULAR MOMENTUM, M.E.ROSE, WILEYACV13232
C     ------------------------------------------------------------------ACV13233
      IMPLICIT REAL*8(D)                                                ACV13234
      COMMON/BKDF/DLOGF(0:2000)                                         ACV13235
      DWR3=0.D0                                                         ACV13236
      IF(M1T+M2T-M3T.NE.0)GOTO 20                                       ACV13237
      DC=DELTA(J1T,J2T,J3T)                                             ACV13238
      IF(DC.EQ.12345D0)GOTO 20                                          ACV13239
      I1=J3T-J2T+M1T                                                    ACV13240
      I2=J3T-J1T-M2T                                                    ACV13241
      I3=J1T+J2T-J3T                                                    ACV13242
      I4=J1T-M1T                                                        ACV13243
      IF(BTEST(I4,0))GOTO 20                                            ACV13244
      I5=J2T+M2T                                                        ACV13245
      IF(BTEST(I5,0))GOTO 20                                            ACV13246
      ITMIN=MAX0(0,-I1,-I2)                                             ACV13247
      ITMAX=MIN0(I3,I4,I5)                                              ACV13248
      IF(ITMIN.GT.ITMAX)GOTO 20                                         ACV13249
      DTOP=(DLOG(DFLOAT(J3T+1))+DC+DLOGF(J1T+M1T)+DLOGF(J1T-M1T)+       ACV13250
     1DLOGF(J2T+M2T)+DLOGF(J2T-M2T)+DLOGF(J3T+M3T)+                     ACV13251
     2DLOGF(J3T-M3T))/DFLOAT(2)                                         ACV13252
      DO 10 IT=ITMIN,ITMAX,2                                            ACV13253
      DBOT=DLOGF(I3-IT)+DLOGF(I4-IT)+DLOGF(I5-IT)+                      ACV13254
     1DLOGF(IT)+DLOGF(I1+IT)+DLOGF(I2+IT)                               ACV13255
      DSUM=DEXP(DTOP-DBOT)                                              ACV13256
      IF(BTEST(IT,1))THEN                                               ACV13257
      DWR3=DWR3-DSUM                                                    ACV13258
      ELSE                                                              ACV13259
      DWR3=DWR3+DSUM                                                    ACV13260
      ENDIF                                                             ACV13261
   10 CONTINUE                                                          ACV13262
   20 RETURN                                                            ACV13263
      END                                                               ACV13264
      FUNCTION DRR3(J1T,J2T,L2T,L1T,J3T,L3T)                            ACV13265
C     ------------------------------------------------------------------ACV13266
C     RACAH COEFFICIENTS FOR R3--TRIANGLE RELATION CHECKED IN DELTA     ACV13267
C     REFERENCES--THE 3-J AND 6-J SYMBOLS, M.ROTENBERG, R.BIVINS,       ACV13268
C                 N.METROPOLIS AND J.K.WOOTEN, MIT PRESS                ACV13269
C     ------------------------------------------------------------------ACV13270
      IMPLICIT REAL*8(D)                                                ACV13271
      COMMON/BKDF/DLOGF(0:2000)                                         ACV13272
      DRR3=0.D0                                                         ACV13273
      DX=DELTA(J1T,J2T,J3T)                                             ACV13274
      IF(DX.EQ.12345D0)GOTO 20                                          ACV13275
      DC=DX                                                             ACV13276
      DX=DELTA(L1T,L2T,J3T)                                             ACV13277
      IF(DX.EQ.12345D0)GOTO 20                                          ACV13278
      DC=DX+DC                                                          ACV13279
      DX=DELTA(L1T,J2T,L3T)                                             ACV13280
      IF(DX.EQ.12345D0)GOTO 20                                          ACV13281
      DC=DX+DC                                                          ACV13282
      DX=DELTA(J1T,L2T,L3T)                                             ACV13283
      IF(DX.EQ.12345D0)GOTO 20                                          ACV13284
      DC=(DX+DC)/2.D0                                                   ACV13285
      I1=J3T+L3T-J1T-L1T                                                ACV13286
      I2=J3T+L3T-J2T-L2T                                                ACV13287
      I3=J1T+J2T+L1T+L2T+2                                              ACV13288
      I4=J1T+J2T-J3T                                                    ACV13289
      I5=L1T+L2T-J3T                                                    ACV13290
      I6=J1T+L2T-L3T                                                    ACV13291
      I7=L1T+J2T-L3T                                                    ACV13292
      ITMIN=MAX0(0,-I1,-I2)                                             ACV13293
      ITMAX=MIN0(I3,I4,I5,I6,I7)                                        ACV13294
      IF(ITMIN.GT.ITMAX)GOTO 20                                         ACV13295
      DO 10 IT=ITMIN,ITMAX,2                                            ACV13296
      DSUM=DEXP(DC+DLOGF(I3-IT)-(DLOGF(I4-IT)+DLOGF(I5-IT)+             ACV13297
     1DLOGF(I6-IT)+DLOGF(I7-IT)+DLOGF(IT)+DLOGF(I1+IT)+DLOGF(I2+IT)))   ACV13298
      IF(BTEST(IT,1))THEN                                               ACV13299
      DRR3=DRR3-DSUM                                                    ACV13300
      ELSE                                                              ACV13301
      DRR3=DRR3+DSUM                                                    ACV13302
      ENDIF                                                             ACV13303
   10 CONTINUE                                                          ACV13304
   20 RETURN                                                            ACV13305
      END                                                               ACV13306
      FUNCTION DJHR3(J1T,J2T,J3T,J4T,J5T,J6T,J7T,J8T,J9T)               ACV13307
C     ------------------------------------------------------------------ACV13308
C     JAHN-HOPE COEFFICIENTS FOR R3--TRIANGLE RELATIONS CHECKED IN DELTAACV13309
C     REFERENCES--ANGULAR MOMENTUM IN QUANTUM MECHANICS, A.R.EDMONDS,   ACV13310
C                 PRINCETON                                             ACV13311
C     ------------------------------------------------------------------ACV13312
      IMPLICIT REAL*8(D)                                                ACV13313
      DJHR3=0.D0                                                        ACV13314
      ITMIN=MAX0(IABS(J1T-J9T),IABS(J2T-J6T),IABS(J4T-J8T))             ACV13315
      ITMAX=MIN0(J1T+J9T,J2T+J6T,J4T+J8T)                               ACV13316
      IF(ITMIN.GT.ITMAX)RETURN                                          ACV13317
      DO 10 IT=ITMIN,ITMAX,2                                            ACV13318
   10 DJHR3=DJHR3+DFLOAT(IT+1)*DRR3(J1T,J9T,J4T,J8T,IT,J7T)*            ACV13319
     1DRR3(J2T,J6T,J8T,J4T,IT,J5T)*DRR3(J1T,J9T,J2T,J6T,IT,J3T)         ACV13320
      DJHR3=DSQRT(DFLOAT((J3T+1)*(J6T+1)*(J7T+1)*(J8T+1)))*DJHR3        ACV13321
      RETURN                                                            ACV13322
      END                                                               ACV13323
      FUNCTION D3JR3(J1T,J2T,J3T,M1T,M2T,M3T)                           ACV13324
C     ------------------------------------------------------------------ACV13325
C     3J COEFFICIENTS FOR R3--TRIANGLE RELATIONS CHECKED IN DELTA       ACV13326
C     REFERENCES--ELEMENTARY THEORY OF ANGULAR MOMENTUM, M.E.ROSE, WILEYACV13327
C     ------------------------------------------------------------------ACV13328
      IMPLICIT REAL*8(D)                                                ACV13329
      COMMON/BKDF/DLOGF(0:2000)                                         ACV13330
      D3JR3=0.D0                                                        ACV13331
      IF(M1T+M2T-M3T.NE.0)GOTO 20                                       ACV13332
      DC=DELTA(J1T,J2T,J3T)                                             ACV13333
      IF(DC.EQ.12345D0)GOTO 20                                          ACV13334
      I1=J3T-J2T+M1T                                                    ACV13335
      I2=J3T-J1T-M2T                                                    ACV13336
      I3=J1T+J2T-J3T                                                    ACV13337
      I4=J1T-M1T                                                        ACV13338
      IF(BTEST(I4,0))GOTO 20                                            ACV13339
      I5=J2T+M2T                                                        ACV13340
      IF(BTEST(I5,0))GOTO 20                                            ACV13341
      ITMIN=MAX0(0,-I1,-I2)                                             ACV13342
      ITMAX=MIN0(I3,I4,I5)                                              ACV13343
      IF(ITMIN.GT.ITMAX)GOTO 20                                         ACV13344
      DTOP=(DC+DLOGF(J1T+M1T)+DLOGF(J1T-M1T)+                           ACV13345
     1DLOGF(J2T+M2T)+DLOGF(J2T-M2T)+DLOGF(J3T+M3T)+                     ACV13346
     2DLOGF(J3T-M3T))/DFLOAT(2)                                         ACV13347
      DO 10 IT=ITMIN,ITMAX,2                                            ACV13348
      DBOT=DLOGF(I3-IT)+DLOGF(I4-IT)+DLOGF(I5-IT)+                      ACV13349
     1DLOGF(IT)+DLOGF(I1+IT)+DLOGF(I2+IT)                               ACV13350
      DSUM=DEXP(DTOP-DBOT)                                              ACV13351
      IF(BTEST(IT,1))THEN                                               ACV13352
      D3JR3=D3JR3-DSUM                                                  ACV13353
      ELSE                                                              ACV13354
      D3JR3=D3JR3+DSUM                                                  ACV13355
      ENDIF                                                             ACV13356
   10 CONTINUE                                                          ACV13357
      IF(BTEST(I1-I2,1))D3JR3=-D3JR3                                    ACV13358
   20 RETURN                                                            ACV13359
      END                                                               ACV13360
      FUNCTION D6JR3(J1T,J2T,J3T,L1T,L2T,L3T)                           ACV13361
C     ------------------------------------------------------------------ACV13362
C     6J COEFFICIENTS FOR R3--TRIANGLE RELATION CHECKED IN DELTA        ACV13363
C     REFERENCES--THE 3-J AND 6-J SYMBOLS, M.ROTENBERG, R.BIVINS,       ACV13364
C                 N.METROPOLIS AND J.K.WOOTEN, MIT PRESS                ACV13365
C     ------------------------------------------------------------------ACV13366
      IMPLICIT REAL*8(D)                                                ACV13367
      COMMON/BKDF/DLOGF(0:2000)                                         ACV13368
      D6JR3=0.D0                                                        ACV13369
      DX=DELTA(J1T,J2T,J3T)                                             ACV13370
      IF(DX.EQ.12345D0)GOTO 20                                          ACV13371
      DC=DX                                                             ACV13372
      DX=DELTA(L1T,L2T,J3T)                                             ACV13373
      IF(DX.EQ.12345D0)GOTO 20                                          ACV13374
      DC=DX+DC                                                          ACV13375
      DX=DELTA(L1T,J2T,L3T)                                             ACV13376
      IF(DX.EQ.12345D0)GOTO 20                                          ACV13377
      DC=DX+DC                                                          ACV13378
      DX=DELTA(J1T,L2T,L3T)                                             ACV13379
      IF(DX.EQ.12345D0)GOTO 20                                          ACV13380
      DC=(DX+DC)/2.D0                                                   ACV13381
      I1=J3T+L3T-J1T-L1T                                                ACV13382
      I2=J3T+L3T-J2T-L2T                                                ACV13383
      I3=J1T+J2T+L1T+L2T+2                                              ACV13384
      I4=J1T+J2T-J3T                                                    ACV13385
      I5=L1T+L2T-J3T                                                    ACV13386
      I6=J1T+L2T-L3T                                                    ACV13387
      I7=L1T+J2T-L3T                                                    ACV13388
      ITMIN=MAX0(0,-I1,-I2)                                             ACV13389
      ITMAX=MIN0(I3,I4,I5,I6,I7)                                        ACV13390
      IF(ITMIN.GT.ITMAX)GOTO 20                                         ACV13391
      DO 10 IT=ITMIN,ITMAX,2                                            ACV13392
      DSUM=DEXP(DC+DLOGF(I3-IT)-(DLOGF(I4-IT)+DLOGF(I5-IT)+             ACV13393
     1DLOGF(I6-IT)+DLOGF(I7-IT)+DLOGF(IT)+DLOGF(I1+IT)+DLOGF(I2+IT)))   ACV13394
      IF(BTEST(IT,1))THEN                                               ACV13395
      D6JR3=D6JR3-DSUM                                                  ACV13396
      ELSE                                                              ACV13397
      D6JR3=D6JR3+DSUM                                                  ACV13398
      ENDIF                                                             ACV13399
   10 CONTINUE                                                          ACV13400
      IF(.NOT.BTEST(I3,1))D6JR3=-D6JR3                                  ACV13401
   20 RETURN                                                            ACV13402
      END                                                               ACV13403
      FUNCTION D9JR3(J1T,J2T,J3T,J4T,J5T,J6T,J7T,J8T,J9T)               ACV13404
C     ------------------------------------------------------------------ACV13405
C     9J COEFFICIENTS FOR R3--TRIANGLE RELATIONS CHECKED IN DELTA       ACV13406
C     REFERENCES--ANGULAR MOMENTUM, D.M.BRINK AND G.R.SATCHLER, OXFORD  ACV13407
C     ------------------------------------------------------------------ACV13408
      IMPLICIT REAL*8(D)                                                ACV13409
      D9JR3=0.D0                                                        ACV13410
      ITMIN=MAX0(IABS(J1T-J9T),IABS(J2T-J6T),IABS(J4T-J8T))             ACV13411
      ITMAX=MIN0(J1T+J9T,J2T+J6T,J4T+J8T)                               ACV13412
      IF(ITMIN.GT.ITMAX)RETURN                                          ACV13413
      DO 10 IT=ITMIN,ITMAX,2                                            ACV13414
   10 D9JR3=D9JR3+DFLOAT(IT+1)*DRR3(J1T,J9T,J4T,J8T,IT,J7T)*            ACV13415
     1DRR3(J2T,J6T,J8T,J4T,IT,J5T)*DRR3(J1T,J9T,J2T,J6T,IT,J3T)         ACV13416
      RETURN                                                            ACV13417
      END                                                               ACV13418
      FUNCTION DELTA(J1T,J2T,J3T)                                       ACV13419
C     ------------------------------------------------------------------ACV13420
C     DELTA FOR R3 ROUTINES--TRIANGLE RELATIONS CHECKED                 ACV13421
C     ------------------------------------------------------------------ACV13422
      IMPLICIT REAL*8(D)                                                ACV13423
      COMMON/BKDF/DLOGF(0:2000)                                         ACV13424
      DELTA=12345.D0                                                    ACV13425
      I1=J1T+J2T-J3T                                                    ACV13426
      IF(BTEST(I1,0))GOTO 10                                            ACV13427
      IF(I1.LT.0)GOTO 10                                                ACV13428
      I2=J2T+J3T-J1T                                                    ACV13429
      IF(I2.LT.0)GOTO 10                                                ACV13430
      I3=J3T+J1T-J2T                                                    ACV13431
      IF(I3.LT.0)GOTO 10                                                ACV13432
      DELTA=DLOGF(I1)+DLOGF(I2)+DLOGF(I3)-DLOGF(J1T+J2T+J3T+2)          ACV13433
   10 RETURN                                                            ACV13434
      END                                                               ACV13435
C********************                                                   ACV13436
