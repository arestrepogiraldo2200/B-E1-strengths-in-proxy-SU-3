C ACVISU3RME.  SU(3) REDUCED MATRIX ELEMENT PACKAGE.  C. BAHRI,           ACVI0000
C 1   J.P. DRAAYER.                                                       ACVI0000
C REF. IN COMP. PHYS. COMMUN. 83 (1994) 59                                ACVI0000
      program hwsgen                                                    ACVI0001
CB * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *ACVI0002
*                                                                      *ACVI0003
*                ***  Highest Weight State Generator ***               *ACVI0004
*                             ** (hwsgen) **                           *ACVI0005
*                                                                      *ACVI0006
CB * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *ACVI0007
C                                                                       ACVI0008
C Author:  Chairul Bahri                                                ACVI0009
C          Department of Physics and Astronomy                          ACVI0010
C          Louisiana State University                                   ACVI0011
C          Baton Rouge LA 70803 USA                                     ACVI0012
C                                                                       ACVI0013
C          e-mail address: phbahr @ lsuvm.sncc.lsu.edu (lsuvm.bitnet)   ACVI0014
C                           bahri @ rouge.phys.lsu.edu                  ACVI0015
C          phone: USA (504)-388-2261                                    ACVI0016
C                     (504)-388-6846                                    ACVI0017
C          fax:   USA (504)-388-5855                                    ACVI0018
C                                                                       ACVI0019
C ----------------------------------------------------------------------ACVI0020
C                                                                       ACVI0021
C Updates: 09/02/90: original ... IBM 3090/600                          ACVI0022
C          11/21/90: first test ... SUCCESSFUL !                        ACVI0023
C          01/13/91: single (LM MU) input.                              ACVI0024
C          01/22/91: multi (LM MU).                                     ACVI0025
C          08/15/91: implement weighted search tree (wst).              ACVI0026
C          11/15/91: SUN 4 work station.                                ACVI0027
C          11/19/91: SUCCESSFUL in SUN!                                 ACVI0028
C          11/13/92: SUCCESSFUL in IBM RS/6000 Model 560, AIX 3.2.3!    ACVI0029
C ----------------------------------------------------------------------ACVI0030
C                                                                       ACVI0031
C General description: PROGRAM 1                                        ACVI0032
C                                                                       ACVI0033
C    Main program for generating the highest weight state for a given   ACVI0034
C    set of SU(3)xSU(M) quantum numbers by applying the SU(3)xSU(M)     ACVI0035
C    raising operators and solving the homogeneous simultaneous linear  ACVI0036
C    equations resulted by the applying the raising operators.          ACVI0037
C                                                                       ACVI0038
C ----------------------------------------------------------------------ACVI0039
C  Reference:                                                           ACVI0040
C    C. Bahri, research notes.                                          ACVI0041
C ----------------------------------------------------------------------ACVI0042
      implicit real*8(d)                                                ACVI0043
C                                                                       ACVI0044
      parameter ( IHWSFILE = 8 )   ! unformatted hws file               ACVI0045
      character*20 cfile,          ! input filename                     ACVI0046
     1             cofile,         ! default hws filename               ACVI0047
     2             ca, cb          ! dummy names                        ACVI0048
C                  infile          ! formatted input file unit for sets ACVI0049
C                                  !   of quantum numbers               ACVI0050
C                  lfile           ! log file for intermediate results  ACVI0051
C                  irme            ! 0/1 for calculating rme later      ACVI0052
C                  ieta            ! oscillator shell number            ACVI0053
      parameter ( MXKETS =  50000, ! max dimension of kets              ACVI0054
     1            MXSOLN = 100000, ! max dimension of dsolmat           ACVI0055
     2            MXTREE = 140000 )! tree size                          ACVI0056
      common / HWSARR /                                                 ACVI0057
     1         dsolmat( MXSOLN ),  ! coefficients of hws for bit states ACVI0058
     2         lbtree(-10:MXTREE ),! binary tree for quantum labels     ACVI0059
     3         kets( MXKETS )      ! bit states                         ACVI0060
      dimension isolmat( 2*MXSOLN )                                     ACVI0061
      equivalence ( isolmat, dsolmat )                                  ACVI0062
C                                                                       ACVI0063
CW**********************************************************************ACVI0064
CW                                                                      ACVI0065
      write(6,'(/a,a/)') ' ENTERING *** HIGHEST WEIGHT STATE GENERATOR',ACVI0066
     1                   ' (HWSGEN) ***'                                ACVI0067
CW                                                                      ACVI0068
CW**********************************************************************ACVI0069
C                                                                       ACVI0070
      call headfile( IHWSFILE, 0, 0, infile, cofile, ca, cb, cfile )    ACVI0071
C+IBM infile = 15                                                       ACVI0072
      call heading( infile, IHWSFILE, lfile, irme )                     ACVI0073
      if(infile.ne.5) open( unit = infile, file = cfile ) ! old file    ACVI0074
      if(lfile.gt.0 .and. lfile.ne.6) open( unit = lfile,               ACVI0075
     1                                      file = 'hws.log' )          ACVI0076
      open( unit = IHWSFILE, file = cofile, form = 'unformatted' ) ! newACVI0077
      call hws( infile, lfile, irme,                                    ACVI0078
     1          dsolmat, kets, lbtree, MXSOLN, MXKETS, MXTREE,          ACVI0079
     2          ikets, isoln, ieta )                                    ACVI0080
C     Dump the required data to external file/tree.                     ACVI0081
      iloc  = lbtree(-10)                                               ACVI0082
      write(6,1020) iloc, ikets, isoln                                  ACVI0083
1020  format(/' TREE (state) STATISTICS:      NODES       kets',        ACVI0084
     1        '    dsolmat'/i36,2i11)                                   ACVI0085
      iloc  = iloc * ( lbtree(-3 ) + 3 )                                ACVI0086
      isoln = 2 * isoln                                                 ACVI0087
      write( IHWSFILE ) ieta                        ! signature         ACVI0088
      call wrdata(0, IHWSFILE, iloc, lbtree, -10)   ! labels            ACVI0089
      call wrdata(0, IHWSFILE, ikets, kets, 1)      ! bit states        ACVI0090
      call wrdata(0, IHWSFILE, isoln, isolmat, 1)   ! coefficients      ACVI0091
      close( unit = IHWSFILE, status = 'keep' )                         ACVI0092
      if(lfile.gt.0 .and. lfile.ne.6) close( unit = lfile,              ACVI0093
     1                                       status = 'keep' )          ACVI0094
      if(infile.ne.5) close( unit = infile, status = 'keep' )           ACVI0095
      stop                                                              ACVI0096
      end                                                               ACVI0097
      subroutine headfile( ihwsfile, iopbfile, irmefile, infile,        ACVI0098
     1                     chfile,   cofile,   crfile,   cfile )        ACVI0099
C                                                                       ACVI0100
C     Read file structure of rme. (Can be omitted.)                     ACVI0101
C                                                                       ACVI0102
      implicit logical(t)                                               ACVI0103
C                   ihwsfile       ! unformatted input file unit for hwsACVI0104
C                   iopbfile       ! unformatted input file unit for opbACVI0105
C                   irmefile       ! unformatted output file unit fr rmeACVI0106
C                   infile         ! formatted input file unit for sets ACVI0107
C                                  !   of quantum numbers               ACVI0108
      character*(*) cfile,         ! input filename                     ACVI0109
     1              chfile,        ! default hws filename               ACVI0110
     2              cofile,        ! default opb filename               ACVI0111
     3              crfile         ! default rme filename               ACVI0112
      character*20  ctemp          ! temporary name                     ACVI0113
C                                                                       ACVI0114
1     format(a)                                                         ACVI0115
2     format(' =>',10i5)                                                ACVI0116
3     format(' =>',4x,a)                                                ACVI0117
      infile = 15                                                       ACVI0118
      cfile = ' '                                                       ACVI0119
C                                                                       ACVI0120
      write(6,1)         ' Enter input file unit: '                     ACVI0121
      read (5,*)         infile                                         ACVI0122
      write(6,2)         infile                                         ACVI0123
      write(6,*)                                                        ACVI0124
      infile = iabs( infile )                                           ACVI0125
      if( infile.eq.0 .or. infile.eq.6 .or.                             ACVI0126
     1    infile.eq.4 .or. infile.eq.ihwsfile .or.                      ACVI0127
     2    infile.eq.iopbfile .or. infile.eq.irmefile ) then             ACVI0128
        write(6,1) ' Sorry, file unit 15 is assigned. '                 ACVI0129
        infile = 15                                                     ACVI0130
      endif                                                             ACVI0131
      if( infile .ne. 5 ) then                                          ACVI0132
        write(6,1)       ' Enter input file name: '                     ACVI0133
        read (5,'(a)')   ctemp                                          ACVI0134
        write(6,3)       ctemp                                          ACVI0135
        call compchr( ctemp )                                           ACVI0136
        if( .not. tchkdt( ctemp ) ) then                                ACVI0137
          chfile = 'hwsw.' // ctemp(1:14)                               ACVI0138
          cofile = 'opbw.' // ctemp(1:14)                               ACVI0139
          crfile = 'rmew.' // ctemp(1:14)                               ACVI0140
          cfile  = 'hwsir.' // ctemp(1:14)                              ACVI0141
        else                                                            ACVI0142
          chfile = 'hwsw.out'                                           ACVI0143
          cofile = 'opbw.out'                                           ACVI0144
          crfile = 'rmew.out'                                           ACVI0145
          cfile  = ctemp                                                ACVI0146
        endif                                                           ACVI0147
      else                                                              ACVI0148
        chfile = 'hwsw.out'                                             ACVI0149
        cofile = 'opbw.out'                                             ACVI0150
        crfile = 'rmew.out'                                             ACVI0151
      endif                                                             ACVI0152
      write(6,3) chfile                                                 ACVI0153
      write(6,3) cofile                                                 ACVI0154
      write(6,3) crfile                                                 ACVI0155
      write(6,3) cfile                                                  ACVI0156
      write(6,*)                                                        ACVI0157
C                                                                       ACVI0158
      return                                                            ACVI0159
C----*end of headfile*--------------------------------------------------ACVI0160
      end                                                               ACVI0161
      subroutine heading( infile, ihwsfile, lfile, irme )               ACVI0162
C                                                                       ACVI0163
C     Read input parameters for rme.                                    ACVI0164
C                                                                       ACVI0165
      implicit logical(t)                                               ACVI0166
      character*20 ctemp                                                ACVI0167
C               lfile              ! log file for intermediate results  ACVI0168
C                                                                       ACVI0169
1     format(a)                                                         ACVI0170
2     format(' =>',10i5)                                                ACVI0171
3     format(' =>',4x,a)                                                ACVI0172
      write(6,1)         ' Enter log file unit number:'                 ACVI0173
      write(6,1)         '       0: none'                               ACVI0174
      write(6,1)         '     1-4: some intermediate results'          ACVI0175
      write(6,1)         '     >=7: all intermediate results'           ACVI0176
      read (5,*,end=999) lfile                                          ACVI0177
      write(6,2)         lfile                                          ACVI0178
      write(6,*)                                                        ACVI0179
      write(6,1)         ' Do you want to calculate rme later?'         ACVI0180
      write(6,1)         '       0: no  (needs one set  of q.n.)'       ACVI0181
      write(6,1)         '       1: yes (needs two sets of q.n.)'       ACVI0182
      read (5,'(18a4)')  ctemp                                          ACVI0183
      write(6,3)         ctemp                                          ACVI0184
      write(6,*)                                                        ACVI0185
      call compchr( ctemp )                                             ACVI0186
      irme = 0                                                          ACVI0187
      if( ctemp(1:1).eq.'y'.or.ctemp(1:1).eq.'Y'.or.ctemp(1:1).eq.'1' ) ACVI0188
     1    irme = 1                                                      ACVI0189
C                                                                       ACVI0190
      lfile = iabs( lfile )                                             ACVI0191
      if( lfile.eq.5 .or. lfile.eq.infile .or. lfile.eq.ihwsfile )      ACVI0192
     1    call error(' HEADING: Fileunit has been assigned')            ACVI0193
C                                                                       ACVI0194
      return                                                            ACVI0195
999   write(6,1) ' I quit. Please enter numbers.'                       ACVI0196
      stop                                                              ACVI0197
C----*end of heading*-------------------------------------------------- ACVI0198
      end                                                               ACVI0199
      LOGICAL FUNCTION TCHKDT( CDOT )                                   ACVI0200
C                                                                       ACVI0201
C     Check whether there is a dot in a string.                         ACVI0202
C                                                                       ACVI0203
      CHARACTER*(*) CDOT                                                ACVI0204
      TCHKDT = .FALSE.                                                  ACVI0205
C     NLEN = 0                                                          ACVI0206
C     DO WHILE( .NOT. EOL )                                             ACVI0207
C       NLEN = NLEN + 1                                                 ACVI0208
      DO NLEN = 1, 20                                                   ACVI0209
        IF( CDOT( NLEN:NLEN ) .EQ. '.' ) THEN                           ACVI0210
          TCHKDT = .TRUE.                                               ACVI0211
          RETURN                                                        ACVI0212
        END IF                                                          ACVI0213
      END DO                                                            ACVI0214
      RETURN                                                            ACVI0215
      END                                                               ACVI0216
C                                                                       ACVI0217
      SUBROUTINE COMPCHR( CHARC )                                       ACVI0218
C                                                                       ACVI0219
C     Delete blanks before the first character.                         ACVI0220
C                                                                       ACVI0221
      CHARACTER*(*) CHARC                                               ACVI0222
      NLEN = 1                                                          ACVI0223
      DO WHILE( CHARC(NLEN:NLEN).EQ.' ' .AND. NLEN.LE.20 )              ACVI0224
        NLEN = NLEN + 1                                                 ACVI0225
      ENDDO                                                             ACVI0226
      CHARC = CHARC( NLEN: )                                            ACVI0227
      RETURN                                                            ACVI0228
      END                                                               ACVI0229
C                                                                       ACVI0230
      SUBROUTINE WRTOUT(LFILE,ICALL,MLTP,IQNUMS,DASH,LITER)             ACVI0231
C                                                                       ACVI0232
C     Write out the heading for each irreps.                            ACVI0233
C                                                                       ACVI0234
      DIMENSION IQNUMS(*)                                               ACVI0235
      CHARACTER DASH(78)                                                ACVI0236
      CHARACTER*(*) LITER                                               ACVI0237
      WRITE(LFILE,'(A)') LITER                                          ACVI0238
      WRITE(LFILE,'(A,I5)') ' STATE NUMBER ', ICALL                     ACVI0239
      WRITE(LFILE,*) DASH                                               ACVI0240
      WRITE(LFILE,1000) (' ',I,I=1,MLTP)                                ACVI0241
      WRITE(LFILE,1010) (IQNUMS(I),I=1,MLTP+2)                          ACVI0242
      WRITE(LFILE,*) DASH                                               ACVI0243
1000  FORMAT(' ',T6,'Lm',T12,'Mu',T15,10(A,2X,'f',I1,'~'))              ACVI0244
1010  FORMAT(' ',12I6)                                                  ACVI0245
      RETURN                                                            ACVI0246
      END                                                               ACVI0247
CB * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *ACVI0248
*                                                                      *ACVI0249
*                ***  HIGHEST WEIGHT STATE GENERATOR ***               *ACVI0250
*                              ** (HWS) **                             *ACVI0251
*                                                                      *ACVI0252
CB * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *ACVI0253
C                                                                       ACVI0254
C Author:  Chairul Bahri                                                ACVI0255
C          Department of Physics and Astronomy                          ACVI0256
C          Louisiana State University                                   ACVI0257
C          Baton Rouge LA 70803 USA                                     ACVI0258
C                                                                       ACVI0259
C          E-MAIL: PHBAHR @ LSUVM.SNCC.LSU.EDU                          ACVI0260
C                  bahri  @ rouge.phys.lsu.edu                          ACVI0261
C          PHONE:  USA (504)-388-2261                                   ACVI0262
C                      (504)-388-6846                                   ACVI0263
C          FAX:    USA (504)-388-5855                                   ACVI0264
C                                                                       ACVI0265
C ----------------------------------------------------------------------ACVI0266
C                                                                       ACVI0267
C Updates: 09/02/90: original ... IBM 3090/600                          ACVI0268
C          11/21/90: first test ... SUCCESSFUL !                        ACVI0269
C          01/13/91: single (LM MU) input.                              ACVI0270
C          01/22/91: multi (LM MU).                                     ACVI0271
C          08/15/91: using weighted search tree (wst).                  ACVI0272
C          11/15/91: SUN 4.                                             ACVI0273
C          12/04/91: in subroutine form ...                             ACVI0274
C          11/14/92: IBM RS/6000-560.                                   ACVI0275
C                                                                       ACVI0276
C ----------------------------------------------------------------------ACVI0277
C                                                                       ACVI0278
C General description: SUB-PROGRAM 1                                    ACVI0279
C                                                                       ACVI0280
C    Main subprogram for generating the highest weight state for a givenACVI0281
C    set of SU(3)xSU(M) quantum numbers by applying the SU(3)xSU(M)     ACVI0282
C    raising operators and solving the homogeneous simultaneous linear  ACVI0283
C    equations resulted from the applying the raising operators.        ACVI0284
C                                                                       ACVI0285
C ----------------------------------------------------------------------ACVI0286
C  Reference:                                                           ACVI0287
C    C. Bahri, research notes.                                          ACVI0288
C ----------------------------------------------------------------------ACVI0289
      subroutine hws( infile, logfile, irme,                            ACVI0290
     ]                dsolmat, kets, lbtree, mxsoln, mxkets, mxtree,    ACVI0291
     ]                ikets, isoln, neta )                              ACVI0292
      implicit real*8(d), logical(t)                                    ACVI0293
C                                                                       ACVI0294
      character dash(78)           ! ---...---                          ACVI0295
C               infile             ! formatted input file unit for sets ACVI0296
C                                  !   of quantum numbers               ACVI0297
C               logfile            ! log file for intermediate results  ACVI0298
C               irme               ! 0/1 for calculating rme later      ACVI0299
      dimension dsolmat( * ),      ! coefficients of hws for bit states ACVI0300
     1          lbtree( -10:* ),   ! binary tree for quantum labels     ACVI0301
     2          kets( * )          ! bit states                         ACVI0302
C               mxsoln             ! max dimension of dsolmat           ACVI0303
C               mxkets             ! max dimension of kets              ACVI0304
C               mxtree             ! max dimension of lbtree            ACVI0305
C               isoln              ! # elements in dsolmat              ACVI0306
C               ikets              ! # elements in kets                 ACVI0307
C               neta               ! oscillator shell number            ACVI0308
      common / HWSCON / lfile,     ! = logfile                          ACVI0309
     1                  ieta,      ! = neta                             ACVI0310
     2                  nwords,    ! # words in one bit state           ACVI0311
     3                  nbits,     ! bit length for # levels in shell   ACVI0312
     4                  mbits      ! length of bits in a word           ACVI0313
C                       mltp       ! internal symmetry multiplicity     ACVI0314
      dimension iqnums( 24 ),      ! array for quantum numbers          ACVI0315
     1          label( 3 )         ! temp array of packed q.n. for tree ACVI0316
      real title( 18 )             ! title ( temporary )                ACVI0317
C                                                                       ACVI0318
      data nbit16 / 65535 /                                             ACVI0319
C     data nbit16 / zffff /        ! for IBM-3090 only                  ACVI0320
C                                                                       ACVI0321
      lfile = logfile                                                   ACVI0322
      tfile = lfile .gt. 0                                              ACVI0323
      trme  = irme .eq. 1                                               ACVI0324
C                                                                       ACVI0325
      write(6,'(a)') ' Enter internal symmetry multiplicity: '          ACVI0326
      write(6,'(a)') ' (e.g. 2 for spin, 4 for spin-isospin) '          ACVI0327
      read (5,*,end=999) mltp                                           ACVI0328
      write(6,'(a,i5)')  ' =>', mltp                                    ACVI0329
7     format(18a4)                                                      ACVI0330
      if( infile.eq.5 ) then                                            ACVI0331
        write(6,*) ' ieta'                                              ACVI0332
        read(infile,*) ieta                                             ACVI0333
      else                                                              ACVI0334
        read(infile,7,end=999) title         ! shell number             ACVI0335
        read(infile,*,end=999) ieta                                     ACVI0336
      end if                                                            ACVI0337
      if(tfile) write(lfile,'(a,i4)') ' Oscillator shell number =',ieta ACVI0338
      neta = ieta                                                       ACVI0339
C                                                                       ACVI0340
C     Start ...                                                         ACVI0341
C                                                                       ACVI0342
      call initlz( mltp ) ! Initialization                              ACVI0343
C                                                                       ACVI0344
      ikets = 0                                                         ACVI0345
      isoln = 0                                                         ACVI0346
      mltpp2 = mltp + 2                                                 ACVI0347
      if( trme ) then                                                   ACVI0348
        mltpmax = mltpp2 + mltpp2          ! bra & ket                  ACVI0349
      else                                                              ACVI0350
        mltpmax = mltpp2                   ! ket only                   ACVI0351
      endif                                                             ACVI0352
C                                                                       ACVI0353
C     Setup lbtree for quantum numbers: 2 data.                         ACVI0354
      nkeys = (mltp + 1)/4 + 1                                          ACVI0355
      n = mxtree/(nkeys+5) - 1                                          ACVI0356
      call tsetlf( lbtree, n, nkeys, 2)                                 ACVI0357
      nkeym1 = nkeys - 1                                                ACVI0358
      n = mod( mltpp2, 4 )                                              ACVI0359
      if( n .eq. 0 ) n = 4                                              ACVI0360
C     Read SU(3)xSU(M) quantum numbers.                                 ACVI0361
      if( infile.eq.5 ) then                                            ACVI0362
        write(6,*) ' Lm Mu f~''s'                                       ACVI0363
      else                                                              ACVI0364
        read(infile,7,end=999) title                                    ACVI0365
      endif                                                             ACVI0366
C                                                                       ACVI0367
      idash = 6 * (mltpp2 + 1)                                          ACVI0368
      do i = 1, idash                                                   ACVI0369
        dash( i ) = '-'                                                 ACVI0370
      end do                                                            ACVI0371
      do i = idash+1, 78                                                ACVI0372
        dash( i ) = ' '                                                 ACVI0373
      end do                                                            ACVI0374
C                                                                       ACVI0375
      icall = 0                                                         ACVI0376
C                                                                       ACVI0377
      do while( .true. )  ! * * * * * * * * * * * * * * * * * * * * * * ACVI0378
C                                                                       ACVI0379
        if( infile.eq.5 ) then                                          ACVI0380
          read(infile,*) (iqnums(i),i=1,mltpmax)                        ACVI0381
        else                                                            ACVI0382
          read(infile,*,end=999) (iqnums(i),i=1,mltpmax)                ACVI0383
        end if                                                          ACVI0384
C                                                                       ACVI0385
        iq = 0                                                          ACVI0386
        kadd = 0                                                        ACVI0387
C                                                                       ACVI0388
        do k = 0, irme                                                  ACVI0389
C         Pack labels: each 4.                                          ACVI0390
          do i = 1, nkeys                                               ACVI0391
            label( i ) = 0                                              ACVI0392
          end do                                                        ACVI0393
          do i = 1, nkeym1                                              ACVI0394
            do j = 1, 4                                                 ACVI0395
              iq = iq + 1                                               ACVI0396
              label( i ) = ior( ishft(label(i),8), iqnums(iq) )         ACVI0397
            end do                                                      ACVI0398
          end do                                                        ACVI0399
          do j = 1, n                                                   ACVI0400
            iq = iq + 1                                                 ACVI0401
            label(nkeys) = ior(ishft(label(nkeys),8),iqnums(iq))        ACVI0402
          end do                                                        ACVI0403
          call tchk( label, lbtree, *110 ) ! Check labels.              ACVI0404
            icall = icall + 1                                           ACVI0405
            if( lbtree(-10) .le. lbtree(-9) ) then                      ACVI0406
              iloc = lbtree(-5)                                         ACVI0407
              lbtree(iloc+lbtree(-2)) = 1  ! equal priority             ACVI0408
              kf = kadd + 1                                             ACVI0409
C             Write out quantum numbers.                                ACVI0410
              call wrtout( 6, icall, mltp, iqnums(kf), dash, ' ')       ACVI0411
              if( tfile ) call wrtout( lfile, icall, mltp, iqnums(kf),  ACVI0412
     1                                 dash, ' ')                       ACVI0413
C             Check quantum numbers.                                    ACVI0414
              do i = kf, kadd+3                                         ACVI0415
                if( iqnums(i) .lt. 0 )                                  ACVI0416
     1            call attn(' HWS: Negative Quantum Numbers.',*110)     ACVI0417
              end do                                                    ACVI0418
              do i = kadd+4, kadd+mltpp2                                ACVI0419
                if( iqnums(i) .lt. 0 )                                  ACVI0420
     1            call attn(' HWS: Negative Quantum Numbers.',*110)     ACVI0421
                if( iqnums(i) .gt. iqnums(i-1) )                        ACVI0422
     1            call attn(' HWS: f2 > f1.',*110)                      ACVI0423
              end do                                                    ACVI0424
C                                                                       ACVI0425
              maxket = mxkets - ikets                                   ACVI0426
              maxsol = mxsoln - isoln                                   ACVI0427
C             Generate basis.                                           ACVI0428
              call genmpb( mltp, iqnums(kf), kets(ikets+1),             ACVI0429
     1                     nbkets, maxket, *999 )                       ACVI0430
C             Apply SU(3)xSU(M) raising operators.                      ACVI0431
              if( nbkets .gt. maxsol )                                  ACVI0432
     1          call attn(' HWS: # bit states too large. Skip.',*110)   ACVI0433
              if( nbkets .gt. 0 ) call appaij( mltp, iqnums(kf+2),      ACVI0434
     1            kets(ikets+1), dsolmat(isoln+1),                      ACVI0435
     2            nbkets, nalpha, maxsol, trme, *999 )                  ACVI0436
C             Save in lbtree.                                           ACVI0437
              if(nbkets.le.nbit16 .and. nalpha.le.nbit16) then          ACVI0438
                iloc = iloc + lbtree(-4) + 1                            ACVI0439
                lbtree(iloc) = ior( ishft(ikets,16), nbkets )           ACVI0440
                lbtree(iloc+1) = ior( ishft(isoln,16), nalpha )         ACVI0441
                call tins( label, lbtree )                              ACVI0442
                ikets = ikets + nbkets*mltp                             ACVI0443
                isoln = isoln + nbkets*nalpha                           ACVI0444
              else                                                      ACVI0445
                write(6,*) 'This IRREP is not saved in lbtree.'         ACVI0446
              endif                                                     ACVI0447
            else                                                        ACVI0448
              call wrtout( 6, icall, mltp, iqnums(kf), dash,            ACVI0449
     1           ' Increase nodes for tree, MXTREE, not processing:')   ACVI0450
            endif                                                       ACVI0451
            write(6,'(a,2i11)') ' Check tree:',icall,lbtree(-10)        ACVI0452
            write(6,'(a,2i11)') '       data:',ikets,isoln              ACVI0453
110       continue                                                      ACVI0454
          kadd = kadd + mltpp2                                          ACVI0455
C                                                                       ACVI0456
        enddo                                                           ACVI0457
C                                                                       ACVI0458
      enddo             ! * * * * * * * * * * * * * * * * * * * * * * * ACVI0459
999   continue                                                          ACVI0460
      return                                                            ACVI0461
C ---*end of hws*-------------------------------------------------------ACVI0462
      end                                                               ACVI0463
      SUBROUTINE CMPRSS(NDIM1,IVAR2,NDIM2,NDIM3,DVAR2,*)                ACVI0464
C                                                                       ACVI0465
C     Compress for nonzero elements of the coefficient only.            ACVI0466
C                                                                       ACVI0467
      IMPLICIT REAL*8(D),LOGICAL(T)                                     ACVI0468
      PARAMETER (DZERO=1.D-12,MXCOMP=5000)                              ACVI0469
      DIMENSION IVAR2(NDIM1,*),DVAR2(*)                                 ACVI0470
C     DIMENSION IVAR2(NDIM1,NDIM2), DVAR2(NDIM2,NDIM3)                  ACVI0471
      COMMON /COMPA / TCOMPA(MXCOMP)                                    ACVI0472
      IF(NDIM2.GT.MXCOMP) THEN                                          ACVI0473
        CALL ATTN(' CMPRSS: Sorry, cannot compress.',*10)               ACVI0474
10      RETURN 1                                                        ACVI0475
      ENDIF                                                             ACVI0476
      NEWDIM=0                                                          ACVI0477
      DO I=1,NDIM2                                                      ACVI0478
        TCOMP=.FALSE.                                                   ACVI0479
        ICOEF=I                                                         ACVI0480
        DO J=1,NDIM3                                                    ACVI0481
          IF(DABS(DVAR2(ICOEF)) .GT. DZERO) THEN                        ACVI0482
            TCOMP=.TRUE.                                                ACVI0483
            GO TO 20                                                    ACVI0484
          ENDIF                                                         ACVI0485
          ICOEF=ICOEF+NDIM2                                             ACVI0486
        ENDDO                                                           ACVI0487
20      CONTINUE                                                        ACVI0488
        IF(TCOMP) THEN                                                  ACVI0489
          NEWDIM=NEWDIM+1                                               ACVI0490
C         Compress kets.                                                ACVI0491
          DO J=1,NDIM1                                                  ACVI0492
            IVAR2(J,NEWDIM)=IVAR2(J,I)                                  ACVI0493
          ENDDO                                                         ACVI0494
        ENDIF                                                           ACVI0495
        TCOMPA(I)=TCOMP                                                 ACVI0496
      ENDDO                                                             ACVI0497
C     Compress coefficients.                                            ACVI0498
      ICOEF =0                                                          ACVI0499
      NWCOEF=0                                                          ACVI0500
      DO J=1,NDIM3                                                      ACVI0501
        DO I=1,NDIM2                                                    ACVI0502
          ICOEF=ICOEF+1                                                 ACVI0503
          IF(TCOMPA(I)) THEN                                            ACVI0504
            NWCOEF=NWCOEF+1                                             ACVI0505
            DVAR2(NWCOEF)=DVAR2(ICOEF)                                  ACVI0506
          ENDIF                                                         ACVI0507
        ENDDO                                                           ACVI0508
CW      WRITE(6,*) 'NEWDIM',DFLOAT(NWCOEF)/DFLOAT(J)                    ACVI0509
      ENDDO                                                             ACVI0510
      NDIM2=NEWDIM                                                      ACVI0511
CW    WRITE(6,*) 'NEWDIM',NEWDIM                                        ACVI0512
      RETURN                                                            ACVI0513
C ---*end of CMPRSS*----------------------------------------------------ACVI0514
      END                                                               ACVI0515
CB----------------------------------------------------------------------ACVI0516
C                                                                       ACVI0517
C                          **************                               ACVI0518
C                          *** INITLZ ***                               ACVI0519
C                          **************                               ACVI0520
C                                                                       ACVI0521
C ----------------------------------------------------------------------ACVI0522
C Author:  Chairul Bahri                                                ACVI0523
C ----------------------------------------------------------------------ACVI0524
C Updates: 07/90 ==> Rewrite of the original (Reske/Lopes'z version.)   ACVI0525
C          A major change: spacelike bit representation so it can apply ACVI0526
C          for SU(M) internal symmetry.                                 ACVI0527
C          07/09/90: generating spacelike many particle basis.          ACVI0528
C          08/04/90: using U(M) counter generator.                      ACVI0529
C          08/18/90: running successfully for HWSGEN !!!                ACVI0530
C          11/04/90: implemented in RMECODE.                            ACVI0531
C          06/25/91: compress for non-zero coefficients of basis only.  ACVI0532
C          11/15/91: SUN 4. (combined BASGEN and AIJOPS.)               ACVI0533
C ----------------------------------------------------------------------ACVI0534
C                                                                       ACVI0535
C    The subroutine INITLZ initializes all single-particle basis of     ACVI0536
C    U(3)xU(M) states for a shell number IETA (max of 10 because of the ACVI0537
C    max number of particles occupying each state ==> 32 bits) and the  ACVI0538
C    coefficients for raising operators Czx, Cxy, and U(M)-Aij.         ACVI0539
C                                                                       ACVI0540
C    M can be any number as long as it satisfies the condition that the ACVI0541
C    total size does not exceed 32 bits, e.g.:                          ACVI0542
C       IETA          N        M                                        ACVI0543
C         0           1       32  > NOT ENOUGH STORAGE TO STORE         ACVI0544
C         1           3       16  > ALL CONFIGURATIONS                  ACVI0545
C         2           6       10                                        ACVI0546
C         .           .        .                                        ACVI0547
C         .           .        .                                        ACVI0548
C         .           .        .                                        ACVI0549
C        10          66        4                                        ACVI0550
C --------------------------------------------------------------------  ACVI0551
C                                                                       ACVI0552
      SUBROUTINE INITLZ( NMAX )                                         ACVI0553
      IMPLICIT REAL*8(D), LOGICAL(T)                                    ACVI0554
C     NMAX   : internal symmetry multiplicity in a state.               ACVI0555
C     NMLVLS : size of significant bits in every word.                  ACVI0556
C                                                                       ACVI0557
      COMMON /HWSCON/ LFILE, IETA, NWORDS, NBITS, MBITS                 ACVI0558
C     LFILE  : log file for intermediate results.                       ACVI0559
C     IETA   : oscillator shell number.                                 ACVI0560
C     NWORDS : number of words of bits in KETS occupied by one bit      ACVI0561
C              state.                                                   ACVI0562
C     NBITS  : length of bits for number of levels in the shell.        ACVI0563
C     MBITS  : length of bits for a word.                               ACVI0564
C                                                                       ACVI0565
      PARAMETER ( NSTSMX = 66, MXBITS = 32, MXCONS = 1024 )             ACVI0566
C     NSTSMX : max number of levels for a given shell IETA.             ACVI0567
C     MXBITS : number of bits in one word.                              ACVI0568
C     MXCONS : max number of U(M) configurations = 2**10                ACVI0569
      COMMON /SUNQN / ICONU3(NSTSMX),IPTRMP(12),NSUN(MXCONS)            ACVI0570
     *               ,ISUN(MXCONS)                                      ACVI0571
C     Spatial U(3) counters for U(N) IKET distribution.                 ACVI0572
C     ICONU3 : U(3) IKETs for all states in a given shell. ==> packed   ACVI0573
C     Internal counters for U(M) IKET distribution.                     ACVI0574
C     NSUN   : binary representations of U(M) configuration.            ACVI0575
C     ISUN   : value of single U(M) configuration. ==> packed           ACVI0576
C     IPTRMP : pointer for number of particles.                         ACVI0577
C     e.g.:                                                             ACVI0578
C     U(2) COUNTERS:                                                    ACVI0579
C              1 for zero particles,                                    ACVI0580
C           2--3 for one particle,                                      ACVI0581
C              4 for two particles                                      ACVI0582
C     DATA IPTRMP/1,2,4,5, 2*0/                                         ACVI0583
C     Binary representation of SU(2) or U(2) IKETs.                     ACVI0584
C     NO.   STATE    2*MS    N1 N2 (BINARY) --> DECIMAL                 ACVI0585
C            1         0      0  0                 0                    ACVI0586
C            2        +1      1  0                 2                    ACVI0587
C            3        -1      0  1                 1                    ACVI0588
C            4         0      1  1                 3                    ACVI0589
C     DATA ISUN/0,Z100,   1,Z101, 12*0/                                 ACVI0590
C     DATA NSUN/0,   2,   1,   3, 12*0/                                 ACVI0591
C                                                                       ACVI0592
C     SU(3) raising operator coefficients.                              ACVI0593
      COMMON /SU3COF/ DCOEFZX(NSTSMX),DCOEFXY(NSTSMX),NMLVLS(3)         ACVI0594
C     DCOEFZX : Czx coefficients.                                       ACVI0595
C     DCOEFXY : Cxy coefficients.                                       ACVI0596
C     NMLVLS  : size of significant bits in every word.                 ACVI0597
C                                                                       ACVI0598
C     Tree for bit patterns.                                            ACVI0599
      PARAMETER ( MXACTS = 100000, MXT = 140000 )                       ACVI0600
C     MXACTS : maximum dimension of distinct bit patterns.              ACVI0601
C     MXT    : max number of nodes in trees.                            ACVI0602
      COMMON /PATTRN/ DACTS(MXACTS),IACTS(MXACTS),IRSPTR(MXACTS)        ACVI0603
     *               ,KETU3(MXACTS,3),KETREE(-10:MXT),IACTMX            ACVI0604
C     DACTS  : buffer linked-list for coefficients.                     ACVI0605
C     IACTS  : linked-list index for coefficients.                      ACVI0606
C     IRSPTR : signature for raising operators.                         ACVI0607
C     KETU3/KETUN : bit pattern as key.                                 ACVI0608
C     KETREE : binary tree for bit patterns.                            ACVI0609
C     IACTMX : last non-empty index.                                    ACVI0610
      CHARACTER*35 BIN, PACK                                            ACVI0611
C     BIN    : Binary representations (in character type).              ACVI0612
C                                                                       ACVI0613
C     Packing Functions == 3 labels                                     ACVI0614
      IPACK(I1,I2,I3)=IOR(I3,ISHFT(IOR(I2,ISHFT(I1,10)),10))            ACVI0615
C                                                                       ACVI0616
      TFILE=LFILE.GT.0                                                  ACVI0617
      IF(IETA.GT.10)CALL ERROR(' INITLZ: Shell Number greater than 10') ACVI0618
C due to the structure of U(3) packing, NQUANTAmax=880 for IETA=10)     ACVI0619
      IF(NMAX.GT.10)CALL ERROR(' INITLZ: Multiplicity greater than 10') ACVI0620
C due to the size of the storage, MXCONS=1024 for NMAX=10)              ACVI0621
C                                                                       ACVI0622
C     Initialize U(3) counters and raising-operator coefficients.       ACVI0623
      NSTS=0                                                            ACVI0624
      IF(TFILE) WRITE(LFILE,'(/A/)') ' ** U( 3 ) space counters:'       ACVI0625
      DO NZ=IETA,0,-1                                                   ACVI0626
        DO NX=IETA-NZ,0,-1                                              ACVI0627
          NY=IETA-NZ-NX                                                 ACVI0628
          NSTS=NSTS+1                                                   ACVI0629
CW----------------------------------------------------------------------ACVI0630
          IF(TFILE) WRITE(LFILE,'(2I8,2I4)') NSTS,NZ,NX,NY              ACVI0631
CW----------------------------------------------------------------------ACVI0632
          ICONU3( NSTS )  = IPACK( NZ, NX, NY )                         ACVI0633
          DCOEFZX( NSTS ) = DSQRT( DFLOAT( (NZ+1)*NX ) )                ACVI0634
          DCOEFXY( NSTS ) = DSQRT( DFLOAT( (NX+1)*NY ) )                ACVI0635
        END DO                                                          ACVI0636
      END DO                                                            ACVI0637
C                                                                       ACVI0638
      MAXBIT = MXBITS - 1                                               ACVI0639
      DO J = 0, MAXBIT                                                  ACVI0640
        IF(BTEST(NSTS,J)) NBITS=J+1                                     ACVI0641
        IF(BTEST(MAXBIT,J)) MBITS=J+1                                   ACVI0642
      END DO                                                            ACVI0643
      MBITS = -MBITS                                                    ACVI0644
      IF((NMAX*NBITS).GT.MXBITS)                                        ACVI0645
     1  CALL ERROR(' INITLZ: Too many configurations')                  ACVI0646
      NWORDS=(NSTS+MXBITS-1)/MXBITS                                     ACVI0647
C                                                                       ACVI0648
C     Pointers for applications or SU(M) raising operators.             ACVI0649
      DO I = 1, NWORDS                                                  ACVI0650
        NFINAL = MIN0( NSTS, MXBITS )                                   ACVI0651
        NMLVLS( I ) = NFINAL - 1                                        ACVI0652
        NSTS = NSTS - NFINAL                                            ACVI0653
      END DO                                                            ACVI0654
C                                                                       ACVI0655
C     Initialize nodes for tree KETREE: 2 data.                         ACVI0656
      NPAIRS = 2 * NWORDS                                               ACVI0657
      KNODES = MXT / (NPAIRS + 5)                                       ACVI0658
      CALL TSETLF( KETREE, KNODES, NPAIRS, 2 )                          ACVI0659
      IACTMX = 0                                                        ACVI0660
C                                                                       ACVI0661
C     Initialize all U(M) counters.                                     ACVI0662
      NMAXM1=NMAX-1                                                     ACVI0663
      NMAXM2=NMAX-2                                                     ACVI0664
      NMAXP2=NMAX+2                                                     ACVI0665
      NSUN(1)=0                                                         ACVI0666
      IPTRMP(1)=1                                                       ACVI0667
      IPTR=1                                                            ACVI0668
      IRUN=1                                                            ACVI0669
      DO I=1,NMAX                                                       ACVI0670
        NSUNTM=0                                                        ACVI0671
        IPTR=IPTR+1                                                     ACVI0672
        IPTRMP(IPTR)=IRUN+1                                             ACVI0673
        DO J=0,I-1                                                      ACVI0674
          NSUNTM=IBSET(NSUNTM,J)                                        ACVI0675
        END DO                                                          ACVI0676
100     CONTINUE                                                        ACVI0677
          IRUN=IRUN+1                                                   ACVI0678
          NSUN(IRUN)=NSUNTM                                             ACVI0679
          DO J=NMAXM2,0,-1                                              ACVI0680
            IF(BTEST(NSUNTM,J))THEN                                     ACVI0681
              NSUNTM=IBCLR(NSUNTM,J)                                    ACVI0682
              NTOT=0                                                    ACVI0683
              DO K=0,J                                                  ACVI0684
                IF(BTEST(NSUNTM,K)) NTOT=NTOT+1                         ACVI0685
              END DO                                                    ACVI0686
              DO K=J+1,NMAXM1                                           ACVI0687
                IF(NTOT.LT.I) THEN                                      ACVI0688
                  NSUNTM=IBSET(NSUNTM,K)                                ACVI0689
                  NTOT=NTOT+1                                           ACVI0690
                ELSE                                                    ACVI0691
                  NSUNTM=IBCLR(NSUNTM,K)                                ACVI0692
                ENDIF                                                   ACVI0693
              END DO                                                    ACVI0694
              IF(NTOT.EQ.I)THEN                                         ACVI0695
                GOTO 100                                                ACVI0696
              ELSE                                                      ACVI0697
                NSUNTM=IBSET(NSUNTM,J)                                  ACVI0698
              ENDIF                                                     ACVI0699
            ENDIF                                                       ACVI0700
          END DO                                                        ACVI0701
      END DO                                                            ACVI0702
      IPTRMP(IPTR+1)=IRUN+1                                             ACVI0703
C     Pack U(M) counters for checking many particle basis.              ACVI0704
      DO I=1,IRUN                                                       ACVI0705
        ISUN(I)=0                                                       ACVI0706
      END DO                                                            ACVI0707
      DO I=2,IRUN                                                       ACVI0708
        IF(BTEST(NSUN(I),0)) ISUN(I)=1                                  ACVI0709
        DO J=1,NMAXM1                                                   ACVI0710
          IF(BTEST(NSUN(I),J))THEN                                      ACVI0711
            ISUN(I)=IOR(ISHFT(ISUN(I),NBITS),1)                         ACVI0712
          ELSE                                                          ACVI0713
            ISUN(I)=ISHFT(ISUN(I),NBITS)                                ACVI0714
          ENDIF                                                         ACVI0715
        END DO                                                          ACVI0716
      END DO                                                            ACVI0717
CW----------------------------------------------------------------------ACVI0718
CW    Write out the configurations                                      ACVI0719
      IF(TFILE) THEN                                                    ACVI0720
      WRITE(LFILE,1000)NMAX,IRUN                                        ACVI0721
 1000 FORMAT(/' ** U(',I2,' ) internal counters:',I7/)                  ACVI0722
      KEMPTY=NBITS+1                                                    ACVI0723
      KM=NMAX*NBITS                                                     ACVI0724
      KMAX=KM+NMAXM1                                                    ACVI0725
      DO I=1,IRUN                                                       ACVI0726
        CALL SETBIN(LFILE,.FALSE.,BIN,NMAX+1,-1,35,35)                  ACVI0727
        CALL DBOVRC(NMAX,NSUN(I))                                       ACVI0728
        IF(KMAX.LE.35)THEN                                              ACVI0729
          CALL SETBIN(LFILE,.FALSE.,PACK,KEMPTY,KM,35,35)               ACVI0730
          CALL DBOVRC(KM,ISUN(I))                                       ACVI0731
        ENDIF                                                           ACVI0732
        WRITE(LFILE,'(2I8,4X,A,A)')I,NSUN(I),BIN,PACK                   ACVI0733
      END DO                                                            ACVI0734
      WRITE(LFILE,'(/A,I1,A/)')                                         ACVI0735
     1    ' Pointer for U(',NMAX,') configurations:'                    ACVI0736
      WRITE(LFILE,'(10I4)')(IPTRMP(I),I=1,IPTR+1)                       ACVI0737
      ENDIF                                                             ACVI0738
CW----------------------------------------------------------------------ACVI0739
C     Setup binary I/O for the purpose in the calling program.          ACVI0740
      CALL SETBIN(LFILE,.TRUE.,BIN,11,-1,35,35)                         ACVI0741
      WRITE(6,*)                                                        ACVI0742
      RETURN                                                            ACVI0743
C ---*end of INITLZ*----------------------------------------------------ACVI0744
      END                                                               ACVI0745
CB----------------------------------------------------------------------ACVI0746
C                                                                       ACVI0747
C                          **************                               ACVI0748
C                          *** GENMPB ***                               ACVI0749
C                          **************                               ACVI0750
C                                                                       ACVI0751
C ----------------------------------------------------------------------ACVI0752
C Author:  Chairul Bahri                                                ACVI0753
C ----------------------------------------------------------------------ACVI0754
C Updates: 07/90 ==> Rewrite of the original (Reske/Lopes'z version.)   ACVI0755
C          A major change: spacelike bit representation so it can apply ACVI0756
C          for SU(M) internal symmetry.                                 ACVI0757
C          07/09/90: generating spacelike many particle basis.          ACVI0758
C          08/04/90: using U(M) counter generator.                      ACVI0759
C          08/18/90: running successfully for HWSGEN !!!                ACVI0760
C          11/04/90: implemented in RMECODE.                            ACVI0761
C          11/15/91: SUN 4.                                             ACVI0762
C ----------------------------------------------------------------------ACVI0763
C                                                                       ACVI0764
C    The subroutine GENMPB generates the many-particle basis of         ACVI0765
C    U(3)xU(M) states for a shell number IETA.) All the possible U(M)   ACVI0766
C    states which may occur for each U(3) irrep are generated for the   ACVI0767
C    possible distribution of Np particles in the shell IETA. Only thoseACVI0768
C    many-body states which belong to the specified input total         ACVI0769
C    U(3)xU(M) irrep are generated.                                     ACVI0770
C                                                                       ACVI0771
C    ... see INITLZ ...                                                 ACVI0772
C                                                                       ACVI0773
C ----------------------------------------------------------------------ACVI0774
C                                                                       ACVI0775
      subroutine genmpb( nmax, iqnums, kets, nkets, mxkets, * )         ACVI0776
      implicit logical(t)                                               ACVI0777
C               nmax               ! internal symmetry multiplicity=M   ACVI0778
      dimension iqnums( * ),       ! SU(3)xU(M) quantum numbers         ACVI0779
     1          kets( nmax, * )    ! many-body bit states of nwords elm.ACVI0780
C               nkets              ! number of bit states generated     ACVI0781
C               mxkets             ! max dimension of kets              ACVI0782
C     common variables from subroutine hws.                             ACVI0783
      common / HWSCON / lfile, ieta, nwords, nbits, mbits               ACVI0784
C     U(3) counters for U(N) distribution and internal counters for U(M)ACVI0785
      parameter ( NSTSMX = 66,     ! max number of levels for a shell   ACVI0786
     1            MXBITS = 32,     ! number of bits in one word         ACVI0787
     2            MXCONS = 1024 )  ! max U(M) configurations = 2**10    ACVI0788
      dimension mdist( NSTSMX ),   ! f symmetry, Young box              ACVI0789
     1          jbot( NSTSMX ),    ! lower                              ACVI0790
     2          jconun( NSTSMX ),  ! conf. counter of U(M) in a U(3) IR ACVI0791
     3          jtop( NSTSMX )     ! upper                              ACVI0792
      common / SUNQN /                                                  ACVI0793
     1         iconu3( NSTSMX ),   ! U(3) labels for all states         ACVI0794
     2         iptrmp( 12 ),       ! pointer for number of particles    ACVI0795
     3         nsun( MXCONS ),     ! binary rep. of U(M) conf.          ACVI0796
     4         isun( MXCONS )      ! packed rep. of U(M) conf.          ACVI0797
C     Packing Functions == 3 labels                                     ACVI0798
      ipack(i1,i2,i3) = ior(i3,ishft(ior(i2,ishft(i1,10)),10))          ACVI0799
C                                                                       ACVI0800
C     Generate First distribution (and reducing the unnecessary ones.)  ACVI0801
      maxket = mxkets / nmax                                            ACVI0802
      tfile  = lfile.gt.0                                               ACVI0803
      nsts   = (ieta+1) * (ieta+2) / 2                                  ACVI0804
      nmaxp2 = nmax + 2                                                 ACVI0805
      nkets  = 0                                                        ACVI0806
      ketsmx = 0                                                        ACVI0807
      np = 0                                                            ACVI0808
      do I = 3, nmaxp2                                                  ACVI0809
        np = np + iqnums( i )                                           ACVI0810
      ENDDO                                                             ACVI0811
C                                                                       ACVI0812
C     Generate initial U(3) distribution for Np particles.              ACVI0813
      ny = (ieta*np - (iqnums(1) + 2*iqnums(2)) ) / 3                   ACVI0814
      nx = ny + iqnums(2)                                               ACVI0815
      nz = nx + iqnums(1)                                               ACVI0816
      nzxy = ipack(nz,nx,ny)                                            ACVI0817
C     nzxy = ipack(ny+iqnums(2)+iqnums(1),ny+iqnums(2),ny)              ACVI0818
CW----------------------------------------------------------------------ACVI0819
CW    Writing the projection quantum numbers.                           ACVI0820
      IF(TFILE) THEN                                                    ACVI0821
      WRITE(LFILE,1000)                                                 ACVI0822
 1000 FORMAT(/' ',T8,'NP',T12,'eps',T16,'MLamT')                        ACVI0823
      WRITE(LFILE,1020) NP,2*NZ-NX-NY,NX-NY                             ACVI0824
      WRITE(LFILE,1010) (' ',I,I=1,NMAX)                                ACVI0825
 1010 FORMAT(' ',T8,'Nz',T13,'Nx',T18,'Ny',T22,10(A,'f',I1,'~',1X))     ACVI0826
      WRITE(LFILE,1020) NZ,NX,NY,(IQNUMS(I),I=3,NMAXP2)                 ACVI0827
 1020 FORMAT(' ',3X,13I5)                                               ACVI0828
      ENDIF                                                             ACVI0829
CW----------------------------------------------------------------------ACVI0830
C                                                                       ACVI0831
C     Generate initial U(M) distribution for Np particles.              ACVI0832
      NF=IQNUMS(3)                                                      ACVI0833
      DO I=4,NMAXP2                                                     ACVI0834
        NF=IOR(ISHFT(NF,NBITS),IQNUMS(I))                               ACVI0835
      END DO                                                            ACVI0836
C     Generate U(N) space Highest Weight State.                         ACVI0837
      DO I=1,NSTS             ! check SUBROUTINE UNXUM                  ACVI0838
        MDIST(I)=0                                                      ACVI0839
      END DO                                                            ACVI0840
      DO  J=3,NMAXP2                                                    ACVI0841
        DO I=1,IQNUMS(J)                                                ACVI0842
          MDIST(I)=MDIST(I)+1                                           ACVI0843
        END DO                                                          ACVI0844
      END DO                                                            ACVI0845
C     MORE=1                                                            ACVI0846
C                                                                       ACVI0847
C     Match weights for U(3) distribution of Np particles.              ACVI0848
C     with input values                                                 ACVI0849
      NSTSM1=NSTS-1                                                     ACVI0850
100   ISUM=0                                                            ACVI0851
        DO I=1,NSTS                                                     ACVI0852
          ISUM=ISUM+MDIST(I)*ICONU3(I)                                  ACVI0853
        END DO                                                          ACVI0854
C       Generate allowed U(3) distribution for Np particles.            ACVI0855
        IF(ISUM.EQ.NZXY)THEN                                            ACVI0856
C         Bounds and counter for state in U(M) distribution.            ACVI0857
          do i = 1, nsts                                                ACVI0858
            jbot( i )   = iptrmp( mdist(i)+1 )                          ACVI0859
            jconun( i ) = jbot( i )                                     ACVI0860
            jtop( i )   = iptrmp( mdist(i)+2 ) - 1                      ACVI0861
          end do                                                        ACVI0862
110       JCONVC=0                                                      ACVI0863
C           Label state according to U(M).                              ACVI0864
            DO I=1,NSTS                                                 ACVI0865
              JCONVC=JCONVC+ISUN(JCONUN(I))                             ACVI0866
            END DO                                                      ACVI0867
C                                                                       ACVI0868
C           Match U(M) IKETs with input values.                         ACVI0869
C           Generate allowed state.                                     ACVI0870
            IF(JCONVC.EQ.NF)THEN                                        ACVI0871
              NKETS=NKETS+1                                             ACVI0872
              KETSMX=KETSMX+1                                           ACVI0873
              if( ketsmx .gt. maxket )                                  ACVI0874
     1          call attn(' GENMPB: Overflow in array KETS',*999)       ACVI0875
              DO IUN=1,NMAX                                             ACVI0876
                K=0                                                     ACVI0877
                IFLAG=0                                                 ACVI0878
                DO J=0,NSTSM1                                           ACVI0879
                  I=MOD(J,MXBITS)                                       ACVI0880
                  IF(IFLAG.EQ.0) KOLD=KETSMX                            ACVI0881
                  IF(BTEST(NSUN(JCONUN(J+1)),IUN-1)) K=IBSET(K,I)       ACVI0882
                  IF(I.EQ.0 .AND. J.NE.0)THEN                           ACVI0883
                    KETS(IUN,KETSMX)=K                                  ACVI0884
                    IFLAG=1                                             ACVI0885
                    KETSMX=KETSMX+1                                     ACVI0886
                    K=0                                                 ACVI0887
                  END IF                                                ACVI0888
                END DO                                                  ACVI0889
                KETS(IUN,KETSMX)=K                                      ACVI0890
                IF(IUN.NE.NMAX) KETSMX=KOLD                             ACVI0891
              END DO                                                    ACVI0892
            END IF                                                      ACVI0893
C           Gett new U(M) distribution if possibilities not exhausted   ACVI0894
C           for U(3) irrep IU3.                                         ACVI0895
            DO IU3=NSTS,1,-1                                            ACVI0896
              JCONUN(IU3)=JCONUN(IU3)+1                                 ACVI0897
              IF(JCONUN(IU3).LE.JTOP(IU3)) GOTO 110                     ACVI0898
              JCONUN(IU3)=JBOT(IU3)                                     ACVI0899
            END DO                                                      ACVI0900
        ENDIF                                                           ACVI0901
C                                                                       ACVI0902
C       Search new U(3) distribution for Np particles.                  ACVI0903
        DO J=NSTSM1,1,-1                                                ACVI0904
          IF(MDIST(J).GE.1)THEN                                         ACVI0905
            MDIST(J)=MDIST(J)-1                                         ACVI0906
            NTOT=0                                                      ACVI0907
            DO K=1,J                                                    ACVI0908
              NTOT=NTOT+MDIST(K)                                        ACVI0909
            END DO                                                      ACVI0910
            DO K=J+1,NSTS                                               ACVI0911
              MDIST(K)=MIN0(NP-NTOT,NMAX)                               ACVI0912
              NTOT=NTOT+MDIST(K)                                        ACVI0913
            END DO                                                      ACVI0914
            IF(NTOT.EQ.NP)THEN                                          ACVI0915
              GOTO 100                                                  ACVI0916
            ELSE                                                        ACVI0917
              MDIST(J)=MDIST(J)+1                                       ACVI0918
            END IF                                                      ACVI0919
          END IF                                                        ACVI0920
        END DO                                                          ACVI0921
C       MORE=0                                                          ACVI0922
C       PRINT *,'GENMPB: All patterns exhausted'                        ACVI0923
CW----------------------------------------------------------------------ACVI0924
      if( nkets .eq. 0 )                                                ACVI0925
     1  write(6,*) ' NO STATES EXIST WITH THESE QUANTUM NUMBERS'        ACVI0926
CW----------------------------------------------------------------------ACVI0927
C                                                                       ACVI0928
      return                                                            ACVI0929
999   return 1                                                          ACVI0930
C ---*end of genmpb*----------------------------------------------------ACVI0931
      end                                                               ACVI0932
CB----------------------------------------------------------------------ACVI0933
C                                                                       ACVI0934
C                          **************                               ACVI0935
C                          *** APPAIJ ***                               ACVI0936
C                          **************                               ACVI0937
C                                                                       ACVI0938
C ----------------------------------------------------------------------ACVI0939
C Author:  Chairul Bahri                                                ACVI0940
C ----------------------------------------------------------------------ACVI0941
C Updates: 07/90 ==> Rewrite of the original (Reske/Lopes'z version.)   ACVI0942
C          A major change: spacelike bit representation so it can apply ACVI0943
C          for SU(M) internal symmetry.                                 ACVI0944
C          08/18/90: running successfully for HWSGEN !!!                ACVI0945
C          11/11/90: implemented in RMECODE.                            ACVI0946
C          06/25/91: compress for non-zero coefficients of basis only.  ACVI0947
C          11/15/91: SUN 4.                                             ACVI0948
C ----------------------------------------------------------------------ACVI0949
C                                                                       ACVI0950
C    The subroutine APPAIJ generates the highest weight (HWS) for       ACVI0951
C    a given number number of particles in a shell IETA after           ACVI0952
C    applying the action of U(3) (Czx,Cxy) and U(M) raising operators   ACVI0953
C    to all possible bit states generated with GENMPB, and then setting ACVI0954
C    up a homogeneous linear system of equations for the coefficients   ACVI0955
C    in the HWS expansion, which is solved by using SPARSE.             ACVI0956
C                                                                       ACVI0957
C ----------------------------------------------------------------------ACVI0958
C                                                                       ACVI0959
      subroutine appaij( nmax, iqnums, kets, dsolmat, nbkets, nsol,     ACVI0960
     ]                   mxnsol, tcmprs, * )                            ACVI0961
      implicit real*8(d),logical(t)                                     ACVI0962
      dimension iqnums( * ),       ! U(M) f~ quanum numbers             ACVI0963
     1          kets(nmax, *),     ! many-body bit states               ACVI0964
     2          dsolmat( * )       ! solution matrix for hws vectors    ACVI0965
C               nmax               ! internal symmetry multiplicity     ACVI0966
C               nbkets             ! # bit states, kets                 ACVI0967
C               nsol               ! # independent solution =           ACVI0968
C                                  !   U(N) > U(3) multiplicity         ACVI0969
C               mxnsol             ! max # solutions to linear system   ACVI0970
C                                  !   of equations                     ACVI0971
C               tcmprs             ! .true.                             ACVI0972
C                                  !   if keep non-zero coefficients    ACVI0973
      common / HWSCON / lfile,     ! logfile                            ACVI0974
     1                  ieta,      ! oscillator shell number            ACVI0975
     2                  nwords,    ! # words in one bit state           ACVI0976
     3                  nbits,     ! bit length for # levels in shell   ACVI0977
     4                  mbits      ! length of bits in a word           ACVI0978
C                                                                       ACVI0979
      parameter ( NSTSMX = 66,     ! max # levels for a given shell     ACVI0980
     1            MXBITS = 32,     ! # bits in a word                   ACVI0981
     2            MASK = MXBITS-1, ! significant bit                    ACVI0982
     3            MXCONS = 1024 )  ! max # U(M) configurations = 2**10  ACVI0983
C                                                                       ACVI0984
C     SU(3) raising operator coefficients.                              ACVI0985
      common / SU3COF /                                                 ACVI0986
     1         dcoefzx( NSTSMX ),  ! Czx coefficients                   ACVI0987
     2         dcoefxy( NSTSMX ),  ! Cxy coefficients                   ACVI0988
     3         nmlvls( 3 )         ! size of significant bits in word   ACVI0989
C                                                                       ACVI0990
C     Tree for bit patterns.                                            ACVI0991
      parameter ( MXACTS = 100000, ! max dim. of distinct bit patterns  ACVI0992
     1            MXTREE = 140000 )! tree size                          ACVI0993
      common / PATTRN /                                                 ACVI0994
     1         dacts( MXACTS ),    ! buffer linked-list for coefficientsACVI0995
     2         iacts( MXACTS ),    ! linked-list index for coefficients ACVI0996
     3         irsptr( MXACTS ),   ! signature for raising operators    ACVI0997
     4         ketu3( MXACTS,3 ),  ! bit pattern as key                 ACVI0998
     5         ketree(-10:MXTREE), ! binary tree for bit patterns       ACVI0999
     6         iactmx              ! last non-empty index               ACVI1000
      dimension ketun( mxacts )                                         ACVI1001
      equivalence ( ketun( 1 ), ketu3( 1,1 ) )                          ACVI1002
C     data izx,ixy /za0000000,zb0000000/! Pointer for raising operators ACVI1003
      data izx,ixy /-1610612736,-1342177280/                            ACVI1004
C                                                                       ACVI1005
C     Tree for coefficients bit patterns after raising operations.      ACVI1006
      parameter ( MXIEQN = 100000, ! max # equations                    ACVI1007
     1            MXIDIS = 100000, ! max # linearly independent eqns    ACVI1008
     2            MXICOE = 100000, ! 2 * max # non-zero coefficients in ACVI1009
     3                             !   the row reduced system of eqns.  ACVI1010
     4            MXIVEC = 10000 ) ! max # non-zero elements allowed forACVI1011
C                                  !   row vector                       ACVI1012
      common / COEFFS /                                                 ACVI1013
     1         deqns( MXIEQN ),    ! buffer linked-list for coefficientsACVI1014
     2         ieqns( 2,MXIEQN ),  ! linked-list index                  ACVI1015
     3         iqtree( -10:MXTREE )! tree for coefficients              ACVI1016
      dimension ketx( 30 ),        ! key                                ACVI1017
     1          ketemp( 6 ),       ! temporary variables for ketx       ACVI1018
     2          ntop( 9 )          ! f~ symmetry -> phase               ACVI1019
C                                                                       ACVI1020
C     Arrays for SPARSE/ROWRED subroutines.                             ACVI1021
      common / ISPARS / ivec(MXIVEC), jvec(MXIVEC), idis(MXIDIS)        ACVI1022
     *                 ,locptr(MXIDIS), icoef(MXICOE)                   ACVI1023
      common / DSPARS / dvec(MXIVEC), dxvec(MXIVEC), dcoef(MXICOE)      ACVI1024
C                                                                       ACVI1025
      tsmfi  = lfile.gt.0          ! small log file                     ACVI1026
      tfile  = lfile.gt.6          ! extended log file                  ACVI1027
      nsts   = mod( (ieta+1)*(ieta+2)/2, MXBITS )                       ACVI1028
C                                                                       ACVI1029
C     Initialize arrays for subroutine ROWRED/SOLVE.                    ACVI1030
      call sparse( MXIVEC, MXICOE, MXIDIS )                             ACVI1031
      idismx = 0                   ! counters for elements in row       ACVI1032
      icoemx = 0                   !                                    ACVI1033
      mxxcoe = mxicoe / 2          ! counters for rows                  ACVI1034
      ihalf  = 1                   !                                    ACVI1035
      jhalf  = mxxcoe + 1          !                                    ACVI1036
      mxxdis = mxidis / 2          !                                    ACVI1037
      khalf  = 1                   !                                    ACVI1038
      lhalf  = mxxdis + 1          !                                    ACVI1039
C                                                                       ACVI1040
C     Initialize tree for coefficients                                  ACVI1041
      nfinal = nmax * nwords                                            ACVI1042
      npairs = 2 * nwords                                               ACVI1043
      if( tfile ) then                                                  ACVI1044
        nfnlm  = nfinal - nmax                                          ACVI1045
        nfnlp1 = nfnlm + 1                                              ACVI1046
      end if                                                            ACVI1047
      mxnodes = MXTREE / (nfinal + 5)                                   ACVI1048
      call tsetlf( iqtree, mxnodes, nfinal, 2 )                         ACVI1049
C                                                                       ACVI1050
C***********************************************************************ACVI1051
C                                                                       ACVI1052
C     APPLYING RAISING OPERATORS Czx,Cxy,Aij                            ACVI1053
C                                                                       ACVI1054
C***********************************************************************ACVI1055
C     WRITE(6,'(//A//)') ' *** ACTION OF RAISING OPERATORS ***'         ACVI1056
C                                                                       ACVI1057
      iqtmx = 0                                                         ACVI1058
      ieqnmx = 0                                                        ACVI1059
C                                                                       ACVI1060
      ketdis = 0                                                        ACVI1061
C                                                                       ACVI1062
C     Loop through all basis kets.                                      ACVI1063
      do 470 iket = 1, nbkets                                           ACVI1064
C       Store each basis ket in KETX for I/O and IQTREE.                ACVI1065
        k=0                                                             ACVI1066
        do 200 j = 1, nwords                                            ACVI1067
          do 200 i = 1, nmax                                            ACVI1068
            k = k + 1                                                   ACVI1069
200         ketx( k )=kets( i, ketdis+j )                               ACVI1070
CW----------------------------------------------------------------------ACVI1071
CW    WRITING STATE OVER WHICH THE ACTION WILL BE CONSIDERED.           ACVI1072
      IF(TFILE) THEN                                                    ACVI1073
      WRITE(LFILE,1000) IKET                                            ACVI1074
      DO 210 IDUS=1,NFNLM                                               ACVI1075
        CALL DBOVRL(KETX(IDUS))                                         ACVI1076
  210   IF(MOD(IDUS,NMAX).EQ.0) WRITE(LFILE,*)                          ACVI1077
      DO 220 IDUS=NFNLP1,NFINAL                                         ACVI1078
  220   CALL DBOVRC(NSTS,KETX(IDUS))                                    ACVI1079
      WRITE(LFILE,'(/A)') ' action of raising operators:'               ACVI1080
      ENDIF                                                             ACVI1081
CW----------------------------------------------------------------------ACVI1082
C       ACTION !!!                                                      ACVI1083
        DO 370 IUN=1,NMAX                                               ACVI1084
C         Store each basis ket in KETX for manipulation.                ACVI1085
          DO 230 J=1,NWORDS                                             ACVI1086
            KETEMP(J)=KETX(J)                                           ACVI1087
230         KETX(J)=KETS(IUN,KETDIS+J)                                  ACVI1088
          DO 240 J=NWORDS+1,NPAIRS                                      ACVI1089
            KETEMP(J)=KETX(J)                                           ACVI1090
240         KETX(J)=0                                                   ACVI1091
C         Check the basis in KETREE.                                    ACVI1092
CB 08/16/91 10:45 ***>                                                  ACVI1093
          CALL TCHK(KETX,KETREE,*290)                                   ACVI1094
          ILOC = KETREE(-5)                   ! new node                ACVI1095
          KETREE( ILOC+KETREE(-2) ) = 1       ! equal priority          ACVI1096
          ILOC = ILOC+KETREE(-4) + 1          ! first data              ACVI1097
CB <*** 08/16/91 17:48                                                  ACVI1098
C         Counters for creation and annihilation places.                ACVI1099
          NZ=MXBITS                                                     ACVI1100
          NX=MXBITS+1                                                   ACVI1101
          NY=MXBITS+2                                                   ACVI1102
C         Action of Czx & Cxy on first U(3) irrep yields zero.          ACVI1103
          IMLVLU3=1                                                     ACVI1104
C         Keeping track of states in between for change of phase.       ACVI1105
          NBTWN=0                                                       ACVI1106
          IF(BTEST(KETX(1),0)) NBTWN=1                                  ACVI1107
C         Run through all U(3) irreps in KETX.                          ACVI1108
          DO 280 I=1,IETA                                               ACVI1109
            DO 270 J=1,I                                                ACVI1110
              IMLVLU3=IMLVLU3+1                                         ACVI1111
C                                                                       ACVI1112
C             Counter for word in ket.                                  ACVI1113
              IZWPTR=ISHFT(NZ,MBITS)                                    ACVI1114
              IXWPTR=ISHFT(NX,MBITS)                                    ACVI1115
              IYWPTR=ISHFT(NY,MBITS)                                    ACVI1116
C                                                                       ACVI1117
              IZWORD=KETX(IZWPTR)                                       ACVI1118
              IXWORD=KETX(IXWPTR)                                       ACVI1119
              IYWORD=KETX(IYWPTR)                                       ACVI1120
C                                                                       ACVI1121
C             Counter for bit in word.                                  ACVI1122
              IZMASK=IAND(NZ,MASK)                                      ACVI1123
              IXMASK=IAND(NX,MASK)                                      ACVI1124
              IYMASK=IAND(NY,MASK)                                      ACVI1125
C                                                                       ACVI1126
C             Test occupances.                                          ACVI1127
              IF(BTEST(IXWORD,IXMASK))THEN                              ACVI1128
C                                                                       ACVI1129
C               ACTION OF Czx                                           ACVI1130
C                                                                       ACVI1131
                IF(.NOT.BTEST(IZWORD,IZMASK))THEN                       ACVI1132
C                 Keep track of the phase.                              ACVI1133
                  NBTWN=NBTWN+1                                         ACVI1134
C                 Destroy!                                              ACVI1135
                  KETX(IXWPTR)=IBCLR(KETX(IXWPTR),IXMASK)               ACVI1136
C                 Create!                                               ACVI1137
                  KETX(IZWPTR)=IBSET(KETX(IZWPTR),IZMASK)               ACVI1138
                  IACTMX=IACTMX+1                                       ACVI1139
                  if( iactmx .gt. mxacts )                              ACVI1140
     1              call attn(' APPAIJ: Increase MXACTS!',*999)         ACVI1141
C                 IF(IACTMX.GT.MXACTS) CALL ERROR('Increase MXACTS!')   ACVI1142
                  IF(BTEST(NBTWN,0))THEN                                ACVI1143
                    DACTS(IACTMX)=DCOEFZX(IMLVLU3)                      ACVI1144
                  ELSE                                                  ACVI1145
                    DACTS(IACTMX)=-DCOEFZX(IMLVLU3)                     ACVI1146
                  ENDIF                                                 ACVI1147
                  DO 250 K=1,NWORDS                                     ACVI1148
250                 KETU3(IACTMX,K)=KETX(K)                             ACVI1149
                  IRSPTR(IACTMX)=IZX     ! Czx pointer.                 ACVI1150
C                 Begin insertion in KETREE for the obtained kets.      ACVI1151
                  IACTS(IACTMX)=KETREE(ILOC+1)                          ACVI1152
                  KETREE(ILOC)=KETREE(ILOC)+1                           ACVI1153
                  KETREE(ILOC+1)=IACTMX                                 ACVI1154
C                 Restore values after insertion in tree.               ACVI1155
                  KETX(IXWPTR)=IXWORD                                   ACVI1156
                  KETX(IZWPTR)=IZWORD                                   ACVI1157
                ENDIF                                                   ACVI1158
C                                                                       ACVI1159
              ELSE                                                      ACVI1160
                IF(BTEST(IZWORD,IZMASK)) NBTWN=NBTWN-1                  ACVI1161
C                                                                       ACVI1162
C               ACTION OF Cxy                                           ACVI1163
C                                                                       ACVI1164
                IF(BTEST(IYWORD,IYMASK))THEN                            ACVI1165
C                 Destroy!                                              ACVI1166
                  KETX(IYWPTR)=IBCLR(KETX(IYWPTR),IYMASK)               ACVI1167
C                 Create!                                               ACVI1168
                  KETX(IXWPTR)=IBSET(KETX(IXWPTR),IXMASK)               ACVI1169
                  IACTMX=IACTMX+1                                       ACVI1170
                  IF(IACTMX.GT.MXACTS)                                  ACVI1171
     1              call attn(' APPAIJ: Increase MXACTS!',*999)         ACVI1172
C    1                CALL ERROR(' APPAIJ: Increase MXACTS!')           ACVI1173
                  DACTS(IACTMX)=DCOEFXY(IMLVLU3+1)                      ACVI1174
                  DO 260 K=1,NWORDS                                     ACVI1175
260                 KETU3(IACTMX,K)=KETX(K)                             ACVI1176
                  IRSPTR(IACTMX)=IXY    ! Cxy pointer                   ACVI1177
C                 Begin insertion in KETREE for the obtained kets.      ACVI1178
                  IACTS(IACTMX)=KETREE(ILOC+1)                          ACVI1179
                  KETREE(ILOC)=KETREE(ILOC)+1                           ACVI1180
                  KETREE(ILOC+1)=IACTMX                                 ACVI1181
C                 Restore values after insertion in tree. .             ACVI1182
                  KETX(IYWPTR)=IYWORD                                   ACVI1183
                  KETX(IXWPTR)=IXWORD                                   ACVI1184
                ENDIF                                                   ACVI1185
              ENDIF                                                     ACVI1186
C             Skip to next U(3) irreps.                                 ACVI1187
              NZ=NZ+1                                                   ACVI1188
              NX=NX+1                                                   ACVI1189
              NY=NY+1                                                   ACVI1190
270         CONTINUE                                                    ACVI1191
C           Keep track of possible change in phase.                     ACVI1192
            IXWORD=KETX(ISHFT(NX,MBITS))                                ACVI1193
            IXMASK=IAND(NX,MASK)                                        ACVI1194
            IF(BTEST(IXWORD,IXMASK)) NBTWN=NBTWN+1                      ACVI1195
C           Change counter for site of possible action.                 ACVI1196
            IMLVLU3=IMLVLU3+1                                           ACVI1197
            NX=NX+1                                                     ACVI1198
            NY=NY+1                                                     ACVI1199
280       CONTINUE                                                      ACVI1200
C                                                                       ACVI1201
CB 08/16/91 10:45 ***>                                                  ACVI1202
          CALL TINS(KETX,KETREE)                                        ACVI1203
          CALL TCHK(KETX,KETREE,*290)                                   ACVI1204
290       ILOC=KETREE(-5)+KETREE(-4)+1                                  ACVI1205
CB <*** 08/17/91 17:47                                                  ACVI1206
C         Restore values after manipulation.                            ACVI1207
          DO 300 J=1,NPAIRS                                             ACVI1208
300         KETX(J)=KETEMP(J)                                           ACVI1209
C         Load all obtained kets for IQTREE insertion.                  ACVI1210
          IVECMX=KETREE(ILOC)                                           ACVI1211
          INEXT=KETREE(ILOC+1)                                          ACVI1212
          DO 360 I=IVECMX,1,-1                                          ACVI1213
            IEQNMX=IEQNMX+1                                             ACVI1214
            IF(IEQNMX.GT.MXIEQN)                                        ACVI1215
     1      call attn(' APPAIJ: Too many eqns. Increase MXIEQN!',*999)  ACVI1216
C    1      CALL ERROR(' APPAIJ: Too many equations. Increase MXIEQN!') ACVI1217
            DEQNS(IEQNMX)=DACTS(INEXT)                                  ACVI1218
            IEQNS(1,IEQNMX)=IKET                                        ACVI1219
C           Load kets.                                                  ACVI1220
            K=IUN                                                       ACVI1221
            DO 310 L=1,NWORDS                                           ACVI1222
              KETEMP(L)=KETX(K)                                         ACVI1223
              KETX(K)=KETU3(INEXT,L)                                    ACVI1224
310         K=K+NMAX                                                    ACVI1225
CW----------------------------------------------------------------------ACVI1226
CW    WRITING COEFFICIENT AND STATE OBTAINED.                           ACVI1227
      IF(TFILE) THEN                                                    ACVI1228
        IF(IRSPTR(INEXT).EQ.IZX) THEN                                   ACVI1229
          WRITE(LFILE,'(A)') '   Czx |'                                 ACVI1230
        ELSE IF(IRSPTR(INEXT).EQ.IXY) THEN                              ACVI1231
          WRITE(LFILE,'(A)') '   Cxy |'                                 ACVI1232
        ENDIF                                                           ACVI1233
      WRITE(LFILE,1100) DEQNS(IEQNMX)                                   ACVI1234
      DO 320 IDUS=1,NFNLM                                               ACVI1235
        CALL DBOVRL(KETX(IDUS))                                         ACVI1236
  320   IF(MOD(IDUS,NMAX).EQ.0) WRITE(LFILE,*)                          ACVI1237
      DO 330 IDUS=NFNLP1,NFINAL                                         ACVI1238
  330   CALL DBOVRC(NSTS,KETX(IDUS))                                    ACVI1239
      ENDIF                                                             ACVI1240
CW----------------------------------------------------------------------ACVI1241
C           Putting the pointer in the last KETX for IQTREE search.     ACVI1242
            KETX(NFINAL)=IEOR(KETX(NFINAL),IRSPTR(INEXT))               ACVI1243
C           Begin insertion in IQTREE for setting system of equations   ACVI1244
C           naming unknown for coefficient.                             ACVI1245
CB 08/16/91 11:30 ***>                                                  ACVI1246
            CALL TCHK(KETX,IQTREE,*340)                                 ACVI1247
C           Counter for different nodes in IQTREE.                      ACVI1248
            IQTMX = IQTMX + 1                                           ACVI1249
            ILOC = IQTREE(-5)                   ! new node              ACVI1250
            IQTREE( ILOC+IQTREE(-2) ) = 1       ! equal priority        ACVI1251
            CALL TINS(KETX,IQTREE)                                      ACVI1252
            CALL TCHK(KETX,IQTREE,*340)                                 ACVI1253
340         ILOC = IQTREE(-5) + IQTREE(-4) + 1  ! first data            ACVI1254
CB <*** 08/16/91 17:48                                                  ACVI1255
C           Pointing to previous coefficient in actual equation and     ACVI1256
C           recording the occurence of this resulting ket.              ACVI1257
            IEQNS(2,IEQNMX)=IQTREE(ILOC+1)                              ACVI1258
            IQTREE(ILOC)=IQTREE(ILOC)+1                                 ACVI1259
            IQTREE(ILOC+1)=IEQNMX                                       ACVI1260
C           Restoring values after insertion in IQTREE.                 ACVI1261
            KETX(NFINAL)=IEOR(KETX(NFINAL),IRSPTR(INEXT))               ACVI1262
            K=IUN                                                       ACVI1263
            DO 350 L=1,NWORDS                                           ACVI1264
              KETX(K)=KETEMP(L)                                         ACVI1265
350         K=K+NMAX                                                    ACVI1266
360       INEXT=IACTS(INEXT)                                            ACVI1267
370     CONTINUE                                                        ACVI1268
C       End application of U(3) raising operators.                      ACVI1269
C                                                                       ACVI1270
C               ACTION OF Aij                                           ACVI1271
C                                                                       ACVI1272
        K=0                                                             ACVI1273
        DO IUN=1,NMAX-1                                                 ACVI1274
          NTOP(IUN)=IQNUMS(IUN)                                         ACVI1275
        ENDDO                                                           ACVI1276
C       Running through the number of words in ket.                     ACVI1277
        DO 460 J=1,NWORDS                                               ACVI1278
C         Skipping the first.                                           ACVI1279
          K=K+1                                                         ACVI1280
          DO 450 IUN=2,NMAX                                             ACVI1281
            K=K+1                                                       ACVI1282
            IUNM1=IUN-1                                                 ACVI1283
            IPWORD=IAND(IEOR(KETX(K-1),KETX(K)),KETX(K))                ACVI1284
            IF(IPWORD .EQ. 0) GOTO 450                                  ACVI1285
C           Storing a pair of kets in KETEMP for KETREE.                ACVI1286
            KTEMP=IUN                                                   ACVI1287
            DO 380 JTEMP=1,NWORDS                                       ACVI1288
              KETEMP(JTEMP)=KETX(KTEMP-1)                               ACVI1289
              KETEMP(JTEMP+NWORDS)=KETX(KTEMP)                          ACVI1290
380         KTEMP=KTEMP+NMAX                                            ACVI1291
            KETEMP(NWORDS)=IOR(KETEMP(NWORDS),ISHFT(J,28))              ACVI1292
C           Checking the action in KETREE.                              ACVI1293
CB 08/16/91 10:51 ***>                                                  ACVI1294
            CALL TCHK(KETEMP,KETREE,*400)                               ACVI1295
            ILOC = KETREE(-5)                 ! new node                ACVI1296
            KETREE( ILOC+KETREE(-2) ) = 1     ! equal priority          ACVI1297
            ILOC = ILOC+KETREE(-4) + 1        ! first data              ACVI1298
CB <*** 08/16/91 17:48                                                  ACVI1299
C           Keeping track of states in between for change of phase.     ACVI1300
            NBTWN=NTOP(IUNM1)                                           ACVI1301
C           Loop through number of U(3) irreps.                         ACVI1302
            DO 390 I=0,NMLVLS(J)                                        ACVI1303
              IF(BTEST(KETX(K-1),I)) THEN                               ACVI1304
                NBTWN=NBTWN-1                                           ACVI1305
                NTOP(IUNM1)=NTOP(IUNM1)-1                               ACVI1306
              ENDIF                                                     ACVI1307
              IF(BTEST(IPWORD,I))THEN                                   ACVI1308
                IACTMX=IACTMX+1                                         ACVI1309
                IF(IACTMX.GT.MXACTS)                                    ACVI1310
     1              call attn(' APPAIJ: Increase MXACTS!',*999)         ACVI1311
C    1              CALL ERROR(' APPAIJ: Increase MXACTS!')             ACVI1312
                IF(BTEST(NBTWN,0))THEN                                  ACVI1313
                  DACTS(IACTMX)=-1.0D0                                  ACVI1314
                ELSE                                                    ACVI1315
                  DACTS(IACTMX)=1.0D0                                   ACVI1316
                ENDIF                                                   ACVI1317
                KETUN(IACTMX)=IBSET(0,I)                                ACVI1318
CB              KETU3(IACTMX,1)=IBSET(0,I)                              ACVI1319
C               Putting Aij pointer.                                    ACVI1320
                IRSPTR(IACTMX)=ISHFT(IUNM1,28)                          ACVI1321
C               Begin insertion in KETREE for the action.               ACVI1322
                IACTS(IACTMX)=KETREE(ILOC+1)                            ACVI1323
                KETREE(ILOC)=KETREE(ILOC)+1                             ACVI1324
                KETREE(ILOC+1)=IACTMX                                   ACVI1325
              ENDIF                                                     ACVI1326
              IF(BTEST(KETX(K),I)) NBTWN=NBTWN+1                        ACVI1327
390         CONTINUE                                                    ACVI1328
C                                                                       ACVI1329
CB 08/16/91 10:54 ***>                                                  ACVI1330
            CALL TINS(KETEMP,KETREE)                                    ACVI1331
            CALL TCHK(KETEMP,KETREE,*400)                               ACVI1332
400         ILOC=KETREE(-5)+KETREE(-4)+1                                ACVI1333
CB <*** 08/16/91 17:49                                                  ACVI1334
C           Loading all obtained kets for IQTREE insertion.             ACVI1335
            IVECMX=KETREE(ILOC)                                         ACVI1336
            INEXT=KETREE(ILOC+1)                                        ACVI1337
            DO 440 I=IVECMX,1,-1                                        ACVI1338
              IEQNMX=IEQNMX+1                                           ACVI1339
              IF(IEQNMX.GT.MXIEQN)                                      ACVI1340
     1          call attn(' APPAIJ: Increase MXIEQN!',*999)             ACVI1341
C             IF(IEQNMX.GT.MXIEQN) CALL ERROR(' APPAIJ: Too many eqn.') ACVI1342
              DEQNS(IEQNMX)=DACTS(INEXT)                                ACVI1343
              IEQNS(1,IEQNMX)=IKET                                      ACVI1344
C             Loading kets.                                             ACVI1345
              KETX(K-1)=IEOR(KETX(K-1),KETUN(INEXT))                    ACVI1346
              KETX(K)=IEOR(KETX(K),KETUN(INEXT))                        ACVI1347
CW----------------------------------------------------------------------ACVI1348
CW    WRITING COEFFICIENT AND STATE OBTAINED.                           ACVI1349
      IF(TFILE) THEN                                                    ACVI1350
      WRITE(LFILE,'(A,2I1,A)') '   C',IUNM1,IUN,' |'                    ACVI1351
      WRITE(LFILE,1100) DEQNS(IEQNMX)                                   ACVI1352
      DO 410 IDUS=1,NFNLM                                               ACVI1353
        CALL DBOVRL(KETX(IDUS))                                         ACVI1354
  410   IF(MOD(IDUS,NMAX).EQ.0) WRITE(LFILE,*)                          ACVI1355
      DO 420 IDUS=NFNLP1,NFINAL                                         ACVI1356
  420   CALL DBOVRC(NSTS,KETX(IDUS))                                    ACVI1357
      ENDIF                                                             ACVI1358
CW----------------------------------------------------------------------ACVI1359
C             Putting the pointer in the last KETX for IQTREE search.   ACVI1360
              KETX(NFINAL)=IEOR(KETX(NFINAL),IRSPTR(INEXT))             ACVI1361
C             Begin insertion in IQTREE for setting system of           ACVI1362
C             equations naming unknown for coefficient.                 ACVI1363
CB 08/16/91 11:40 ***>                                                  ACVI1364
              CALL TCHK(KETX,IQTREE,*430)                               ACVI1365
C             Counter for different nodes in IQTREE.                    ACVI1366
              IQTMX = IQTMX + 1                                         ACVI1367
              ILOC = IQTREE(-5)                 ! new node              ACVI1368
              IQTREE( ILOC+IQTREE(-2) ) = 1     ! equal priority        ACVI1369
              CALL TINS(KETX,IQTREE)                                    ACVI1370
              CALL TCHK(KETX,IQTREE,*430)                               ACVI1371
430           ILOC = IQTREE(-5) + IQTREE(-4) + 1! first data            ACVI1372
CB <*** 08/16/91 17:49                                                  ACVI1373
C             Pointing to previous coefficient in actual equation       ACVI1374
C             and recording the occurence of this resulting ket.        ACVI1375
              IEQNS(2,IEQNMX)=IQTREE(ILOC+1)                            ACVI1376
              IQTREE(ILOC)=IQTREE(ILOC)+1                               ACVI1377
              IQTREE(ILOC+1)=IEQNMX                                     ACVI1378
C             Restoring values after insertion in IQTREE.               ACVI1379
              KETX(NFINAL)=IEOR(KETX(NFINAL),IRSPTR(INEXT))             ACVI1380
              KETX(K-1)=IEOR(KETX(K-1),KETUN(INEXT))                    ACVI1381
              KETX(K)=IEOR(KETX(K),KETUN(INEXT))                        ACVI1382
440         INEXT=IACTS(INEXT)                                          ACVI1383
C                                                                       ACVI1384
450       CONTINUE                                                      ACVI1385
460     CONTINUE                                                        ACVI1386
C       End application of Aij operator.                                ACVI1387
C                                                                       ACVI1388
        KETDIS=KETDIS+NWORDS                                            ACVI1389
470   CONTINUE                                                          ACVI1390
C                                                                       ACVI1391
C ****** LOOKING AT LINEAR HOMOGENEOUS SYSTEM OF EQUATIONS *************ACVI1392
C                                                                       ACVI1393
      WRITE(6,'(/A,I10)') ' NUMBER OF BASIS OBTAINED:',IQTMX            ACVI1394
      IF(TSMFI)                                                         ACVI1395
     1 WRITE(LFILE,'(/A,I10)') ' Number of basis obtained:',IQTMX       ACVI1396
      DO 490 I=1,IQTMX                                                  ACVI1397
        CALL TOUT(0,0,I,I,1,IQTREE)                                     ACVI1398
CB 08/16/91 11:40 ***>                                                  ACVI1399
        ILOC=IQTREE(-5)+IQTREE(-4)+1                                    ACVI1400
CB <*** 08/16/91 17:49                                                  ACVI1401
        IVECMX=IQTREE(ILOC)                                             ACVI1402
        INEXT=IQTREE(ILOC+1)                                            ACVI1403
C                                                                       ACVI1404
        DO 480 J=IVECMX,1,-1                                            ACVI1405
          IVEC(J)=IEQNS(1,INEXT)                                        ACVI1406
          DVEC(J)=DEQNS(INEXT)                                          ACVI1407
480     INEXT=IEQNS(2,INEXT)                                            ACVI1408
CW----------------------------------------------------------------------ACVI1409
CW    WRITING COEFFICIENT AND STATE OBTAINED.                           ACVI1410
      IF(TFILE) THEN                                                    ACVI1411
      WRITE(LFILE,'(/A,I10)') ' LINEAR EQUATION NUMBER:',I              ACVI1412
      WRITE(LFILE,1200)(IVEC(J),DVEC(J),J=1,IVECMX)                     ACVI1413
      ENDIF                                                             ACVI1414
CW----------------------------------------------------------------------ACVI1415
C                                                                       ACVI1416
        CALL ROWRED(IVEC, DVEC, JVEC, DXVEC,                            ACVI1417
     1              ICOEF(IHALF),DCOEF(IHALF),IDIS(KHALF),              ACVI1418
     2              ICOEF(JHALF),DCOEF(JHALF),IDIS(LHALF),IVECMX,*999)  ACVI1419
C       Changing sides in arrays for coefficients.                      ACVI1420
        IF(IVECMX.GT.0)THEN                                             ACVI1421
          KAUX=IHALF                                                    ACVI1422
          IHALF=JHALF                                                   ACVI1423
          JHALF=KAUX                                                    ACVI1424
          KAUX=KHALF                                                    ACVI1425
          KHALF=LHALF                                                   ACVI1426
          LHALF=KAUX                                                    ACVI1427
        ENDIF                                                           ACVI1428
490   CONTINUE                                                          ACVI1429
C                                                                       ACVI1430
C***********************************************************************ACVI1431
C                                                                       ACVI1432
C     SOLVE REDUCED MATRIX FOR COEFFICIENTS                             ACVI1433
C                                                                       ACVI1434
C***********************************************************************ACVI1435
      CALL GENSOL(ICOEF(IHALF),DCOEF(IHALF),IDIS(KHALF),DSOLMAT,NBKETS, ACVI1436
     1            LOCPTR,NSOL,MXNSOL,*999)                              ACVI1437
C                                                                       ACVI1438
C ****** WRITING SOLUTION **********************************************ACVI1439
      IF(NSOL.EQ.0)THEN                                                 ACVI1440
        WRITE(6,*) 'THERE IS ONLY TRIVIAL SOLUTION'                     ACVI1441
      ELSE                                                              ACVI1442
        WRITE(6,*) '     NUMBER OF SOLUTIONS:',NSOL                     ACVI1443
        WRITE(6,*) '         NUMBER OF BASIS:',NBKETS                   ACVI1444
        IF(TSMFI) THEN                                                  ACVI1445
          WRITE(LFILE,*) '         Number of basis:',NBKETS             ACVI1446
          WRITE(LFILE,*) '     Number of solutions:',NSOL               ACVI1447
          NSTART=1                                                      ACVI1448
          NSTOP=NBKETS                                                  ACVI1449
          DO 500 J=1,NSOL                                               ACVI1450
            WRITE(LFILE,'(/A,I5)') ' SOLUTION NUMBER',J                 ACVI1451
            WRITE(LFILE,1300) (DSOLMAT(I),I=NSTART,NSTOP)               ACVI1452
            NSTART=NSTOP+1                                              ACVI1453
            NSTOP=NSTOP+NBKETS                                          ACVI1454
500       CONTINUE                                                      ACVI1455
        ENDIF                                                           ACVI1456
C                                                                       ACVI1457
C       Compress for non-zero coefficients of the basis only.           ACVI1458
            IF(TCMPRS) THEN                                             ACVI1459
        CALL CMPRSS(NMAX,KETS,NBKETS,NSOL,DSOLMAT,*520)                 ACVI1460
        WRITE(6,*) '     NUMBER OF NEW BASIS:',NBKETS                   ACVI1461
        IF(TSMFI) THEN                                                  ACVI1462
          WRITE(LFILE,'(/A/)') ' COMPRESS !!!'                          ACVI1463
          WRITE(LFILE,*) '     Number of new basis:',NBKETS             ACVI1464
          NSTART=1                                                      ACVI1465
          NSTOP=NBKETS                                                  ACVI1466
          DO 510 J=1,NSOL                                               ACVI1467
            WRITE(LFILE,'(/A,I5)') ' SOLUTION NUMBER',J                 ACVI1468
            WRITE(LFILE,1300) (DSOLMAT(I),I=NSTART,NSTOP)               ACVI1469
            NSTART=NSTOP+1                                              ACVI1470
            NSTOP=NSTOP+NBKETS                                          ACVI1471
510       CONTINUE                                                      ACVI1472
        ENDIF                                                           ACVI1473
            ENDIF                                                       ACVI1474
      ENDIF                                                             ACVI1475
C                                                                       ACVI1476
520   RETURN                                                            ACVI1477
999   return 1                                                          ACVI1478
C980  WRITE(6,*) ' *** BOMB *** : TOO MANY NODES IN TREE'               ACVI1479
C     STOP                                                              ACVI1480
1000  FORMAT(/' CONSIDERING STATE ',I5/)                                ACVI1481
1100  FORMAT(F19.7)                                                     ACVI1482
1200  FORMAT(' ',T8,'STATE =',I7,T40,'COEFF =',F19.7)                   ACVI1483
1300  FORMAT(10(4(F13.7,5X)/))                                          ACVI1484
1400  FORMAT(10(4(1PE20.12)/))                                          ACVI1485
C ---*end of APPAIJ*----------------------------------------------------ACVI1486
      END                                                               ACVI1487
C ----------------------------------------------------------------------ACVI1488
C                                                                       ACVI1489
C                          **************                               ACVI1490
C                          *** SPARSE ***                               ACVI1491
C                          **************                               ACVI1492
C                                                                       ACVI1493
C The subroutine SPARSE row reduces a linear system of equations        ACVI1494
C produced by HWSGEN to obtain the coefficients for the                 ACVI1495
C expansion of the U(3)xU(4) highest weight state in terms              ACVI1496
C of the SETGMN basis vectors. The Gauss-Jordan method is used          ACVI1497
C for the row reduction, and is carried out as each equation is         ACVI1498
C passed in sequence to the subroutine. Whenever the number of          ACVI1499
C equations (NUMLINEQ) is less than the number of U(3)xU(4) basis       ACVI1500
C vectors (NBKETS), the solution to the system of equations is          ACVI1501
C obtained by setting the last NBKETS-NUMLINEQ coefficients each        ACVI1502
C equal to one in turn, with the remainder zero, and then               ACVI1503
C Gramm-Schmidt ortho-normalizing the resulting solution vectors.       ACVI1504
C                                                                       ACVI1505
C --------------------------------------------------------------------  ACVI1506
C                                                                       ACVI1507
C This subprogram has the following submodules.                         ACVI1508
C                                                                       ACVI1509
C   SPARSE : Initialize counters for later use.                         ACVI1510
C   ROWRED : Row reduce the system of equations as each row is          ACVI1511
C            passed in sequence to the subroutine.                      ACVI1512
C   GENSOL : Solve the row reduced system of equations and              ACVI1513
C            Gramm-Schmidt ortho-normalize the solution vectors.        ACVI1514
C                                                                       ACVI1515
C --------------------------------------------------------------------  ACVI1516
C                                                                       ACVI1517
C A call to SPARSE transfers the following variables.                   ACVI1518
C                                                                       ACVI1519
C   IVEC   : IKET values of non-zero elements in current row vector.    ACVI1520
C   DVEC   : Actual values of non-zero elements in current row vector.  ACVI1521
C   JVEC   : Temporary working array for storage of IVEC.               ACVI1522
C   DXVEC  : Temporary working array for row reduction of VEC.          ACVI1523
C   MXIVEC : Dim of IVEC, DVEC, JVEC, and DXVEC in main program.        ACVI1524
C   MXICOE : Dim of ICOEF, DCOEF, JCOEF, and DXCOEF in main program.    ACVI1525
C   MXIDIS : Dim of IDIS, JDIS, LOCPTR, and DSOLMAT in main progam.     ACVI1526
C   LOCPTR : Locate pointer for indicating independent solutions        ACVI1527
C            in DSOLMAT.                                                ACVI1528
C   ZEROMX : Variable which stores smallest difference between          ACVI1529
C            terms in row reduction process.                            ACVI1530
C   SAVEMN : Variable which stores largest difference between           ACVI1531
C            terms in row reduction process.                            ACVI1532
C                                                                       ACVI1533
C --------------------------------------------------------------------  ACVI1534
C                                                                       ACVI1535
C A call to ROWRED transfers the following variables.                   ACVI1536
C                                                                       ACVI1537
C   ICOEF  : IKET values of non-zero elements in previously row         ACVI1538
C            reduced system of equations.                               ACVI1539
C   DCOEF  : Actual values of non-zero elements in previously row       ACVI1540
C            reduced system of equations.                               ACVI1541
C   IDIS   : Index positions in ICOEF and DCOEF of last non-zero        ACVI1542
C            elements of row vectors.                                   ACVI1543
C   JCOEF  : Temporary working array for storage of ICOEF.              ACVI1544
C   DXCOEF : Temporary working array for current row reduction          ACVI1545
C            of DCOEF.                                                  ACVI1546
C   JDIS   : Temporary working array for storage of IDIS.               ACVI1547
C   IVECMX : Number of non-zero coefficients in currently               ACVI1548
C            passed vector DVEC.                                        ACVI1549
C                                                                       ACVI1550
C --------------------------------------------------------------------  ACVI1551
C                                                                       ACVI1552
C A call to GENSOL transfers the following variables.                   ACVI1553
C                                                                       ACVI1554
C   ICOEF   : Described above.                                          ACVI1555
C   DCOEF   : Described above.                                          ACVI1556
C   IDIS    : Described above.                                          ACVI1557
C   DSOLMAT : Solution matrix of HWS vectors on output.                 ACVI1558
C   IDIM    : Dimension of solution vectors.                            ACVI1559
C   NSOL    : Number of HWS solution vectors.                           ACVI1560
C   MXNSOL  : Max dimension of DSOLMAT in main program.                 ACVI1561
C                                                                       ACVI1562
C ----------------------------------------------------------------------ACVI1563
C                                                                       ACVI1564
      SUBROUTINE SPARSE( MXVEC, MXCOE, MXDIS )                          ACVI1565
      IMPLICIT REAL*8 (D)                                               ACVI1566
C                                                                       ACVI1567
      COMMON /SAVSMX/ MXIVEC, MXICOE, MXIDIS                            ACVI1568
      COMMON /SAVSPR/ DZEROMX, DSAVEMN, MXXCOE, ICOEMX, IDISMX          ACVI1569
      PARAMETER(DTOLER=1.0D-04)                                         ACVI1570
      DIMENSION ICOEF(*),IDIS(*),IVEC(*),JCOEF(*),JDIS(*),JVEC(*)       ACVI1571
      DIMENSION DCOEF(*),DVEC(*),DXCOEF(*),DXVEC(*),DSOLMAT(IDIM,*)     ACVI1572
C     Auxilliary storage.                                               ACVI1573
      DIMENSION LOCPTR(*)                                               ACVI1574
C     Initial values for internal counters.                             ACVI1575
      MXIVEC=MXVEC                                                      ACVI1576
      MXICOE=MXCOE                                                      ACVI1577
      MXIDIS=MXDIS                                                      ACVI1578
      MXXCOE=MXICOE/2                                                   ACVI1579
      ICOEMX=0                                                          ACVI1580
      IDISMX=0                                                          ACVI1581
      DZEROMX=0.0D0                                                     ACVI1582
      DSAVEMN=1.0D+30                                                   ACVI1583
      RETURN                                                            ACVI1584
C                                                                       ACVI1585
C ----------------------------------------------------------------------ACVI1586
C     **************                                                    ACVI1587
C     *** ROWRED ***                                                    ACVI1588
C     **************                                                    ACVI1589
C ----------------------------------------------------------------------ACVI1590
C                                                                       ACVI1591
      ENTRY ROWRED(IVEC, DVEC, JVEC, DXVEC,                             ACVI1592
     1             ICOEF,DCOEF,IDIS,JCOEF,DXCOEF,JDIS,IVECMX,*)         ACVI1593
C     Reduce new vector and insert into solution matrix.                ACVI1594
      IF(IDISMX.GT.0)THEN                                               ACVI1595
C     (1) Reduce new vector by running through previous rows.           ACVI1596
      IROW=0                                                            ACVI1597
110     IROW=IROW+1                                                     ACVI1598
C       Value of IDIS(IROW+1) gives the position of the last            ACVI1599
C       element in the array ICOEF(*) that belongs to a                 ACVI1600
C       given row (IROW).                                               ACVI1601
        IROWMN=IDIS(IROW)+1                                             ACVI1602
        ICURR=IROWMN                                                    ACVI1603
        IFIRST=ICOEF(ICURR)                                             ACVI1604
C                                                                       ACVI1605
C       Run through elements of input vector and compare with           ACVI1606
C       current vector in reduced system.                               ACVI1607
        JINPUT=1                                                        ACVI1608
120     IF(IVEC(JINPUT).GT.IFIRST)THEN                                  ACVI1609
C         (a)If non-zero component considered in input vector           ACVI1610
C            appears after first non-zero element in current            ACVI1611
C            previous row, then go to the next row.                     ACVI1612
          JINPUT=IVECMX+1                                               ACVI1613
          JVECMX=IVECMX                                                 ACVI1614
        ELSE IF(IVEC(JINPUT).EQ.IFIRST)THEN                             ACVI1615
C         (b)If non_zero component considered in vector                 ACVI1616
C            appears in the same place as first non_zero element        ACVI1617
C            in actual row, then perform reduction process.             ACVI1618
C                                                                       ACVI1619
C         Do not change previous elements in vector.                    ACVI1620
          DCONST=DVEC(JINPUT)                                           ACVI1621
          MARK=JINPUT                                                   ACVI1622
          JVECMX=JINPUT-1                                               ACVI1623
          ICURR=ICURR+1                                                 ACVI1624
          JINPUT=JINPUT+1                                               ACVI1625
          IROWMX=IDIS(IROW+1)                                           ACVI1626
C         Actual reduction process.                                     ACVI1627
C         IFLAG = 1 stops reduction process.                            ACVI1628
          IFLAG=0                                                       ACVI1629
140       IF(ICURR.LE.IROWMX)THEN                                       ACVI1630
            IF(JINPUT.LE.IVECMX)THEN                                    ACVI1631
              IF(ICOEF(ICURR).LT.IVEC(JINPUT))THEN                      ACVI1632
C             Assign a new element to current vector                    ACVI1633
C             and shift to next element in current row.                 ACVI1634
                JVECMX=JVECMX+1                                         ACVI1635
                JVEC(JVECMX)=ICOEF(ICURR)                               ACVI1636
                DXVEC(JVECMX)=-DCONST*DCOEF(ICURR)                      ACVI1637
                ICURR=ICURR+1                                           ACVI1638
              ELSE IF(ICOEF(ICURR).EQ.IVEC(JINPUT))THEN                 ACVI1639
C             Alter an element in current vector and                    ACVI1640
C             shift to next element in both rows.                       ACVI1641
                DEX=DVEC(JINPUT)-DCONST*DCOEF(ICURR)                    ACVI1642
                DABSEX=DABS(DEX)                                        ACVI1643
                IF(DABSEX.LT.DTOLER)THEN                                ACVI1644
                  IF(DABSEX.GT.DZEROMX) DZEROMX=DABSEX                  ACVI1645
                ELSE                                                    ACVI1646
                  IF(DABSEX.LT.DSAVEMN) DSAVEMN=DABSEX                  ACVI1647
                  JVECMX=JVECMX+1                                       ACVI1648
                  JVEC(JVECMX)=IVEC(JINPUT)                             ACVI1649
                  DXVEC(JVECMX)=DEX                                     ACVI1650
                END IF                                                  ACVI1651
                ICURR=ICURR+1                                           ACVI1652
                JINPUT=JINPUT+1                                         ACVI1653
              ELSE IF(ICOEF(ICURR).GT.IVEC(JINPUT))THEN                 ACVI1654
C             Retain element which already exists in                    ACVI1655
C             current vector.                                           ACVI1656
                JVECMX=JVECMX+1                                         ACVI1657
                JVEC(JVECMX)=IVEC(JINPUT)                               ACVI1658
                DXVEC(JVECMX)=DVEC(JINPUT)                              ACVI1659
                JINPUT=JINPUT+1                                         ACVI1660
              END IF                                                    ACVI1661
            ELSE                                                        ACVI1662
              IFLAG=1                                                   ACVI1663
C             If no more non-zero elements found in current             ACVI1664
C             vector, then assign as many new values as                 ACVI1665
C             remaining elements in row of input vector.                ACVI1666
              DO 150 KI=ICURR,IROWMX                                    ACVI1667
                JVECMX=JVECMX+1                                         ACVI1668
                JVEC(JVECMX)=ICOEF(KI)                                  ACVI1669
                DXVEC(JVECMX)=-DCONST*DCOEF(KI)                         ACVI1670
150           CONTINUE                                                  ACVI1671
              JINPUT=IVECMX                                             ACVI1672
            END IF                                                      ACVI1673
          ELSE                                                          ACVI1674
            IFLAG=1                                                     ACVI1675
            IF(JINPUT.LE.IVECMX)THEN                                    ACVI1676
C           If no more non-zero elements found in current               ACVI1677
C           row, then leave remaining elements in input                 ACVI1678
C           vector as is.                                               ACVI1679
              DO 160 KJ=JINPUT,IVECMX                                   ACVI1680
                JVECMX=JVECMX+1                                         ACVI1681
                JVEC(JVECMX)=IVEC(KJ)                                   ACVI1682
                DXVEC(JVECMX)=DVEC(KJ)                                  ACVI1683
160           CONTINUE                                                  ACVI1684
            ELSE                                                        ACVI1685
              IF(JVECMX.EQ.0)THEN                                       ACVI1686
C             Vector was a linear combination of rows.                  ACVI1687
                IROW=IDISMX                                             ACVI1688
                JINPUT=IVECMX                                           ACVI1689
              END IF                                                    ACVI1690
            END IF                                                      ACVI1691
          END IF                                                        ACVI1692
          IF(IFLAG.EQ.0) GO TO 140                                      ACVI1693
          IF(JVECMX.GT.MXIVEC)                                          ACVI1694
     1      CALL ATTN(' ROWRED: Overflow in number of rows.',*999)      ACVI1695
C    1      CALL ERROR(' ROWRED: ALLOWED NUMBER OF ROWS EXCEEDED')      ACVI1696
C         Vector reduction completed                                    ACVI1697
          IVECMX=JVECMX                                                 ACVI1698
          DO 170 KM=MARK,IVECMX                                         ACVI1699
            IVEC(KM)=JVEC(KM)                                           ACVI1700
            DVEC(KM)=DXVEC(KM)                                          ACVI1701
170       CONTINUE                                                      ACVI1702
        ELSE IF(IVEC(JINPUT).LT.IFIRST)THEN                             ACVI1703
C         (c)If non-zero component of input vector appears before       ACVI1704
C            first non-zero compenent of current previous vector,       ACVI1705
C            then skip to next component.                               ACVI1706
          JINPUT=JINPUT+1                                               ACVI1707
        END IF                                                          ACVI1708
        IF(JINPUT.LE.IVECMX)GO TO 120                                   ACVI1709
        IF(IROW.LT.IDISMX) GO TO 110                                    ACVI1710
      END IF                                                            ACVI1711
C                                                                       ACVI1712
C     Reduce previous rows in reduced system of equations and           ACVI1713
C     append reduced vector as new row.                                 ACVI1714
      IF(IVECMX.GT.0)THEN                                               ACVI1715
        IF(JVECMX.GT.MXIVEC)                                            ACVI1716
     1    CALL ATTN(' ROWRED: Overflow in number of rows.',*999)        ACVI1717
C    1    CALL ERROR(' ROWRED: ALLOWED NUMBER OF ROWS EXCEEDED')        ACVI1718
C                                                                       ACVI1719
C       Set first non-zero element in vector equal to 1.                ACVI1720
        DCONST=DVEC(1)                                                  ACVI1721
        DO 180 I=1,IVECMX                                               ACVI1722
180       DVEC(I)=DVEC(I)/DCONST                                        ACVI1723
C                                                                       ACVI1724
C       (2) Now reduce all rows existing previously to input vector.    ACVI1725
        JDIS(1)=0                                                       ACVI1726
        JCOEMX=0                                                        ACVI1727
        IF(IDISMX.GT.0)THEN                                             ACVI1728
          DO 210 IROW=1,IDISMX                                          ACVI1729
            IROWMN=IDIS(IROW)+1                                         ACVI1730
            IROWMX=IDIS(IROW+1)                                         ACVI1731
            JINPUT=1                                                    ACVI1732
            IFIRST=IVEC(JINPUT)                                         ACVI1733
            ICURR=IROWMN                                                ACVI1734
            IFLAG=0                                                     ACVI1735
220         IF(ICOEF(ICURR).LT.IFIRST)THEN                              ACVI1736
C             Cannot reduce current row with given input vector.        ACVI1737
C             Save values of current row because JCOEF(*) will          ACVI1738
C             become ICOEF(*) in the next call to ROWRED.               ACVI1739
230           JCOEMX=JCOEMX+1                                           ACVI1740
              JCOEF(JCOEMX)=ICOEF(ICURR)                                ACVI1741
              DXCOEF(JCOEMX)=DCOEF(ICURR)                               ACVI1742
              ICURR=ICURR+1                                             ACVI1743
              IF(ICOEF(ICURR).LT.IFIRST.AND.ICURR.LE.IROWMX)            ACVI1744
     1          GO TO 230                                               ACVI1745
              IF(ICURR.GT.IROWMX) IFLAG=1                               ACVI1746
            ELSE IF(ICOEF(ICURR).EQ.IFIRST)THEN                         ACVI1747
C             Make reduction in row considered.                         ACVI1748
              DCONST=DCOEF(ICURR)                                       ACVI1749
              ICURR=ICURR+1                                             ACVI1750
              JINPUT=JINPUT+1                                           ACVI1751
C             Run through elements of current previous vector           ACVI1752
C             and compare with elements of input vector.                ACVI1753
240           IF(ICURR.LE.IROWMX)THEN                                   ACVI1754
                IF(JINPUT.LE.IVECMX)THEN                                ACVI1755
                  IF(ICOEF(ICURR).LT.IVEC(JINPUT))THEN                  ACVI1756
C                   Do not change element in row and save.              ACVI1757
250                 JCOEMX=JCOEMX+1                                     ACVI1758
                    JCOEF(JCOEMX)=ICOEF(ICURR)                          ACVI1759
                    DXCOEF(JCOEMX)=DCOEF(ICURR)                         ACVI1760
                    ICURR=ICURR+1                                       ACVI1761
                    IF(ICOEF(ICURR).LT.IVEC(JINPUT)                     ACVI1762
     1                 .AND.ICURR.LE.IROWMX) GO TO 250                  ACVI1763
                  ELSE IF(ICOEF(ICURR).EQ.IVEC(JINPUT))THEN             ACVI1764
C                   Change value of element in row and shift            ACVI1765
C                   to next element in both rows.                       ACVI1766
                    DEX=DCOEF(ICURR)-DCONST*DVEC(JINPUT)                ACVI1767
                    DABSEX=DABS(DEX)                                    ACVI1768
                    IF(DABSEX.LT.DTOLER)THEN                            ACVI1769
                      IF(DABSEX.GT.DZEROMX) DZEROMX=DABSEX              ACVI1770
                    ELSE                                                ACVI1771
                      IF(DABSEX.LT.DSAVEMN) DSAVEMN=DABSEX              ACVI1772
                      JCOEMX=JCOEMX+1                                   ACVI1773
                      JCOEF(JCOEMX)=ICOEF(ICURR)                        ACVI1774
                      DXCOEF(JCOEMX)=DEX                                ACVI1775
                    END IF                                              ACVI1776
                    ICURR=ICURR+1                                       ACVI1777
                    JINPUT=JINPUT+1                                     ACVI1778
                  ELSE IF(ICOEF(ICURR).GT.IVEC(JINPUT))THEN             ACVI1779
C                   Assign a new element into row.                      ACVI1780
260                 JCOEMX=JCOEMX+1                                     ACVI1781
                    JCOEF(JCOEMX)=IVEC(JINPUT)                          ACVI1782
                    DXCOEF(JCOEMX)=-DCONST*DVEC(JINPUT)                 ACVI1783
                    JINPUT=JINPUT+1                                     ACVI1784
                    IF(ICOEF(ICURR).GT.IVEC(JINPUT)                     ACVI1785
     1                 .AND.JINPUT.LE.IVECMX) GO TO 260                 ACVI1786
                  END IF                                                ACVI1787
                ELSE                                                    ACVI1788
                  IFLAG=1                                               ACVI1789
C                 If all elements in vector exhausted, then do          ACVI1790
C                 not change remaining elements in actual row.          ACVI1791
                  DO 270 KI=ICURR,IROWMX                                ACVI1792
                    JCOEMX=JCOEMX+1                                     ACVI1793
                    JCOEF(JCOEMX)=ICOEF(KI)                             ACVI1794
270                 DXCOEF(JCOEMX)=DCOEF(KI)                            ACVI1795
                END IF                                                  ACVI1796
              ELSE                                                      ACVI1797
                IFLAG=1                                                 ACVI1798
                IF(JINPUT.LE.IVECMX)THEN                                ACVI1799
C                 If all elements in current row exhausted,             ACVI1800
C                 then assign as many new elements in this row          ACVI1801
C                 as elements not yet considered in input               ACVI1802
C                 vector.                                               ACVI1803
                  DO 280 KJ=JINPUT,IVECMX                               ACVI1804
                    JCOEMX=JCOEMX+1                                     ACVI1805
                    JCOEF(JCOEMX)=IVEC(KJ)                              ACVI1806
280                 DXCOEF(JCOEMX)=-DCONST*DVEC(KJ)                     ACVI1807
                ELSE                                                    ACVI1808
                  IF(JCOEMX.LE.JDIS(IROW))                              ACVI1809
     1              CALL ATTN(' ROWRED: Internal check fails.',*999)    ACVI1810
C    1              CALL ERROR(' ROWRED: INTERNAL CHECK FAILS')         ACVI1811
                ENDIF                                                   ACVI1812
              ENDIF                                                     ACVI1813
              IF(IFLAG.EQ.0)GO TO 240                                   ACVI1814
            ELSE IF(ICOEF(ICURR).GT.IFIRST)THEN                         ACVI1815
C             Leave elements in current previous row unchanged.         ACVI1816
              IFLAG=1                                                   ACVI1817
              DO 290 KI=ICURR,IROWMX                                    ACVI1818
                JCOEMX=JCOEMX+1                                         ACVI1819
                JCOEF(JCOEMX)=ICOEF(KI)                                 ACVI1820
290             DXCOEF(JCOEMX)=DCOEF(KI)                                ACVI1821
            END IF                                                      ACVI1822
C           Repeat until vector or current row is exhausted.            ACVI1823
            IF(IFLAG.EQ.0) GO TO 220                                    ACVI1824
            JDIS(IROW+1)=JCOEMX                                         ACVI1825
210       CONTINUE                                                      ACVI1826
        END IF                                                          ACVI1827
C     (3) Append processed input vector to end of matrix.               ACVI1828
        IF(JCOEMX+IVECMX.GT.MXXCOE) CALL ATTN(' ROWRED: '//             ACVI1829
     1    'Overflow in arrays JCOEF and DXCOEF.',*999)                  ACVI1830
C    1    CALL ERROR(' ROWRED: OVERFLOW IN ARRAYS JCOEF AND DXCOEF')    ACVI1831
        DO 310 INEW=1,IVECMX                                            ACVI1832
          JCOEMX=JCOEMX+1                                               ACVI1833
          JCOEF(JCOEMX)=IVEC(INEW)                                      ACVI1834
  310     DXCOEF(JCOEMX)=DVEC(INEW)                                     ACVI1835
        JDISMX=IDISMX+1                                                 ACVI1836
CW----------------------------------------------------------------------ACVI1837
CW    Show me all the rows after complete reduction.                    ACVI1838
C     WRITE(6,*)                                                        ACVI1839
C     WRITE(6,*) 'SHOW ALL THE ROWS AFTER REDUCTION'                    ACVI1840
C     WRITE(6,*) 'JCOEMX = ',JCOEMX                                     ACVI1841
C     DO 325,I=1,JCOEMX                                                 ACVI1842
C325    WRITE(6,*) 'I = ',I,' JCOEF = ',JCOEF(I),' DXCOEF = ',DXCOEF(I) ACVI1843
C     WRITE(6,*)                                                        ACVI1844
CW----------------------------------------------------------------------ACVI1845
        IF(JDISMX.GT.MXIDIS)                                            ACVI1846
     1    CALL ATTN(' ROWRED: Overflow in array JDIS.',*999)            ACVI1847
C    1    CALL ERROR(' ROWRED: OVERFLOW IN ARRAY JDIS')                 ACVI1848
        JDIS(JDISMX+1)=JCOEMX                                           ACVI1849
        ICOEMX=JCOEMX                                                   ACVI1850
        IDISMX=JDISMX                                                   ACVI1851
      END IF                                                            ACVI1852
      RETURN                                                            ACVI1853
C                                                                       ACVI1854
C ----------------------------------------------------------------------ACVI1855
C     **************                                                    ACVI1856
C     *** GENSOL ***                                                    ACVI1857
C     **************                                                    ACVI1858
C ----------------------------------------------------------------------ACVI1859
C                                                                       ACVI1860
      ENTRY GENSOL(ICOEF,DCOEF,IDIS,DSOLMAT,IDIM,LOCPTR,NSOL,MXNSOL,*)  ACVI1861
C                                                                       ACVI1862
C     Generate orthonormal solution vectors and store them in           ACVI1863
C     DSOLMAT.                                                          ACVI1864
CW----------------------------------------------------------------------ACVI1865
      WRITE(6,*)                                                        ACVI1866
      WRITE(6,*) 'NUMBER OF ROWS AND COLUMNS IN SOLUTION MATRIX'        ACVI1867
      WRITE(6,*)   ' IDISMX = ',IDISMX,' IDIM = ',IDIM                  ACVI1868
CW----------------------------------------------------------------------ACVI1869
      DO 410 I=1,IDIM                                                   ACVI1870
410     LOCPTR(I)=0                                                     ACVI1871
      IF(IDISMX.GT.0)THEN                                               ACVI1872
C       Using LOCPTR(ICOL)=IROW to identify the row (IROW) in which     ACVI1873
C       column ICOL is the first non-zero entry. Values of ICOL         ACVI1874
C       which represent independent variables have LOCPTR(ICOL)=0.      ACVI1875
        DO 420 IROW=1,IDISMX                                            ACVI1876
          IROWMN=IDIS(IROW)+1                                           ACVI1877
          ICOL=ICOEF(IROWMN)                                            ACVI1878
420       LOCPTR(ICOL)=IROW                                             ACVI1879
      END IF                                                            ACVI1880
C     Determining how many independent parameters there are and         ACVI1881
C     identifying them with negative numbers.                           ACVI1882
      NSOL=0                                                            ACVI1883
      DO 430 ICOL=1,IDIM                                                ACVI1884
        IF(LOCPTR(ICOL).EQ.0)THEN                                       ACVI1885
          NSOL=NSOL+1                                                   ACVI1886
          LOCPTR(ICOL)=-NSOL                                            ACVI1887
        END IF                                                          ACVI1888
430   CONTINUE                                                          ACVI1889
      IF(NSOL.NE.IDIM-IDISMX) THEN                                      ACVI1890
        WRITE(6,'(//A,2I11)') ' NSOL(??)',NSOL,IDIM-IDISMX              ACVI1891
        CALL ERROR(' GENSOL: ERROR IN PROCESSING INDEPENDENT VARIABLES')ACVI1892
      ENDIF                                                             ACVI1893
      IF(IDIM*NSOL.GT.MXNSOL)                                           ACVI1894
     1  CALL ATTN(' GENSOL: Overflow in array DSOLMAT.',*999)           ACVI1895
C    1  CALL ERROR(' GENSOL: OVERFLOW IN ARRAY DSOLMAT')                ACVI1896
      IF(NSOL.NE.0)THEN                                                 ACVI1897
        DO 440 J=1,NSOL                                                 ACVI1898
          DO 440 I=1,IDIM                                               ACVI1899
            DSOLMAT(I,J)=0.0D0                                          ACVI1900
440     CONTINUE                                                        ACVI1901
C       Each column in DSOLMAT(IVECTOR,ISOL) represents an              ACVI1902
C       independent solution vector. Each ISOL column of DSOLMAT        ACVI1903
C       is determined by assigning one of the IDIM-IDISMX independent   ACVI1904
C       variables a value of 1 and the remainder a value of 0.          ACVI1905
        DO 460 ICOL=1,IDIM                                              ACVI1906
          JROW=LOCPTR(ICOL)                                             ACVI1907
          IF(JROW.LT.0)THEN                                             ACVI1908
            DSOLMAT(ICOL,-JROW)=1.0D0                                   ACVI1909
          ELSE IF(JROW.GT.0)THEN                                        ACVI1910
C           Insert non-zero elements of reduced row into solution       ACVI1911
C           matrix. Observe that first element in row was already       ACVI1912
C           considered.                                                 ACVI1913
            IROWMN=IDIS(JROW)+2                                         ACVI1914
            IROWMX=IDIS(JROW+1)                                         ACVI1915
            IF(IROWMN.LE.IROWMX)THEN                                    ACVI1916
              DO 470 K=IROWMN,IROWMX                                    ACVI1917
                KCOL=ICOEF(K)                                           ACVI1918
                JROW=LOCPTR(KCOL)                                       ACVI1919
C               ICOEF(KCOL) must correspond to one independent variable.ACVI1920
C               IF(JROW.GE.0)                                           ACVI1921
C    1            CALL ERROR(' GENSOL: '//                              ACVI1922
C    2              'ERROR IN CONSTRUCTING MATRIX SOLUTION')            ACVI1923
                IF(JROW.GE.0) CALL ATTN(' GENSOL: '//                   ACVI1924
     1            'Error in constructing matrix solution.',*999)        ACVI1925
                DSOLMAT(ICOL,-JROW)=-DCOEF(K)                           ACVI1926
470           CONTINUE                                                  ACVI1927
            END IF                                                      ACVI1928
          END IF                                                        ACVI1929
460     CONTINUE                                                        ACVI1930
C                                                                       ACVI1931
C       (2) Gramm-Schmidt ortho-normalize the solution basis vectors.   ACVI1932
        DO 510,J=1,NSOL                                                 ACVI1933
C        Extracting from the next column an independent vector          ACVI1934
C        orthogonal to all the previous ones.                           ACVI1935
         IF(J.GT.1)THEN                                                 ACVI1936
           DO 520 K=1,J-1                                               ACVI1937
             DOT=0.D0                                                   ACVI1938
             DO 530 I=1,IDIM                                            ACVI1939
530            DOT=DOT+DSOLMAT(I,K)*DSOLMAT(I,J)                        ACVI1940
             DO 540 I=1,IDIM                                            ACVI1941
540            DSOLMAT(I,J)=DSOLMAT(I,J)-DOT*DSOLMAT(I,K)               ACVI1942
520        CONTINUE                                                     ACVI1943
         END IF                                                         ACVI1944
C        Normalizing vector.                                            ACVI1945
         DRENORM=0.D0                                                   ACVI1946
         DO 550 I=1,IDIM                                                ACVI1947
550        DRENORM=DRENORM+DSOLMAT(I,J)**2                              ACVI1948
         DRENORM=DSQRT(DRENORM)                                         ACVI1949
         DO 560 I=IDIM,1,-1                                             ACVI1950
560        DSOLMAT(I,J)=DSOLMAT(I,J)/DRENORM                            ACVI1951
510     CONTINUE                                                        ACVI1952
      END IF                                                            ACVI1953
      RETURN                                                            ACVI1954
999   RETURN 1                                                          ACVI1955
C990  STOP                                                              ACVI1956
      END                                                               ACVI1957
CB----------------------------------------------------------------------ACVI1958
C                          ***************                              ACVI1959
C                          **  UTILITY ***                              ACVI1960
C                          ***************                              ACVI1961
C ----------------------------------------------------------------------ACVI1962
      SUBROUTINE SETBIN(LFILX,TWRITX,BINX,NEMPTX,NFIRSTX,NMAXX,NTOTALX) ACVI1963
C     Generating the binary I/O.                                        ACVI1964
C                                                                       ACVI1965
      COMMON /SAVARG/ LFILE, TWRITE, NEMPTY, NFIRST, NMAX, NTOTAL       ACVI1966
      COMMON /SAVCHR/ BIN                                               ACVI1967
      COMMON /SAVINT/ IPM, NMAXP1, NMAX2, NMAXP2, NMAX3                 ACVI1968
      LOGICAL TWRITE, TWRITX                                            ACVI1969
      CHARACTER*(*) BINX                                                ACVI1970
      CHARACTER BIN*115                                                 ACVI1971
      LFILE = LFILX                                                     ACVI1972
      TWRITE = TWRITX                                                   ACVI1973
      NEMPTY = NEMPTX                                                   ACVI1974
      NFIRST = NFIRSTX                                                  ACVI1975
      NMAX = NMAXX                                                      ACVI1976
      NTOTAL = NTOTALX                                                  ACVI1977
      NMAXP1=NMAX+6                                                     ACVI1978
      NMAX2=NMAX+NMAX+5                                                 ACVI1979
      NMAXP2=NMAX2+6                                                    ACVI1980
      NMAX3=3*NMAX+10                                                   ACVI1981
      IF(NFIRST.LE.0)THEN                                               ACVI1982
        IPM=1                                                           ACVI1983
      ELSE                                                              ACVI1984
        IPM=-1                                                          ACVI1985
      ENDIF                                                             ACVI1986
      DO J=1,115                                                        ACVI1987
        BIN(J:J)=' '                                                    ACVI1988
      ENDDO                                                             ACVI1989
      RETURN                                                            ACVI1990
C                                                                       ACVI1991
      ENTRY DBOVRL(ILEFT,IL,IR,IRGHT)                                   ACVI1992
C     Converting a number from decimal to binary. (1a)                  ACVI1993
C-----------------------------------------------------------------------ACVI1994
C     Reading bits in number.                                           ACVI1995
      I=IPM                                                             ACVI1996
      DO J=1,NMAX                                ! LEFT                 ACVI1997
        IF(MOD(J,NEMPTY).NE.0)THEN                                      ACVI1998
           IF(BTEST(ILEFT,NFIRST+I))THEN                                ACVI1999
              BIN(J:J)='1'                                              ACVI2000
           ELSE                                                         ACVI2001
              BIN(J:J)='0'                                              ACVI2002
           ENDIF                                                        ACVI2003
           I=I+IPM                                                      ACVI2004
        ENDIF                                                           ACVI2005
      ENDDO                                                             ACVI2006
                  IF(NMAX.EQ.NTOTAL) THEN                               ACVI2007
C     Writing number in binary form.                                    ACVI2008
      IF(TWRITE) THEN                                                   ACVI2009
C       WRITE(LFILE,100) ILEFT,BIN                                      ACVI2010
        write(lfile,100) ileft, (bin(i:i), i=1, j)                      ACVI2011
      ELSE                                                              ACVI2012
        DO J = 1, NTOTAL                                                ACVI2013
          BINX(J:J) = BIN(J:J)                                          ACVI2014
        END DO                                                          ACVI2015
      END IF                                                            ACVI2016
                  ELSE                                                  ACVI2017
      I=IPM                                                             ACVI2018
      NMIN=MOD(NMAXP1-1,NEMPTY)                                         ACVI2019
      DO J=NMAXP1,NMAX2                          ! MIDDLE               ACVI2020
        IF(MOD(J,NEMPTY).NE.NMIN)THEN                                   ACVI2021
           NF=NFIRST+I                                                  ACVI2022
           IF(BTEST(IL,NF))THEN                                         ACVI2023
              IF(BTEST(IR,NF))THEN                                      ACVI2024
                 BIN(J:J)='I'                                           ACVI2025
              ELSE                                                      ACVI2026
                 BIN(J:J)='V'                                           ACVI2027
              ENDIF                                                     ACVI2028
           ELSE                                                         ACVI2029
              IF(BTEST(IR,NF))THEN                                      ACVI2030
                 BIN(J:J)='A'                                           ACVI2031
              ELSE                                                      ACVI2032
                 BIN(J:J)='-'                                           ACVI2033
              ENDIF                                                     ACVI2034
           ENDIF                                                        ACVI2035
           I=I+IPM                                                      ACVI2036
        ENDIF                                                           ACVI2037
      ENDDO                                                             ACVI2038
      I=IPM                                                             ACVI2039
      NMIN=MOD(NMAXP2-1,NEMPTY)                                         ACVI2040
      DO J=NMAXP2,NMAX3                          ! RIGHT                ACVI2041
        IF(MOD(J,NEMPTY).NE.NMIN)THEN                                   ACVI2042
           IF(BTEST(IRGHT,NFIRST+I))THEN                                ACVI2043
              BIN(J:J)='1'                                              ACVI2044
           ELSE                                                         ACVI2045
              BIN(J:J)='0'                                              ACVI2046
           ENDIF                                                        ACVI2047
           I=I+IPM                                                      ACVI2048
        ENDIF                                                           ACVI2049
      ENDDO                                                             ACVI2050
C     Writing number in binary form.                                    ACVI2051
C     IF(TWRITE) WRITE(LFILE,'(T2,A)') BIN                              ACVI2052
      if( twrite ) write(lfile,'(t2,115a)') (bin(i:i), i=1, j)          ACVI2053
                  ENDIF                                                 ACVI2054
      RETURN                                                            ACVI2055
C                                                                       ACVI2056
      ENTRY DBOVRP(NBITS,ILEFT,IL,IR,IRGHT)                             ACVI2057
C     Converting a number from decimal to binary. (1b)                  ACVI2058
C-----------------------------------------------------------------------ACVI2059
C     Reading bits in number.                                           ACVI2060
      I=IPM                                                             ACVI2061
      J=0                                                               ACVI2062
      DO WHILE(J.LE.NMAX .AND. IABS(I).LE.NBITS)   ! LEFT               ACVI2063
        J=J+1                                                           ACVI2064
        IF(MOD(J,NEMPTY).NE.0) THEN                                     ACVI2065
           IF(BTEST(ILEFT,NFIRST+I))THEN                                ACVI2066
              BIN(J:J)='1'                                              ACVI2067
           ELSE                                                         ACVI2068
              BIN(J:J)='0'                                              ACVI2069
           ENDIF                                                        ACVI2070
           I=I+IPM                                                      ACVI2071
        ENDIF                                                           ACVI2072
      ENDDO                                                             ACVI2073
      I=IPM                                                             ACVI2074
      J=NMAXP1-1                                                        ACVI2075
      NMIN=MOD(J,NEMPTY)                                                ACVI2076
      DO WHILE(J.LE.NMAX2 .AND. IABS(I).LE.NBITS)   ! MIDDLE            ACVI2077
        J=J+1                                                           ACVI2078
        IF(MOD(J,NEMPTY).NE.NMIN) THEN                                  ACVI2079
           NF=NFIRST+I                                                  ACVI2080
           IF(BTEST(IL,NF))THEN                                         ACVI2081
              IF(BTEST(IR,NF))THEN                                      ACVI2082
                 BIN(J:J)='I'                                           ACVI2083
              ELSE                                                      ACVI2084
                 BIN(J:J)='A'                                           ACVI2085
              ENDIF                                                     ACVI2086
           ELSE                                                         ACVI2087
              IF(BTEST(IR,NF))THEN                                      ACVI2088
                 BIN(J:J)='V'                                           ACVI2089
              ELSE                                                      ACVI2090
                 BIN(J:J)='-'                                           ACVI2091
              ENDIF                                                     ACVI2092
           ENDIF                                                        ACVI2093
           I=I+IPM                                                      ACVI2094
        ENDIF                                                           ACVI2095
      ENDDO                                                             ACVI2096
      I=IPM                                                             ACVI2097
      J=NMAXP2-1                                                        ACVI2098
      NMIN=MOD(J,NEMPTY)                                                ACVI2099
      DO WHILE(J.LE.NMAX3 .AND. IABS(I).LE.NBITS)   ! RIGHT             ACVI2100
        J=J+1                                                           ACVI2101
        IF(MOD(J,NEMPTY).NE.NMIN) THEN                                  ACVI2102
           IF(BTEST(IRGHT,NFIRST+I))THEN                                ACVI2103
              BIN(J:J)='1'                                              ACVI2104
           ELSE                                                         ACVI2105
              BIN(J:J)='0'                                              ACVI2106
           ENDIF                                                        ACVI2107
           I=I+IPM                                                      ACVI2108
        ENDIF                                                           ACVI2109
      ENDDO                                                             ACVI2110
C     Writing number in binary form.                                    ACVI2111
      IF(TWRITE) THEN                                                   ACVI2112
C       WRITE(LFILE,'(T2,A)') BIN                                       ACVI2113
        write(lfile,'(t2,115a)') (bin(i:i), i=1, j)                     ACVI2114
      ELSE                                                              ACVI2115
        DO J = 1, NTOTAL                                                ACVI2116
          BINX(J:J) = BIN(J:J)                                          ACVI2117
        END DO                                                          ACVI2118
      END IF                                                            ACVI2119
      RETURN                                                            ACVI2120
C                                                                       ACVI2121
      ENTRY DBOVRC(NBITS,ILEFT,IL,IR,IRGHT)                             ACVI2122
C     Converting a number from decimal to binary. (2)                   ACVI2123
C-----------------------------------------------------------------------ACVI2124
C     Reading bits in number.                                           ACVI2125
      I=IPM                                                             ACVI2126
      J=0                                                               ACVI2127
      DO WHILE(J.LE.NMAX .AND. IABS(I).LE.NBITS)   ! LEFT               ACVI2128
        J=J+1                                                           ACVI2129
        IF(MOD(J,NEMPTY).NE.0) THEN                                     ACVI2130
           IF(BTEST(ILEFT,NFIRST+I))THEN                                ACVI2131
              BIN(J:J)='1'                                              ACVI2132
           ELSE                                                         ACVI2133
              BIN(J:J)='0'                                              ACVI2134
           ENDIF                                                        ACVI2135
           I=I+IPM                                                      ACVI2136
        ENDIF                                                           ACVI2137
      ENDDO                                                             ACVI2138
                  IF(NMAX.EQ.NTOTAL) THEN                               ACVI2139
C     Writing number in binary form.                                    ACVI2140
      IF(TWRITE) THEN                                                   ACVI2141
C       WRITE(LFILE,100) ILEFT,BIN                                      ACVI2142
        write(lfile,100) ileft, (bin(i:i), i=1, j)                      ACVI2143
      ELSE                                                              ACVI2144
        DO J = 1, NTOTAL                                                ACVI2145
          BINX(J:J) = BIN(J:J)                                          ACVI2146
        END DO                                                          ACVI2147
      END IF                                                            ACVI2148
                  ELSE                                                  ACVI2149
      NMAXP1=NBITS+6                                                    ACVI2150
      NMAX2=NBITS+NBITS+5                                               ACVI2151
      NMAXP2=NMAX2+6                                                    ACVI2152
      NMAX3=3*NBITS+10                                                  ACVI2153
      I=IPM                                                             ACVI2154
      J=NMAXP1-1                                                        ACVI2155
      NMIN=MOD(J,NEMPTY)                                                ACVI2156
      DO WHILE(J.LE.NMAX2 .AND. IABS(I).LE.NBITS)   ! MIDDLE            ACVI2157
        J=J+1                                                           ACVI2158
        IF(MOD(J,NEMPTY).NE.NMIN) THEN                                  ACVI2159
           NF=NFIRST+I                                                  ACVI2160
           IF(BTEST(IL,NF))THEN                                         ACVI2161
              IF(BTEST(IR,NF))THEN                                      ACVI2162
                 BIN(J:J)='I'                                           ACVI2163
              ELSE                                                      ACVI2164
                 BIN(J:J)='A'                                           ACVI2165
              ENDIF                                                     ACVI2166
           ELSE                                                         ACVI2167
              IF(BTEST(IR,NF))THEN                                      ACVI2168
                 BIN(J:J)='V'                                           ACVI2169
              ELSE                                                      ACVI2170
                 BIN(J:J)='-'                                           ACVI2171
              ENDIF                                                     ACVI2172
           ENDIF                                                        ACVI2173
           I=I+IPM                                                      ACVI2174
        ENDIF                                                           ACVI2175
      ENDDO                                                             ACVI2176
      I=IPM                                                             ACVI2177
      J=NMAXP2-1                                                        ACVI2178
      NMIN=MOD(J,NEMPTY)                                                ACVI2179
      DO WHILE(J.LE.NMAX3 .AND. IABS(I).LE.NBITS)   ! RIGHT             ACVI2180
        J=J+1                                                           ACVI2181
        IF(MOD(J,NEMPTY).NE.NMIN) THEN                                  ACVI2182
           IF(BTEST(IRGHT,NFIRST+I))THEN                                ACVI2183
              BIN(J:J)='1'                                              ACVI2184
           ELSE                                                         ACVI2185
              BIN(J:J)='0'                                              ACVI2186
           ENDIF                                                        ACVI2187
           I=I+IPM                                                      ACVI2188
        ENDIF                                                           ACVI2189
      ENDDO                                                             ACVI2190
C     Writing number in binary form.                                    ACVI2191
C     IF(TWRITE) WRITE(LFILE,'(T2,A)') BIN                              ACVI2192
      if( twrite ) write(lfile,'(t2,115a)') (bin(i:i), i=1, j)          ACVI2193
                  ENDIF                                                 ACVI2194
C100  FORMAT(' ',I23,2X,A)                                              ACVI2195
100   format(' ',i23,2x,115a)                                           ACVI2196
      RETURN                                                            ACVI2197
C ---*end of SETBIN*----------------------------------------------------ACVI2198
      END                                                               ACVI2199
C                                                                       ACVI2200
      SUBROUTINE ATTN(LITER,*)                                          ACVI2201
C     Warning because of gross error.                                   ACVI2202
      CHARACTER*(*) LITER                                               ACVI2203
      WRITE(6,'(//A,A)') ' ***** ATTENTION ',LITER                      ACVI2204
      RETURN 1                                                          ACVI2205
C     RETURN                                                            ACVI2206
C                                                                       ACVI2207
      ENTRY ERROR(LITER)                                                ACVI2208
C     Exit program because of gross error.                              ACVI2209
      WRITE(6,'(//A,A)') ' ***** FATAL ERROR ',LITER                    ACVI2210
      STOP                                                              ACVI2211
C ---*end of ATTN*------------------------------------------------------ACVI2212
      END                                                               ACVI2213
C                                                                       ACVI2214
      SUBROUTINE WRDATA(IOTYPE,IOFILE,NUMI,IRAY,NF)                     ACVI2215
C                                                                       ACVI2216
C     ORIGINAL: JPD ('PHDRYR.RESKE2.FORT(IODATA)')                      ACVI2217
C     MODIFIED: CB (2/15/91) FOR GENERAL PURPOSE                        ACVI2218
C                                                                       ACVI2219
CB    IMPLICIT REAL*8(A-H,O-Z)                                          ACVI2220
      DIMENSION IRAY(NF:1)                                              ACVI2221
      IF (IOTYPE.EQ.0) THEN                                             ACVI2222
         WRITE(IOFILE)NUMI                                              ACVI2223
         WRITE(IOFILE)(IRAY(I),I=NF,0)                                  ACVI2224
      ELSE                                                              ACVI2225
         READ(IOFILE)NUMI                                               ACVI2226
         READ(IOFILE)(IRAY(I),I=NF,0)                                   ACVI2227
      ENDIF                                                             ACVI2228
C     SPLIT: LRECL=1724 & BLKSIZE=32760 --> 430 4 BYTE SETS PER RECORD  ACVI2229
C            (BUFFER SIZE IS 32760 SO SPLIT INTO 19 RECORDS WITH        ACVI2230
C             4 BYTES PER RECORD PLUS 4 BYTES PER BLOCK OVERHEAD:       ACVI2231
C             19*(1720+4)+4=32760 FOR OPTIMUM BUFFER UTILIZATION)       ACVI2232
      NRUN=NUMI                                                         ACVI2233
      NREC=NRUN/430                                                     ACVI2234
      NBEG=1                                                            ACVI2235
      NEND=430                                                          ACVI2236
C     I/O FULL BLOCKS, INTEGER ARRAY                                    ACVI2237
      DO 10 NR=1,NREC                                                   ACVI2238
         IF (IOTYPE.EQ.0) THEN                                          ACVI2239
            WRITE(IOFILE)(IRAY(N),N=NBEG,NEND)                          ACVI2240
         ELSE                                                           ACVI2241
            READ(IOFILE)(IRAY(N),N=NBEG,NEND)                           ACVI2242
         ENDIF                                                          ACVI2243
         NBEG=NEND+1                                                    ACVI2244
10       NEND=NEND+430                                                  ACVI2245
C     I/O RESIDUAL                                                      ACVI2246
      IF (NBEG.LE.NRUN) THEN                                            ACVI2247
         IF (IOTYPE.EQ.0) THEN                                          ACVI2248
            WRITE(IOFILE)(IRAY(N),N=NBEG,NRUN)                          ACVI2249
         ELSE                                                           ACVI2250
            READ(IOFILE)(IRAY(N),N=NBEG,NRUN)                           ACVI2251
         ENDIF                                                          ACVI2252
      ENDIF                                                             ACVI2253
      RETURN                                                            ACVI2254
      END                                                               ACVI2255
C ----------------------------------------------------------------------ACVI2256
C                                                                       ACVI2257
C                       ***********************                         ACVI2258
C                       ***   WST PACKAGE   ***                         ACVI2259
C                       ***********************                         ACVI2260
C                                                                       ACVI2261
C                -------------------------------------                  ACVI2262
C                                                                       ACVI2263
C             *******************************************               ACVI2264
C             ***    WEIGHTED SEARCH TREE ROUTINES    ***               ACVI2265
C             ***                for                  ***               ACVI2266
C             ***  SCIENTIFIC (FORTRAN) APPLICATIONS  ***               ACVI2267
C             *******************************************               ACVI2268
C                                                                       ACVI2269
C ----------------------------------------------------------------------ACVI2270
C                                                                       ACVI2271
C Authors: Soon Park, C. Bahri, J. P. Draayer and S.-Q. Zheng           ACVI2272
C          Department of Physics (Computer Science)                     ACVI2273
C          Louisiana State University                                   ACVI2274
C          Baton Rouge LA                                               ACVI2275
C          USA 70803-4001                                               ACVI2276
C                                                                       ACVI2277
C                 BITNET:  PHDRYR @ LSUMVS or LSUVM                     ACVI2278
C                 TELEX:   559184                                       ACVI2279
C                 PHONE:   USA-504-388-2261                             ACVI2280
C                 FAX:     USA-504-388-5855                             ACVI2281
C                                                                       ACVI2282
C Version: 1.1    LSU (07/01/91)                                        ACVI2283
C Version: 2.1    LSU (04/01/94)                                        ACVI2284
C                                                                       ACVI2285
C ----------------------------------------------------------------------ACVI2286
C                                                                       ACVI2287
C Updates: 07/90 Original from a FORTRAN code written by Soon Park      ACVI2288
C          04/94 Inclusion of TMRG subroutine by C. Bahri               ACVI2289
C                                                                       ACVI2290
C ----------------------------------------------------------------------ACVI2291
C                                                                       ACVI2292
C General comments on the package:                                      ACVI2293
C                                                                       ACVI2294
C   This package is written in FORTRAN since most scientific programs   ACVI2295
C   require FORTRAN compatibility and for it the existing scientific    ACVI2296
C   subroutine libraries are the most extensive and efficient.          ACVI2297
C                                                                       ACVI2298
C   The tree is a linear array ID(-10:*) consisting of eleven integers, ACVI2299
C   ID(-10:0), that specifies the structure of the tree, see TSET for   ACVI2300
C   documentation on this, and node information starting with ID(1:1).  ACVI2301
C   The latter includes the key(s), data, balance factor and priority,  ACVI2302
C   as well as the left and right child pointers. See the documentation ACVI2303
C   on each subroutine for further details.                             ACVI2304
C                                                                       ACVI2305
C ----------------------------------------------------------------------ACVI2306
C                                                                       ACVI2307
C This numerical database package consists of 7 different subroutines:  ACVI2308
C                                                                       ACVI2309
C   1. TSET    -->  initializes a storage area for use as a binary tree ACVI2310
C       TSETLL -->  ... entry in TSET for a tree with a linked-list     ACVI2311
C       TSETLF -->  ... entry in TSET for a tree without a linked-list  ACVI2312
C   2. TCHK    -->  performs the search operation on the data structure ACVI2313
C   3. TADD    -->  add a new element to an existing WST data structure ACVI2314
C   4. TINS    -->  called by TADD to insert new node into existing WST ACVI2315
C   5. TDEL    -->  called by TADD to delete low priority node in a WST ACVI2316
C   6. TOUT    -->  generates output information on specific tree nodes ACVI2317
C   7. TMRG    -->  merges two trees that have the same structure       ACVI2318
C                                                                       ACVI2319
C The four routines TSET, TCHK, TADD and TOUT subject to user control   ACVI2320
C while TINS and TDEL are only used by TADD and not otherwise needed.   ACVI2321
C The type of information stored in the buffer may vary from applicaton ACVI2322
C to application. This can be achieved by editing TADD. For example,    ACVI2323
C if variable length integer data is to be stored, then BUFFER must be  ACVI2324
C replaced by, INTGER, an integer array throughout. Likewise if BUFFER  ACVI2325
C holds double precision or complex data, the statement that defines    ACVI2326
C BUFFER in TADD must be changed accordingly. The program can also be   ACVI2327
C used in other ways, for example, each node can refer to more than a   ACVI2328
C single buffer. If this is done, both TSET and TADD must be modified   ACVI2329
C in a rather obvious way to add the required multiple link-list and    ACVI2330
C buffer arrays.                                                        ACVI2331
C                                                                       ACVI2332
C ----------------------------------------------------------------------ACVI2333
C                                                                       ACVI2334
C                           **************                              ACVI2335
C                           ***  TSET  ***                              ACVI2336
C                           **************                              ACVI2337
C                                                                       ACVI2338
C The subroutine TSET must be called before inserting the first item    ACVI2339
C into the tree. This call fixes the first 11 values of the array ID:   ACVI2340
C                                                                       ACVI2341
C      ID(-10): Number of nodes currently in the tree    ==> ID(-10)    ACVI2342
C      ID(-9) : Maximum nodes in in tree                 ==> MXNODE     ACVI2343
C      ID(-8) : Parent node of ID(-7), see next entry    ==> NF         ACVI2344
C      ID(-7) : Node to be balanced                      ==> NA         ACVI2345
C      ID(-6) : Parent node pointer                      ==> NQ         ACVI2346
C      ID(-5) : Current node pointer                     ==> Determined ACVI2347
C      ID(-4) : Number of integers assigned the key      ==> NKEY       ACVI2348
C      ID(-3) : Number of integers for key and data      ==> NSUM       ACVI2349
C      ID(-2) : Position of the priority in a node       ==> NPR        ACVI2350
C      ID(-1) : Position of the next available node      ==> Determined ACVI2351
C      ID( 0) : Root pointer (-1 for empty tree)         ==> Determined ACVI2352
C                                                                       ACVI2353
C The subroutine TSET also initializes all entries of the array LLBUFF  ACVI2354
C where the link-list information for BUFFER is stored. In particular,  ACVI2355
C the first three which refer to the free space are assigned values:    ACVI2356
C                                                                       ACVI2357
C      LLBUFF(-2) : Tail of free space                                  ACVI2358
C      LLBUFF(-1) : Head of free space                                  ACVI2359
C      LLBUFF( 0) : Size of free space                                  ACVI2360
C                                                                       ACVI2361
C ----------------------------------------------------------------------ACVI2362
C                                                                       ACVI2363
      SUBROUTINE TSET                                                   ACVI2364
C                                                                       ACVI2365
      INTEGER ID(-10:*),LLBUFF(-2:*)                                    ACVI2366
C                                                                       ACVI2367
      ENTRY TSETLL(ID,MXNODE,NKEY,NDAT,LLBUFF,MXBUFF)                   ACVI2368
C                                                                       ACVI2369
C Initialize link-list array values                                     ACVI2370
C                                                                       ACVI2371
      DO 10 J=1,MXBUFF-1                                                ACVI2372
10       LLBUFF(J)=J+1                                                  ACVI2373
      LLBUFF(MXBUFF)=-1                                                 ACVI2374
C                                                                       ACVI2375
C Initialize free space pointer/counter                                 ACVI2376
C                                                                       ACVI2377
      LLBUFF(-2)=MXBUFF                                                 ACVI2378
      LLBUFF(-1)=1                                                      ACVI2379
      LLBUFF(0)=MXBUFF                                                  ACVI2380
C                                                                       ACVI2381
      ENTRY TSETLF(ID,MXNODE,NKEY,NDAT)                                 ACVI2382
C                                                                       ACVI2383
C Initialize WST parameters                                             ACVI2384
C                                                                       ACVI2385
      NSUM=NKEY+NDAT                                                    ACVI2386
      IMAX=MXNODE*(NSUM+3)                                              ACVI2387
      ID(-10)=0                                                         ACVI2388
      ID(-9)=MXNODE                                                     ACVI2389
      ID(-8)=-1                                                         ACVI2390
      ID(-7)=-1                                                         ACVI2391
      ID(-6)=-1                                                         ACVI2392
      ID(-5)=0                                                          ACVI2393
      ID(-4)=NKEY                                                       ACVI2394
      ID(-3)=NSUM                                                       ACVI2395
      ID(-2)=NSUM+1                                                     ACVI2396
      ID(-1)=0                                                          ACVI2397
      ID(0)=-1                                                          ACVI2398
      DO 30 I=1,IMAX                                                    ACVI2399
30       ID(I)=0                                                        ACVI2400
      RETURN                                                            ACVI2401
      END                                                               ACVI2402
C ----------------------------------------------------------------------ACVI2403
C                                                                       ACVI2404
C                            ************                               ACVI2405
C                            *** TCHK ***                               ACVI2406
C                            ************                               ACVI2407
C                                                                       ACVI2408
C The subroutine TCHK constructs a weighted search tree or locates a    ACVI2409
C specific node in a weighted search tree. The subroutine uses two      ACVI2410
C arrays, NEWKEY and ID. A call to TCHK(,,) generates a search of ID forACVI2411
C the key stored in NEWKEY. If the search is successful the priority of ACVI2412
C the node will be increased. If not successful, a normal return is     ACVI2413
C generated and the position, ID(-5) where a new node can be inserted   ACVI2414
C is given. After this subroutine is called, the position(s) where the  ACVI2415
C key(s) of the new node are to be assigned are ID(-5)+1, ID(-5)+2,     ACVI2416
C ID(-5)+3, ... and the positions of the location and size of the new   ACVI2417
C incoming data are ID(IPOS)+1 and ID(IPOS)+2, respectively, where      ACVI2418
C IPOS = ID(-5)+ID(-4) and ID(-4) is the number of integer words set    ACVI2419
C aside for the key.                                                    ACVI2420
C                                                                       ACVI2421
C   RETURN 1 --> Successful search (key already in the tree,            ACVI2422
C                ID(-5) is set to point to the located node).           ACVI2423
C                The priority of the located node is updated.           ACVI2424
C                                                                       ACVI2425
C ----------------------------------------------------------------------ACVI2426
C                                                                       ACVI2427
      SUBROUTINE TCHK(NEWKEY,ID,*)                                      ACVI2428
C                                                                       ACVI2429
      INTEGER NEWKEY(*),ID(-10:*)                                       ACVI2430
      LOGICAL FLAG                                                      ACVI2431
C                                                                       ACVI2432
C Integer data for a DEC system                                         ACVI2433
      DATA NF0,N2,N3/1073741823,-2147483648,-1073741824/                ACVI2434
C Hexadecimal data for a IBM system                                     ACVI2435
C                                                                       ACVI2436
C      DATA NF0,N2,N3/Z3FFFFFFF,Z80000000,ZC0000000/                     ACVI2437
      NKEY=ID(-4)                                                       ACVI2438
      NPR=ID(-2)                                                        ACVI2439
      NLC=NPR+1                                                         ACVI2440
      NRC=NPR+2                                                         ACVI2441
C                                                                       ACVI2442
C Special case (empty tree)                                             ACVI2443
C                                                                       ACVI2444
      IF (ID(0).EQ.-1) THEN                                             ACVI2445
         ID(-5)=0                                                       ACVI2446
         RETURN                                                         ACVI2447
      ENDIF                                                             ACVI2448
C                                                                       ACVI2449
C Normal case (non-empty tree)                                          ACVI2450
C                                                                       ACVI2451
      NF=-1                                                             ACVI2452
      NA=ID(0)                                                          ACVI2453
      NP=ID(0)                                                          ACVI2454
      NQ=-1                                                             ACVI2455
100   IF (NP.NE.-1) THEN                                                ACVI2456
C                                                                       ACVI2457
C Check if the balance factor of NP is 0 or not                         ACVI2458
C                                                                       ACVI2459
         IF (ID(NP+NPR).LT.0) THEN                                      ACVI2460
            NA=NP                                                       ACVI2461
            NF=NQ                                                       ACVI2462
         ENDIF                                                          ACVI2463
         DO 101 I=1,NKEY                                                ACVI2464
            IF (NEWKEY(I).LT.ID(NP+I)) THEN                             ACVI2465
               NQ=NP                                                    ACVI2466
               NO=ID(NP+NLC)                                            ACVI2467
               IF (NO.EQ.-1) THEN                                       ACVI2468
                  NP=-1                                                 ACVI2469
               ELSEIF (ID(NO+NRC).EQ.NP) THEN                           ACVI2470
                  IF (ISHFT(ID(NP+NPR),-30).EQ.2) THEN                  ACVI2471
                     NP=NO                                              ACVI2472
                  ELSE                                                  ACVI2473
                     NP=-1                                              ACVI2474
                  ENDIF                                                 ACVI2475
               ELSE                                                     ACVI2476
                  NP=NO                                                 ACVI2477
               ENDIF                                                    ACVI2478
               GOTO 100                                                 ACVI2479
            ELSEIF (NEWKEY(I).GT.ID(NP+I)) THEN                         ACVI2480
               NQ=NP                                                    ACVI2481
               NO=ID(NP+NLC)                                            ACVI2482
               IF (NO.EQ.-1) THEN                                       ACVI2483
                  NP=-1                                                 ACVI2484
               ELSEIF (ID(NO+NRC).EQ.NP) THEN                           ACVI2485
                  IF (ISHFT(ID(NP+NPR),-30).EQ.2) THEN                  ACVI2486
                     NP=-1                                              ACVI2487
                  ELSE                                                  ACVI2488
                     NP=NO                                              ACVI2489
                  ENDIF                                                 ACVI2490
               ELSE                                                     ACVI2491
                  NP=ID(NO+NRC)                                         ACVI2492
               ENDIF                                                    ACVI2493
               GOTO 100                                                 ACVI2494
            ENDIF                                                       ACVI2495
101      CONTINUE                                                       ACVI2496
C                                                                       ACVI2497
C If a match is found, increase the node priority, reconstruct the      ACVI2498
C linear array, set the pointer to the node location, and RETURN 1      ACVI2499
C                                                                       ACVI2500
         NA=NP                                                          ACVI2501
         FLAG=.TRUE.                                                    ACVI2502
         ITEM=IAND(ID(NP+NPR),NF0)                                      ACVI2503
         ITEMP=ITEM+ISHFT(ISHFT(ITEM,24),-16)                           ACVI2504
C                                                                       ACVI2505
C Saturation condition for the priority                                 ACVI2506
C                                                                       ACVI2507
         IF (ITEMP.LE.NF0) THEN                                         ACVI2508
            ID(NP+NPR)=ITEMP+IAND(ID(NP+NPR),N3)                        ACVI2509
            ITEMPR=ITEMP                                                ACVI2510
         ELSE                                                           ACVI2511
            ITEMPR=ITEM                                                 ACVI2512
         ENDIF                                                          ACVI2513
         NL=ID(-1)                                                      ACVI2514
102      NB=NA+NA+NRC                                                   ACVI2515
         IF (NB.LT.NL) THEN                                             ACVI2516
            NC=NB+NRC                                                   ACVI2517
            NBPR=IAND(ID(NB+NPR),NF0)                                   ACVI2518
            NCPR=IAND(ID(NC+NPR),NF0)                                   ACVI2519
            IF ((NBPR.LE.NCPR).OR.(NC.GE.NL)) THEN                      ACVI2520
               IF (NBPR.LT.ITEMPR) THEN                                 ACVI2521
                  IF (FLAG) THEN                                        ACVI2522
                     NARC=ID(NA+NRC)                                    ACVI2523
                     NALC=ID(NA+NLC)                                    ACVI2524
                     IF (NARC.NE.-1) THEN                               ACVI2525
                        IF (ID(NARC+NLC).EQ.NA) THEN                    ACVI2526
                           ID(NARC+NLC)=NL                              ACVI2527
                        ELSE                                            ACVI2528
                           IF (ID(ID(NARC+NRC)+NLC).EQ.NA) THEN         ACVI2529
                              ID(ID(NARC+NRC)+NLC)=NL                   ACVI2530
                           ELSE                                         ACVI2531
                              ID(ID(NARC+NLC)+NRC)=NL                   ACVI2532
                           ENDIF                                        ACVI2533
                        ENDIF                                           ACVI2534
                     ENDIF                                              ACVI2535
                     IF (NALC.NE.-1) THEN                               ACVI2536
                        IF (ID(NALC+NRC).EQ.NA) THEN                    ACVI2537
                           ID(NALC+NRC)=NL                              ACVI2538
                        ELSE                                            ACVI2539
                           ID(ID(NALC+NRC)+NRC)=NL                      ACVI2540
                        ENDIF                                           ACVI2541
                     ENDIF                                              ACVI2542
                     DO 105 I=1,NRC                                     ACVI2543
105                     ID(NL+I)=ID(NA+I)                               ACVI2544
                     FLAG=.FALSE.                                       ACVI2545
                  ENDIF                                                 ACVI2546
                  NBRC=ID(NB+NRC)                                       ACVI2547
                  NBLC=ID(NB+NLC)                                       ACVI2548
                  IF (NBRC.EQ.-1) THEN                                  ACVI2549
                     ID(0)=NA                                           ACVI2550
                  ELSE                                                  ACVI2551
                     IF (ID(NBRC+NLC).EQ.NB) THEN                       ACVI2552
                        ID(NBRC+NLC)=NA                                 ACVI2553
                     ELSE                                               ACVI2554
                        IF (ID(ID(NBRC+NRC)+NLC).EQ.NB) THEN            ACVI2555
                           ID(ID(NBRC+NRC)+NLC)=NA                      ACVI2556
                        ELSE                                            ACVI2557
                           ID(ID(NBRC+NLC)+NRC)=NA                      ACVI2558
                        ENDIF                                           ACVI2559
                     ENDIF                                              ACVI2560
                  ENDIF                                                 ACVI2561
                  IF (NBLC.NE.-1) THEN                                  ACVI2562
                     IF (ID(NBLC+NRC).EQ.NB) THEN                       ACVI2563
                        ID(NBLC+NRC)=NA                                 ACVI2564
                     ELSE                                               ACVI2565
                        ID(ID(NBLC+NRC)+NRC)=NA                         ACVI2566
                     ENDIF                                              ACVI2567
                  ENDIF                                                 ACVI2568
                  DO 109 I=1,NRC                                        ACVI2569
                     ID(NA+I)=ID(NB+I)                                  ACVI2570
109               CONTINUE                                              ACVI2571
                  NA=NB                                                 ACVI2572
                  GOTO 102                                              ACVI2573
               ENDIF                                                    ACVI2574
            ELSE                                                        ACVI2575
               IF (NCPR.LT.ITEMPR) THEN                                 ACVI2576
                  IF (FLAG) THEN                                        ACVI2577
                     NARC=ID(NA+NRC)                                    ACVI2578
                     NALC=ID(NA+NLC)                                    ACVI2579
                     IF (NARC.NE.-1) THEN                               ACVI2580
                        IF (ID(NARC+NLC).EQ.NA) THEN                    ACVI2581
                           ID(NARC+NLC)=NL                              ACVI2582
                        ELSE                                            ACVI2583
                           IF (ID(ID(NARC+NRC)+NLC).EQ.NA) THEN         ACVI2584
                              ID(ID(NARC+NRC)+NLC)=NL                   ACVI2585
                           ELSE                                         ACVI2586
                              ID(ID(NARC+NLC)+NRC)=NL                   ACVI2587
                           ENDIF                                        ACVI2588
                        ENDIF                                           ACVI2589
                     ENDIF                                              ACVI2590
                     IF (NALC.NE.-1) THEN                               ACVI2591
                        IF (ID(NALC+NRC).EQ.NA) THEN                    ACVI2592
                           ID(NALC+NRC)=NL                              ACVI2593
                        ELSE                                            ACVI2594
                           ID(ID(NALC+NRC)+NRC)=NL                      ACVI2595
                        ENDIF                                           ACVI2596
                     ENDIF                                              ACVI2597
                     DO 112 I=1,NRC                                     ACVI2598
112                     ID(NL+I)=ID(NA+I)                               ACVI2599
                     FLAG=.FALSE.                                       ACVI2600
                  ENDIF                                                 ACVI2601
                  NCRC=ID(NC+NRC)                                       ACVI2602
                  NCLC=ID(NC+NLC)                                       ACVI2603
                  IF (NCRC.EQ.-1) THEN                                  ACVI2604
                     ID(0)=NA                                           ACVI2605
                  ELSE                                                  ACVI2606
                     IF (ID(NCRC+NLC).EQ.NC) THEN                       ACVI2607
                        ID(NCRC+NLC)=NA                                 ACVI2608
                     ELSE                                               ACVI2609
                        IF (ID(ID(NCRC+NRC)+NLC).EQ.NC) THEN            ACVI2610
                           ID(ID(NCRC+NRC)+NLC)=NA                      ACVI2611
                        ELSE                                            ACVI2612
                           ID(ID(NCRC+NLC)+NRC)=NA                      ACVI2613
                        ENDIF                                           ACVI2614
                     ENDIF                                              ACVI2615
                  ENDIF                                                 ACVI2616
                  IF (NCLC.NE.-1) THEN                                  ACVI2617
                     IF (ID(NCLC+NRC).EQ.NC) THEN                       ACVI2618
                        ID(NCLC+NRC)=NA                                 ACVI2619
                     ELSE                                               ACVI2620
                        ID(ID(NCLC+NRC)+NRC)=NA                         ACVI2621
                     ENDIF                                              ACVI2622
                  ENDIF                                                 ACVI2623
                  DO 116 I=1,NRC                                        ACVI2624
                     ID(NA+I)=ID(NC+I)                                  ACVI2625
116               CONTINUE                                              ACVI2626
                  NA=NC                                                 ACVI2627
                  GOTO 102                                              ACVI2628
               ENDIF                                                    ACVI2629
            ENDIF                                                       ACVI2630
         ENDIF                                                          ACVI2631
118      IF (FLAG) GOTO 124                                             ACVI2632
         NLRC=ID(NL+NRC)                                                ACVI2633
         NLLC=ID(NL+NLC)                                                ACVI2634
         IF (NLRC.EQ.-1) THEN                                           ACVI2635
            ID(0)=NA                                                    ACVI2636
         ELSE                                                           ACVI2637
            IF (ID(NLRC+NLC).EQ.NL) THEN                                ACVI2638
               ID(NLRC+NLC)=NA                                          ACVI2639
            ELSE                                                        ACVI2640
               IF (ID(ID(NLRC+NRC)+NLC).EQ.NL) THEN                     ACVI2641
                  ID(ID(NLRC+NRC)+NLC)=NA                               ACVI2642
               ELSE                                                     ACVI2643
                  ID(ID(NLRC+NLC)+NRC)=NA                               ACVI2644
               ENDIF                                                    ACVI2645
            ENDIF                                                       ACVI2646
         ENDIF                                                          ACVI2647
         IF (NLLC.NE.-1) THEN                                           ACVI2648
            IF (ID(NLLC+NRC).EQ.NL) THEN                                ACVI2649
               ID(NLLC+NRC)=NA                                          ACVI2650
            ELSE                                                        ACVI2651
               ID(ID(NLLC+NRC)+NRC)=NA                                  ACVI2652
            ENDIF                                                       ACVI2653
         ENDIF                                                          ACVI2654
         DO 120 I=1,NRC                                                 ACVI2655
            ID(NA+I)=ID(NL+I)                                           ACVI2656
            ID(NL+I)=0                                                  ACVI2657
120      CONTINUE                                                       ACVI2658
124      ID(-5)=NA                                                      ACVI2659
         RETURN 1                                                       ACVI2660
      ENDIF                                                             ACVI2661
      ID(-8)=NF                                                         ACVI2662
      ID(-7)=NA                                                         ACVI2663
      ID(-6)=NQ                                                         ACVI2664
      ID(-5)=ID(-1)+NRC                                                 ACVI2665
      RETURN                                                            ACVI2666
      END                                                               ACVI2667
C ----------------------------------------------------------------------ACVI2668
C                                                                       ACVI2669
C                            ************                               ACVI2670
C                            *** TADD ***                               ACVI2671
C                            ************                               ACVI2672
C                                                                       ACVI2673
C The subroutine TADD allows the user to add elements to an existing    ACVI2674
C database provided space is available. If it is not and the priority   ACVI2675
C of the incoming element is greater than the lowest priority element   ACVI2676
C currently in the database that element is eliminated to make room for ACVI2677
C the new one. The last scenario is repeated as necessary to accommodateACVI2678
C the new element. The logic of the program is the following:           ACVI2679
C                                                                       ACVI2680
C Step 1) check if either the tree or buffer is full                    ACVI2681
C         if yes, check if the priority of incoming item is higher      ACVI2682
C                 than that of the lowest priority node in the tree     ACVI2683
C                 if yes, delete lowest priority node and goto Step 1)  ACVI2684
C                 else RETURN                                           ACVI2685
C         else    goto Step 2)                                          ACVI2686
C Step 2) insert the incoming item into the database                    ACVI2687
C                                                                       ACVI2688
C ----------------------------------------------------------------------ACVI2689
C                                                                       ACVI2690
      SUBROUTINE TADD(NEWKEY,NEWDAT,BULOAD,NOSIZE,NPBASE,ID,BUFFER,     ACVI2691
     &                LLBUFF)                                           ACVI2692
C                                                                       ACVI2693
      INTEGER NEWKEY(*),NEWDAT(*),ID(-10:*),LLBUFF(-2:*)                ACVI2694
      REAL    BULOAD(*),BUFFER(*)                     ! *** storage typeACVI2695
C     REAL*8  BULOAD(*),BUFFER(*)                     ! *** storage typeACVI2696
      DATA    NF0/1073741823/                                           ACVI2697
C                                                                       ACVI2698
C ...define the number of keys in a node                                ACVI2699
C                                                                       ACVI2700
      NKEY=ID(-4)                                                       ACVI2701
C                                                                       ACVI2702
C ...calculate current priority                                         ACVI2703
C                                                                       ACVI2704
      NPCURR=NPBASE+ISHFT(NPBASE,8)                                     ACVI2705
C                                                                       ACVI2706
C ...check to see if the tree and buffer can receive the new item       ACVI2707
C                                                                       ACVI2708
      IF ((ID(-10).LE.ID(-9)).AND.(LLBUFF(0).GE.NOSIZE)) THEN           ACVI2709
         GOTO 300                                                       ACVI2710
C                                                                       ACVI2711
C ...check if priority of new item is greater than lowest priority      ACVI2712
C                                                                       ACVI2713
      ELSE                                                              ACVI2714
100      IF (NPCURR.GT.IAND(ID(ID(-2)),NF0)) THEN                       ACVI2715
C                                                                       ACVI2716
C Rearrange LLBUFF array which holds free space list information        ACVI2717
C                                                                       ACVI2718
            IX=ID(NKEY+1)                                               ACVI2719
            IF (LLBUFF(0).EQ.0) THEN                                    ACVI2720
               LLBUFF(0)=ID(NKEY+2)                                     ACVI2721
               LLBUFF(-1)=IX                                            ACVI2722
            ELSE                                                        ACVI2723
               LLBUFF(0)=LLBUFF(0)+ID(NKEY+2)                           ACVI2724
               LLBUFF(LLBUFF(-2))=IX                                    ACVI2725
            ENDIF                                                       ACVI2726
            LLBUFF(-2)=IX                                               ACVI2727
            DO 200 I=1,ID(NKEY+2)-1                                     ACVI2728
               IX=LLBUFF(IX)                                            ACVI2729
200            LLBUFF(-2)=IX                                            ACVI2730
C                                                                       ACVI2731
            CALL TDEL(ID)                                               ACVI2732
C                                                                       ACVI2733
            IF (LLBUFF(0).LT.NOSIZE) GOTO 100                           ACVI2734
            CALL TCHK(NEWKEY,ID)                                        ACVI2735
            GOTO 300                                                    ACVI2736
C                                                                       ACVI2737
C ...doing nothing when new priority is less than lowest priority       ACVI2738
C                                                                       ACVI2739
         ELSE                                                           ACVI2740
            RETURN                                                      ACVI2741
         ENDIF                                                          ACVI2742
      ENDIF                                                             ACVI2743
300   CONTINUE                                                          ACVI2744
C                                                                       ACVI2745
C Load the key, data and its priority into the tree                     ACVI2746
C                                                                       ACVI2747
C ...last key position of the current node                              ACVI2748
C                                                                       ACVI2749
      IFIND=ID(-5)+NKEY                                                 ACVI2750
C                                                                       ACVI2751
C ...load integer data into tree after the last key                     ACVI2752
C                                                                       ACVI2753
      NEWDAT(1)=LLBUFF(-1)                                              ACVI2754
      NEWDAT(2)=NOSIZE                                                  ACVI2755
      NDAT=ID(-3)-ID(-4)                                                ACVI2756
      DO 400 IKK=1,NDAT                                                 ACVI2757
400   ID(IFIND+IKK)=NEWDAT(IKK)                                         ACVI2758
C                                                                       ACVI2759
C ...load real data into the BUFFER and update LLBUFF                   ACVI2760
C                                                                       ACVI2761
      DO 500 IKK=1,NOSIZE                                               ACVI2762
         BUFFER(LLBUFF(-1))=BULOAD(IKK)                                 ACVI2763
500      LLBUFF(-1)=LLBUFF(LLBUFF(-1))                                  ACVI2764
C                                                                       ACVI2765
C ...reduce size of free space                                          ACVI2766
C                                                                       ACVI2767
      LLBUFF(0)=LLBUFF(0)-NOSIZE                                        ACVI2768
C                                                                       ACVI2769
C ...load the priority                                                  ACVI2770
C                                                                       ACVI2771
      IFIND=ID(-5)+ID(-2)                                               ACVI2772
      ID(IFIND)=NPCURR                                                  ACVI2773
C                                                                       ACVI2774
C ...call TINS to complete (balance tree and restructure heap)          ACVI2775
C                                                                       ACVI2776
      CALL TINS(NEWKEY,ID)                                              ACVI2777
      RETURN                                                            ACVI2778
      END                                                               ACVI2779
C ----------------------------------------------------------------------ACVI2780
C                                                                       ACVI2781
C                            ************                               ACVI2782
C                            *** TINS ***                               ACVI2783
C                            ************                               ACVI2784
C                                                                       ACVI2785
C The subroutine TINS is used to insert a new item into a WST. A call   ACVI2786
C to TINS should follow a call to TCHK since TCHK tell where to insert  ACVI2787
C the new item. After the item is placed in the WST, TINS rebalances theACVI2788
C the tree and reconstruct the priority heap.                           ACVI2789
C                                                                       ACVI2790
C ----------------------------------------------------------------------ACVI2791
C                                                                       ACVI2792
      SUBROUTINE TINS(NEWKEY,ID)                                        ACVI2793
C                                                                       ACVI2794
      INTEGER NEWKEY(*),ID(-10:*)                                       ACVI2795
C                                                                       ACVI2796
C Integer data for a DEC system                                         ACVI2797
      DATA NF0,N2,N3/1073741823,-2147483648,-1073741824/                ACVI2798
C Hexadecimal data for a IBM  system                                    ACVI2799
C                                                                       ACVI2800
C     DATA NF0,N2,N3/Z3FFFFFFF,Z80000000,ZC0000000/                     ACVI2801
      ID(-10)=ID(-10)+1                                                 ACVI2802
      NKEY=ID(-4)                                                       ACVI2803
      NPR=ID(-2)                                                        ACVI2804
      NLC=NPR+1                                                         ACVI2805
      NRC=NPR+2                                                         ACVI2806
C                                                                       ACVI2807
C Special case (empty tree)                                             ACVI2808
C                                                                       ACVI2809
      IF (ID(0).EQ.-1) THEN                                             ACVI2810
         ID(-1)=NRC                                                     ACVI2811
         ID(0)=0                                                        ACVI2812
C                                                                       ACVI2813
C The first node position                                               ACVI2814
C                                                                       ACVI2815
         ID(NRC)=-1                                                     ACVI2816
         ID(NLC)=-1                                                     ACVI2817
         DO 10 I=1,NKEY                                                 ACVI2818
10          ID(I)=NEWKEY(I)                                             ACVI2819
         RETURN                                                         ACVI2820
      ENDIF                                                             ACVI2821
C                                                                       ACVI2822
C Normal case (non-empty tree)                                          ACVI2823
C                                                                       ACVI2824
      NY=ID(-5)                                                         ACVI2825
      NQ=ID(-6)                                                         ACVI2826
      NA=ID(-7)                                                         ACVI2827
      NF=ID(-8)                                                         ACVI2828
C                                                                       ACVI2829
C Set the pointer to indicate the position of the new node and load key ACVI2830
C                                                                       ACVI2831
      DO 130 I=1,NKEY                                                   ACVI2832
130      ID(NY+I)=NEWKEY(I)                                             ACVI2833
      ID(NY+NLC)=-1                                                     ACVI2834
      IF (ID(NQ+NLC).EQ.-1) THEN                                        ACVI2835
         ID(NY+NRC)=NQ                                                  ACVI2836
         ID(NQ+NLC)=NY                                                  ACVI2837
      ELSE                                                              ACVI2838
         DO 131 I=1,NKEY                                                ACVI2839
            IF (NEWKEY(I).LT.ID(NQ+I)) THEN                             ACVI2840
               ID(NY+NRC)=ID(NQ+NLC)                                    ACVI2841
               ID(NQ+NLC)=NY                                            ACVI2842
               ID(NQ+NPR)=IAND(ID(NQ+NPR),NF0)                          ACVI2843
               GOTO 999                                                 ACVI2844
            ELSEIF (NEWKEY(I).GT.ID(NQ+I)) THEN                         ACVI2845
               ID(NY+NRC)=NQ                                            ACVI2846
               ID(ID(NQ+NLC)+NRC)=NY                                    ACVI2847
               ID(NQ+NPR)=IAND(ID(NQ+NPR),NF0)                          ACVI2848
               GOTO 999                                                 ACVI2849
            ENDIF                                                       ACVI2850
131      CONTINUE                                                       ACVI2851
      ENDIF                                                             ACVI2852
      DO 161 I=1,NKEY                                                   ACVI2853
         IF (NEWKEY(I).LT.ID(NA+I)) THEN                                ACVI2854
            NP=ID(NA+NLC)                                               ACVI2855
            NB=NP                                                       ACVI2856
            IND=N2                                                      ACVI2857
            GOTO 200                                                    ACVI2858
         ELSEIF (NEWKEY(I).GT.ID(NA+I)) THEN                            ACVI2859
            NO=ID(NA+NLC)                                               ACVI2860
            IF (ID(NO+NRC).EQ.NA) THEN                                  ACVI2861
               NP=NO                                                    ACVI2862
            ELSE                                                        ACVI2863
               NP=ID(NO+NRC)                                            ACVI2864
            ENDIF                                                       ACVI2865
            NB=NP                                                       ACVI2866
            IND=N3                                                      ACVI2867
            GOTO 200                                                    ACVI2868
         ENDIF                                                          ACVI2869
161   CONTINUE                                                          ACVI2870
200   IF (NP.NE.NY) THEN                                                ACVI2871
         DO 211 I=1,NKEY                                                ACVI2872
            IF (NEWKEY(I).LT.ID(NP+I)) THEN                             ACVI2873
               ID(NP+NPR)=IAND(ID(NP+NPR),NF0)+N2                       ACVI2874
               NP=ID(NP+NLC)                                            ACVI2875
               GOTO 200                                                 ACVI2876
            ELSEIF (NEWKEY(I).GT.ID(NP+I)) THEN                         ACVI2877
               ID(NP+NPR)=IAND(ID(NP+NPR),NF0)+N3                       ACVI2878
               NO=ID(NP+NLC)                                            ACVI2879
               IF (ID(NO+NRC).EQ.NP) THEN                               ACVI2880
                  NP=NO                                                 ACVI2881
               ELSE                                                     ACVI2882
                  NP=ID(NO+NRC)                                         ACVI2883
               ENDIF                                                    ACVI2884
               GOTO 200                                                 ACVI2885
            ENDIF                                                       ACVI2886
211      CONTINUE                                                       ACVI2887
      ENDIF                                                             ACVI2888
      IF (ID(NA+NPR).GE.0) THEN                                         ACVI2889
         ID(NA+NPR)=IAND(ID(NA+NPR),NF0)+IND                            ACVI2890
         GOTO 999                                                       ACVI2891
      ELSEIF (IAND(ID(NA+NPR),N3).NE.IND) THEN                          ACVI2892
         ID(NA+NPR)=IAND(ID(NA+NPR),NF0)                                ACVI2893
         GOTO 999                                                       ACVI2894
      ENDIF                                                             ACVI2895
      IF (IND.EQ.N2) THEN                                               ACVI2896
         IF (ISHFT(ID(NB+NPR),-30).EQ.2) THEN                           ACVI2897
            IF (ID(ID(NA+NLC)+NRC).EQ.NA) THEN                          ACVI2898
               ID(NB+NRC)=ID(NA+NRC)                                    ACVI2899
               ID(ID(NB+NLC)+NRC)=NA                                    ACVI2900
               ID(NA+NLC)=-1                                            ACVI2901
               ID(NA+NRC)=NB                                            ACVI2902
            ELSE                                                        ACVI2903
               NC=ID(ID(NB+NLC)+NRC)                                    ACVI2904
               ID(NC+NRC)=ID(NB+NRC)                                    ACVI2905
               ID(NB+NRC)=ID(NA+NRC)                                    ACVI2906
               ID(ID(NB+NLC)+NRC)=NA                                    ACVI2907
               ID(NA+NLC)=NC                                            ACVI2908
               ID(NA+NRC)=NB                                            ACVI2909
            ENDIF                                                       ACVI2910
            ID(NA+NPR)=IAND(ID(NA+NPR),NF0)                             ACVI2911
            ID(NB+NPR)=IAND(ID(NB+NPR),NF0)                             ACVI2912
         ELSE                                                           ACVI2913
            IF (ID(ID(NB+NLC)+NRC).EQ.NB) THEN                          ACVI2914
               NC=ID(NB+NLC)                                            ACVI2915
            ELSE                                                        ACVI2916
               NC=ID(ID(NB+NLC)+NRC)                                    ACVI2917
            ENDIF                                                       ACVI2918
            IF (ID(NC+NPR).GE.0) THEN                                   ACVI2919
               ID(NC+NRC)=ID(NA+NRC)                                    ACVI2920
               ID(NA+NLC)=-1                                            ACVI2921
               ID(NA+NRC)=NC                                            ACVI2922
               ID(NB+NLC)=-1                                            ACVI2923
               ID(NC+NLC)=NB                                            ACVI2924
               ID(NA+NPR)=IAND(ID(NA+NPR),NF0)                          ACVI2925
               ID(NB+NPR)=IAND(ID(NB+NPR),NF0)                          ACVI2926
            ELSE                                                        ACVI2927
               IF (ID(ID(NC+NLC)+NRC).EQ.NC) THEN                       ACVI2928
                  IF (ISHFT(ID(NC+NPR),-30).EQ.2) THEN                  ACVI2929
                     ID(NA+NLC)=ID(NB+NRC)                              ACVI2930
                     ID(ID(NB+NLC)+NRC)=ID(NC+NLC)                      ACVI2931
                     ID(ID(NC+NLC)+NRC)=NB                              ACVI2932
                  ELSE                                                  ACVI2933
                     ID(NA+NLC)=ID(NC+NLC)                              ACVI2934
                     ID(ID(NA+NLC)+NRC)=ID(NB+NRC)                      ACVI2935
                     ID(ID(NB+NLC)+NRC)=NB                              ACVI2936
                  ENDIF                                                 ACVI2937
               ELSE                                                     ACVI2938
                  ID(NA+NLC)=ID(ID(NC+NLC)+NRC)                         ACVI2939
                  ID(ID(NA+NLC)+NRC)=ID(NB+NRC)                         ACVI2940
                  ID(ID(NB+NLC)+NRC)=ID(NC+NLC)                         ACVI2941
                  ID(ID(NC+NLC)+NRC)=NB                                 ACVI2942
               ENDIF                                                    ACVI2943
               ID(NC+NLC)=NB                                            ACVI2944
               ID(NC+NRC)=ID(NA+NRC)                                    ACVI2945
               ID(NA+NRC)=NC                                            ACVI2946
               ID(NB+NRC)=NA                                            ACVI2947
               IF (ISHFT(ID(NC+NPR),-30).EQ.2) THEN                     ACVI2948
                  ID(NA+NPR)=IAND(ID(NA+NPR),NF0)+N3                    ACVI2949
                  ID(NB+NPR)=IAND(ID(NB+NPR),NF0)                       ACVI2950
               ELSE                                                     ACVI2951
                  ID(NA+NPR)=IAND(ID(NA+NPR),NF0)                       ACVI2952
                  ID(NB+NPR)=IAND(ID(NB+NPR),NF0)+N2                    ACVI2953
               ENDIF                                                    ACVI2954
               ID(NC+NPR)=IAND(ID(NC+NPR),NF0)                          ACVI2955
            ENDIF                                                       ACVI2956
            NB=NC                                                       ACVI2957
         ENDIF                                                          ACVI2958
      ELSE                                                              ACVI2959
         IF (ISHFT(ID(NB+NPR),-30).EQ.3) THEN                           ACVI2960
            IF (ID(ID(NA+NLC)+NRC).EQ.NA) THEN                          ACVI2961
               ID(NA+NLC)=-1                                            ACVI2962
               ID(NB+NRC)=ID(NA+NRC)                                    ACVI2963
               ID(NA+NRC)=ID(NB+NLC)                                    ACVI2964
               ID(NB+NLC)=NA                                            ACVI2965
            ELSE                                                        ACVI2966
               NC=ID(NB+NLC)                                            ACVI2967
               ID(NB+NLC)=NA                                            ACVI2968
               ID(NB+NRC)=ID(NA+NRC)                                    ACVI2969
               ID(ID(NA+NLC)+NRC)=NC                                    ACVI2970
               ID(NA+NRC)=ID(NC+NRC)                                    ACVI2971
               ID(NC+NRC)=NA                                            ACVI2972
            ENDIF                                                       ACVI2973
            ID(NA+NPR)=IAND(ID(NA+NPR),NF0)                             ACVI2974
            ID(NB+NPR)=IAND(ID(NB+NPR),NF0)                             ACVI2975
         ELSE                                                           ACVI2976
            NC=ID(NB+NLC)                                               ACVI2977
            IF (ID(NC+NPR).GE.0) THEN                                   ACVI2978
               ID(NC+NRC)=ID(NA+NRC)                                    ACVI2979
               ID(NC+NLC)=NA                                            ACVI2980
               ID(NA+NLC)=-1                                            ACVI2981
               ID(NA+NRC)=NB                                            ACVI2982
               ID(NB+NLC)=-1                                            ACVI2983
               ID(NB+NRC)=NC                                            ACVI2984
               ID(NA+NPR)=IAND(ID(NA+NPR),NF0)                          ACVI2985
               ID(NB+NPR)=IAND(ID(NB+NPR),NF0)                          ACVI2986
            ELSE                                                        ACVI2987
               IF (ID(ID(NC+NLC)+NRC).EQ.NC) THEN                       ACVI2988
                  IF (ISHFT(ID(NC+NPR),-30).EQ.3) THEN                  ACVI2989
                     ID(ID(NA+NLC)+NRC)=NA                              ACVI2990
                     ID(ID(NC+NLC)+NRC)=ID(NC+NRC)                      ACVI2991
                     ID(NB+NLC)=ID(NC+NLC)                              ACVI2992
                  ELSE                                                  ACVI2993
                     ID(ID(NA+NLC)+NRC)=ID(NC+NLC)                      ACVI2994
                     ID(ID(NC+NLC)+NRC)=NA                              ACVI2995
                     ID(NB+NLC)=ID(NC+NRC)                              ACVI2996
                  ENDIF                                                 ACVI2997
               ELSE                                                     ACVI2998
                  ID(ID(NA+NLC)+NRC)=ID(NC+NLC)                         ACVI2999
                  ID(NB+NLC)=ID(ID(NC+NLC)+NRC)                         ACVI3000
                  ID(ID(NB+NLC)+NRC)=ID(NC+NRC)                         ACVI3001
                  ID(ID(NC+NLC)+NRC)=NA                                 ACVI3002
               ENDIF                                                    ACVI3003
               ID(NC+NLC)=NA                                            ACVI3004
               ID(NC+NRC)=ID(NA+NRC)                                    ACVI3005
               ID(NA+NRC)=NB                                            ACVI3006
               ID(NB+NRC)=NC                                            ACVI3007
               IF (ISHFT(ID(NC+NPR),-30).EQ.3) THEN                     ACVI3008
                  ID(NA+NPR)=IAND(ID(NA+NPR),NF0)+N2                    ACVI3009
                  ID(NB+NPR)=IAND(ID(NB+NPR),NF0)                       ACVI3010
               ELSE                                                     ACVI3011
                  ID(NA+NPR)=IAND(ID(NA+NPR),NF0)                       ACVI3012
                  ID(NB+NPR)=IAND(ID(NB+NPR),NF0)+N3                    ACVI3013
               ENDIF                                                    ACVI3014
               ID(NC+NPR)=IAND(ID(NC+NPR),NF0)                          ACVI3015
            ENDIF                                                       ACVI3016
            NB=NC                                                       ACVI3017
         ENDIF                                                          ACVI3018
      ENDIF                                                             ACVI3019
      IF (NF.EQ.-1) THEN                                                ACVI3020
         ID(0)=NB                                                       ACVI3021
      ELSEIF (NA.EQ.ID(NF+NLC)) THEN                                    ACVI3022
         ID(NF+NLC)=NB                                                  ACVI3023
      ELSE                                                              ACVI3024
         ID(ID(NF+NLC)+NRC)=NB                                          ACVI3025
      ENDIF                                                             ACVI3026
999   CONTINUE                                                          ACVI3027
C                                                                       ACVI3028
C Reconstruct the priority array                                        ACVI3029
C                                                                       ACVI3030
      NL=ID(-1)                                                         ACVI3031
      NRC2=NRC+NRC                                                      ACVI3032
1000  IF (MOD(NL,NRC2).EQ.0) THEN                                       ACVI3033
         NA=NL/2-NRC                                                    ACVI3034
      ELSE                                                              ACVI3035
         NA=(NL-NRC)/2                                                  ACVI3036
      ENDIF                                                             ACVI3037
      IF (IAND(ID(NA+NPR),NF0).GT.IAND(ID(NY+NPR),NF0)) THEN            ACVI3038
         NARC=ID(NA+NRC)                                                ACVI3039
         NALC=ID(NA+NLC)                                                ACVI3040
         IF (NARC.EQ.-1) THEN                                           ACVI3041
            ID(0)=NL                                                    ACVI3042
         ELSE                                                           ACVI3043
            IF (ID(NARC+NLC).EQ.NA) THEN                                ACVI3044
               ID(NARC+NLC)=NL                                          ACVI3045
            ELSE                                                        ACVI3046
               IF (ID(ID(NARC+NRC)+NLC).EQ.NA) THEN                     ACVI3047
                  ID(ID(NARC+NRC)+NLC)=NL                               ACVI3048
               ELSE                                                     ACVI3049
                  ID(ID(NARC+NLC)+NRC)=NL                               ACVI3050
               ENDIF                                                    ACVI3051
            ENDIF                                                       ACVI3052
         ENDIF                                                          ACVI3053
         IF (NALC.NE.-1) THEN                                           ACVI3054
            IF (ID(NALC+NRC).EQ.NA) THEN                                ACVI3055
               ID(NALC+NRC)=NL                                          ACVI3056
            ELSE                                                        ACVI3057
               ID(ID(NALC+NRC)+NRC)=NL                                  ACVI3058
            ENDIF                                                       ACVI3059
         ENDIF                                                          ACVI3060
         DO 1400 I=1,NRC                                                ACVI3061
            ID(NL+I)=ID(NA+I)                                           ACVI3062
1400     CONTINUE                                                       ACVI3063
         NL=NA                                                          ACVI3064
         IF (NL.GT.0) GOTO 1000                                         ACVI3065
      ENDIF                                                             ACVI3066
      NYRC=ID(NY+NRC)                                                   ACVI3067
      NYLC=ID(NY+NLC)                                                   ACVI3068
      IF (NYRC.EQ.-1) THEN                                              ACVI3069
         ID(0)=NL                                                       ACVI3070
      ELSE                                                              ACVI3071
         IF (ID(NYRC+NLC).EQ.NY) THEN                                   ACVI3072
            ID(NYRC+NLC)=NL                                             ACVI3073
         ELSE                                                           ACVI3074
            IF (ID(ID(NYRC+NRC)+NLC).EQ.NY) THEN                        ACVI3075
               ID(ID(NYRC+NRC)+NLC)=NL                                  ACVI3076
            ELSE                                                        ACVI3077
               ID(ID(NYRC+NLC)+NRC)=NL                                  ACVI3078
            ENDIF                                                       ACVI3079
         ENDIF                                                          ACVI3080
      ENDIF                                                             ACVI3081
      IF (NYLC.NE.-1) THEN                                              ACVI3082
         IF (ID(NYLC+NRC).EQ.NY) THEN                                   ACVI3083
            ID(NYLC+NRC)=NL                                             ACVI3084
         ELSE                                                           ACVI3085
            ID(ID(NYLC+NRC)+NRC)=NL                                     ACVI3086
         ENDIF                                                          ACVI3087
      ENDIF                                                             ACVI3088
      DO 7000 I=1,NRC                                                   ACVI3089
         ID(NL+I)=ID(NY+I)                                              ACVI3090
7000  CONTINUE                                                          ACVI3091
      ID(-1)=ID(-1)+NRC                                                 ACVI3092
      ID(-5)=NL                                                         ACVI3093
      RETURN                                                            ACVI3094
      END                                                               ACVI3095
C ----------------------------------------------------------------------ACVI3096
C                                                                       ACVI3097
C                            ************                               ACVI3098
C                            *** TDEL ***                               ACVI3099
C                            ************                               ACVI3100
C                                                                       ACVI3101
C The subroutine TDEL deletes a node from a weighted search tree,       ACVI3102
C keeping the binary tree balanced. The node deleted from the tree      ACVI3103
C is the first element of the linear array ID. Since the linear array   ACVI3104
C also has a heap structure, the root of a heap, which is the first     ACVI3105
C element of an array, is the lowest priority element (minheap). The    ACVI3106
C subroutine, TDEL, also manages the free space in the array BUFFER.    ACVI3107
C The pointer array LLBUFF is used for this purpose.                    ACVI3108
C                                                                       ACVI3109
C ----------------------------------------------------------------------ACVI3110
C                                                                       ACVI3111
      SUBROUTINE TDEL(ID)                                               ACVI3112
C                                                                       ACVI3113
      INTEGER ID(-10:*)                                                 ACVI3114
      LOGICAL FLAG                                                      ACVI3115
C                                                                       ACVI3116
C Integer data for a DEC system                                         ACVI3117
      DATA NF0,N2,N3/1073741823,-2147483648,-1073741824/                ACVI3118
C Hexadecimal data for a IBM system                                     ACVI3119
C                                                                       ACVI3120
C      DATA NF0,N2,N3/Z3FFFFFFF,Z80000000,ZC0000000/                     ACVI3121
C                                                                       ACVI3122
C Initialize some integer constants                                     ACVI3123
C                                                                       ACVI3124
      NKEY=ID(-4)                                                       ACVI3125
      NPR=ID(-2)                                                        ACVI3126
      NLC=NPR+1                                                         ACVI3127
      NRC=NPR+2                                                         ACVI3128
C                                                                       ACVI3129
C NA keeps track of most recent node with BF(0)                         ACVI3130
C NF is the parent of NA                                                ACVI3131
C NQ follows NP through the tree                                        ACVI3132
C                                                                       ACVI3133
      NF=-1                                                             ACVI3134
      NA=ID(0)                                                          ACVI3135
      NP=ID(0)                                                          ACVI3136
      NQ=-1                                                             ACVI3137
      NR=-1                                                             ACVI3138
100   IF (NP.NE.-1) THEN                                                ACVI3139
C                                                                       ACVI3140
C Looking for the last node with BF(0) that is not a leaf               ACVI3141
C                                                                       ACVI3142
         IF ((ID(NP+NPR).GE.0).AND.(ID(NP+NLC).NE.-1)) THEN             ACVI3143
            NA=NP                                                       ACVI3144
            NF=NQ                                                       ACVI3145
         ELSEIF (NQ.NE.-1) THEN                                         ACVI3146
            IF (FLAG) THEN                                              ACVI3147
               NO=ID(NQ+NLC)                                            ACVI3148
               IF (ID(NO+NRC).EQ.NQ) THEN                               ACVI3149
                  NRCHILD=NO                                            ACVI3150
               ELSE                                                     ACVI3151
                  NRCHILD=ID(NO+NRC)                                    ACVI3152
               ENDIF                                                    ACVI3153
               IF ((ISHFT(ID(NQ+NPR),-30).EQ.3)                         ACVI3154
     *             .AND.(ID(NRCHILD+NPR).GE.0)) THEN                    ACVI3155
                  NA=NQ                                                 ACVI3156
                  NF=NR                                                 ACVI3157
               ENDIF                                                    ACVI3158
            ELSE                                                        ACVI3159
               IF ((ISHFT(ID(NQ+NPR),-30).EQ.2)                         ACVI3160
     *             .AND.(ID(ID(NQ+NLC)+NPR).GE.0)) THEN                 ACVI3161
                  NA=NQ                                                 ACVI3162
                  NF=NR                                                 ACVI3163
               ENDIF                                                    ACVI3164
            ENDIF                                                       ACVI3165
         ENDIF                                                          ACVI3166
         DO 101 I=1,NKEY                                                ACVI3167
            IF (ID(I).LT.ID(NP+I)) THEN                                 ACVI3168
               NR=NQ                                                    ACVI3169
               NQ=NP                                                    ACVI3170
               NP=ID(NP+NLC)                                            ACVI3171
               FLAG=.TRUE.                                              ACVI3172
               GOTO 100                                                 ACVI3173
            ELSEIF (ID(I).GT.ID(NP+I)) THEN                             ACVI3174
               NR=NQ                                                    ACVI3175
               NQ=NP                                                    ACVI3176
               NO=ID(NP+NLC)                                            ACVI3177
               IF (ID(NO+NRC).EQ.NP) THEN                               ACVI3178
                  NP=NO                                                 ACVI3179
               ELSE                                                     ACVI3180
                  NP=ID(NO+NRC)                                         ACVI3181
               ENDIF                                                    ACVI3182
               FLAG=.FALSE.                                             ACVI3183
               GOTO 100                                                 ACVI3184
            ENDIF                                                       ACVI3185
101      CONTINUE                                                       ACVI3186
C                                                                       ACVI3187
C A match is found                                                      ACVI3188
C                                                                       ACVI3189
         GOTO 130                                                       ACVI3190
      ENDIF                                                             ACVI3191
C                                                                       ACVI3192
C A match is not found                                                  ACVI3193
C                                                                       ACVI3194
      WRITE(6,*) 'TREE IS EMPTY.'                                       ACVI3195
      GOTO 9999                                                         ACVI3196
130   CONTINUE                                                          ACVI3197
      ID(-10)=ID(-10)-1                                                 ACVI3198
C                                                                       ACVI3199
C Matched node does not have a child                                    ACVI3200
C                                                                       ACVI3201
      IF (ID(NP+NLC).EQ.-1) THEN                                        ACVI3202
         IF (NQ.EQ.-1) THEN                                             ACVI3203
            ID(0)=-1                                                    ACVI3204
            ID(-1)=0                                                    ACVI3205
            GOTO 9999                                                   ACVI3206
         ELSEIF (ID(NP+NRC).EQ.NQ) THEN                                 ACVI3207
            IF (ID(NQ+NLC).EQ.NP) THEN                                  ACVI3208
               ID(NQ+NLC)=-1                                            ACVI3209
            ELSE                                                        ACVI3210
               ID(ID(NQ+NLC)+NRC)=NQ                                    ACVI3211
            ENDIF                                                       ACVI3212
         ELSE                                                           ACVI3213
            ID(NQ+NLC)=ID(NP+NRC)                                       ACVI3214
         ENDIF                                                          ACVI3215
C                                                                       ACVI3216
C Matched node has only one child                                       ACVI3217
C                                                                       ACVI3218
      ELSEIF (ID(ID(NP+NLC)+NRC).EQ.NP) THEN                            ACVI3219
         NO=ID(NP+NLC)                                                  ACVI3220
         IF (ID(NP+NRC).EQ.-1) THEN                                     ACVI3221
            ID(0)=ID(NO)                                                ACVI3222
            ID(NO+NRC)=-1                                               ACVI3223
            GOTO 999                                                    ACVI3224
         ELSE                                                           ACVI3225
            IF (ID(NQ+NLC).EQ.NP) THEN                                  ACVI3226
               ID(NQ+NLC)=NO                                            ACVI3227
            ELSE                                                        ACVI3228
               ID(ID(NQ+NLC)+NRC)=NO                                    ACVI3229
            ENDIF                                                       ACVI3230
            ID(NO+NRC)=ID(NP+NRC)                                       ACVI3231
         ENDIF                                                          ACVI3232
C                                                                       ACVI3233
C Match node has both children                                          ACVI3234
C                                                                       ACVI3235
      ELSE                                                              ACVI3236
         NI=NQ                                                          ACVI3237
         NH=NP                                                          ACVI3238
         NG=ID(NP+NLC)                                                  ACVI3239
         IF (ID(NP+NPR).GE.0) THEN                                      ACVI3240
            NA=NH                                                       ACVI3241
            NF=NI                                                       ACVI3242
         ELSEIF ((ISHFT(ID(NP+NPR),-30).EQ.3).AND.                      ACVI3243
     *           (ID(ID(NG+NRC)+NPR).GE.0)) THEN                        ACVI3244
            NA=NH                                                       ACVI3245
            NF=NI                                                       ACVI3246
         ENDIF                                                          ACVI3247
         IF (ID(NG+NLC).EQ.-1) THEN                                     ACVI3248
            ID(ID(NG+NRC)+NRC)=NG                                       ACVI3249
            ID(NG+NLC)=ID(NG+NRC)                                       ACVI3250
            ID(NG+NRC)=ID(NP+NRC)                                       ACVI3251
         ELSEIF (ID(ID(NG+NLC)+NRC).EQ.NG) THEN                         ACVI3252
           IF (ISHFT(ID(NG+NPR),-30).EQ.2) THEN                         ACVI3253
               ID(ID(NG+NLC)+NRC)=ID(NG+NRC)                            ACVI3254
               ID(ID(NG+NRC)+NRC)=NG                                    ACVI3255
               ID(NG+NRC)=ID(NP+NRC)                                    ACVI3256
            ELSE                                                        ACVI3257
               NH=NG                                                    ACVI3258
               NG=ID(NG+NLC)                                            ACVI3259
               ID(NH+NLC)=-1                                            ACVI3260
               ID(ID(NH+NRC)+NRC)=NG                                    ACVI3261
               ID(NG+NLC)=NH                                            ACVI3262
               ID(NG+NRC)=ID(NP+NRC)                                    ACVI3263
            ENDIF                                                       ACVI3264
         ELSE                                                           ACVI3265
150         NI=NH                                                       ACVI3266
            NH=NG                                                       ACVI3267
            NG=ID(ID(NG+NLC)+NRC)                                       ACVI3268
            IF(ID(NH+NPR).GE.0) THEN                                    ACVI3269
               NA=NH                                                    ACVI3270
               NF=NI                                                    ACVI3271
            ELSEIF ((ISHFT(ID(NH+NPR),-30).EQ.2)                        ACVI3272
     *          .AND.(ID(ID(NH+NLC)+NPR).GE.0)) THEN                    ACVI3273
               NA=NH                                                    ACVI3274
               NF=NI                                                    ACVI3275
            ENDIF                                                       ACVI3276
            IF (ID(NG+NLC).EQ.-1) THEN                                  ACVI3277
               ID(ID(NH+NLC)+NRC)=NH                                    ACVI3278
            ELSEIF (ID(ID(NG+NLC)+NRC).EQ.NG) THEN                      ACVI3279
               IF (ISHFT(ID(NG+NPR),-30).EQ.2) THEN                     ACVI3280
                  ID(ID(NH+NLC)+NRC)=ID(NG+NLC)                         ACVI3281
                  ID(ID(NG+NLC)+NRC)=NH                                 ACVI3282
               ELSE                                                     ACVI3283
                  NH=NG                                                 ACVI3284
                  NG=ID(NG+NLC)                                         ACVI3285
                  ID(NH+NLC)=-1                                         ACVI3286
               ENDIF                                                    ACVI3287
            ELSE                                                        ACVI3288
               GOTO 150                                                 ACVI3289
            ENDIF                                                       ACVI3290
            ID(ID(ID(NP+NLC)+NRC)+NRC)=NG                               ACVI3291
            ID(NG+NLC)=ID(NP+NLC)                                       ACVI3292
            ID(NG+NRC)=ID(NP+NRC)                                       ACVI3293
         ENDIF                                                          ACVI3294
         IF (NQ.EQ.-1) THEN                                             ACVI3295
            ID(0)=NG                                                    ACVI3296
         ELSEIF (ID(NQ+NLC).EQ.NP) THEN                                 ACVI3297
            ID(NQ+NLC)=NG                                               ACVI3298
         ELSE                                                           ACVI3299
            ID(ID(NQ+NLC)+NRC)=NG                                       ACVI3300
         ENDIF                                                          ACVI3301
         ID(NG+NPR)=IAND(ID(NG+NPR),NF0)+IAND(ID(NP+NPR),N3)            ACVI3302
         IF (NA.EQ.NP) NA=NG                                            ACVI3303
         IF (NF.EQ.NP) NF=NG                                            ACVI3304
         DO 160 I=1,NKEY                                                ACVI3305
160         ID(I)=ID(NG+I)                                              ACVI3306
      ENDIF                                                             ACVI3307
C                                                                       ACVI3308
C Balance the tree                                                      ACVI3309
C                                                                       ACVI3310
500   CONTINUE                                                          ACVI3311
      IF (ID(NA+NLC).EQ.-1) THEN                                        ACVI3312
         ID(NA+NPR)=IAND(ID(NA+NPR),NF0)                                ACVI3313
         GOTO 999                                                       ACVI3314
      ENDIF                                                             ACVI3315
      DO 180 I=1,NKEY                                                   ACVI3316
         IF (ID(I).LT.ID(NA+I)) THEN                                    ACVI3317
            GOTO 190                                                    ACVI3318
         ELSEIF (ID(I).GT.ID(NA+I)) THEN                                ACVI3319
            NO=ID(NA+NLC)                                               ACVI3320
            IF (NO.EQ.-1) THEN                                          ACVI3321
               NP=-1                                                    ACVI3322
            ELSEIF (ID(NO+NRC).EQ.NA) THEN                              ACVI3323
               IF (ISHFT(ID(NA+NPR),-30).EQ.2) THEN                     ACVI3324
                  NP=-1                                                 ACVI3325
               ELSE                                                     ACVI3326
                  NP=NO                                                 ACVI3327
               ENDIF                                                    ACVI3328
            ELSE                                                        ACVI3329
               NP=ID(NO+NRC)                                            ACVI3330
            ENDIF                                                       ACVI3331
            NB=NO                                                       ACVI3332
            IND=N2                                                      ACVI3333
            GOTO 300                                                    ACVI3334
         ENDIF                                                          ACVI3335
180   CONTINUE                                                          ACVI3336
190   NO=ID(NA+NLC)                                                     ACVI3337
      IF (NO.EQ.-1) THEN                                                ACVI3338
         NP=-1                                                          ACVI3339
      ELSEIF (ID(NO+NRC).EQ.NA) THEN                                    ACVI3340
         IF (ISHFT(ID(NA+NPR),-30).EQ.2) THEN                           ACVI3341
            NP=NO                                                       ACVI3342
         ELSE                                                           ACVI3343
            NP=-1                                                       ACVI3344
         ENDIF                                                          ACVI3345
      ELSE                                                              ACVI3346
         NP=NO                                                          ACVI3347
      ENDIF                                                             ACVI3348
      IF (ID(NO+NRC).EQ.NA) THEN                                        ACVI3349
         NB=NO                                                          ACVI3350
      ELSE                                                              ACVI3351
         NB=ID(NO+NRC)                                                  ACVI3352
      ENDIF                                                             ACVI3353
      IND=N3                                                            ACVI3354
300   IF (ID(NA+NPR).GE.0) THEN                                         ACVI3355
         ID(NA+NPR)=IAND(ID(NA+NPR),NF0)+IND                            ACVI3356
         GOTO 600                                                       ACVI3357
      ELSEIF (IAND(ID(NA+NPR),N3).NE.IND) THEN                          ACVI3358
         ID(NA+NPR)=IAND(ID(NA+NPR),NF0)                                ACVI3359
         GOTO 600                                                       ACVI3360
      ENDIF                                                             ACVI3361
C                                                                       ACVI3362
C Rotate                                                                ACVI3363
C                                                                       ACVI3364
      IF (IND.EQ.N2) THEN                                               ACVI3365
         IF (ISHFT(ID(NB+NPR),-30).EQ.3) THEN                           ACVI3366
            IF (ID(ID(NB+NLC)+NRC).EQ.NB) THEN                          ACVI3367
               NC=ID(NB+NLC)                                            ACVI3368
            ELSE                                                        ACVI3369
               NC=ID(ID(NB+NLC)+NRC)                                    ACVI3370
            ENDIF                                                       ACVI3371
            IF (ID(NC+NPR).GE.0) THEN                                   ACVI3372
               IF (ID(NB+NRC).EQ.NA) THEN                               ACVI3373
                  ID(NC+NRC)=ID(NA+NRC)                                 ACVI3374
                  ID(NA+NLC)=-1                                         ACVI3375
                  ID(NA+NRC)=NC                                         ACVI3376
                  ID(NB+NLC)=-1                                         ACVI3377
                  ID(NC+NLC)=NB                                         ACVI3378
                  ID(NA+NPR)=IAND(ID(NA+NPR),NF0)                       ACVI3379
                  ID(NB+NPR)=IAND(ID(NB+NPR),NF0)                       ACVI3380
               ELSE                                                     ACVI3381
                  IDCL=ID(NC+NLC)                                       ACVI3382
                  ID(NA+NLC)=ID(IDCL+NRC)                               ACVI3383
                  ID(ID(NA+NLC)+NRC)=ID(NB+NRC)                         ACVI3384
                  ID(ID(NB+NLC)+NRC)=IDCL                               ACVI3385
                  ID(IDCL+NRC)=NB                                       ACVI3386
                  ID(NC+NLC)=NB                                         ACVI3387
                  ID(NC+NRC)=ID(NA+NRC)                                 ACVI3388
                  ID(NA+NRC)=NC                                         ACVI3389
                  ID(NB+NRC)=NA                                         ACVI3390
                  ID(NA+NPR)=IAND(ID(NA+NPR),NF0)                       ACVI3391
                  ID(NB+NPR)=IAND(ID(NB+NPR),NF0)                       ACVI3392
               ENDIF                                                    ACVI3393
            ELSE                                                        ACVI3394
               IDCL=ID(NC+NLC)                                          ACVI3395
               IF (ID(IDCL+NRC).EQ.NC) THEN                             ACVI3396
                  IF (ISHFT(ID(NC+NPR),-30).EQ.2) THEN                  ACVI3397
                     ID(NA+NLC)=ID(NB+NRC)                              ACVI3398
                     ID(ID(NB+NLC)+NRC)=IDCL                            ACVI3399
                     ID(IDCL+NRC)=NB                                    ACVI3400
                  ELSE                                                  ACVI3401
                     ID(NA+NLC)=IDCL                                    ACVI3402
                     ID(ID(NA+NLC)+NRC)=ID(NB+NRC)                      ACVI3403
                     ID(ID(NB+NLC)+NRC)=NB                              ACVI3404
                  ENDIF                                                 ACVI3405
               ELSE                                                     ACVI3406
                  ID(NA+NLC)=ID(IDCL+NRC)                               ACVI3407
                  ID(ID(NA+NLC)+NRC)=ID(NB+NRC)                         ACVI3408
                  ID(ID(NB+NLC)+NRC)=IDCL                               ACVI3409
                  ID(IDCL+NRC)=NB                                       ACVI3410
               ENDIF                                                    ACVI3411
               ID(NC+NLC)=NB                                            ACVI3412
               ID(NC+NRC)=ID(NA+NRC)                                    ACVI3413
               ID(NA+NRC)=NC                                            ACVI3414
               ID(NB+NRC)=NA                                            ACVI3415
               IF (ISHFT(ID(NC+NPR),-30).EQ.2) THEN                     ACVI3416
                  ID(NA+NPR)=IAND(ID(NA+NPR),NF0)+N3                    ACVI3417
                  ID(NB+NPR)=IAND(ID(NB+NPR),NF0)                       ACVI3418
               ELSE                                                     ACVI3419
                  ID(NA+NPR)=IAND(ID(NA+NPR),NF0)                       ACVI3420
                  ID(NB+NPR)=IAND(ID(NB+NPR),NF0)+N2                    ACVI3421
               ENDIF                                                    ACVI3422
               ID(NC+NPR)=IAND(ID(NC+NPR),NF0)                          ACVI3423
            ENDIF                                                       ACVI3424
            NB=NC                                                       ACVI3425
         ELSE                                                           ACVI3426
            IF (ID(NB+NRC).EQ.NA) THEN                                  ACVI3427
               IF (ID(NB+NPR).GE.0) THEN                                ACVI3428
                  NC=ID(ID(NB+NLC)+NRC)                                 ACVI3429
                  ID(NA+NLC)=NC                                         ACVI3430
                  ID(NC+NRC)=NA                                         ACVI3431
               ELSE                                                     ACVI3432
                  ID(NA+NLC)=-1                                         ACVI3433
               ENDIF                                                    ACVI3434
               ID(NB+NRC)=ID(NA+NRC)                                    ACVI3435
               ID(ID(NB+NLC)+NRC)=NA                                    ACVI3436
               ID(NA+NRC)=NB                                            ACVI3437
            ELSE                                                        ACVI3438
               NC=ID(ID(NB+NLC)+NRC)                                    ACVI3439
               ID(NC+NRC)=ID(NB+NRC)                                    ACVI3440
               ID(NB+NRC)=ID(NA+NRC)                                    ACVI3441
               ID(ID(NB+NLC)+NRC)=NA                                    ACVI3442
               ID(NA+NLC)=NC                                            ACVI3443
               ID(NA+NRC)=NB                                            ACVI3444
            ENDIF                                                       ACVI3445
            IF (ID(NB+NPR).GE.0) THEN                                   ACVI3446
               ID(NA+NPR)=IAND(ID(NA+NPR),NF0)+N2                       ACVI3447
               ID(NB+NPR)=IAND(ID(NB+NPR),NF0)+N3                       ACVI3448
            ELSE                                                        ACVI3449
               ID(NA+NPR)=IAND(ID(NA+NPR),NF0)                          ACVI3450
               ID(NB+NPR)=IAND(ID(NB+NPR),NF0)                          ACVI3451
            ENDIF                                                       ACVI3452
         ENDIF                                                          ACVI3453
      ELSE                                                              ACVI3454
         IF (ISHFT(ID(NB+NPR),-30).EQ.2) THEN                           ACVI3455
            NC=ID(NB+NLC)                                               ACVI3456
            IF (ID(NC+NPR).GE.0) THEN                                   ACVI3457
               IF (ID(NA+NLC).EQ.NB) THEN                               ACVI3458
                  ID(NC+NRC)=ID(NA+NRC)                                 ACVI3459
                  ID(NC+NLC)=NA                                         ACVI3460
                  ID(NA+NLC)=-1                                         ACVI3461
                  ID(NA+NRC)=NB                                         ACVI3462
                  ID(NB+NLC)=-1                                         ACVI3463
                  ID(NB+NRC)=NC                                         ACVI3464
                  ID(NA+NPR)=IAND(ID(NA+NPR),NF0)                       ACVI3465
                  ID(NB+NPR)=IAND(ID(NB+NPR),NF0)                       ACVI3466
               ELSE                                                     ACVI3467
                  IDCL=ID(NC+NLC)                                       ACVI3468
                  ID(ID(NA+NLC)+NRC)=IDCL                               ACVI3469
                  ID(NB+NLC)=ID(IDCL+NRC)                               ACVI3470
                  ID(ID(NB+NLC)+NRC)=ID(NC+NRC)                         ACVI3471
                  ID(IDCL+NRC)=NA                                       ACVI3472
                  ID(NC+NLC)=NA                                         ACVI3473
                  ID(NC+NRC)=ID(NA+NRC)                                 ACVI3474
                  ID(NA+NRC)=NB                                         ACVI3475
                  ID(NB+NRC)=NC                                         ACVI3476
                  ID(NA+NPR)=IAND(ID(NA+NPR),NF0)                       ACVI3477
                  ID(NB+NPR)=IAND(ID(NB+NPR),NF0)                       ACVI3478
               ENDIF                                                    ACVI3479
            ELSE                                                        ACVI3480
               IDCL=ID(NC+NLC)                                          ACVI3481
               IF (ID(IDCL+NRC).EQ.NC) THEN                             ACVI3482
                  IF (ISHFT(ID(NC+NPR),-30).EQ.3) THEN                  ACVI3483
                     ID(ID(NA+NLC)+NRC)=NA                              ACVI3484
                     ID(IDCL+NRC)=ID(NC+NRC)                            ACVI3485
                     ID(NB+NLC)=IDCL                                    ACVI3486
                  ELSE                                                  ACVI3487
                     ID(ID(NA+NLC)+NRC)=IDCL                            ACVI3488
                     ID(IDCL+NRC)=NA                                    ACVI3489
                     ID(NB+NLC)=ID(NC+NRC)                              ACVI3490
                  ENDIF                                                 ACVI3491
               ELSE                                                     ACVI3492
                  ID(ID(NA+NLC)+NRC)=IDCL                               ACVI3493
                  ID(NB+NLC)=ID(IDCL+NRC)                               ACVI3494
                  ID(ID(NB+NLC)+NRC)=ID(NC+NRC)                         ACVI3495
                  ID(IDCL+NRC)=NA                                       ACVI3496
               ENDIF                                                    ACVI3497
               ID(NC+NLC)=NA                                            ACVI3498
               ID(NC+NRC)=ID(NA+NRC)                                    ACVI3499
               ID(NA+NRC)=NB                                            ACVI3500
               ID(NB+NRC)=NC                                            ACVI3501
               IF (ISHFT(ID(NC+NPR),-30).EQ.3) THEN                     ACVI3502
                  ID(NA+NPR)=IAND(ID(NA+NPR),NF0)+N2                    ACVI3503
                  ID(NB+NPR)=IAND(ID(NB+NPR),NF0)                       ACVI3504
               ELSE                                                     ACVI3505
                  ID(NA+NPR)=IAND(ID(NA+NPR),NF0)                       ACVI3506
                  ID(NB+NPR)=IAND(ID(NB+NPR),NF0)+N3                    ACVI3507
               ENDIF                                                    ACVI3508
               ID(NC+NPR)=IAND(ID(NC+NPR),NF0)                          ACVI3509
            ENDIF                                                       ACVI3510
            NB=NC                                                       ACVI3511
         ELSE                                                           ACVI3512
            IF (ID(NA+NLC).EQ.NB) THEN                                  ACVI3513
               IF (ID(NB+NPR).GE.0) THEN                                ACVI3514
                  NC=ID(NB+NLC)                                         ACVI3515
                  ID(NB+NRC)=ID(NA+NRC)                                 ACVI3516
                  ID(NA+NLC)=NC                                         ACVI3517
                  ID(NA+NRC)=ID(NC+NRC)                                 ACVI3518
                  ID(NC+NRC)=NA                                         ACVI3519
               ELSE                                                     ACVI3520
                  ID(NA+NLC)=-1                                         ACVI3521
                  ID(NB+NRC)=ID(NA+NRC)                                 ACVI3522
                  ID(NA+NRC)=ID(NB+NLC)                                 ACVI3523
               ENDIF                                                    ACVI3524
               ID(NB+NLC)=NA                                            ACVI3525
            ELSE                                                        ACVI3526
               NC=ID(NB+NLC)                                            ACVI3527
               ID(NB+NLC)=NA                                            ACVI3528
               ID(NB+NRC)=ID(NA+NRC)                                    ACVI3529
               ID(ID(NA+NLC)+NRC)=NC                                    ACVI3530
               ID(NA+NRC)=ID(NC+NRC)                                    ACVI3531
               ID(NC+NRC)=NA                                            ACVI3532
            ENDIF                                                       ACVI3533
            IF (ID(NB+NPR).GE.0) THEN                                   ACVI3534
               ID(NA+NPR)=IAND(ID(NA+NPR),NF0)+N3                       ACVI3535
               ID(NB+NPR)=IAND(ID(NB+NPR),NF0)+N2                       ACVI3536
            ELSE                                                        ACVI3537
               ID(NA+NPR)=IAND(ID(NA+NPR),NF0)                          ACVI3538
               ID(NB+NPR)=IAND(ID(NB+NPR),NF0)                          ACVI3539
            ENDIF                                                       ACVI3540
         ENDIF                                                          ACVI3541
      ENDIF                                                             ACVI3542
      IF (NF.EQ.-1) THEN                                                ACVI3543
         ID(0)=NB                                                       ACVI3544
      ELSEIF (NA.EQ.ID(NF+NLC)) THEN                                    ACVI3545
         ID(NF+NLC)=NB                                                  ACVI3546
      ELSE                                                              ACVI3547
         ID(ID(NF+NLC)+NRC)=NB                                          ACVI3548
      ENDIF                                                             ACVI3549
600   NF=NA                                                             ACVI3550
      NA=NP                                                             ACVI3551
      IF (NP.NE.-1) GOTO 500                                            ACVI3552
999   CONTINUE                                                          ACVI3553
C                                                                       ACVI3554
C Reconstruct the priority array                                        ACVI3555
C                                                                       ACVI3556
      NA=0                                                              ACVI3557
      NL=ID(-1)-NRC                                                     ACVI3558
      NLPR=IAND(ID(NL+NPR),NF0)                                         ACVI3559
1000  NB=NA+NA+NRC                                                      ACVI3560
      IF (NB.LT.NL)THEN                                                 ACVI3561
         NC=NB+NRC                                                      ACVI3562
         NBPR=IAND(ID(NB+NPR),NF0)                                      ACVI3563
         NCPR=IAND(ID(NC+NPR),NF0)                                      ACVI3564
         IF ((NBPR.LE.NCPR).OR.(NC.GT.NL)) THEN                         ACVI3565
            IF (NBPR.LT.NLPR) THEN                                      ACVI3566
               NBRC=ID(NB+NRC)                                          ACVI3567
               NBLC=ID(NB+NLC)                                          ACVI3568
               IF (NBRC.EQ.-1) THEN                                     ACVI3569
                  ID(0)=NA                                              ACVI3570
               ELSE                                                     ACVI3571
                  IF (ID(NBRC+NLC).EQ.NB) THEN                          ACVI3572
                     ID(NBRC+NLC)=NA                                    ACVI3573
                  ELSE                                                  ACVI3574
                     IF (ID(ID(NBRC+NRC)+NLC).EQ.NB) THEN               ACVI3575
                        ID(ID(NBRC+NRC)+NLC)=NA                         ACVI3576
                     ELSE                                               ACVI3577
                        ID(ID(NBRC+NLC)+NRC)=NA                         ACVI3578
                     ENDIF                                              ACVI3579
                  ENDIF                                                 ACVI3580
               ENDIF                                                    ACVI3581
               IF (NBLC.NE.-1) THEN                                     ACVI3582
                  IF (ID(NBLC+NRC).EQ.NB) THEN                          ACVI3583
                     ID(NBLC+NRC)=NA                                    ACVI3584
                  ELSE                                                  ACVI3585
                     ID(ID(NBLC+NRC)+NRC)=NA                            ACVI3586
                  ENDIF                                                 ACVI3587
               ENDIF                                                    ACVI3588
               DO 1400 I=1,NRC                                          ACVI3589
                  ID(NA+I)=ID(NB+I)                                     ACVI3590
1400           CONTINUE                                                 ACVI3591
               NA=NB                                                    ACVI3592
               GOTO 1000                                                ACVI3593
            ENDIF                                                       ACVI3594
         ELSE                                                           ACVI3595
            IF (NCPR.LT.NLPR) THEN                                      ACVI3596
               NCRC=ID(NC+NRC)                                          ACVI3597
               NCLC=ID(NC+NLC)                                          ACVI3598
               IF (NCRC.EQ.-1) THEN                                     ACVI3599
                  ID(0)=NA                                              ACVI3600
               ELSE                                                     ACVI3601
                  IF (ID(NCRC+NLC).EQ.NC) THEN                          ACVI3602
                     ID(NCRC+NLC)=NA                                    ACVI3603
                  ELSE                                                  ACVI3604
                     IF (ID(ID(NCRC+NRC)+NLC).EQ.NC) THEN               ACVI3605
                        ID(ID(NCRC+NRC)+NLC)=NA                         ACVI3606
                     ELSE                                               ACVI3607
                        ID(ID(NCRC+NLC)+NRC)=NA                         ACVI3608
                     ENDIF                                              ACVI3609
                  ENDIF                                                 ACVI3610
               ENDIF                                                    ACVI3611
               IF (NCLC.NE.-1) THEN                                     ACVI3612
                  IF (ID(NCLC+NRC).EQ.NC) THEN                          ACVI3613
                     ID(NCLC+NRC)=NA                                    ACVI3614
                  ELSE                                                  ACVI3615
                     ID(ID(NCLC+NRC)+NRC)=NA                            ACVI3616
                  ENDIF                                                 ACVI3617
               ENDIF                                                    ACVI3618
               DO 2400 I=1,NRC                                          ACVI3619
                  ID(NA+I)=ID(NC+I)                                     ACVI3620
2400           CONTINUE                                                 ACVI3621
               NA=NC                                                    ACVI3622
               GOTO 1000                                                ACVI3623
            ENDIF                                                       ACVI3624
         ENDIF                                                          ACVI3625
      ENDIF                                                             ACVI3626
      NLRC=ID(NL+NRC)                                                   ACVI3627
      NLLC=ID(NL+NLC)                                                   ACVI3628
      IF (NLRC.EQ.-1) THEN                                              ACVI3629
         ID(0)=NA                                                       ACVI3630
      ELSE                                                              ACVI3631
         IF (ID(NLRC+NLC).EQ.NL) THEN                                   ACVI3632
            ID(NLRC+NLC)=NA                                             ACVI3633
         ELSE                                                           ACVI3634
            IF (ID(ID(NLRC+NRC)+NLC).EQ.NL) THEN                        ACVI3635
               ID(ID(NLRC+NRC)+NLC)=NA                                  ACVI3636
            ELSE                                                        ACVI3637
               ID(ID(NLRC+NLC)+NRC)=NA                                  ACVI3638
            ENDIF                                                       ACVI3639
         ENDIF                                                          ACVI3640
      ENDIF                                                             ACVI3641
      IF (NLLC.NE.-1) THEN                                              ACVI3642
         IF (ID(NLLC+NRC).EQ.NL) THEN                                   ACVI3643
            ID(NLLC+NRC)=NA                                             ACVI3644
         ELSE                                                           ACVI3645
            ID(ID(NLLC+NRC)+NRC)=NA                                     ACVI3646
         ENDIF                                                          ACVI3647
      ENDIF                                                             ACVI3648
      DO 7000 I=1,NRC                                                   ACVI3649
         ID(NA+I)=ID(NL+I)                                              ACVI3650
         ID(NL+I)=0                                                     ACVI3651
7000  CONTINUE                                                          ACVI3652
      ID(-1)=NL                                                         ACVI3653
9999  RETURN                                                            ACVI3654
      END                                                               ACVI3655
C ----------------------------------------------------------------------ACVI3656
C                                                                       ACVI3657
C                           **************                              ACVI3658
C                           ***  TOUT  ***                              ACVI3659
C                           **************                              ACVI3660
C                                                                       ACVI3661
C The subroutine TOUT traverses the weighted search tree in order to    ACVI3662
C specific nodes and sends the outputs to an output devise designated   ACVI3663
C by NFILE. The traversal is done in ascending order if NAD=1 and in    ACVI3664
C descending order otherwise. A description of the arguments follows:   ACVI3665
C                                                                       ACVI3666
C   NFILE = File number assigned to output device where the results are ACVI3667
C           to be written. In the special case NFILE = 0, TOUT does not ACVI3668
C           output anything, however after upon returning ID(-5) points ACVI3669
C           to the current node position.                               ACVI3670
C     NAD = 1 if the traversal is to be in ascending order, otherwise   ACVI3671
C           it is done in descending order.                             ACVI3672
C     MIN = Number of the first node that is to be retrieved.           ACVI3673
C     MAX = Number of the last node that is to be retrieved.            ACVI3674
C   NSTEP = Step size; nodes that are integer multiples of NSTEP beyond ACVI3675
C           MIN up to MAX will be retrieved.                            ACVI3676
C      ID = Name of the tree array that is to be traversed.             ACVI3677
C                                                                       ACVI3678
C An internal array called STACK is used to store information about the ACVI3679
C path traversed. It is dimensioned 32 since this is the maximum height ACVI3680
C the tree can have.  This limit is set by the fact that 2**32-1 is the ACVI3681
C largest integer pointer that can be used on a 32 bit machine.         ACVI3682
C                                                                       ACVI3683
C ----------------------------------------------------------------------ACVI3684
C                                                                       ACVI3685
      SUBROUTINE TOUT(NFILE,NAD,MIN,MAX,NSTEP,ID)                       ACVI3686
C                                                                       ACVI3687
      INTEGER ID(-10:*),STACK(32)                                       ACVI3688
      NSUM=ID(-3)                                                       ACVI3689
      NLC=NSUM+2                                                        ACVI3690
      NRC=NLC+1                                                         ACVI3691
      I=0                                                               ACVI3692
      M=0                                                               ACVI3693
      NUM=MIN                                                           ACVI3694
      NODE=ID(0)                                                        ACVI3695
      IF (NAD.EQ.1) THEN                                                ACVI3696
200      IF (NODE.NE.-1) THEN                                           ACVI3697
            I=I+1                                                       ACVI3698
            STACK(I)=NODE                                               ACVI3699
            NO=ID(NODE+NLC)                                             ACVI3700
            IF (NO.EQ.-1) THEN                                          ACVI3701
               NODE=-1                                                  ACVI3702
            ELSEIF (ID(NO+NRC).EQ.NODE) THEN                            ACVI3703
               IF (ISHFT(ID(NODE+NPR),-30).EQ.2) THEN                   ACVI3704
                  NODE=NO                                               ACVI3705
               ELSE                                                     ACVI3706
                  NODE=-1                                               ACVI3707
               ENDIF                                                    ACVI3708
            ELSE                                                        ACVI3709
               NODE=NO                                                  ACVI3710
            ENDIF                                                       ACVI3711
            GOTO 200                                                    ACVI3712
300         M=M+1                                                       ACVI3713
            ID(-5)=NODE                                                 ACVI3714
            IF ((M.EQ.NUM).OR.(M.EQ.MAX)) THEN                          ACVI3715
               IF (NFILE.GT.0) WRITE(NFILE,1000)                        ACVI3716
     *           M,(ID(NODE+II),II=1,NSUM)                              ACVI3717
               IF (M.EQ.MAX) GOTO 999                                   ACVI3718
               NUM=NUM+NSTEP                                            ACVI3719
            ENDIF                                                       ACVI3720
            NO=ID(NODE+NLC)                                             ACVI3721
            IF (NO.EQ.-1) THEN                                          ACVI3722
               NODE=-1                                                  ACVI3723
            ELSEIF (ID(NO+NRC).EQ.NODE) THEN                            ACVI3724
               IF (ISHFT(ID(NODE+NPR),-30).EQ.2) THEN                   ACVI3725
                  NODE=-1                                               ACVI3726
               ELSE                                                     ACVI3727
                  NODE=NO                                               ACVI3728
               ENDIF                                                    ACVI3729
            ELSE                                                        ACVI3730
               NODE=ID(NO+NRC)                                          ACVI3731
            ENDIF                                                       ACVI3732
            GOTO 200                                                    ACVI3733
         ENDIF                                                          ACVI3734
         IF (I.NE.0) THEN                                               ACVI3735
            NODE=STACK(I)                                               ACVI3736
            I=I-1                                                       ACVI3737
            GOTO 300                                                    ACVI3738
         ENDIF                                                          ACVI3739
      ELSE                                                              ACVI3740
400      IF (NODE.NE.-1) THEN                                           ACVI3741
            I=I+1                                                       ACVI3742
            STACK(I)=NODE                                               ACVI3743
            NO=ID(NODE+NLC)                                             ACVI3744
            IF (NO.EQ.-1) THEN                                          ACVI3745
               NODE=-1                                                  ACVI3746
            ELSEIF (ID(NO+NRC).EQ.NODE) THEN                            ACVI3747
               IF (ISHFT(ID(NODE+NPR),-30).EQ.2) THEN                   ACVI3748
                  NODE=-1                                               ACVI3749
               ELSE                                                     ACVI3750
                  NODE=NO                                               ACVI3751
               ENDIF                                                    ACVI3752
            ELSE                                                        ACVI3753
               NODE=ID(NO+NRC)                                          ACVI3754
            ENDIF                                                       ACVI3755
            GOTO 400                                                    ACVI3756
500         M=M+1                                                       ACVI3757
            ID(-5)=NODE                                                 ACVI3758
            IF ((M.EQ.NUM).OR.(M.EQ.MAX)) THEN                          ACVI3759
               IF (NFILE.GT.0) WRITE(NFILE,1000)                        ACVI3760
     *           M,(ID(NODE+II),II=1,NSUM)                              ACVI3761
               IF (M.EQ.MAX) GOTO 999                                   ACVI3762
               NUM=NUM+NSTEP                                            ACVI3763
            ENDIF                                                       ACVI3764
            NO=ID(NODE+NLC)                                             ACVI3765
            IF (NO.EQ.-1) THEN                                          ACVI3766
               NODE=-1                                                  ACVI3767
            ELSEIF (ID(NO+NRC).EQ.NODE) THEN                            ACVI3768
               IF (ISHFT(ID(NODE+NPR),-30).EQ.2) THEN                   ACVI3769
                  NODE=NO                                               ACVI3770
               ELSE                                                     ACVI3771
                  NODE=-1                                               ACVI3772
               ENDIF                                                    ACVI3773
            ELSE                                                        ACVI3774
               NODE=NO                                                  ACVI3775
            ENDIF                                                       ACVI3776
            GOTO 400                                                    ACVI3777
         ENDIF                                                          ACVI3778
         IF (I.NE.0) THEN                                               ACVI3779
            NODE=STACK(I)                                               ACVI3780
            I=I-1                                                       ACVI3781
            GOTO 500                                                    ACVI3782
         ENDIF                                                          ACVI3783
      ENDIF                                                             ACVI3784
      WRITE(NFILE,*)                                                    ACVI3785
      WRITE(NFILE,*)'  WARNING:  TOUT SEARCH HAS GONE OUT OF BOUNDS!'   ACVI3786
999   RETURN                                                            ACVI3787
1000  FORMAT(3X,I6,':',2X,14I5/(12X,14I5))                              ACVI3788
      END                                                               ACVI3789
C ----------------------------------------------------------------------ACVI3790
C                                                                       ACVI3791
C                            ************                               ACVI3792
C                            *** TMRG ***                               ACVI3793
C                            ************                               ACVI3794
C                                                                       ACVI3795
C The subroutine TMRG allows the user to merge two trees that have the  ACVI3796
C same structure. Elements of the first tree are stored in the second.  ACVI3797
C If the merge results in an overflow condition, notification is given  ACVI3798
C and the lowest priority item is eliminated in favor of the incoming   ACVI3799
C element with higher priority. If the trees have a different structure,ACVI3800
C the merge is automatically aborted (alternate return) after giving an ACVI3801
C appropriate warning message. If the element already exists in the     ACVI3802
C second tree, the priority is set to the maximum of the two.           ACVI3803
C                                                                       ACVI3804
C   RETURN 1 --> Unsuccessful match (the trees have different           ACVI3805
C                structures) so merge is aborted.                       ACVI3806
C                                                                       ACVI3807
C ----------------------------------------------------------------------ACVI3808
C                                                                       ACVI3809
      SUBROUTINE TMRG(ID1,BUFF1,LLBF1,ID2,BUFF2,LLBF2,*)                ACVI3810
C                                                                       ACVI3811
      INTEGER ID1(-10:*),LLBF1(-2:*),ID2(-10:*),LLBF2(-2:*)             ACVI3812
      REAL    BUFF1(*),BUFF2(*)                       ! *** storage typeACVI3813
C     REAL*8  BUFF1(*),BUFF2(*)                       ! *** storage typeACVI3814
      PARAMETER (MAXDAT=1000)                                           ACVI3815
      INTEGER NDAT(2)                                                   ACVI3816
      REAL    BULOAD(MAXDAT)                          ! *** temp storageACVI3817
C     REAL*8  BULOAD(MAXDAT)                          ! *** temp storageACVI3818
C     DATA    NF0/Z3FFFFFFF/                                            ACVI3819
      DATA    NF0/1073741823/                                           ACVI3820
C                                                                       ACVI3821
C Check the structure of the trees for compatibility                    ACVI3822
C                                                                       ACVI3823
      IF (ID1(-4).NE.ID2(-4).OR.ID1(-3).NE.ID2(-3)) THEN                ACVI3824
         WRITE(6,*) ' TREE 1 AND 2 HAVE DIFFERENT STRUCTURES'           ACVI3825
         RETURN 1                                                       ACVI3826
      ENDIF                                                             ACVI3827
      IF (ID1(-10)+ID2(-10).GE.ID2(-9).OR.LLBF1(-1).GE.LLBF2(0)) THEN   ACVI3828
         WRITE(6,*) ' SIZE OF THE INCOMING NODE IS LARGER THAN MAXIMUM,'ACVI3829
         WRITE(6,*) ' ONE OR MORE LOW PRIORITY NODES WILL BE DELETED.'  ACVI3830
      ENDIF                                                             ACVI3831
C                                                                       ACVI3832
C List out information from tree 1                                      ACVI3833
C                                                                       ACVI3834
      DO INODE1=1,ID1(-10)                                              ACVI3835
         CALL TOUT(0,0,INODE1,INODE1,1,ID1)                             ACVI3836
C                                                                       ACVI3837
C ...find the key and data indices in first tree                        ACVI3838
C                                                                       ACVI3839
         IPNODE=ID1(-5)                                                 ACVI3840
         IPKEY=IPNODE+1                                                 ACVI3841
         IPRIOR=ID1(IPNODE+ID1(-2))                                     ACVI3842
C                                                                       ACVI3843
C Check whether or not the incoming node is new                         ACVI3844
C                                                                       ACVI3845
         CALL TCHK(ID1(IPKEY),ID2,*100)                                 ACVI3846
C                                                                       ACVI3847
C ...find the data information from first tree                          ACVI3848
C                                                                       ACVI3849
         IPDATA=IPNODE+ID1(-4)+1                                        ACVI3850
         IFIND=ID1(IPDATA)                                              ACVI3851
         NOSIZE=ID1(IPDATA+1)                                           ACVI3852
C                                                                       ACVI3853
C ...abort if the size is too big                                       ACVI3854
C                                                                       ACVI3855
         IF (NOSIZE.GT.MAXDAT) THEN                                     ACVI3856
            WRITE(6,*) ' SIZE OF TEMPORARY BUFFER BULOAD IS TOO SMALL.' ACVI3857
            WRITE(6,*) ' INCREASE MAXDAT IN SUBROUTINE TMRG AND'        ACVI3858
            WRITE(6,*) ' CONTINUE.'                                     ACVI3859
            STOP                                                        ACVI3860
         ENDIF                                                          ACVI3861
C                                                                       ACVI3862
         DO IKK=1,NOSIZE                                                ACVI3863
            BULOAD(IKK)=BUFF1(IFIND)                                    ACVI3864
            IFIND=LLBF1(IFIND)                                          ACVI3865
         ENDDO                                                          ACVI3866
C                                                                       ACVI3867
C Store key, data, and information in second tree                       ACVI3868
C                                                                       ACVI3869
         CALL TADD(ID1(IPKEY),NDAT,BULOAD,NOSIZE,IPRIOR,ID2,BUFF2,LLBF2)ACVI3870
         GOTO 200                                                       ACVI3871
C                                                                       ACVI3872
C Update the priority if the incoming node is old                       ACVI3873
C                                                                       ACVI3874
100      IFIND=ID2(-5)+ID2(-2)                                          ACVI3875
         ID2(IFIND)=MAX0(ID2(IFIND),IPRIOR)                             ACVI3876
C                                                                       ACVI3877
200      CONTINUE                                                       ACVI3878
      ENDDO                                                             ACVI3879
      RETURN                                                            ACVI3880
      END                                                               ACVI3881
C********************                                                   ACVI3882
