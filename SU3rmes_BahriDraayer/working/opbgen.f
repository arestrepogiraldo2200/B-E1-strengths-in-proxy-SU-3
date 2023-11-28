      program opbgen                                                    ACVI3884
CB * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *ACVI3885
*                                                                      *ACVI3886
*                ***  Operator Basis Generator   ***                   *ACVI3887
*                            ** (opbgen) **                            *ACVI3888
*                                                                      *ACVI3889
CB * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *ACVI3890
C                                                                       ACVI3891
C Author:  Chairul Bahri                                                ACVI3892
C          Department of Physics and Astronomy                          ACVI3893
C          Louisiana State University                                   ACVI3894
C          Baton Rouge LA 70803 USA                                     ACVI3895
C                                                                       ACVI3896
C          e-mail address: phbahr @ lsuvm.sncc.lsu.edu (lsuvm.bitnet)   ACVI3897
C                           bahri @ rouge.phys.lsu.edu                  ACVI3898
C          phone: USA (504)-388-2261                                    ACVI3899
C                     (504)-388-6846                                    ACVI3900
C          fax:   USA (504)-388-5855                                    ACVI3901
C                                                                       ACVI3902
C ----------------------------------------------------------------------ACVI3903
C                                                                       ACVI3904
C Updates: 04/09/90: original ... IBM 3090/                             ACVI3905
C          11/14/92: IBM RS/6000-560.                                   ACVI3906
C ----------------------------------------------------------------------ACVI3907
C                                                                       ACVI3908
C General description: PROGRAM 2                                        ACVI3909
C                                                                       ACVI3910
C    Main program for generating all coupled-SU(3) tensor operator basisACVI3911
C    in the terms of creation and annihilation operators, a+'s and a's. ACVI3912
C                                                                       ACVI3913
C       T(lm mu) eps 2j 2m =                                            ACVI3914
C                                                                       ACVI3915
C          SUM <(lm1 mu1) ep1 2j1;(mu2 lm2) -ep2 2j2 || (lm mu) eps 2j> ACVI3916
C                                                                       ACVI3917
C              < 2j1 2m1;2j2 -2m2 | 2j 2m >                             ACVI3918
C                  (p-r)               (lm1 mu1)            (lm2 mu2)   ACVI3919
C              (-1)      (a+ a+ ... a+)          (a a ... a)            ACVI3920
C                                                                       ACVI3921
C ----------------------------------------------------------------------ACVI3922
C  Reference:                                                           ACVI3923
C    C. Bahri, research notes.                                          ACVI3924
C    D. Braunschweig, Comp.Phys.Comp. 14, 109 (1978)                    ACVI3925
C    M.F.O'Reilly, J.Math.Phys. 23, 2033 (1982)                         ACVI3926
C ----------------------------------------------------------------------ACVI3927
      implicit real*8(d)                                                ACVI3928
C                                                                       ACVI3929
      parameter ( IOPBFILE = 9 )   ! unformatted opb file               ACVI3930
      character*20 cofile          ! default opb filename               ACVI3931
      parameter ( MXWCOF = 100000, ! max dimension of dwsu3             ACVI3932
     1            MXNUMB =   1000, ! max number of (lm mu) pairs lmmua  ACVI3933
     2            MXTREE = 500000 )! tree size                          ACVI3934
      common / DSU3OP /                                                 ACVI3935
     1         dwsu3( MXWCOF )     ! coefficients for SU(3) operators   ACVI3936
      common / TREES /                                                  ACVI3937
     1         mltree(-10:MXTREE ),! labels for uncoupled tensors       ACVI3938
     2         lbtree(-10:MXTREE ),! labels for coupled tensors (saved) ACVI3939
     3         lctree(-10:MXTREE ) ! labels for uncoupled tensors (svd) ACVI3940
      dimension iwsu3( 2*MXWCOF )                                       ACVI3941
      equivalence ( iwsu3, dwsu3 )                                      ACVI3942
C                                                                       ACVI3943
CW**********************************************************************ACVI3944
CW                                                                      ACVI3945
      write(6,'(/a,a/)') ' ENTERING *** OPERATOR BASIS GENERATOR',      ACVI3946
     1                   ' (OPBGEN) ***'                                ACVI3947
CW                                                                      ACVI3948
CW**********************************************************************ACVI3949
C                                                                       ACVI3950
      call heading( IOPBFILE, cofile, lfile )                           ACVI3951
      if(lfile.gt.0 .and. lfile.ne.6) open( unit = lfile,               ACVI3952
     1                                      file = 'opb.log' )          ACVI3953
      open( unit = IOPBFILE, file = cofile, form = 'unformatted',       ACVI3954
     1      status = 'new' )                                            ACVI3955
      call opb( lfile, mltree, lbtree, lctree, MXTREE,                  ACVI3956
     1          dwsu3, MXWCOF, mxwsu3, ieta )                           ACVI3957
C     Dump the required data to external file/tree.                     ACVI3958
      iloc  = lbtree(-10)                                               ACVI3959
      ilocc = lctree(-10)                                               ACVI3960
      write(6,1000) iloc, ilocc, mxwsu3                                 ACVI3961
1000  format(/' TREE STATISTICS: lbtree NODES = ',i10/                  ACVI3962
     1        '                  lctree       = ',i10/                  ACVI3963
     2        ' DATA STATISTICS:  dwsu3       = ',i10/)                 ACVI3964
      write( IOPBFILE ) ieta       ! signature                          ACVI3965
      iloc  = iloc *  ( lbtree(-3 ) + 3 )                               ACVI3966
      ilocc = ilocc * ( lctree(-3 ) + 3 )                               ACVI3967
      mxwsu3 = 2 * mxwsu3                                               ACVI3968
      call wrdata(0, IOPBFILE, iloc,  lbtree, -10 )                     ACVI3969
      call wrdata(0, IOPBFILE, ilocc, lctree, -10 )                     ACVI3970
      call wrdata(0, IOPBFILE, mxwsu3, dwsu3,   1 )                     ACVI3971
      close( unit = IOPBFILE, status = 'keep' )                         ACVI3972
      if(lfile.gt.0 .and. lfile.ne.6) close( unit = lfile,              ACVI3973
     1                                       status = 'keep' )          ACVI3974
      stop                                                              ACVI3975
C ---*end of opbgen*----------------------------------------------------ACVI3976
      end                                                               ACVI3977
      subroutine heading( iopbfile, cofile, lfile )                     ACVI3978
C                                                                       ACVI3979
C     Read file structure of opb. (Can be omitted.)                     ACVI3980
C                                                                       ACVI3981
      implicit logical(t)                                               ACVI3982
      character*(*) cofile         ! default opb filename               ACVI3983
C                   lfile          ! log file for intermediate results  ACVI3984
      character*20  ctemp          ! temporary filename                 ACVI3985
C                                                                       ACVI3986
1     format(a)                                                         ACVI3987
2     format(' =>',10i5)                                                ACVI3988
3     format(' =>',4x,a)                                                ACVI3989
      lfile = 0                                                         ACVI3990
C                                                                       ACVI3991
      write(6,1)         ' Enter output file name: '                    ACVI3992
      read (5,'(18a4)')  ctemp                                          ACVI3993
      cofile = 'opbw.' // ctemp(1:14)                                   ACVI3994
      write(6,3)         cofile                                         ACVI3995
      write(6,*)                                                        ACVI3996
      write(6,1)         ' Enter log file unit number:'                 ACVI3997
      write(6,1)         '       0: none'                               ACVI3998
      write(6,1)         '     >=1: some intermediate results'          ACVI3999
      read (5,*,end=999) lfile                                          ACVI4000
      write(6,2)         lfile                                          ACVI4001
      write(6,*)                                                        ACVI4002
C                                                                       ACVI4003
      lfile = iabs( lfile )                                             ACVI4004
      if( lfile.eq.5 .or. lfile.eq.iopbfile )                           ACVI4005
     1  call error(' HEADING:The file unit has been assigned.')         ACVI4006
C                                                                       ACVI4007
      return                                                            ACVI4008
999   write(6,1) ' I quit. Please enter numbers.'                       ACVI4009
      stop                                                              ACVI4010
C----*end of heading*-------------------------------------------------- ACVI4011
      end                                                               ACVI4012
CB * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *ACVI4013
*                                                                      *ACVI4014
*                ***  OPERATOR BASIS GENERATOR ***                     *ACVI4015
*                           ** (OPB) **                                *ACVI4016
*                                                                      *ACVI4017
CB * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *ACVI4018
C                                                                       ACVI4019
C Author:  Chairul Bahri                                                ACVI4020
C          Department of Physics and Astronomy                          ACVI4021
C          Louisiana State University                                   ACVI4022
C          Baton Rouge LA 70803 USA                                     ACVI4023
C                                                                       ACVI4024
C          e-mail: phbahr @ lsuvm.sncc.lsu.edu (lsuvm.bitnet)           ACVI4025
C                   bahri @ rouge.phys.lsu.edu                          ACVI4026
C          phone: USA (504)-388-2261                                    ACVI4027
C                     (504)-388-6846                                    ACVI4028
C          fax:   USA (504)-388-5855                                    ACVI4029
C                                                                       ACVI4030
C ----------------------------------------------------------------------ACVI4031
C                                                                       ACVI4032
C Updates: 04/09/90: original ... IBM 3090/                             ACVI4033
C          11/14/92: IBM RS/6000-560.                                   ACVI4034
C ----------------------------------------------------------------------ACVI4035
C                                                                       ACVI4036
C General description: SUB-PROGRAM 2                                    ACVI4037
C                                                                       ACVI4038
C    Main subprogram for generating all coupled-SU(3) tensor operator   ACVI4039
C    in the terms of creation and annihilation operators, a+'s and a's. ACVI4040
C                                                                       ACVI4041
C       T(lm mu) eps 2j 2m =                                            ACVI4042
C                                                                       ACVI4043
C          SUM <(lm1 mu1) ep1 2j1;(mu2 lm2) -ep2 2j2 || (lm mu) eps 2j> ACVI4044
C                                                                       ACVI4045
C              < 2j1 2m1;2j2 -2m2 | 2j 2m >                             ACVI4046
C                  (p-r)               (lm1 mu1)            (lm2 mu2)   ACVI4047
C              (-1)      (a+ a+ ... a+)          (a a ... a)            ACVI4048
C                                                                       ACVI4049
C ----------------------------------------------------------------------ACVI4050
C  Reference:                                                           ACVI4051
C    C. Bahri, research notes.                                          ACVI4052
C    D. Braunschweig, Comp.Phys.Comp. 14, 109 (1978)                    ACVI4053
C    M.F.O'Reilly, J.Math.Phys. 23, 2033 (1982)                         ACVI4054
C ----------------------------------------------------------------------ACVI4055
C                                                                       ACVI4056
      subroutine opb( logfile, mltree, lbtree, lctree, mxtree,          ACVI4057
     ]                dwigsu3, mxwcof, nx, neta )                       ACVI4058
      implicit real*8(d), logical(t)                                    ACVI4059
      character*4 cflag, cphase                                         ACVI4060
C                                                                       ACVI4061
C               logfile            ! log file for intermediate results  ACVI4062
      dimension mltree(-10:* ),    ! binary tree for uncoupled tensors  ACVI4063
     1          lbtree(-10:* ),    ! binary tree for coupled tensors    ACVI4064
     2          lctree(-10:* ),    ! binary tree for uncoupled tensors  ACVI4065
     3          dwigsu3( * )       ! coefficients of SU(3) couplings    ACVI4066
C               mxtree             ! max dimension of trees             ACVI4067
C               mxwcof             ! max dimension of dwigsu3           ACVI4068
C               neta               ! oscillator shell number            ACVI4069
      common / OPBCON / lfile,     ! logfile                            ACVI4070
     1                  ieta,      ! oscillator shell number            ACVI4071
     2                  nbits,     ! bit length for # levels in shell   ACVI4072
     3                  nbitx,     ! significant bits                   ACVI4073
     4                  cflag,     ! flag for new SU(3) irreps          ACVI4074
     5                  cphase,    ! SU(3) phase                        ACVI4075
     6                  tfile      ! yes, write to log file             ACVI4076
      common / INDICS / ix,        ! mltree index for dwigsu3           ACVI4077
     1                  iy,        ! lbtree index for dwigsu3           ACVI4078
     2                  maxy       ! max index for dwigsu3              ACVI4079
      parameter ( MXCOEF = 100000, ! max dimension of dwigsu3           ACVI4080
     1            MXNUMB =   1000 )! max number of (lm mu) pairs lmmua  ACVI4081
      common / ISU3OP /                                                 ACVI4082
     1         next( MXCOEF ),     ! linked-list index for mltree       ACVI4083
     2         levels( MXCOEF ),   ! U(3) levels                        ACVI4084
     3         nxy( MXCOEF )       ! (nx,ny) values for mltree          ACVI4085
      common / SU3PTR /                                                 ACVI4086
     1         lmmua( 2, MXNUMB ), ! array of packed (lm,mu) of SU(3)   ACVI4087
     2         iptr( MXNUMB )      ! level pointers for generations     ACVI4088
C     Packing functions                                                 ACVI4089
      data nbit8 / 255 /           ! zff                                ACVI4090
      ipack(i,j,k,l) = ior(l,ishft(ior(k,ishft(ior(j,ishft(i,           ACVI4091
     >                 8)),8)),8))                                      ACVI4092
      iupac(index,iover,ibits) = iand( ishft(index,iover),ibits )       ACVI4093
C                                                                       ACVI4094
      lfile = logfile                                                   ACVI4095
      tfile = lfile .gt. 0                                              ACVI4096
      if( mxwcof .lt. MXCOEF )                                          ACVI4097
     1  call attn(' OPB: Less dimension for coefficients!')             ACVI4098
      cflag  = '    '                                                   ACVI4099
      cphase = '    '                                                   ACVI4100
C                                                                       ACVI4101
      write(6,'(a)') ' Oscillator shell number (ieta)'                  ACVI4102
      read (5,*) ieta                                                   ACVI4103
      write(6,'(a,i5)') ' =>', ieta                                     ACVI4104
      write(6,*)                                                        ACVI4105
      neta = ieta                                                       ACVI4106
      if( tfile ) write(lfile,'(a,i4)')' Oscillator shell number =',ietaACVI4107
      write(6,'(a)') ' Operator structure (#a+,#a)'                     ACVI4108
      read (5,*) n1, n2                                                 ACVI4109
      write(6,'(a,2i5)') ' =>', n1, n2                                  ACVI4110
      write(6,*)                                                        ACVI4111
      nx = max0( n1, n2 )                                               ACVI4112
      tactive = ieta .le. 2 .and. nx .ge. 2                             ACVI4113
C                                                                       ACVI4114
      call blocks                         ! initialize binomials        ACVI4115
      nlevel = (ieta+1) * (ieta+2) / 2                                  ACVI4116
      if( MXCOEF .lt. nlevel ) call error(' OPB: Too few MXCOEF.')      ACVI4117
C                                                                       ACVI4118
C     Initialize trees and data.                                        ACVI4119
C                                                                       ACVI4120
      inodes = mxtree / 7 - 1                                           ACVI4121
      call tsetlf( mltree, inodes, 2, 2 )   ! 2 keys, 2 data            ACVI4122
      inodes = mxtree / 9 - 1                                           ACVI4123
      call tsetlf( lbtree, inodes, 4, 2 )   ! 4 keys, 2 data            ACVI4124
      inodes = mxtree / 6 - 1                                           ACVI4125
      call tsetlf( lctree, inodes, 2, 1 )   ! 2 keys, 1 data            ACVI4126
C                                                                       ACVI4127
      do i = 1, mxwcof                                                  ACVI4128
        dwigsu3( i ) = 0.d0                                             ACVI4129
      end do                                                            ACVI4130
C                                                                       ACVI4131
C     Occupance numbers in bit forms.                                   ACVI4132
      do i = 0, 31                                                      ACVI4133
        if( btest(nlevel,i) ) nbits=i+1   ! check # bits in the level   ACVI4134
      end do                                                            ACVI4135
      nbitx = ishft( 1, nbits ) - 1                                     ACVI4136
C                                                                       ACVI4137
C     ... Start ...                                                     ACVI4138
C                                                                       ACVI4139
      call su3fp( ieta, nx, lmmua, iptr, MXNUMB )                       ACVI4140
C     Initialize a+ operator                                            ACVI4141
      call a1op( mltree, next, levels, nxy, dwigsu3 )                   ACVI4142
C                                                                       ACVI4143
      call optlm( mltree, lbtree, lctree, next, levels, nxy, dwigsu3,   ACVI4144
     1  mxwcof, 1, ieta, 0,-1, 0, ieta, 0)! a+ a                        ACVI4145
      call optlm( mltree, lbtree, lctree, next, levels, nxy, dwigsu3,   ACVI4146
     1  mxwcof, 1, ieta, 0, 1, ieta, 0, 0)! a+ a+                       ACVI4147
      do i = 2, nx - 1                    ! a+ ... a+                   ACVI4148
        ibot = iptr( i )                                                ACVI4149
        itop = iptr( i+1 ) - 1                                          ACVI4150
        do j = ibot, itop                                               ACVI4151
          lmmu = lmmua( 2, j )                                          ACVI4152
          lm   = iupac( lmmu, -16, nbit8 )                              ACVI4153
          mu   = iupac( lmmu,  -8, nbit8 )                              ACVI4154
          call optlm( mltree, lbtree, lctree, next, levels, nxy,        ACVI4155
     1      dwigsu3, mxwcof, 1, ieta, 0, i, lm, mu, j )                 ACVI4156
        end do                                                          ACVI4157
      end do                                                            ACVI4158
C                                                                       ACVI4159
      if( tactive ) then                                                ACVI4160
        if( n1 .eq. 0 ) then                                            ACVI4161
          ibot = 1                                                      ACVI4162
          itop = 1                                                      ACVI4163
        else                                                            ACVI4164
          ibot = iptr( n1 )                                             ACVI4165
          itop = iptr( n1+1 ) - 1                                       ACVI4166
        end if                                                          ACVI4167
        if( n2 .eq. 0 ) then                                            ACVI4168
          jbot = 1                                                      ACVI4169
          jtop = 1                                                      ACVI4170
        else                                                            ACVI4171
          jbot = iptr( n2 )                                             ACVI4172
          jtop = iptr( n2+1 ) - 1                                       ACVI4173
        end if                                                          ACVI4174
        do i = ibot, itop                                               ACVI4175
          lmmu = lmmua( 2, i )                                          ACVI4176
          if( n1 .eq. 0 ) lmmu = 0                                      ACVI4177
          lm1  = iupac( lmmu, -16, nbit8 )                              ACVI4178
          mu1  = iupac( lmmu,  -8, nbit8 )                              ACVI4179
          do j = jbot, jtop                                             ACVI4180
            lmmu = lmmua( 2, j )                                        ACVI4181
            if( n2 .eq. 0 ) lmmu = 0                                    ACVI4182
            lm2  = iupac( lmmu, -16, nbit8 )                            ACVI4183
            mu2  = iupac( lmmu,  -8, nbit8 )                            ACVI4184
            call optlm( mltree, lbtree, lctree, next, levels, nxy,      ACVI4185
     1        dwigsu3, mxwcof, n1, lm1, mu1,-n2, mu2, lm2, 0 )          ACVI4186
          end do                                                        ACVI4187
        end do                                                          ACVI4188
      end if                                                            ACVI4189
C                                                                       ACVI4190
      nx = maxy                                                         ACVI4191
      return                                                            ACVI4192
C ---*end of opb*------------------------------------------------------ ACVI4193
      end                                                               ACVI4194
C ----------------------------------------------------------------------ACVI4195
C                                                                       ACVI4196
C                           *************                               ACVI4197
C                           *** A1OP  ***                               ACVI4198
C                           *************                               ACVI4199
C                                                                       ACVI4200
C ----------------------------------------------------------------------ACVI4201
C Author:  Chairul Bahri (LSU)                                          ACVI4202
C References:  -*-                                                      ACVI4203
C ----------------------------------------------------------------------ACVI4204
C Updates: 11/90 ==> original (IBM 3090)                                ACVI4205
C          11/28/90: implemented in RMECODE.                            ACVI4206
C          04/09/90: phase included.                                    ACVI4207
C ----------------------------------------------------------------------ACVI4208
C                                                                       ACVI4209
C    The subroutine A1OP generates the matrix elements of oscillator    ACVI4210
C    creation operator a+ between various SU(3) basis. This is the      ACVI4211
C    building block for any SU(3) tensors.                              ACVI4212
C                                                                       ACVI4213
C ----------------------------------------------------------------------ACVI4214
C                                                                       ACVI4215
      subroutine a1op( mltree, next, levels, nxy, dwigsu3 )             ACVI4216
      implicit real*8(d), logical(t)                                    ACVI4217
      character*4 cflag, cphase                                         ACVI4218
      parameter ( D1 = 1.d0 )                                           ACVI4219
C                                                                       ACVI4220
      dimension mltree(-10:* ),    ! binary tree for uncoupled tensors  ACVI4221
     1          next( * ),         ! linked list index for mltree       ACVI4222
     2          levels( * ),       ! U(3) levels                        ACVI4223
     3          nxy( * ),          ! (nx,ny) values for mltree          ACVI4224
     4          dwigsu3( * )       ! = 1.0, cfp                         ACVI4225
      common / OPBCON / lfile,     ! logfile                            ACVI4226
     1                  ieta,      ! oscillator shell number            ACVI4227
     2                  nbits,     ! bit length for # levels in shell   ACVI4228
     3                  nbitx,     ! significant bits                   ACVI4229
     4                  cflag,     ! flag for new SU(3) irreps          ACVI4230
     5                  cphase,    ! SU(3) phase                        ACVI4231
     6                  tfile      ! yes, write to log file             ACVI4232
      common / INDICS / ix,        ! mltree index for dwigsu3           ACVI4233
     1                  iy,        ! lbtree index for dwigsu3           ACVI4234
     2                  maxy       ! max index for dwigsu3              ACVI4235
      dimension label( 2 )                                              ACVI4236
C     Packing functions                                                 ACVI4237
      data nbit8 / 255 /           ! zff                                ACVI4238
      ipack( i, j, k, l ) = ior(l,ishft(ior(k,ishft(ior(j,ishft(i,      ACVI4239
     1       8)),8)),8))                                                ACVI4240
      iupac( index, iover, ibits) = iand( ishft( index, iover), ibits ) ACVI4241
C                                                                       ACVI4242
C a+(ieta 0)                                                            ACVI4243
      ix = 0                                                            ACVI4244
CW----------------------------------------------------------------------ACVI4245
      write(6,100) ieta                                                 ACVI4246
      if( tfile ) write(lfile,100) ieta                                 ACVI4247
100   format(/' operator 1: a+(',i3,'  0)'/)                            ACVI4248
CW----------------------------------------------------------------------ACVI4249
      cflag(2:2) = '#'                                                  ACVI4250
      do nz = ieta, 0, -1                                               ACVI4251
        ie3 = 3*nz - ieta                                               ACVI4252
        j3t = ieta - nz                                                 ACVI4253
CW----------------------------------------------------------------------ACVI4254
      if( tfile ) write(lfile,'(a,2i3)') '   >> eps  2Lam:',ie3,j3t     ACVI4255
CW----------------------------------------------------------------------ACVI4256
        label( 1 ) = ipack(1, ieta, 0, j3t)                             ACVI4257
        label( 2 ) = ipack(1, 0, 0, nz)                                 ACVI4258
        do nx = j3t, 0, -1                                              ACVI4259
          ix  = ix + 1                                                  ACVI4260
          ny  = ieta - (nz + nx)                                        ACVI4261
          m3t = nx - ny                                                 ACVI4262
          dwigsu3( ix ) = D1                                            ACVI4263
          call tchk( label, mltree, *110 )                              ACVI4264
            cflag(3:3) = '*'                                            ACVI4265
            iloc = mltree(-5 )                                          ACVI4266
            mltree( iloc + mltree(-2) ) = 1                             ACVI4267
            call tins( label, mltree )                                  ACVI4268
            call tchk( label, mltree, *110 )                            ACVI4269
110       iloc = mltree(-5 ) + mltree(-4 ) + 1                          ACVI4270
          next( ix )       = mltree( iloc+1 )                           ACVI4271
          mltree( iloc )   = mltree( iloc ) + 1                         ACVI4272
          mltree( iloc+1 ) = ix                                         ACVI4273
          levels( ix )     = ix                                         ACVI4274
          nxy( ix )        = ipack( 0, 0, nx, ny )                      ACVI4275
CW----------------------------------------------------------------------ACVI4276
CW    Write coefficients.                                               ACVI4277
      if( tfile ) then                                                  ACVI4278
        if( ix .eq. 1 ) then                                            ACVI4279
          write(lfile,120) cflag, ix, j3t, m3t, ix, dwigsu3( ix )       ACVI4280
        else                                                            ACVI4281
          write(lfile,130) cflag, ix, j3t, m3t, ix, dwigsu3( ix )       ACVI4282
        end if                                                          ACVI4283
      end if                                                            ACVI4284
120   format(18x,'2Lam',10x,'2M',23x,'<  ;  |  >'/                      ACVI4285
     1       a,i6,2i12,12x,' a+(',i2,')',f14.7)                         ACVI4286
130   format(a,i6,2i12,12x,' a+(',i2,')',f14.7)                         ACVI4287
CW----------------------------------------------------------------------ACVI4288
          cflag(3:3) = ' '                                              ACVI4289
        end do                                                          ACVI4290
      end do                                                            ACVI4291
      maxy = ix                                                         ACVI4292
      return                                                            ACVI4293
C ---*end of a1op*------------------------------------------------------ACVI4294
      end                                                               ACVI4295
C ----------------------------------------------------------------------ACVI4296
C                                                                       ACVI4297
C                           *************                               ACVI4298
C                           *** OPTLM ***                               ACVI4299
C                           *************                               ACVI4300
C                                                                       ACVI4301
C ----------------------------------------------------------------------ACVI4302
C Author:  Chairul Bahri (LSU)                                          ACVI4303
C References:  -*-                                                      ACVI4304
C ----------------------------------------------------------------------ACVI4305
C Updates: 11/90 ==> original (IBM 3090)                                ACVI4306
C          11/28/90: implemented in RMECODE.                            ACVI4307
C          04/09/90: phase included.                                    ACVI4308
C          11/17/92: IBM RS/6000-560.                                   ACVI4309
C ----------------------------------------------------------------------ACVI4310
C                                                                       ACVI4311
C    The subroutine OPTLM generates the matrix elements of coupled SU(3)ACVI4312
C    tensor operators t1xt2 (lm mu) between various SU(3) basis, recur- ACVI4313
C    sively.                                                            ACVI4314
C                                                                       ACVI4315
C ----------------------------------------------------------------------ACVI4316
C                                                                       ACVI4317
      subroutine optlm( mltree, lbtree, lctree, next, levels, nxy,      ACVI4318
     ]           dwigsu3, mxwcof, n1, lm1, mu1, n2a, lm2, mu2, icfp )   ACVI4319
      implicit real*8(d), logical(t)                                    ACVI4320
      character*4 cflag, cphase                                         ACVI4321
      parameter ( DZERO = 1.d-12 )                                      ACVI4322
C                                                                       ACVI4323
      dimension mltree(-10:* ),    ! binary tree for uncoupled tensors  ACVI4324
     1          lbtree(-10:* ),    ! binary tree for coupled tensors    ACVI4325
     2          lctree(-10:* ),    ! binary tree for uncoupledtensors   ACVI4326
     3          next( * ),         ! linked list index for mltree       ACVI4327
     4          levels( * ),       ! U(3) levels                        ACVI4328
     5          nxy( * ),          ! (nx,ny) values for mltree          ACVI4329
     6          dwigsu3( * )       ! coefficients for coupled SU(3)     ACVI4330
C               mxwcof             ! max # operator basis               ACVI4331
C               n1                 ! # a+ in t1                         ACVI4332
C               lm1, mu1           ! SU(3) tensor of t1                 ACVI4333
C               n2a                ! # a+ (or a) in t2                  ACVI4334
C               lm2, mu2           ! SU(3) tensor of t2                 ACVI4335
C               icfp               ! fractional parentage index for     ACVI4336
C                                  !   more than a+a+                   ACVI4337
      common / OPBCON / lfile,     ! logfile                            ACVI4338
     1                  ieta,      ! oscillator shell number            ACVI4339
     2                  nbits,     ! bit length for # levels in shell   ACVI4340
     3                  nbitx,     ! significant bits                   ACVI4341
     4                  cflag,     ! flag for new SU(3) irreps          ACVI4342
     5                  cphase,    ! SU(3) phase                        ACVI4343
     6                  tfile      ! yes, write to log file             ACVI4344
      character*4 cann             ! annihilation operator              ACVI4345
      common / INDICS / ix,        ! mltree index for dwigsu3           ACVI4346
     1                  iy,        ! lbtree index for dwigsu3           ACVI4347
     2                  maxy       ! max index for dwigsu3              ACVI4348
      dimension label( 4 ),        ! labels for trees                   ACVI4349
     1          lvla( 8 )          ! packed occupance number, max 4     ACVI4350
C                                                                       ACVI4351
C     SU(3) package                                                     ACVI4352
C                                                                       ACVI4353
      parameter ( NCE1 = 9, NCE2 = 13244, NCEX = 42 )                   ACVI4354
      parameter ( NCW1 = 9, NCW2 = 42, NCW3 = 9030 )                    ACVI4355
      parameter ( KIMAX1 = 3*NCE2, KIMAX2 = 3*NCW3 )                    ACVI4356
      parameter ( NCW22 = NCW2*NCW2 )                                   ACVI4357
C     common / SU3EXT / tj2ta( NCE2 )                                   ACVI4358
      common / SU3EXI / j1ta( NCE2 ), j2ta( NCE2 ), iea( NCE2 )         ACVI4359
      common / SU3I   / j1smax( NCW22 ), j1tmax( NCW22 ),               ACVI4360
     1                  j2smax( NCW2 ), j2tmax( NCW2 ),  indmat( NCW22 )ACVI4361
      common / SU3EXD / dewu3( KIMAX1 )                                 ACVI4362
      common / SU3D   / dwu3( KIMAX2 )                                  ACVI4363
C                                                                       ACVI4364
C     Packing functions                                                 ACVI4365
      data nbit4, nbit8 / 16, 255 /! zf, zff                            ACVI4366
      ipack( i, j, k, l ) = ior( l, ishft( ior( k, ishft( ior( j,       ACVI4367
     1        ishft( i, 8 )), 8 )), 8 ))                                ACVI4368
      ipack5( i, j, k, l, m ) = ior( m, ishft( ior( l, ishft( ior( k,   ACVI4369
     1        ishft( ior( j, ishft( i, 4 )), 8 )), 8 )), 8 ))           ACVI4370
      iupac(index,iover,ibits) = iand( ishft( index, iover), ibits )    ACVI4371
C                                                                       ACVI4372
C t(lm1 mu1) x t(lm2 mu2)                                               ACVI4373
C n1           n2                                                       ACVI4374
C              n2a<0 for annihilation operators.                        ACVI4375
      iz = 1                                                            ACVI4376
      n2 = iabs( n2a )                                                  ACVI4377
      if( n1 .lt. 0 ) call error(' OPTLM: Creation ops. preferable.')   ACVI4378
      ntot   = n1 + n2                                                  ACVI4379
      ndbits = iabs( n1 - n2 ) * nbits                                  ACVI4380
      nmbits = min0( n1, n2 )  * nbits                                  ACVI4381
      n1eta  = n1 * ieta                                                ACVI4382
      n2eta  = n2 * ieta                                                ACVI4383
      if( n2a .lt. 0 ) iphsu3 = 2 * (mu2 - lm2)   ! SU(3) phase: phi    ACVI4384
      if( n2a .lt. 0 ) then                                             ACVI4385
        cann = ') a('                                                   ACVI4386
      else                                                              ACVI4387
        cann = ')a+('                                                   ACVI4388
      end if                                                            ACVI4389
CW----------------------------------------------------------------------ACVI4390
      write(6,100) ntot, n1, lm1, mu1, n2, lm2, mu2                     ACVI4391
      if( tfile ) write(lfile,100) ntot, n1, lm1, mu1, n2, lm2, mu2     ACVI4392
100   format(//' operator ',i1,': ',i1,'(',2i3,') x ',i1,'(',2i3,')'/)  ACVI4393
CW----------------------------------------------------------------------ACVI4394
      label1 = ipack( lm1, mu1, lm2, mu2 )                              ACVI4395
      label( 4 ) = 0                                                    ACVI4396
      lam = lm1 + mu1 + lm2 + mu2                                       ACVI4397
C     Run over all possible SU(3) coupled tensors                       ACVI4398
      do  220 lm = 0, lam                                               ACVI4399
       do 220 mu = 0, lam                                               ACVI4400
        call u3mult( lm1, mu1, lm2, mu2, lm, mu, mult, *220 )           ACVI4401
CW----------------------------------------------------------------------ACVI4402
      write(6,110) lm, mu                                               ACVI4403
      if( tfile ) write(lfile,110) lm, mu                               ACVI4404
110   format(' SU(3) *** (',2i3,') ***')                                ACVI4405
CW----------------------------------------------------------------------ACVI4406
C                                                                       ACVI4407
C                           SU(3) > SU(2) x U(1)                        ACVI4408
C                   <   t1   ;   t2   |   T  > v for all isoscalars     ACVI4409
          call xewu3( lm1,mu1, lm2,mu2, lm,mu, 1, NEC,                  ACVI4410
     1      kromax,indmax, dewu3,j1ta,j2ta,iea, NCE1,NCE2, KIMAX1 )     ACVI4411
          ie3top = 2 * lm + mu                                          ACVI4412
          do ip = 0, lm                                                 ACVI4413
            do iq = 0, mu                                               ACVI4414
              ie3 = ie3top - 3*(ip+iq)                                  ACVI4415
              j3t = mu + ip - iq                                        ACVI4416
              label2 = ipack5( icfp, mult, lm, mu, j3t )                ACVI4417
CW----------------------------------------------------------------------ACVI4418
      if( tfile ) write(lfile,'(a,2i3)') '   >> eps  2Lam:',ie3,j3t     ACVI4419
CW----------------------------------------------------------------------ACVI4420
              call xwu3( lm1,mu1, lm2,mu2, lm,mu, ie3,j3t, NEC, dewu3,  ACVI4421
     1          kromax,indmax, dwu3,j1smax,j1tmax,j2smax,j2tmax,        ACVI4422
     2          iesmax,ie2max, indmat, NCW1,NCW2,NCW3, KIMAX2 )         ACVI4423
              do ies = 1, iesmax                                        ACVI4424
                ie2 = ie2max - 3*(iesmax-ies)                           ACVI4425
                ie1 = ie3 - ie2                                         ACVI4426
                nz1 = (n1eta + ie1) / 3      ! ie1=3*nz1-n1*ieta        ACVI4427
                if( n2a .lt. 0 ) then                                   ACVI4428
                  iphu1 = iphsu3 + ie2                                  ACVI4429
                  ie2   = -ie2                                          ACVI4430
                end if                                                  ACVI4431
                nz2 = (n2eta + ie2) / 3                                 ACVI4432
C                                                                       ACVI4433
C                           SU(2) > U(1)                                ACVI4434
C                                                                       ACVI4435
C               Run over all J2T of t2                                  ACVI4436
                do 210 j2s = 1, j2smax( ies )                           ACVI4437
                  j2t = j2tmax(ies) - 2*(j2s-1)                         ACVI4438
                  iesj2s = ies + NCW2*(j2s-1)                           ACVI4439
                  if( n2a.ge.0 ) then                                   ACVI4440
                    label( 1 ) = ipack( 1, lm2, mu2, j2t )              ACVI4441
                  else                                                  ACVI4442
                    label( 1 ) = ipack( 1, mu2, lm2, j2t )              ACVI4443
                  end if                                                ACVI4444
                  label( 2 ) = ipack( n2, 0, 0, nz2 )                   ACVI4445
C                 Load t2 from mltree                                   ACVI4446
                  call tchk( label, mltree, *120 )                      ACVI4447
CW                  write(lfile,*) 't2 not found'                       ACVI4448
                    go to 210                                           ACVI4449
120               iloc   = mltree(-5 ) + mltree(-4 ) + 1                ACVI4450
                  ivecmx = mltree( iloc )                               ACVI4451
                  inext  = mltree( iloc + 1 )                           ACVI4452
                  do iv = ivecmx, 1, -1                                 ACVI4453
                    lvl2  = levels( inext )                             ACVI4454
                    index = nxy( inext )                                ACVI4455
                    nx    = iupac( index, -8, nbit8 )                   ACVI4456
                    ny    = iupac( index,  0, nbit8 )                   ACVI4457
                    dwig2 = dwigsu3( inext )                            ACVI4458
                    m2t   = nx - ny          ! m = (nx-ny)/2            ACVI4459
                    if( n2a .lt. 0 ) then                               ACVI4460
                      iphase = (3*m2t + iphu1) / 6                      ACVI4461
                      if( btest(iphase,0) ) then                        ACVI4462
                        dwig2  = -dwig2                                 ACVI4463
                        cphase = ' (-)'                                 ACVI4464
                      end if                                            ACVI4465
                      m2t = -m2t                                        ACVI4466
                    end if                                              ACVI4467
CW----------------------------------------------------------------------ACVI4468
      lvls = lvl2                                                       ACVI4469
      do ilvl = 5, 8                                                    ACVI4470
        lvla( ilvl ) = iand( lvls, nbitx )                              ACVI4471
        lvls = ishft( lvls,-nbits )                                     ACVI4472
      end do                                                            ACVI4473
CW----------------------------------------------------------------------ACVI4474
C                   Run over all J1T of t1                              ACVI4475
                    do 200 j1s = 1, j1smax( iesj2s )                    ACVI4476
                      j1t = j1tmax(iesj2s) - 2*(j1s-1)                  ACVI4477
                      ind = (indmat(iesj2s) - j1t) / 2                  ACVI4478
                      label( 1 ) = ipack( 1, lm1, mu1, j1t )            ACVI4479
                      label( 2 ) = ipack( n1, 0, 0, nz1 )               ACVI4480
C                     Load t1 from mltree                               ACVI4481
                      call tchk( label, mltree, *130 )                  ACVI4482
CW                      write(lfile,*) 't1 not found'                   ACVI4483
                        go to 200                                       ACVI4484
130                   iloc   = mltree(-5 ) + mltree(-4 ) + 1            ACVI4485
                      jvecmx = mltree( iloc )                           ACVI4486
                      jnext  = mltree( iloc + 1 )                       ACVI4487
                      do jv = jvecmx, 1, -1                             ACVI4488
                        lvl1  = levels( jnext )                         ACVI4489
                        index = nxy( jnext )                            ACVI4490
                        nx    = iupac( index, -8, nbit8 )               ACVI4491
                        ny    = iupac( index,  0, nbit8 )               ACVI4492
                        dwig1 = dwig2 * dwigsu3( jnext )                ACVI4493
                        m1t   = nx - ny                                 ACVI4494
CW----------------------------------------------------------------------ACVI4495
      lvls = lvl1                                                       ACVI4496
      do ilvl = 1, 4                                                    ACVI4497
        lvla( ilvl ) = iand( lvls, nbitx )                              ACVI4498
        lvls = ishft( lvls,-nbits )                                     ACVI4499
      end do                                                            ACVI4500
CW----------------------------------------------------------------------ACVI4501
C                                                                       ACVI4502
C                       SU(3) coupling                                  ACVI4503
C                                                                       ACVI4504
                        m3t = m1t + m2t                                 ACVI4505
                        if( iabs(m3t)   .le. j3t .and.                  ACVI4506
     1                      dabs(dwig1) .ge. DZERO ) then               ACVI4507
C                         Set key for trees                             ACVI4508
                          label( 1 ) = label1                           ACVI4509
                          label( 2 ) = label2                           ACVI4510
                          jflag = 0                                     ACVI4511
                          if( n2a .ge. 0 ) then                         ACVI4512
                            lvls  = lvl1                                ACVI4513
                            index = lvl2                                ACVI4514
                            lvls  = ishft( lvls, ndbits )               ACVI4515
                            if( lvls .le. index ) then                  ACVI4516
                              jflag = 1                                 ACVI4517
                              lvls  = ior( index, ishft(lvls,nmbits) )  ACVI4518
                              label( 3 ) = lvls                         ACVI4519
                              label( 4 ) = 0                            ACVI4520
                            end if                                      ACVI4521
                          else                                          ACVI4522
                            jflag = 1                                   ACVI4523
                            label( 3 ) = lvl1                           ACVI4524
                            label( 4 ) = lvl2                           ACVI4525
                          end if                                        ACVI4526
C                         Load the trees                                ACVI4527
                          if( jflag .eq. 1 ) then                       ACVI4528
                            ix = ix + 1                                 ACVI4529
                            dwig1 = dwig1*dwr3(j1t,j2t,j3t,m1t,m2t,m3t) ACVI4530
C   lbtree                                                              ACVI4531
                            call tchk( label, lbtree, *140 )            ACVI4532
                              if( lbtree(-10).gt.lbtree(-9) ) go to 999 ACVI4533
                              cflag(2:2) = '#'                          ACVI4534
                              iloc = lbtree(-5 )                        ACVI4535
                              if( n2a .ge. 0 ) then                     ACVI4536
                                kprty = ipack(127,0,0,0) - 1            ACVI4537
                              else                                      ACVI4538
                                kprty = 0                               ACVI4539
                              end if                                    ACVI4540
                              lbtree( iloc+lbtree(-2) )   = kprty       ACVI4541
                              lbtree( iloc+lbtree(-4)+1 ) = maxy + 1    ACVI4542
                              lbtree( iloc+lbtree(-4)+2 ) = kromax      ACVI4543
                              call tins( label, lbtree )                ACVI4544
                              call tchk( label, lbtree, *140 )          ACVI4545
140                         iloc  = lbtree(-5 ) + lbtree(-4 ) + 1       ACVI4546
                            iy    = lbtree( iloc ) - 1                  ACVI4547
                            minro = kromax * (ind-1)                    ACVI4548
                            maxro = minro + kromax                      ACVI4549
                            miny  = iy                                  ACVI4550
                            if( iy+mult .gt. mxwcof )                   ACVI4551
     1                        call attn(' OPTLM: Too many opb', *999 )  ACVI4552
                            do indro = minro+1, maxro                   ACVI4553
                              iy = iy + 1                               ACVI4554
                              dwigsu3( iy ) = dwigsu3( iy ) +           ACVI4555
     1                                        dwig1 * dwu3( indro )     ACVI4556
                            end do ! indro                              ACVI4557
                            maxy = max0( maxy, iy )                     ACVI4558
C   lctree                                                              ACVI4559
                          if( n2a .ge. 0 ) then                         ACVI4560
                            label( 1 ) = label2                         ACVI4561
                            label( 2 ) = lvls                           ACVI4562
                            call tchk( label, lctree, *160 )            ACVI4563
                              if( lctree(-10).gt.lctree(-9) ) go to 999 ACVI4564
                              iloc = lctree(-5 )                        ACVI4565
                              lctree( iloc+lctree(-2) ) = kromax        ACVI4566
                              call tins( label, lctree )                ACVI4567
                              call tchk( label, lctree, *150 )          ACVI4568
150                           iloc = lctree(-5 ) + lctree(-4 ) + 1      ACVI4569
                              lctree( iloc ) = iy                       ACVI4570
160                         label( 2 ) = ipack( ntot, 0, 0, nz1+nz2 )   ACVI4571
C   mltree                                                              ACVI4572
                            call tchk( label, mltree, *170 )            ACVI4573
                              if( mltree(-10).gt.mltree(-9) ) go to 999 ACVI4574
                              cflag(3:3) = '*'                          ACVI4575
                              iloc = mltree(-5 )                        ACVI4576
                              mltree( iloc+mltree(-2) ) = kromax        ACVI4577
                              call tins( label, mltree )                ACVI4578
                              call tchk( label, mltree, *170 )          ACVI4579
170                         iloc = mltree(-5 ) + mltree(-4 ) + 1        ACVI4580
                            next( ix ) = mltree( iloc + 1 )             ACVI4581
                            mltree( iloc ) = mltree( iloc ) + 1         ACVI4582
                            mltree( iloc + 1 ) = ix                     ACVI4583
                            levels( ix ) = lvls                         ACVI4584
                            nxy( ix )    = nxy(inext) + nxy(jnext)      ACVI4585
                          end if                                        ACVI4586
CW----------------------------------------------------------------------ACVI4587
CW    Write coefficients.                                               ACVI4588
      if( tfile ) then                                                  ACVI4589
        if( iz .eq. 1 ) then                                            ACVI4590
          write(lfile,180) cflag, miny+1, j1t,j2t,j3t, m1t,m2t,m3t,     ACVI4591
     1      (lvla(i1),i1=1,4), cann, (lvla(i2),i2=5,8), ind, cphase,    ACVI4592
     2      jnext,inext,('  ',dwu3(minro+i),dwigsu3(miny+i),i=1,mult)   ACVI4593
          iz = 0                                                        ACVI4594
        else                                                            ACVI4595
          write(lfile,190) cflag, miny+1, j1t,j2t,j3t, m1t,m2t,m3t,     ACVI4596
     1      lvla, ind, cphase, jnext, inext,                            ACVI4597
     2      ('  ',dwu3(minro+i),dwigsu3(miny+i),i=1,mult)               ACVI4598
        end if                                                          ACVI4599
      end if                                                            ACVI4600
180   format(13x,'2L1 L2 L3',3x,'2M1 M2 M3',26x,'ind phi',13x,          ACVI4601
     1                                '<  ;  ||  >',4x,'<  ;  |  >'/    ACVI4602
     2       a,i6,3x,3i3,3x,3i3,' a+(',4i2,a,4i2,')',i4,a,2i4,a,2f14.7, ACVI4603
     3       10(a/77x,2f14.7))                                          ACVI4604
190   format(a,i6,3x,3i3,3x,3i3,4x,    4i2,4x,4i2,1x,i4,a,2i4,a,2f14.7, ACVI4605
     3       10(a/77x,2f14.7))                                          ACVI4606
CW----------------------------------------------------------------------ACVI4607
                          cflag = '    '                                ACVI4608
                        end if                                          ACVI4609
                        end if                                          ACVI4610
                        jnext = next( jnext )                           ACVI4611
                      end do ! jv                                       ACVI4612
200                 continue ! j1s                                      ACVI4613
                    cphase = '    '                                     ACVI4614
                    inext = next( inext )                               ACVI4615
                  end do ! iv                                           ACVI4616
210             continue ! j2s                                          ACVI4617
              end do ! ies                                              ACVI4618
            end do ! iq                                                 ACVI4619
          end do ! ip                                                   ACVI4620
220   continue ! mu, lm                                                 ACVI4621
      return                                                            ACVI4622
999   write(6,*) ' OPTLM: Tree overloaded.'                             ACVI4623
      stop                                                              ACVI4624
C ---*end of optlm*-----------------------------------------------------ACVI4625
      end                                                               ACVI4626
C ----------------------------------------------------------------------ACVI4627
C                                                                       ACVI4628
      subroutine su3fp( ieta, nx, lmmua, iptr, mxnumb )                 ACVI4629
C                                                                       ACVI4630
C ----------------------------------------------------------------------ACVI4631
C Author: Chairul Bahri (LSU 11/90 ... original)                        ACVI4632
C ----------------------------------------------------------------------ACVI4633
C     Sub-program to generate all SU(3) quantum numbers of coupled a+'s.ACVI4634
C ----------------------------------------------------------------------ACVI4635
C                                                                       ACVI4636
C               ieta              ! shell number                        ACVI4637
C               nx                ! generation level                    ACVI4638
      dimension lmmua( 2, * ),    ! array of packed (lm,mu) of SU(3)    ACVI4639
     1          iptr( * )         ! level pointers for generations      ACVI4640
C               mxnumb            ! max # (lm mu) pairs lmmua           ACVI4641
C     Packing functions                                                 ACVI4642
      data nbit8 / 255 /           ! zff                                ACVI4643
      ipack( i, j, k, l ) = ior(l,ishft(ior(k,ishft(ior(j,ishft(i,      ACVI4644
     1       8)),8)),8))                                                ACVI4645
      iupac( index, iover, ibits) = iand( ishft( index, iover), ibits ) ACVI4646
C                                                                       ACVI4647
C     Initialize lmmua                                                  ACVI4648
      do i = 1, mxnumb-2                                                ACVI4649
        lmmua( 1, i ) = 0                                               ACVI4650
        lmmua( 2, i ) = 0                                               ACVI4651
      end do                                                            ACVI4652
C                                                                       ACVI4653
C     First generation                                                  ACVI4654
      lmi = ieta                                                        ACVI4655
      mui = 0                                                           ACVI4656
      lmmua( 1, 1 ) = ipack( lmi, mui, 0, 0 )                           ACVI4657
      lmmua( 2, 1 ) = ipack( 1, lmi, mui, 0 )                           ACVI4658
      iptr( 1 ) = 1                                                     ACVI4659
      l = 1                                                             ACVI4660
      j = l                                                             ACVI4661
CW    write(6,110) 1                                                    ACVI4662
CW    write(6,120) 1, 0,0, ieta,0, 1,ieta,0                             ACVI4663
C     Next generations                                                  ACVI4664
      do i = 2, nx                                                      ACVI4665
CW      write(6,110) i                                                  ACVI4666
        do k = iptr( i-1 ), j                                           ACVI4667
          lmmuk = lmmua( 2, k )                                         ACVI4668
          lmk   = iupac( lmmuk, -16, nbit8 )                            ACVI4669
          muk   = iupac( lmmuk,  -8, nbit8 )                            ACVI4670
          lam   = lmi + mui + lmk + muk                                 ACVI4671
          do lm = 0, lam                                                ACVI4672
            do mu = 0, lam                                              ACVI4673
              call u3mult( lmi, mui, lmk, muk, lm, mu, kromax, *100 )   ACVI4674
                l = l + 1                                               ACVI4675
                if( l .gt. mxnumb-2 )                                   ACVI4676
     1            call attn(' SU3FP: Too many ops.',*999 )              ACVI4677
                lmmua( 1, l ) = ipack( lmi, mui, lmk, muk )             ACVI4678
                lmmua( 2, l ) = ipack( kromax, lm, mu, 0 )              ACVI4679
CW              write(6,120) l, lmi,mui, lmk,muk, kromax,lm,mu          ACVI4680
100           continue                                                  ACVI4681
            end do                                                      ACVI4682
          end do                                                        ACVI4683
        end do                                                          ACVI4684
        iptr( i ) = j + 1                                               ACVI4685
        j = l                                                           ACVI4686
      end do                                                            ACVI4687
      iptr( nx + 1 ) = l + 1                                            ACVI4688
CW    write(6,'(a,i3,a,i4)') ' last pointer(',nx+1,') ',iptr(nx+1)      ACVI4689
110   format(' Generation ',i3)                                         ACVI4690
120   format(8x,' # ',i3,' (',2i3,')x(',2i3,') .. ',i2,'(',2i3,')')     ACVI4691
999   return                                                            ACVI4692
C ---*end of su3fp*-----------------------------------------------------ACVI4693
      end                                                               ACVI4694
CB----------------------------------------------------------------------ACVI4695
C                          ***************                              ACVI4696
C                          **  UTILITY ***                              ACVI4697
C                          ***************                              ACVI4698
C ----------------------------------------------------------------------ACVI4699
      SUBROUTINE SETBIN(LFILX,TWRITX,BINX,NEMPTX,NFIRSTX,NMAXX,NTOTALX) ACVI4700
C     Generating the binary I/O.                                        ACVI4701
C                                                                       ACVI4702
      COMMON /SAVARG/ LFILE, TWRITE, NEMPTY, NFIRST, NMAX, NTOTAL       ACVI4703
      COMMON /SAVCHR/ BIN                                               ACVI4704
      COMMON /SAVINT/ IPM, NMAXP1, NMAX2, NMAXP2, NMAX3                 ACVI4705
      LOGICAL TWRITE, TWRITX                                            ACVI4706
      CHARACTER*(*) BINX                                                ACVI4707
      CHARACTER BIN*115                                                 ACVI4708
      LFILE = LFILX                                                     ACVI4709
      TWRITE = TWRITX                                                   ACVI4710
      NEMPTY = NEMPTX                                                   ACVI4711
      NFIRST = NFIRSTX                                                  ACVI4712
      NMAX = NMAXX                                                      ACVI4713
      NTOTAL = NTOTALX                                                  ACVI4714
      NMAXP1=NMAX+6                                                     ACVI4715
      NMAX2=NMAX+NMAX+5                                                 ACVI4716
      NMAXP2=NMAX2+6                                                    ACVI4717
      NMAX3=3*NMAX+10                                                   ACVI4718
      IF(NFIRST.LE.0)THEN                                               ACVI4719
        IPM=1                                                           ACVI4720
      ELSE                                                              ACVI4721
        IPM=-1                                                          ACVI4722
      ENDIF                                                             ACVI4723
      DO J=1,115                                                        ACVI4724
        BIN(J:J)=' '                                                    ACVI4725
      ENDDO                                                             ACVI4726
      RETURN                                                            ACVI4727
C                                                                       ACVI4728
      ENTRY DBOVRL(ILEFT,IL,IR,IRGHT)                                   ACVI4729
C     Converting a number from decimal to binary. (1a)                  ACVI4730
C-----------------------------------------------------------------------ACVI4731
C     Reading bits in number.                                           ACVI4732
      I=IPM                                                             ACVI4733
      DO J=1,NMAX                                ! LEFT                 ACVI4734
        IF(MOD(J,NEMPTY).NE.0)THEN                                      ACVI4735
           IF(BTEST(ILEFT,NFIRST+I))THEN                                ACVI4736
              BIN(J:J)='1'                                              ACVI4737
           ELSE                                                         ACVI4738
              BIN(J:J)='0'                                              ACVI4739
           ENDIF                                                        ACVI4740
           I=I+IPM                                                      ACVI4741
        ENDIF                                                           ACVI4742
      ENDDO                                                             ACVI4743
                  IF(NMAX.EQ.NTOTAL) THEN                               ACVI4744
C     Writing number in binary form.                                    ACVI4745
      IF(TWRITE) THEN                                                   ACVI4746
C       WRITE(LFILE,100) ILEFT,BIN                                      ACVI4747
        write(lfile,100) ileft, (bin(i:i), i=1, j)                      ACVI4748
      ELSE                                                              ACVI4749
        DO J = 1, NTOTAL                                                ACVI4750
          BINX(J:J) = BIN(J:J)                                          ACVI4751
        END DO                                                          ACVI4752
      END IF                                                            ACVI4753
                  ELSE                                                  ACVI4754
      I=IPM                                                             ACVI4755
      NMIN=MOD(NMAXP1-1,NEMPTY)                                         ACVI4756
      DO J=NMAXP1,NMAX2                          ! MIDDLE               ACVI4757
        IF(MOD(J,NEMPTY).NE.NMIN)THEN                                   ACVI4758
           NF=NFIRST+I                                                  ACVI4759
           IF(BTEST(IL,NF))THEN                                         ACVI4760
              IF(BTEST(IR,NF))THEN                                      ACVI4761
                 BIN(J:J)='I'                                           ACVI4762
              ELSE                                                      ACVI4763
                 BIN(J:J)='V'                                           ACVI4764
              ENDIF                                                     ACVI4765
           ELSE                                                         ACVI4766
              IF(BTEST(IR,NF))THEN                                      ACVI4767
                 BIN(J:J)='A'                                           ACVI4768
              ELSE                                                      ACVI4769
                 BIN(J:J)='-'                                           ACVI4770
              ENDIF                                                     ACVI4771
           ENDIF                                                        ACVI4772
           I=I+IPM                                                      ACVI4773
        ENDIF                                                           ACVI4774
      ENDDO                                                             ACVI4775
      I=IPM                                                             ACVI4776
      NMIN=MOD(NMAXP2-1,NEMPTY)                                         ACVI4777
      DO J=NMAXP2,NMAX3                          ! RIGHT                ACVI4778
        IF(MOD(J,NEMPTY).NE.NMIN)THEN                                   ACVI4779
           IF(BTEST(IRGHT,NFIRST+I))THEN                                ACVI4780
              BIN(J:J)='1'                                              ACVI4781
           ELSE                                                         ACVI4782
              BIN(J:J)='0'                                              ACVI4783
           ENDIF                                                        ACVI4784
           I=I+IPM                                                      ACVI4785
        ENDIF                                                           ACVI4786
      ENDDO                                                             ACVI4787
C     Writing number in binary form.                                    ACVI4788
C     IF(TWRITE) WRITE(LFILE,'(T2,A)') BIN                              ACVI4789
      if( twrite ) write(lfile,'(t2,115a)') (bin(i:i), i=1, j)          ACVI4790
                  ENDIF                                                 ACVI4791
      RETURN                                                            ACVI4792
C                                                                       ACVI4793
      ENTRY DBOVRP(NBITS,ILEFT,IL,IR,IRGHT)                             ACVI4794
C     Converting a number from decimal to binary. (1b)                  ACVI4795
C-----------------------------------------------------------------------ACVI4796
C     Reading bits in number.                                           ACVI4797
      I=IPM                                                             ACVI4798
      J=0                                                               ACVI4799
      DO WHILE(J.LE.NMAX .AND. IABS(I).LE.NBITS)   ! LEFT               ACVI4800
        J=J+1                                                           ACVI4801
        IF(MOD(J,NEMPTY).NE.0) THEN                                     ACVI4802
           IF(BTEST(ILEFT,NFIRST+I))THEN                                ACVI4803
              BIN(J:J)='1'                                              ACVI4804
           ELSE                                                         ACVI4805
              BIN(J:J)='0'                                              ACVI4806
           ENDIF                                                        ACVI4807
           I=I+IPM                                                      ACVI4808
        ENDIF                                                           ACVI4809
      ENDDO                                                             ACVI4810
      I=IPM                                                             ACVI4811
      J=NMAXP1-1                                                        ACVI4812
      NMIN=MOD(J,NEMPTY)                                                ACVI4813
      DO WHILE(J.LE.NMAX2 .AND. IABS(I).LE.NBITS)   ! MIDDLE            ACVI4814
        J=J+1                                                           ACVI4815
        IF(MOD(J,NEMPTY).NE.NMIN) THEN                                  ACVI4816
           NF=NFIRST+I                                                  ACVI4817
           IF(BTEST(IL,NF))THEN                                         ACVI4818
              IF(BTEST(IR,NF))THEN                                      ACVI4819
                 BIN(J:J)='I'                                           ACVI4820
              ELSE                                                      ACVI4821
                 BIN(J:J)='A'                                           ACVI4822
              ENDIF                                                     ACVI4823
           ELSE                                                         ACVI4824
              IF(BTEST(IR,NF))THEN                                      ACVI4825
                 BIN(J:J)='V'                                           ACVI4826
              ELSE                                                      ACVI4827
                 BIN(J:J)='-'                                           ACVI4828
              ENDIF                                                     ACVI4829
           ENDIF                                                        ACVI4830
           I=I+IPM                                                      ACVI4831
        ENDIF                                                           ACVI4832
      ENDDO                                                             ACVI4833
      I=IPM                                                             ACVI4834
      J=NMAXP2-1                                                        ACVI4835
      NMIN=MOD(J,NEMPTY)                                                ACVI4836
      DO WHILE(J.LE.NMAX3 .AND. IABS(I).LE.NBITS)   ! RIGHT             ACVI4837
        J=J+1                                                           ACVI4838
        IF(MOD(J,NEMPTY).NE.NMIN) THEN                                  ACVI4839
           IF(BTEST(IRGHT,NFIRST+I))THEN                                ACVI4840
              BIN(J:J)='1'                                              ACVI4841
           ELSE                                                         ACVI4842
              BIN(J:J)='0'                                              ACVI4843
           ENDIF                                                        ACVI4844
           I=I+IPM                                                      ACVI4845
        ENDIF                                                           ACVI4846
      ENDDO                                                             ACVI4847
C     Writing number in binary form.                                    ACVI4848
      IF(TWRITE) THEN                                                   ACVI4849
C       WRITE(LFILE,'(T2,A)') BIN                                       ACVI4850
        write(lfile,'(t2,115a)') (bin(i:i), i=1, j)                     ACVI4851
      ELSE                                                              ACVI4852
        DO J = 1, NTOTAL                                                ACVI4853
          BINX(J:J) = BIN(J:J)                                          ACVI4854
        END DO                                                          ACVI4855
      END IF                                                            ACVI4856
      RETURN                                                            ACVI4857
C                                                                       ACVI4858
      ENTRY DBOVRC(NBITS,ILEFT,IL,IR,IRGHT)                             ACVI4859
C     Converting a number from decimal to binary. (2)                   ACVI4860
C-----------------------------------------------------------------------ACVI4861
C     Reading bits in number.                                           ACVI4862
      I=IPM                                                             ACVI4863
      J=0                                                               ACVI4864
      DO WHILE(J.LE.NMAX .AND. IABS(I).LE.NBITS)   ! LEFT               ACVI4865
        J=J+1                                                           ACVI4866
        IF(MOD(J,NEMPTY).NE.0) THEN                                     ACVI4867
           IF(BTEST(ILEFT,NFIRST+I))THEN                                ACVI4868
              BIN(J:J)='1'                                              ACVI4869
           ELSE                                                         ACVI4870
              BIN(J:J)='0'                                              ACVI4871
           ENDIF                                                        ACVI4872
           I=I+IPM                                                      ACVI4873
        ENDIF                                                           ACVI4874
      ENDDO                                                             ACVI4875
                  IF(NMAX.EQ.NTOTAL) THEN                               ACVI4876
C     Writing number in binary form.                                    ACVI4877
      IF(TWRITE) THEN                                                   ACVI4878
C       WRITE(LFILE,100) ILEFT,BIN                                      ACVI4879
        write(lfile,100) ileft, (bin(i:i), i=1, j)                      ACVI4880
      ELSE                                                              ACVI4881
        DO J = 1, NTOTAL                                                ACVI4882
          BINX(J:J) = BIN(J:J)                                          ACVI4883
        END DO                                                          ACVI4884
      END IF                                                            ACVI4885
                  ELSE                                                  ACVI4886
      NMAXP1=NBITS+6                                                    ACVI4887
      NMAX2=NBITS+NBITS+5                                               ACVI4888
      NMAXP2=NMAX2+6                                                    ACVI4889
      NMAX3=3*NBITS+10                                                  ACVI4890
      I=IPM                                                             ACVI4891
      J=NMAXP1-1                                                        ACVI4892
      NMIN=MOD(J,NEMPTY)                                                ACVI4893
      DO WHILE(J.LE.NMAX2 .AND. IABS(I).LE.NBITS)   ! MIDDLE            ACVI4894
        J=J+1                                                           ACVI4895
        IF(MOD(J,NEMPTY).NE.NMIN) THEN                                  ACVI4896
           NF=NFIRST+I                                                  ACVI4897
           IF(BTEST(IL,NF))THEN                                         ACVI4898
              IF(BTEST(IR,NF))THEN                                      ACVI4899
                 BIN(J:J)='I'                                           ACVI4900
              ELSE                                                      ACVI4901
                 BIN(J:J)='A'                                           ACVI4902
              ENDIF                                                     ACVI4903
           ELSE                                                         ACVI4904
              IF(BTEST(IR,NF))THEN                                      ACVI4905
                 BIN(J:J)='V'                                           ACVI4906
              ELSE                                                      ACVI4907
                 BIN(J:J)='-'                                           ACVI4908
              ENDIF                                                     ACVI4909
           ENDIF                                                        ACVI4910
           I=I+IPM                                                      ACVI4911
        ENDIF                                                           ACVI4912
      ENDDO                                                             ACVI4913
      I=IPM                                                             ACVI4914
      J=NMAXP2-1                                                        ACVI4915
      NMIN=MOD(J,NEMPTY)                                                ACVI4916
      DO WHILE(J.LE.NMAX3 .AND. IABS(I).LE.NBITS)   ! RIGHT             ACVI4917
        J=J+1                                                           ACVI4918
        IF(MOD(J,NEMPTY).NE.NMIN) THEN                                  ACVI4919
           IF(BTEST(IRGHT,NFIRST+I))THEN                                ACVI4920
              BIN(J:J)='1'                                              ACVI4921
           ELSE                                                         ACVI4922
              BIN(J:J)='0'                                              ACVI4923
           ENDIF                                                        ACVI4924
           I=I+IPM                                                      ACVI4925
        ENDIF                                                           ACVI4926
      ENDDO                                                             ACVI4927
C     Writing number in binary form.                                    ACVI4928
C     IF(TWRITE) WRITE(LFILE,'(T2,A)') BIN                              ACVI4929
      if( twrite ) write(lfile,'(t2,115a)') (bin(i:i), i=1, j)          ACVI4930
                  ENDIF                                                 ACVI4931
C100  FORMAT(' ',I23,2X,A)                                              ACVI4932
100   format(' ',i23,2x,115a)                                           ACVI4933
      RETURN                                                            ACVI4934
C ---*end of SETBIN*----------------------------------------------------ACVI4935
      END                                                               ACVI4936
C                                                                       ACVI4937
      SUBROUTINE ATTN(LITER,*)                                          ACVI4938
C     Warning because of gross error.                                   ACVI4939
      CHARACTER*(*) LITER                                               ACVI4940
      WRITE(6,'(//A,A)') ' ***** ATTENTION ',LITER                      ACVI4941
      RETURN 1                                                          ACVI4942
C     RETURN                                                            ACVI4943
C                                                                       ACVI4944
      ENTRY ERROR(LITER)                                                ACVI4945
C     Exit program because of gross error.                              ACVI4946
      WRITE(6,'(//A,A)') ' ***** FATAL ERROR ',LITER                    ACVI4947
      STOP                                                              ACVI4948
C ---*end of ATTN*------------------------------------------------------ACVI4949
      END                                                               ACVI4950
C                                                                       ACVI4951
      SUBROUTINE WRDATA(IOTYPE,IOFILE,NUMI,IRAY,NF)                     ACVI4952
C                                                                       ACVI4953
C     ORIGINAL: JPD ('PHDRYR.RESKE2.FORT(IODATA)')                      ACVI4954
C     MODIFIED: CB (2/15/91) FOR GENERAL PURPOSE                        ACVI4955
C                                                                       ACVI4956
CB    IMPLICIT REAL*8(A-H,O-Z)                                          ACVI4957
      DIMENSION IRAY(NF:1)                                              ACVI4958
      IF (IOTYPE.EQ.0) THEN                                             ACVI4959
         WRITE(IOFILE)NUMI                                              ACVI4960
         WRITE(IOFILE)(IRAY(I),I=NF,0)                                  ACVI4961
      ELSE                                                              ACVI4962
         READ(IOFILE)NUMI                                               ACVI4963
         READ(IOFILE)(IRAY(I),I=NF,0)                                   ACVI4964
      ENDIF                                                             ACVI4965
C     SPLIT: LRECL=1724 & BLKSIZE=32760 --> 430 4 BYTE SETS PER RECORD  ACVI4966
C            (BUFFER SIZE IS 32760 SO SPLIT INTO 19 RECORDS WITH        ACVI4967
C             4 BYTES PER RECORD PLUS 4 BYTES PER BLOCK OVERHEAD:       ACVI4968
C             19*(1720+4)+4=32760 FOR OPTIMUM BUFFER UTILIZATION)       ACVI4969
      NRUN=NUMI                                                         ACVI4970
      NREC=NRUN/430                                                     ACVI4971
      NBEG=1                                                            ACVI4972
      NEND=430                                                          ACVI4973
C     I/O FULL BLOCKS, INTEGER ARRAY                                    ACVI4974
      DO 10 NR=1,NREC                                                   ACVI4975
         IF (IOTYPE.EQ.0) THEN                                          ACVI4976
            WRITE(IOFILE)(IRAY(N),N=NBEG,NEND)                          ACVI4977
         ELSE                                                           ACVI4978
            READ(IOFILE)(IRAY(N),N=NBEG,NEND)                           ACVI4979
         ENDIF                                                          ACVI4980
         NBEG=NEND+1                                                    ACVI4981
10       NEND=NEND+430                                                  ACVI4982
C     I/O RESIDUAL                                                      ACVI4983
      IF (NBEG.LE.NRUN) THEN                                            ACVI4984
         IF (IOTYPE.EQ.0) THEN                                          ACVI4985
            WRITE(IOFILE)(IRAY(N),N=NBEG,NRUN)                          ACVI4986
         ELSE                                                           ACVI4987
            READ(IOFILE)(IRAY(N),N=NBEG,NRUN)                           ACVI4988
         ENDIF                                                          ACVI4989
      ENDIF                                                             ACVI4990
      RETURN                                                            ACVI4991
      END                                                               ACVI4992
C ----------------------------------------------------------------------ACVI4993
C                                                                       ACVI4994
C                       ***********************                         ACVI4995
C                       ***   WST PACKAGE   ***                         ACVI4996
C                       ***********************                         ACVI4997
C                                                                       ACVI4998
C                -------------------------------------                  ACVI4999
C                                                                       ACVI5000
C             *******************************************               ACVI5001
C             ***    WEIGHTED SEARCH TREE ROUTINES    ***               ACVI5002
C             ***                for                  ***               ACVI5003
C             ***  SCIENTIFIC (FORTRAN) APPLICATIONS  ***               ACVI5004
C             *******************************************               ACVI5005
C                                                                       ACVI5006
C ----------------------------------------------------------------------ACVI5007
C                                                                       ACVI5008
C Authors: Soon Park, C. Bahri, J. P. Draayer and S.-Q. Zheng           ACVI5009
C          Department of Physics (Computer Science)                     ACVI5010
C          Louisiana State University                                   ACVI5011
C          Baton Rouge LA                                               ACVI5012
C          USA 70803-4001                                               ACVI5013
C                                                                       ACVI5014
C                 BITNET:  PHDRYR @ LSUMVS or LSUVM                     ACVI5015
C                 TELEX:   559184                                       ACVI5016
C                 PHONE:   USA-504-388-2261                             ACVI5017
C                 FAX:     USA-504-388-5855                             ACVI5018
C                                                                       ACVI5019
C Version: 1.1    LSU (07/01/91)                                        ACVI5020
C Version: 2.1    LSU (04/01/94)                                        ACVI5021
C                                                                       ACVI5022
C ----------------------------------------------------------------------ACVI5023
C                                                                       ACVI5024
C Updates: 07/90 Original from a FORTRAN code written by Soon Park      ACVI5025
C          04/94 Inclusion of TMRG subroutine by C. Bahri               ACVI5026
C                                                                       ACVI5027
C ----------------------------------------------------------------------ACVI5028
C                                                                       ACVI5029
C General comments on the package:                                      ACVI5030
C                                                                       ACVI5031
C   This package is written in FORTRAN since most scientific programs   ACVI5032
C   require FORTRAN compatibility and for it the existing scientific    ACVI5033
C   subroutine libraries are the most extensive and efficient.          ACVI5034
C                                                                       ACVI5035
C   The tree is a linear array ID(-10:*) consisting of eleven integers, ACVI5036
C   ID(-10:0), that specifies the structure of the tree, see TSET for   ACVI5037
C   documentation on this, and node information starting with ID(1:1).  ACVI5038
C   The latter includes the key(s), data, balance factor and priority,  ACVI5039
C   as well as the left and right child pointers. See the documentation ACVI5040
C   on each subroutine for further details.                             ACVI5041
C                                                                       ACVI5042
C ----------------------------------------------------------------------ACVI5043
C                                                                       ACVI5044
C This numerical database package consists of 7 different subroutines:  ACVI5045
C                                                                       ACVI5046
C   1. TSET    -->  initializes a storage area for use as a binary tree ACVI5047
C       TSETLL -->  ... entry in TSET for a tree with a linked-list     ACVI5048
C       TSETLF -->  ... entry in TSET for a tree without a linked-list  ACVI5049
C   2. TCHK    -->  performs the search operation on the data structure ACVI5050
C   3. TADD    -->  add a new element to an existing WST data structure ACVI5051
C   4. TINS    -->  called by TADD to insert new node into existing WST ACVI5052
C   5. TDEL    -->  called by TADD to delete low priority node in a WST ACVI5053
C   6. TOUT    -->  generates output information on specific tree nodes ACVI5054
C   7. TMRG    -->  merges two trees that have the same structure       ACVI5055
C                                                                       ACVI5056
C The four routines TSET, TCHK, TADD and TOUT subject to user control   ACVI5057
C while TINS and TDEL are only used by TADD and not otherwise needed.   ACVI5058
C The type of information stored in the buffer may vary from applicaton ACVI5059
C to application. This can be achieved by editing TADD. For example,    ACVI5060
C if variable length integer data is to be stored, then BUFFER must be  ACVI5061
C replaced by, INTGER, an integer array throughout. Likewise if BUFFER  ACVI5062
C holds double precision or complex data, the statement that defines    ACVI5063
C BUFFER in TADD must be changed accordingly. The program can also be   ACVI5064
C used in other ways, for example, each node can refer to more than a   ACVI5065
C single buffer. If this is done, both TSET and TADD must be modified   ACVI5066
C in a rather obvious way to add the required multiple link-list and    ACVI5067
C buffer arrays.                                                        ACVI5068
C                                                                       ACVI5069
C ----------------------------------------------------------------------ACVI5070
C                                                                       ACVI5071
C                           **************                              ACVI5072
C                           ***  TSET  ***                              ACVI5073
C                           **************                              ACVI5074
C                                                                       ACVI5075
C The subroutine TSET must be called before inserting the first item    ACVI5076
C into the tree. This call fixes the first 11 values of the array ID:   ACVI5077
C                                                                       ACVI5078
C      ID(-10): Number of nodes currently in the tree    ==> ID(-10)    ACVI5079
C      ID(-9) : Maximum nodes in in tree                 ==> MXNODE     ACVI5080
C      ID(-8) : Parent node of ID(-7), see next entry    ==> NF         ACVI5081
C      ID(-7) : Node to be balanced                      ==> NA         ACVI5082
C      ID(-6) : Parent node pointer                      ==> NQ         ACVI5083
C      ID(-5) : Current node pointer                     ==> Determined ACVI5084
C      ID(-4) : Number of integers assigned the key      ==> NKEY       ACVI5085
C      ID(-3) : Number of integers for key and data      ==> NSUM       ACVI5086
C      ID(-2) : Position of the priority in a node       ==> NPR        ACVI5087
C      ID(-1) : Position of the next available node      ==> Determined ACVI5088
C      ID( 0) : Root pointer (-1 for empty tree)         ==> Determined ACVI5089
C                                                                       ACVI5090
C The subroutine TSET also initializes all entries of the array LLBUFF  ACVI5091
C where the link-list information for BUFFER is stored. In particular,  ACVI5092
C the first three which refer to the free space are assigned values:    ACVI5093
C                                                                       ACVI5094
C      LLBUFF(-2) : Tail of free space                                  ACVI5095
C      LLBUFF(-1) : Head of free space                                  ACVI5096
C      LLBUFF( 0) : Size of free space                                  ACVI5097
C                                                                       ACVI5098
C ----------------------------------------------------------------------ACVI5099
C                                                                       ACVI5100
      SUBROUTINE TSET                                                   ACVI5101
C                                                                       ACVI5102
      INTEGER ID(-10:*),LLBUFF(-2:*)                                    ACVI5103
C                                                                       ACVI5104
      ENTRY TSETLL(ID,MXNODE,NKEY,NDAT,LLBUFF,MXBUFF)                   ACVI5105
C                                                                       ACVI5106
C Initialize link-list array values                                     ACVI5107
C                                                                       ACVI5108
      DO 10 J=1,MXBUFF-1                                                ACVI5109
10       LLBUFF(J)=J+1                                                  ACVI5110
      LLBUFF(MXBUFF)=-1                                                 ACVI5111
C                                                                       ACVI5112
C Initialize free space pointer/counter                                 ACVI5113
C                                                                       ACVI5114
      LLBUFF(-2)=MXBUFF                                                 ACVI5115
      LLBUFF(-1)=1                                                      ACVI5116
      LLBUFF(0)=MXBUFF                                                  ACVI5117
C                                                                       ACVI5118
      ENTRY TSETLF(ID,MXNODE,NKEY,NDAT)                                 ACVI5119
C                                                                       ACVI5120
C Initialize WST parameters                                             ACVI5121
C                                                                       ACVI5122
      NSUM=NKEY+NDAT                                                    ACVI5123
      IMAX=MXNODE*(NSUM+3)                                              ACVI5124
      ID(-10)=0                                                         ACVI5125
      ID(-9)=MXNODE                                                     ACVI5126
      ID(-8)=-1                                                         ACVI5127
      ID(-7)=-1                                                         ACVI5128
      ID(-6)=-1                                                         ACVI5129
      ID(-5)=0                                                          ACVI5130
      ID(-4)=NKEY                                                       ACVI5131
      ID(-3)=NSUM                                                       ACVI5132
      ID(-2)=NSUM+1                                                     ACVI5133
      ID(-1)=0                                                          ACVI5134
      ID(0)=-1                                                          ACVI5135
      DO 30 I=1,IMAX                                                    ACVI5136
30       ID(I)=0                                                        ACVI5137
      RETURN                                                            ACVI5138
      END                                                               ACVI5139
C ----------------------------------------------------------------------ACVI5140
C                                                                       ACVI5141
C                            ************                               ACVI5142
C                            *** TCHK ***                               ACVI5143
C                            ************                               ACVI5144
C                                                                       ACVI5145
C The subroutine TCHK constructs a weighted search tree or locates a    ACVI5146
C specific node in a weighted search tree. The subroutine uses two      ACVI5147
C arrays, NEWKEY and ID. A call to TCHK(,,) generates a search of ID forACVI5148
C the key stored in NEWKEY. If the search is successful the priority of ACVI5149
C the node will be increased. If not successful, a normal return is     ACVI5150
C generated and the position, ID(-5) where a new node can be inserted   ACVI5151
C is given. After this subroutine is called, the position(s) where the  ACVI5152
C key(s) of the new node are to be assigned are ID(-5)+1, ID(-5)+2,     ACVI5153
C ID(-5)+3, ... and the positions of the location and size of the new   ACVI5154
C incoming data are ID(IPOS)+1 and ID(IPOS)+2, respectively, where      ACVI5155
C IPOS = ID(-5)+ID(-4) and ID(-4) is the number of integer words set    ACVI5156
C aside for the key.                                                    ACVI5157
C                                                                       ACVI5158
C   RETURN 1 --> Successful search (key already in the tree,            ACVI5159
C                ID(-5) is set to point to the located node).           ACVI5160
C                The priority of the located node is updated.           ACVI5161
C                                                                       ACVI5162
C ----------------------------------------------------------------------ACVI5163
C                                                                       ACVI5164
      SUBROUTINE TCHK(NEWKEY,ID,*)                                      ACVI5165
C                                                                       ACVI5166
      INTEGER NEWKEY(*),ID(-10:*)                                       ACVI5167
      LOGICAL FLAG                                                      ACVI5168
C                                                                       ACVI5169
C Integer data for a DEC system                                         ACVI5170
      DATA NF0,N2,N3/1073741823,-2147483648,-1073741824/                ACVI5171
C Hexadecimal data for a IBM system                                     ACVI5172
C                                                                       ACVI5173
C      DATA NF0,N2,N3/Z3FFFFFFF,Z80000000,ZC0000000/                     ACVI5174
      NKEY=ID(-4)                                                       ACVI5175
      NPR=ID(-2)                                                        ACVI5176
      NLC=NPR+1                                                         ACVI5177
      NRC=NPR+2                                                         ACVI5178
C                                                                       ACVI5179
C Special case (empty tree)                                             ACVI5180
C                                                                       ACVI5181
      IF (ID(0).EQ.-1) THEN                                             ACVI5182
         ID(-5)=0                                                       ACVI5183
         RETURN                                                         ACVI5184
      ENDIF                                                             ACVI5185
C                                                                       ACVI5186
C Normal case (non-empty tree)                                          ACVI5187
C                                                                       ACVI5188
      NF=-1                                                             ACVI5189
      NA=ID(0)                                                          ACVI5190
      NP=ID(0)                                                          ACVI5191
      NQ=-1                                                             ACVI5192
100   IF (NP.NE.-1) THEN                                                ACVI5193
C                                                                       ACVI5194
C Check if the balance factor of NP is 0 or not                         ACVI5195
C                                                                       ACVI5196
         IF (ID(NP+NPR).LT.0) THEN                                      ACVI5197
            NA=NP                                                       ACVI5198
            NF=NQ                                                       ACVI5199
         ENDIF                                                          ACVI5200
         DO 101 I=1,NKEY                                                ACVI5201
            IF (NEWKEY(I).LT.ID(NP+I)) THEN                             ACVI5202
               NQ=NP                                                    ACVI5203
               NO=ID(NP+NLC)                                            ACVI5204
               IF (NO.EQ.-1) THEN                                       ACVI5205
                  NP=-1                                                 ACVI5206
               ELSEIF (ID(NO+NRC).EQ.NP) THEN                           ACVI5207
                  IF (ISHFT(ID(NP+NPR),-30).EQ.2) THEN                  ACVI5208
                     NP=NO                                              ACVI5209
                  ELSE                                                  ACVI5210
                     NP=-1                                              ACVI5211
                  ENDIF                                                 ACVI5212
               ELSE                                                     ACVI5213
                  NP=NO                                                 ACVI5214
               ENDIF                                                    ACVI5215
               GOTO 100                                                 ACVI5216
            ELSEIF (NEWKEY(I).GT.ID(NP+I)) THEN                         ACVI5217
               NQ=NP                                                    ACVI5218
               NO=ID(NP+NLC)                                            ACVI5219
               IF (NO.EQ.-1) THEN                                       ACVI5220
                  NP=-1                                                 ACVI5221
               ELSEIF (ID(NO+NRC).EQ.NP) THEN                           ACVI5222
                  IF (ISHFT(ID(NP+NPR),-30).EQ.2) THEN                  ACVI5223
                     NP=-1                                              ACVI5224
                  ELSE                                                  ACVI5225
                     NP=NO                                              ACVI5226
                  ENDIF                                                 ACVI5227
               ELSE                                                     ACVI5228
                  NP=ID(NO+NRC)                                         ACVI5229
               ENDIF                                                    ACVI5230
               GOTO 100                                                 ACVI5231
            ENDIF                                                       ACVI5232
101      CONTINUE                                                       ACVI5233
C                                                                       ACVI5234
C If a match is found, increase the node priority, reconstruct the      ACVI5235
C linear array, set the pointer to the node location, and RETURN 1      ACVI5236
C                                                                       ACVI5237
         NA=NP                                                          ACVI5238
         FLAG=.TRUE.                                                    ACVI5239
         ITEM=IAND(ID(NP+NPR),NF0)                                      ACVI5240
         ITEMP=ITEM+ISHFT(ISHFT(ITEM,24),-16)                           ACVI5241
C                                                                       ACVI5242
C Saturation condition for the priority                                 ACVI5243
C                                                                       ACVI5244
         IF (ITEMP.LE.NF0) THEN                                         ACVI5245
            ID(NP+NPR)=ITEMP+IAND(ID(NP+NPR),N3)                        ACVI5246
            ITEMPR=ITEMP                                                ACVI5247
         ELSE                                                           ACVI5248
            ITEMPR=ITEM                                                 ACVI5249
         ENDIF                                                          ACVI5250
         NL=ID(-1)                                                      ACVI5251
102      NB=NA+NA+NRC                                                   ACVI5252
         IF (NB.LT.NL) THEN                                             ACVI5253
            NC=NB+NRC                                                   ACVI5254
            NBPR=IAND(ID(NB+NPR),NF0)                                   ACVI5255
            NCPR=IAND(ID(NC+NPR),NF0)                                   ACVI5256
            IF ((NBPR.LE.NCPR).OR.(NC.GE.NL)) THEN                      ACVI5257
               IF (NBPR.LT.ITEMPR) THEN                                 ACVI5258
                  IF (FLAG) THEN                                        ACVI5259
                     NARC=ID(NA+NRC)                                    ACVI5260
                     NALC=ID(NA+NLC)                                    ACVI5261
                     IF (NARC.NE.-1) THEN                               ACVI5262
                        IF (ID(NARC+NLC).EQ.NA) THEN                    ACVI5263
                           ID(NARC+NLC)=NL                              ACVI5264
                        ELSE                                            ACVI5265
                           IF (ID(ID(NARC+NRC)+NLC).EQ.NA) THEN         ACVI5266
                              ID(ID(NARC+NRC)+NLC)=NL                   ACVI5267
                           ELSE                                         ACVI5268
                              ID(ID(NARC+NLC)+NRC)=NL                   ACVI5269
                           ENDIF                                        ACVI5270
                        ENDIF                                           ACVI5271
                     ENDIF                                              ACVI5272
                     IF (NALC.NE.-1) THEN                               ACVI5273
                        IF (ID(NALC+NRC).EQ.NA) THEN                    ACVI5274
                           ID(NALC+NRC)=NL                              ACVI5275
                        ELSE                                            ACVI5276
                           ID(ID(NALC+NRC)+NRC)=NL                      ACVI5277
                        ENDIF                                           ACVI5278
                     ENDIF                                              ACVI5279
                     DO 105 I=1,NRC                                     ACVI5280
105                     ID(NL+I)=ID(NA+I)                               ACVI5281
                     FLAG=.FALSE.                                       ACVI5282
                  ENDIF                                                 ACVI5283
                  NBRC=ID(NB+NRC)                                       ACVI5284
                  NBLC=ID(NB+NLC)                                       ACVI5285
                  IF (NBRC.EQ.-1) THEN                                  ACVI5286
                     ID(0)=NA                                           ACVI5287
                  ELSE                                                  ACVI5288
                     IF (ID(NBRC+NLC).EQ.NB) THEN                       ACVI5289
                        ID(NBRC+NLC)=NA                                 ACVI5290
                     ELSE                                               ACVI5291
                        IF (ID(ID(NBRC+NRC)+NLC).EQ.NB) THEN            ACVI5292
                           ID(ID(NBRC+NRC)+NLC)=NA                      ACVI5293
                        ELSE                                            ACVI5294
                           ID(ID(NBRC+NLC)+NRC)=NA                      ACVI5295
                        ENDIF                                           ACVI5296
                     ENDIF                                              ACVI5297
                  ENDIF                                                 ACVI5298
                  IF (NBLC.NE.-1) THEN                                  ACVI5299
                     IF (ID(NBLC+NRC).EQ.NB) THEN                       ACVI5300
                        ID(NBLC+NRC)=NA                                 ACVI5301
                     ELSE                                               ACVI5302
                        ID(ID(NBLC+NRC)+NRC)=NA                         ACVI5303
                     ENDIF                                              ACVI5304
                  ENDIF                                                 ACVI5305
                  DO 109 I=1,NRC                                        ACVI5306
                     ID(NA+I)=ID(NB+I)                                  ACVI5307
109               CONTINUE                                              ACVI5308
                  NA=NB                                                 ACVI5309
                  GOTO 102                                              ACVI5310
               ENDIF                                                    ACVI5311
            ELSE                                                        ACVI5312
               IF (NCPR.LT.ITEMPR) THEN                                 ACVI5313
                  IF (FLAG) THEN                                        ACVI5314
                     NARC=ID(NA+NRC)                                    ACVI5315
                     NALC=ID(NA+NLC)                                    ACVI5316
                     IF (NARC.NE.-1) THEN                               ACVI5317
                        IF (ID(NARC+NLC).EQ.NA) THEN                    ACVI5318
                           ID(NARC+NLC)=NL                              ACVI5319
                        ELSE                                            ACVI5320
                           IF (ID(ID(NARC+NRC)+NLC).EQ.NA) THEN         ACVI5321
                              ID(ID(NARC+NRC)+NLC)=NL                   ACVI5322
                           ELSE                                         ACVI5323
                              ID(ID(NARC+NLC)+NRC)=NL                   ACVI5324
                           ENDIF                                        ACVI5325
                        ENDIF                                           ACVI5326
                     ENDIF                                              ACVI5327
                     IF (NALC.NE.-1) THEN                               ACVI5328
                        IF (ID(NALC+NRC).EQ.NA) THEN                    ACVI5329
                           ID(NALC+NRC)=NL                              ACVI5330
                        ELSE                                            ACVI5331
                           ID(ID(NALC+NRC)+NRC)=NL                      ACVI5332
                        ENDIF                                           ACVI5333
                     ENDIF                                              ACVI5334
                     DO 112 I=1,NRC                                     ACVI5335
112                     ID(NL+I)=ID(NA+I)                               ACVI5336
                     FLAG=.FALSE.                                       ACVI5337
                  ENDIF                                                 ACVI5338
                  NCRC=ID(NC+NRC)                                       ACVI5339
                  NCLC=ID(NC+NLC)                                       ACVI5340
                  IF (NCRC.EQ.-1) THEN                                  ACVI5341
                     ID(0)=NA                                           ACVI5342
                  ELSE                                                  ACVI5343
                     IF (ID(NCRC+NLC).EQ.NC) THEN                       ACVI5344
                        ID(NCRC+NLC)=NA                                 ACVI5345
                     ELSE                                               ACVI5346
                        IF (ID(ID(NCRC+NRC)+NLC).EQ.NC) THEN            ACVI5347
                           ID(ID(NCRC+NRC)+NLC)=NA                      ACVI5348
                        ELSE                                            ACVI5349
                           ID(ID(NCRC+NLC)+NRC)=NA                      ACVI5350
                        ENDIF                                           ACVI5351
                     ENDIF                                              ACVI5352
                  ENDIF                                                 ACVI5353
                  IF (NCLC.NE.-1) THEN                                  ACVI5354
                     IF (ID(NCLC+NRC).EQ.NC) THEN                       ACVI5355
                        ID(NCLC+NRC)=NA                                 ACVI5356
                     ELSE                                               ACVI5357
                        ID(ID(NCLC+NRC)+NRC)=NA                         ACVI5358
                     ENDIF                                              ACVI5359
                  ENDIF                                                 ACVI5360
                  DO 116 I=1,NRC                                        ACVI5361
                     ID(NA+I)=ID(NC+I)                                  ACVI5362
116               CONTINUE                                              ACVI5363
                  NA=NC                                                 ACVI5364
                  GOTO 102                                              ACVI5365
               ENDIF                                                    ACVI5366
            ENDIF                                                       ACVI5367
         ENDIF                                                          ACVI5368
118      IF (FLAG) GOTO 124                                             ACVI5369
         NLRC=ID(NL+NRC)                                                ACVI5370
         NLLC=ID(NL+NLC)                                                ACVI5371
         IF (NLRC.EQ.-1) THEN                                           ACVI5372
            ID(0)=NA                                                    ACVI5373
         ELSE                                                           ACVI5374
            IF (ID(NLRC+NLC).EQ.NL) THEN                                ACVI5375
               ID(NLRC+NLC)=NA                                          ACVI5376
            ELSE                                                        ACVI5377
               IF (ID(ID(NLRC+NRC)+NLC).EQ.NL) THEN                     ACVI5378
                  ID(ID(NLRC+NRC)+NLC)=NA                               ACVI5379
               ELSE                                                     ACVI5380
                  ID(ID(NLRC+NLC)+NRC)=NA                               ACVI5381
               ENDIF                                                    ACVI5382
            ENDIF                                                       ACVI5383
         ENDIF                                                          ACVI5384
         IF (NLLC.NE.-1) THEN                                           ACVI5385
            IF (ID(NLLC+NRC).EQ.NL) THEN                                ACVI5386
               ID(NLLC+NRC)=NA                                          ACVI5387
            ELSE                                                        ACVI5388
               ID(ID(NLLC+NRC)+NRC)=NA                                  ACVI5389
            ENDIF                                                       ACVI5390
         ENDIF                                                          ACVI5391
         DO 120 I=1,NRC                                                 ACVI5392
            ID(NA+I)=ID(NL+I)                                           ACVI5393
            ID(NL+I)=0                                                  ACVI5394
120      CONTINUE                                                       ACVI5395
124      ID(-5)=NA                                                      ACVI5396
         RETURN 1                                                       ACVI5397
      ENDIF                                                             ACVI5398
      ID(-8)=NF                                                         ACVI5399
      ID(-7)=NA                                                         ACVI5400
      ID(-6)=NQ                                                         ACVI5401
      ID(-5)=ID(-1)+NRC                                                 ACVI5402
      RETURN                                                            ACVI5403
      END                                                               ACVI5404
C ----------------------------------------------------------------------ACVI5405
C                                                                       ACVI5406
C                            ************                               ACVI5407
C                            *** TADD ***                               ACVI5408
C                            ************                               ACVI5409
C                                                                       ACVI5410
C The subroutine TADD allows the user to add elements to an existing    ACVI5411
C database provided space is available. If it is not and the priority   ACVI5412
C of the incoming element is greater than the lowest priority element   ACVI5413
C currently in the database that element is eliminated to make room for ACVI5414
C the new one. The last scenario is repeated as necessary to accommodateACVI5415
C the new element. The logic of the program is the following:           ACVI5416
C                                                                       ACVI5417
C Step 1) check if either the tree or buffer is full                    ACVI5418
C         if yes, check if the priority of incoming item is higher      ACVI5419
C                 than that of the lowest priority node in the tree     ACVI5420
C                 if yes, delete lowest priority node and goto Step 1)  ACVI5421
C                 else RETURN                                           ACVI5422
C         else    goto Step 2)                                          ACVI5423
C Step 2) insert the incoming item into the database                    ACVI5424
C                                                                       ACVI5425
C ----------------------------------------------------------------------ACVI5426
C                                                                       ACVI5427
      SUBROUTINE TADD(NEWKEY,NEWDAT,BULOAD,NOSIZE,NPBASE,ID,BUFFER,     ACVI5428
     &                LLBUFF)                                           ACVI5429
C                                                                       ACVI5430
      INTEGER NEWKEY(*),NEWDAT(*),ID(-10:*),LLBUFF(-2:*)                ACVI5431
      REAL    BULOAD(*),BUFFER(*)                     ! *** storage typeACVI5432
C     REAL*8  BULOAD(*),BUFFER(*)                     ! *** storage typeACVI5433
      DATA    NF0/1073741823/                                           ACVI5434
C                                                                       ACVI5435
C ...define the number of keys in a node                                ACVI5436
C                                                                       ACVI5437
      NKEY=ID(-4)                                                       ACVI5438
C                                                                       ACVI5439
C ...calculate current priority                                         ACVI5440
C                                                                       ACVI5441
      NPCURR=NPBASE+ISHFT(NPBASE,8)                                     ACVI5442
C                                                                       ACVI5443
C ...check to see if the tree and buffer can receive the new item       ACVI5444
C                                                                       ACVI5445
      IF ((ID(-10).LE.ID(-9)).AND.(LLBUFF(0).GE.NOSIZE)) THEN           ACVI5446
         GOTO 300                                                       ACVI5447
C                                                                       ACVI5448
C ...check if priority of new item is greater than lowest priority      ACVI5449
C                                                                       ACVI5450
      ELSE                                                              ACVI5451
100      IF (NPCURR.GT.IAND(ID(ID(-2)),NF0)) THEN                       ACVI5452
C                                                                       ACVI5453
C Rearrange LLBUFF array which holds free space list information        ACVI5454
C                                                                       ACVI5455
            IX=ID(NKEY+1)                                               ACVI5456
            IF (LLBUFF(0).EQ.0) THEN                                    ACVI5457
               LLBUFF(0)=ID(NKEY+2)                                     ACVI5458
               LLBUFF(-1)=IX                                            ACVI5459
            ELSE                                                        ACVI5460
               LLBUFF(0)=LLBUFF(0)+ID(NKEY+2)                           ACVI5461
               LLBUFF(LLBUFF(-2))=IX                                    ACVI5462
            ENDIF                                                       ACVI5463
            LLBUFF(-2)=IX                                               ACVI5464
            DO 200 I=1,ID(NKEY+2)-1                                     ACVI5465
               IX=LLBUFF(IX)                                            ACVI5466
200            LLBUFF(-2)=IX                                            ACVI5467
C                                                                       ACVI5468
            CALL TDEL(ID)                                               ACVI5469
C                                                                       ACVI5470
            IF (LLBUFF(0).LT.NOSIZE) GOTO 100                           ACVI5471
            CALL TCHK(NEWKEY,ID)                                        ACVI5472
            GOTO 300                                                    ACVI5473
C                                                                       ACVI5474
C ...doing nothing when new priority is less than lowest priority       ACVI5475
C                                                                       ACVI5476
         ELSE                                                           ACVI5477
            RETURN                                                      ACVI5478
         ENDIF                                                          ACVI5479
      ENDIF                                                             ACVI5480
300   CONTINUE                                                          ACVI5481
C                                                                       ACVI5482
C Load the key, data and its priority into the tree                     ACVI5483
C                                                                       ACVI5484
C ...last key position of the current node                              ACVI5485
C                                                                       ACVI5486
      IFIND=ID(-5)+NKEY                                                 ACVI5487
C                                                                       ACVI5488
C ...load integer data into tree after the last key                     ACVI5489
C                                                                       ACVI5490
      NEWDAT(1)=LLBUFF(-1)                                              ACVI5491
      NEWDAT(2)=NOSIZE                                                  ACVI5492
      NDAT=ID(-3)-ID(-4)                                                ACVI5493
      DO 400 IKK=1,NDAT                                                 ACVI5494
400   ID(IFIND+IKK)=NEWDAT(IKK)                                         ACVI5495
C                                                                       ACVI5496
C ...load real data into the BUFFER and update LLBUFF                   ACVI5497
C                                                                       ACVI5498
      DO 500 IKK=1,NOSIZE                                               ACVI5499
         BUFFER(LLBUFF(-1))=BULOAD(IKK)                                 ACVI5500
500      LLBUFF(-1)=LLBUFF(LLBUFF(-1))                                  ACVI5501
C                                                                       ACVI5502
C ...reduce size of free space                                          ACVI5503
C                                                                       ACVI5504
      LLBUFF(0)=LLBUFF(0)-NOSIZE                                        ACVI5505
C                                                                       ACVI5506
C ...load the priority                                                  ACVI5507
C                                                                       ACVI5508
      IFIND=ID(-5)+ID(-2)                                               ACVI5509
      ID(IFIND)=NPCURR                                                  ACVI5510
C                                                                       ACVI5511
C ...call TINS to complete (balance tree and restructure heap)          ACVI5512
C                                                                       ACVI5513
      CALL TINS(NEWKEY,ID)                                              ACVI5514
      RETURN                                                            ACVI5515
      END                                                               ACVI5516
C ----------------------------------------------------------------------ACVI5517
C                                                                       ACVI5518
C                            ************                               ACVI5519
C                            *** TINS ***                               ACVI5520
C                            ************                               ACVI5521
C                                                                       ACVI5522
C The subroutine TINS is used to insert a new item into a WST. A call   ACVI5523
C to TINS should follow a call to TCHK since TCHK tell where to insert  ACVI5524
C the new item. After the item is placed in the WST, TINS rebalances theACVI5525
C the tree and reconstruct the priority heap.                           ACVI5526
C                                                                       ACVI5527
C ----------------------------------------------------------------------ACVI5528
C                                                                       ACVI5529
      SUBROUTINE TINS(NEWKEY,ID)                                        ACVI5530
C                                                                       ACVI5531
      INTEGER NEWKEY(*),ID(-10:*)                                       ACVI5532
C                                                                       ACVI5533
C Integer data for a DEC system                                         ACVI5534
      DATA NF0,N2,N3/1073741823,-2147483648,-1073741824/                ACVI5535
C Hexadecimal data for a IBM  system                                    ACVI5536
C                                                                       ACVI5537
C      DATA NF0,N2,N3/Z3FFFFFFF,Z80000000,ZC0000000/                     ACVI5538
      ID(-10)=ID(-10)+1                                                 ACVI5539
      NKEY=ID(-4)                                                       ACVI5540
      NPR=ID(-2)                                                        ACVI5541
      NLC=NPR+1                                                         ACVI5542
      NRC=NPR+2                                                         ACVI5543
C                                                                       ACVI5544
C Special case (empty tree)                                             ACVI5545
C                                                                       ACVI5546
      IF (ID(0).EQ.-1) THEN                                             ACVI5547
         ID(-1)=NRC                                                     ACVI5548
         ID(0)=0                                                        ACVI5549
C                                                                       ACVI5550
C The first node position                                               ACVI5551
C                                                                       ACVI5552
         ID(NRC)=-1                                                     ACVI5553
         ID(NLC)=-1                                                     ACVI5554
         DO 10 I=1,NKEY                                                 ACVI5555
10          ID(I)=NEWKEY(I)                                             ACVI5556
         RETURN                                                         ACVI5557
      ENDIF                                                             ACVI5558
C                                                                       ACVI5559
C Normal case (non-empty tree)                                          ACVI5560
C                                                                       ACVI5561
      NY=ID(-5)                                                         ACVI5562
      NQ=ID(-6)                                                         ACVI5563
      NA=ID(-7)                                                         ACVI5564
      NF=ID(-8)                                                         ACVI5565
C                                                                       ACVI5566
C Set the pointer to indicate the position of the new node and load key ACVI5567
C                                                                       ACVI5568
      DO 130 I=1,NKEY                                                   ACVI5569
130      ID(NY+I)=NEWKEY(I)                                             ACVI5570
      ID(NY+NLC)=-1                                                     ACVI5571
      IF (ID(NQ+NLC).EQ.-1) THEN                                        ACVI5572
         ID(NY+NRC)=NQ                                                  ACVI5573
         ID(NQ+NLC)=NY                                                  ACVI5574
      ELSE                                                              ACVI5575
         DO 131 I=1,NKEY                                                ACVI5576
            IF (NEWKEY(I).LT.ID(NQ+I)) THEN                             ACVI5577
               ID(NY+NRC)=ID(NQ+NLC)                                    ACVI5578
               ID(NQ+NLC)=NY                                            ACVI5579
               ID(NQ+NPR)=IAND(ID(NQ+NPR),NF0)                          ACVI5580
               GOTO 999                                                 ACVI5581
            ELSEIF (NEWKEY(I).GT.ID(NQ+I)) THEN                         ACVI5582
               ID(NY+NRC)=NQ                                            ACVI5583
               ID(ID(NQ+NLC)+NRC)=NY                                    ACVI5584
               ID(NQ+NPR)=IAND(ID(NQ+NPR),NF0)                          ACVI5585
               GOTO 999                                                 ACVI5586
            ENDIF                                                       ACVI5587
131      CONTINUE                                                       ACVI5588
      ENDIF                                                             ACVI5589
      DO 161 I=1,NKEY                                                   ACVI5590
         IF (NEWKEY(I).LT.ID(NA+I)) THEN                                ACVI5591
            NP=ID(NA+NLC)                                               ACVI5592
            NB=NP                                                       ACVI5593
            IND=N2                                                      ACVI5594
            GOTO 200                                                    ACVI5595
         ELSEIF (NEWKEY(I).GT.ID(NA+I)) THEN                            ACVI5596
            NO=ID(NA+NLC)                                               ACVI5597
            IF (ID(NO+NRC).EQ.NA) THEN                                  ACVI5598
               NP=NO                                                    ACVI5599
            ELSE                                                        ACVI5600
               NP=ID(NO+NRC)                                            ACVI5601
            ENDIF                                                       ACVI5602
            NB=NP                                                       ACVI5603
            IND=N3                                                      ACVI5604
            GOTO 200                                                    ACVI5605
         ENDIF                                                          ACVI5606
161   CONTINUE                                                          ACVI5607
200   IF (NP.NE.NY) THEN                                                ACVI5608
         DO 211 I=1,NKEY                                                ACVI5609
            IF (NEWKEY(I).LT.ID(NP+I)) THEN                             ACVI5610
               ID(NP+NPR)=IAND(ID(NP+NPR),NF0)+N2                       ACVI5611
               NP=ID(NP+NLC)                                            ACVI5612
               GOTO 200                                                 ACVI5613
            ELSEIF (NEWKEY(I).GT.ID(NP+I)) THEN                         ACVI5614
               ID(NP+NPR)=IAND(ID(NP+NPR),NF0)+N3                       ACVI5615
               NO=ID(NP+NLC)                                            ACVI5616
               IF (ID(NO+NRC).EQ.NP) THEN                               ACVI5617
                  NP=NO                                                 ACVI5618
               ELSE                                                     ACVI5619
                  NP=ID(NO+NRC)                                         ACVI5620
               ENDIF                                                    ACVI5621
               GOTO 200                                                 ACVI5622
            ENDIF                                                       ACVI5623
211      CONTINUE                                                       ACVI5624
      ENDIF                                                             ACVI5625
      IF (ID(NA+NPR).GE.0) THEN                                         ACVI5626
         ID(NA+NPR)=IAND(ID(NA+NPR),NF0)+IND                            ACVI5627
         GOTO 999                                                       ACVI5628
      ELSEIF (IAND(ID(NA+NPR),N3).NE.IND) THEN                          ACVI5629
         ID(NA+NPR)=IAND(ID(NA+NPR),NF0)                                ACVI5630
         GOTO 999                                                       ACVI5631
      ENDIF                                                             ACVI5632
      IF (IND.EQ.N2) THEN                                               ACVI5633
         IF (ISHFT(ID(NB+NPR),-30).EQ.2) THEN                           ACVI5634
            IF (ID(ID(NA+NLC)+NRC).EQ.NA) THEN                          ACVI5635
               ID(NB+NRC)=ID(NA+NRC)                                    ACVI5636
               ID(ID(NB+NLC)+NRC)=NA                                    ACVI5637
               ID(NA+NLC)=-1                                            ACVI5638
               ID(NA+NRC)=NB                                            ACVI5639
            ELSE                                                        ACVI5640
               NC=ID(ID(NB+NLC)+NRC)                                    ACVI5641
               ID(NC+NRC)=ID(NB+NRC)                                    ACVI5642
               ID(NB+NRC)=ID(NA+NRC)                                    ACVI5643
               ID(ID(NB+NLC)+NRC)=NA                                    ACVI5644
               ID(NA+NLC)=NC                                            ACVI5645
               ID(NA+NRC)=NB                                            ACVI5646
            ENDIF                                                       ACVI5647
            ID(NA+NPR)=IAND(ID(NA+NPR),NF0)                             ACVI5648
            ID(NB+NPR)=IAND(ID(NB+NPR),NF0)                             ACVI5649
         ELSE                                                           ACVI5650
            IF (ID(ID(NB+NLC)+NRC).EQ.NB) THEN                          ACVI5651
               NC=ID(NB+NLC)                                            ACVI5652
            ELSE                                                        ACVI5653
               NC=ID(ID(NB+NLC)+NRC)                                    ACVI5654
            ENDIF                                                       ACVI5655
            IF (ID(NC+NPR).GE.0) THEN                                   ACVI5656
               ID(NC+NRC)=ID(NA+NRC)                                    ACVI5657
               ID(NA+NLC)=-1                                            ACVI5658
               ID(NA+NRC)=NC                                            ACVI5659
               ID(NB+NLC)=-1                                            ACVI5660
               ID(NC+NLC)=NB                                            ACVI5661
               ID(NA+NPR)=IAND(ID(NA+NPR),NF0)                          ACVI5662
               ID(NB+NPR)=IAND(ID(NB+NPR),NF0)                          ACVI5663
            ELSE                                                        ACVI5664
               IF (ID(ID(NC+NLC)+NRC).EQ.NC) THEN                       ACVI5665
                  IF (ISHFT(ID(NC+NPR),-30).EQ.2) THEN                  ACVI5666
                     ID(NA+NLC)=ID(NB+NRC)                              ACVI5667
                     ID(ID(NB+NLC)+NRC)=ID(NC+NLC)                      ACVI5668
                     ID(ID(NC+NLC)+NRC)=NB                              ACVI5669
                  ELSE                                                  ACVI5670
                     ID(NA+NLC)=ID(NC+NLC)                              ACVI5671
                     ID(ID(NA+NLC)+NRC)=ID(NB+NRC)                      ACVI5672
                     ID(ID(NB+NLC)+NRC)=NB                              ACVI5673
                  ENDIF                                                 ACVI5674
               ELSE                                                     ACVI5675
                  ID(NA+NLC)=ID(ID(NC+NLC)+NRC)                         ACVI5676
                  ID(ID(NA+NLC)+NRC)=ID(NB+NRC)                         ACVI5677
                  ID(ID(NB+NLC)+NRC)=ID(NC+NLC)                         ACVI5678
                  ID(ID(NC+NLC)+NRC)=NB                                 ACVI5679
               ENDIF                                                    ACVI5680
               ID(NC+NLC)=NB                                            ACVI5681
               ID(NC+NRC)=ID(NA+NRC)                                    ACVI5682
               ID(NA+NRC)=NC                                            ACVI5683
               ID(NB+NRC)=NA                                            ACVI5684
               IF (ISHFT(ID(NC+NPR),-30).EQ.2) THEN                     ACVI5685
                  ID(NA+NPR)=IAND(ID(NA+NPR),NF0)+N3                    ACVI5686
                  ID(NB+NPR)=IAND(ID(NB+NPR),NF0)                       ACVI5687
               ELSE                                                     ACVI5688
                  ID(NA+NPR)=IAND(ID(NA+NPR),NF0)                       ACVI5689
                  ID(NB+NPR)=IAND(ID(NB+NPR),NF0)+N2                    ACVI5690
               ENDIF                                                    ACVI5691
               ID(NC+NPR)=IAND(ID(NC+NPR),NF0)                          ACVI5692
            ENDIF                                                       ACVI5693
            NB=NC                                                       ACVI5694
         ENDIF                                                          ACVI5695
      ELSE                                                              ACVI5696
         IF (ISHFT(ID(NB+NPR),-30).EQ.3) THEN                           ACVI5697
            IF (ID(ID(NA+NLC)+NRC).EQ.NA) THEN                          ACVI5698
               ID(NA+NLC)=-1                                            ACVI5699
               ID(NB+NRC)=ID(NA+NRC)                                    ACVI5700
               ID(NA+NRC)=ID(NB+NLC)                                    ACVI5701
               ID(NB+NLC)=NA                                            ACVI5702
            ELSE                                                        ACVI5703
               NC=ID(NB+NLC)                                            ACVI5704
               ID(NB+NLC)=NA                                            ACVI5705
               ID(NB+NRC)=ID(NA+NRC)                                    ACVI5706
               ID(ID(NA+NLC)+NRC)=NC                                    ACVI5707
               ID(NA+NRC)=ID(NC+NRC)                                    ACVI5708
               ID(NC+NRC)=NA                                            ACVI5709
            ENDIF                                                       ACVI5710
            ID(NA+NPR)=IAND(ID(NA+NPR),NF0)                             ACVI5711
            ID(NB+NPR)=IAND(ID(NB+NPR),NF0)                             ACVI5712
         ELSE                                                           ACVI5713
            NC=ID(NB+NLC)                                               ACVI5714
            IF (ID(NC+NPR).GE.0) THEN                                   ACVI5715
               ID(NC+NRC)=ID(NA+NRC)                                    ACVI5716
               ID(NC+NLC)=NA                                            ACVI5717
               ID(NA+NLC)=-1                                            ACVI5718
               ID(NA+NRC)=NB                                            ACVI5719
               ID(NB+NLC)=-1                                            ACVI5720
               ID(NB+NRC)=NC                                            ACVI5721
               ID(NA+NPR)=IAND(ID(NA+NPR),NF0)                          ACVI5722
               ID(NB+NPR)=IAND(ID(NB+NPR),NF0)                          ACVI5723
            ELSE                                                        ACVI5724
               IF (ID(ID(NC+NLC)+NRC).EQ.NC) THEN                       ACVI5725
                  IF (ISHFT(ID(NC+NPR),-30).EQ.3) THEN                  ACVI5726
                     ID(ID(NA+NLC)+NRC)=NA                              ACVI5727
                     ID(ID(NC+NLC)+NRC)=ID(NC+NRC)                      ACVI5728
                     ID(NB+NLC)=ID(NC+NLC)                              ACVI5729
                  ELSE                                                  ACVI5730
                     ID(ID(NA+NLC)+NRC)=ID(NC+NLC)                      ACVI5731
                     ID(ID(NC+NLC)+NRC)=NA                              ACVI5732
                     ID(NB+NLC)=ID(NC+NRC)                              ACVI5733
                  ENDIF                                                 ACVI5734
               ELSE                                                     ACVI5735
                  ID(ID(NA+NLC)+NRC)=ID(NC+NLC)                         ACVI5736
                  ID(NB+NLC)=ID(ID(NC+NLC)+NRC)                         ACVI5737
                  ID(ID(NB+NLC)+NRC)=ID(NC+NRC)                         ACVI5738
                  ID(ID(NC+NLC)+NRC)=NA                                 ACVI5739
               ENDIF                                                    ACVI5740
               ID(NC+NLC)=NA                                            ACVI5741
               ID(NC+NRC)=ID(NA+NRC)                                    ACVI5742
               ID(NA+NRC)=NB                                            ACVI5743
               ID(NB+NRC)=NC                                            ACVI5744
               IF (ISHFT(ID(NC+NPR),-30).EQ.3) THEN                     ACVI5745
                  ID(NA+NPR)=IAND(ID(NA+NPR),NF0)+N2                    ACVI5746
                  ID(NB+NPR)=IAND(ID(NB+NPR),NF0)                       ACVI5747
               ELSE                                                     ACVI5748
                  ID(NA+NPR)=IAND(ID(NA+NPR),NF0)                       ACVI5749
                  ID(NB+NPR)=IAND(ID(NB+NPR),NF0)+N3                    ACVI5750
               ENDIF                                                    ACVI5751
               ID(NC+NPR)=IAND(ID(NC+NPR),NF0)                          ACVI5752
            ENDIF                                                       ACVI5753
            NB=NC                                                       ACVI5754
         ENDIF                                                          ACVI5755
      ENDIF                                                             ACVI5756
      IF (NF.EQ.-1) THEN                                                ACVI5757
         ID(0)=NB                                                       ACVI5758
      ELSEIF (NA.EQ.ID(NF+NLC)) THEN                                    ACVI5759
         ID(NF+NLC)=NB                                                  ACVI5760
      ELSE                                                              ACVI5761
         ID(ID(NF+NLC)+NRC)=NB                                          ACVI5762
      ENDIF                                                             ACVI5763
999   CONTINUE                                                          ACVI5764
C                                                                       ACVI5765
C Reconstruct the priority array                                        ACVI5766
C                                                                       ACVI5767
      NL=ID(-1)                                                         ACVI5768
      NRC2=NRC+NRC                                                      ACVI5769
1000  IF (MOD(NL,NRC2).EQ.0) THEN                                       ACVI5770
         NA=NL/2-NRC                                                    ACVI5771
      ELSE                                                              ACVI5772
         NA=(NL-NRC)/2                                                  ACVI5773
      ENDIF                                                             ACVI5774
      IF (IAND(ID(NA+NPR),NF0).GT.IAND(ID(NY+NPR),NF0)) THEN            ACVI5775
         NARC=ID(NA+NRC)                                                ACVI5776
         NALC=ID(NA+NLC)                                                ACVI5777
         IF (NARC.EQ.-1) THEN                                           ACVI5778
            ID(0)=NL                                                    ACVI5779
         ELSE                                                           ACVI5780
            IF (ID(NARC+NLC).EQ.NA) THEN                                ACVI5781
               ID(NARC+NLC)=NL                                          ACVI5782
            ELSE                                                        ACVI5783
               IF (ID(ID(NARC+NRC)+NLC).EQ.NA) THEN                     ACVI5784
                  ID(ID(NARC+NRC)+NLC)=NL                               ACVI5785
               ELSE                                                     ACVI5786
                  ID(ID(NARC+NLC)+NRC)=NL                               ACVI5787
               ENDIF                                                    ACVI5788
            ENDIF                                                       ACVI5789
         ENDIF                                                          ACVI5790
         IF (NALC.NE.-1) THEN                                           ACVI5791
            IF (ID(NALC+NRC).EQ.NA) THEN                                ACVI5792
               ID(NALC+NRC)=NL                                          ACVI5793
            ELSE                                                        ACVI5794
               ID(ID(NALC+NRC)+NRC)=NL                                  ACVI5795
            ENDIF                                                       ACVI5796
         ENDIF                                                          ACVI5797
         DO 1400 I=1,NRC                                                ACVI5798
            ID(NL+I)=ID(NA+I)                                           ACVI5799
1400     CONTINUE                                                       ACVI5800
         NL=NA                                                          ACVI5801
         IF (NL.GT.0) GOTO 1000                                         ACVI5802
      ENDIF                                                             ACVI5803
      NYRC=ID(NY+NRC)                                                   ACVI5804
      NYLC=ID(NY+NLC)                                                   ACVI5805
      IF (NYRC.EQ.-1) THEN                                              ACVI5806
         ID(0)=NL                                                       ACVI5807
      ELSE                                                              ACVI5808
         IF (ID(NYRC+NLC).EQ.NY) THEN                                   ACVI5809
            ID(NYRC+NLC)=NL                                             ACVI5810
         ELSE                                                           ACVI5811
            IF (ID(ID(NYRC+NRC)+NLC).EQ.NY) THEN                        ACVI5812
               ID(ID(NYRC+NRC)+NLC)=NL                                  ACVI5813
            ELSE                                                        ACVI5814
               ID(ID(NYRC+NLC)+NRC)=NL                                  ACVI5815
            ENDIF                                                       ACVI5816
         ENDIF                                                          ACVI5817
      ENDIF                                                             ACVI5818
      IF (NYLC.NE.-1) THEN                                              ACVI5819
         IF (ID(NYLC+NRC).EQ.NY) THEN                                   ACVI5820
            ID(NYLC+NRC)=NL                                             ACVI5821
         ELSE                                                           ACVI5822
            ID(ID(NYLC+NRC)+NRC)=NL                                     ACVI5823
         ENDIF                                                          ACVI5824
      ENDIF                                                             ACVI5825
      DO 7000 I=1,NRC                                                   ACVI5826
         ID(NL+I)=ID(NY+I)                                              ACVI5827
7000  CONTINUE                                                          ACVI5828
      ID(-1)=ID(-1)+NRC                                                 ACVI5829
      ID(-5)=NL                                                         ACVI5830
      RETURN                                                            ACVI5831
      END                                                               ACVI5832
C ----------------------------------------------------------------------ACVI5833
C                                                                       ACVI5834
C                            ************                               ACVI5835
C                            *** TDEL ***                               ACVI5836
C                            ************                               ACVI5837
C                                                                       ACVI5838
C The subroutine TDEL deletes a node from a weighted search tree,       ACVI5839
C keeping the binary tree balanced. The node deleted from the tree      ACVI5840
C is the first element of the linear array ID. Since the linear array   ACVI5841
C also has a heap structure, the root of a heap, which is the first     ACVI5842
C element of an array, is the lowest priority element (minheap). The    ACVI5843
C subroutine, TDEL, also manages the free space in the array BUFFER.    ACVI5844
C The pointer array LLBUFF is used for this purpose.                    ACVI5845
C                                                                       ACVI5846
C ----------------------------------------------------------------------ACVI5847
C                                                                       ACVI5848
      SUBROUTINE TDEL(ID)                                               ACVI5849
C                                                                       ACVI5850
      INTEGER ID(-10:*)                                                 ACVI5851
      LOGICAL FLAG                                                      ACVI5852
C                                                                       ACVI5853
C Integer data for a DEC system                                         ACVI5854
      DATA NF0,N2,N3/1073741823,-2147483648,-1073741824/                ACVI5855
C Hexadecimal data for a IBM system                                     ACVI5856
C                                                                       ACVI5857
C      DATA NF0,N2,N3/Z3FFFFFFF,Z80000000,ZC0000000/                     ACVI5858
C                                                                       ACVI5859
C Initialize some integer constants                                     ACVI5860
C                                                                       ACVI5861
      NKEY=ID(-4)                                                       ACVI5862
      NPR=ID(-2)                                                        ACVI5863
      NLC=NPR+1                                                         ACVI5864
      NRC=NPR+2                                                         ACVI5865
C                                                                       ACVI5866
C NA keeps track of most recent node with BF(0)                         ACVI5867
C NF is the parent of NA                                                ACVI5868
C NQ follows NP through the tree                                        ACVI5869
C                                                                       ACVI5870
      NF=-1                                                             ACVI5871
      NA=ID(0)                                                          ACVI5872
      NP=ID(0)                                                          ACVI5873
      NQ=-1                                                             ACVI5874
      NR=-1                                                             ACVI5875
100   IF (NP.NE.-1) THEN                                                ACVI5876
C                                                                       ACVI5877
C Looking for the last node with BF(0) that is not a leaf               ACVI5878
C                                                                       ACVI5879
         IF ((ID(NP+NPR).GE.0).AND.(ID(NP+NLC).NE.-1)) THEN             ACVI5880
            NA=NP                                                       ACVI5881
            NF=NQ                                                       ACVI5882
         ELSEIF (NQ.NE.-1) THEN                                         ACVI5883
            IF (FLAG) THEN                                              ACVI5884
               NO=ID(NQ+NLC)                                            ACVI5885
               IF (ID(NO+NRC).EQ.NQ) THEN                               ACVI5886
                  NRCHILD=NO                                            ACVI5887
               ELSE                                                     ACVI5888
                  NRCHILD=ID(NO+NRC)                                    ACVI5889
               ENDIF                                                    ACVI5890
               IF ((ISHFT(ID(NQ+NPR),-30).EQ.3)                         ACVI5891
     *             .AND.(ID(NRCHILD+NPR).GE.0)) THEN                    ACVI5892
                  NA=NQ                                                 ACVI5893
                  NF=NR                                                 ACVI5894
               ENDIF                                                    ACVI5895
            ELSE                                                        ACVI5896
               IF ((ISHFT(ID(NQ+NPR),-30).EQ.2)                         ACVI5897
     *             .AND.(ID(ID(NQ+NLC)+NPR).GE.0)) THEN                 ACVI5898
                  NA=NQ                                                 ACVI5899
                  NF=NR                                                 ACVI5900
               ENDIF                                                    ACVI5901
            ENDIF                                                       ACVI5902
         ENDIF                                                          ACVI5903
         DO 101 I=1,NKEY                                                ACVI5904
            IF (ID(I).LT.ID(NP+I)) THEN                                 ACVI5905
               NR=NQ                                                    ACVI5906
               NQ=NP                                                    ACVI5907
               NP=ID(NP+NLC)                                            ACVI5908
               FLAG=.TRUE.                                              ACVI5909
               GOTO 100                                                 ACVI5910
            ELSEIF (ID(I).GT.ID(NP+I)) THEN                             ACVI5911
               NR=NQ                                                    ACVI5912
               NQ=NP                                                    ACVI5913
               NO=ID(NP+NLC)                                            ACVI5914
               IF (ID(NO+NRC).EQ.NP) THEN                               ACVI5915
                  NP=NO                                                 ACVI5916
               ELSE                                                     ACVI5917
                  NP=ID(NO+NRC)                                         ACVI5918
               ENDIF                                                    ACVI5919
               FLAG=.FALSE.                                             ACVI5920
               GOTO 100                                                 ACVI5921
            ENDIF                                                       ACVI5922
101      CONTINUE                                                       ACVI5923
C                                                                       ACVI5924
C A match is found                                                      ACVI5925
C                                                                       ACVI5926
         GOTO 130                                                       ACVI5927
      ENDIF                                                             ACVI5928
C                                                                       ACVI5929
C A match is not found                                                  ACVI5930
C                                                                       ACVI5931
      WRITE(6,*) 'TREE IS EMPTY.'                                       ACVI5932
      GOTO 9999                                                         ACVI5933
130   CONTINUE                                                          ACVI5934
      ID(-10)=ID(-10)-1                                                 ACVI5935
C                                                                       ACVI5936
C Matched node does not have a child                                    ACVI5937
C                                                                       ACVI5938
      IF (ID(NP+NLC).EQ.-1) THEN                                        ACVI5939
         IF (NQ.EQ.-1) THEN                                             ACVI5940
            ID(0)=-1                                                    ACVI5941
            ID(-1)=0                                                    ACVI5942
            GOTO 9999                                                   ACVI5943
         ELSEIF (ID(NP+NRC).EQ.NQ) THEN                                 ACVI5944
            IF (ID(NQ+NLC).EQ.NP) THEN                                  ACVI5945
               ID(NQ+NLC)=-1                                            ACVI5946
            ELSE                                                        ACVI5947
               ID(ID(NQ+NLC)+NRC)=NQ                                    ACVI5948
            ENDIF                                                       ACVI5949
         ELSE                                                           ACVI5950
            ID(NQ+NLC)=ID(NP+NRC)                                       ACVI5951
         ENDIF                                                          ACVI5952
C                                                                       ACVI5953
C Matched node has only one child                                       ACVI5954
C                                                                       ACVI5955
      ELSEIF (ID(ID(NP+NLC)+NRC).EQ.NP) THEN                            ACVI5956
         NO=ID(NP+NLC)                                                  ACVI5957
         IF (ID(NP+NRC).EQ.-1) THEN                                     ACVI5958
            ID(0)=ID(NO)                                                ACVI5959
            ID(NO+NRC)=-1                                               ACVI5960
            GOTO 999                                                    ACVI5961
         ELSE                                                           ACVI5962
            IF (ID(NQ+NLC).EQ.NP) THEN                                  ACVI5963
               ID(NQ+NLC)=NO                                            ACVI5964
            ELSE                                                        ACVI5965
               ID(ID(NQ+NLC)+NRC)=NO                                    ACVI5966
            ENDIF                                                       ACVI5967
            ID(NO+NRC)=ID(NP+NRC)                                       ACVI5968
         ENDIF                                                          ACVI5969
C                                                                       ACVI5970
C Match node has both children                                          ACVI5971
C                                                                       ACVI5972
      ELSE                                                              ACVI5973
         NI=NQ                                                          ACVI5974
         NH=NP                                                          ACVI5975
         NG=ID(NP+NLC)                                                  ACVI5976
         IF (ID(NP+NPR).GE.0) THEN                                      ACVI5977
            NA=NH                                                       ACVI5978
            NF=NI                                                       ACVI5979
         ELSEIF ((ISHFT(ID(NP+NPR),-30).EQ.3).AND.                      ACVI5980
     *           (ID(ID(NG+NRC)+NPR).GE.0)) THEN                        ACVI5981
            NA=NH                                                       ACVI5982
            NF=NI                                                       ACVI5983
         ENDIF                                                          ACVI5984
         IF (ID(NG+NLC).EQ.-1) THEN                                     ACVI5985
            ID(ID(NG+NRC)+NRC)=NG                                       ACVI5986
            ID(NG+NLC)=ID(NG+NRC)                                       ACVI5987
            ID(NG+NRC)=ID(NP+NRC)                                       ACVI5988
         ELSEIF (ID(ID(NG+NLC)+NRC).EQ.NG) THEN                         ACVI5989
           IF (ISHFT(ID(NG+NPR),-30).EQ.2) THEN                         ACVI5990
               ID(ID(NG+NLC)+NRC)=ID(NG+NRC)                            ACVI5991
               ID(ID(NG+NRC)+NRC)=NG                                    ACVI5992
               ID(NG+NRC)=ID(NP+NRC)                                    ACVI5993
            ELSE                                                        ACVI5994
               NH=NG                                                    ACVI5995
               NG=ID(NG+NLC)                                            ACVI5996
               ID(NH+NLC)=-1                                            ACVI5997
               ID(ID(NH+NRC)+NRC)=NG                                    ACVI5998
               ID(NG+NLC)=NH                                            ACVI5999
               ID(NG+NRC)=ID(NP+NRC)                                    ACVI6000
            ENDIF                                                       ACVI6001
         ELSE                                                           ACVI6002
150         NI=NH                                                       ACVI6003
            NH=NG                                                       ACVI6004
            NG=ID(ID(NG+NLC)+NRC)                                       ACVI6005
            IF(ID(NH+NPR).GE.0) THEN                                    ACVI6006
               NA=NH                                                    ACVI6007
               NF=NI                                                    ACVI6008
            ELSEIF ((ISHFT(ID(NH+NPR),-30).EQ.2)                        ACVI6009
     *          .AND.(ID(ID(NH+NLC)+NPR).GE.0)) THEN                    ACVI6010
               NA=NH                                                    ACVI6011
               NF=NI                                                    ACVI6012
            ENDIF                                                       ACVI6013
            IF (ID(NG+NLC).EQ.-1) THEN                                  ACVI6014
               ID(ID(NH+NLC)+NRC)=NH                                    ACVI6015
            ELSEIF (ID(ID(NG+NLC)+NRC).EQ.NG) THEN                      ACVI6016
               IF (ISHFT(ID(NG+NPR),-30).EQ.2) THEN                     ACVI6017
                  ID(ID(NH+NLC)+NRC)=ID(NG+NLC)                         ACVI6018
                  ID(ID(NG+NLC)+NRC)=NH                                 ACVI6019
               ELSE                                                     ACVI6020
                  NH=NG                                                 ACVI6021
                  NG=ID(NG+NLC)                                         ACVI6022
                  ID(NH+NLC)=-1                                         ACVI6023
               ENDIF                                                    ACVI6024
            ELSE                                                        ACVI6025
               GOTO 150                                                 ACVI6026
            ENDIF                                                       ACVI6027
            ID(ID(ID(NP+NLC)+NRC)+NRC)=NG                               ACVI6028
            ID(NG+NLC)=ID(NP+NLC)                                       ACVI6029
            ID(NG+NRC)=ID(NP+NRC)                                       ACVI6030
         ENDIF                                                          ACVI6031
         IF (NQ.EQ.-1) THEN                                             ACVI6032
            ID(0)=NG                                                    ACVI6033
         ELSEIF (ID(NQ+NLC).EQ.NP) THEN                                 ACVI6034
            ID(NQ+NLC)=NG                                               ACVI6035
         ELSE                                                           ACVI6036
            ID(ID(NQ+NLC)+NRC)=NG                                       ACVI6037
         ENDIF                                                          ACVI6038
         ID(NG+NPR)=IAND(ID(NG+NPR),NF0)+IAND(ID(NP+NPR),N3)            ACVI6039
         IF (NA.EQ.NP) NA=NG                                            ACVI6040
         IF (NF.EQ.NP) NF=NG                                            ACVI6041
         DO 160 I=1,NKEY                                                ACVI6042
160         ID(I)=ID(NG+I)                                              ACVI6043
      ENDIF                                                             ACVI6044
C                                                                       ACVI6045
C Balance the tree                                                      ACVI6046
C                                                                       ACVI6047
500   CONTINUE                                                          ACVI6048
      IF (ID(NA+NLC).EQ.-1) THEN                                        ACVI6049
         ID(NA+NPR)=IAND(ID(NA+NPR),NF0)                                ACVI6050
         GOTO 999                                                       ACVI6051
      ENDIF                                                             ACVI6052
      DO 180 I=1,NKEY                                                   ACVI6053
         IF (ID(I).LT.ID(NA+I)) THEN                                    ACVI6054
            GOTO 190                                                    ACVI6055
         ELSEIF (ID(I).GT.ID(NA+I)) THEN                                ACVI6056
            NO=ID(NA+NLC)                                               ACVI6057
            IF (NO.EQ.-1) THEN                                          ACVI6058
               NP=-1                                                    ACVI6059
            ELSEIF (ID(NO+NRC).EQ.NA) THEN                              ACVI6060
               IF (ISHFT(ID(NA+NPR),-30).EQ.2) THEN                     ACVI6061
                  NP=-1                                                 ACVI6062
               ELSE                                                     ACVI6063
                  NP=NO                                                 ACVI6064
               ENDIF                                                    ACVI6065
            ELSE                                                        ACVI6066
               NP=ID(NO+NRC)                                            ACVI6067
            ENDIF                                                       ACVI6068
            NB=NO                                                       ACVI6069
            IND=N2                                                      ACVI6070
            GOTO 300                                                    ACVI6071
         ENDIF                                                          ACVI6072
180   CONTINUE                                                          ACVI6073
190   NO=ID(NA+NLC)                                                     ACVI6074
      IF (NO.EQ.-1) THEN                                                ACVI6075
         NP=-1                                                          ACVI6076
      ELSEIF (ID(NO+NRC).EQ.NA) THEN                                    ACVI6077
         IF (ISHFT(ID(NA+NPR),-30).EQ.2) THEN                           ACVI6078
            NP=NO                                                       ACVI6079
         ELSE                                                           ACVI6080
            NP=-1                                                       ACVI6081
         ENDIF                                                          ACVI6082
      ELSE                                                              ACVI6083
         NP=NO                                                          ACVI6084
      ENDIF                                                             ACVI6085
      IF (ID(NO+NRC).EQ.NA) THEN                                        ACVI6086
         NB=NO                                                          ACVI6087
      ELSE                                                              ACVI6088
         NB=ID(NO+NRC)                                                  ACVI6089
      ENDIF                                                             ACVI6090
      IND=N3                                                            ACVI6091
300   IF (ID(NA+NPR).GE.0) THEN                                         ACVI6092
         ID(NA+NPR)=IAND(ID(NA+NPR),NF0)+IND                            ACVI6093
         GOTO 600                                                       ACVI6094
      ELSEIF (IAND(ID(NA+NPR),N3).NE.IND) THEN                          ACVI6095
         ID(NA+NPR)=IAND(ID(NA+NPR),NF0)                                ACVI6096
         GOTO 600                                                       ACVI6097
      ENDIF                                                             ACVI6098
C                                                                       ACVI6099
C Rotate                                                                ACVI6100
C                                                                       ACVI6101
      IF (IND.EQ.N2) THEN                                               ACVI6102
         IF (ISHFT(ID(NB+NPR),-30).EQ.3) THEN                           ACVI6103
            IF (ID(ID(NB+NLC)+NRC).EQ.NB) THEN                          ACVI6104
               NC=ID(NB+NLC)                                            ACVI6105
            ELSE                                                        ACVI6106
               NC=ID(ID(NB+NLC)+NRC)                                    ACVI6107
            ENDIF                                                       ACVI6108
            IF (ID(NC+NPR).GE.0) THEN                                   ACVI6109
               IF (ID(NB+NRC).EQ.NA) THEN                               ACVI6110
                  ID(NC+NRC)=ID(NA+NRC)                                 ACVI6111
                  ID(NA+NLC)=-1                                         ACVI6112
                  ID(NA+NRC)=NC                                         ACVI6113
                  ID(NB+NLC)=-1                                         ACVI6114
                  ID(NC+NLC)=NB                                         ACVI6115
                  ID(NA+NPR)=IAND(ID(NA+NPR),NF0)                       ACVI6116
                  ID(NB+NPR)=IAND(ID(NB+NPR),NF0)                       ACVI6117
               ELSE                                                     ACVI6118
                  IDCL=ID(NC+NLC)                                       ACVI6119
                  ID(NA+NLC)=ID(IDCL+NRC)                               ACVI6120
                  ID(ID(NA+NLC)+NRC)=ID(NB+NRC)                         ACVI6121
                  ID(ID(NB+NLC)+NRC)=IDCL                               ACVI6122
                  ID(IDCL+NRC)=NB                                       ACVI6123
                  ID(NC+NLC)=NB                                         ACVI6124
                  ID(NC+NRC)=ID(NA+NRC)                                 ACVI6125
                  ID(NA+NRC)=NC                                         ACVI6126
                  ID(NB+NRC)=NA                                         ACVI6127
                  ID(NA+NPR)=IAND(ID(NA+NPR),NF0)                       ACVI6128
                  ID(NB+NPR)=IAND(ID(NB+NPR),NF0)                       ACVI6129
               ENDIF                                                    ACVI6130
            ELSE                                                        ACVI6131
               IDCL=ID(NC+NLC)                                          ACVI6132
               IF (ID(IDCL+NRC).EQ.NC) THEN                             ACVI6133
                  IF (ISHFT(ID(NC+NPR),-30).EQ.2) THEN                  ACVI6134
                     ID(NA+NLC)=ID(NB+NRC)                              ACVI6135
                     ID(ID(NB+NLC)+NRC)=IDCL                            ACVI6136
                     ID(IDCL+NRC)=NB                                    ACVI6137
                  ELSE                                                  ACVI6138
                     ID(NA+NLC)=IDCL                                    ACVI6139
                     ID(ID(NA+NLC)+NRC)=ID(NB+NRC)                      ACVI6140
                     ID(ID(NB+NLC)+NRC)=NB                              ACVI6141
                  ENDIF                                                 ACVI6142
               ELSE                                                     ACVI6143
                  ID(NA+NLC)=ID(IDCL+NRC)                               ACVI6144
                  ID(ID(NA+NLC)+NRC)=ID(NB+NRC)                         ACVI6145
                  ID(ID(NB+NLC)+NRC)=IDCL                               ACVI6146
                  ID(IDCL+NRC)=NB                                       ACVI6147
               ENDIF                                                    ACVI6148
               ID(NC+NLC)=NB                                            ACVI6149
               ID(NC+NRC)=ID(NA+NRC)                                    ACVI6150
               ID(NA+NRC)=NC                                            ACVI6151
               ID(NB+NRC)=NA                                            ACVI6152
               IF (ISHFT(ID(NC+NPR),-30).EQ.2) THEN                     ACVI6153
                  ID(NA+NPR)=IAND(ID(NA+NPR),NF0)+N3                    ACVI6154
                  ID(NB+NPR)=IAND(ID(NB+NPR),NF0)                       ACVI6155
               ELSE                                                     ACVI6156
                  ID(NA+NPR)=IAND(ID(NA+NPR),NF0)                       ACVI6157
                  ID(NB+NPR)=IAND(ID(NB+NPR),NF0)+N2                    ACVI6158
               ENDIF                                                    ACVI6159
               ID(NC+NPR)=IAND(ID(NC+NPR),NF0)                          ACVI6160
            ENDIF                                                       ACVI6161
            NB=NC                                                       ACVI6162
         ELSE                                                           ACVI6163
            IF (ID(NB+NRC).EQ.NA) THEN                                  ACVI6164
               IF (ID(NB+NPR).GE.0) THEN                                ACVI6165
                  NC=ID(ID(NB+NLC)+NRC)                                 ACVI6166
                  ID(NA+NLC)=NC                                         ACVI6167
                  ID(NC+NRC)=NA                                         ACVI6168
               ELSE                                                     ACVI6169
                  ID(NA+NLC)=-1                                         ACVI6170
               ENDIF                                                    ACVI6171
               ID(NB+NRC)=ID(NA+NRC)                                    ACVI6172
               ID(ID(NB+NLC)+NRC)=NA                                    ACVI6173
               ID(NA+NRC)=NB                                            ACVI6174
            ELSE                                                        ACVI6175
               NC=ID(ID(NB+NLC)+NRC)                                    ACVI6176
               ID(NC+NRC)=ID(NB+NRC)                                    ACVI6177
               ID(NB+NRC)=ID(NA+NRC)                                    ACVI6178
               ID(ID(NB+NLC)+NRC)=NA                                    ACVI6179
               ID(NA+NLC)=NC                                            ACVI6180
               ID(NA+NRC)=NB                                            ACVI6181
            ENDIF                                                       ACVI6182
            IF (ID(NB+NPR).GE.0) THEN                                   ACVI6183
               ID(NA+NPR)=IAND(ID(NA+NPR),NF0)+N2                       ACVI6184
               ID(NB+NPR)=IAND(ID(NB+NPR),NF0)+N3                       ACVI6185
            ELSE                                                        ACVI6186
               ID(NA+NPR)=IAND(ID(NA+NPR),NF0)                          ACVI6187
               ID(NB+NPR)=IAND(ID(NB+NPR),NF0)                          ACVI6188
            ENDIF                                                       ACVI6189
         ENDIF                                                          ACVI6190
      ELSE                                                              ACVI6191
         IF (ISHFT(ID(NB+NPR),-30).EQ.2) THEN                           ACVI6192
            NC=ID(NB+NLC)                                               ACVI6193
            IF (ID(NC+NPR).GE.0) THEN                                   ACVI6194
               IF (ID(NA+NLC).EQ.NB) THEN                               ACVI6195
                  ID(NC+NRC)=ID(NA+NRC)                                 ACVI6196
                  ID(NC+NLC)=NA                                         ACVI6197
                  ID(NA+NLC)=-1                                         ACVI6198
                  ID(NA+NRC)=NB                                         ACVI6199
                  ID(NB+NLC)=-1                                         ACVI6200
                  ID(NB+NRC)=NC                                         ACVI6201
                  ID(NA+NPR)=IAND(ID(NA+NPR),NF0)                       ACVI6202
                  ID(NB+NPR)=IAND(ID(NB+NPR),NF0)                       ACVI6203
               ELSE                                                     ACVI6204
                  IDCL=ID(NC+NLC)                                       ACVI6205
                  ID(ID(NA+NLC)+NRC)=IDCL                               ACVI6206
                  ID(NB+NLC)=ID(IDCL+NRC)                               ACVI6207
                  ID(ID(NB+NLC)+NRC)=ID(NC+NRC)                         ACVI6208
                  ID(IDCL+NRC)=NA                                       ACVI6209
                  ID(NC+NLC)=NA                                         ACVI6210
                  ID(NC+NRC)=ID(NA+NRC)                                 ACVI6211
                  ID(NA+NRC)=NB                                         ACVI6212
                  ID(NB+NRC)=NC                                         ACVI6213
                  ID(NA+NPR)=IAND(ID(NA+NPR),NF0)                       ACVI6214
                  ID(NB+NPR)=IAND(ID(NB+NPR),NF0)                       ACVI6215
               ENDIF                                                    ACVI6216
            ELSE                                                        ACVI6217
               IDCL=ID(NC+NLC)                                          ACVI6218
               IF (ID(IDCL+NRC).EQ.NC) THEN                             ACVI6219
                  IF (ISHFT(ID(NC+NPR),-30).EQ.3) THEN                  ACVI6220
                     ID(ID(NA+NLC)+NRC)=NA                              ACVI6221
                     ID(IDCL+NRC)=ID(NC+NRC)                            ACVI6222
                     ID(NB+NLC)=IDCL                                    ACVI6223
                  ELSE                                                  ACVI6224
                     ID(ID(NA+NLC)+NRC)=IDCL                            ACVI6225
                     ID(IDCL+NRC)=NA                                    ACVI6226
                     ID(NB+NLC)=ID(NC+NRC)                              ACVI6227
                  ENDIF                                                 ACVI6228
               ELSE                                                     ACVI6229
                  ID(ID(NA+NLC)+NRC)=IDCL                               ACVI6230
                  ID(NB+NLC)=ID(IDCL+NRC)                               ACVI6231
                  ID(ID(NB+NLC)+NRC)=ID(NC+NRC)                         ACVI6232
                  ID(IDCL+NRC)=NA                                       ACVI6233
               ENDIF                                                    ACVI6234
               ID(NC+NLC)=NA                                            ACVI6235
               ID(NC+NRC)=ID(NA+NRC)                                    ACVI6236
               ID(NA+NRC)=NB                                            ACVI6237
               ID(NB+NRC)=NC                                            ACVI6238
               IF (ISHFT(ID(NC+NPR),-30).EQ.3) THEN                     ACVI6239
                  ID(NA+NPR)=IAND(ID(NA+NPR),NF0)+N2                    ACVI6240
                  ID(NB+NPR)=IAND(ID(NB+NPR),NF0)                       ACVI6241
               ELSE                                                     ACVI6242
                  ID(NA+NPR)=IAND(ID(NA+NPR),NF0)                       ACVI6243
                  ID(NB+NPR)=IAND(ID(NB+NPR),NF0)+N3                    ACVI6244
               ENDIF                                                    ACVI6245
               ID(NC+NPR)=IAND(ID(NC+NPR),NF0)                          ACVI6246
            ENDIF                                                       ACVI6247
            NB=NC                                                       ACVI6248
         ELSE                                                           ACVI6249
            IF (ID(NA+NLC).EQ.NB) THEN                                  ACVI6250
               IF (ID(NB+NPR).GE.0) THEN                                ACVI6251
                  NC=ID(NB+NLC)                                         ACVI6252
                  ID(NB+NRC)=ID(NA+NRC)                                 ACVI6253
                  ID(NA+NLC)=NC                                         ACVI6254
                  ID(NA+NRC)=ID(NC+NRC)                                 ACVI6255
                  ID(NC+NRC)=NA                                         ACVI6256
               ELSE                                                     ACVI6257
                  ID(NA+NLC)=-1                                         ACVI6258
                  ID(NB+NRC)=ID(NA+NRC)                                 ACVI6259
                  ID(NA+NRC)=ID(NB+NLC)                                 ACVI6260
               ENDIF                                                    ACVI6261
               ID(NB+NLC)=NA                                            ACVI6262
            ELSE                                                        ACVI6263
               NC=ID(NB+NLC)                                            ACVI6264
               ID(NB+NLC)=NA                                            ACVI6265
               ID(NB+NRC)=ID(NA+NRC)                                    ACVI6266
               ID(ID(NA+NLC)+NRC)=NC                                    ACVI6267
               ID(NA+NRC)=ID(NC+NRC)                                    ACVI6268
               ID(NC+NRC)=NA                                            ACVI6269
            ENDIF                                                       ACVI6270
            IF (ID(NB+NPR).GE.0) THEN                                   ACVI6271
               ID(NA+NPR)=IAND(ID(NA+NPR),NF0)+N3                       ACVI6272
               ID(NB+NPR)=IAND(ID(NB+NPR),NF0)+N2                       ACVI6273
            ELSE                                                        ACVI6274
               ID(NA+NPR)=IAND(ID(NA+NPR),NF0)                          ACVI6275
               ID(NB+NPR)=IAND(ID(NB+NPR),NF0)                          ACVI6276
            ENDIF                                                       ACVI6277
         ENDIF                                                          ACVI6278
      ENDIF                                                             ACVI6279
      IF (NF.EQ.-1) THEN                                                ACVI6280
         ID(0)=NB                                                       ACVI6281
      ELSEIF (NA.EQ.ID(NF+NLC)) THEN                                    ACVI6282
         ID(NF+NLC)=NB                                                  ACVI6283
      ELSE                                                              ACVI6284
         ID(ID(NF+NLC)+NRC)=NB                                          ACVI6285
      ENDIF                                                             ACVI6286
600   NF=NA                                                             ACVI6287
      NA=NP                                                             ACVI6288
      IF (NP.NE.-1) GOTO 500                                            ACVI6289
999   CONTINUE                                                          ACVI6290
C                                                                       ACVI6291
C Reconstruct the priority array                                        ACVI6292
C                                                                       ACVI6293
      NA=0                                                              ACVI6294
      NL=ID(-1)-NRC                                                     ACVI6295
      NLPR=IAND(ID(NL+NPR),NF0)                                         ACVI6296
1000  NB=NA+NA+NRC                                                      ACVI6297
      IF (NB.LT.NL)THEN                                                 ACVI6298
         NC=NB+NRC                                                      ACVI6299
         NBPR=IAND(ID(NB+NPR),NF0)                                      ACVI6300
         NCPR=IAND(ID(NC+NPR),NF0)                                      ACVI6301
         IF ((NBPR.LE.NCPR).OR.(NC.GT.NL)) THEN                         ACVI6302
            IF (NBPR.LT.NLPR) THEN                                      ACVI6303
               NBRC=ID(NB+NRC)                                          ACVI6304
               NBLC=ID(NB+NLC)                                          ACVI6305
               IF (NBRC.EQ.-1) THEN                                     ACVI6306
                  ID(0)=NA                                              ACVI6307
               ELSE                                                     ACVI6308
                  IF (ID(NBRC+NLC).EQ.NB) THEN                          ACVI6309
                     ID(NBRC+NLC)=NA                                    ACVI6310
                  ELSE                                                  ACVI6311
                     IF (ID(ID(NBRC+NRC)+NLC).EQ.NB) THEN               ACVI6312
                        ID(ID(NBRC+NRC)+NLC)=NA                         ACVI6313
                     ELSE                                               ACVI6314
                        ID(ID(NBRC+NLC)+NRC)=NA                         ACVI6315
                     ENDIF                                              ACVI6316
                  ENDIF                                                 ACVI6317
               ENDIF                                                    ACVI6318
               IF (NBLC.NE.-1) THEN                                     ACVI6319
                  IF (ID(NBLC+NRC).EQ.NB) THEN                          ACVI6320
                     ID(NBLC+NRC)=NA                                    ACVI6321
                  ELSE                                                  ACVI6322
                     ID(ID(NBLC+NRC)+NRC)=NA                            ACVI6323
                  ENDIF                                                 ACVI6324
               ENDIF                                                    ACVI6325
               DO 1400 I=1,NRC                                          ACVI6326
                  ID(NA+I)=ID(NB+I)                                     ACVI6327
1400           CONTINUE                                                 ACVI6328
               NA=NB                                                    ACVI6329
               GOTO 1000                                                ACVI6330
            ENDIF                                                       ACVI6331
         ELSE                                                           ACVI6332
            IF (NCPR.LT.NLPR) THEN                                      ACVI6333
               NCRC=ID(NC+NRC)                                          ACVI6334
               NCLC=ID(NC+NLC)                                          ACVI6335
               IF (NCRC.EQ.-1) THEN                                     ACVI6336
                  ID(0)=NA                                              ACVI6337
               ELSE                                                     ACVI6338
                  IF (ID(NCRC+NLC).EQ.NC) THEN                          ACVI6339
                     ID(NCRC+NLC)=NA                                    ACVI6340
                  ELSE                                                  ACVI6341
                     IF (ID(ID(NCRC+NRC)+NLC).EQ.NC) THEN               ACVI6342
                        ID(ID(NCRC+NRC)+NLC)=NA                         ACVI6343
                     ELSE                                               ACVI6344
                        ID(ID(NCRC+NLC)+NRC)=NA                         ACVI6345
                     ENDIF                                              ACVI6346
                  ENDIF                                                 ACVI6347
               ENDIF                                                    ACVI6348
               IF (NCLC.NE.-1) THEN                                     ACVI6349
                  IF (ID(NCLC+NRC).EQ.NC) THEN                          ACVI6350
                     ID(NCLC+NRC)=NA                                    ACVI6351
                  ELSE                                                  ACVI6352
                     ID(ID(NCLC+NRC)+NRC)=NA                            ACVI6353
                  ENDIF                                                 ACVI6354
               ENDIF                                                    ACVI6355
               DO 2400 I=1,NRC                                          ACVI6356
                  ID(NA+I)=ID(NC+I)                                     ACVI6357
2400           CONTINUE                                                 ACVI6358
               NA=NC                                                    ACVI6359
               GOTO 1000                                                ACVI6360
            ENDIF                                                       ACVI6361
         ENDIF                                                          ACVI6362
      ENDIF                                                             ACVI6363
      NLRC=ID(NL+NRC)                                                   ACVI6364
      NLLC=ID(NL+NLC)                                                   ACVI6365
      IF (NLRC.EQ.-1) THEN                                              ACVI6366
         ID(0)=NA                                                       ACVI6367
      ELSE                                                              ACVI6368
         IF (ID(NLRC+NLC).EQ.NL) THEN                                   ACVI6369
            ID(NLRC+NLC)=NA                                             ACVI6370
         ELSE                                                           ACVI6371
            IF (ID(ID(NLRC+NRC)+NLC).EQ.NL) THEN                        ACVI6372
               ID(ID(NLRC+NRC)+NLC)=NA                                  ACVI6373
            ELSE                                                        ACVI6374
               ID(ID(NLRC+NLC)+NRC)=NA                                  ACVI6375
            ENDIF                                                       ACVI6376
         ENDIF                                                          ACVI6377
      ENDIF                                                             ACVI6378
      IF (NLLC.NE.-1) THEN                                              ACVI6379
         IF (ID(NLLC+NRC).EQ.NL) THEN                                   ACVI6380
            ID(NLLC+NRC)=NA                                             ACVI6381
         ELSE                                                           ACVI6382
            ID(ID(NLLC+NRC)+NRC)=NA                                     ACVI6383
         ENDIF                                                          ACVI6384
      ENDIF                                                             ACVI6385
      DO 7000 I=1,NRC                                                   ACVI6386
         ID(NA+I)=ID(NL+I)                                              ACVI6387
         ID(NL+I)=0                                                     ACVI6388
7000  CONTINUE                                                          ACVI6389
      ID(-1)=NL                                                         ACVI6390
9999  RETURN                                                            ACVI6391
      END                                                               ACVI6392
C ----------------------------------------------------------------------ACVI6393
C                                                                       ACVI6394
C                           **************                              ACVI6395
C                           ***  TOUT  ***                              ACVI6396
C                           **************                              ACVI6397
C                                                                       ACVI6398
C The subroutine TOUT traverses the weighted search tree in order to    ACVI6399
C specific nodes and sends the outputs to an output devise designated   ACVI6400
C by NFILE. The traversal is done in ascending order if NAD=1 and in    ACVI6401
C descending order otherwise. A description of the arguments follows:   ACVI6402
C                                                                       ACVI6403
C   NFILE = File number assigned to output device where the results are ACVI6404
C           to be written. In the special case NFILE = 0, TOUT does not ACVI6405
C           output anything, however after upon returning ID(-5) points ACVI6406
C           to the current node position.                               ACVI6407
C     NAD = 1 if the traversal is to be in ascending order, otherwise   ACVI6408
C           it is done in descending order.                             ACVI6409
C     MIN = Number of the first node that is to be retrieved.           ACVI6410
C     MAX = Number of the last node that is to be retrieved.            ACVI6411
C   NSTEP = Step size; nodes that are integer multiples of NSTEP beyond ACVI6412
C           MIN up to MAX will be retrieved.                            ACVI6413
C      ID = Name of the tree array that is to be traversed.             ACVI6414
C                                                                       ACVI6415
C An internal array called STACK is used to store information about the ACVI6416
C path traversed. It is dimensioned 32 since this is the maximum height ACVI6417
C the tree can have.  This limit is set by the fact that 2**32-1 is the ACVI6418
C largest integer pointer that can be used on a 32 bit machine.         ACVI6419
C                                                                       ACVI6420
C ----------------------------------------------------------------------ACVI6421
C                                                                       ACVI6422
      SUBROUTINE TOUT(NFILE,NAD,MIN,MAX,NSTEP,ID)                       ACVI6423
C                                                                       ACVI6424
      INTEGER ID(-10:*),STACK(32)                                       ACVI6425
      NSUM=ID(-3)                                                       ACVI6426
      NLC=NSUM+2                                                        ACVI6427
      NRC=NLC+1                                                         ACVI6428
      I=0                                                               ACVI6429
      M=0                                                               ACVI6430
      NUM=MIN                                                           ACVI6431
      NODE=ID(0)                                                        ACVI6432
      IF (NAD.EQ.1) THEN                                                ACVI6433
200      IF (NODE.NE.-1) THEN                                           ACVI6434
            I=I+1                                                       ACVI6435
            STACK(I)=NODE                                               ACVI6436
            NO=ID(NODE+NLC)                                             ACVI6437
            IF (NO.EQ.-1) THEN                                          ACVI6438
               NODE=-1                                                  ACVI6439
            ELSEIF (ID(NO+NRC).EQ.NODE) THEN                            ACVI6440
               IF (ISHFT(ID(NODE+NPR),-30).EQ.2) THEN                   ACVI6441
                  NODE=NO                                               ACVI6442
               ELSE                                                     ACVI6443
                  NODE=-1                                               ACVI6444
               ENDIF                                                    ACVI6445
            ELSE                                                        ACVI6446
               NODE=NO                                                  ACVI6447
            ENDIF                                                       ACVI6448
            GOTO 200                                                    ACVI6449
300         M=M+1                                                       ACVI6450
            ID(-5)=NODE                                                 ACVI6451
            IF ((M.EQ.NUM).OR.(M.EQ.MAX)) THEN                          ACVI6452
               IF (NFILE.GT.0) WRITE(NFILE,1000)                        ACVI6453
     *           M,(ID(NODE+II),II=1,NSUM)                              ACVI6454
               IF (M.EQ.MAX) GOTO 999                                   ACVI6455
               NUM=NUM+NSTEP                                            ACVI6456
            ENDIF                                                       ACVI6457
            NO=ID(NODE+NLC)                                             ACVI6458
            IF (NO.EQ.-1) THEN                                          ACVI6459
               NODE=-1                                                  ACVI6460
            ELSEIF (ID(NO+NRC).EQ.NODE) THEN                            ACVI6461
               IF (ISHFT(ID(NODE+NPR),-30).EQ.2) THEN                   ACVI6462
                  NODE=-1                                               ACVI6463
               ELSE                                                     ACVI6464
                  NODE=NO                                               ACVI6465
               ENDIF                                                    ACVI6466
            ELSE                                                        ACVI6467
               NODE=ID(NO+NRC)                                          ACVI6468
            ENDIF                                                       ACVI6469
            GOTO 200                                                    ACVI6470
         ENDIF                                                          ACVI6471
         IF (I.NE.0) THEN                                               ACVI6472
            NODE=STACK(I)                                               ACVI6473
            I=I-1                                                       ACVI6474
            GOTO 300                                                    ACVI6475
         ENDIF                                                          ACVI6476
      ELSE                                                              ACVI6477
400      IF (NODE.NE.-1) THEN                                           ACVI6478
            I=I+1                                                       ACVI6479
            STACK(I)=NODE                                               ACVI6480
            NO=ID(NODE+NLC)                                             ACVI6481
            IF (NO.EQ.-1) THEN                                          ACVI6482
               NODE=-1                                                  ACVI6483
            ELSEIF (ID(NO+NRC).EQ.NODE) THEN                            ACVI6484
               IF (ISHFT(ID(NODE+NPR),-30).EQ.2) THEN                   ACVI6485
                  NODE=-1                                               ACVI6486
               ELSE                                                     ACVI6487
                  NODE=NO                                               ACVI6488
               ENDIF                                                    ACVI6489
            ELSE                                                        ACVI6490
               NODE=ID(NO+NRC)                                          ACVI6491
            ENDIF                                                       ACVI6492
            GOTO 400                                                    ACVI6493
500         M=M+1                                                       ACVI6494
            ID(-5)=NODE                                                 ACVI6495
            IF ((M.EQ.NUM).OR.(M.EQ.MAX)) THEN                          ACVI6496
               IF (NFILE.GT.0) WRITE(NFILE,1000)                        ACVI6497
     *           M,(ID(NODE+II),II=1,NSUM)                              ACVI6498
               IF (M.EQ.MAX) GOTO 999                                   ACVI6499
               NUM=NUM+NSTEP                                            ACVI6500
            ENDIF                                                       ACVI6501
            NO=ID(NODE+NLC)                                             ACVI6502
            IF (NO.EQ.-1) THEN                                          ACVI6503
               NODE=-1                                                  ACVI6504
            ELSEIF (ID(NO+NRC).EQ.NODE) THEN                            ACVI6505
               IF (ISHFT(ID(NODE+NPR),-30).EQ.2) THEN                   ACVI6506
                  NODE=NO                                               ACVI6507
               ELSE                                                     ACVI6508
                  NODE=-1                                               ACVI6509
               ENDIF                                                    ACVI6510
            ELSE                                                        ACVI6511
               NODE=NO                                                  ACVI6512
            ENDIF                                                       ACVI6513
            GOTO 400                                                    ACVI6514
         ENDIF                                                          ACVI6515
         IF (I.NE.0) THEN                                               ACVI6516
            NODE=STACK(I)                                               ACVI6517
            I=I-1                                                       ACVI6518
            GOTO 500                                                    ACVI6519
         ENDIF                                                          ACVI6520
      ENDIF                                                             ACVI6521
      WRITE(NFILE,*)                                                    ACVI6522
      WRITE(NFILE,*)'  WARNING:  TOUT SEARCH HAS GONE OUT OF BOUNDS!'   ACVI6523
999   RETURN                                                            ACVI6524
1000  FORMAT(3X,I6,':',2X,14I5/(12X,14I5))                              ACVI6525
      END                                                               ACVI6526
C ----------------------------------------------------------------------ACVI6527
C                                                                       ACVI6528
C                            ************                               ACVI6529
C                            *** TMRG ***                               ACVI6530
C                            ************                               ACVI6531
C                                                                       ACVI6532
C The subroutine TMRG allows the user to merge two trees that have the  ACVI6533
C same structure. Elements of the first tree are stored in the second.  ACVI6534
C If the merge results in an overflow condition, notification is given  ACVI6535
C and the lowest priority item is eliminated in favor of the incoming   ACVI6536
C element with higher priority. If the trees have a different structure,ACVI6537
C the merge is automatically aborted (alternate return) after giving an ACVI6538
C appropriate warning message. If the element already exists in the     ACVI6539
C second tree, the priority is set to the maximum of the two.           ACVI6540
C                                                                       ACVI6541
C   RETURN 1 --> Unsuccessful match (the trees have different           ACVI6542
C                structures) so merge is aborted.                       ACVI6543
C                                                                       ACVI6544
C ----------------------------------------------------------------------ACVI6545
C                                                                       ACVI6546
      SUBROUTINE TMRG(ID1,BUFF1,LLBF1,ID2,BUFF2,LLBF2,*)                ACVI6547
C                                                                       ACVI6548
      INTEGER ID1(-10:*),LLBF1(-2:*),ID2(-10:*),LLBF2(-2:*)             ACVI6549
      REAL    BUFF1(*),BUFF2(*)                       ! *** storage typeACVI6550
C     REAL*8  BUFF1(*),BUFF2(*)                       ! *** storage typeACVI6551
      PARAMETER (MAXDAT=1000)                                           ACVI6552
      INTEGER NDAT(2)                                                   ACVI6553
      REAL    BULOAD(MAXDAT)                          ! *** temp storageACVI6554
C     REAL*8  BULOAD(MAXDAT)                          ! *** temp storageACVI6555
C     DATA    NF0/Z3FFFFFFF/                                            ACVI6556
      DATA    NF0/1073741823/                                           ACVI6557
C                                                                       ACVI6558
C Check the structure of the trees for compatibility                    ACVI6559
C                                                                       ACVI6560
      IF (ID1(-4).NE.ID2(-4).OR.ID1(-3).NE.ID2(-3)) THEN                ACVI6561
         WRITE(6,*) ' TREE 1 AND 2 HAVE DIFFERENT STRUCTURES'           ACVI6562
         RETURN 1                                                       ACVI6563
      ENDIF                                                             ACVI6564
      IF (ID1(-10)+ID2(-10).GE.ID2(-9).OR.LLBF1(-1).GE.LLBF2(0)) THEN   ACVI6565
         WRITE(6,*) ' SIZE OF THE INCOMING NODE IS LARGER THAN MAXIMUM,'ACVI6566
         WRITE(6,*) ' ONE OR MORE LOW PRIORITY NODES WILL BE DELETED.'  ACVI6567
      ENDIF                                                             ACVI6568
C                                                                       ACVI6569
C List out information from tree 1                                      ACVI6570
C                                                                       ACVI6571
      DO INODE1=1,ID1(-10)                                              ACVI6572
         CALL TOUT(0,0,INODE1,INODE1,1,ID1)                             ACVI6573
C                                                                       ACVI6574
C ...find the key and data indices in first tree                        ACVI6575
C                                                                       ACVI6576
         IPNODE=ID1(-5)                                                 ACVI6577
         IPKEY=IPNODE+1                                                 ACVI6578
         IPRIOR=ID1(IPNODE+ID1(-2))                                     ACVI6579
C                                                                       ACVI6580
C Check whether or not the incoming node is new                         ACVI6581
C                                                                       ACVI6582
         CALL TCHK(ID1(IPKEY),ID2,*100)                                 ACVI6583
C                                                                       ACVI6584
C ...find the data information from first tree                          ACVI6585
C                                                                       ACVI6586
         IPDATA=IPNODE+ID1(-4)+1                                        ACVI6587
         IFIND=ID1(IPDATA)                                              ACVI6588
         NOSIZE=ID1(IPDATA+1)                                           ACVI6589
C                                                                       ACVI6590
C ...abort if the size is too big                                       ACVI6591
C                                                                       ACVI6592
         IF (NOSIZE.GT.MAXDAT) THEN                                     ACVI6593
            WRITE(6,*) ' SIZE OF TEMPORARY BUFFER BULOAD IS TOO SMALL.' ACVI6594
            WRITE(6,*) ' INCREASE MAXDAT IN SUBROUTINE TMRG AND'        ACVI6595
            WRITE(6,*) ' CONTINUE.'                                     ACVI6596
            STOP                                                        ACVI6597
         ENDIF                                                          ACVI6598
C                                                                       ACVI6599
         DO IKK=1,NOSIZE                                                ACVI6600
            BULOAD(IKK)=BUFF1(IFIND)                                    ACVI6601
            IFIND=LLBF1(IFIND)                                          ACVI6602
         ENDDO                                                          ACVI6603
C                                                                       ACVI6604
C Store key, data, and information in second tree                       ACVI6605
C                                                                       ACVI6606
         CALL TADD(ID1(IPKEY),NDAT,BULOAD,NOSIZE,IPRIOR,ID2,BUFF2,LLBF2)ACVI6607
         GOTO 200                                                       ACVI6608
C                                                                       ACVI6609
C Update the priority if the incoming node is old                       ACVI6610
C                                                                       ACVI6611
100      IFIND=ID2(-5)+ID2(-2)                                          ACVI6612
         ID2(IFIND)=MAX0(ID2(IFIND),IPRIOR)                             ACVI6613
C                                                                       ACVI6614
200      CONTINUE                                                       ACVI6615
      ENDDO                                                             ACVI6616
      RETURN                                                            ACVI6617
      END                                                               ACVI6618
C ----------------------------------------------------------------------ACVI6619
C                                                                       ACVI6620
C                       ***********************                         ACVI6621
C                       ***   SU3 PACKAGE   ***                         ACVI6622
C                       ***********************                         ACVI6623
C                                                                       ACVI6624
C ----------------------------------------------------------------------ACVI6625
      SUBROUTINE BLOCKS                                                 ACVI6626
C     ------------------------------------------------------------------ACVI6627
C     BINOMIAL COEFFICIENTS AND FACTORIALS ***** SEE COMMENT BELOW *****ACVI6628
C     ------------------------------------------------------------------ACVI6629
C     UPDATE/MOD: (MTS,06-76)  H.SATO             EXPANDED RANGE        ACVI6630
C                 (LSU,11-78)  J.P.DRAAYER        LOG FACTORIALS        ACVI6631
C                 (LSU,05-80)  J.P.DRAAYER        LOG BINOMIALS         ACVI6632
C                 (LSU,08-81)  J.P.DRAAYER        EXPANDED RANGE        ACVI6633
C                 (LSU,01-83)  J.P.DRAAYER        EXTENDED PRECISION    ACVI6634
C                 (LSU,03-83)  J.P.DRAAYER        D,Q MIX & LOG CUTS    ACVI6635
C                 (LSU,11-84)  J.P.DRAAYER        MODIFIED FACTORIALS   ACVI6636
C                 (LSU,01-88)  J.P.DRAAYER        DLOGF RANGE/INDEX     ACVI6637
C                 (LSU,10-89)  J.P.DRAAYER        BINOMIAL INVERSES     ACVI6638
C                 (LSU,11-89)  J.P.DRAAYER        POWERS OF TWO ARRAY   ACVI6639
C                                                                       ACVI6640
C     BKDB--BINOMIAL (BINO) COEFFICIENTS (EXPANDED 6-76,8-81,1-83)      ACVI6641
C           SCALE: BINO(I,J)=DBINO(I*(I+1)/2+J+1)                       ACVI6642
C           RANGE: BINO(0,0)=DBINO(1) TO BINO(128,128)=DBINO(8385)      ACVI6643
C           ADDED: 2**I = DTWOS(I) WHERE I=-128,128                     ACVI6644
C     BKQB--BINOMIAL (BINO) COEFFICIENTS (EXPANDED 6-76,8-81,1-83)      ACVI6645
C           SCALE: BINO(I,J)=QBINO(I*(I+1)/2+J+1)                       ACVI6646
C           RANGE: BINO(0,0)=QBINO(1) TO BINO(192,192)=QBINO(18721)     ACVI6647
C           ADDED: 2**I = QTWOS(I) WHERE I=-192,192                     ACVI6648
C     BKDF--LOG FACTORIALS (FACT) (INSERTED 11-78, MODIFIED 01-88)      ACVI6649
C           SCALE: LNFACT(I)=DLOGF(2*I)                                 ACVI6650
C           RANGE: LNFACT(0)=DLOGF(0) TO LNFACT(1000)=DLOGF(2000)       ACVI6651
C                                                                       ACVI6652
C          ********************************************************     ACVI6653
C          **  BLOCKS INPUT MUST BE PREGENERATED USING SU3GENBK  **     ACVI6654
C          ********************************************************     ACVI6655
C                                                                       ACVI6656
C     ------------------------------------------------------------------ACVI6657
      IMPLICIT REAL*8(D),REAL*16(Q)                                     ACVI6658
      COMMON/BKDB/DBINO(8385),DBINV(8385),DTWOS(-128:128)               ACVI6659
      COMMON/BKQB/QBINO(18721),QBINV(18721),QTWOS(-192:192)             ACVI6660
      COMMON/BKDF/DLOGF(0:2000)                                         ACVI6661
      READ(4)DLOGF                                                      ACVI6662
      READ(4)DBINO                                                      ACVI6663
      READ(4)DBINV                                                      ACVI6664
      READ(4)QBINO                                                      ACVI6665
      READ(4)QBINV                                                      ACVI6666
      READ(4)DTWOS                                                      ACVI6667
      READ(4)QTWOS                                                      ACVI6668
      RETURN                                                            ACVI6669
      END                                                               ACVI6670
      FUNCTION MULTU3(LX1,MX1,LX2,MX2,LX3,MX3)                          ACVI6671
C     ------------------------------------------------------------------ACVI6672
C     MULTIPLICITY IN U3 COUPLING (SEE MULTTEST FOR VARIOUS VERSIONS)   ACVI6673
C     ... FASTER THAN THE DRAAYER ORIGINAL AND THE MILLENER KAS FUNCTIONACVI6674
C     ------------------------------------------------------------------ACVI6675
C     UPDATE/MOD: (LSU,11-78)  J.P.DRAAYER        ORIGINAL VERSION      ACVI6676
C                 (BNL,06-87)  J.MILLENER         MILLENER VERSION      ACVI6677
C                 (LSU,10-89)  J.P.DRAAYER        PROFILER OPTIMIZED    ACVI6678
C     ------------------------------------------------------------------ACVI6679
      MULTU3=0                                                          ACVI6680
      NX=LX1+LX2-LX3-MX1-MX2+MX3                                        ACVI6681
      MX=NX/3                                                           ACVI6682
      IF(3*MX.NE.NX)RETURN                                              ACVI6683
      IF(MX.GE.0)THEN                                                   ACVI6684
         L1=LX1                                                         ACVI6685
         L2=LX2                                                         ACVI6686
         L3=LX3                                                         ACVI6687
         M1=MX1                                                         ACVI6688
         M2=MX2                                                         ACVI6689
         M3=MX3                                                         ACVI6690
      ELSE                                                              ACVI6691
         L1=MX1                                                         ACVI6692
         L2=MX2                                                         ACVI6693
         L3=MX3                                                         ACVI6694
         M1=LX1                                                         ACVI6695
         M2=LX2                                                         ACVI6696
         M3=LX3                                                         ACVI6697
         MX=-MX                                                         ACVI6698
      ENDIF                                                             ACVI6699
      NX=MX+M1+M2-M3                                                    ACVI6700
      MU=MIN0(L1-MX,M2)                                                 ACVI6701
      IF(MU.LT.0)RETURN                                                 ACVI6702
      NU=MIN0(L2-MX,M1)                                                 ACVI6703
      IF(NU.LT.0)RETURN                                                 ACVI6704
      MULTU3=MAX0(MIN0(NX,NU)-MAX0(NX-MU,0)+1,0)                        ACVI6705
      RETURN                                                            ACVI6706
      END                                                               ACVI6707
      SUBROUTINE U3MULT(LX1,MX1,LX2,MX2,LX3,MX3,MULTU3,*)               ACVI6708
C     ------------------------------------------------------------------ACVI6709
C     MULTIPLICITY IN U3 COUPLING (SEE MULTTEST FOR VARIOUS VERSIONS)   ACVI6710
C     ... MULTU3 IS FASTER THAN THE DRAAYER ORIGINAL AND MILLENER KAS   ACVI6711
C     ... SUBROUTINE FORM EXECUTES FASTEST ... RETURN 1 BRANCH OPTION   ACVI6712
C     ------------------------------------------------------------------ACVI6713
C     UPDATE/MOD: (LSU,11-78)  J.P.DRAAYER     ORIGINAL FUNCTION FORM   ACVI6714
C                 (BNL,06-87)  J.MILLENER      MILLENER FUNCTION ... KASACVI6715
C                 (LSU,10-89)  J.P.DRAAYER     PROFILER OPTIMIZED FORM  ACVI6716
C                 (LSU,11-89)  J.P.DRAAYER     SUBROUTINE FORM FOR SPEEDACVI6717
C     ------------------------------------------------------------------ACVI6718
      MULTU3=0                                                          ACVI6719
      NX=LX1+LX2-LX3-MX1-MX2+MX3                                        ACVI6720
      MX=NX/3                                                           ACVI6721
      IF(3*MX.NE.NX)RETURN 1                                            ACVI6722
      IF(MX.GE.0)THEN                                                   ACVI6723
         L1=LX1                                                         ACVI6724
         L2=LX2                                                         ACVI6725
         L3=LX3                                                         ACVI6726
         M1=MX1                                                         ACVI6727
         M2=MX2                                                         ACVI6728
         M3=MX3                                                         ACVI6729
      ELSE                                                              ACVI6730
         L1=MX1                                                         ACVI6731
         L2=MX2                                                         ACVI6732
         L3=MX3                                                         ACVI6733
         M1=LX1                                                         ACVI6734
         M2=LX2                                                         ACVI6735
         M3=LX3                                                         ACVI6736
         MX=-MX                                                         ACVI6737
      ENDIF                                                             ACVI6738
      NX=MX+M1+M2-M3                                                    ACVI6739
      MU=MIN0(L1-MX,M2)                                                 ACVI6740
      IF(MU.LT.0)RETURN 1                                               ACVI6741
      NU=MIN0(L2-MX,M1)                                                 ACVI6742
      IF(NU.LT.0)RETURN 1                                               ACVI6743
      MULTU3=MAX0(MIN0(NX,NU)-MAX0(NX-MU,0)+1,0)                        ACVI6744
      IF(MULTU3.NE.0)RETURN                                             ACVI6745
      RETURN 1                                                          ACVI6746
      END                                                               ACVI6747
      SUBROUTINE XEWU3S(INC,LAM1,MU1,LAM2,MU2,NEC,NNC,KR0A,KR0B,DEWU3P, ACVI6748
     1J1TA,IAB,ICD,INDMAX,DEWU3,KR0MAX)                                 ACVI6749
C     ------------------------------------------------------------------ACVI6750
C     SUPPORT ROUTINE FOR XEWU3S (X PREFIX FOR 6-81 VERSION)            ACVI6751
C     ------------------------------------------------------------------ACVI6752
C     UPDATE/MOD: (LSU,06-81)  J.P.DRAAYER        INDEXING DEWU3        ACVI6753
C                                                                       ACVI6754
C     PARAMETERS--                                                      ACVI6755
C       INC=0:<(LAM1,MU1)????;(LAM2,MU2)????::KR0(LAM3,MU3)HIGH>        ACVI6756
C             FROM <(LAM1,MU1)????;(LAM2,MU2)HIGH::KR0(LAM3,MU3)HIGH>   ACVI6757
C       INC=1:<(LAM1,MU1)????;(LAM2,MU2)????::KR0(LAM3,MU3)HIGH>        ACVI6758
C             FROM <(LAM1,MU1)HIGH;(LAM2,MU2)????::KR0(LAM3,MU3)HIGH>   ACVI6759
C     ------------------------------------------------------------------ACVI6760
      IMPLICIT REAL*8(D)                                                ACVI6761
      DIMENSION DEWU3(1),J1TA(1),DEWU3P(1),IAB(1),ICD(1)                ACVI6762
      INDEX(J1TD,LAM1,J1T,J2TD,LAM2,J2T)=1+J2TD*(J2TD+1)*(3*J1TD+J2TD+5)ACVI6763
     1/6+(J1TD+1)*(LAM2+J2TD-J2T)/2+(LAM1+J1TD-J1T)/2                   ACVI6764
      INDQ=-KR0MAX                                                      ACVI6765
      IF(INC.EQ.1)INDQ=INDQ+(INDMAX-NNC)*KR0MAX                         ACVI6766
      INDPQ=-KR0MAX                                                     ACVI6767
      DO 10 IIQ=1,NNC                                                   ACVI6768
      INDQ=INDQ+KR0MAX                                                  ACVI6769
      INDPQ=INDPQ+KR0MAX                                                ACVI6770
      DO 10 KR0=KR0A,KR0B                                               ACVI6771
      KI=KR0+INDQ                                                       ACVI6772
      KIP=KR0+INDPQ                                                     ACVI6773
   10 DEWU3(KI)=DEWU3P(KIP)                                             ACVI6774
      IF(NEC.EQ.0)RETURN                                                ACVI6775
      IF(INC.EQ.1)GO TO 15                                              ACVI6776
      L1=LAM1                                                           ACVI6777
      M1=MU1                                                            ACVI6778
      L2=LAM2                                                           ACVI6779
      M2=MU2                                                            ACVI6780
      GO TO 20                                                          ACVI6781
   15 L1=LAM2                                                           ACVI6782
      M1=MU2                                                            ACVI6783
      L2=LAM1                                                           ACVI6784
      M2=MU1                                                            ACVI6785
   20 LL1=L1+1                                                          ACVI6786
      MM1=M1+1                                                          ACVI6787
      LL2=L2+1                                                          ACVI6788
      MM2=M2+1                                                          ACVI6789
      LM1=LL1+MM1                                                       ACVI6790
      LM2=LL2+MM2                                                       ACVI6791
      DO 70 J2TD=1,NEC                                                  ACVI6792
      J1TD=NEC-J2TD                                                     ACVI6793
      J2D=J2TD-1                                                        ACVI6794
      J1D=J1TD+1                                                        ACVI6795
      IIQ2A=J2TD-M2                                                     ACVI6796
      IF(IIQ2A.LT.0)IIQ2A=0                                             ACVI6797
      IIQ2A=IIQ2A+1                                                     ACVI6798
      IIQ2B=J2TD                                                        ACVI6799
      IF(L2.LT.IIQ2B)IIQ2B=L2                                           ACVI6800
      IIQ2B=IIQ2B+1                                                     ACVI6801
      IIQ1A=J1TD-M1                                                     ACVI6802
      IF(IIQ1A.LT.0)IIQ1A=0                                             ACVI6803
      IIQ1A=IIQ1A+1                                                     ACVI6804
      IIQ1B=J1TD                                                        ACVI6805
      IF(L1.LT.IIQ1B)IIQ1B=L1                                           ACVI6806
      IIQ1B=IIQ1B+1                                                     ACVI6807
      DO 70 IIQ2=IIQ2A,IIQ2B                                            ACVI6808
      IQ2=IIQ2-1                                                        ACVI6809
      IP2=J2TD-IQ2                                                      ACVI6810
      J2T=L2+IP2-IQ2                                                    ACVI6811
      JJ2T=J2T+1                                                        ACVI6812
      IQ=-1                                                             ACVI6813
      IP=-1                                                             ACVI6814
      IF(IP2.EQ.0)GO TO 25                                              ACVI6815
      IQ2P=IQ2                                                          ACVI6816
      IP2P=IP2-1                                                        ACVI6817
      IQ2D=0                                                            ACVI6818
      IP2D=1                                                            ACVI6819
      IF(INC.EQ.1)IP=1                                                  ACVI6820
      NM=IP2*(M2-IP2P)*(LL2+IP2)                                        ACVI6821
      DN=DFLOAT(J2T)                                                    ACVI6822
      GO TO 30                                                          ACVI6823
   25 IQ2P=IQ2-1                                                        ACVI6824
      IP2P=IP2                                                          ACVI6825
      IQ2D=1                                                            ACVI6826
      IP2D=0                                                            ACVI6827
      IF(INC.EQ.0)IQ=1                                                  ACVI6828
      NM=IQ2*(L2-IQ2P)*(LM2-IQ2)                                        ACVI6829
      DN=DFLOAT(J2T+2)                                                  ACVI6830
   30 J2TP=L2+IP2P-IQ2P                                                 ACVI6831
      JTA=J2TD-IQ2P                                                     ACVI6832
      JTB=NNC-JTA                                                       ACVI6833
      NQD=NEC-IQ2P                                                      ACVI6834
      DO 70 IIQ1=IIQ1A,IIQ1B                                            ACVI6835
      IQ1=IIQ1-1                                                        ACVI6836
      IP1=J1TD-IQ1                                                      ACVI6837
      J1T=L1+IP1-IQ1                                                    ACVI6838
      IF(INC.EQ.0)IND=INDEX(J1TD,L1,J1T,J2TD,L2,J2T)                    ACVI6839
      IF(INC.EQ.1)IND=INDEX(J2TD,L2,J2T,J1TD,L1,J1T)                    ACVI6840
      IF(J1TA(IND).LT.0)GO TO 70                                        ACVI6841
      IF(IP1.EQ.M1)GO TO 50                                             ACVI6842
      IQ1P=IQ1                                                          ACVI6843
      IP1P=IP1+1                                                        ACVI6844
      J1TP=J1T+1                                                        ACVI6845
      IF(INC.EQ.0)GO TO 35                                              ACVI6846
      INDP=INDEX(J2D,L2,J2TP,J1D,L1,J1TP)                               ACVI6847
      IF(J1TA(INDP).LT.0)GO TO 50                                       ACVI6848
      J123=JTB-IQ1                                                      ACVI6849
      GO TO 40                                                          ACVI6850
   35 INDP=INDEX(J1D,L1,J1TP,J2D,L2,J2TP)                               ACVI6851
      IF(J1TA(INDP).LT.0)GO TO 50                                       ACVI6852
      J123=JTA+IQ1                                                      ACVI6853
   40 IF(IP2D.EQ.1)I=IAB(J123)                                          ACVI6854
      IF(IQ2D.EQ.1)I=ICD(NQD-IQ1)                                       ACVI6855
      I=JJ2T*IP1P*(MM1-IP1P)*(LL1+IP1P)*I                               ACVI6856
      DC=DSQRT(DFLOAT(I)/(DFLOAT((J1T+2)*J1TP*NM)*DN))                  ACVI6857
      IF(IQ.LT.0)DC=-DC                                                 ACVI6858
      INDQ=(IND-1)*KR0MAX                                               ACVI6859
      INDPQ=(INDP-1)*KR0MAX                                             ACVI6860
      DO 45 KR0=KR0A,KR0B                                               ACVI6861
      KI=KR0+INDQ                                                       ACVI6862
      KIP=KR0+INDPQ                                                     ACVI6863
   45 DEWU3(KI)=DC*DEWU3(KIP)                                           ACVI6864
   50 IF(IQ1.EQ.L1)GO TO 70                                             ACVI6865
      IQ1P=IQ1+1                                                        ACVI6866
      IP1P=IP1                                                          ACVI6867
      J1TP=J1T-1                                                        ACVI6868
      IF(INC.EQ.0)GO TO 55                                              ACVI6869
      INDP=INDEX(J2D,L2,J2TP,J1D,L1,J1TP)                               ACVI6870
      IF(J1TA(INDP).LT.0)GO TO 70                                       ACVI6871
      J123=JTB-IQ1                                                      ACVI6872
      GO TO 60                                                          ACVI6873
   55 INDP=INDEX(J1D,L1,J1TP,J2D,L2,J2TP)                               ACVI6874
      IF(J1TA(INDP).LT.0)GO TO 70                                       ACVI6875
      J123=JTA+IQ1                                                      ACVI6876
   60 IF(IP2D.EQ.1)I=ICD(NQD-IQ1)                                       ACVI6877
      IF(IQ2D.EQ.1)I=IAB(J123)                                          ACVI6878
      I=JJ2T*IQ1P*(LL1-IQ1P)*(LM1-IQ1P)*I                               ACVI6879
      DC=DSQRT(DFLOAT(I)/(DFLOAT((J1TP+2)*J1T*NM)*DN))                  ACVI6880
      IF(IP.LT.0)DC=-DC                                                 ACVI6881
      INDQ=(IND-1)*KR0MAX                                               ACVI6882
      INDPQ=(INDP-1)*KR0MAX                                             ACVI6883
      DO 65 KR0=KR0A,KR0B                                               ACVI6884
      KI=KR0+INDQ                                                       ACVI6885
      KIP=KR0+INDPQ                                                     ACVI6886
   65 DEWU3(KI)=DEWU3(KI)+DC*DEWU3(KIP)                                 ACVI6887
   70 CONTINUE                                                          ACVI6888
      RETURN                                                            ACVI6889
      END                                                               ACVI6890
      SUBROUTINE XEWU3(LAM1X,MU1X,LAM2X,MU2X,LAM3X,MU3X,I3,NEC,KR0MAX,  ACVI6891
     1INDMAX,DEWU3,J1TA,J2TA,IEA,N1,N2,KIMAX1)                          ACVI6892
C     ------------------------------------------------------------------ACVI6893
C     EXTREMAL WIGNER COEFFICIENTS FOR U3 (X PREFIX FOR 6-81 VERSION)   ACVI6894
C     ------------------------------------------------------------------ACVI6895
C     UPDATE/MOD: (LSU,05-80)  J.P.DRAAYER        LOG BINOMIALS         ACVI6896
C                 (LSU,06-81)  J.P.DRAAYER        INDEXING DEWU3        ACVI6897
C                 (LSU,03-83)  J.P.DRAAYER        SPACE SAVING MEASURE  ACVI6898
C                 (LSU,02-87)  J.P.DRAAYER        OVERFLOW CORRECTION   ACVI6899
C                 (LSU,10-89)  J.P.DRAAYER        ZERO OUT RELOCATED    ACVI6900
C                                                                       ACVI6901
C     REFERENCES--J.P.DRAAYER AND Y.AKIYAMA, J.MATH.PHYS.14(1973)1904   ACVI6902
C                 K.T.HECHT, NUCL.PHYS.62(1965)1                        ACVI6903
C     PARAMETERS--(I3) : (1)=GHW, (0)=GLW                               ACVI6904
C       EXTERNAL--N1=MAX(KR0MAX)                                        ACVI6905
C                 N2=MAX(INDMAX)=NX*(NX+1)*(NX+2)/6, NX=MAX(LAM2+MU2+1) ACVI6906
C*                KIMAX1=MAX(KR0MAX*INDMAX)                             ACVI6907
C       INTERNAL--X1=ABS(N1*NX)                                         ACVI6908
C                 X2=ABS(NX)                                            ACVI6909
C     EXTENSIONS--CHANGE EXTERNAL PARAMETERS IN CALL                    ACVI6910
C                 ADJUST INTERNAL PARAMETERS BELOW                      ACVI6911
C*    DIMENSIONS--DEWU3(N1*N2->KIMAX1),J1TA(N2),J2TA(N2),IEA(N2),       ACVI6912
C                 DEWU3P(X1),DZ(X2),J1TAP(X2),IAB(X2),ICD(X2)           ACVI6913
C       COMMENTS--ASSUME MAX N1=9,NX=42,N2=13244                        ACVI6914
C                        SET X1=378,X2=42                               ACVI6915
C                 DZ ARRAY ADDED FOR CORRECTING THE OVERFLOW PROBLEM    ACVI6916
C     ------------------------------------------------------------------ACVI6917
      IMPLICIT REAL*8(D)                                                ACVI6918
      COMMON/BKDB/DBINO(8385),DBINV(8385),DTWOS(-128:128)               ACVI6919
      COMMON/BKDF/DLOGF(0:2000)                                         ACVI6920
      DIMENSION DEWU3(1),J1TA(1),J2TA(1),IEA(1),                        ACVI6921
     1          DEWU3P(378),DZ(42),J1TAP(42),IAB(42),ICD(42)            ACVI6922
      INDEX(J1TD,LAM1,J1T,J2TD,LAM2,J2T)=1+J2TD*(J2TD+1)*(3*J1TD+J2TD+5)ACVI6923
     1/6+(J1TD+1)*(LAM2+J2TD-J2T)/2+(LAM1+J1TD-J1T)/2                   ACVI6924
C     DIMENSION CHECKS (LSU,6-81)-START                                 ACVI6925
      IF(N1.GT.9)GO TO 200                                              ACVI6926
      NX=XLAM2+XMU2+1                                                   ACVI6927
      IF(NX.GT.42)GO TO 210                                             ACVI6928
C     DIMENSION CHECKS (LSU,6-81)-START                                 ACVI6929
      KR0MAX=MULTU3(LAM1X,MU1X,LAM2X,MU2X,LAM3X,MU3X)                   ACVI6930
      IF(KR0MAX.EQ.0)RETURN                                             ACVI6931
      IF(I3.EQ.1)GO TO 10                                               ACVI6932
      LAM1=MU1X                                                         ACVI6933
      LAM2=MU2X                                                         ACVI6934
      LAM3=MU3X                                                         ACVI6935
      MU1=LAM1X                                                         ACVI6936
      MU2=LAM2X                                                         ACVI6937
      MU3=LAM3X                                                         ACVI6938
      GO TO 15                                                          ACVI6939
   10 LAM1=LAM1X                                                        ACVI6940
      LAM2=LAM2X                                                        ACVI6941
      LAM3=LAM3X                                                        ACVI6942
      MU1=MU1X                                                          ACVI6943
      MU2=MU2X                                                          ACVI6944
      MU3=MU3X                                                          ACVI6945
   15 NEC=(LAM1+LAM2-LAM3+2*(MU1+MU2-MU3))/3                            ACVI6946
      IAH=(LAM2+LAM3-LAM1-NEC)/2                                        ACVI6947
      IBH=(LAM3+LAM1-LAM2+NEC+2)/2                                      ACVI6948
      ICH=(LAM1+LAM2-LAM3-NEC)/2                                        ACVI6949
      IDH=(LAM1+LAM2+LAM3-NEC+2)/2                                      ACVI6950
      DO 20 I=1,NEC                                                     ACVI6951
      IAB(I)=(IAH+I)*(IBH-I)                                            ACVI6952
   20 ICD(I)=(ICH+I)*(IDH+I)                                            ACVI6953
      NCDMAX=MULTHY(LAM1,MU1,LAM2,MU2,LAM3,MU3)                         ACVI6954
      NEC=NEC-NCDMAX                                                    ACVI6955
      LAM2=LAM2-NCDMAX                                                  ACVI6956
      MU2=MU2-NCDMAX                                                    ACVI6957
      NCDMIN=1                                                          ACVI6958
   25 IF(NCDMIN.EQ.NCDMAX)GO TO 30                                      ACVI6959
      IF(MULTU3(LAM1,MU1,LAM2+1,MU2+1,LAM3,MU3).GT.0)GO TO 30           ACVI6960
      NEC=NEC+1                                                         ACVI6961
      LAM2=LAM2+1                                                       ACVI6962
      MU2=MU2+1                                                         ACVI6963
      NCDMIN=NCDMIN+1                                                   ACVI6964
      GO TO 25                                                          ACVI6965
C     DIMENSION MODIFICATION (LSU,6-81)-START                           ACVI6966
   30 NNCMAX=NEC+NCDMAX-NCDMIN+2                                        ACVI6967
      KITEST=KR0MAX*(NNCMAX)*(NNCMAX+1)*(NNCMAX+2)/6                    ACVI6968
      IF(KITEST.GT.KIMAX1)GO TO 220                                     ACVI6969
C     DIMENSION MODIFICATION (LSU,6-81)--STOP                           ACVI6970
      DO I=1,KITEST                                                     ACVI6971
         DEWU3(I)=0.D0                                                  ACVI6972
      ENDDO                                                             ACVI6973
      LL1=LAM1+1                                                        ACVI6974
      MM1=MU1+1                                                         ACVI6975
      LL2=LAM2+1                                                        ACVI6976
      MM2=MU2+1                                                         ACVI6977
      IA1=2*LAM1+4*MU1                                                  ACVI6978
      IB1=4*LAM1+2*MU1                                                  ACVI6979
      IC1=IB1-IA1                                                       ACVI6980
      IA2=2*LAM2+4*MU2                                                  ACVI6981
      IB2=4*LAM2+2*MU2                                                  ACVI6982
      IC2=IB2-IA2                                                       ACVI6983
      IS1=LL1+MM1                                                       ACVI6984
      IS2=LL2+MM2                                                       ACVI6985
      ISS=MM1+LAM2+MU2-NEC                                              ACVI6986
      IE3=-(LAM3+2*MU3)                                                 ACVI6987
      IEH=-(LAM2+2*MU2+3)                                               ACVI6988
      KR0CNT=0                                                          ACVI6989
      DO 135 NCD=NCDMIN,NCDMAX                                          ACVI6990
      NEC=NEC+1                                                         ACVI6991
      LAM2=LAM2+1                                                       ACVI6992
      MU2=MU2+1                                                         ACVI6993
      NNC=NEC+1                                                         ACVI6994
      INDMAX=NNC*(NNC+1)*(NNC+2)/6                                      ACVI6995
      IA2=IA2+6                                                         ACVI6996
      IB2=IB2+6                                                         ACVI6997
      IS2=IS2+2                                                         ACVI6998
      ISS=ISS+1                                                         ACVI6999
      IEH=IEH-3                                                         ACVI7000
      LL2=LAM2+1                                                        ACVI7001
      MM2=MU2+1                                                         ACVI7002
      LN1=LAM1+NEC                                                      ACVI7003
      LN2=LAM2+NEC                                                      ACVI7004
      INN=NEC*NNC/2                                                     ACVI7005
      IF(NCD.EQ.NCDMIN)GO TO 40                                         ACVI7006
      DO 35 I=1,KITEST                                                  ACVI7007
   35 DEWU3(I)=0.D0                                                     ACVI7008
   40 DO 45 IND=1,INDMAX                                                ACVI7009
      IEA(IND)=-1000                                                    ACVI7010
      J2TA(IND)=-1000                                                   ACVI7011
   45 J1TA(IND)=-1000                                                   ACVI7012
      IE2=IEH                                                           ACVI7013
      I=1000                                                            ACVI7014
      DO 55 IIE=1,NNC                                                   ACVI7015
      IE2=IE2+3                                                         ACVI7016
      IE1=IE3-IE2                                                       ACVI7017
      J2TD=IIE-1                                                        ACVI7018
      J1TD=NNC-IIE                                                      ACVI7019
      JJ2TA=IA2-IE2                                                     ACVI7020
      JJ2TB=IB2+IE2                                                     ACVI7021
      IF(JJ2TB.LT.JJ2TA)JJ2TA=JJ2TB                                     ACVI7022
      JJ2TA=JJ2TA/3+1                                                   ACVI7023
      JJ2TB=JJ2TA-IABS(IC2-IE2)/3                                       ACVI7024
      JJ1TA=IA1-IE1                                                     ACVI7025
      JJ1TB=IB1+IE1                                                     ACVI7026
      IF(JJ1TB.LT.JJ1TA)JJ1TA=JJ1TB                                     ACVI7027
      JJ1TA=JJ1TA/3+1                                                   ACVI7028
      JJ1TB=JJ1TA-IABS(IC1-IE1)/3                                       ACVI7029
      J=0                                                               ACVI7030
      DO 50 JJ2T=1,JJ2TB,2                                              ACVI7031
      J2T=JJ2TA-JJ2T                                                    ACVI7032
      L=IABS(J2T-LAM3)                                                  ACVI7033
      M=J2T+LAM3                                                        ACVI7034
      DO 50 JJ1T=1,JJ1TB,2                                              ACVI7035
      J1T=JJ1TA-JJ1T                                                    ACVI7036
      IF(J1T.LT.L)GO TO 50                                              ACVI7037
      IF(J1T.GT.M)GO TO 50                                              ACVI7038
      IND=INDEX(J1TD,LAM1,J1T,J2TD,LAM2,J2T)                            ACVI7039
      IEA(IND)=IE2                                                      ACVI7040
      J2TA(IND)=J2T                                                     ACVI7041
      J1TA(IND)=J1T                                                     ACVI7042
      J=J+1                                                             ACVI7043
   50 CONTINUE                                                          ACVI7044
      IF(J.LT.I)I=J                                                     ACVI7045
   55 CONTINUE                                                          ACVI7046
      IF(I.EQ.0)GO TO 135                                               ACVI7047
      IF(KR0CNT.EQ.0)GO TO 80                                           ACVI7048
C     GENERATE <(LAM1,MU1)????;(LAM2,MU2)HIGH::KR0(LAM3,MU3)HIGH>       ACVI7049
C     FROM <(LAM1,MU1)????;(LAM2-1,MU2-1)HIGH::KR0(LAM3,MU3)HIGH>       ACVI7050
      INDQ=-KR0MAX                                                      ACVI7051
      DO 75 IND=1,NNC                                                   ACVI7052
      INDQ=INDQ+KR0MAX                                                  ACVI7053
      J1T=J1TA(IND)                                                     ACVI7054
      IF(J1T.LT.0)GO TO 75                                              ACVI7055
      IQ1=(LN1-J1T)/2                                                   ACVI7056
      IF(IQ1.LE.0)GO TO 65                                              ACVI7057
      J1TP=J1T+1                                                        ACVI7058
      INDP=(LN1-J1TP-1)/2+1                                             ACVI7059
      IF(J1TAP(INDP).LT.0)GO TO 65                                      ACVI7060
      I=IAB(IQ1)*IQ1*(LL1-IQ1)*(IS1-IQ1)                                ACVI7061
      DC=-DSQRT(DFLOAT(I)/DFLOAT((J1T+2)*J1TP))                         ACVI7062
      INDPQ=(INDP-1)*KR0MAX                                             ACVI7063
      DO 60 KR0=1,KR0CNT                                                ACVI7064
      KI=KR0+INDQ                                                       ACVI7065
      KIP=KR0+INDPQ                                                     ACVI7066
   60 DEWU3(KI)=DC*DEWU3P(KIP)                                          ACVI7067
   65 IP1=NEC-IQ1                                                       ACVI7068
      IF(IP1.LE.0)GO TO 75                                              ACVI7069
      J1TP=J1T-1                                                        ACVI7070
      INDP=(LN1-J1TP-1)/2+1                                             ACVI7071
      IF(J1TAP(INDP).LT.0)GO TO 75                                      ACVI7072
      I=ICD(NNC-IND)*IP1*(MM1-IP1)*(LL1+IP1)                            ACVI7073
      DC=DSQRT(DFLOAT(I)/DFLOAT((J1TP+2)*J1T))                          ACVI7074
      INDPQ=(INDP-1)*KR0MAX                                             ACVI7075
      DO 70 KR0=1,KR0CNT                                                ACVI7076
      KI=KR0+INDQ                                                       ACVI7077
      KIP=KR0+INDPQ                                                     ACVI7078
   70 DEWU3(KI)=DEWU3(KI)+DC*DEWU3P(KIP)                                ACVI7079
   75 CONTINUE                                                          ACVI7080
      INC=0                                                             ACVI7081
      IF(KR0CNT.EQ.KR0MAX)GO TO 125                                     ACVI7082
C     EVALUATE <(LAM1,MU1)HIGH;(LAM2,MU2)????::KR0(LAM3,MU3)HIGH>       ACVI7083
C     WITH (LAM2,MU2) A MINIMUM FOR KR0=KR0CNT                          ACVI7084
   80 KR0CNT=KR0CNT+1                                                   ACVI7085
      I=0                                                               ACVI7086
      IND=INDEX(0,LAM1,LAM1,NEC,LAM2,LN2)-1                             ACVI7087
      INDQ=-KR0MAX                                                      ACVI7088
      DO 85 IIQ2=1,NNC                                                  ACVI7089
      INDQ=INDQ+KR0MAX                                                  ACVI7090
      IND=IND+1                                                         ACVI7091
      KI=KR0CNT+INDQ                                                    ACVI7092
      DEWU3P(KI)=0.D0                                                   ACVI7093
      IF(J1TA(IND).LT.0)GO TO 85                                        ACVI7094
      I=I+1                                                             ACVI7095
      IIQ2B=IIQ2                                                        ACVI7096
   85 CONTINUE                                                          ACVI7097
C                                                                       ACVI7098
C     *****MODIFIED TO AVOID OVERFLOW (LSU,5-80)-START                  ACVI7099
C      ****DLOGB CHANGED TO SAVE SPACE (LSU,3-83)****                   ACVI7100
C      ****FURTHER OVERFLOW CORRECTION (LSU,2-87)****                   ACVI7101
C                                                                       ACVI7102
      IIQ2A=IIQ2B-I+1                                                   ACVI7103
      IQ2B=IIQ2B-1                                                      ACVI7104
      INDQ=(IIQ2A-2)*KR0MAX                                             ACVI7105
      IZ=0                                                              ACVI7106
      DO 115 IIQ2=IIQ2A,IIQ2B                                           ACVI7107
      IZ=IZ+1                                                           ACVI7108
      DZ(IZ)=1.D0                                                       ACVI7109
      INDQ=INDQ+KR0MAX                                                  ACVI7110
      L=LL2-IIQ2                                                        ACVI7111
C                                                                       ACVI7112
C --> START NUMERATOR PRODUCT LOOP                                      ACVI7113
C                                                                       ACVI7114
      IX=L-IIQ2+NNC+1                                                   ACVI7115
      IF(IX.EQ.0)GO TO 120                                              ACVI7116
      IY=IABS(IX)                                                       ACVI7117
      IN=IX/IY                                                          ACVI7118
      DN=DLOG(DFLOAT(IY))                                               ACVI7119
      IF(IIQ2A.EQ.IIQ2B)GO TO 95                                        ACVI7120
      DO 90 I=IIQ2A,IQ2B                                                ACVI7121
      J=NNC-I                                                           ACVI7122
      IF(I.LT.IIQ2)THEN                                                 ACVI7123
      K=IAB(J)*ICD(J)*(IS2-I)                                           ACVI7124
      ELSE                                                              ACVI7125
      K=MM2-J                                                           ACVI7126
      ENDIF                                                             ACVI7127
      IF(K.EQ.0)GO TO 120                                               ACVI7128
      IF(K.LT.0)IN=-IN                                                  ACVI7129
   90 DN=DLOG(DFLOAT(IABS(K)))+DN                                       ACVI7130
   95 DN=DN+DLOG(DBINO(INN+IIQ2))                                       ACVI7131
C                                                                       ACVI7132
C --> END NUMERATOR PRODUCT LOOP & START DENOMINATOR PRODUCT LOOP       ACVI7133
C                                                                       ACVI7134
      ID=1                                                              ACVI7135
      DD=0.D0                                                           ACVI7136
      DO 100 I=1,NNC                                                    ACVI7137
      IX=I+L                                                            ACVI7138
      IF(IX.LT.0)ID=-ID                                                 ACVI7139
  100 DD=DLOG(DFLOAT(I+L))+DD                                           ACVI7140
C                                                                       ACVI7141
C --> END DENOMINATOR PRODUCT LOOP & START INNER PRODUCT/SUM LOOP       ACVI7142
C                                                                       ACVI7143
      IP2=NNC-IIQ2                                                      ACVI7144
C                                                                       ACVI7145
C     MULTIPLY BY SMALL NUMBER --> DEXP(-172) LIMIT FOR IBM SYSTEMS     ACVI7146
C                                                                       ACVI7147
      DZ(IZ)=DEXP(-DMIN1(DLOGF(2*IP2),172.D0))                          ACVI7148
      IIP2=IP2+1                                                        ACVI7149
      M=IP2*IIP2/2                                                      ACVI7150
      DS=0.D0                                                           ACVI7151
      DO 110 I=1,IIP2                                                   ACVI7152
      DC=DZ(IZ)*DBINO(I+M)                                              ACVI7153
      IF(IIP2.EQ.1)GO TO 110                                            ACVI7154
      DO 105 J=1,IP2                                                    ACVI7155
      IF(J.LT.I)THEN                                                    ACVI7156
         K=(J+L)*(ISS+J)                                                ACVI7157
      ELSE                                                              ACVI7158
         K=IAB(J)                                                       ACVI7159
      ENDIF                                                             ACVI7160
  105 DC=DFLOAT(K)*DC                                                   ACVI7161
  110 DS=DS+DC                                                          ACVI7162
C                                                                       ACVI7163
C --> END INNER PRODUCT/SUM LOOP & ASSIGN unnormalized DEWU3P VALUE     ACVI7164
C                                                                       ACVI7165
      IF(2*(IP2/2).NE.IP2)DS=-DS                                        ACVI7166
      KI=KR0CNT+INDQ                                                    ACVI7167
  115 DEWU3P(KI)=DFLOAT(IN*ID)*DS*DEXP((DN-DD)/2.D0)                    ACVI7168
C                                                                       ACVI7169
C --> START renormalization PROCEDURE                                   ACVI7170
C                                                                       ACVI7171
      DMIN=1.D0                                                         ACVI7172
      IZ=0                                                              ACVI7173
      DO 117 IIQ2=IIQ2A,IIQ2B                                           ACVI7174
      IZ=IZ+1                                                           ACVI7175
  117 DMIN=DMIN1(DMIN,DZ(IZ))                                           ACVI7176
      INDQ=(IIQ2A-2)*KR0MAX                                             ACVI7177
      IZ=0                                                              ACVI7178
      DO 118 IIQ2=IIQ2A,IIQ2B                                           ACVI7179
      IZ=IZ+1                                                           ACVI7180
      INDQ=INDQ+KR0MAX                                                  ACVI7181
      KI=KR0CNT+INDQ                                                    ACVI7182
  118 DEWU3P(KI)=(DMIN/DZ(IZ))*DEWU3P(KI)                               ACVI7183
C                                                                       ACVI7184
C     *****MODIFIED TO AVOID OVERFLOW (LSU,5-80)--STOP                  ACVI7185
C                                                                       ACVI7186
      KR0A=KR0CNT                                                       ACVI7187
      KR0B=KR0CNT                                                       ACVI7188
C     GENERATE <(LAM1,MU1)????;(LAM2,MU2)????::KR0(LAM3,MU3)HIGH>       ACVI7189
C     FROM <(LAM1,MU1)HIGH;(LAM2,MU2)????::KR0(LAM3,MU3)HIGH>           ACVI7190
      CALL XEWU3S(1,LAM1,MU1,LAM2,MU2,NEC,NNC,KR0A,KR0B,DEWU3P,J1TA,    ACVI7191
     1IAB,ICD,INDMAX,DEWU3,KR0MAX)                                      ACVI7192
      INC=1                                                             ACVI7193
      GO TO 125                                                         ACVI7194
  120 KR0CNT=KR0CNT-1                                                   ACVI7195
      WRITE(1,195)XLAM1,XMU1,XLAM2,XMU2,XLAM3,XMU3,KR0MAX,KR0CNT        ACVI7196
  125 IF(KR0CNT.EQ.0)GO TO 135                                          ACVI7197
      INDQ=-KR0MAX                                                      ACVI7198
      DO 130 IND=1,NNC                                                  ACVI7199
      INDQ=INDQ+KR0MAX                                                  ACVI7200
      J1TAP(IND)=J1TA(IND)                                              ACVI7201
      DO 130 KR0=1,KR0CNT                                               ACVI7202
      KI=KR0+INDQ                                                       ACVI7203
  130 DEWU3P(KI)=DEWU3(KI)                                              ACVI7204
  135 CONTINUE                                                          ACVI7205
      IF(KR0CNT.EQ.0)RETURN                                             ACVI7206
      KR0A=1                                                            ACVI7207
      KR0B=KR0CNT-INC                                                   ACVI7208
      IF(KR0B.EQ.0)GO TO 140                                            ACVI7209
C     GENERATE <(LAM1,MU1)????;(LAM2,MU2)????::KR0(LAM3,MU3)HIGH>       ACVI7210
C     FROM <(LAM1,MU1)????;(LAM2,MU2)HIGH::KR0(LAM3,MU3)HIGH>           ACVI7211
      CALL XEWU3S(0,LAM1,MU1,LAM2,MU2,NEC,NNC,KR0A,KR0B,DEWU3P,J1TA,    ACVI7212
     1IAB,ICD,INDMAX,DEWU3,KR0MAX)                                      ACVI7213
  140 CONTINUE                                                          ACVI7214
C     RENORMALIZE VIA LARGEST ELEMENT TO AVOID OVERFLOW (LSU,5-80)-STARTACVI7215
      DO 142 KR0=1,KR0CNT                                               ACVI7216
      DC=1.D0                                                           ACVI7217
      INDQ=-KR0MAX                                                      ACVI7218
      DO 141 IND=1,INDMAX                                               ACVI7219
      INDQ=INDQ+KR0MAX                                                  ACVI7220
      KI=KR0+INDQ                                                       ACVI7221
  141 DC=DMAX1(DC,DABS(DEWU3(KI)))                                      ACVI7222
      INDQ=-KR0MAX                                                      ACVI7223
      DO 142 IND=1,INDMAX                                               ACVI7224
      INDQ=INDQ+KR0MAX                                                  ACVI7225
      KI=KR0+INDQ                                                       ACVI7226
  142 DEWU3(KI)=DEWU3(KI)/DC                                            ACVI7227
C     RENORMALIZE VIA LARGEST ELEMENT TO AVOID OVERFLOW (LSU,5-80)--STOPACVI7228
C     ORTHONORMALIZATION OF SOLUTIONS                                   ACVI7229
      DO 165 KR0=1,KR0CNT                                               ACVI7230
      KR0PA=1                                                           ACVI7231
      IF(INC.EQ.1.AND.KR0.EQ.KR0CNT)KR0PA=KR0CNT                        ACVI7232
      DO 165 KR0P=KR0PA,KR0                                             ACVI7233
      DN=0.D0                                                           ACVI7234
      INDQ=-KR0MAX                                                      ACVI7235
      DO 145 IND=1,INDMAX                                               ACVI7236
      INDQ=INDQ+KR0MAX                                                  ACVI7237
      KI=KR0+INDQ                                                       ACVI7238
      KIP=KR0P+INDQ                                                     ACVI7239
  145 DN=DN+DEWU3(KI)*DEWU3(KIP)                                        ACVI7240
      IF(KR0P.EQ.KR0)GO TO 155                                          ACVI7241
      IF(DABS(DN).LT.1.D-12)GO TO 165                                   ACVI7242
      INDQ=-KR0MAX                                                      ACVI7243
      DO 150 IND=1,INDMAX                                               ACVI7244
      INDQ=INDQ+KR0MAX                                                  ACVI7245
      KI=KR0+INDQ                                                       ACVI7246
      KIP=KR0P+INDQ                                                     ACVI7247
  150 DEWU3(KI)=DEWU3(KI)-DN*DEWU3(KIP)                                 ACVI7248
      GO TO 165                                                         ACVI7249
  155 DN=1.D0/DSQRT(DN)                                                 ACVI7250
      INDQ=-KR0MAX                                                      ACVI7251
      DO 160 IND=1,INDMAX                                               ACVI7252
      INDQ=INDQ+KR0MAX                                                  ACVI7253
      KI=KR0+INDQ                                                       ACVI7254
  160 DEWU3(KI)=DN*DEWU3(KI)                                            ACVI7255
  165 CONTINUE                                                          ACVI7256
C     SET PHASE CONVENTION (K.T.HECHT, NUCL.PHYS.62(1965)1)             ACVI7257
      IE2=IE3+(LAM1+2*MU1)                                              ACVI7258
      IPH=2*(LAM1+LAM2-LAM3+MU1+MU2-MU3+KR0MAX)                         ACVI7259
      INDQ=-KR0MAX                                                      ACVI7260
      DO 180 IND=1,INDMAX                                               ACVI7261
      INDQ=INDQ+KR0MAX                                                  ACVI7262
      IF(IEA(IND).NE.IE2)GO TO 180                                      ACVI7263
      I=IPH+J1TA(IND)+J2TA(IND)-LAM3                                    ACVI7264
      DO 175 KR0=1,KR0MAX                                               ACVI7265
      I=I-2                                                             ACVI7266
      KI=KR0+INDQ                                                       ACVI7267
      J=I                                                               ACVI7268
      IF(DEWU3(KI).LT.0.D0)J=J-2                                        ACVI7269
      IF(4*(J/4).EQ.J)GO TO 175                                         ACVI7270
      INDPQ=-KR0MAX                                                     ACVI7271
      DO 170 INDP=1,INDMAX                                              ACVI7272
      INDPQ=INDPQ+KR0MAX                                                ACVI7273
      KIP=KR0+INDPQ                                                     ACVI7274
  170 DEWU3(KIP)=-DEWU3(KIP)                                            ACVI7275
  175 CONTINUE                                                          ACVI7276
      GO TO 185                                                         ACVI7277
  180 CONTINUE                                                          ACVI7278
  185 IF(I3.EQ.1)RETURN                                                 ACVI7279
      INDQ=-KR0MAX                                                      ACVI7280
      DO 190 IND=1,INDMAX                                               ACVI7281
      INDQ=INDQ+KR0MAX                                                  ACVI7282
      IEA(IND)=-IEA(IND)                                                ACVI7283
      I=IPH+J1TA(IND)+J2TA(IND)-LAM3                                    ACVI7284
      DO 190 KR0=1,KR0MAX                                               ACVI7285
      I=I-2                                                             ACVI7286
      KI=KR0+INDQ                                                       ACVI7287
      IF(4*(I/4).NE.I)DEWU3(KI)=-DEWU3(KI)                              ACVI7288
  190 CONTINUE                                                          ACVI7289
  195 FORMAT(28H *****U3 COUPLING ERROR*****,3X,3(4X,2I3),3X,           ACVI7290
     115HRO(ABSOLUTE) = ,I2,3X,15HRO(RELATIVE) = ,I2,4X,                ACVI7291
     226H*****REPORT TO AUTHOR*****)                                    ACVI7292
      RETURN                                                            ACVI7293
C     DIMENSION CHECKS (LSU,6-81)-START                                 ACVI7294
  200 WRITE(6,205)N1                                                    ACVI7295
  205 FORMAT(36H ***** XEWU3 DIMENSION OVERFLOW: N1=,I10)               ACVI7296
      GO TO 230                                                         ACVI7297
  210 WRITE(6,215)NX                                                    ACVI7298
  215 FORMAT(36H ***** XEWU3 DIMENSION OVERFLOW: NX=,I10)               ACVI7299
      GO TO 230                                                         ACVI7300
  220 WRITE(6,225)KITEST,KIMAX1                                         ACVI7301
  225 FORMAT(40H ***** XEWU3 DIMENSION OVERFLOW: KITEST=,I10,5X,        ACVI7302
     17HKIMAX1=,I10)                                                    ACVI7303
  230 STOP                                                              ACVI7304
C     DIMENSION CHECKS (LSU,6-81)--STOP                                 ACVI7305
      END                                                               ACVI7306
      SUBROUTINE XWU3(LAM1,MU1,LAM2,MU2,LAM3,MU3,IE,JT,NEC,DEWU3,       ACVI7307
     1KR0MAX,INDMAX,DWU3,J1SMAX,J1TMAX,J2SMAX,J2TMAX,IESMAX,IE2MAX,     ACVI7308
     2INDMAT,N1,N2,N3,KIMAX2)                                           ACVI7309
C     ------------------------------------------------------------------ACVI7310
C     WIGNER COEFFICIENTS FOR U3 (X PREFIX FOR 6-81 VERSION)            ACVI7311
C     ------------------------------------------------------------------ACVI7312
C     UPDATE/MOD: (LSU,06-81)  J.P.DRAAYER        INDEXING OF DEWU3     ACVI7313
C                 (LSU,11-89)  J.P.DRAAYER        DWU3 ZERO-OUT RANGE   ACVI7314
C                                                                       ACVI7315
C     REFERENCES--J.P.DRAAYER AND Y.AKIYAMA, J.MATH.PHYS.14(1973)1904   ACVI7316
C                 K.T.HECHT, NUCL.PHYS.62(1965)1                        ACVI7317
C     PARAMETERS--SEE ALSO XEWU3                                        ACVI7318
C       EXTERNAL--N1=MAX(KR0MAX)                                        ACVI7319
C                 N2=MAX(IESMAX) SAFE TO SET N2=MAX(LAM2+MU2+1)         ACVI7320
C                 N3=MAX(DIM(LAM2,MU2))                                 ACVI7321
C                 NA=MAX(INDMAX)=NX*(NX+1)*(NX+2)/6, NX=MAX(LAM2+MU2+1) ACVI7322
C                 NB=MAX(J2SMAX) SAFE TO SET NB=MAX(LAM2+MU2+1)         ACVI7323
C*                KIMAX1=MAX(KR0MAX*INDMAX) (SEE XEWU3)                 ACVI7324
C*                KIMAX2=MAX(KR0MAX*DIM(LAM2,MU2))                      ACVI7325
C       INTERNAL--X1=ABS(N1*N3)                                         ACVI7326
C                 X2=ABS(N2)                                            ACVI7327
C                 X3=ABS(N2*NB)                                         ACVI7328
C     EXTENSIONS--CHANGE EXTERNAL PARAMETERS IN CALL STATEMENT          ACVI7329
C                 ADJUST INTERNAL PARAMETERS BELOW                      ACVI7330
C*    DIMENSIONS--DEWU3(N1*NA),DWU3(N1*N3),J1SMAX(N2*NB),               ACVI7331
C                 J1TMAX(N2*NB),J2SMAX(N2),J2TMAX(N2),INDMAT(N2*NB),    ACVI7332
C                 DWU3P(X1),J2TMAP(X2),INDMAP(X3)                       ACVI7333
C*      COMMENTS--USE N1*NA->KIMAX1,N1*N3->KIMAX2                       ACVI7334
C                 ASSUME MAX N1=9,N2=42,N3=9030                         ACVI7335
C                        SET X1=27090,X2=42,X3=1764 (X1=3*N3,FIXED)     ACVI7336
C     ------------------------------------------------------------------ACVI7337
      IMPLICIT REAL*8(D),INTEGER(X)                                     ACVI7338
C     IMPLICIT INTEGER(X)                                               ACVI7339
      DIMENSION DEWU3(1),DWU3(1),J1SMAX(1),                             ACVI7340
     1          J1TMAX(1),J2SMAX(1),J2TMAX(1),INDMAT(1),                ACVI7341
     2          DWU3P(27090),J2TMAP(42),INDMAP(1764)                    ACVI7342
      INDEX(J1TD,LAM1,J1T,J2TD,LAM2,J2T)=1+J2TD*(J2TD+1)*(3*J1TD+J2TD+5)ACVI7343
     1/6+(J1TD+1)*(LAM2+J2TD-J2T)/2+(LAM1+J1TD-J1T)/2                   ACVI7344
      IDM(LAM,MU)=(LAM+1)*(MU+1)*(LAM+MU+2)/2                           ACVI7345
C     DIMENSION CHECKS (LSU,6-81)-START                                 ACVI7346
      IF(N1.GT.9)GO TO 130                                              ACVI7347
      IF(N2.GT.42)GO TO 140                                             ACVI7348
      IF(N3.GT.9030)GO TO 150                                           ACVI7349
      IDTEST=KR0MAX*IDM(LAM2,MU2)                                       ACVI7350
      IF(IDTEST.GT.KIMAX2.OR.IDTEST.GT.27090)GO TO 160                  ACVI7351
C     DIMENSION CHECKS (LSU,6-81)--STOP                                 ACVI7352
      LL1=LAM1+1                                                        ACVI7353
      MM1=MU1+1                                                         ACVI7354
      LL2=LAM2+1                                                        ACVI7355
      MM2=MU2+1                                                         ACVI7356
      LL3=LAM3+1                                                        ACVI7357
      MM3=MU3+1                                                         ACVI7358
      LM1=LAM1+MU1                                                      ACVI7359
      LM2=LAM2+MU2                                                      ACVI7360
      LLMM1=LL1+MM1                                                     ACVI7361
      LLMM2=LL2+MM2                                                     ACVI7362
      LLMM3=LL3+MM3                                                     ACVI7363
      JJTD=(IE+LAM3+2*MU3)/3+1                                          ACVI7364
      IP=(JJTD+JT-LL3)/2                                                ACVI7365
      NCC=NEC-1                                                         ACVI7366
      INC=1                                                             ACVI7367
      IQ3=0                                                             ACVI7368
      IP3=-1                                                            ACVI7369
      J3T=LAM3+IP3-IQ3                                                  ACVI7370
      DO 125 JJ3TD=1,JJTD                                               ACVI7371
C     DO 10 N=1,KIMAX2                                                  ACVI7372
      DO 10 N=1,IDTEST                                                  ACVI7373
   10 DWU3(N)=0.D0                                                      ACVI7374
      NCC=NCC+1                                                         ACVI7375
      IF(IP3.EQ.IP)INC=0                                                ACVI7376
      IF(INC.EQ.1)GO TO 15                                              ACVI7377
      IQ3=IQ3+1                                                         ACVI7378
      J3T=J3T-1                                                         ACVI7379
      NM=(LL3-IQ3)*IQ3*(LLMM3-IQ3)                                      ACVI7380
      GO TO 20                                                          ACVI7381
   15 IP3=IP3+1                                                         ACVI7382
      J3T=J3T+1                                                         ACVI7383
      NM=(MM3-IP3)*IP3*(LL3+IP3)                                        ACVI7384
   20 JJ2TDA=NCC-LM1                                                    ACVI7385
      IF(JJ2TDA.LT.0)JJ2TDA=0                                           ACVI7386
      JJ2TDA=JJ2TDA+1                                                   ACVI7387
      JJ2TDB=LM2                                                        ACVI7388
      IF(NCC.LT.JJ2TDB)JJ2TDB=NCC                                       ACVI7389
      JJ2TDB=JJ2TDB+1                                                   ACVI7390
      JJ2TDC=JJ2TDA                                                     ACVI7391
      IND=0                                                             ACVI7392
      IES=0                                                             ACVI7393
      DO 115 JJ2TD=JJ2TDA,JJ2TDB                                        ACVI7394
      J2TD=JJ2TD-1                                                      ACVI7395
      J1TD=NCC-J2TD                                                     ACVI7396
      IES=IES+1                                                         ACVI7397
      IIQ2A=J2TD-MU2                                                    ACVI7398
      IF(IIQ2A.LT.0)IIQ2A=0                                             ACVI7399
      IIQ2A=IIQ2A+1                                                     ACVI7400
      IIQ2B=J2TD                                                        ACVI7401
      IF(LAM2.LT.IIQ2B)IIQ2B=LAM2                                       ACVI7402
      IIQ2B=IIQ2B+1                                                     ACVI7403
      IIQ1A=J1TD-MU1                                                    ACVI7404
      IF(IIQ1A.LT.0)IIQ1A=0                                             ACVI7405
      IIQ1A=IIQ1A+1                                                     ACVI7406
      IIQ1B=J1TD                                                        ACVI7407
      IF(LAM1.LT.IIQ1B)IIQ1B=LAM1                                       ACVI7408
      IIQ1B=IIQ1B+1                                                     ACVI7409
      J2S=0                                                             ACVI7410
      DO 105 IIQ2=IIQ2A,IIQ2B                                           ACVI7411
      IQ2=IIQ2-1                                                        ACVI7412
      IP2=J2TD-IQ2                                                      ACVI7413
      J2T=LAM2+IP2-IQ2                                                  ACVI7414
      J23S=J2T+J3T                                                      ACVI7415
      J23D=J3T-J2T                                                      ACVI7416
      J23H=IABS(J23D)                                                   ACVI7417
      J1S=0                                                             ACVI7418
      DO 100 IIQ1=IIQ1A,IIQ1B                                           ACVI7419
      IQ1=IIQ1-1                                                        ACVI7420
      IP1=J1TD-IQ1                                                      ACVI7421
      J1T=LAM1+IP1-IQ1                                                  ACVI7422
      IF(J1T.LT.J23H)GO TO 100                                          ACVI7423
      IF(J1T.GT.J23S)GO TO 100                                          ACVI7424
      J1TS=J1T                                                          ACVI7425
      J2TS=J2T                                                          ACVI7426
      INDQ=IND*KR0MAX                                                   ACVI7427
      IND=IND+1                                                         ACVI7428
      J1S=J1S+1                                                         ACVI7429
      IF(JJ3TD.EQ.1)GO TO 90                                            ACVI7430
      JA=(J23S-J1T)/2                                                   ACVI7431
      JJA=JA+1                                                          ACVI7432
      JB=(J23D+J1T)/2                                                   ACVI7433
      JJB=JB+1                                                          ACVI7434
      JC=(J1T+J23S)/2+1                                                 ACVI7435
      JJC=JC+1                                                          ACVI7436
      JD=(J1T-J23D)/2+1                                                 ACVI7437
      JJD=JD-1                                                          ACVI7438
      IESP=J2TD-JJ2TDP                                                  ACVI7439
      DO 85 I=1,4                                                       ACVI7440
      IF(I.EQ.1)IESP=IESP+1                                             ACVI7441
      IF(I.EQ.3)IESP=IESP+1                                             ACVI7442
      IF(IESP.LT.1)GO TO 85                                             ACVI7443
      IF(IESP.GT.IESMAX)GO TO 85                                        ACVI7444
      GO TO (25,35,45,55),I                                             ACVI7445
   25 J2TP=J2T+1                                                        ACVI7446
      J1TP=J1T                                                          ACVI7447
      IF(J1TP.LT.IABS(J2TP-J3TP))GO TO 85                               ACVI7448
      IF(J1TP.GT.J2TP+J3TP)GO TO 85                                     ACVI7449
      M=IQ2                                                             ACVI7450
      IF(M.EQ.0)GO TO 85                                                ACVI7451
      N=LL2-M                                                           ACVI7452
      N=(LLMM2-M)*N                                                     ACVI7453
      J12TP=J2T+1                                                       ACVI7454
      IF(INC.EQ.1)GO TO 30                                              ACVI7455
      IAB=JJA                                                           ACVI7456
      ICD=JJC                                                           ACVI7457
      IPH=1                                                             ACVI7458
      GO TO 65                                                          ACVI7459
   30 IAB=JB                                                            ACVI7460
      ICD=JD                                                            ACVI7461
      IPH=-1                                                            ACVI7462
      GO TO 65                                                          ACVI7463
   35 J2TP=J2T-1                                                        ACVI7464
      J1TP=J1T                                                          ACVI7465
      IF(J1TP.LT.IABS(J2TP-J3TP))GO TO 85                               ACVI7466
      IF(J1TP.GT.J2TP+J3TP)GO TO 85                                     ACVI7467
      M=IP2                                                             ACVI7468
      IF(M.EQ.0)GO TO 85                                                ACVI7469
      N=MM2-M                                                           ACVI7470
      N=(LLMM2-N)*N                                                     ACVI7471
      J12TP=J2T                                                         ACVI7472
      IF(INC.EQ.1)GO TO 40                                              ACVI7473
      IAB=JJB                                                           ACVI7474
      ICD=JJD                                                           ACVI7475
      IPH=1                                                             ACVI7476
      GO TO 65                                                          ACVI7477
   40 IAB=JA                                                            ACVI7478
      ICD=JC                                                            ACVI7479
      IPH=1                                                             ACVI7480
      GO TO 65                                                          ACVI7481
   45 J2TP=J2T                                                          ACVI7482
      J1TP=J1T+1                                                        ACVI7483
      IF(J1TP.LT.IABS(J2TP-J3TP))GO TO 85                               ACVI7484
      IF(J1TP.GT.J2TP+J3TP)GO TO 85                                     ACVI7485
      M=IQ1                                                             ACVI7486
      IF(M.EQ.0)GO TO 85                                                ACVI7487
      N=LL1-M                                                           ACVI7488
      N=(LLMM1-M)*N                                                     ACVI7489
      J12TP=J1T+1                                                       ACVI7490
      IF(INC.EQ.1)GO TO 50                                              ACVI7491
      IAB=JJB                                                           ACVI7492
      ICD=JJC                                                           ACVI7493
      IPH=1                                                             ACVI7494
      GO TO 65                                                          ACVI7495
   50 IAB=JA                                                            ACVI7496
      ICD=JD                                                            ACVI7497
      IPH=1                                                             ACVI7498
      GO TO 65                                                          ACVI7499
   55 J2TP=J2T                                                          ACVI7500
      J1TP=J1T-1                                                        ACVI7501
      IF(J1TP.LT.IABS(J2TP-J3TP))GO TO 85                               ACVI7502
      IF(J1TP.GT.J2TP+J3TP)GO TO 85                                     ACVI7503
      M=IP1                                                             ACVI7504
      IF(M.EQ.0)GO TO 85                                                ACVI7505
      N=MM1-M                                                           ACVI7506
      N=(LLMM1-N)*N                                                     ACVI7507
      J12TP=J1T                                                         ACVI7508
      IF(INC.EQ.1)GO TO 60                                              ACVI7509
      IAB=JJA                                                           ACVI7510
      ICD=JJD                                                           ACVI7511
      IPH=-1                                                            ACVI7512
      GO TO 65                                                          ACVI7513
   60 IAB=JB                                                            ACVI7514
      ICD=JC                                                            ACVI7515
      IPH=1                                                             ACVI7516
   65 IF(J12TP.GT.0)GO TO 70                                            ACVI7517
      IF(INC.EQ.1)IAB=1                                                 ACVI7518
      IF(INC.EQ.0)ICD=1                                                 ACVI7519
      DC=1.D0                                                           ACVI7520
      GO TO 75                                                          ACVI7521
   70 DC=DFLOAT(J12TP*(J12TP+1))                                        ACVI7522
   75 J2SP=(J2TMAP(IESP)-J2TP+2)/2                                      ACVI7523
      INDP=(INDMAP(IESP+(J2SP-1)*N2)-J1TP)/2                            ACVI7524
      DC=DSQRT(DFLOAT(IAB*ICD*M*N)/(DFLOAT(NM)*DC))                     ACVI7525
      IF(IPH.LT.0)DC=-DC                                                ACVI7526
      INDPQ=(INDP-1)*KR0MAX                                             ACVI7527
      DO 80 KR0=1,KR0MAX                                                ACVI7528
      KI=KR0+INDQ                                                       ACVI7529
      KIP=KR0+INDPQ                                                     ACVI7530
   80 DWU3(KI)=DWU3(KI)+DC*DWU3P(KIP)                                   ACVI7531
   85 CONTINUE                                                          ACVI7532
      GO TO 100                                                         ACVI7533
   90 INDPQ=(INDEX(J1TD,LAM1,J1T,J2TD,LAM2,J2T)-1)*KR0MAX               ACVI7534
      DO 95 KR0=1,KR0MAX                                                ACVI7535
      KI=KR0+INDQ                                                       ACVI7536
      KIP=KR0+INDPQ                                                     ACVI7537
   95 DWU3(KI)=DEWU3(KIP)                                               ACVI7538
  100 CONTINUE                                                          ACVI7539
      IF(J1S.EQ.0)GO TO 105                                             ACVI7540
      IESJ2S=IES+J2S*N2                                                 ACVI7541
      J2S=J2S+1                                                         ACVI7542
      J1SMAX(IESJ2S)=J1S                                                ACVI7543
      J1TMAX(IESJ2S)=J1TS+2*(J1S-1)                                     ACVI7544
      INDMAT(IESJ2S)=2*IND+J1TS                                         ACVI7545
  105 CONTINUE                                                          ACVI7546
      IF(J2S.NE.0)GO TO 110                                             ACVI7547
      IES=IES-1                                                         ACVI7548
      IF(IES.EQ.0)JJ2TDC=JJ2TDC+1                                       ACVI7549
      GO TO 115                                                         ACVI7550
  110 J2SMAX(IES)=J2S                                                   ACVI7551
      J2TMAX(IES)=J2TS+2*(J2S-1)                                        ACVI7552
  115 CONTINUE                                                          ACVI7553
      IESMAX=IES                                                        ACVI7554
      IF(JJ3TD.EQ.JJTD)GO TO 125                                        ACVI7555
      J3TP=J3T                                                          ACVI7556
      JJ2TDP=JJ2TDC                                                     ACVI7557
      IND=0                                                             ACVI7558
      DO 120 IES=1,IESMAX                                               ACVI7559
      J2TMAP(IES)=J2TMAX(IES)                                           ACVI7560
      J2SB=J2SMAX(IES)                                                  ACVI7561
      J2SQ=-N2                                                          ACVI7562
      DO 120 J2S=1,J2SB                                                 ACVI7563
      J2SQ=J2SQ+N2                                                      ACVI7564
      IESJ2S=IES+J2SQ                                                   ACVI7565
      INDMAP(IESJ2S)=INDMAT(IESJ2S)                                     ACVI7566
      J1SB=J1SMAX(IESJ2S)                                               ACVI7567
      DO 120 J1S=1,J1SB                                                 ACVI7568
      INDQ=IND*KR0MAX                                                   ACVI7569
      IND=IND+1                                                         ACVI7570
      DO 120 KR0=1,KR0MAX                                               ACVI7571
      KI=KR0+INDQ                                                       ACVI7572
  120 DWU3P(KI)=DWU3(KI)                                                ACVI7573
  125 CONTINUE                                                          ACVI7574
      INDMAX=IND                                                        ACVI7575
      IE2MAX=-(LAM2+2*MU2)+3*(JJ2TDC-1)+3*(IES-1)                       ACVI7576
      RETURN                                                            ACVI7577
C     DIMENSION CHECKS (LSU,6-81)-START                                 ACVI7578
  130 WRITE(6,135)N1                                                    ACVI7579
  135 FORMAT(35H *****XWU3 DIMENSION OVERFLOW:  N1=,I10)                ACVI7580
      GO TO 170                                                         ACVI7581
  140 WRITE(6,145)N2                                                    ACVI7582
  145 FORMAT(35H *****XWU3 DIMENSION OVERFLOW:  N2=,I10)                ACVI7583
      GO TO 170                                                         ACVI7584
  150 WRITE(6,155)N3                                                    ACVI7585
  155 FORMAT(35H *****XWU3 DIMENSION OVERFLOW:  N3=,I10)                ACVI7586
      GO TO 170                                                         ACVI7587
  160 WRITE(6,165)IDTEST                                                ACVI7588
  165 FORMAT(39H *****XWU3 DIMENSION OVERFLOW:  IDTEST=,I10)            ACVI7589
  170 STOP                                                              ACVI7590
C     DIMENSION CHECKS (LSU,6-81)--STOP                                 ACVI7591
      END                                                               ACVI7592
      FUNCTION MULTHY(L1,M1,L2,M2,L3,M3)                                ACVI7593
C     ------------------------------------------------------------------ACVI7594
C     MULTIPLICITY (THEORY) IN U3 COUPLING                              ACVI7595
C     ------------------------------------------------------------------ACVI7596
      DIMENSION IX(6)                                                   ACVI7597
      MULTHY=0                                                          ACVI7598
      IX(1)=L1+L2-L3+2*(M1+M2-M3)                                       ACVI7599
      IX(2)=M1+M2-M3+2*(L1+L2-L3)                                       ACVI7600
      IX(3)=2*L2+M2+M1-L1-M3+L3                                         ACVI7601
      IX(4)=2*M2+L2+L1-M1-L3+M3                                         ACVI7602
      IX(5)=L3+M2-L1+2*(M3+L2-M1)                                       ACVI7603
      IX(6)=M3+L2-M1+2*(L3+M2-L1)                                       ACVI7604
      IXMIN=1000                                                        ACVI7605
      DO 10 I=1,6                                                       ACVI7606
      IXDB3=IX(I)/3                                                     ACVI7607
      IF(3*IXDB3.LT.IX(I))RETURN                                        ACVI7608
      IF(IXDB3.LT.IXMIN)IXMIN=IXDB3                                     ACVI7609
   10 CONTINUE                                                          ACVI7610
      IF(IXMIN.LT.0)RETURN                                              ACVI7611
      MULTHY=MIN0(IXMIN,L2,M2)+1                                        ACVI7612
      RETURN                                                            ACVI7613
      END                                                               ACVI7614
C     ------------------------------------------------------------------ACVI7615
C                                                                       ACVI7616
C     AUTHOR: ORIGINAL CODE BY J. P. DRAAYER (U. OF MICHIGAN, 1970-1974)ACVI7617
C                                                                       ACVI7618
C     UPDATE: 02/25/79 --> MODIFIED MTS R3PACK (LOG FACTORIALS INSERTED)ACVI7619
C             01/10/88 --> MODIFIED FOR ENHANCED APPLICATIONS (V/P WORK)ACVI7620
C                          1) DLOGF REDEFINED & RANGE EXTENDED IN BLOCKSACVI7621
C                          2) DLOGF IMPLEMENTED TO AVOID XFLOWS IN DELTAACVI7622
C                          3) BTEST INSERTED TO SIMPLIFY & IMPROVE CODESACVI7623
C                                                                       ACVI7624
C     NOTICE: 1) BLKNEW MUST BE CALLED (ONCE) BEFORE USING THE PROGRAMS.ACVI7625
C             2) THE RANGE OF DLOGF CAN BE EXTENDED FOR HIGHER J VALUES.ACVI7626
C             3) A SIMILAR PACKAGE FOR VECTOR APPLICATIONS IS AVAILABLE.ACVI7627
C                                                                       ACVI7628
C     ------------------------------------------------------------------ACVI7629
      FUNCTION DWR3(J1T,J2T,J3T,M1T,M2T,M3T)                            ACVI7630
C     ------------------------------------------------------------------ACVI7631
C     WIGNER COEFFICIENTS FOR R3--TRIANGLE RELATIONS CHECKED IN DELTA   ACVI7632
C     REFERENCES--ELEMENTARY THEORY OF ANGULAR MOMENTUM, M.E.ROSE, WILEYACVI7633
C     ------------------------------------------------------------------ACVI7634
      IMPLICIT REAL*8(D)                                                ACVI7635
      COMMON/BKDF/DLOGF(0:2000)                                         ACVI7636
      DWR3=0.D0                                                         ACVI7637
      IF(M1T+M2T-M3T.NE.0)GOTO 20                                       ACVI7638
      DC=DELTA(J1T,J2T,J3T)                                             ACVI7639
      IF(DC.EQ.12345D0)GOTO 20                                          ACVI7640
      I1=J3T-J2T+M1T                                                    ACVI7641
      I2=J3T-J1T-M2T                                                    ACVI7642
      I3=J1T+J2T-J3T                                                    ACVI7643
      I4=J1T-M1T                                                        ACVI7644
      IF(BTEST(I4,0))GOTO 20                                            ACVI7645
      I5=J2T+M2T                                                        ACVI7646
      IF(BTEST(I5,0))GOTO 20                                            ACVI7647
      ITMIN=MAX0(0,-I1,-I2)                                             ACVI7648
      ITMAX=MIN0(I3,I4,I5)                                              ACVI7649
      IF(ITMIN.GT.ITMAX)GOTO 20                                         ACVI7650
      DTOP=(DLOG(DFLOAT(J3T+1))+DC+DLOGF(J1T+M1T)+DLOGF(J1T-M1T)+       ACVI7651
     1DLOGF(J2T+M2T)+DLOGF(J2T-M2T)+DLOGF(J3T+M3T)+                     ACVI7652
     2DLOGF(J3T-M3T))/DFLOAT(2)                                         ACVI7653
      DO 10 IT=ITMIN,ITMAX,2                                            ACVI7654
      DBOT=DLOGF(I3-IT)+DLOGF(I4-IT)+DLOGF(I5-IT)+                      ACVI7655
     1DLOGF(IT)+DLOGF(I1+IT)+DLOGF(I2+IT)                               ACVI7656
      DSUM=DEXP(DTOP-DBOT)                                              ACVI7657
      IF(BTEST(IT,1))THEN                                               ACVI7658
      DWR3=DWR3-DSUM                                                    ACVI7659
      ELSE                                                              ACVI7660
      DWR3=DWR3+DSUM                                                    ACVI7661
      ENDIF                                                             ACVI7662
   10 CONTINUE                                                          ACVI7663
   20 RETURN                                                            ACVI7664
      END                                                               ACVI7665
      FUNCTION DRR3(J1T,J2T,L2T,L1T,J3T,L3T)                            ACVI7666
C     ------------------------------------------------------------------ACVI7667
C     RACAH COEFFICIENTS FOR R3--TRIANGLE RELATION CHECKED IN DELTA     ACVI7668
C     REFERENCES--THE 3-J AND 6-J SYMBOLS, M.ROTENBERG, R.BIVINS,       ACVI7669
C                 N.METROPOLIS AND J.K.WOOTEN, MIT PRESS                ACVI7670
C     ------------------------------------------------------------------ACVI7671
      IMPLICIT REAL*8(D)                                                ACVI7672
      COMMON/BKDF/DLOGF(0:2000)                                         ACVI7673
      DRR3=0.D0                                                         ACVI7674
      DX=DELTA(J1T,J2T,J3T)                                             ACVI7675
      IF(DX.EQ.12345D0)GOTO 20                                          ACVI7676
      DC=DX                                                             ACVI7677
      DX=DELTA(L1T,L2T,J3T)                                             ACVI7678
      IF(DX.EQ.12345D0)GOTO 20                                          ACVI7679
      DC=DX+DC                                                          ACVI7680
      DX=DELTA(L1T,J2T,L3T)                                             ACVI7681
      IF(DX.EQ.12345D0)GOTO 20                                          ACVI7682
      DC=DX+DC                                                          ACVI7683
      DX=DELTA(J1T,L2T,L3T)                                             ACVI7684
      IF(DX.EQ.12345D0)GOTO 20                                          ACVI7685
      DC=(DX+DC)/2.D0                                                   ACVI7686
      I1=J3T+L3T-J1T-L1T                                                ACVI7687
      I2=J3T+L3T-J2T-L2T                                                ACVI7688
      I3=J1T+J2T+L1T+L2T+2                                              ACVI7689
      I4=J1T+J2T-J3T                                                    ACVI7690
      I5=L1T+L2T-J3T                                                    ACVI7691
      I6=J1T+L2T-L3T                                                    ACVI7692
      I7=L1T+J2T-L3T                                                    ACVI7693
      ITMIN=MAX0(0,-I1,-I2)                                             ACVI7694
      ITMAX=MIN0(I3,I4,I5,I6,I7)                                        ACVI7695
      IF(ITMIN.GT.ITMAX)GOTO 20                                         ACVI7696
      DO 10 IT=ITMIN,ITMAX,2                                            ACVI7697
      DSUM=DEXP(DC+DLOGF(I3-IT)-(DLOGF(I4-IT)+DLOGF(I5-IT)+             ACVI7698
     1DLOGF(I6-IT)+DLOGF(I7-IT)+DLOGF(IT)+DLOGF(I1+IT)+DLOGF(I2+IT)))   ACVI7699
      IF(BTEST(IT,1))THEN                                               ACVI7700
      DRR3=DRR3-DSUM                                                    ACVI7701
      ELSE                                                              ACVI7702
      DRR3=DRR3+DSUM                                                    ACVI7703
      ENDIF                                                             ACVI7704
   10 CONTINUE                                                          ACVI7705
   20 RETURN                                                            ACVI7706
      END                                                               ACVI7707
      FUNCTION DJHR3(J1T,J2T,J3T,J4T,J5T,J6T,J7T,J8T,J9T)               ACVI7708
C     ------------------------------------------------------------------ACVI7709
C     JAHN-HOPE COEFFICIENTS FOR R3--TRIANGLE RELATIONS CHECKED IN DELTAACVI7710
C     REFERENCES--ANGULAR MOMENTUM IN QUANTUM MECHANICS, A.R.EDMONDS,   ACVI7711
C                 PRINCETON                                             ACVI7712
C     ------------------------------------------------------------------ACVI7713
      IMPLICIT REAL*8(D)                                                ACVI7714
      DJHR3=0.D0                                                        ACVI7715
      ITMIN=MAX0(IABS(J1T-J9T),IABS(J2T-J6T),IABS(J4T-J8T))             ACVI7716
      ITMAX=MIN0(J1T+J9T,J2T+J6T,J4T+J8T)                               ACVI7717
      IF(ITMIN.GT.ITMAX)RETURN                                          ACVI7718
      DO 10 IT=ITMIN,ITMAX,2                                            ACVI7719
   10 DJHR3=DJHR3+DFLOAT(IT+1)*DRR3(J1T,J9T,J4T,J8T,IT,J7T)*            ACVI7720
     1DRR3(J2T,J6T,J8T,J4T,IT,J5T)*DRR3(J1T,J9T,J2T,J6T,IT,J3T)         ACVI7721
      DJHR3=DSQRT(DFLOAT((J3T+1)*(J6T+1)*(J7T+1)*(J8T+1)))*DJHR3        ACVI7722
      RETURN                                                            ACVI7723
      END                                                               ACVI7724
      FUNCTION D3JR3(J1T,J2T,J3T,M1T,M2T,M3T)                           ACVI7725
C     ------------------------------------------------------------------ACVI7726
C     3J COEFFICIENTS FOR R3--TRIANGLE RELATIONS CHECKED IN DELTA       ACVI7727
C     REFERENCES--ELEMENTARY THEORY OF ANGULAR MOMENTUM, M.E.ROSE, WILEYACVI7728
C     ------------------------------------------------------------------ACVI7729
      IMPLICIT REAL*8(D)                                                ACVI7730
      COMMON/BKDF/DLOGF(0:2000)                                         ACVI7731
      D3JR3=0.D0                                                        ACVI7732
      IF(M1T+M2T-M3T.NE.0)GOTO 20                                       ACVI7733
      DC=DELTA(J1T,J2T,J3T)                                             ACVI7734
      IF(DC.EQ.12345D0)GOTO 20                                          ACVI7735
      I1=J3T-J2T+M1T                                                    ACVI7736
      I2=J3T-J1T-M2T                                                    ACVI7737
      I3=J1T+J2T-J3T                                                    ACVI7738
      I4=J1T-M1T                                                        ACVI7739
      IF(BTEST(I4,0))GOTO 20                                            ACVI7740
      I5=J2T+M2T                                                        ACVI7741
      IF(BTEST(I5,0))GOTO 20                                            ACVI7742
      ITMIN=MAX0(0,-I1,-I2)                                             ACVI7743
      ITMAX=MIN0(I3,I4,I5)                                              ACVI7744
      IF(ITMIN.GT.ITMAX)GOTO 20                                         ACVI7745
      DTOP=(DC+DLOGF(J1T+M1T)+DLOGF(J1T-M1T)+                           ACVI7746
     1DLOGF(J2T+M2T)+DLOGF(J2T-M2T)+DLOGF(J3T+M3T)+                     ACVI7747
     2DLOGF(J3T-M3T))/DFLOAT(2)                                         ACVI7748
      DO 10 IT=ITMIN,ITMAX,2                                            ACVI7749
      DBOT=DLOGF(I3-IT)+DLOGF(I4-IT)+DLOGF(I5-IT)+                      ACVI7750
     1DLOGF(IT)+DLOGF(I1+IT)+DLOGF(I2+IT)                               ACVI7751
      DSUM=DEXP(DTOP-DBOT)                                              ACVI7752
      IF(BTEST(IT,1))THEN                                               ACVI7753
      D3JR3=D3JR3-DSUM                                                  ACVI7754
      ELSE                                                              ACVI7755
      D3JR3=D3JR3+DSUM                                                  ACVI7756
      ENDIF                                                             ACVI7757
   10 CONTINUE                                                          ACVI7758
      IF(BTEST(I1-I2,1))D3JR3=-D3JR3                                    ACVI7759
   20 RETURN                                                            ACVI7760
      END                                                               ACVI7761
      FUNCTION D6JR3(J1T,J2T,J3T,L1T,L2T,L3T)                           ACVI7762
C     ------------------------------------------------------------------ACVI7763
C     6J COEFFICIENTS FOR R3--TRIANGLE RELATION CHECKED IN DELTA        ACVI7764
C     REFERENCES--THE 3-J AND 6-J SYMBOLS, M.ROTENBERG, R.BIVINS,       ACVI7765
C                 N.METROPOLIS AND J.K.WOOTEN, MIT PRESS                ACVI7766
C     ------------------------------------------------------------------ACVI7767
      IMPLICIT REAL*8(D)                                                ACVI7768
      COMMON/BKDF/DLOGF(0:2000)                                         ACVI7769
      D6JR3=0.D0                                                        ACVI7770
      DX=DELTA(J1T,J2T,J3T)                                             ACVI7771
      IF(DX.EQ.12345D0)GOTO 20                                          ACVI7772
      DC=DX                                                             ACVI7773
      DX=DELTA(L1T,L2T,J3T)                                             ACVI7774
      IF(DX.EQ.12345D0)GOTO 20                                          ACVI7775
      DC=DX+DC                                                          ACVI7776
      DX=DELTA(L1T,J2T,L3T)                                             ACVI7777
      IF(DX.EQ.12345D0)GOTO 20                                          ACVI7778
      DC=DX+DC                                                          ACVI7779
      DX=DELTA(J1T,L2T,L3T)                                             ACVI7780
      IF(DX.EQ.12345D0)GOTO 20                                          ACVI7781
      DC=(DX+DC)/2.D0                                                   ACVI7782
      I1=J3T+L3T-J1T-L1T                                                ACVI7783
      I2=J3T+L3T-J2T-L2T                                                ACVI7784
      I3=J1T+J2T+L1T+L2T+2                                              ACVI7785
      I4=J1T+J2T-J3T                                                    ACVI7786
      I5=L1T+L2T-J3T                                                    ACVI7787
      I6=J1T+L2T-L3T                                                    ACVI7788
      I7=L1T+J2T-L3T                                                    ACVI7789
      ITMIN=MAX0(0,-I1,-I2)                                             ACVI7790
      ITMAX=MIN0(I3,I4,I5,I6,I7)                                        ACVI7791
      IF(ITMIN.GT.ITMAX)GOTO 20                                         ACVI7792
      DO 10 IT=ITMIN,ITMAX,2                                            ACVI7793
      DSUM=DEXP(DC+DLOGF(I3-IT)-(DLOGF(I4-IT)+DLOGF(I5-IT)+             ACVI7794
     1DLOGF(I6-IT)+DLOGF(I7-IT)+DLOGF(IT)+DLOGF(I1+IT)+DLOGF(I2+IT)))   ACVI7795
      IF(BTEST(IT,1))THEN                                               ACVI7796
      D6JR3=D6JR3-DSUM                                                  ACVI7797
      ELSE                                                              ACVI7798
      D6JR3=D6JR3+DSUM                                                  ACVI7799
      ENDIF                                                             ACVI7800
   10 CONTINUE                                                          ACVI7801
      IF(.NOT.BTEST(I3,1))D6JR3=-D6JR3                                  ACVI7802
   20 RETURN                                                            ACVI7803
      END                                                               ACVI7804
      FUNCTION D9JR3(J1T,J2T,J3T,J4T,J5T,J6T,J7T,J8T,J9T)               ACVI7805
C     ------------------------------------------------------------------ACVI7806
C     9J COEFFICIENTS FOR R3--TRIANGLE RELATIONS CHECKED IN DELTA       ACVI7807
C     REFERENCES--ANGULAR MOMENTUM, D.M.BRINK AND G.R.SATCHLER, OXFORD  ACVI7808
C     ------------------------------------------------------------------ACVI7809
      IMPLICIT REAL*8(D)                                                ACVI7810
      D9JR3=0.D0                                                        ACVI7811
      ITMIN=MAX0(IABS(J1T-J9T),IABS(J2T-J6T),IABS(J4T-J8T))             ACVI7812
      ITMAX=MIN0(J1T+J9T,J2T+J6T,J4T+J8T)                               ACVI7813
      IF(ITMIN.GT.ITMAX)RETURN                                          ACVI7814
      DO 10 IT=ITMIN,ITMAX,2                                            ACVI7815
   10 D9JR3=D9JR3+DFLOAT(IT+1)*DRR3(J1T,J9T,J4T,J8T,IT,J7T)*            ACVI7816
     1DRR3(J2T,J6T,J8T,J4T,IT,J5T)*DRR3(J1T,J9T,J2T,J6T,IT,J3T)         ACVI7817
      RETURN                                                            ACVI7818
      END                                                               ACVI7819
      FUNCTION DELTA(J1T,J2T,J3T)                                       ACVI7820
C     ------------------------------------------------------------------ACVI7821
C     DELTA FOR R3 ROUTINES--TRIANGLE RELATIONS CHECKED                 ACVI7822
C     ------------------------------------------------------------------ACVI7823
      IMPLICIT REAL*8(D)                                                ACVI7824
      COMMON/BKDF/DLOGF(0:2000)                                         ACVI7825
      DELTA=12345.D0                                                    ACVI7826
      I1=J1T+J2T-J3T                                                    ACVI7827
      IF(BTEST(I1,0))GOTO 10                                            ACVI7828
      IF(I1.LT.0)GOTO 10                                                ACVI7829
      I2=J2T+J3T-J1T                                                    ACVI7830
      IF(I2.LT.0)GOTO 10                                                ACVI7831
      I3=J3T+J1T-J2T                                                    ACVI7832
      IF(I3.LT.0)GOTO 10                                                ACVI7833
      DELTA=DLOGF(I1)+DLOGF(I2)+DLOGF(I3)-DLOGF(J1T+J2T+J3T+2)          ACVI7834
   10 RETURN                                                            ACVI7835
      END                                                               ACVI7836
C***************************                                            ACVI7837
       
