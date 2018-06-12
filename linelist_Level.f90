program linelist_trove
 implicit none

 integer isym,k,jf,ji,ilevel,ilevelf,ileveli,info,j0,i
 integer nfiles, jmax, ilines
 integer nlevels,itermi,itermf,J_t,imis,vi,vf

 real(8),allocatable :: energies(:)
 integer,allocatable :: J(:),NN(:),Nlines(:),numj(:),v(:)

 logical :: energeyfile_do = .false. 

 real(8) :: acoef,abscoef,energy,energyf,energyi,tranfreq,linestr,ZPE,coeff

 character(50) intfilename(400),enrfilename,enrfilenameout,intfilenameout,dscrfilename,filename
 character(1)  :: ch_t1,PQR
 character(200) :: ch_t200
 character(400) :: form
 integer        :: unit_f,ilevelmax,nmodes,nsym,nclasses
 character(len=3),allocatable   :: symm(:)
 integer   :: gns,DeltaV

 type numT
      !
      integer,pointer   :: nn(:)
      integer,pointer   :: sym(:)
      integer,pointer   :: J(:)
      real(8),pointer   :: energy(:)
      !
 end type numT

 type(numT),allocatable :: level(:,:)

 real(8),parameter :: planck=6.6260693d-27,avogno=6.0221415d+23,vellgt=2.99792458d+10,boltz=1.380658d-16

    !read number of intensity files
    read*,nfiles
    !
    if (nfiles>400) then 
      print('("Too many files (>400):",i)'),nfiles
      stop 'Too many files'
    endif
    !
    !read intensities filenames
    do i = 1,nfiles
     read*,intfilename(i)
    enddo
    !
    !read energies filename
    read*,enrfilename
    !
    !read intensity filename - out
     read*,intfilenameout
    !
    !read energies filename - out 
     read*,enrfilenameout

    !read nmodes, Nsym, nclasses
    read*,gns
    !
    ! read the rnage for DeltaV
    read*,deltaV
    !
    nsym = 2
    !

    !  
    !compute constants
    !ln2=log(2.0)
    !pi=acos(-1.0)
    !beta=planck*vellgt/(boltz*temp)
    !cmcoef=avogno/(8.0*pi*vellgt)
    !dpwcoef=sqrt(2.0*ln2*boltz*avogno)/vellgt
    !
    !read the energy file 
    open(unit=11,file=trim(enrfilename))
    !
    allocate(numj(0:1000),stat=info); 
    if (info/=0) stop 'error: numj is out of memory'
    !
    i = 0
    numj = 0
    jmax = 0
    !
    ! skip 5 lines from fort.8
    !
    do i=1,5
      read(11,"(a200)") ch_t200
    enddo
    !
    i = 0
    ZPE = 0
    !
    do
       !
       read(11,*,end=11) vi,j_t,energy
       !
       !print*,v,j,ilevel,energy
       !
       i = i + 1
       !
       numj(j_t) = numj(j_t) + 1
       !
       jmax = max(jmax,j_t)
       !
       if (j_t==0.and.vi==0) ZPE = energy
       !
       cycle
11     exit
    end do
    !
    jmax = jmax + 1
    !
    print("(' jmax = ',i4)"),jmax
    !
    nlevels = i
    print("(' nlevels = ',i)"),nlevels
    !
    rewind(11)
    !
    allocate(energies(nlevels),v(nlevels),J(nlevels),NN(0:Jmax),stat=info); 
    if (info/=0) stop 'error: energies,J is out of memory'
    !
    if ( all(trim(enrfilenameout)/=(/'none','NONE','null','NULL'/)) ) then  
      !
      energeyfile_do = .true.
      !
      open(unit=12,file=trim(enrfilenameout),action='write',status='replace')
      !
    endif
    !
    NN(:) = 0
    !
    J_t =  0
    !
    print("('Reading energies...')")
    !
    ilevelmax = 0
    !
    if (jmax>1000) stop 'Jmax>1000'
    !
    do i=1,5
      read(11,"(a200)") ch_t200
    enddo
    !
    do i = 1,nlevels
       !
       read(11,*) v(i),J(i),energies(i)
       !
       Ji     = J(i)
       !
       NN(Ji) = NN(Ji) + 1
       J_t = Ji
       !
       if (energeyfile_do) then
         write(12,"(i12,1x,f12.6,1x,i6,1x,i7,2x,i4)") i,energies(i)-ZPE,gns*(2*J(i)+1),j(i),v(i)
       endif
       !
    end do
    !
    !print("('Maximal number of levels in a block = ',i9)"), ilevelmax
    !
    !print("(' Energy ranges:')")
    !
    !do i = 0,jmax
    !    if ( NN(i)/=0 ) then
    !         print("(i4,2x,'[',i8,'...',i8,'] = ',i8)"),i,NN(i),NN(i)-NN(i)+1
    !      else
    !         print("(i4,2x,'[',i8,'...',i8,'] = ',i8,' tot = ',i8)"),i,NN(i),NN(i)-NN(i)+1,NN(i)-NN(i)+1
    !    endif
    !enddo
    !
    close(11)
    if (energeyfile_do) close(12)
    !
    !deallocate(nvib,lvib,tauvib,vvib,symvib)
    !
    print("(' N of levels =  ',i)"),nlevels
    !
    ! finish if the intensity output is not requested
    !
    if (trim(intfilenameout)=="none".or.trim(intfilenameout)=="NONE") stop 
    !
    print("('Generate the Transition file ...')")
    !
    open(unit=13,file=trim(intfilenameout),action='write',status='replace')
    !
    !start loop over all transitions
    !
    allocate(Nlines(0:jmax))
    !
    ilines = 0
    Nlines = 0
    !
    do i = 1,nfiles
      !
      open(unit=1,file=trim(intfilename(i)))
      !
      print("( 'process ',a,3x,i)"), intfilename(i),ilines
      !
      ilines = 0
      imis = 0
      !
      do k=1,11
        read(1,"(a200)") ch_t200
      enddo
      !
      do
         !   read new line from intensities file
         !
         read(1,"(1x,a1,1x,i3,3x,i3,2x,i3,f10.2,f14.6,e14.5)",end=20) PQR,ji,vf,vi,energyi,tranfreq,acoef
         !
         select case (PQR)
          case ("R") 
            Jf = Ji+1  
          case ("P") 
            Jf = Ji-1  
          case ("Q") 
            Jf = Ji  
          case default
            stop 'Wrong PQR'
         end select

!R(  0)    0 -  0    374.21     -0.605582   6.77594D-08   1.00000D+00   1.70838D+00
         !
         Nlines(max(Ji,Jf))= Nlines(max(Ji,Jf)) + 1
         !
         if (jf>jmax.or.ji>jmax) cycle
         !
         if (vf-vi > deltaV) cycle
         !
         itermi = 0
         !
         do k = 1,nlevels
           !
           if (ji ==J(k) .and. vi == v(k)) then 
             !
             if (abs(energyi-energies(k))>5e-2) then 
               write(6,"('Wrong energy: ',2i,2f13.6)") Ji,vi,energyi,energies(k)
             endif
             !
             imis = imis + 1
             !
             itermi = k
             !
             exit
             !
           endif
           !
         enddo
         !
         itermf = 0
         !
         do k = 1,nlevels
           !
           if (jf ==J(k) .and. vf == v(k)) then 
             !
             if (abs(energyi-tranfreq-energies(k))>5e-2) then 
               write(6,"('Wrong energy: ',2i,2f13.6)") Jf,vf,energyf,energies(k)
               stop 'Wrong energy found'
             endif
             !
             imis = imis + 1
             !
             itermf = k
             !
             exit
             !
           endif
           !
         enddo
         !
         if (itermi==0.or.itermf==0) then 
           !
           print("( 'no match for ',a1,3i3,2f14.6 )"), PQR,ji,vf,vi,energyi,tranfreq
           !stop 'cannot find the correspondence in the energy-file'
           !
           itermf = -1 ; itermi = -1
           !
           cycle
           !
         endif
         !
         tranfreq = energies(itermf)-energies(itermi)
         !
         unit_f = 13
         !
         ilines = ilines + 1
         !
         write(unit_f,"(2i12,1x,es10.4,1x,f15.6)"),itermf,itermi,acoef,tranfreq
         !
         cycle
      20  exit
      enddo
      !
      if (imis /= 0) print("( 'missassigned idescr = ',i)"), imis
      !
    enddo 
    !
    do i = 0,jmax
     print("(i4,2x,i12)"), i,Nlines(i)
    enddo
    !
    deallocate(energies,J,v)
    !
    print("( 'Done!')")
    !
    close(1)
    close(13)


end program linelist_trove
