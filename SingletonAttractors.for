ccccccc....this program is to find the single attractors of a network.....ccccc
      program main
      implicit none
      integer i,j,k,ii,jj,iii,kk,kk1
      integer n
      parameter(n=23)

      real*8 matrix(n,n)
      integer gm(n,n),rm(n,n)
	integer newcount,ncount
      
	integer peak(n)
	integer last,time

	integer,allocatable:: bfold(:,:),bfnew(:,:)
      integer,allocatable:: bfnew1(:,:),bfnew2(:,:)
      
	integer newcount1,ncount1,ncount2
	integer dtraj(n)
	integer kst
	parameter (kst=1000)
      
	integer nbas
	real*8 t1,t2
      
	open(2,file="cputime.dat",status="unknown")
	open(3,file="basin.dat",status="unknown")

      
	open(4,file="Thelper_cell23.dat",status="old")   
      do i=1,n    !Input the network matrix
	   read(4,100)matrix(i,1:n)
	enddo
      
         gm=0
         rm=0
         do i=1,n
            do j=1,n
               if(matrix(i,j)>0.1)then
                 gm(i,j)=1
               else if(matrix(i,j)<-0.1)then
                 rm(i,j)=1
               endif
           enddo
        enddo

	  call cpu_time(t1)

  
        call checkrm(peak,rm,n)

        ncount=1
	  ALLOCATE(bfnew1(ncount,n))
	  do i=1,n
	     bfnew1(1,i)=-10
	  enddo

	  last=0
	  time=0
	  do while(last.eq.0)
	       time=time+1
	       ALLOCATE(bfold(ncount,n),bfnew(2*ncount,n))
	       bfold(1:ncount,1:n)=bfnew1(1:ncount,1:n)

             DEALLOCATE (bfnew1)

             call bfinal(bfnew,newcount,bfold,ncount,peak,gm,rm,n)


             ALLOCATE(bfnew1(newcount,n))
	       bfnew1=0
	       last=1
	       newcount1=0
	       do i=1,newcount
	          ncount1=0
	          do j=1,n
	             if (bfnew(i,j)<-5)then
			       last=0
	               ncount1=ncount1+1
	               goto 110
			     endif
	          enddo

	         if(ncount1.eq.0)then
	           nbas=nbas+1
	           dtraj(1:n)=bfnew(i,1:n)
	           write(3,100)real(dtraj(1:n))
	         endif
            
110	         if(ncount1.ge.1)then
	           newcount1=newcount1+1
                 bfnew1(newcount1,1:n)=bfnew(i,1:n)
	         endif
	      enddo

	   
c           write(*,*)time,newcount1
	      ncount=newcount1

	      DEALLOCATE(bfold,bfnew)

	      if (ncount>500000)goto 111	       
	  enddo

	  if (ncount.le.500000)then
	     DEALLOCATE (bfnew1)
	     goto 130
	  endif

c	  write(*,*)bfnew1(1,1:n)
c       write(*,*)ncount
      
111     ncount2=ncount
        ALLOCATE(bfnew2(ncount,n))
        bfnew2(1:ncount2,1:n)=bfnew1(1:ncount2,1:n)
        DEALLOCATE (bfnew1)
        kk=0
	  do while(kk<ncount2)

           if(mod(kk,10000).eq.0)write(*,*)kk
	     kk1=kk
	     kk=kk+kst
	     ncount=kst

	     if(kk.ge.ncount2)then
	        kk=ncount2
	        ncount=ncount2-kk1
	     endif

	     ALLOCATE(bfnew1(ncount,n))
	     bfnew1(1:ncount,1:n)=bfnew2((kk1+1):kk,1:n)

	     last=0
	     time=0
	     do while(last.eq.0)
	        time=time+1
	        ALLOCATE(bfold(ncount,n),bfnew(2*ncount,n))
	        bfold(1:ncount,1:n)=bfnew1(1:ncount,1:n)

              DEALLOCATE (bfnew1)

              call bfinal(bfnew,newcount,bfold,ncount,peak,gm,rm,n)

c	         write(*,*)time,newcount

              ALLOCATE(bfnew1(newcount,n))
              bfnew1=0

	        last=1
	        newcount1=0
	        do i=1,newcount
	            ncount1=0
	           do j=1,n
	              if (bfnew(i,j)<-5)then
			         last=0
	                 ncount1=ncount1+1
	                 goto 120
			      endif
	           enddo

	           if(ncount1.eq.0)then
	              nbas=nbas+1
	              dtraj(1:n)=bfnew(i,1:n)
	             write(3,100)real(dtraj(1:n))
	           endif
            
120	           if(ncount1.ge.1)then
	             newcount1=newcount1+1
                   bfnew1(newcount1,1:n)=bfnew(i,1:n)
	           endif
	       enddo

           
	       ncount=newcount1

	       DEALLOCATE(bfold,bfnew)

	     enddo

           DEALLOCATE (bfnew1)
	  enddo

        DEALLOCATE (bfnew2)

130	  call cpu_time(t2)

	  write(2,100)t2-t1,real(nbas)
	  write(*,*)t2-t1,nbas


100   format(1X,999f16.8)
      end program main

ccccccccccccccccccc.....the rules 1-7.........ccccccccccccccccccccccccccccccccccc
      subroutine rule(basnew,basold,gm,rm,n,tf)
      implicit none
      integer i,j,k
      integer n
      integer basold(n),basnew(n)
      integer gm(n,n),rm(n,n)
      
	integer tfp,tf
	integer logi,pan
      
	tfp=1
      basnew=basold
      do i=1,n
         !rule 1£ºif S_i=1 and r_{ji}=1, then S_j=0;
         do j=1,n
            if((j.ne.i).and.(basold(j)*rm(i,j).eq.1))then
	         if (basnew(i).eq.1)then
	            tfp=0
	            goto 10
	         else
                  basnew(i)=0
	         endif
            endif
         enddo

        if(basold(i)>0)then
        !rule 2£ºif S_i=1 and r_{ij}=1, then S_j=0;
           do j=1,n
              if(j.ne.i.and.(rm(i,j)>0))then
	           if (basnew(j).eq.1)then
	              tfp=0
	              goto 10
	           else
                    basnew(j)=0
	           endif
              endif
           enddo

        !rule 3£º if S_i=1, r_{ii}=1, g_{ij_0}=1 and \sum_{j\ne j_0}({S_j}g_{ij})=0, then S_{j_0}=1;
           if(rm(i,i)>0)then
             logi=0
             do j=1,n
                if(j.ne.i)then
                  logi=logi+basold(j)*gm(i,j)

                  if (basold(j)*gm(i,j)>0)then
                    k=j
                  endif
                endif
             enddo

             if(logi.eq.1)then
	          if (basnew(k).eq.0)then
	            tfp=0
	            goto 10
	          else
                  basnew(k)=1
	          endif
             endif
           endif

        endif
      
         !rule 5£ºif r_{ii}=1 and \sum_{j\ne i ({S_j}g_{ij})=0, then S_i=0;
         if (rm(i,i)>0)then
            pan=1
            do j=1,n
               if (j.ne.i)then
                  pan=pan.and.(gm(i,j)*basold(j).eq.0)
               endif
            enddo

           if(pan>0)then
	        if (basnew(i).eq.1)then
	            tfp=0
	            goto 10
	        else
                  basnew(i)=0
	        endif
           endif
         endif

           pan=1      
           do j=1,n
              if (j.ne.i)then
                 pan=pan.and.(rm(i,j)*basold(j).eq.0)
              endif
           enddo

        !rule 6£º if g_{ii}=1 and \sum_{j\ne i} ({S_j}r_{ij})=0, then S_i=1;

	  !rule 7£ºif \sum_{j\ne i} ({S_j}r_{ij})=0, and there is a node j_0 such that S_{j_0}g_{ij_0}=1, then S_i=1.

	   logi=0
	   do j=1,n
	      if((j.ne.i).and.(basold(j)*gm(i,j).eq.1))then
	            logi=logi+1
	      endif
	   enddo
         logi=logi+gm(i,i)

	   if((logi>0).and.(pan>0))then
	      if (basnew(i).eq.0)then
	          tfp=0
	          goto 10
	      else
                basnew(i)=1
	      endif
	   endif

	   !rule 4£ºif S_i=0 and \sum_{j\ne i} ({S_j}r_{ij})=0, g_{ij_0}=1 and \sum_{j\ne j_0}({S_j}g_{ij})=0, then S_{j_0}=0;
	   if((basold(i).eq.0).and.(pan>0))then
	       do j=1,n
	          if((j.ne.i).and.(gm(i,j)>0))then
	             if(basnew(j).eq.1)then
	                tfp=0
	                goto 10
	             else
	                basnew(j)=0
	             endif
	          endif
	       enddo
	    endif
	           
      enddo
      
10    tf=tfp
      end subroutine rule


cccccccccc......step 2...cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine checkrm(peak,rm,n)
	implicit none
	integer i,j,k
	integer n
	integer rm(n,n)
      
	integer guodu,igd
      integer check(n)
	integer peak(n)
      
	do i=1,n
	   check(i)=0
	   do j=1,n
	      if((j.ne.i).and.(rm(j,i)>0))then
		     check(i)=check(i)+1
		  endif

	      if((j.ne.i).and.(rm(i,j)>0))then
		     check(i)=check(i)+1
		  endif
	   enddo
	enddo
      
      peak=check

	end subroutine checkrm


	subroutine posi(posit,peak,numi,n)
	implicit none
	integer i
      integer n
	integer numi
	integer peak(n)
	integer posit
      
	posit=0
	do i=1,n
	   if(numi.eq.peak(i))then
	     posit=i
	   endif
	enddo

	end subroutine posi





ccccccc......step 3.......ccccccccccccccccccccccccccccccccccc
	subroutine step(basin2,basoldintial,peakcheck,gm,rm,n)
	implicit none
	integer i,j,k,ii
      integer n
	integer gm(n,n),rm(n,n)
	integer peakcheck
	integer basoldintial(n),basold(n),basnew(n)

	integer tf,tch
	integer basin2(2,n)

      basin2=0
	do ii=1,2
         basold=basoldintial
	   if(ii.eq.1)then
	     basold(peakcheck)=1
	   else
	     basold(peakcheck)=0
	   endif
	
      tf=1
	do while(tf.eq.1)
	   basnew=0
	   call rule(basnew,basold,gm,rm,n,tf)
      
	   if(tf.eq.0)then
	      basin2(ii,1:n)=10
	   else
	      tch=0
	      do i=1,n
		     if(basold(i).ne.basnew(i))then
			   tch=tch+1
	         endif
	      enddo

	      if(tch.eq.0)then
	         tf=0
               basin2(ii,1:n)=basnew(1:n)
	      endif
	   endif
	        
	   basold=basnew
	enddo
	enddo
      
	endsubroutine step
      
	subroutine bfinal(bfnew,newcount,bfold,ncount,peak,gm,rm,n)
	implicit none
	integer i,j,k
	integer n,ncount
	integer gm(n,n),rm(n,n)
	integer peak(n)
	integer bfold(ncount,n)

	integer basin2(2,n),basoldintial(n)
	integer negcount
	integer positold
	integer peakcheck
      
	integer bfnew(2*ncount,n)
	integer newcount
      
	bfnew=0
	newcount=0
      do i=1,ncount
         basoldintial(1:n)=bfold(i,1:n)

	   positold=-1
         peakcheck=0
	   do j=1,n
	      if((basoldintial(j)<-5).and.(positold<peak(j)))then
	         positold=peak(j)
	         peakcheck=j
	      endif	      
	   enddo


	    call step(basin2,basoldintial,peakcheck,gm,rm,n)
	      do k=1,2
		     if(basin2(k,1)<5)then
                  newcount=newcount+1
                  bfnew(newcount,1:n)=basin2(k,1:n)
	         endif
	      enddo
	enddo

      end subroutine bfinal


