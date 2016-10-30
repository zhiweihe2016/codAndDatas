ccccccc....this program is to find preimages of any state.....ccccc
      program main
	implicit none
	integer i,j,k
	integer n
	parameter (n=11)
      real*8 matrix(n,n)
      integer gm(n,n),rm(n,n)

	integer tf,tfckeck
	integer state(n)
	integer staold(n),stanew(n)

	integer ncount,newcount
	integer last
	integer prefinalnew(1000,n),prefinalold(1000,n)


      open(2,file="preimage.dat",status="unknown")

	open(1,file="Budding_yeast11.dat",status="old")	
      do i=1,n    ! input the network
         read(1,100)matrix(i,1:n)
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
      
	state=0
	do i=1,n
	   state(i)=0
	   staold(i)=-100
	enddo

	state(2)=1
	state(5)=1
	state(9)=1


      tf=1
	do while(tf.eq.1)
	   call rule(stanew,staold,state,gm,rm,n,tf)
	   if(tf.eq.0)then
	     stanew(1:n)=100
	     goto 101
	   endif
	      
         tfckeck=0
	   do i=1,n
	      if (staold(i).ne.stanew(i))then
	         tfckeck=tfckeck+1
	      endif
	   enddo

	   if(tfckeck.eq.0)then
	      tf=0
	   endif		       
	   staold=stanew
	enddo

c	write(*,*)stanew

101	prefinalold(1,1:n)=stanew(1:n)
      ncount=1
      last=0
	do while(last.eq.0)
	   call prefinal(prefinalnew,newcount,prefinalold,
     $           ncount,state,gm,rm,n)
	   write(*,*)newcount

	   last=1
	   do i=1,newcount
	      do j=1,n
	         if (prefinalnew(i,j)<-50)then
	            last=0
	            goto 102
	         endif
	     enddo
	   enddo

102      prefinalold=prefinalnew
         ncount=newcount
	enddo

	do i=1,newcount
	   write(2,100)real(prefinalnew(i,1:n))
	enddo

	do i=1,newcount
	   write(*,*)prefinalnew(i,1:n)
	enddo

	
100   format(1X,100f16.8)
      end program main


      subroutine rule(stanew,staold,state,gm,rm,n,tf)
      implicit none
      integer i,j,k
      integer n
      integer gm(n,n),rm(n,n)
	integer state(n)
      
	integer tfp,tf
	integer logi,panduan,nodeone
	integer nodepeak
	integer stanew(n),staold(n)
      
	stanew=staold
      tfp=1
	do i=1,n      
	   !rule 1：if S_i(t+1)=1
	   if (state(i).eq.1)then  

	      !rule 1.1: $r_{ij}=1$, then $S_j=0$;
	      do j=1,n
	         if ((i.ne.j).and.(rm(i,j)>0))then  
	             if(stanew(j).eq.1)then
	               tfp=0
	               goto 1001
	             else
	               stanew(j)=0
	             endif
	         endif
	      enddo
	   endif
	enddo

	do i=1,n
         if (state(i).eq.1)then  
	      panduan=0
	      nodeone=0
	      do j=1,n
	          if((j.ne.i).and.(stanew(j)*gm(i,j).ne.0))then
	             panduan=panduan+1
	              nodeone=j
	          endif
	      enddo
         
            !rule 1.2: $r_{ii}=1$, $g_{ij_0}=1$, $\sum_{j\ne i} ({S_j}r_{ij})=0$ and $\sum_{j\ne j_0}({S_j}g_{ij})=0$, then $S_{j_0}=1$;
	       if ((panduan.eq.1).and.(rm(i,i)>0))then
	          if(stanew(nodeone).eq.0)then
	             tfp=0
	             goto 1001
	           else
	             stanew(nodeone)=1
	           endif
	        endif

	       !rule 1.4: if $g_{ii}=0$, $g_{ij_0}=1$, $\sum_{j\ne i} ({S_j}r_{ij})=0$, $\sum_{j\ne j_0}({S_j}g_{ij})=0$ and $S_i=0$, then $S_{j_0}=1$;
	  if ((panduan.eq.1).and.(gm(i,i).eq.0).and.(stanew(i).eq.0))then
	          if(stanew(nodeone).eq.0)then
	             tfp=0
	             goto 1001
	           else
	             stanew(nodeone)=1
	           endif
	        endif

	        

             !rule 1.3: if $g_{ii}=0$, $r_{ii}=0$, $\sum_{j\ne i} ({S_j}r_{ij})=0$ and $\sum_{j\ne i} ({S_j}g_{ij})=0$, then $S_i=1$;
             if ((panduan.eq.0).and.(gm(i,i)*rm(i,i).eq.0))then
	             if(stanew(i).eq.0)then
	               tfp=0
	               goto 1001
	             else
	               stanew(i)=1
	             endif
	       endif

             
             !规则1.5: if $r_{ii}=1$ and $\sum_{j\ne i}({S_j}g_{ij})=0$, then the current state $S(t)$ is not a pre-image.
             if ((panduan.eq.0).and.(rm(i,i)>0))then
	           tfp=0
	           goto 1001
	       endif
		 				       
	   endif   !endif_state1
	enddo


	   !rule 2：若S_i(t+1)=0，且i不受有效抑制作用:  规则2.1,若j促进i，则S_j(t)=0;
	   !                                           规则2.2，若g_ii=0且r_ii=0,则S_i(t)=0;
	   !                                           规则2.3，若g_ii=1，则矛盾
	   !       若S_i(t+1)=0，且只有一个S_j(t)*rm(i,j).ne.0: 规则2.4，若对i的有效促进作用非零，或是i有自促进作用，则S_j(t)=1;   
	   
	   


	 !rule 2.1：if $\sum_{j\ne i} ({S_j}r_{ij})=0$ and there is a node $j_0$ such that $g_{ij_0}=1$, then $S_{j_0}=0$;
	 do i=1,n                              
         if (state(i).eq.0)then  !if_state2.1
	      logi=1
	      do j=1,n
	         if (j.ne.i)then
	            logi=logi.and.(stanew(j)*rm(i,j).eq.0)
	         endif
	      enddo

	      if(logi.eq.1)then   !endif_logi2.1 

	        do j=1,n
	           if((j.ne.i).and.(gm(i,j)>0))then        
	             if(stanew(j).eq.1)then
	               tfp=0

	               goto 1001
	             else
	               stanew(j)=0
	             endif
	           endif
	        enddo

	      endif    !endif_logi2.1

	    endif    !endif_state2.1
	enddo



       !rule 2.2-2.5      
	do i=1,n
	   if (state(i).eq.0)then  !if_state2.2
	      panduan=0
	      nodepeak=0
	      do j=1,n
	         if ((j.ne.i).and.(stanew(j)*rm(i,j).ne.0))then
	            panduan=panduan+1
	            nodepeak=j
	         endif
	      enddo

	      if(panduan.eq.0)then   !if_logi2.2
               ! rule 2.2: if $g_{ii}=0$, $r_{ii}=0$ and $\sum_{j\ne i} ({S_j}r_{ij})=0$, then $S_i=0$;
	        if((rm(i,i).eq.0).and.(gm(i,i).eq.0))then   
	           if(stanew(i).eq.1)then            
	               tfp=0
	               goto 1001
	           else
	               stanew(i)=0
	           endif
	        endif

    
               !rule 2.5：if $g_{ii}=1$ and $\sum_{j\ne i} ({S_j}r_{ij})=0$, then the current state $S(t)$ is not a pre-image;
              if(gm(i,i)>0)then        
	           tfp=0
	            goto 1001
	        endif

	      endif    !endif_logi2.2


            !rule 2.3 and 2.4：if $r_{ij_0}=1$, and $\sum_{j\ne j_0} ({S_j}r_{ij})=0$,and $g_{ii}=1$ or there is a node $j_0^{'}$ such that $S_jg_{ij_0^{'}}=1$, then $S_{j_0}=1$;
            if(panduan.eq.1)then   !if_logi2.4
	        logi=0
	        do j=1,n
	           if ((j.ne.i).and.(stanew(j)*gm(i,j)>0))then
	              logi=logi+1
	           endif
	        enddo
              logi=logi+gm(i,i)

	        if (logi>0)then
                 if(stanew(nodepeak).eq.0)then            
	               tfp=0
	               goto 1001
	           else
	               stanew(nodepeak)=1
	           endif
	        endif

	      endif !endif_logi2.4
	        
	   endif    !endif_state2.2

	enddo

	   

1001	tf=tfp            	        	       	          
      end subroutine rule

      subroutine step(pre2,staoldinitial,state,peakcheck,gm,rm,n)
	implicit none
	integer i,j,k,ii
      integer n
      integer gm(n,n),rm(n,n)
	integer peakcheck

	integer state(n)
	integer staoldinitial(n),staold(n),stanew(n)
      
	integer tf,tch
	integer pre2(2,n)

      pre2=0
	do ii=1,2
	   staold(1:n)=staoldinitial(1:n)
	   if(ii.eq.1)then
	     staold(peakcheck)=1
	   else
	     staold(peakcheck)=0
	   endif
        
	   tf=1
	   do while(tf.eq.1)
	      stanew=0
	      call rule(stanew,staold,state,gm,rm,n,tf)

c          if (ii.eq.1)write(*,*)stanew

            if(tf.eq.0)then
	         pre2(ii,1:n)=100
	      else
	         tch=0
	         do i=1,n
			    if(staold(i).ne.stanew(i))then
	               tch=tch+1
	            endif
	         enddo

	         if(tch.eq.0)then
	           tf=0
	           pre2(ii,1:n)=stanew(1:n)
	         endif
	      endif

	      staold=stanew
	   enddo
	enddo
      
	endsubroutine step
      

	subroutine prefinal(prefinalnew,newcount,prefinalold,
     $           ncount,state,gm,rm,n)
	implicit none
	integer i,j,k
	integer n,ncount,newcount
	integer gm(n,n),rm(n,n)
	integer state(n)
	integer prefinalold(1000,n)
      
	integer staoldinitial(n)
	integer peakcheck
	integer pre2(2,n)

	integer negcount
	integer prefinalnew(1000,n)

	prefinalnew=0
	newcount=0
	do i=1,ncount
	   staoldinitial(1:n)=prefinalold(i,1:n)
	   negcount=0
         do j=1,n
	      if(staoldinitial(j)<-10)then
	        negcount=negcount+1
              peakcheck=j
	        goto 2001
	      endif
	   enddo

2001	   if (negcount.eq.0)then
            newcount=newcount+1
            prefinalnew(newcount,1:n)=staoldinitial(1:n)
	   else
	      call step(pre2,staoldinitial,state,peakcheck,gm,rm,n)

            do k=1,2
c	         write(*,*)k*10,pre2(k,1:n)

	         if(pre2(k,1)<50)then
	           newcount=newcount+1
                 prefinalnew(newcount,1:n)=pre2(k,1:n)
	         endif
	      enddo
	   endif
	enddo

	end subroutine prefinal


	        
