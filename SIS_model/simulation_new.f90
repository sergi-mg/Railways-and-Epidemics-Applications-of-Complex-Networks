!Author: Sergi Martínez Galindo
program simulations_network2
    implicit none
    integer :: endemic, i, ios
    double precision :: final_t
    character*64 :: file_name, file_name_s
    !parametres
    character*32 sN,sgamma,slambda, sseed, sE, sN_sim
    double precision :: gamma, lambda
    integer :: N, seed, E, N_sim

    integer, dimension(:), allocatable :: P, D, D_sim
    common/parameters/lambda 

    !Parameters reading

    !Order to execute the program:
    !.\name.exe N gamma lambda seed E
    !double precision (2,3), integer (1,4,5,6)

    call getarg(1 , sN)
    call getarg(2 , sgamma)
    call getarg(3 , slambda)
    call getarg(4, sseed)
    call getarg(5, sE)
    call getarg(6, sN_sim)


    READ (sN,*) N
    READ (sgamma,*) gamma
    READ (slambda,*) lambda
    READ (sseed, *) seed
    READ (sE, *) E
    READ (sN_sim, *) N_sim

    !write(*,*)lambda
    write(*,*)"Start program.",lambda,gamma,N,seed

    !create the file for (N,lambda)
    
    file_name_s="../data_results_fortran/results_simulation.dat"
    open(15, file=trim(file_name_s))

    !read the network
    allocate(P(1:N+1), D(1:2*E), D_sim(1:2*E))
    !number of neighbours
    write(file_name,"(A,I0,A)")"../data_networks/P",seed,".dat"
    open(11, file=trim(file_name))
    do i = 1, N+1
        read(11, *, iostat=ios) P(i)
        if (ios /= 0) then
            write(*,*) "Error reading P at i=", i
            stop
        end if
    end do
    close(11)

    !neighbours
    write(file_name,"(A,I0,A)")"../data_networks/D",seed,".dat"
    open(12, file=trim(file_name))
    do i = 1, 2*E
        read(12, *, iostat=ios) D(i)
        if (ios /= 0) then
            write(*,*) "Error reading P at i=", i
            stop
        end if
    end do
    close(12)

    !indexes on D symmetric
    write(file_name,"(A,I0,A)")"../data_networks/D_sim",seed,".dat"
    open(13, file=trim(file_name))
    do i = 1, 2*E
        read(13, *, iostat=ios) D_sim(i)
        if (ios /= 0) then
            write(*,*) "Error reading P at i=", i
            stop
        end if
    end do
    close(13)

    !write(*,*)"Start simulation loop"
    do i=1,N_sim
        !simulation
        !write(*,*)"New simulation",i
        call init_genrand(i)
        call gillespie_SIS(N,E,P,D,D_sim,final_t,endemic)
        !save the results
        write(15,"(e25.15,3X,I0)")final_t, endemic
    enddo
    !write(*,*) "End of simulations"

    close(15)
    write(*,*)"End of program reached.",lambda,seed,i

end program

subroutine add_link_i(node_s_i,node_susceptible,index_j,E,D_sim,D_act,active,Na)
    implicit none 
    !inputs
    integer :: node_s_i,node_susceptible,E,Na,index_j
    integer :: D_sim(1:2*E),D_act(1:2*E),active(1:E,1:2,1:2)
    !subrotuine
    integer :: index_k

    Na=Na+1
    !add the index of the active link to the s to i neighbourhood 
    D_act(index_j)=Na 
    !add the index of the active link to the susceptible node 
    index_k=D_sim(index_j)
    D_act(index_k)=Na
    !add the active link at the end of active (nodes)
    active(Na,1,1)=node_s_i
    active(Na,2,1)=node_susceptible
    !add the active link at the end of active (D_index)
    active(Na,1,2)=index_j
    active(Na,2,2)=index_k

    return
end subroutine

subroutine add_link_r(node_infected,node_i_s,index_j,E,D_sim,D_act,active,Na)
    implicit none 
    !inputs
    integer :: node_infected,node_i_s,E,Na,index_j
    integer :: D_sim(1:2*E),D_act(1:2*E),active(1:E,1:2,1:2)
    !subrotuine
    integer :: index_k

    Na=Na+1
    !add the index of the active link to the i to s neighbourhood 
    D_act(index_j)=Na 
    !add the index of the active link to the infected node
    index_k=D_sim(index_j)
    D_act(index_k)=Na
    !add the link at the end of active (nodes)
    active(Na,1,1)=node_infected
    active(Na,2,1)=node_i_s
    !add the link at the end of active (D_index)
    active(Na,1,2)=index_k
    active(Na,2,2)=index_j

    return
end subroutine

subroutine remove_link(link,index_j,E,D_sim,D_act,active,Na)
    implicit none 
    !inputs
    integer :: E,Na,index_j,link
    integer :: D_sim(1:2*E),D_act(1:2*E),active(1:E,1:2,1:2)
    !subrotuine
    integer :: node_1,node_2,k_1,k_2,index_k

    !move the last link to the position of the futurely deleted link 
    active(link,:,1)=active(Na,:,1)
    active(link,:,2)=active(Na,:,2)
    !read nodes from link which changes position in active
    node_1=active(link,1,1)
    node_2=active(link,2,1)
    !read D_index from link which changes position in active
    k_1=active(link,1,2)
    k_2=active(link,2,2)

    !update node 1
    D_act(k_1)=link

    !update node 2
    D_act(k_2)=link

    !delete the last link in active
    active(Na,:,1)=[0,0]
    active(Na,:,2)=[0,0]
    Na=Na-1

    !delete the index in the modified node neighbourhood (deleted link, input)
    D_act(index_j)=0

    !delete the index in the neighbour node neighbourhood (deleted link)
    index_k=D_sim(index_j)
    D_act(index_k)=0

    return
end subroutine


subroutine recovery(index_i_s,N,E,p,D,D_sim,D_act,active,infected,Na,Ni)
    implicit none 
    !inputs
    integer :: N,E,Na,Ni,index_i_s
    integer :: p(1:N+1),D(1:2*E),D_sim(1:2*E),D_act(1:2*E),active(1:E,1:2,1:2),infected(1:N)
    !subrotuine
    integer :: j, node_i_s, j_i, j_f, link_j, n_node

    !node modified
    node_i_s=infected(index_i_s)

    !infected vector update
    infected(index_i_s)=infected(Ni)
    infected(Ni)=0
    Ni=Ni-1


    !delete/add active links

    j_i=p(node_i_s)
    j_f=p(node_i_s+1)-1

    do j=j_i,j_f
        !delete links
        link_j=D_act(j)
        n_node=D(j)
        if (link_j.ne.0) then   
            call remove_link(link_j,j,E,D_sim,D_act,active,Na)
        !add links
        else
            call add_link_r(n_node,node_i_s,j,E,D_sim,D_act,active,Na)
        endif
    enddo

    return
end subroutine

subroutine infection(link_i,N,E,p,D,D_sim,D_act,active,infected,infected_end,Na,Ni,N_end)
    implicit none 
    !inputs
    integer :: link_i,N,E,Na,Ni,N_end
    integer :: p(1:N+1),D(1:2*E),D_sim(1:2*E),D_act(1:2*E),active(1:E,1:2,1:2),infected(1:N),infected_end(1:N)
    !subrotuine
    integer :: j, node_s_i, j_i, j_f, link_j, n_node

    !node modified
    node_s_i=active(link_i,2,1)

    !infected vector update
    Ni=Ni+1
    infected(Ni)=node_s_i
    if (infected_end(node_s_i).eq.0) then
        infected_end(node_s_i)=1
        N_end=N_end+1
    endif


    !delete/add active links
    j_i=p(node_s_i)
    j_f=p(node_s_i+1)-1

    do j=j_i,j_f
        !delete links
        link_j=D_act(j)
        n_node=D(j)
        if (link_j.ne.0) then   
            call remove_link(link_j,j,E,D_sim,D_act,active,Na)
        !add links
        else
            call add_link_i(node_s_i,n_node,j,E,D_sim,D_act,active,Na)
        endif
    enddo


    return
end subroutine

subroutine gillespie_step(p,D,D_sim,D_act,active,infected,Ni,Na,E,N,tau,infected_end,N_end)
    implicit none 
    !input
    integer :: N,E,Ni,Na,N_end
    integer :: p(1:N+1),D(1:2*E),D_sim(1:2*E),D_act(1:2*E),active(1:E,1:2,1:2),infected(1:N),infected_end(1:N)
    !subroutine
    double precision :: genrand_real3, tau, u1, u2, r2, a_i, a_r, lambda, a0, u3
    integer :: index_i_s, link_s_i

    common/parameters/lambda

    !1. calculate a_i and a_r
    a_i=Na*lambda
    a_r=Ni

    !2. calculate a0
    a0=a_i+a_r  

    !3. time?
    u1=genrand_real3()
    tau=-log(u1)/a0

    !4. reaction?
    u2=genrand_real3()
    r2=(a_i+a_r)*u2

    if (r2.lt.a_i) then
        !infection
        !5. which link?
        u3=genrand_real3()
        link_s_i=ceiling(u3*Na)
        call infection(link_s_i,N,E,p,D,D_sim,D_act,active,infected,infected_end,Na,Ni,N_end)
    else
        !recovery
        !5. which node?
        u3=genrand_real3()
        index_i_s=ceiling(u3*Ni)
        call recovery(index_i_s,N,E,p,D,D_sim,D_act,active,infected,Na,Ni,N_end)
    endif 

    return
end subroutine


subroutine gillespie_SIS(N,E,p,D,D_sim,final_t,endemic)
    implicit none 
    !input
    integer :: N,E
    integer :: p(1:N+1),D(1:2*E),D_sim(1:2*E)
    !subroutine
    integer :: D_act(1:2*E),active(1:E,1:2,1:2),infected(1:N),infected_end(1:N)
    integer :: Na,Ni,node_selected,k_i,k_f,k,index_l,i,N_end,neighbour_node,k_ns
    double precision :: genrand_real3,u,t,tau
    !output
    double precision :: final_t
    integer :: endemic

    !1. initialize supplementary vectors and matrices
    D_act=0
    active=0
    infected=0
    infected_end=0
    endemic=0
    Na=0
    Ni=0
    N_end=0

    !1. select a random initial infected node
    k_ns=3 
    do while (k_ns.lt.4) 
        u=genrand_real3()
        node_selected=ceiling(u*N)
        k_i=p(node_selected)
        k_f=p(node_selected+1)
        k_ns=k_f-k_i
    enddo


    !add to infected
    infected(1)=node_selected
    infected_end(node_selected)=1
    Ni=Ni+1
    N_end=N_end+1

    !add all its active links
    k_i=p(node_selected)
    k_f=p(node_selected+1)-1

    do k=k_i,k_f
        !neighbour
        neighbour_node=D(k)
        Na=Na+1

        !selected's neighbourhood
        D_act(k)=Na

        !neighbour's neighbourhood
        index_l=D_sim(k)
        D_act(index_l)=Na

        !new active link (nodes)
        active(Na,:,1)=[node_selected,neighbour_node]

        !new active link (D_index)
        active(Na,:,2)=[k, index_l]
    enddo

    !now we start with the time loop
    !set initial time to 0
    t=0
    do i=1,10**9
        call gillespie_step(p,D,D_sim,D_act,active,infected,Ni,Na,E,N,tau,infected_end,N_end)
        t=t+tau
        if (Ni.eq.0) then
            exit
        elseif (N_end.gt.N/2) then
            endemic=1
            exit
        endif
        
    enddo
    if (i.eq.10**9) then
        write(*,*) "Not enough iterations"
    endif
    final_t=t
    return
end subroutine



!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
!mt19937ar

subroutine init_genrand(s)
      integer s
      integer N
      integer DONE
      integer ALLBIT_MASK
      parameter (N=624)
      parameter (DONE=123456789)
      integer mti,initialized
      integer mt(0:N-1)
      common /mt_state1/ mti,initialized
      common /mt_state2/ mt
      common /mt_mask1/ ALLBIT_MASK
!
      call mt_initln
      mt(0)=iand(s,ALLBIT_MASK)
      do 100 mti=1,N-1
        mt(mti)=1812433253*&
               ieor(mt(mti-1),ishft(mt(mti-1),-30))+mti
        mt(mti)=iand(mt(mti),ALLBIT_MASK)
  100 continue
      initialized=DONE
!
      return
      end
!-----------------------------------------------------------------------
!     initialize by an array with array-length
!     init_key is the array for initializing keys
!     key_length is its length
!-----------------------------------------------------------------------
      subroutine init_by_array(init_key,key_length)
      integer init_key(0:*)
      integer key_length
      integer N
      integer ALLBIT_MASK
      integer TOPBIT_MASK
      parameter (N=624)
      integer i,j,k
      integer mt(0:N-1)
      common /mt_state2/ mt
      common /mt_mask1/ ALLBIT_MASK
      common /mt_mask2/ TOPBIT_MASK
!
      call init_genrand(19650218)
      i=1
      j=0
      do 100 k=max(N,key_length),1,-1
        mt(i)=ieor(mt(i),ieor(mt(i-1),ishft(mt(i-1),-30))*1664525)&
                +init_key(j)+j
        mt(i)=iand(mt(i),ALLBIT_MASK)
        i=i+1
        j=j+1
        if(i.ge.N)then
          mt(0)=mt(N-1)
          i=1
        endif
        if(j.ge.key_length)then
          j=0
        endif
  100 continue
      do 200 k=N-1,1,-1
        mt(i)=ieor(mt(i),ieor(mt(i-1),ishft(mt(i-1),-30))*1566083941)-i
        mt(i)=iand(mt(i),ALLBIT_MASK)
        i=i+1
        if(i.ge.N)then
          mt(0)=mt(N-1)
          i=1
        endif
  200 continue
      mt(0)=TOPBIT_MASK
!
      return
      end
!-----------------------------------------------------------------------
!     generates a random number on [0,0xffffffff]-interval
!-----------------------------------------------------------------------
      function genrand_int32()
      integer genrand_int32
      integer N,M
      integer DONE
      integer UPPER_MASK,LOWER_MASK,MATRIX_A
      integer T1_MASK,T2_MASK
      parameter (N=624)
      parameter (M=397)
      parameter (DONE=123456789)
      integer mti,initialized
      integer mt(0:N-1)
      integer y,kk
      integer mag01(0:1)
      common /mt_state1/ mti,initialized
      common /mt_state2/ mt
      common /mt_mask3/ UPPER_MASK,LOWER_MASK,MATRIX_A,T1_MASK,T2_MASK
      common /mt_mag01/ mag01
!
      if(initialized.ne.DONE)then
        call init_genrand(21641)
      endif
!
      if(mti.ge.N)then
        do 100 kk=0,N-M-1
          y=ior(iand(mt(kk),UPPER_MASK),iand(mt(kk+1),LOWER_MASK))
          mt(kk)=ieor(ieor(mt(kk+M),ishft(y,-1)),mag01(iand(y,1)))
  100   continue
        do 200 kk=N-M,N-1-1
          y=ior(iand(mt(kk),UPPER_MASK),iand(mt(kk+1),LOWER_MASK))
          mt(kk)=ieor(ieor(mt(kk+(M-N)),ishft(y,-1)),mag01(iand(y,1)))
  200   continue
        y=ior(iand(mt(N-1),UPPER_MASK),iand(mt(0),LOWER_MASK))
        mt(kk)=ieor(ieor(mt(M-1),ishft(y,-1)),mag01(iand(y,1)))
        mti=0
      endif
!
      y=mt(mti)
      mti=mti+1
!
      y=ieor(y,ishft(y,-11))
      y=ieor(y,iand(ishft(y,7),T1_MASK))
      y=ieor(y,iand(ishft(y,15),T2_MASK))
      y=ieor(y,ishft(y,-18))
!
      genrand_int32=y
      return
      end
!-----------------------------------------------------------------------
!     generates a random number on [0,0x7fffffff]-interval
!-----------------------------------------------------------------------
      function genrand_int31()
      integer genrand_int31
      integer genrand_int32
      genrand_int31=int(ishft(genrand_int32(),-1))
      return
      end
!-----------------------------------------------------------------------
!     generates a random number on [0,1]-real-interval
!-----------------------------------------------------------------------
      function genrand_real1()
      double precision genrand_real1,r
      integer genrand_int32
      r=dble(genrand_int32())
      if(r.lt.0.d0)r=r+2.d0**32
      genrand_real1=r/4294967295.d0
      return
      end
!-----------------------------------------------------------------------
!     generates a random number on [0,1)-real-interval
!-----------------------------------------------------------------------
      function genrand_real2()
      double precision genrand_real2,r
      integer genrand_int32
      r=dble(genrand_int32())
      if(r.lt.0.d0)r=r+2.d0**32
      genrand_real2=r/4294967296.d0
      return
      end
!-----------------------------------------------------------------------
!     generates a random number on (0,1)-real-interval
!-----------------------------------------------------------------------
      function genrand_real3()
      double precision genrand_real3,r
      integer genrand_int32
      r=dble(genrand_int32())
      if(r.lt.0.d0)r=r+2.d0**32
      genrand_real3=(r+0.5d0)/4294967296.d0
      return
      end
!-----------------------------------------------------------------------
!     generates a random number on [0,1) with 53-bit resolution
!-----------------------------------------------------------------------
      function genrand_res53()
      double precision genrand_res53
      integer genrand_int32
      double precision a,b
      a=dble(ishft(genrand_int32(),-5))
      b=dble(ishft(genrand_int32(),-6))
      if(a.lt.0.d0)a=a+2.d0**32
      if(b.lt.0.d0)b=b+2.d0**32
      genrand_res53=(a*67108864.d0+b)/9007199254740992.d0
      return
      end
!-----------------------------------------------------------------------
!     initialize large number (over 32-bit constant number)
!-----------------------------------------------------------------------
      subroutine mt_initln
      integer ALLBIT_MASK
      integer TOPBIT_MASK
      integer UPPER_MASK,LOWER_MASK,MATRIX_A,T1_MASK,T2_MASK
      integer mag01(0:1)
      common /mt_mask1/ ALLBIT_MASK
      common /mt_mask2/ TOPBIT_MASK
      common /mt_mask3/ UPPER_MASK,LOWER_MASK,MATRIX_A,T1_MASK,T2_MASK
      common /mt_mag01/ mag01
!    TOPBIT_MASK = Z'80000000'
!    ALLBIT_MASK = Z'ffffffff'
!    UPPER_MASK  = Z'80000000'
!    LOWER_MASK  = Z'7fffffff'
!    MATRIX_A    = Z'9908b0df'
!    T1_MASK     = Z'9d2c5680'
!    T2_MASK     = Z'efc60000'
      TOPBIT_MASK=1073741824
      TOPBIT_MASK=ishft(TOPBIT_MASK,1)
      ALLBIT_MASK=2147483647
      ALLBIT_MASK=ior(ALLBIT_MASK,TOPBIT_MASK)
      UPPER_MASK=TOPBIT_MASK
      LOWER_MASK=2147483647
      MATRIX_A=419999967
      MATRIX_A=ior(MATRIX_A,TOPBIT_MASK)
      T1_MASK=489444992
      T1_MASK=ior(T1_MASK,TOPBIT_MASK)
      T2_MASK=1875247104
      T2_MASK=ior(T2_MASK,TOPBIT_MASK)
      mag01(0)=0
      mag01(1)=MATRIX_A
      return
      end

