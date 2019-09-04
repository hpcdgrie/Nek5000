c-----------------------------------------------------------------------
      subroutine properties_co2(h,rho,mu,alpha,beta,T)
c beta = - drho/dh

      real*8 h,rho,mu,alpha,beta,T
      real*8 a_T(5,7),a_rho(5,7),a_mu(5,7)
      real*8 a_alpha(5,7),a_beta(5,7)
      real*8 h0,h1,h2,h3,h4,h5

      h = h - 232114.0 ! reference point at 8 mpa, 301.15k

c     150 K, 280 K, 305K, 310K, 400K, 600K
      h0 = -552879.0
      h1 = -294409.0
      h2 = -211471.0
      h3 = -124847
      h4 = 47347.0
      h5 = 276505.0

      a_T(1,1) = 107.667
      a_T(1,2) = -3.00948E-3
      a_T(1,3) = -1.63250E-8
      a_T(1,4) = -4.24059E-14
      a_T(1,5) = -6.54849E-20
      a_T(1,6) = -5.63292E-26
      a_T(1,7) = -2.06318E-32
  
      a_T(2,1) = 714.511
      a_T(2,2) = 0.0115313
      a_T(2,3) = 1.25421E-7
      a_T(2,4) = 6.83923E-13
      a_T(2,5) = 2.01027E-18
      a_T(2,6) = 3.09332E-24
      a_T(2,7) = 1.96713E-30

      a_T(3,1) = 716.287
      a_T(3,2) = 0.0135318
      a_T(3,3) = 1.92477E-7
      a_T(3,4) = 1.48708E-12
      a_T(3,5) = 6.51092E-18
      a_T(3,6) = 1.51638E-23
      a_T(3,7) = 1.45133E-29
 
      a_T(4,1) = 363.556
      a_T(4,2) = 7.13355E-4
      a_T(4,3) = 1.465E-9
      a_T(4,4) = -6.57488E-15
      a_T(4,5) = 1.72399E-20
      a_T(4,6) = 2.78569E-26
      a_T(4,7) = -9.16921E-31

      a_T(5,1) = 363.534
      a_T(5,2) = 7.14428E-4
      a_T(5,3) = 1.45066E-9
      a_T(5,4) = -6.52713E-15
      a_T(5,5) = 1.73124E-20
      a_T(5,6) = -2.69883E-26
      a_T(5,7) = 1.88425E-32  

      a_rho(1,1) = -1426.43
      a_rho(1,2) = -0.0206929
      a_rho(1,3) = -8.33343E-8
      a_rho(1,4) = -2.04659E-13
      a_rho(1,5) = -2.91442E-19
      a_rho(1,6) = -2.23045E-25
      a_rho(1,7) = -7.06877E-32

      a_rho(2,1) = 9233.81
      a_rho(2,2) = 0.231914
      a_rho(2,3) = 2.3851E-6
      a_rho(2,4) = 1.25694E-11
      a_rho(2,5) = 3.67100E-17
      a_rho(2,6) = 5.67385E-23
      a_rho(2,7) = 3.6372E-29
  
      a_rho(3,1) = 6614.25
      a_rho(3,2) = 0.244217
      a_rho(3,3) = 3.82971E-6
      a_rho(3,4) = 3.15155E-11
      a_rho(3,5) = 1.44333E-16
      a_rho(3,6) = 3.46696E-22
      a_rho(3,7) = 3.40351E-28

      a_rho(4,1) = 149.492
      a_rho(4,2) = -6.66837E-4
      a_rho(4,3) = 3.24725E-9
      a_rho(4,4) = -1.40349E-14
      a_rho(4,5) = 5.70357E-20
      a_rho(4,6) = -3.23439E-25
      a_rho(4,7) = -1.73987E-30

      a_rho(5,1) = 149.349
      a_rho(5,2) = -6.5753E-4
      a_rho(5,3) = 3.01580E-9
      a_rho(5,4) = -1.14388E-14
      a_rho(5,5) = 3.10613E-20
      a_rho(5,6) = -5.14415E-26
      a_rho(5,7) = 3.81899E-32

      a_mu(1,1) = 0.0123092
      a_mu(1,2) = 1.95612E-7
      a_mu(1,3) = 1.29383E-12
      a_mu(1,4) = 4.54965E-18
      a_mu(1,5) = 9.00169E-24
      a_mu(1,6) = 9.51635E-30
      a_mu(1,7) = 4.22977E-36
  
      a_mu(2,1) = 0.00126419
      a_mu(2,2) = 2.97986E-08
      a_mu(2,3) = 2.90930E-13
      a_mu(2,4) = 1.49024E-18
      a_mu(2,5) = 4.28494E-24
      a_mu(2,6) = 6.56018E-30
      a_mu(2,7) = 4.20210E-36
  
      a_mu(3,1) = 4.84806E-4
      a_mu(3,2) = 1.71068E-08
      a_mu(3,3) = 2.58375E-13
      a_mu(3,4) = 2.04412E-18
      a_mu(3,5) = 8.93188E-24
      a_mu(3,6) = 2.02554E-29
      a_mu(3,7) = 1.86029E-35

      a_mu(4,1) = 2.02437E-5
      a_mu(4,2) = 1.59633E-11
      a_mu(4,3) = 1.74130E-16
      a_mu(4,4) = -1.01322E-21
      a_mu(4,5) = 4.42989E-27
      a_mu(4,6) = -1.99263E-32
      a_mu(4,7) = -9.53787E-39
 
      a_mu(5,1) = 2.02331E-5
      a_mu(5,2) = 1.66431E-11
      a_mu(5,3) = 1.56597E-16
      a_mu(5,4) = -7.80107E-22
      a_mu(5,5) = 2.26872E-27
      a_mu(5,6) = -3.83642E-33
      a_mu(5,7) = 2.86644E-39

      a_alpha(1,1) = 0.0000462205
      a_alpha(1,2) = 1.90974E-9
      a_alpha(1,3) = 1.72367E-14
      a_alpha(1,4) = 6.52294E-20
      a_alpha(1,5) = 1.34204E-25
      a_alpha(1,6) = 1.43153E-31
      a_alpha(1,7) = 6.17044E-38

      a_alpha(2,1) = 0.0306225
      a_alpha(2,2) = 7.96386E-7
      a_alpha(2,3) = 8.53845E-12
      a_alpha(2,4) = 4.83780E-17
      a_alpha(2,5) = 1.52923E-22
      a_alpha(2,6) = 2.55810E-28
      a_alpha(2,7) = 1.76992E-34

      a_alpha(3,1) = -0.00174308
      a_alpha(3,2) = -6.74465E-08
      a_alpha(3,3) = -1.07161E-12
      a_alpha(3,4) = -9.01347E-18
      a_alpha(3,5) = -4.24279E-23
      a_alpha(3,6) = -1.05939E-28
      a_alpha(3,7) = -1.09352E-34
  
      a_alpha(4,1) = 0.0000199554
      a_alpha(4,2) = 9.20474E-11
      a_alpha(4,3) = -1.27884E-16
      a_alpha(4,4) = -2.55165E-22
      a_alpha(4,5) = -1.69803E-27
      a_alpha(4,6) = 1.37300E-31
      a_alpha(4,7) = 8.79488E-37

      a_alpha(5,1) = 2.03651E-5
      a_alpha(5,2) = 7.39909E-11
      a_alpha(5,3) = 1.51522E-16
      a_alpha(5,4) = -1.97392E-21
      a_alpha(5,5) = 8.47239E-27
      a_alpha(5,6) = -1.74011E-32
      a_alpha(5,7) = 1.41211E-38
      
      a_beta(1,1) = 0.0408315
      a_beta(1,2) = 4.61529E-7
      a_beta(1,3) = 2.39881E-12
      a_beta(1,4) = 6.88361E-18
      a_beta(1,5) = 1.13411E-23
      a_beta(1,6) = 1.01057E-29
      a_beta(1,7) = 3.79171E-36

      a_beta(2,1) = -1.99041
      a_beta(2,2) = -4.66546E-5
      a_beta(2,3) = -4.52483E-10
      a_beta(2,4) = -2.33282E-15
      a_beta(2,5) = -6.75036E-21
      a_beta(2,6) = -1.03995E-26
      a_beta(2,7) = -6.66509E-33

      a_beta(3,1) = 0.311298
      a_beta(3,2) = 1.2605E-5
      a_beta(3,3) = 2.11710E-10
      a_beta(3,4) = 1.87729E-15
      a_beta(3,5) = 9.27157E-21
      a_beta(3,6) = 2.41291E-26
      a_beta(3,7) = 2.57938E-32
  
      a_beta(4,1) = 0.000666531
      a_beta(4,2) = -6.49321E-9
      a_beta(4,3) = 4.35547E-14
      a_beta(4,4) = -2.13127E-19
      a_beta(4,5) = 8.88253E-25
      a_beta(4,6) = -3.05950E-30
      a_beta(4,7) = -6.04744E-35

      a_beta(5,1) = 0.000664008
      a_beta(5,2) = -6.33472E-9
      a_beta(5,3) = 3.99329E-14
      a_beta(5,4) = -1.76700E-19
      a_beta(5,5) = 5.19251E-25
      a_beta(5,6) = -8.96889E-31
      a_beta(5,7) = 6.81624E-37

c if enthalpy is too small or too large, break.  
c      if (h.lt.h0) then
c       if (nid.eq.0) write(6,*) "h too small.Temperature is below 150K"
c        call exitt()
c      endif

c      if (h.gt.h5) then
c        if (nid.eq.0) write(6,*) "h too large.Temperature is above 600K"
c        call exitt()
c      endif
c=====================================================================================
      if ((h.lt.h1).and.(h.ge.h0)) then
          T = 0.0
          do i = 1,7
            T = T + a_T(1,i)*h**dble(i-1)
          enddo

          rho = 0.0
          do i = 1,7
            rho = rho + a_rho(1,i)*h**dble(i-1)
          enddo

          mu = 0.0
          do i = 1,7
            mu = mu + a_mu(1,i)*h**dble(i-1)
          enddo

          alpha = 0.0
          do i = 1,7
           alpha = alpha + a_alpha(1,i)*h**dble(i-1) 
          enddo
	   
         beta = 0.0
          do i = 1,7
           beta = beta + a_beta(1,i)*h**dble(i-1) 
         enddo
c         do i = 2,7
c           beta = beta + dble(i-1)*a_rho(1,i)*h**dble(i-2)
c         enddo
      elseif ((h.lt.h2).and.(h.ge.h1)) then
          T = 0.0
          do i = 1,7
            T = T + a_T(2,i)*h**dble(i-1)
          enddo

          rho = 0.0
          do i = 1,7
            rho = rho + a_rho(2,i)*h**dble(i-1)
          enddo

          mu = 0.0
          do i = 1,7
            mu = mu + a_mu(2,i)*h**dble(i-1)
          enddo

          alpha = 0.0
          do i = 1,7
           alpha = alpha + a_alpha(2,i)*h**dble(i-1) 
          enddo
	   
         beta = 0.0
          do i = 1,7
           beta = beta + a_beta(2,i)*h**dble(i-1) 
         enddo
c         do i = 2,7
c           beta = beta + dble(i-1)*a_rho(2,i)*h**dble(i-2)
c         enddo
      elseif ((h.lt.h3).and.(h.ge.h2)) then
          T = 0.0
          do i = 1,7
            T = T + a_T(3,i)*h**dble(i-1)
          enddo

          rho = 0.0
          do i = 1,7
            rho = rho + a_rho(3,i)*h**dble(i-1)
          enddo

          mu = 0.0
          do i = 1,7
            mu = mu + a_mu(3,i)*h**dble(i-1)
          enddo

          alpha = 0.0
          do i = 1,7
           alpha = alpha + a_alpha(3,i)*h**dble(i-1) 
          enddo
	   
         beta = 0.0
          do i = 1,7
           beta = beta + a_beta(3,i)*h**dble(i-1) 
         enddo
c         do i = 2,7
c           beta = beta + dble(i-1)*a_rho(3,i)*h**dble(i-2)
c         enddo
      elseif ((h.lt.h4).and.(h.ge.h3)) then
          T = 0.0
          do i = 1,7
            T = T + a_T(4,i)*h**dble(i-1)
          enddo

          rho = 0.0
          do i = 1,7
            rho = rho + a_rho(4,i)*h**dble(i-1)
          enddo

          mu = 0.0
          do i = 1,7
            mu = mu + a_mu(4,i)*h**dble(i-1)
          enddo

          alpha = 0.0
          do i = 1,7
           alpha = alpha + a_alpha(4,i)*h**dble(i-1) 
          enddo
	   
         beta = 0.0
          do i = 1,7
           beta = beta + a_beta(4,i)*h**dble(i-1) 
         enddo
c         do i = 2,7
c           beta = beta + dble(i-1)*a_rho(4,i)*h**dble(i-2)
c         enddo
      elseif ((h.lt.h5).and.(h.ge.h4)) then
          T = 0.0
          do i = 1,7
            T = T + a_T(5,i)*h**dble(i-1)
          enddo

          rho = 0.0
          do i = 1,7
            rho = rho + a_rho(5,i)*h**dble(i-1)
          enddo

          mu = 0.0
          do i = 1,7
            mu = mu + a_mu(5,i)*h**dble(i-1)
          enddo 
          alpha = 0.0
          do i = 1,7
           alpha = alpha + a_alpha(5,i)*h**dble(i-1) 
          enddo
	   
         beta = 0.0
          do i = 1,7
           beta = beta + a_beta(5,i)*h**dble(i-1) 
         enddo
c         do i = 2,7
c           beta = beta + dble(i-1)*a_rho(5,i)*h**dble(i-2)
c         enddo
      endif
c      beta = -1.0*beta

      return
      end
c-----------------------------------------------------------------------