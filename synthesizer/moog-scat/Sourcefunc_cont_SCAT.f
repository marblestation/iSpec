      subroutine sourcefunc_cont_scat
c********************************************************************************************************************
c     Calculates the *CONTINUUM* source function which incorporates both a
c     scattering and an absorption component.  The relevant quantities are S_cont (source function),
c     J_cont (the mean intensity), and Flux_cont (the total flux).  Employs the short characteristics methodology
c     for the solution of radiative transfer. Note that angle integration weights are computed in the 
c     Ang_Weight_SCAT subroutine.
c********************************************************************************************************************
      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Linex.com'
      include 'Scat.com'

c***  Local Arrays/Variables
      real*8  densx(100), alo(100), ood(100)
      real*8  Thomson_cont(100), Scat_opacity(100), Therm_opacity(100)
      real*8  eta(100), etat(100), dtau1_cont(100)
      real*8  dtau(100), WOP(100), W1(100), W0(100)
      real*8  J_cont_OLD(100), J_cont_moog(100), S_cont_moog(100)
      integer IT

c*** DETERMINE: rhox for MODELS that do not contain these values.  CHANGE below coding as calculation of rhox is
c    complicated.
      do i=1, ntau
         if (modtype .eq. 'KURTYPE   ' .or.
     >       modtype .eq. 'KURUCZ    ' .or.
     >       modtype .eq. 'WEB2MARC  ' .or.
     >       modtype .eq. 'WEBMARCS  ') cycle
         rhox(i) = tauref(i)/(kapref(i)/rho(i))
      enddo   

c*** SET-UP: delta_tau variable (name: dtau1_cont).
      do i=1, ntau
        ood(i)        = kaplam(i) / rho(i)
        if (i .eq. 1)  cycle
        densx(i)      = 0.5 * abs((rhox(i) - rhox(i-1)))
        dtau1_cont(i) = densx(i) * (ood(i-1) + ood(i))
        if (i .eq. ntau) cycle
      enddo

c*** SET-UP: Scat_opacity, Therm_opacity, and Thomson terms.  
c    Note the S_cont determination includes ONLY continuum-associated quantities.
c    As is standard, Therm_opacity(i) = kaplam(i) - Scat_opacity(i).
      do i=1, ntau
         Scat_opacity(i)  = kaplamsca(i)
         Therm_opacity(i) = kaplamabs(i)
         Thomson_cont(i)  = Scat_opacity(i)/kaplam(i)
      enddo

c*** ACCELERATION OF CONVERGENCE: To accelerate the convergence of the solution of the RTE, it is necessary to employ the technique of 
c    accelerated lambda iteration (ALI).  The steps below set-up the lambda operator.  
c    Convergence Requirements (original values: Converg_iter = 1.e-5 and Max_iter = 25):
      Converg_iter = 2.E-3
      Max_iter     = 65
c     Gamma controls the amount of acceleration.
c     Gamma = 0 : No Acceleration
c           = 1 : Full Acceleration
c           > 1 : Acceleration with damping  
      gamma = 1.
      if (gamma .eq. 0.) then
         do i=1, ntau
           alo(i) = 1.
         enddo
      else
        do i=1, ntau
          if (i .eq. 1) then
            dtau_alo = dtau1_cont(2)
          else if (i .eq. ntau) then
            dtau_alo = 0.5 * dtau1_cont(ntau)                                    | APPLY necessary additonal damping for the inner boundary
          else
            dtau_alo = 0.5 * (dtau1_cont(i) + dtau1_cont(i+1))                   | COMMENT: sum (/integral) of 1/mu * dmu = Sum of mu * dmu = 0.5
          endif
          dtau_alo     = dtau_alo * Thomson_cont(i) / (2.*gamma)                 | EMPLOY the damping factor of 2. * gamma (alternate: dtau_alo = dtau_alo * 1. /)
          exp_dtau_alo = exp(-dtau_alo)                                          | CREATE lambda operator with above-listed THOMSON factor
          if (exp_dtau_alo .ne. 1.) then
            alo(i) = dtau_alo / (1. - exp(-dtau_alo))
          else
            alo(i) = 1.
          endif
        enddo
      endif

c***  SET-UP: Fundamental Quantities (e.g., emissivity from planck function)
      Const1       = 1.43878858E08
      Const2       = 3.972610376E08
      H_cont       = Const1/wave
      H_cont_core  = H_cont/t(ntau)
      B_core       = Const2/((exp(H_cont_core)-1.)*wave*wave*wave)
      H_cont_core1 = H_cont/ t(ntau-1)
      B_core1      = Const2/((exp(H_cont_core1)-1.)*wave*wave*wave)
      DB_core      = B_core-B_core1

c***  SET-UP: Critical parameters such as the Planck Function (B), the Initial Mean Intensity (J_cont_OLD), etc.
      do i=1, ntau
        H_cont_t      = H_cont / t(i)
        B_planck(i)   = Const2 / ((exp(H_cont_t)-1.)*wave*wave*wave)
        J_cont_OLD(i) = B_planck(i)
        eta(i)        = B_planck(i) * Therm_opacity(i)
        etat(i)       = Scat_opacity(i) * J_cont_OLD(i)
      enddo
 
c     ---------------------------------
c***  START of the main ITERATION loop
c     ---------------------------------
      IT   = 0
      do
        IT = IT + 1
        do i=1,ntau
           J_cont_OLD(i)  = J_cont(i)
           J_cont(i)      = 0.
           S_cont(i)      = (eta(i) + etat(i)) / kaplam(i)
        enddo
c     -------------------------
c***  START of the ANGLE loop
c     -------------------------
        Flux_cont  = 0.  
        do j=1, mmu
          XI = 0.                                                          | NO INCIDENT RADIATION
          do i=2, ntau                                                     | BEGIN Depth Loop: Integration Inwards
            dtau(i)   = dtau1_cont(i) / mu(j)
            exptau    = dexp(-dtau(i))
            WOP(i)    = (exptau - 1.) / dtau(i)
            W0(i)     = 1. + WOP(i)
            W1(i)     = -exptau - WOP(i)
            XI        = XI*exptau + W0(i)*S_cont(i) + W1(i)*S_cont(i-1)
            J_cont(i) = J_cont(i) + 0.5 * wtmu(j) * XI
          enddo                                                            | End Depth Loop  
          XI = 0.5 * wtmu(j) * (B_core+(DB_core*mu(j)/dtau1_cont(ntau)))   | PREPARE Inner Boundary: multiplied by 2/dtau
          J_cont(ntau) = J_cont(ntau) + XI                                 | SET Inner Boundary Condition
          do i=ntau-1, 1, -1                                               | BEGIN Depth Loop: Integration Outwards
            dtau(i)   = dtau1_cont(i+1) / mu(j)
            exptau    = dexp(-dtau(i))
            WOP(i)    = (exptau - 1.) / dtau(i)
            W0(i)     = 1. + WOP(i)
            W1(i)     = -exptau - WOP(i)
            XI        = XI*exptau + W0(i)*S_cont(i) + W1(i)*S_cont(i+1)
            J_cont(i) = J_cont(i) + 0.5 * wtmu(j) * XI
          enddo                                                             | End Depth Loop 
          Flux_cont    = Flux_cont + XI * 0.5 * wtmu(j) * mu(j)             | CALCULATE Emergent Flux (CRITICAL OUTPUT QUANTITY!)
c     ------------------
c***  END of ANGLE Loop
c     ------------------
        enddo

c***  ACCELERATION OF CONVERGENCE: Calculate new etat and prepare to exit iteration loop (note ALI applied)
        q_max = 0.
        i_max = 0.
        do i=1, ntau
          etat_new = Scat_opacity(i) * J_cont(i)                              | SET-UP the etat_new quantity
          DE       = etat_new - etat(i)                                       | DETERMINE new value of delta(etat)
          DE       = DE * alo(i)                                              | APPLY acceleration to delta(etat)
          etat_new = etat(i) + DE                                             | CALCULATE the "accerlerated" etat_new quantity 
          QNA      = (etat_new - etat(i)) / etat_new                          | COMPUTE convergence_quantity_1 
          Q        = ABS(etat_new - etat(i)) / etat_new                       | COMPUTE convergence_quantity_2
          if (Q .gt. q_max) then
            q_max   = Q
            q_maxNA = QNA
            i_max   = i
          endif
          etat(i) = etat_new                                                   | SET-UP and pass the iterated etat quantity to next loop
        enddo
        if (q_max .le. Converg_iter .or. Max_iter .eq. 1) exit                 | DETERMINE if convergence requirements met and EXIT routine (IMPT. STEP!)
        if (IT .ge. Max_iter) then                                             | NOTIFY user if the maximum number of iterations exceeded
           write (0,'(A,F12.5,1X,1p7e10.2,1X,I3,1X,I3)')           
     >          'Maximum number of iterations exceeded : wave, q_max= ',
     >          wave, q_max, i_max, ntau
           exit                                                                | EXIT routine if maximum number of iterations exceeded
        endif
c     -----------------------
c***  END of ITERATION loop
c     -----------------------
      enddo

c***  CONVERT variables to MOOG/Edmonds format/units (if desired).
       do i=1,ntau
          S_cont_moog(i) = 2.9977518E26*(1/(wave**2))
     >                *S_cont(i)
          J_cont_moog(i) = 2.9977518E26*(1/(wave**2))
     >                *J_cont(i)
       enddo
       Flux_cont_moog = 2.9977518E26*(1/(wave**2))
     >             *Flux_cont

c     ---------------------
c***  END of ROUTINE
c     ---------------------           
      return
      end

c***  Tabulation of the flux due to pure absorption (not requisite).
c      Flux_pureabs = Flux_pureabs + B_planck(1)
c      Flux_pureabs_moog =  2.9977518E26*(1/(wave**2))
c     >             *Flux_pureabs

