c      program test
c      pi = acos(-1.0)
c      icon = 0
c      call  sigma_en_bound(icon,nuc_id,ei,q2,p_ini,th_i,phi_i,
c     & sigma_en_df1,sigma_en_st,sigma_en_free,sigma_en_lc,eta,ifl)
c
c      icon = 1
c      nuc_id = 1
c      q2 = 2.0
c      ei = 5.76
c      p_ini = 0.7
c      th_i  = 170.0*pi/180.0
c      phi_i = 0.0*pi/180.0 
c      
c      call  sigma_en_bound(icon,nuc_id,ei,q2,p_ini,th_i,phi_i,
c    &     sigma_en_df1,sigma_en_st,sigma_en_free,sigma_en_lc,eta,ifl)
c
c      print *,p_ini,sigma_en_df1,sigma_en_st,sigma_en_free,sigma_en_lc,
c     &        eta
c      end
      
      subroutine sigma_en_bound(icon,nuc_id,ei,q2,p_ini,th_i,phi_i,
     & sigma_en_df1,sigma_en_st,sigma_en_free,sigma_en_lc,eta,ifl)
************************************************************************
*     icon  = 0 - initialization, 1 - calculation
*     nuc_id = 1 - proton -1 - neutron
*     ei = initial electron energy GeV
*     q2 = Q2
*     p_ini = inital momenutm of the nucleon GeV/c
*     th_i = angle of initial nucleon vs q - in radians
*     phi_i = azimuthal angle of initial nucleon in the Ref frame with (qe) = (zx)
*     sigma_en_df1  - de-forest CC1
*     sigma_en_St   - for stational nucleon with Ei and theta_e
*     sigma_en_free - for free nucleon with the same initial momentum
*     sigma_en_lc   - light-cone approximation
*     eta           -   virtuality parameter Eq.(30) Phys. Rev C08, 035202, 2018
*                       https://arxiv.org/abs/1805.04639
*     ifl           - 0 - ok, 1 - kinematically forbidden
*      
*      
*     FIU, July 2019
*     Miami, FL
***************************************************************************      
      common/parms/pipi,pmpm,dmdm
      common/par/pi,pm,pmp,pmn,dm,eb
      real       pi,pm,pmp,pmn,dm,eb
      common/electron/eii,esc,ue
      real            ei,esc,ue,eii
      common/ELECTRON_en/ei_i,esc_i,ue_i
      common/PHOTON_en/Q2_i,q0_i,qv_i
      common/MASS_en/pm_i,dm_i,pi_i 


      
      dimension q(3),p2(3),p1(3),s1(3),s2(3)
      real, dimension(2,2) :: jN_re,jN_im
      complex, dimension(0:3,2,2) :: jN,jN_b,jN_off_vn,jN_ins_lf
      complex, dimension(4,4) :: gamma_pl,gamma_mn,prod
      common/gamma_lf/gamma_pl,gamma_mn
      real v_l,v_t,v_tl,v_tt
      real w_l,w_t,w_tl,w_tt
      complex, dimension(0:3,0:3) :: WN_s   , WN_s_b
      complex, dimension(0:3,0:3) :: WN_sum , WN_sum_b

      sigma_en_df1  = 0.0
      sigma_en_free = 0.0
      sigma_en_lc  = 0.0

      noff = 0

      if(icon.eq.0)then
      pi  = acos(-1.0)
      pi_i = pi
        pm  = 0.938279         
        pm_i = pm
        dm  = 1.875628
        dm_i = dm
        pmp = pm
        pmn = 0.939565  
	eb  = 0.00221
        pipi = pi
        pmpm = pm
        dmdm = dm
        alpha = 1.0/137.0

*        a = 12.0
*        z = 6.0
*       call nuclear_masses(a,z,am)


        
        nuc = 0                 ! initialization
      call j_N_matrix_s_onof(nuc,q2,q,p2,s2,p1,s1,jN,jN_off_vn,jN_ins_lf
     &       ,e_off,ppl_off,eta,noff)

      return
      endif
      
        nuc  = 1             ! the knock-out nucleon is proton(1), neutron(-1)
        nuc = nuc_id

      
*     ei   = 5.76          ! initial electron energy
        eii = ei
        pr   = p_ini !0.5       ! recoil nucleon momentum GeV/C
        thr  = pi - th_i      !30.*pi/180.   ! polar angle of the recoil nucleon
        phir = pi + phi_i     !pi  ! azimuthal angle of the recoil nucleon
        prx = pr*sin(thr)*cos(phir)
        pry = pr*sin(thr)*sin(phir)
        prz = pr*cos(thr)
        
*        q2 = 4.0
        ei_i = ei
        q2_i = q2

        print *,ei,pr,thr,phir,q2
        
************************************************************************
*       definition q0 by the QE condition of q+d --> pf + pr scattering
************************************************************************
	ifl    = 0
 	call kine_q0(q2,pr,thr,q0,ifl)
	if(ifl.eq.1)return !goto 2
 	esc  = ei - q0         !scattered electron energy
	arg = sqrt(q2/4.0/ei/esc)
	if(arg**2.gt.1.0)then
        ifl =  1                 !goto 2
        return
        endif
      
	ue  = 2.0*asin(arg)
        ue_i = ue
        
 	qv  = sqrt(q2 + q0**2)
        q0_i = q0
        qv_i = qv
        
        x = q2/2.0/pm/q0
 
        q(1)  = 0.0; q(2)  = 0.0;  q(3) = qv

****************************************************************
*         Definition of Electron Current's structure functions
*         in the reference frame in which q||z
*****************************************************************
	cs_th_eq = (ei - esc*cos(ue))/qv
	sn_th_eq =       esc*sin(ue) /qv
	phi_eq = 0.0
	th_eq = acos(cs_th_eq)

	cs_th_escq = (ei*cos(ue)-esc)/qv
	sn_th_escq =  ei*sin(ue)/qv
	phi_escq   = 0.0
	if((th_eq + ue).ge.pi)then
	phi_escq   = pi
	endif

         v_l  = q2**2/qv**4
         v_t  = q2/2.0/qv**2 + tan(ue/2.0)**2
         v_tt = q2/2.0/qv**2
         v_tl = q2/qv**2*sqrt(q2/qv**2 + tan(ue/2.0)**2)

****************************************************************
*         Definition of Structure Functions of nucleon
*         in the reference frame in which q||z
*****************************************************************
        s1(1) = 0.0;  s1(2) =  0.0; s1(3) = 1.0  ! initial nucleon's polarization axis
        p1(1) = -prx; p1(2) = -pry; p1(3) = -prz ! initial nucleon's momentum        
        ein = sqrt(pm**2 + prx**2 + pry**2 + prz**2)
        
        s2(1) = 0.0; s2(2) = 0.0; s2(3) = 1.0 ! final nucleon's polarization axis       
        p2(1) = p1(1)+q(1); p2(2) = p1(2)+q(2); p2(3) = p1(3)+q(3) !final nucleon momentum
        pf = sqrt(p2(1)**2 + p2(2)**2 + p2(3)**2)
        ef = sqrt(pm**2+pf**2)
*     nuc = 1
       jN_ins_lf = 0.0 
       jN_off_vn = 0.0 

      call j_N_matrix_s_onof(nuc,q2,q,p2,s2,p1,s1,jN,jN_off_vn,jN_ins_lf
     &                                          ,e_off,ppl_off,eta,noff)

*      print *,p1
      jN_b = jN + jN_ins_lf

      WN_sum   = (0.0,0.0)
      WN_sum_b = (0.0,0.0)
      
      do ks2 = 1,2
      do ks1 = 1,1
         do mu = 0,3
         do nu = 0,3
         WN_s(mu,nu) =  conjg(jN(mu,ks2,ks1))*jN(nu,ks2,ks1)     
         WN_s_b(mu,nu)= conjg(jN_b(mu,ks2,ks1))*jN_b(nu,ks2,ks1)     
         enddo
         enddo
      WN_sum = WN_sum + WN_s
      WN_sum_b = WN_sum_b + WN_s_b
      enddo
      enddo

*      print *,WN_sum
         
         V_L_N = qv**4/Q2**2 *
     &     (real(WN_sum(0,0)) - 2.0*q0/qv*real(WN_sum(0,3))
     &                        + q0**2/qv**2*real(WN_sum(3,3)))
         V_T_N =   Real(WN_sum(2,2)) +  Real(WN_sum(1,1))
         V_TT_N =  Real(WN_sum(2,2)) -  Real(WN_sum(1,1))
         V_TL_N = 2.0*qv**2/Q2*
     &      (q0/qv* Real(WN_sum(3,1)) - Real(WN_sum(0,1)))

*         a_un = (1.0/2.0)*  ! initial spin average
         a_un = 1.0 *   
     &    (v_l*V_L_N + v_t*V_t_N  + v_tl*V_TL_N + v_tt*V_TT_N)
*         print *, "aun",a_un

***************************************************************
*     Bound nucleon strucuture function calculations
***************************************************************         
         
         V_L_N_b = qv**4/Q2**2 *
     &     (real(WN_sum_b(0,0)) - 2.0*q0/qv*real(WN_sum_b(0,3))
     &                        + q0**2/qv**2*real(WN_sum_b(3,3)))
         V_T_N_b  =   Real(WN_sum_b(2,2)) +  Real(WN_sum_b(1,1))
         V_TT_N_b =   Real(WN_sum_b(2,2)) -  Real(WN_sum_b(1,1))
         V_TL_N_b = 2.0*qv**2/Q2*
     &      (q0/qv* Real(WN_sum_b(3,1)) - Real(WN_sum_b(0,1)))

         a_un_b = 1.0 *  
     &   (v_l*V_L_N_b + v_t*V_t_N_b  + v_tl*V_TL_N_b + v_tt*V_TT_N_b)
*         print *, "aun",a_un
         

         sigma_mott = alpha**2*cos(ue/2.0)**2/(4.0*ei**2*sin(ue/2.0)**4)
     &   *0.389385*1000.*1000.  ! nanobar
*         sigma_en0 = sigma_mott/(2.0*dm*ef)*a_un
         sigma_en_free  = sigma_mott/(2.0*ein*2.0*ef)*a_un
         sigma_en_lc    = sigma_mott/(2.0*(dm/2.)*2.0*ef)*a_un_b


         px = -prx
         py = -pry
         pz = -prz
*         nuc = 1
         nsc = 1
         sigma_en_st  = sigma_en_stat(nuc,ei,esc,ue,q2,q0)
         sigma_en_df1 = sigma_en(px,py,pz,nuc,nsc)

*         sigma_en_df1 = s_en_df
*         sigma_en_lc = sigma_en_bound
         
*        print *,"sigma",sigma_en_free,s_en_df,s_en_st,sigma_en_bound,eta
*********************************************************************
  2    continue
       return
       end


      real function sigma_en_stat(nuc,ei,esc,ue,q2,q0)
      implicit none     
      real q2,q0,tau,ef
      integer nuc
      common/par/pi,pm,pmp,pmn,dm,eb
      real       pi,pm,pmp,pmn,dm,eb,alpha
      real       ei,esc,ue
      real  sigma_mott,G_E,G_M,F1s,F2s
      real ph_fc,w2,w1,sen0 
      alpha = 1.0/137.0
      sigma_mott = alpha**2*cos(ue/2.0)**2/(4.0*ei**2*sin(ue/2.0)**4)
      Ef = pm + q0

      
      tau = q2/4.0/pm**2
      G_E = F1s(nuc,q2) - tau * F2s(nuc,q2)     
      G_M = F1s(nuc,q2) + F2s(nuc,q2)

      ph_fc  = 1.0/(2.0*pm)/(2.0*Ef)
      ph_fc  = pm/Ef
      w2 = (G_E**2 + tau * G_M**2)/(1.0+tau)
      w1 = tau*G_M**2
      sen0 = sigma_mott*ph_fc*(w2 + 2.0*tan(ue/2)**2 *w1)
      sigma_en_stat = sen0*0.389385*1000.*1000.
      return
      end

	subroutine kine_q0(q2,pr,thr,q0,ifl)
	implicit none     
	real               q2,pr,thr,q0
	integer ifl
        common/par/pi,pm,pmp,pmn,dm,eb
        real       pi,pm,pmp,pmn,dm,eb
        common/electron/eil,erl,ue
        real            eil,erl,ue
	integer i_test
	real er,arm_a,arm_b,arm_g
	real A,B,C,x1,x2
	real test

 	ifl   = 1
	er    = sqrt(pm**2 + pr**2)
   	arm_a = -q2 + dm**2 - 2.0*er*dm
	arm_b = 2.0*(dm-er)
	arm_g = -2.0*pr*cos(thr)

*	if(arm_g.eq.0.0)then
	if(arm_g**2.le.0.00001)then
	q0    = - arm_a/arm_b
	if(q0.lt.0.0.or.q0.gt.eil)return
	ifl   = 0.0
	return
	endif
	A     = arm_b**2 - arm_g**2
	B     = 2.0*arm_a*arm_b
	C     = arm_a**2 - q2*arm_g**2
	if(A.eq.0.0)then
	q0 = -C/B
	if(q0.lt.0.or.q0.gt.eil)return
 	else
   	call  squar_eq(a,b,c,x1,x2,i_test)
	if(i_test.eq.1)return
*	write(6,*)'xx',thr*180.0/pi,a,b,c,x1,x2
    	q0 = x2
	test  = (arm_a + q0*arm_b)/arm_g  
	if(test.lt.0.0.or.q0.lt.0.0.or.q0.gt.eil)then
     	q0 = x1
	test  = (arm_a + q0*arm_b)/arm_g  
	if(test.lt.0.0.or.q0.lt.0.0.or.q0.gt.eil)return
	endif
	endif
	ifl   = 0	
 	return
 	end



        SUBROUTINE SQUAR_EQ(A,B,C,X1,X2,I_TEST)
	implicit none
	real a,b,c,x1,x2,dskm
	integer i_test
        I_TEST = 0
        DSKM  = B**2 - 4.0*A*C
        IF(DSKM.LT.0.0)THEN
                IF(DSKM.GT.-1.0E-4)THEN
                DSKM = 0.0
                GOTO 111
                ELSE
                I_TEST = 1
                ENDIF
        RETURN
        ENDIF
111     X1 = ( -B + SQRT(DSKM) ) / 2.0 / A
        X2 = ( -B - SQRT(DSKM) ) / 2.0 / A
        RETURN
        END




      
 
      subroutine j_N_matrix_s_onof(nuc,q2,q,p2,s2,p1,s1,jN,jN_off_vn,
     &                             jN_ins_lf,e_off,ppl_off,eta,noff) 
********************************************************************************************
* Calculates electromagnetic current of nucleon
* for initial nucleon polarization s1 and final s2
* up - means polarization along s1 or s2
* down - polarization opposite to s1 and s2
*  <s2,p2|J^\mu|p1,s1>  = \bar u(p_2,s_2)\Gamma^\mu u(p_1,s_1) - jN_re/im(s2,s1), where
*  1,1 = u2,u1, 1,2=u2,d1,  2,1=d2,u1, 2,2=d2,d1
*
* and polarization axises for initial lepton is defined by s1 
*     and for final lepton s2
* here \Gamma_\mu = \gamma_mu F1 + i/(2M) F2 q^\nu \sigma_{\mu\nu}
*
* nuc - (0) initializes the \gamma and \sigma matrices
* nuc - (1) - proton (-1) neutron
* q2 - Q2
* q  - q(3)  three momentum of the virtual photon
* p1 - p1(3) three momentum of the initial nucleon
* s1 - s1(3) polarization of initial nucleon
* p2 - p2(3) three momentum of the final nucleon 
* s2 - s2(3) polarization of final nucleon
* jN(0:3,2f,2i) - complex out - returns re and im part of the current
*               0:3 are mu components,  2f - final   nucleon's spin 1-up 2-down
*                                       2i - initial nucleon's spin 1-up 2-down
*up and down means parallel and aniparallel to s1 (for initial) and s2 (for final)
*! common/parms/pi,pm,dm should be provided 
*
*jN_off_vn off-shell current ubar Gamma gamma_0 u (e_off-e_on)/2pm
*in virtual nucleon approximation
*
*jN_ins_lf instantaneous current  ubar Gamma gamma_mn u (p+_off-p+_on)/4pm
*in Light-Front approximation
*      
* ics - defines the case for form-factor parameterization
* inside the subroutine
* if the code is used within other code then ics should
* be passed by common/
*
*     z axis is defined along q
*      
* Misak Sargsian 
* March 2007,FIU, Miami
* March-April 2008,FIU, Miami
* error in Kelly's parameterization is corrected
* June 2010, Miami
* Light-Front approximation case is added
* November 2017 Miami
****************************************************************************************
      implicit none
      complex, dimension(0:3,4,4):: gamma,bgamma
      complex, dimension(4,4) :: gamma_pl,gamma_mn
      complex, dimension(0:3,0:3,4,4) :: sigma
      complex, dimension(4):: u1,uc1,u2,uc2,u2bar,ul0,ulx,uly,ulz,u1bar
      complex, dimension(4)::g0u1,g3u1,gmnu1
      complex, dimension(0:3,4):: ucgamma,ucmat,u2bar_bgamma,gammau
      real,    dimension(3) :: q,p2,p1
      real,    dimension(3)   :: s2,s1
      complex, dimension(0:3,2,2):: jN,jN_off_vn,jN_ins_lf

      complex c,cm,j0,jx,jy,jz
      real pi,pm,dm
      real p1x,p1y,p1z,e1,p2x,p2y,p2z,e2,qx,qy,qz,q0,q2,qv
      real pf_pl, qpl
      real e_on,e_off,ppl_on,ppl_off,f_off
      real alpha_i,alpha_s,pt2
      real eta
      integer nuc,is1,isp1,is2,isp2      
      real norm_test
      common/gammamatices_s/gamma,sigma
      common/gamma_lf/gamma_pl,gamma_mn
      common/parms/pi,pm,dm
      common/formfcase/ics
      integer ics,noff

      jN = cmplx(0.0,0.0)
      jN_off_vn = cmplx(0.0,0.0)
      jN_ins_lf = cmplx(0.0,0.0)

      if(nuc.eq.0)then

      ics = 1  ! defines the parameterization of nucleon form-factors
               ! 1 - SLAC, 2 - Kelly, 3-Bodek, BB, Arrington
               ! 4 - Kelly with different Gen ( =0)
               ! 5 - Kelly with different Gen (2Gen_Galster)
      call gammamatricess  ! initializes gamma and sigma matrices

      return
      endif

      p1x = p1(1)
      p1y = p1(2)
      p1z = p1(3)
      e1  = sqrt(pm**2 + p1x**2 + p1y**2 + p1z**2)
      e_on = e1
      e_off = 2.0*pm - e_on
      alpha_s = (e_on + p1z)/(dm/2.0)
      alpha_i = 2.0 - alpha_s
      pt2 =  p1x**2 + p1y**2
      ppl_on = (pm**2 +  pt2)/(dm/2.*alpha_i)

      
      p2x = p2(1)
      p2y = p2(2)
      p2z = p2(3)
      e2  = sqrt(pm**2 + p2x**2 + p2y**2 + p2z**2)
      
      pf_pl = e2 + p2z
      
      
      qx = q(1)
      qy = q(2)
      qz = q(3)
      q0 = sqrt(qx**2+qy**2+qz**2-q2)
      qpl = q0 + qz

      ppl_off = pf_pl - qpl
      ppl_off = dm - (pm**2 + pt2)/(dm/2.0*alpha_s)
      

      


      do is1 = 1,2
                     isp1 =   1
         if(is1.eq.2)isp1 =  -1
      do is2 = 1,2
                     isp2 =  1
         if(is2.eq.2)isp2 = -1
      
*******************************************************
*      Spinor of initial nucleon
*******************************************************
      call  us_s(pm,p1x,p1y,p1z,e1,s1,isp1,u1,uc1)
*      print *,"us",u1
      call vec_mat_prods_inv(gamma,u1,gammau)  
      g0u1 = gammau(0,1:4)
      g3u1 = gammau(3,1:4)
      gmnu1 = g0u1 - g3u1
*      write(24,*)gmnu1
*      write(24,*)"*************"
*      write(24,*)g0u1
*      norm_test  = u1bar(1)*u1(1) + u1bar(2)*u1(2) + 
*     &             u1bar(3)*u1(3) + u1bar(4)*u1(4)
*      print *,"norm",norm_test
*      call  u_s(pm,p1x,p1y,p1z,e1,isp1,u1,uc1)
*      print *,u1
*******************************************************
*      Spinor of final nucleon
*******************************************************
      call  us_s(pm,p2x,p2y,p2z,e2,s2,isp2,u2,uc2)
*      print *,"s",s
*      call  u_s(pm,p2x,p2y,p2z,e2,isp2,u2,uc2)
*      a =   matmul(gamma(1,1:4,1:4),gamma(1,1:4,1:4))
*      b =   matmul(gamma(1,1:4,1:4),gamma(1,1:4,1:4))
*      call vec_vec_sprod(uc1,u1,c)
********************************************************
*     Calculating u2bar
********************************************************
      call vec_mat_prods(uc2,gamma,ucgamma)  
      u2bar = ucgamma(0,1:4)
*      print *,u2bar
******************************************************************
*      Calculating nucleon vertex
******************************************************************
      call nucleon_vertex(nuc,q0,qx,qy,qz,bgamma)

*      print *, "bgamma",bgamma(0,1:4,1:4)
*****************************************************************
*      Calculating \bar u2 \Gamma part
*****************************************************************
      call vec_mat_prods(u2bar,bgamma,u2bar_bgamma)  
      
      ul0 = u2bar_bgamma(0,1:4)
      ulx = u2bar_bgamma(1,1:4)
      uly = u2bar_bgamma(2,1:4)
      ulz = u2bar_bgamma(3,1:4)

*      call vec_vec_sprod(ul0,u2,c)
*      print *, "c",c,2.*e2
*      call vec_vec_sprod(uc2,u1,c)
*      print *, "c",c

****************************************************************
*      Completing the calculation of current components
****************************************************************
      call vec_vec_sprods(ul0,u1,j0)
      call vec_vec_sprods(ulx,u1,jx)
      call vec_vec_sprods(uly,u1,jy)
      call vec_vec_sprods(ulz,u1,jz)

      jN(0,is2,is1) = j0
      jN(1,is2,is1) = jx
      jN(2,is2,is1) = jy
      jN(3,is2,is1) = jz


****************************************************************
*     Completing the calculation of off shell current components
*     in VN approximation
****************************************************************
      call vec_vec_sprods(ul0,g0u1,j0)
      call vec_vec_sprods(ulx,g0u1,jx)
      call vec_vec_sprods(uly,g0u1,jy)
      call vec_vec_sprods(ulz,g0u1,jz)

      f_off = (e_off-e_on)/(2.0*pm)
*      print *,"f_off", f_off,e_on,e_off,pm

      jN_off_vn(0,is2,is1) = j0*f_off
      jN_off_vn(1,is2,is1) = jx*f_off
      jN_off_vn(2,is2,is1) = jy*f_off
      jN_off_vn(3,is2,is1) = jz*f_off

****************************************************************
*     Completing the calculation of the instantaneous current components
*     in Light Front approximation
****************************************************************
      call vec_vec_sprods(ul0,gmnu1,j0)
      call vec_vec_sprods(ulx,gmnu1,jx)
      call vec_vec_sprods(uly,gmnu1,jy)
      call vec_vec_sprods(ulz,gmnu1,jz)

  
      f_off = (ppl_off-ppl_on)/(4.0*pm)
*      f_off = (1.0/4.0/pm)*
*     &        (1.0/dm)*(dm**2 - 4.0*(pm**2+pt2)/alpha_i/(2.0-alpha_i))

      eta = (1.0/Q2)*(4.0*(pm**2+pt2)/alpha_i/(2.0-alpha_i)-dm**2)
      
      jN_ins_lf(0,is2,is1) = j0*f_off
      jN_ins_lf(1,is2,is1) = jx*f_off
      jN_ins_lf(2,is2,is1) = jy*f_off
      jN_ins_lf(3,is2,is1) = jz*f_off
    

*      print *, "j0",j0
      enddo
      enddo

      return
      end

      subroutine vec_vec_sprods(a,b,c)
      complex, dimension(4) :: a, b
      complex c
      c = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)+a(4)*b(4)
      return
      end

      subroutine vec_mat_prods(a,b,pr)
******************************************
*     Calculates gamma^mu u
******************************************      
      complex, dimension(4) :: a
      complex, dimension(0:3,4,4) :: b
      complex, dimension(0:3,4):: pr
      pr = cmplx(0.0,0.0)
      do mu = 0,3
      do j = 1,4
      do i = 1,4
      pr(mu,j)  = pr(mu,j) + a(i)*b(mu,i,j)
      enddo
      enddo
      enddo
*      print *,"*",a
*      print *,"***", pr(0,1:4)
      return
      end

      subroutine vec_mat_prods_inv(bm,a,pr)
      complex, dimension(4) :: a
      complex, dimension(0:3,4,4) :: bm
      complex, dimension(0:3,4):: pr
      pr = cmplx(0.0,0.0)
      do mu = 0,3
      do j = 1,4
      do i = 1,4
      pr(mu,j)  = pr(mu,j) + a(i)*bm(mu,j,i)
      enddo
      enddo
      enddo
*      print *,"*",a
*      print *,"***", pr(0,1:4)
      return
      end




      subroutine nucleon_vertex(nuc,q0,qx,qy,qz,bgamma)
***************************************************************************
*  This subroutine calculates the elastic electromagnetc vertex of nucleon
*   \Gamma_\mu = \gamma_mu F1 + i/(2M) F2 q^\nu \sigma_{\mu\nu}
* 
*   nuc = 1 proton, -1 neutron
*   q0,qx,qy,qz - energy and momentum components of virtual photon
*   bgamma(0:3,4,4)  \gamma_m  4x4 matrix for \mu = 0-3
*
*
* Misak Sargsian
* March 2007, FIU, Miami
******************************************************************************
      implicit none
      common/parms/pi,pm,dm
      common/gammamatices_s/gamma,sigma
      complex, dimension(0:3,4,4):: gamma,bgamma,sigmaq
      complex, dimension(0:3,0:3,4,4) :: sigma
      complex fdir_1,fdir_2
      real pi,pm,dm
      real q0,qx,qy,qz,qv2,q2,F1s,F2s
      integer nuc,mu,j,i
      qv2 = qx**2 + qy**2 + qz**2
      q2  = qv2 - q0**2
      do mu = 0,3
      do j = 1,4
      do i = 1,4
      sigmaq(mu,i,j) = q0*sigma(mu,0,i,j) - qx*sigma(mu,1,i,j)
     &                -qy*sigma(mu,2,i,j) - qz*sigma(mu,3,i,j)

      enddo
      enddo
      enddo
     
      fdir_1 = cmplx(F1s(nuc,q2), 0.0)
      fdir_2 = cmplx(0.0       , 1.0/(2.0*pm)*F2s(nuc,q2))

      bgamma =fdir_1*gamma + fdir_2*sigmaq
*      print *,"sigma",fdir_2,q2,pm
      return
      end





     
**********************************************************
*  Dirac Nucleon Form Factors for Nucleon
**********************************************************
      function F1s(nuc,q2)
*      change
      if(nuc.eq. 1)F1s = F1ps(q2)
      if(nuc.eq.-1)F1s = F1ns(q2)
      return
      end

      function F2s(nuc,q2)
*       change
      if(nuc.eq. 1)F2s = F2ps(q2)
      if(nuc.eq.-1)F2s = F2ns(q2)
      return
      end


      function F1ps(Q2)                                                  
      common/parms/pi,pm,dm
      common/formfcase/ics
      call proton_formfactors(Q2,GEp,GMp,ics)
      tau = Q2/4./pm**2             
      F1ps = (tau*GMp+GEp)/(1.+tau)                                        
      return                                                            
      end                                                            
C     ----------------------------------------                          
      function F2ps(Q2)                                                  
      common/parms/pi,pm,dm
      common/formfcase/ics
      call proton_formfactors(Q2,GEp,GMp,ics)
      tau = Q2/4./pm**2                                                   
      F2ps = (GMp-GEp)/(1.+tau)    
      return                                             
      end                                                            
C     ---------------------------------------------                     
      function F1ns(Q2)                                                  
      common/parms/pi,pm,dm
      common/formfcase/ics
      call neutron_formfactors(Q2,GEn,GMn,ics)
      tau = Q2/4./pm**2                                         
      F1ns = (tau*GMn+GEn)/(1.+tau)                                        
      return                                                            
      end                                                            
C     ---------------------------------------------                     
      function F2ns(Q2)                                                  
      common/parms/pi,pm,dm
      common/formfcase/ics
      call neutron_formfactors(Q2,GEn,GMn,ics)
      tau = Q2/4./pm**2                                                   
      F2ns = (GMn-GEn)/(1.+tau)                                            
      return                                                            
      end                                                            
C     -------------------------------------                             
*******************************************************
*   Charge and Magnetic Form Factor Parameterizations
*   Standard parameterization (one I used before)
*******************************************************
      subroutine proton_formfactors(q2,GEp,GMp,ics)
      common/parms/pi,pm,dm
      tau = Q2/4./pm**2                           

      if(ics.eq.1)then
************************************************************
*  SLAC Parameterization
************************************************************
      GEp = Gs(Q2)                                                         
      GMp = Gs(Q2)*2.79

      elseif(ics.eq.2.or.ics.ge.4)then
************************************************************
*    J.J. Kelly's Parameterization
*    From Phys. Rev. C70, 068202, 2004
************************************************************
      GEp = (1.-0.24*tau)/(1.+10.98*tau+12.82*tau**2+21.97*tau**3)
      GMp = 2.79*(1.+0.12*tau)/(1.+10.97*tau+18.86*tau**2+6.55*tau**3)

      elseif(ics.eq.3)then
************************************************************
*    Bradford, Bodek, Budd, Arrington
*    Hep-ex/0602017
************************************************************
      GEp = (1.-0.0578*tau)/(1.+11.1*tau+13.6*tau**2+33.0*tau**3)
      GMp = 2.79*(1.+0.15*tau)/(1.+11.1*tau+19.6*tau**2+7.54*tau**3)
***************************************************************************      
      endif

      return
      end

      subroutine neutron_formfactors(q2,GEn,GMn,ics)
      common/parms/pi,pm,dm
      tau = Q2/4./pm**2    

      if(ics.eq.1)then
************************************************************
*  SLAC Parameterization
************************************************************
      GEn =  1.91*Gs(Q2)*tau/(1.+5.6*tau)                                   
      GMn = -1.91*Gs(Q2) 

      elseif(ics.eq.2.or.ics.ge.4)then
************************************************************
*    J.J. Kelly's Parameterization
*    From Phys. Rev. C70, 068202, 2004
************************************************************
      GEn0=  1.7*tau/(1.+3.3*tau)*Gs(Q2)                                  
      GMn = -1.91*(1.+2.33*tau)/(1.+14.72*tau+24.20*tau**2+84.1*tau**3)

      GEn = Gen0
      if(ics.eq.4)GEn = 0.0
      if(ics.eq.5)GEn = 2.0*GEn0

      elseif(ics.eq.3)then   
************************************************************
*    Bradford, Bodek, Budd, Arrington
*    Hep-ex/0602017
************************************************************
      GEn = (1.25*tau+1.3*tau**2)/(1.-9.86*tau+305.*tau**2-758.0*tau**3
     &                                                    +802.0*tau**4)

      GMn = -1.91*(1.+1.81*tau)/(1.+14.1*tau+20.7*tau**2+68.7*tau**3)
      endif

      return
      end



      function Gs(Q2)                                                    
      Gs=1./(1.+Q2/0.71)**2                                               
      return                                                            
      end                                                            





      subroutine us_s(m,px,py,pz,e,s,isp,us,usc)
*********************************************************
*     calculates Dirac spinor for
*     polarization axis defined by s
*     INPUTS
*     s - polarization axis three vector
*     isp = 1/-1 polarizationis along/opposite to s
*     OUTPUTS
*     us - dirac spinor, usc - complex conjucated us
*
*********************************************************
      implicit none
      complex, dimension(4)   :: us, usc
      real,    dimension(3)   :: s
      real,    dimension(0:3) :: s4
      real m,px,py,pz,e, s_dot_p
      complex p_pl,p_mn
      complex s_0,s_pl,s_mn,s_z
      integer isp
      real eplm,norm

      p_pl = cmplx(px,py)
      p_mn = cmplx(px,-py)
*      e = sqrt(m**2 + px**2 + py**2 + pz**2)
      eplm = e+m

      s_dot_p = s(1)*px + s(2)*py + s(3)*pz

      s4(0) = s_dot_p/m
      s4(1) = s(1) + s_dot_p/(m*(e+m)) * px
      s4(2) = s(2) + s_dot_p/(m*(e+m)) * py
      s4(3) = s(3) + s_dot_p/(m*(e+m)) * pz

      s_0   = cmplx(s4(0),0.0)
      s_pl  = cmplx(s4(1),s4(2))
      s_mn  = cmplx(s4(1),-s4(2))
      s_z   = cmplx(s4(3),0.0)
      norm = sqrt(2.0/(1.0+s(3)))*sqrt(eplm)/2.0
      if(isp.eq.1)then
      us(1)=norm*(1.0  + s_z - s_0*pz/eplm)
      us(2)=norm*(s_pl - s_0*p_pl/eplm)
      us(3)=norm*(s_0 + (1.0-s_z)*pz/eplm-s_mn*p_pl/eplm)
      us(4)=norm*(-s_pl*pz/eplm + (1.0 + s_z)*p_pl/eplm)

      usc(1)=norm*(1.0  + s_z - s_0*pz/eplm)
      usc(2)=norm*(s_mn - s_0*p_mn/eplm)
      usc(3)=norm*(s_0 +(1.0-s_z)*pz/eplm-s_pl*p_mn/eplm)
      usc(4)=norm*(-s_mn*pz/eplm + (1.0+s_z)*p_mn/eplm)


      elseif(isp.eq.-1.0)then
      us(1)=norm*(-s_mn + s_0*p_mn/eplm)
      us(2)=norm*(1.0 + s_z - s_0*pz/eplm)
      us(3)=norm*((1.0 +s_z)*p_mn/eplm  - s_mn*pz/eplm)
      us(4)=norm*(-s_0 + s_pl*p_mn/eplm - (1.0-s_z)*pz/eplm)

      usc(1)=norm*(-s_pl + s_0*p_pl/eplm)
      usc(2)=norm*(1.0 + s_z - s_0*pz/eplm)
      usc(3)=norm*((1.0 +s_z)*p_pl/eplm - s_pl*pz/eplm)
      usc(4)=norm*(-s_0 +s_mn*p_pl/eplm - (1.0-s_z)*pz/eplm)
      endif

      return
      end



      subroutine u_s(m,px,py,pz,e,isp,u,uc)
***************************************************************
*     This subroutine calculates Dirac Spinor
*     for polarization axis defined by z
*   m - mass
*   px,py,pz,e - three components of momenta and energy
*   isp - 1 spin - up  -1 spin down
*
*   u   - spinor
*   uc  - complex conjugate spinor
*
*
*
* Misak Sargsian
* March 2007, FIU, Miami
*
***************************************************************
      implicit none
      complex, dimension(4):: u,uc
      real m,px,py,pz,e,fct
      integer isp
      fct = sqrt(e+m)
      if(isp.eq.1)then
      u(1) = fct*cmplx(1.0,0.0)
      u(2) =     cmplx(0.0,0.0)
      u(3) = fct/(e+m)*cmplx(pz,0)
      u(4) = fct/(e+m)*cmplx(px,py)

      uc(1) = fct*cmplx(1.0,0.0)
      uc(2) =     cmplx(0.0,0.0)
      uc(3) = fct/(e+m)*cmplx(pz,0)
      uc(4) = fct/(e+m)*cmplx(px,-py)

      elseif(isp.eq.-1)then
      u(1) =     cmplx(0.0,0.0)
      u(2) = fct*cmplx(1.0,0.0)
      u(3) = fct/(e+m)*cmplx(px,-py)
      u(4) = fct/(e+m)*cmplx(-pz,0.0)

      uc(1) =     cmplx(0.0,0.0)
      uc(2) = fct*cmplx(1.0,0.0)
      uc(3) = fct/(e+m)*cmplx(px,py)
      uc(4) = fct/(e+m)*cmplx(-pz,0.0)
      endif
      return
      end


      


      subroutine e_vector(e_vec)
      implicit none
      complex, dimension(-1:1,3) :: e_vec
      complex im,one,zero
      zero = cmplx(0.0,0.0)
      one  = cmplx(1.0,0.0)      
      im   = cmplx(0.0,1.0)
      e_vec(-1,1:3) =  1.0/sqrt(2.0)*(/one, -im  ,zero/)
      e_vec(0,1:3)  =  1.0/sqrt(2.0)*(/zero, zero, one/)
      e_vec(1,1:3)  = -1.0/sqrt(2.0)*(/one,  im  ,zero/)
      return
      end




      subroutine clebsh_gordan(i1,mi1,i2,mi2,i,m,coeff)
*****************************************************************
*      Calculates the Clebsh-Gordan Coefficients
*      spins and projections defined as multiplied by 2
*      for example -1/2 is -1,  -3/2 is -3
*
*      Misak Sargsian
*      March 2007, FIU, Miami
*     
******************************************************************
      implicit none 
      integer i1,mi1,i2,mi2,i,m
      real coeff
      coeff = 0.0
      if(i1.eq.2.and.i2.eq.1.and.i.eq.3)then
      coeff = 0.0
      if(mi1.eq. 2.and.mi2.eq. 1.and.m.eq. 3)coeff = 1.0
      if(mi1.eq. 2.and.mi2.eq.-1.and.m.eq. 1)coeff = sqrt(1.0/3.0)
      if(mi1.eq. 0.and.mi2.eq. 1.and.m.eq. 1)coeff = sqrt(2.0/3.0)
      if(mi1.eq. 0.and.mi2.eq.-1.and.m.eq.-1)coeff = sqrt(2.0/3.0)
      if(mi1.eq.-2.and.mi2.eq. 1.and.m.eq.-1)coeff = sqrt(1.0/3.0)
      if(mi1.eq.-2.and.mi2.eq.-1.and.m.eq.-3)coeff = 1.0
      endif
      return
      end
      

      subroutine gammamatricess
*********************************************************************************
*     This subroutine calculates Dirac gamma and sigma 
*     matrices.
*
*    sigma_{\mu,\nu} = {1\over 2}[\gamma_\mu\gamma_nu - \gamma_\nu\gamma_\mu}
*    gamma(0:3,4,4) 0:3 - are the values of \mu
*    (4,4) components are such
*     a(4,4)
*     a11 a12 a13 a14
*     a21 a22 a23 a24
*     a31 a32 a33 a34
*     a41 a42 a43 a44      
*     gamma_pl(4,4) = gamma_0 + gamma3
*     gamma_mn(4,4) = gamma_0 - gamma3
*     gamma_pl and gamma_mn mattrices are added
*     March 2007, FIU, Miami
*********************************************************************************
      complex, dimension(0:3,4,4) :: gamma
      complex, dimension(4,4) :: gamma_pl,gamma_mn
      complex, dimension(0:3,0:3,4,4) :: sigma
      common/gammamatices_s/gamma,sigma
      common/gamma_lf/gamma_pl,gamma_mn
      complex fct
***********************************************
*  gamma_0
***********************************************
      gamma(0,1,1) = cmplx(1.,0.) 
      gamma(0,1,2) = cmplx(0.,0.)
      gamma(0,1,3) = cmplx(0.,0.)
      gamma(0,1,4) = cmplx(0.,0.)

      gamma(0,2,1) = cmplx(0.,0.) 
      gamma(0,2,2) = cmplx(1.,0.)
      gamma(0,2,3) = cmplx(0.,0.)
      gamma(0,2,4) = cmplx(0.,0.)

      gamma(0,3,1) = cmplx(0.,0.) 
      gamma(0,3,2) = cmplx(0.,0.)
      gamma(0,3,3) = cmplx(-1.,0.)
      gamma(0,3,4) = cmplx(0.,0.)

      gamma(0,4,1) = cmplx(0.,0.) 
      gamma(0,4,2) = cmplx(0.,0.)
      gamma(0,4,3) = cmplx(0.,0.)
      gamma(0,4,4) = cmplx(-1.,0.)
***********************************************
*  gamma_x
***********************************************
      gamma(1,1,1) = cmplx(0.,0.) 
      gamma(1,1,2) = cmplx(0.,0.)
      gamma(1,1,3) = cmplx(0.,0.)
      gamma(1,1,4) = cmplx(1.,0.)

      gamma(1,2,1) = cmplx(0.,0.) 
      gamma(1,2,2) = cmplx(0.,0.)
      gamma(1,2,3) = cmplx(1.,0.)
      gamma(1,2,4) = cmplx(0.,0.)

      gamma(1,3,1) = cmplx(0.,0.) 
      gamma(1,3,2) = cmplx(-1.,0.)
      gamma(1,3,3) = cmplx(0.,0.)
      gamma(1,3,4) = cmplx(0.,0.)

      gamma(1,4,1) = cmplx(-1.,0.) 
      gamma(1,4,2) = cmplx(0.,0.)
      gamma(1,4,3) = cmplx(0.,0.)
      gamma(1,4,4) = cmplx(0.,0.)
***********************************************
*  gamma_y
***********************************************
      gamma(2,1,1) = cmplx(0.,0.) 
      gamma(2,1,2) = cmplx(0.,0.)
      gamma(2,1,3) = cmplx(0.,0.)
      gamma(2,1,4) = cmplx(0.,-1.)

      gamma(2,2,1) = cmplx(0.,0.) 
      gamma(2,2,2) = cmplx(0.,0.)
      gamma(2,2,3) = cmplx(0.,1.)
      gamma(2,2,4) = cmplx(0.,0.)

      gamma(2,3,1) = cmplx(0.,0.) 
      gamma(2,3,2) = cmplx(0.,1.)
      gamma(2,3,3) = cmplx(0.,0.)
      gamma(2,3,4) = cmplx(0.,0.)

      gamma(2,4,1) = cmplx(0.,-1.) 
      gamma(2,4,2) = cmplx(0.,0.)
      gamma(2,4,3) = cmplx(0.,0.)
      gamma(2,4,4) = cmplx(0.,0.)

***********************************************
*  gamma_y
***********************************************
      gamma(3,1,1) = cmplx(0.,0.) 
      gamma(3,1,2) = cmplx(0.,0.)
      gamma(3,1,3) = cmplx(1.,0.)
      gamma(3,1,4) = cmplx(0.,0.)

      gamma(3,2,1) = cmplx(0.,0.) 
      gamma(3,2,2) = cmplx(0.,0.)
      gamma(3,2,3) = cmplx(0.,0.)
      gamma(3,2,4) = cmplx(-1.,0.)

      gamma(3,3,1) = cmplx(-1.,0.) 
      gamma(3,3,2) = cmplx(0.,0.)
      gamma(3,3,3) = cmplx(0.,0.)
      gamma(3,3,4) = cmplx(0.,0.)

      gamma(3,4,1) = cmplx(0.,0.) 
      gamma(3,4,2) = cmplx(1.,0.)
      gamma(3,4,3) = cmplx(0.,0.)
      gamma(3,4,4) = cmplx(0.,0.)

      fct = cmplx(0.0, 0.5)
      do mu = 0,3
      do nu = 0,3
      sigma(mu,nu,1:4,1:4) = fct*
     &  (matmul(gamma(mu,1:4,1:4),gamma(nu,1:4,1:4)) - 
     &   matmul(gamma(nu,1:4,1:4),gamma(mu,1:4,1:4)) )
      enddo
      enddo
      gamma_pl(1:4,1:4) = gamma(0,1:4,1:4) + gamma(3,1:4,1:4)
      gamma_mn(1:4,1:4) = gamma(0,1:4,1:4) - gamma(3,1:4,1:4)
      return
      end


      
************************************************************************************
*     This subroutine calculates nuclear masses though the Weizeker formula
************************************************************************************      
      
       subroutine nuclear_masses(a,z,am)
       implicit none
       real pmp,pmn
       real a,z,am
       real a1,a2,a3,a4,del_a
       integer in,ip
       pmp = 0.938272
       pmn = 0.939565
       a1 = 15.68
       a2 = 18.56
       a3 = 0.717
       a4 = 28.1
       del_a = 0.0
       
       in = (a-z)/2.0 
       ip = z/2.0
       if(float(in).lt.(a-z)/2.0.and.float(ip).lt.z/2.0)then
       del_a = 34.0
       elseif(float(in)+0.01.gt.(a-z)/2.0.and.float(ip)+0.01.gt.z/2.0)
     & then
       del_a = -34.0
       endif
*       print *,"del_a",del_a,a1

       am = z*pmp*1000. + (a-z)*pmn*1000. - a1*a 
     &      + a2 * a**(2.0/3.0) 
     &      + a3 *z**2/a**(1.0/3.0) + a4*(z - (a-z))**2/a 
     &      + del_a*a**(-3.0/4.0)
       am = am/1000.
       return
       end

	function sigma_en(px,py,pz,nuc,nsc)
	COMMON/ELECTRON_en/EI,ER,UE
      	COMMON/PHOTON_en/Q2,Q0,QV
	COMMON/MASS_en/PM,DM,PI 
	CHARACTER*7 NUCL
	PI = ACOS(-1.)
	PM = 0.938279
	DM = 1.875
	N  = nsc  !models
**	N = 2   !models
**	N = 3
*************************************************
* 3-Momentum of knocked out proton before the 
* rescattering
*************************************************
        pfx = px
        pfy = py
	pfz = pz + qv
	pf = sqrt(pfx**2+pfy**2+pfz**2)
*	pix = px
*	piy = py
*	piz = pfx 
	
	arg = pfz/pf
	if(arg.gt. 1.0)arg =  1.
	if(arg.lt.-1.0)arg = -1.
	thf   = acos(arg)
	sn_arg = sin(thf)

	phif = 0.0
	if(sn_arg.ne.0.0)then
	arg = pfx/pf/sn_arg
	if(arg.gt. 1.0)arg =  1.
	if(arg.lt.-1.0)arg = -1.
	phif = acos(arg)
	if(pfy/sn_arg.lt.0.0)phif = 2.0*pi - phif
	endif	

 	GAMMA_F = thf
	PHI_G   = phif
	P       = pf
***********************************************************
	NUCL    = 'PROTON '
	if(nuc.eq.-1)NUCL    = 'NEUTRON'
c*****************************************************************
c	Choice the model for G_en
c*****************************************************************
	IF(N.EQ.1)THEN
	CALL G_en_DF_1(PF,GAMMA_F,PHI_G,NUCL,W_C,W_T,W_I,W_S,G_en)
	ELSEIF(N.EQ.2)THEN
	CALL G_en_DF_2(PF,GAMMA_F,PHI_G,NUCL,W_C,W_T,W_I,W_S,G_en)
	ELSEIF(N.EQ.3.OR.N.EQ.5)THEN
	EF     = sqrt(pm**2 + pf**2)
	AL_F   = ( EF - PF*COS(GAMMA_F)) /PM
        PT     =        PF*SIN(GAMMA_F)
      	AL_q   = ( Q0 - QV ) / PM 
      	AL     = AL_F - AL_Q 
        DEL_M  = AL * ((PM**2 + PT**2)/AL + (PM**2 + PT**2)/(2.0 - AL) 
     &                                                    - DM**2/2.0)         
	IF(N.EQ.3)THEN
	CALL G_en_lc(PF,GAMMA_F,PHI_G,AL,AL_f,AL_q,DEL_M,
     &	                            NUCL,0,W_C,W_T,W_I,W_S,G_en)
	ELSEIF(N.EQ.5)THEN
	CALL G_en_FS(PF,GAMMA_F,PHI_G,AL,AL_f,AL_q,DEL_M,
     &	                            NUCL,0,W_C,W_T,W_I,W_S,G_en)
	ENDIF
	ELSEIF(N.EQ.4)THEN
	CALL G_en_S(PF,GAMMA_F,PHI_G,NUCL,W_C,W_T,W_I,W_S,G_en)
	ENDIF
	sigma_en =  G_EN                     
	RETURN
	END
********************************************************************************
*      Here are included the following  approximations for G_en                *
*  G_en_DF_1(PFIN,GAMMA,PHI,ANUN,W_c,W_t,W_i,W_s,G_en)                         *
*  G_en_DF_2(PFIN,GAMMA,PHI,ANUN,W_c,W_t,W_i,W_s,G_en)                         *
*  G_en_curr(PFIN,GAMMA,PHI,ANUN,W_c,W_t,W_i,W_s,G_en)                         *
*  G_en_FS  (PFIN,GAMMA,PHI,AL,AL_f,AL_q,DEL_M,ANUN,ILDA,W_c,W_t,W_i,W_s,G_en) *
*  G_en_lc  (PFIN,GAMMA,PHI,AL,AL_f,AL_q,DEL_M,ANUN,ILDA,W_c,W_t,W_i,W_s,G_en) *
*  G_en_S   (PFIN,GAMMA,PHI,ANUN,W_c,W_t,W_i,W_s,G_en)                         *
*                                                                              *
********************************************************************************

      	SUBROUTINE G_en_DF_1(PFIN,GAMMA,PHI,ANUN,W_c,W_t,W_i,W_s,G_en)
c****************************************************************
c	The following electromagnetic vertex for nucleon        *
c       have been taken ( De Forest approximation)              *
c            < J = U(G(F1+F2)-(P+P')F2/2M)U >                   *
c**************************************************************** 
	CHARACTER*7 ANUN
	COMMON/ELECTRON_en/EI,ER,UE
      	COMMON/PHOTON_en/Q2,Q0,QV
	COMMON/MASS_en/PM,DM,PI
	common/ggg/gm,TN
	A_c  = Q2**2 / QV**4
	A_t  = Q2/2.0/QV**2 + TAN(UE/2.0)**2
	A_i  = Q2/QV**2 * SQRT(Q2/QV**2 + TAN(UE/2.0)**2)
	A_s  = Q2/QV**2*COS(PHI)**2 + TAN(UE/2.0)**2
	A_tt = Q2/2.0/QV**2
	pin2 = PFIN**2 - 2.*PFIN*QV*COS(GAMMA) + QV**2
	if(pin2.le.0.0)pin2 = 0.0
      	PIN  = SQRT(pin2)
	print *,"piiin",pin
        EFIN = SQRT(PM**2 + PFIN**2)
	EIN  = SQRT(PM**2 + PIN**2)
	Q2_A = QV**2 - (EFIN-EIN)**2  
	IF(ANUN.EQ.'PROTON ')THEN
c************************************
c	Proton contribution         *
c************************************
	F1 = F1P(Q2)
	F2 = F2P(Q2)
	ELSEIF(ANUN.EQ.'NEUTRON')THEN
c************************************
c	Neutron contribution        *
c************************************
	F1 = F1N(Q2)
	F2 = F2N(Q2)
	ENDIF
c*********************************************************
c       Structure functions via de Forest formalism      *
c*********************************************************
	W_c = 1.0/(4.0*EIN*EFIN) * ( 
     &  (EIN+EFIN)**2 * (F1**2 + Q2_A/4.0/PM**2*F2**2)
     &                              - QV**2 * (F1+F2)**2)
c
	W_t = 1.0/(2.0*EIN*EFIN) * Q2_A * (F1+F2)**2
c
	W_s = 1.0/(EIN*EFIN) * PFIN**2*SIN(GAMMA)**2
     &                       * (F1**2 + Q2_A/4.0/PM**2*F2**2)
c
	W_i = - PFIN*SIN(GAMMA)/(EIN*EFIN) * (EIN+EFIN) 
     &                         * (F1**2 + Q2_A/4.0/PM**2*F2**2)
c*********************************************************
c       Elementary Cross section                         *
c*********************************************************
	GM   = GMOTT(EI,UE)
	G_en = GM * (A_c*W_c + A_t*W_t + A_i*W_i*COS(PHI) + A_s*W_s)
	RETURN
	END
c
      	SUBROUTINE G_en_DF_2(PFIN,GAMMA,PHI,ANUN,W_c,W_t,W_i,W_s,G_en)
c****************************************************************
c	The following electromagnetic vertex for nucleon        *
c       have been taken (De Forest approximation)               *
c            < J = U(G*F1-1/2*(Gq-qG)*F2/2M)U >                 *
c**************************************************************** 
	CHARACTER*7 ANUN
	COMMON/ELECTRON_en/EI,ER,UE
      	COMMON/PHOTON_en/Q2,Q0,QV
	COMMON/MASS_en/PM,DM,PI
	A_c  = Q2**2 / QV**4
	A_t  = Q2/2.0/QV**2 + TAN(UE/2.0)**2
	A_i  = Q2/QV**2 * SQRT(Q2/QV**2 + TAN(UE/2.0)**2)
	A_s  = Q2/QV**2*COS(PHI)**2 + TAN(UE/2.0)**2
	A_tt = Q2/2.0/QV**2
	pin2 = PFIN**2 - 2.*PFIN*QV*COS(GAMMA) + QV**2
	if(pin2.le.0.0)pin2 = 0.0
      	PIN  = SQRT(pin2)
        EFIN = SQRT(PM**2 + PFIN**2)
	EIN  = SQRT(PM**2 + PIN**2)
	Q2_A = QV**2 - (EFIN-EIN)**2  
      	PfP = EIN*EFIN - PFIN**2 + PFIN*QV*COS(GAMMA)      
	PfQ  = EFIN*Q0  - PFIN*QV*COS(GAMMA)
	PQ  = EIN *Q0  - PFIN*QV*COS(GAMMA) + QV**2
	IF(ANUN.EQ.'PROTON')THEN
c************************************
c	Proton contribution         *
c************************************
	F1 = F1P(Q2)
	F2 = F2P(Q2)
	ELSEIF(ANUN.EQ.'NEUTRON')THEN
c************************************
c	Neutron contribution        *
c************************************
	F1 = F1N(Q2)
	F2 = F2N(Q2)
	ENDIF
c*********************************************************
c       Structure functions via de Forest formalism      *
c*********************************************************
	W_c = 1.0/(EIN*EFIN) * (
     &	      F1**2*(EIN*EFIN + 1.0/2.0*(PM**2-PfP)) - 
     &        F1*F2*QV**2/2.0 +
     &        F2**2/4.0/PM**2 * (Q0*(EFIN*PQ+EIN*PfQ) +
     &	      Q2*EIN*EFIN - PQ*PfQ - QV**2/2.0*(PfP+PM**2)))
c
	W_t = 1.0/(EIN*EFIN) * (-F1**2*(PM**2-PfP) - 
     &        F1*F2*(PfQ-PQ) +
     &        F2**2/4.0/PM**2*(Q2*(PfP+PM**2) + 2.0*PQ*PfQ))
c
	W_s = PFIN**2*SIN(GAMMA)**2/(EIN*EFIN) * 
     &        ( F1**2 + Q2/4.0/PM**2*F2**2 )
c	
	W_i = - PFIN*SIN(GAMMA)/(EIN*EFIN) * ( (EFIN+EIN)*F1**2 +
     &          F2**2/4.0/PM**2 * ( Q2*(EIN+EFIN) + Q0*(PQ+PfQ)))
c*********************************************************
c       Elementary Cross section                         *
c*********************************************************
	GM   = GMOTT(EI,UE)
	G_en = GM * (A_c*W_c + A_t*W_t + A_i*W_i*COS(PHI) + A_s*W_s)
	RETURN
	END

      	SUBROUTINE G_en_curr(PFIN,GAMMA,PHI,ANUN,W_c,W_t,W_i,W_s,G_en)
c****************************************************************
c	The following electromagnetic vertex for nucleon        *
c       have been taken                                         *
c            < J = U(G*F1-1/2*(Gq-qG)*F2/2M)U >                 *
c       The current conservation not  used and                  *
c	the 0 and Z componenta calculated intependently         *
c**************************************************************** 
	CHARACTER*7 ANUN
	COMMON/ELECTRON_en/EI,ER,UE
      	COMMON/PHOTON_en/Q2,Q0,QV
	COMMON/MASS_en/PM,DM,PI
	A_c  = Q2**2 / QV**4
	A_t  = Q2/2.0/QV**2 + TAN(UE/2.0)**2
	A_i  = Q2/QV**2 * SQRT(Q2/QV**2 + TAN(UE/2.0)**2)
	A_s  = Q2/QV**2*COS(PHI)**2 + TAN(UE/2.0)**2
	A_tt = Q2/2.0/QV**2
c***********************************************************
c       De Forest off shell approximation                  *
c***********************************************************
      	PIN  = SQRT(PFIN**2 - 2.*PFIN*QV*COS(GAMMA) + QV**2)
        EFIN = SQRT(PM**2 + PFIN**2)
	EIN  = SQRT(PM**2 + PIN**2)
	Q2_A = QV**2 - (EFIN-EIN)**2  
      	PfP = EIN*EFIN - PFIN**2 + PFIN*QV*COS(GAMMA)      
	PfQ  = EFIN*Q0  - PFIN*QV*COS(GAMMA)
	PQ  = EIN *Q0  - PFIN*QV*COS(GAMMA) + QV**2
	IF(ANUN.EQ.'PROTON')THEN
c************************************
c	Proton contribution         *
c************************************
	F1 = F1P(Q2)
	F2 = F2P(Q2)
	ELSEIF(ANUN.EQ.'NEUTRON')THEN
c************************************
c	Neutron contribution        *
c************************************
	F1 = F1N(Q2)
	F2 = F2N(Q2)
	ENDIF
c*********************************************************
c       Structure functions                              *
c*********************************************************
	AA_C = 8.0*(EFIN*EIN*(Q2**2/QV**4) + 
     &         Q0/QV**2*(EFIN*PQ+EIN*PfQ)*Q2/QV**2 +
     &         Q0**2/QV**4 * PfQ*PQ - 1.0/2.0*(PfP-PM**2)*Q2/QV**2 )
	AB_C = 4.0 * (PfQ-PQ) * Q2/QV**2
	AG_C = 2.0*Q2**2/QV**4 *( Q0*(EFIN*PQ+EIN*PfQ) + Q2*EFIN*EIN -
     &         QV**2/2.0*(PM**2+PfP) - PQ*PfQ )

	W_c = 1.0/(8.0*EIN*EFIN) * QV**4/Q2**2 * 
     &            (F1**2*AA_C + F1*F2*AB_C + F2**2/PM**2*AG_C)

	W_t = 1.0/(EIN*EFIN) * (-F1**2*(PM**2-PfP) - 
     &        F1*F2*(PfQ-PQ) +
     &        F2**2/4.0/PM**2*(Q2*(PfP+PM**2) + 2.0*PQ*PfQ))

	W_s = PFIN**2*SIN(GAMMA)**2/(EIN*EFIN) * 
     &        ( F1**2 + Q2/4.0/PM**2*F2**2 )

	AA_I = - 4.0*PFIN*SIN(GAMMA) * 
     &           ( Q2/QV**2*(EIN+EFIN) + Q0/QV**2*(PQ+PfQ))
	AG_I = - Q2*PFIN*SIN(GAMMA) *
     &           ( Q2/QV**2*(EIN+EFIN) + Q0/QV**2*(PQ+PfQ))

	W_i = 1.0/(4.0*EIN*EFIN) * QV**2/Q2 *
     &        ( F1**2*AA_I + F2**2/PM**2*AG_I)
c*********************************************************
c       Elementary Cross section                         *
c*********************************************************
	GM   = GMOTT(EI,UE)
	G_en = GM * (A_c*W_c + A_t*W_t + A_i*W_i*COS(PHI) + A_s*W_s)
	RETURN
	END


      	SUBROUTINE G_en_FS(PFIN,GAMMA,PHI,AL,AL_f,AL_q,DEL_M,ANUN,
     &                ILDA,W_c,W_t,W_i,W_s,G_en)
c****************************************************************
c	The following electromagnetic vertex for nucleon        *
c       have been taken                                         *
c            < J = U(G*F1-1/2*(Gq-qG)*F2/2M)U >                 *
c       The current conservation used according of light-cone   *
c	formalism by which the    j+ = -(q+/q_)j_  (FS)         *
c**************************************************************** 
	CHARACTER*7 ANUN
	COMMON/ELECTRON_en/EI,ER,UE
      	COMMON/PHOTON_en/Q2,Q0,QV
	COMMON/MASS_en/PM,DM,PI
     	COMMON/SCALAR/PFD,PFD2,RHMA,PFQ                                         
      	EFIN = SQRT(PM**2 + PFIN**2)                                          
        EFIN = SQRT(PM**2 + PFIN**2)
      	COSD = Q0/QV                                                      
      	SIND = SQRT(Q2)/QV                                                
      	PFD  = EFIN-PFIN*COS(GAMMA)*COSD-PFIN*SIN(GAMMA)*SIND*COS(PHI)
      	PFD2 = PFD**2                                                       
	PFQ  = EFIN*Q0  - PFIN*QV*COS(GAMMA)
	RHMA = DEL_M
	PT   = PFIN * SIN(GAMMA)
	IF(ANUN.EQ.'PROTON')THEN
c************************************
c	Proton contribution         *
c************************************
	F1 = F1P(Q2)
	F2 = F2P(Q2)
	ELSEIF(ANUN.EQ.'NEUTRON')THEN
c************************************
c	Neutron contribution        *
c************************************
	F1 = F1N(Q2)
	F2 = F2N(Q2)
	ENDIF
      	GM=GMOTT(EI,UE)                                                   


	IF(ILDA.EQ.0)THEN
      	S = S1(AL,AL_f,AL_q,F1,F2) + S2(AL,AL_f,AL_q,F1,F2)
 
     	C1B = 1./3. * (QV**2/Q2 * S1(AL,AL_f,AL_q,F1,F2) +                    
     &                       2. * C1(AL,AL_f,AL_q,PT,F1,F2)) 
      	C2B = 1./3. * (QV**2/Q2 * S2(AL,AL_f,AL_q,F1,F2) +                    
     &                       2. * C2(AL,AL_f,AL_q,PT,F1,F2)) 

      	CB = C1B + C2B    
	ELSE
      	S = S1(AL,AL_f,AL_q,F1,F2) 
 
     	CB = 1./3. * (QV**2/Q2 * S1(AL,AL_f,AL_q,F1,F2) +                    
     &                       2. * C1(AL,AL_f,AL_q,PT,F1,F2)) 
	ENDIF                                             
c************************************
c	G_en                        *
c************************************
      	TNG2 = (SIN(UE/2.)/COS(UE/2.))**2                                   
      	G_en = GM/2./EFIN/(4.*PM)*(S + 2.*TNG2 *CB)                          
      	RETURN                                                            
      	END                                                               
C---------------------------------------------------                    
      FUNCTION S1(AL,ALF,AQ,F1,F2)                                        
      COMMON/PHOTON_en/GP2,GP0,GPV
      COMMON/MASS_en/PM,DM_M,PI                       
      COMMON/SCALAR/PFD,PFD2,RHMA,PFQ                                         
      S11 = F1**2*((PFD+1./GPV*(RHMA/2.))**2-GP2/GPV**2*RHMA/4.)    
      S12 = F2**2/4./PM**2*GP2*PFD2                                 
      S1  = 8.*(S11+S12)                                                  
      RETURN                                                            
      END                                                               
C.............................                                          
      FUNCTION S2(AL,ALF,AQ,F1,F2)                                        
      COMMON/PHOTON_en/GP2,GP0,GPV
      COMMON/MASS_en/PM,DM_M,PI                       
      COMMON/SCALAR/PFD,PFD2,RHMA,PFQ                                         
      S1I = F1**2*GP2*ALF*PM/GPV**2                                  
      S2I = F2*F1*(1.-GP0/GPV)*GP2/GPV                         
      S3I = F2**2/2./PM**2*(1.-GP0/GPV)*GP2*PFD                      
      S2  = 2.*RHMA/AL/PM*(S1I - S2I + S3I)     
      RETURN                                                            
      END                                                               
C.............................                                          
      FUNCTION C1(AL,ALF,AQ,PT,F1,F2)                                        
      COMMON/PHOTON_en/GP2,GP0,GPV
      COMMON/MASS_en/PM,DM_M,PI                       
      COMMON/SCALAR/PFD,PFD2,RHMA,PFQ                                         
      C1I = 4.*F1**2*(-ALF*AL/AQ**2*GP2+PT**2-PFQ)                   
      C2I = 6.*F1*F2*GP2                                       
      C3I = F2**2/2./PM**2*(GP2*(6.*PM**2+PFQ+2.*PT**2)-             
     *2.*(AL+ALF)/AQ*PFQ*GP2-2.*ALF**2/AQ**2*GP2**2+4.*PFQ**2)          
      C1  = ( C1I + C2I + C3I )        
      RETURN                                                            
      END                                                               
C......................                                                 
      FUNCTION C2(AL,ALF,AQ,PT,F1,F2)                                        
      COMMON/PHOTON_en/GP2,GP0,GPV
      COMMON/MASS_en/PM,DM_M,PI                       
      COMMON/SCALAR/PFD,PFD2,RHMA,PFQ                                         
      C1I = ALF*F1**2+AQ*F1*F2                           
      C2I = F2**2/4./PM**2*(2.*PFQ*AQ+GP2*ALF/2.)                    
      C2  = 2.0 * RHMA/AL * (C1I +C2I)
      RETURN                                                            
      END                                                               
C------------------------------------------                             

      	SUBROUTINE G_en_lc(PFIN,GAMMA,PHI,AL,AL_f,AL_q,DEL_M,ANUN,
     &                ILDA,W_c,W_t,W_i,W_s,G_en)
c****************************************************************
c	The following electromagnetic vertex for nucleon        *
c       have been taken                                         *
c            < J = U(G*F1-1/2*(Gq-qG)*F2/2M)U >                 *
c       The current conservation used according of light-cone   *
c	formalism by which the    j+ = -(q+/q_)j_               *
c**************************************************************** 
	CHARACTER*7 ANUN
	COMMON/ELECTRON_en/EI,ER,UE
      	COMMON/PHOTON_en/Q2,Q0,QV
	COMMON/MASS_en/PM,DM,PI
	A_c  = Q2**2 / QV**4
	A_t  = Q2/2.0/QV**2 + TAN(UE/2.0)**2
	A_i  = Q2/QV**2 * SQRT(Q2/QV**2 + TAN(UE/2.0)**2)
	A_s  = Q2/QV**2*COS(PHI)**2 + TAN(UE/2.0)**2
	A_tt = Q2/2.0/QV**2
	IF(ILDA.NE.0)DEL_M = 0.0
        EFIN = SQRT(PM**2 + PFIN**2)
	FK   = 8.0 * PM * EFIN
      	PfP  = (EFIN-Q0)*EFIN - PFIN**2 + PFIN*QV*COS(GAMMA) +     
     &         1.0/2.0 * (DEL_M/AL) * AL_f
	PfQ  = EFIN*Q0  - PFIN*QV*COS(GAMMA)
	PQ   = (EFIN-Q0)*Q0  - PFIN*QV*COS(GAMMA) + QV**2 + 
     &         1.0/2.0 * (DEL_M/AL) * AL_q
	IF(ANUN.EQ.'PROTON')THEN
c************************************
c	Proton contribution         *
c************************************
	F1 = F1P(Q2)
	F2 = F2P(Q2)
	ELSEIF(ANUN.EQ.'NEUTRON')THEN
c************************************
c	Neutron contribution        *
c************************************
	F1 = F1N(Q2)
	F2 = F2N(Q2)
	ENDIF
c*********************************************************
c       Structure functions                              *
c*********************************************************
	W_c = 8.0/FK * QV**2 * ( F1**2 * AL_f*AL/AL_q**2 - 
     &        F1*F2/2.0 * (AL_f-AL)/AL_q + 
     &        F2**2/4.0/PM**2 * ( -1.0/2.0 * (PM**2+PfP) + 
     &        Q2*AL_f*AL/AL_q**2 + PQ*AL_f/AL_q + PfQ*AL/AL_q))

	W_t = 8.0/FK * (-F1**2*(PM**2-PfP) - F1*F2*(PfQ-PQ) + 
     &        F2**2/4.0/PM**2 * (Q2*(PfP+PM**2) + 2.0*PQ*PfQ))

	W_s = 8.0*PFIN**2*SIN(GAMMA)**2/FK * (F1**2+Q2/4.0/PM**2*F2**2)

	W_i = 8.0/FK*QV*PFIN*SIN(GAMMA) * (F1**2*(AL_f+AL)/AL_q + 
     &        F2**2/4.0/PM**2*(Q2*(AL_f+AL)/AL_q + PQ + PfQ))
c*********************************************************
c       Elementary Cross section                         *
c*********************************************************
	GM   = GMOTT(EI,UE)
	G_en = GM * (A_c*W_c + A_t*W_t + A_i*W_i*COS(PHI) + A_s*W_s)
	RETURN
	END
C     ---------------------------------------------------               
      	SUBROUTINE G_en_S(PFIN,GAMMA,PHI,ANUN,W_c,W_t,W_i,W_s,G_en)
c****************************************************************
c	The following electromagnetic vertex for nucleon        *
c       have been taken (Saclay approximation)                  *
c            < J = U(G(F1+F2)-(P+P')F2/2M)U >                   *
c**************************************************************** 
	CHARACTER*7 ANUN
	COMMON/ELECTRON_en/EI,ER,UE
      	COMMON/PHOTON_en/Q2,Q0,QV
	COMMON/MASS_en/PM,DM,PI
	common/ggg/gm,TN
	A_c  = Q2**2 / QV**4
	A_t  = Q2/2.0/QV**2 + TAN(UE/2.0)**2
	A_i  = Q2/QV**2 * SQRT(Q2/QV**2 + TAN(UE/2.0)**2)
	A_s  = Q2/QV**2*COS(PHI)**2 + TAN(UE/2.0)**2
	A_tt = Q2/2.0/QV**2
      	PIN  = SQRT(PFIN**2 - 2.*PFIN*QV*COS(GAMMA) + QV**2)
        EFIN = SQRT(PM**2 + PFIN**2)
	EIN  = EFIN - Q0
	PMA = SQRT(EIN**2 - PIN**2)
	CS1  = COS(GAMMA) * (EI - ER*COS(UE)) / QV
	SN1  = SIN(GAMMA) *       ER*SIN(UE)  / QV
	CS2  = COS(GAMMA) * (EI*COS(UE) - ER) / QV
	SN2  = SIN(GAMMA) *  EI*SIN(UE)       / QV
	A3   = PFIN**2 * SN1 * SN2 / 2.0 
	A1   = (EFIN - PFIN*CS1)*(EFIN - PFIN*CS2) + A3
	A2   = (EFIN - PFIN*CS1)*PFIN*SN2 + PFIN*SN1*(EFIN-PFIN*CS2)  
	PfQ  = EFIN*Q0  - PFIN*QV*COS(GAMMA)
	IF(ANUN.EQ.'PROTON ')THEN
c************************************
c	Proton contribution         *
c************************************
	F1 = F1P(Q2)
	F2 = F2P(Q2)
	ELSEIF(ANUN.EQ.'NEUTRON')THEN
c************************************
c	Neutron contribution        *
c************************************
	F1 = F1N(Q2)
	F2 = F2N(Q2)
	ENDIF
c*********************************************************
c       Structure functions via de Forest formalism      *
c*********************************************************
	FAC   = 1.0 / (1.+ Q2/(PM+PMA)**2)
	FORM1 = (F1 - Q2/4.0/PM**2*F2)**2 
	FORM2 =  Q2/(PM+PMA)**2*(F1 + F2)**2
	FORM  =  form1 + form2
	ABC   = EFIN**2 - PFIN**2*COS(GAMMA)**2

	W2_IN = FAC/(EIN*EFIN) *
     &  (EFIN**2 - 2.*EFIN*PFIN*COS(GAMMA)*Q0/QV +
     &  PFIN**2/QV**2*(COS(GAMMA)**2*Q0**2+SIN(GAMMA)**2/2.*Q2))*FORM	

        V_T   = FAC/(EIN*EFIN) * ((ABC - PM**2)*FORM1 +
     &          (ABC + PM**2 + 2.*PfQ**2/Q2)*FORM2)
c
	W_i = 1.0 / (Q2/QV**2*SQRT(Q2/QV**2+TAN(UE/2.0)**2)) *
     &        FAC / (EIN*EFIN*COS(UE/2.0)**2) * (-A2*FORM)
c
	W_s = 1.0/(Q2/2.0/QV**2)*FAC/(EIN*EFIN*COS(UE/2.0)**2)*A3*FORM
c
	W_t = V_T - W_s
c
	W_c   =  QV**4/Q2**2 * (W2_in - Q2/2.0/QV**2*V_T)

c*********************************************************
c       Elementary Cross section                         *
c*********************************************************
	GM   = GMOTT(EI,UE)
	G_en = GM * (A_c*W_c + A_t*W_t + A_i*W_i*COS(PHI) + A_s*W_s)
	RETURN
	END


      FUNCTION F1P(Q2)                                                  
      COMMON/MASS_en/PM,DM,PI
      GMP=G(Q2)*2.79                                                    
      GEP=G(Q2)                                                         
      TAU=Q2/4./PM**2                                                   
      F1P=(TAU*GMP+GEP)/(1.+TAU)                                        
      RETURN                                                            
      END                                                               
C     ----------------------------------------                          
      FUNCTION F2P(Q2)                                                  
      COMMON/MASS_en/PM,DM,PI
      TAU=Q2/4./PM**2                                                   
      GMP=G(Q2)*2.79                                                    
      GEP=G(Q2)                                                         
      F2P=(GMP-GEP)/(1.+TAU)                                            
      RETURN                                                            
      END                                                               
C     ---------------------------------------------                     
      FUNCTION F1N(Q2)                                                  
      COMMON/MASS_en/PM,DM,PI
      TAU=Q2/4./PM**2                                                   
      GEN=1.91*G(Q2)*TAU/(1.+5.6*TAU)                                   
      GMN=-1.91*G(Q2)                                                   
      F1N=(TAU*GMN+GEN)/(1.+TAU)                                        
      RETURN                                                            
      END                                                               
C     ---------------------------------------------                     
      FUNCTION F2N(Q2)                                                  
      COMMON/MASS_en/PM,DM,PI
      TAU=Q2/4./PM**2                                                   
      GEN=1.91*G(Q2)*TAU/(1.+5.6*TAU)                                   
      GMN=-1.91*G(Q2)                                                   
      F2N=(GMN-GEN)/(1.+TAU)                                            
      RETURN                                                            
      END                                                               
C     -------------------------------------                             
      FUNCTION G(Q2)                                                    
      G=1./(1.+Q2/0.71)**2                                               
      RETURN                                                            
      END                                                               

 

C==== MOTT FACTOR=====                                                  
      FUNCTION GMOTT(EI,UE)                                             
c******************************************
c     Mott cross section in nanobarn 
c******************************************
      G1=(1./137)**2*COS(UE/2.)**2                                      
      G2=4.*EI**2*SIN(UE/2.)**4                                         
      GMOTT=G1/G2*0.389385*1000.*1000.                                  
      RETURN                                                            
      END                                                               

