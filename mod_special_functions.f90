! -----------------------------------------------------------------------------
! This file is derived from the original SPECIAL_FUNCTIONS source distributed at:
!   https://people.math.sc.edu/Burkardt/f_src/special_functions/special_functions.f90
!
! Original authors:
!   Shanjie Zhang, Jianming Jin (FORTRAN77 source routines)
!   John Burkardt (FORTRAN90 packaging/distribution)
!
! This local file has been modified for Modern Fortran style and maintenance.
! -----------------------------------------------------------------------------

module mod_special_functions

  use, intrinsic :: iso_fortran_env, only: real64, int32, error_unit

  implicit none
  private

  public :: special_functions

  type :: special_functions
  contains
    procedure, nopass, public :: airya, airyb, airyzo, ajyik, aswfa, &
      aswfb, bernoa, bernob, beta, bjndd, &
      cbk, cchg, cerf, cerror, cerzo, &
      cfc, cfs, cgama, ch12n, chgm, &
      chgu, chgubi, chguit, chgul, chgus, &
      cik01, ciklv, cikna, ciknb, cikva, &
      cikvb, cisia, cisib, cjk, cjy01, &
      cjylv, cjyna, cjynb, cjyva, cjyvb, &
      clpmn, clpn, clqmn, clqn, comelp, &
      cpbdn, cpdla, cpdsa, cpsi, csphik, &
      csphjy, cv0, cva1, cva2, cvf, &
      cvql, cvqm, cy01, cyzo, dvla, &
      dvsa, e1xa, e1xb, e1z, eix, &
      elit, elit3, envj, enxa, enxb, &
      error, eulera, eulerb, fcoef, fcs, &
      fcszo, ffk, gaih, gam0, gamma, &
      gmn, herzo, hygfx, hygfz, ik01a, &
      ik01b, ikna, iknb, ikv, incob, &
      incog, itairy, itika, itikb, itjya, &
      itjyb, itsh0, itsl0, itth0, ittika, &
      ittikb, ittjya, ittjyb, jdzo, jelp, &
      jy01a, jy01b, jyna, jynb, jyndd, &
      jyv, jyzo, klvna, klvnb, klvnzo, &
      kmn, lagzo, lamn, lamv, legzo, &
      lgama, lpmn, lpmns, lpmv, lpn, &
      lpni, lqmn, lqmns, lqna, lqnb, &
      msta1, msta2, mtu0, mtu12, othpl, &
      pbdv, pbvv, pbwa, psi, qstar, &
      r8_gamma_log, rctj, rcty, refine, rmn1, &
      rmn2l, rmn2so, rmn2sp, rswfo, rswfp, &
      scka, sckb, sdmn, segv, sphi, &
      sphj, sphk, sphy, stvh0, stvh1, &
      stvhv, stvl0, stvl1, stvlv, timestamp, &
      vvla, vvsa
  end type special_functions

contains

!> @brief Computes Airy functions and their derivatives.
!> @return The values of Ai(x), Bi(x), Ai'(x), Bi'(x).
!>
!> @param x [in] The argument of the Airy function.
!> @param ai [out] The value of Ai(x).
!> @param bi [out] The value of Bi(x).
!> @param ad [out] The value of Ai'(x).
!> @param bd [out] The value of Bi'(x).
subroutine airya ( x, ai, bi, ad, bd )

!*****************************************************************************80
!
!! AIRYA computes Airy functions and their derivatives.
!
!  Licensing:
!
!    The original FORTRAN77 version of this routine is copyrighted by 
!    Shanjie Zhang and Jianming Jin.  However, they give permission to 
!    incorporate this routine into a user program that the copyright 
!    is acknowledged.
!
!  Modified:
!
!    30 June 2012
!
!  Author:
!
!    Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!       
!  Parameters:
!
!    Input, real(real64) X, the argument of the Airy function.
!
!    Output, real(real64) AI, BI, AD, BD, the values of Ai(x), Bi(x),
!    Ai'(x), Bi'(x).
!
  implicit none

  real(real64), intent(inout) :: ad
  real(real64), intent(inout) :: ai
  real(real64), intent(inout) :: bd
  real(real64), intent(inout) :: bi
  real(real64) c1
  real(real64) c2
  real(real64) pir
  real(real64) sr3
  real(real64) vi1
  real(real64) vi2
  real(real64) vj1
  real(real64) vj2
  real(real64) vk1
  real(real64) vk2
  real(real64) vy1
  real(real64) vy2
  real(real64), intent(in) :: x
  real(real64) xa
  real(real64) xq
  real(real64) z

  xa = abs ( x )
  pir = 0.318309886183891e+00_real64
  c1 = 0.355028053887817e+00_real64
  c2 = 0.258819403792807e+00_real64
  sr3 = 1.732050807568877e+00_real64
  z = xa ** 1.5e+00_real64 / 1.5e+00_real64
  xq = sqrt ( xa )

  call ajyik ( z, vj1, vj2, vy1, vy2, vi1, vi2, vk1, vk2 )

  if ( x == 0.0e+00_real64 ) then
    ai = c1
    bi = sr3 * c1
    ad = - c2
    bd = sr3 * c2
  else if ( 0.0e+00_real64 < x ) then
    ai = pir * xq / sr3 * vk1
    bi = xq * ( pir * vk1 + 2.0e+00_real64 / sr3 * vi1 )
    ad = - xa / sr3 * pir * vk2
    bd = xa * ( pir * vk2 + 2.0e+00_real64 / sr3 * vi2 )
  else
    ai = 0.5e+00_real64 * xq * ( vj1 - vy1 / sr3 )
    bi = - 0.5e+00_real64 * xq * ( vj1 / sr3 + vy1 )
    ad = 0.5e+00_real64 * xa * ( vj2 + vy2 / sr3 )
    bd = 0.5e+00_real64 * xa * ( vj2 / sr3 - vy2 )
  end if

  return
end subroutine airya
!> @brief subroutine airyb.
!> @return None.
!>
!> @param x [in] Argument x.
!> @param ai [inout] Argument ai.
!> @param bi [inout] Argument bi.
!> @param ad [inout] Argument ad.
!> @param bd [inout] Argument bd.
pure subroutine airyb ( x, ai, bi, ad, bd )

!*****************************************************************************80
!
!! AIRYB computes Airy functions and their derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    02 June 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) X, argument of Airy function.
!
!    Output, real(real64) AI, Ai(x).
!
!    Output, real(real64) BI, Bi(x).
!
!    Output, real(real64) AD, Ai'(x).
!
!    Output, real(real64) BD, Bi'(x).
!
  implicit none

  real(real64), intent(inout) :: ad
  real(real64), intent(inout) :: ai
  real(real64), intent(inout) :: bd
  real(real64), intent(inout) :: bi
  real(real64) c1
  real(real64) c2
  real(real64) ck(41)
  real(real64) df
  real(real64) dg
  real(real64) dk(41)
  real(real64) eps
  real(real64) fx
  real(real64) gx
  integer(int32) k
  integer(int32) km
  real(real64) pi
  real(real64) r
  real(real64) rp
  real(real64) sad
  real(real64) sai
  real(real64) sbd
  real(real64) sbi
  real(real64) sda
  real(real64) sdb
  real(real64) sr3
  real(real64) ssa
  real(real64) ssb
  real(real64), intent(in) :: x
  real(real64) xa
  real(real64) xar
  real(real64) xcs
  real(real64) xe
  real(real64) xf
  real(real64) xm
  real(real64) xp1
  real(real64) xq
  real(real64) xr1
  real(real64) xr2
  real(real64) xss

  eps = 1.0e-15_real64
  pi = 3.141592653589793e+00_real64
  c1 = 0.355028053887817e+00_real64
  c2 = 0.258819403792807e+00_real64
  sr3 = 1.732050807568877e+00_real64
  xa = abs ( x )
  xq = sqrt ( xa )

  if ( x <= 0.0e+00_real64 ) then
    xm = 8.0e+00_real64
  else
    xm = 5.0e+00_real64
  end if

  if ( x == 0.0e+00_real64 ) then
    ai = c1
    bi = sr3 * c1
    ad = -c2
    bd = sr3 * c2
    return
  end if

  if ( xa <= xm ) then

    fx = 1.0e+00_real64
    r = 1.0e+00_real64
    do k = 1, 40
      r = r * x / ( 3.0e+00_real64 * k ) * x / ( 3.0e+00_real64 * k - 1.0e+00_real64 ) * x
      fx = fx + r
      if ( abs ( r ) < abs ( fx ) * eps ) then
        exit
      end if
    end do

    gx = x
    r = x
    do k = 1, 40
      r = r * x / ( 3.0e+00_real64 * k ) * x / ( 3.0e+00_real64 * k + 1.0e+00_real64 ) * x
      gx = gx + r
      if ( abs ( r ) < abs ( gx ) * eps ) then
        exit
      end if
    end do

    ai = c1 * fx - c2 * gx
    bi = sr3 * ( c1 * fx + c2 * gx )
    df = 0.5e+00_real64 * x * x
    r = df
    do k = 1, 40
      r = r * x / ( 3.0e+00_real64 * k ) * x / ( 3.0e+00_real64 * k + 2.0e+00_real64 ) * x
      df = df + r 
      if ( abs ( r ) < abs ( df ) * eps ) then
        exit
      end if
    end do

    dg = 1.0e+00_real64
    r = 1.0e+00_real64
    do k = 1, 40
      r = r * x / ( 3.0e+00_real64 * k ) * x / ( 3.0e+00_real64 * k - 2.0e+00_real64 ) * x
      dg = dg + r
      if ( abs ( r ) < abs ( dg ) * eps ) then
        exit
      end if
    end do

    ad = c1 * df - c2 * dg
    bd = sr3 * ( c1 * df + c2 * dg )

  else

    xe = xa * xq / 1.5e+00_real64
    xr1 = 1.0e+00_real64 / xe
    xar = 1.0e+00_real64 / xq
    xf = sqrt ( xar )
    rp = 0.5641895835477563e+00_real64
    r = 1.0e+00_real64
    do k = 1, 40
      r = r * ( 6.0e+00_real64 * k - 1.0e+00_real64 ) &
        / 216.0e+00_real64 * ( 6.0e+00_real64 * k - 3.0e+00_real64 ) &
        / k * ( 6.0e+00_real64 * k - 5.0e+00_real64 ) / ( 2.0e+00_real64 * k - 1.0e+00_real64 )
      ck(k) = r
      dk(k) = - ( 6.0e+00_real64 * k + 1.0e+00_real64 ) / ( 6.0e+00_real64 * k - 1.0e+00_real64 ) * ck(k)
    end do

    km = int ( 24.5e+00_real64 - xa )

    if ( xa < 6.0e+00_real64 ) then
      km = 14
    end if

    if ( 15.0e+00_real64 < xa ) then
      km = 10
    end if

    if ( 0.0e+00_real64 < x ) then
      sai = 1.0e+00_real64
      sad = 1.0e+00_real64
      r = 1.0e+00_real64
      do k = 1, km
        r = - r * xr1
        sai = sai + ck(k) * r
        sad = sad + dk(k) * r
      end do
      sbi = 1.0e+00_real64
      sbd = 1.0e+00_real64
      r = 1.0e+00_real64
      do k = 1, km
        r = r * xr1
        sbi = sbi + ck(k) * r
        sbd = sbd + dk(k) * r
      end do
      xp1 = exp ( - xe )
      ai = 0.5e+00_real64 * rp * xf * xp1 * sai
      bi = rp * xf / xp1 * sbi
      ad = -0.5e+00_real64 * rp / xf * xp1 * sad
      bd = rp / xf / xp1 * sbd
    else
      xcs = cos ( xe + pi / 4.0e+00_real64 )
      xss = sin ( xe + pi / 4.0e+00_real64 )
      ssa = 1.0e+00_real64
      sda = 1.0e+00_real64
      r = 1.0e+00_real64
      xr2 = 1.0e+00_real64 / ( xe * xe )
      do k = 1, km
        r = - r * xr2
        ssa = ssa + ck(2*k) * r
        sda = sda + dk(2*k) * r
      end do
      ssb = ck(1) * xr1
      sdb = dk(1) * xr1
      r = xr1
      do k = 1, km
        r = - r * xr2
        ssb = ssb + ck(2*k+1) * r
        sdb = sdb + dk(2*k+1) * r
      end do
      ai = rp * xf * ( xss * ssa - xcs * ssb )
      bi = rp * xf * ( xcs * ssa + xss * ssb )
      ad = -rp / xf * ( xcs * sda + xss * sdb )
      bd =  rp / xf * ( xss * sda - xcs * sdb )
    end if

  end if

  return
end subroutine airyb
!> @brief subroutine airyzo.
!> @return None.
!>
!> @param nt [in] Argument nt.
!> @param kf [in] Argument kf.
!> @param xa [inout] Argument xa.
!> @param xb [inout] Argument xb.
!> @param xc [inout] Argument xc.
!> @param xd [inout] Argument xd.
subroutine airyzo ( nt, kf, xa, xb, xc, xd )

!*****************************************************************************80
!
!! AIRYZO computes the first NT zeros of Ai(x) and Ai'(x).
!
!   Discussion:
!
!    Compute the first NT zeros of Airy functions Ai(x) and Ai'(x), 
!    a and a', and the associated values of Ai(a') and Ai'(a); and 
!    the first NT zeros of Airy functions Bi(x) and Bi'(x), b and
!    b', and the associated values of Bi(b') and Bi'(b).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    14 March 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) NT, the number of zeros.
!
!    Input, integer(int32) KF, the function code.
!    1 for Ai(x) and Ai'(x);
!    2 for Bi(x) and Bi'(x).
!
!    Output, real(real64) XA(m), a, the m-th zero of Ai(x) or
!    b, the m-th zero of Bi(x).
!
!    Output, real(real64) XB(m), a', the m-th zero of Ai'(x) or
!    b', the m-th zero of Bi'(x).
!
!    Output, real(real64) XC(m), Ai(a') or Bi(b').
!
!    Output, real(real64) XD(m), Ai'(a) or Bi'(b)
!
  implicit none

  integer(int32), intent(in) :: nt

  real(real64) ad
  real(real64) ai
  real(real64) bd
  real(real64) bi
  integer(int32) i
  integer(int32), intent(in) :: kf
  real(real64) pi
  real(real64) rt
  real(real64) rt0
  real(real64) u
  real(real64) u1
  real(real64) x
  real(real64), intent(inout) :: xa(nt)
  real(real64), intent(inout) :: xb(nt)
  real(real64), intent(inout) :: xc(nt)
  real(real64), intent(inout) :: xd(nt)

  pi = 3.141592653589793e+00_real64

  do i = 1, nt

    if (kf == 1) then
      u = 3.0e+00_real64 * pi * ( 4.0e+00_real64 * i - 1 ) / 8.0e+00_real64
      u1 = 1.0e+00_real64 / ( u * u )
      rt0 = - ( u * u ) ** ( 1.0_real64 / 3.0_real64 ) &
        * (((( -15.5902e+00_real64 * u1 + 0.929844e+00_real64 ) * u1 &
        - 0.138889e+00_real64 ) * u1 + 0.10416667e+00_real64 ) * u1 + 1.0e+00_real64 )
    else if ( kf == 2 ) then
      if ( i == 1 ) then
        rt0 = -1.17371e+00_real64
      else
        u = 3.0e+00_real64 * pi * ( 4.0e+00_real64 * i - 3.0e+00_real64 ) / 8.0e+00_real64
        u1 = 1.0e+00_real64 / ( u * u )
        rt0 = - ( u * u ) ** ( 1.0e+00_real64 / 3.0e+00_real64 ) &
          * (((( -15.5902e+00_real64 * u1 + 0.929844e+00_real64 ) * u1 &
          - 0.138889e+00_real64 ) * u1 + 0.10416667e+00_real64 ) * u1 + 1.0e+00_real64 )
      end if
    end if

    do

      x = rt0
      call airyb ( x, ai, bi, ad, bd )

      if ( kf == 1 ) then
        rt = rt0 - ai / ad
      else
        rt = rt0 - bi / bd
      end if

      if ( abs ( ( rt - rt0 ) / rt ) <= 1.0e-09_real64 ) then
        exit
      end if
      rt0 = rt

    end do

    xa(i) = rt
    if ( kf == 1 ) then
      xd(i) = ad
    else
      xd(i) = bd
    end if

  end do

  do i = 1, nt

    if ( kf == 1 ) then
      if ( i == 1 ) then
        rt0 = -1.01879e+00_real64
      else
        u = 3.0e+00_real64 * pi * ( 4.0e+00_real64 * i - 3.0e+00_real64 ) / 8.0e+00_real64
        u1 = 1.0e+00_real64 / ( u * u )
        rt0 = - ( u * u ) ** ( 1.0e+00_real64 / 3.0e+00_real64 ) &
          * (((( 15.0168e+00_real64 * u1 - 0.873954e+00_real64 ) &
          * u1 + 0.121528e+00_real64 ) * u1 - 0.145833e+00_real64 ) * u1 + 1.0e+00_real64 )
      end if
    else if ( kf == 2 ) then
      if ( i == 1 ) then
        rt0 = -2.29444e+00_real64
      else
        u = 3.0e+00_real64 * pi * ( 4.0e+00_real64 * i - 1.0e+00_real64 ) / 8.0e+00_real64
        u1 = 1.0e+00_real64 / ( u * u )
        rt0 = - ( u * u ) ** ( 1.0e+00_real64 / 3.0e+00_real64 ) &
          * (((( 15.0168e+00_real64 * u1 - 0.873954e+00_real64 ) &
          * u1 + 0.121528e+00_real64 ) * u1 - 0.145833e+00_real64 ) * u1 + 1.0e+00_real64 )
      end if
    end if

    do

      x = rt0
      call airyb ( x, ai, bi, ad, bd )

      if ( kf == 1 ) then
        rt = rt0 - ad / ( ai * x )
      else
        rt = rt0 - bd / ( bi * x )
      end if

      if ( abs ( ( rt - rt0 ) / rt ) <= 1.0e-09_real64 ) then
        exit
      end if

      rt0 = rt

    end do
    
    xb(i) = rt
    if ( kf == 1 ) then
      xc(i) = ai
    else
      xc(i) = bi
    end if

  end do

  return
end subroutine airyzo
!> @brief subroutine ajyik.
!> @return None.
!>
!> @param x [in] Argument x.
!> @param vj1 [inout] Argument vj1.
!> @param vj2 [inout] Argument vj2.
!> @param vy1 [inout] Argument vy1.
!> @param vy2 [inout] Argument vy2.
!> @param vi1 [inout] Argument vi1.
!> @param vi2 [inout] Argument vi2.
!> @param vk1 [inout] Argument vk1.
!> @param vk2 [inout] Argument vk2.
pure subroutine ajyik ( x, vj1, vj2, vy1, vy2, vi1, vi2, vk1, vk2 )

!*****************************************************************************80
!
!! AJYIK computes Bessel functions Jv(x), Yv(x), Iv(x), Kv(x).
!
!  Discussion: 
!
!    Compute Bessel functions Jv(x) and Yv(x), and modified Bessel functions 
!    Iv(x) and Kv(x), and their derivatives with v = 1/3, 2/3._real64
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    31 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) X, the argument.  X should not be zero.
!
!    Output, real(real64) VJ1, VJ2, VY1, VY2, VI1, VI2, VK1, VK2,
!    the values of J1/3(x), J2/3(x), Y1/3(x), Y2/3(x), I1/3(x), I2/3(x),
!    K1/3(x), K2/3(x).
!
  implicit none

  real(real64) a0
  real(real64) b0
  real(real64) c0
  real(real64) ck
  real(real64) gn
  real(real64) gn1
  real(real64) gn2
  real(real64) gp1
  real(real64) gp2
  integer(int32) k
  integer(int32) k0
  integer(int32) l
  real(real64) pi
  real(real64) pv1
  real(real64) pv2
  real(real64) px
  real(real64) qx
  real(real64) r
  real(real64) rp
  real(real64) rp2
  real(real64) rq
  real(real64) sk
  real(real64) sum
  real(real64) uj1
  real(real64) uj2
  real(real64) uu0
  real(real64), intent(inout) :: vi1
  real(real64), intent(inout) :: vi2
  real(real64) vil
  real(real64), intent(inout) :: vj1
  real(real64), intent(inout) :: vj2
  real(real64) vjl
  real(real64), intent(inout) :: vk1
  real(real64), intent(inout) :: vk2
  real(real64) vl
  real(real64) vsl
  real(real64) vv
  real(real64) vv0
  real(real64), intent(inout) :: vy1
  real(real64), intent(inout) :: vy2
  real(real64), intent(in) :: x
  real(real64) x2
  real(real64) xk

  if ( x == 0.0e+00_real64 ) then
    vj1 = 0.0e+00_real64
    vj2 = 0.0e+00_real64
    vy1 = -1.0e+300_real64
    vy2 = 1.0e+300_real64
    vi1 = 0.0e+00_real64
    vi2 = 0.0e+00_real64
    vk1 = -1.0e+300_real64
    vk2 = -1.0e+300_real64
    return
  end if

  pi = 3.141592653589793e+00_real64
  rp2 = 0.63661977236758e+00_real64
  gp1 = 0.892979511569249e+00_real64
  gp2 = 0.902745292950934e+00_real64
  gn1 = 1.3541179394264e+00_real64
  gn2 = 2.678938534707747e+00_real64
  vv0 = 0.444444444444444e+00_real64
  uu0 = 1.1547005383793e+00_real64
  x2 = x * x

  if ( x < 35.0e+00_real64 ) then
    k0 = 12
  else if ( x < 50.0e+00_real64 ) then
    k0 = 10
  else
    k0 = 8
  end if

  if ( x <= 12.0e+00_real64 ) then

    do l = 1, 2
      vl = l / 3.0e+00_real64
      vjl = 1.0e+00_real64
      r = 1.0e+00_real64
      do k = 1, 40
        r = -0.25e+00_real64 * r * x2 / ( k * ( k + vl ) )
        vjl = vjl + r
        if ( abs ( r ) < 1.0e-15_real64 ) then
          exit
        end if
      end do

      a0 = ( 0.5e+00_real64 * x ) ** vl
      if ( l == 1 ) then
        vj1 = a0 / gp1 * vjl
      else
        vj2 = a0 / gp2 * vjl
      end if

    end do

  else

    do l = 1, 2

      vv = vv0 * l * l
      px = 1.0e+00_real64
      rp = 1.0e+00_real64

      do k = 1, k0
        rp = - 0.78125e-02_real64 * rp &
          * ( vv - ( 4.0e+00_real64 * k - 3.0e+00_real64 ) ** 2 ) &
          * ( vv - ( 4.0e+00_real64 * k - 1.0e+00_real64 ) ** 2 ) &
          / ( k * ( 2.0e+00_real64 * k - 1.0e+00_real64 ) * x2 )
        px = px + rp
      end do

      qx = 1.0e+00_real64
      rq = 1.0e+00_real64
      do k = 1, k0
        rq = - 0.78125e-02_real64 * rq &
          * ( vv - ( 4.0e+00_real64 * k - 1.0e+00_real64 ) ** 2 ) &
          * ( vv - ( 4.0e+00_real64 * k + 1.0e+00_real64 ) ** 2 ) &
          / ( k * ( 2.0e+00_real64 * k + 1.0e+00_real64 ) * x2 )
        qx = qx + rq
      end do

      qx = 0.125e+00_real64 * ( vv - 1.0e+00_real64 ) * qx / x
      xk = x - ( 0.5e+00_real64 * l / 3.0e+00_real64 + 0.25e+00_real64 ) * pi
      a0 = sqrt ( rp2 / x )
      ck = cos ( xk )
      sk = sin ( xk )
      if ( l == 1) then
        vj1 = a0 * ( px * ck - qx * sk )
        vy1 = a0 * ( px * sk + qx * ck )
      else
        vj2 = a0 * ( px * ck - qx * sk )
        vy2 = a0 * ( px * sk + qx * ck )
      end if

    end do

  end if

  if ( x <= 12.0e+00_real64 ) then

    do l = 1, 2

      vl = l / 3.0e+00_real64
      vjl = 1.0e+00_real64
      r = 1.0e+00_real64
      do k = 1, 40
        r = -0.25e+00_real64 * r * x2 / ( k * ( k - vl ) )
        vjl = vjl + r
        if ( abs ( r ) < 1.0e-15_real64 ) then
          exit
        end if
      end do

      b0 = ( 2.0e+00_real64 / x ) ** vl
      if ( l == 1 ) then
        uj1 = b0 * vjl / gn1
      else
         uj2 = b0 * vjl / gn2
      end if

    end do

    pv1 = pi / 3.0e+00_real64
    pv2 = pi / 1.5e+00_real64
    vy1 = uu0 * ( vj1 * cos ( pv1 ) - uj1 )
    vy2 = uu0 * ( vj2 * cos ( pv2 ) - uj2 )

  end if

  if ( x <= 18.0e+00_real64 ) then

    do l = 1, 2
      vl = l / 3.0e+00_real64
      vil = 1.0e+00_real64
      r = 1.0e+00_real64
      do k = 1, 40
        r = 0.25e+00_real64 * r * x2 / ( k * ( k + vl ) )
        vil = vil + r
        if ( abs ( r ) < 1.0e-15_real64 ) then
          exit
        end if
      end do

      a0 = ( 0.5e+00_real64 * x ) ** vl

      if ( l == 1 ) then
        vi1 = a0 / gp1 * vil
      else
        vi2 = a0 / gp2 * vil
      end if

    end do

  else

    c0 = exp ( x ) / sqrt ( 2.0e+00_real64 * pi * x )

    do l = 1, 2
      vv = vv0 * l * l
      vsl = 1.0e+00_real64
      r = 1.0e+00_real64
      do k = 1, k0
        r = - 0.125e+00_real64 * r &
          * ( vv - ( 2.0e+00_real64 * k - 1.0e+00_real64 ) ** 2 ) / ( k * x )
        vsl = vsl + r
      end do
      if ( l == 1 ) then
        vi1 = c0 * vsl
      else
        vi2 = c0 * vsl
      end if
    end do

  end if

  if ( x <= 9.0e+00_real64 ) then

    do l = 1, 2
      vl = l / 3.0e+00_real64
      if ( l == 1 ) then
        gn = gn1
      else
        gn = gn2
      end if
      a0 = ( 2.0e+00_real64 / x ) ** vl / gn
      sum = 1.0e+00_real64
      r = 1.0e+00_real64
      do k = 1, 60
        r = 0.25e+00_real64 * r * x2 / ( k * ( k - vl ) )
        sum = sum + r
        if ( abs ( r ) < 1.0e-15_real64 ) then
          exit
        end if
      end do

      if ( l == 1 ) then
        vk1 = 0.5e+00_real64 * uu0 * pi * ( sum * a0 - vi1 )
      else
        vk2 = 0.5e+00_real64 * uu0 * pi * ( sum * a0 - vi2 )
      end if

    end do

  else

    c0 = exp ( - x ) * sqrt ( 0.5e+00_real64 * pi / x )

    do l = 1, 2
      vv = vv0 * l * l
      sum = 1.0e+00_real64
      r = 1.0e+00_real64
      do k = 1, k0
        r = 0.125e+00_real64 * r * ( vv - ( 2.0e+00_real64 * k - 1.0e+00_real64 ) ** 2 ) / ( k * x )
        sum = sum + r
      end do
      if ( l == 1 ) then
        vk1 = c0 * sum
      else
        vk2 = c0 * sum
      end if
    end do

  end if

  return
end subroutine ajyik
!> @brief subroutine aswfa.
!> @return None.
!>
!> @param m [in] Argument m.
!> @param n [in] Argument n.
!> @param c [in] Argument c.
!> @param x [inout] Argument x.
!> @param kd [in] Argument kd.
!> @param cv [in] Argument cv.
!> @param s1f [inout] Argument s1f.
!> @param s1d [inout] Argument s1d.
subroutine aswfa ( m, n, c, x, kd, cv, s1f, s1d )

!*****************************************************************************80
!
!! ASWFA: prolate and oblate spheroidal angular functions of the first kind.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    13 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) M, the mode parameter.
!
!    Input, integer(int32) N, the mode parameter, with N = M, M+1, ...
!
!    Input, real(real64) C, the spheroidal parameter.
!
!    Input, real(real64) X, the argument of the angular function.
!    |X| < 1.0_real64.
!
!    Input, integer(int32) KD, the function code.
!    1, the prolate function.
!    -1, the oblate function.
!
!    Input, real(real64) CV, the characteristic value.
!
!    Output, real(real64) S1F, S1D, the angular function of the first
!    kind and its derivative.
!
  implicit none

  real(real64) a0
  real(real64), intent(in) :: c
  real(real64) ck(200)
  real(real64), intent(in) :: cv
  real(real64) d0
  real(real64) d1
  real(real64) df(200)
  real(real64) eps
  integer(int32) ip
  integer(int32) k
  integer(int32), intent(in) :: kd
  integer(int32), intent(in) :: m
  integer(int32), intent(in) :: n
  integer(int32) nm
  integer(int32) nm2
  real(real64) r
  real(real64), intent(inout) :: s1d
  real(real64), intent(inout) :: s1f
  real(real64) su1
  real(real64) su2
  real(real64), intent(inout) :: x
  real(real64) x0
  real(real64) x1

  eps = 1.0e-14_real64
  x0 = x
  x = abs ( x )

  if ( n - m == 2 * int ( ( n - m ) / 2 ) ) then
    ip = 0
  else
    ip = 1
  end if

  nm = 10 + int ( ( n - m ) / 2 + c )
  nm2 = nm / 2 - 2 
  call sdmn ( m, n, c, cv, kd, df )
  call sckb ( m, n, c, df, ck )
  x1 = 1.0e+00_real64 - x * x

  if ( m == 0 .and. x1 == 0.0e+00_real64 ) then
    a0 = 1.0e+00_real64
  else
    a0 = x1 ** ( 0.5e+00_real64 * m )
  end if

  su1 = ck(1)
  do k = 1, nm2
    r = ck(k+1) * x1 ** k
    su1 = su1 + r
    if ( 10 <= k .and. abs ( r / su1 ) < eps ) then
      exit
    end if
  end do

  s1f = a0 * x ** ip * su1

  if ( x == 1.0e+00_real64 ) then

    if ( m == 0 ) then
      s1d = ip * ck(1) - 2.0e+00_real64 * ck(2)
    else if ( m == 1 ) then
      s1d = -1.0e+100_real64
    else if ( m == 2 ) then
      s1d = -2.0e+00_real64 * ck(1)
    else if ( 3 <= m ) then
      s1d = 0.0e+00_real64
    end if

  else

    d0 = ip - m / x1 * x ** ( ip + 1.0e+00_real64 )
    d1 = -2.0e+00_real64 * a0 * x ** ( ip + 1.0e+00_real64 )
    su2 = ck(2)
    do k = 2, nm2
      r = k * ck(k+1) * x1 ** ( k - 1.0e+00_real64 )
      su2 = su2 + r
      if ( 10 <= k .and. abs ( r / su2 ) < eps ) then
        exit
      end if
    end do

    s1d = d0 * a0 * su1 + d1 * su2

  end if

  if ( x0 < 0.0e+00_real64 ) then
    if ( ip == 0 ) then
      s1d = -s1d
    else if ( ip == 1 ) then
      s1f = -s1f
    end if
  end if

  x = x0

  return
end subroutine aswfa
!> @brief subroutine aswfb.
!> @return None.
!>
!> @param m [in] Argument m.
!> @param n [in] Argument n.
!> @param c [in] Argument c.
!> @param x [in] Argument x.
!> @param kd [in] Argument kd.
!> @param cv [in] Argument cv.
!> @param s1f [inout] Argument s1f.
!> @param s1d [inout] Argument s1d.
subroutine aswfb ( m, n, c, x, kd, cv, s1f, s1d )

!*****************************************************************************80
!
!! ASWFB: prolate and oblate spheroidal angular functions of the first kind.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    20 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) M, the mode parameter, m = 0, 1, 2, ...
!
!    Input, integer(int32) N, mode parameter, N = M, M+1, M+2, ...
!
!    Input, real(real64) C, the spheroidal parameter.
!
!    Input, real(real64) X, the argument, with |X| < 1.0_real64.
!
!    Input, integer(int32) KD, the function code.
!    1, the prolate function.
!    -1, the oblate function.
!
!    Input, real(real64) CV, the characteristic value.
!
!    Output, real(real64) S1F, S1D, the angular function of the first
!    kind and its derivative.
!
  implicit none

  real(real64), intent(in) :: c
  real(real64), intent(in) :: cv
  real(real64) df(200)
  real(real64) eps
  integer(int32) ip
  integer(int32) k
  integer(int32), intent(in) :: kd
  integer(int32), intent(in) :: m
  integer(int32) mk
  integer(int32), intent(in) :: n
  integer(int32) nm
  integer(int32) nm2
  real(real64) pd(0:251)
  real(real64) pm(0:251)
  real(real64), intent(inout) :: s1d
  real(real64), intent(inout) :: s1f
  real(real64) su1
  real(real64) sw
  real(real64), intent(in) :: x

  eps = 1.0e-14_real64

  if ( n - m == 2 * int ( ( n - m ) / 2 ) ) then
    ip = 0
  else
    ip = 1
  end if

  nm = 25 + int ( ( n - m ) / 2 + c )
  nm2 = 2 * nm + m
  call sdmn ( m, n, c, cv, kd, df )
  call lpmns ( m, nm2, x, pm, pd )
  su1 = 0.0e+00_real64
  sw = 0.0e+00_real64
  do k = 1, nm
    mk = m + 2 * ( k - 1 ) + ip
    su1 = su1 + df(k) * pm(mk)
    if ( abs ( sw - su1 ) < abs ( su1 ) * eps ) then
      exit
    end if
    sw = su1
  end do

  s1f = ( -1.0e+00_real64 ) ** m * su1

  su1 = 0.0e+00_real64
  sw = 0.0e+00_real64
  do k = 1, nm
    mk = m + 2 * ( k - 1 ) + ip
    su1 = su1 + df(k) * pd(mk)
    if ( abs ( sw - su1 ) < abs ( su1 ) * eps ) then
      exit
    end if
    sw = su1
  end do

  s1d = ( -1.0e+00_real64 ) ** m * su1

  return
end subroutine aswfb
!> @brief subroutine bernoa.
!> @return None.
!>
!> @param n [in] Argument n.
!> @param bn [inout] Argument bn.
pure subroutine bernoa ( n, bn )

!*****************************************************************************80
!
!! BERNOA computes the Bernoulli number Bn.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    11 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) N, the index.
!
!    Output, real(real64) BN, the value of the N-th Bernoulli number.
!
  implicit none

  integer(int32), intent(in) :: n

  real(real64), intent(inout) :: bn(0:n)
  integer(int32) j
  integer(int32) k
  integer(int32) m
  real(real64) r
  real(real64) s

  bn(0) = 1.0e+00_real64
  bn(1) = -0.5e+00_real64

  do m = 2, n
    s = - ( 1.0e+00_real64 / ( m + 1.0e+00_real64 ) - 0.5e+00_real64 )
    do k = 2, m - 1
      r = 1.0e+00_real64
      do j = 2, k
        r = r * ( j + m - k ) / j
      end do
    s = s - r * bn(k)
   end do
   bn(m) = s
  end do

  do m = 3, n, 2
    bn(m) = 0.0e+00_real64
  end do

  return
end subroutine bernoa
!> @brief subroutine bernob.
!> @return None.
!>
!> @param n [in] Argument n.
!> @param bn [inout] Argument bn.
pure subroutine bernob ( n, bn )

!*****************************************************************************80
!
!! BERNOB computes the Bernoulli number Bn.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    11 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) N, the index.
!
!    Output, real(real64) BN, the value of the N-th Bernoulli number.
!
  implicit none

  integer(int32), intent(in) :: n

  real(real64), intent(inout) :: bn(0:n)
  integer(int32) k
  integer(int32) m
  real(real64) r1
  real(real64) r2
  real(real64) s
  real(real64) tpi

  tpi = 6.283185307179586e+00_real64
  bn(0) = 1.0e+00_real64
  bn(1) = -0.5e+00_real64
  bn(2) = 1.0e+00_real64 / 6.0e+00_real64
  r1 = ( 2.0e+00_real64 / tpi )**2

  do m = 4, n, 2

    r1 = - r1 * ( m - 1 ) * m / ( tpi * tpi )
    r2 = 1.0e+00_real64
 
    do k = 2, 10000
      s = ( 1.0e+00_real64 / k ) ** m
      r2 = r2 + s
      if ( s < 1.0e-15_real64 ) then
        exit
      end if
    end do

    bn(m) = r1 * r2

  end do

  return
end subroutine bernob
!> @brief subroutine beta.
!> @return None.
!>
!> @param p [in] Argument p.
!> @param q [in] Argument q.
!> @param bt [inout] Argument bt.
subroutine beta ( p, q, bt )

!*****************************************************************************80
!
!! BETA computes the Beta function B(p,q).
!
!  Licensing:
!
!    The original FORTRAN77 version of this routine is copyrighted by 
!    Shanjie Zhang and Jianming Jin.  However, they give permission to 
!    incorporate this routine into a user program that the copyright 
!    is acknowledged.
!
!  Modified:
!
!    12 March 2012
!
!  Author:
!
!    Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45
!
!  Parameters:
!
!    Input, real(real64) P, Q, the parameters.
!    0 < P, 0 < Q.
!
!    Output, real(real64) BT, the value of B(P,Q).
!
  implicit none

  real(real64), intent(inout) :: bt
  real(real64) gp
  real(real64) gpq
  real(real64) gq
  real(real64), intent(in) :: p
  real(real64) ppq
  real(real64), intent(in) :: q

  call gamma ( p, gp )
  call gamma ( q, gq )
  ppq = p + q
  call gamma ( ppq, gpq )
  bt = gp * gq / gpq

  return
end subroutine beta
!> @brief subroutine bjndd.
!> @return None.
!>
!> @param n [in] Argument n.
!> @param x [in] Argument x.
!> @param bj [inout] Argument bj.
!> @param dj [inout] Argument dj.
!> @param fj [inout] Argument fj.
pure subroutine bjndd ( n, x, bj, dj, fj )

!*****************************************************************************80
!
!! BJNDD computes Bessel functions Jn(x) and first and second derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    11 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) N, the order.
!
!    Input, real(real64) X, the argument.
!
!    Output, real(real64) BJ(N+1), DJ(N+1), FJ(N+1), the values of 
!    Jn(x), Jn'(x) and Jn''(x) in the last entries.
!
  implicit none

  integer(int32), intent(in) :: n

  real(real64), intent(inout) :: bj(n+1)
  real(real64) bs
  real(real64), intent(inout) :: dj(n+1)
  real(real64) f
  real(real64) f0
  real(real64) f1
  real(real64), intent(inout) :: fj(n+1)
  integer(int32) k
  integer(int32) m
  integer(int32) mt
  integer(int32) nt
  real(real64), intent(in) :: x

  do nt = 1, 900
    mt = int ( 0.5e+00_real64 * log10 ( 6.28e+00_real64 * nt ) &
      - nt * log10 ( 1.36e+00_real64 * abs ( x ) / nt ) )
    if ( 20 < mt ) then
      exit
    end if
  end do

  m = nt
  bs = 0.0e+00_real64
  f0 = 0.0e+00_real64
  f1 = 1.0e-35_real64
  do k = m, 0, -1
    f = 2.0e+00_real64 * ( k + 1.0e+00_real64 ) * f1 / x - f0
    if ( k <= n ) then
      bj(k+1) = f
    end if
    if ( k == 2 * int ( k / 2 ) ) then
      bs = bs + 2.0e+00_real64 * f
    end if
    f0 = f1
    f1 = f
  end do

  do k = 0, n
    bj(k+1) = bj(k+1) / ( bs - f )
  end do

  dj(1) = -bj(2)
  fj(1) = -1.0e+00_real64 * bj(1) - dj(1) / x
  do k = 1, n
    dj(k+1) = bj(k) - k * bj(k+1) / x
    fj(k+1) = ( k * k / ( x * x ) - 1.0e+00_real64 ) * bj(k+1) - dj(k+1) / x
  end do

  return
end subroutine bjndd
!> @brief subroutine cbk.
!> @return None.
!>
!> @param m [in] Argument m.
!> @param n [in] Argument n.
!> @param c [in] Argument c.
!> @param cv [in] Argument cv.
!> @param qt [in] Argument qt.
!> @param ck [in] Argument ck.
!> @param bk [inout] Argument bk.
pure subroutine cbk ( m, n, c, cv, qt, ck, bk )

!*****************************************************************************80
!
!! CBK computes coefficients for oblate radial functions with small argument.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    20 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) M, the mode parameter;  M = 0, 1, 2, ...
!
!    Input, integer(int32) N, mode parameter, N = M, M + 1, M + 2, ...
!
!    Input, real(real64) C, spheroidal parameter.
!
!    Input, real(real64) CV, the characteristic value.
!
!    Input, real(real64) QT, ?
!
!    Input, real(real64) CK(*), ?
!
!    Output, real(real64) BK(*), the coefficients.
!
  implicit none

  real(real64), intent(inout) :: bk(200)
  real(real64), intent(in) :: c
  real(real64), intent(in) :: ck(200)
  real(real64), intent(in) :: cv
  real(real64) eps
  integer(int32) i
  integer(int32) i1
  integer(int32) ip
  integer(int32) j
  integer(int32) k
  integer(int32), intent(in) :: m
  integer(int32), intent(in) :: n
  integer(int32) n2
  integer(int32) nm
  real(real64), intent(in) :: qt
  real(real64) r1
  real(real64) s1
  real(real64) sw
  real(real64) t
  real(real64) u(200)
  real(real64) v(200)
  real(real64) w(200)

  eps = 1.0e-14_real64
  if ( n - m == 2 * int ( ( n - m ) / 2 ) ) then
    ip = 0
  else
    ip = 1
  end if
  nm = 25 + int ( 0.5e+00_real64 * ( n - m ) + c )
  u(1) = 0.0e+00_real64
  n2 = nm - 2
  do j = 2, n2
    u(j) = c * c
  end do

  do j = 1, n2
    v(j) = ( 2.0e+00_real64 * j - 1.0e+00_real64 - ip ) &
      * ( 2.0e+00_real64 * ( j - m ) - ip ) + m * ( m - 1.0e+00_real64 ) - cv
  end do

  do j = 1, nm - 1
    w(j) = ( 2.0e+00_real64 * j - ip ) * ( 2.0e+00_real64 * j + 1.0e+00_real64 - ip )
  end do

  if ( ip == 0 ) then

    do k = 0, n2 - 1

      s1 = 0.0e+00_real64
      sw = 0.0e+00_real64
      i1 = k - m + 1

      do i = i1, nm
        if ( 0 <= i ) then
          r1 = 1.0e+00_real64
          do j = 1, k
            r1 = r1 * ( i + m - j ) / j
          end do
          s1 = s1 + ck(i+1) * ( 2.0e+00_real64 * i + m ) * r1
          if ( abs ( s1 - sw ) < abs ( s1 ) * eps ) then
            exit
          end if
          sw = s1
        end if
      end do

      bk(k+1) = qt * s1

    end do

  else if ( ip == 1 ) then

    do k = 0, n2 - 1

      s1 = 0.0e+00_real64
      sw = 0.0e+00_real64
      i1 = k - m + 1

      do i = i1, nm

        if ( 0 <= i ) then

          r1 = 1.0e+00_real64
          do j = 1, k
            r1 = r1 * ( i + m - j ) / j
          end do

          if ( 0 < i ) then
            s1 = s1 + ck(i) * ( 2.0e+00_real64 * i + m - 1 ) * r1
          end if
          s1 = s1 - ck(i+1) * ( 2.0e+00_real64 * i + m ) * r1
          if ( abs ( s1 - sw ) < abs ( s1 ) * eps ) then
            exit
          end if
          sw = s1

        end if

      end do

      bk(k+1) = qt * s1

    end do

  end if

  w(1) = w(1) / v(1)
  bk(1) = bk(1) / v(1)
  do k = 2, n2
    t = v(k) - w(k-1) * u(k)
    w(k) = w(k) / t
    bk(k) = ( bk(k) - bk(k-1) * u(k) ) / t
  end do

  do k = n2 - 1, 1, -1
    bk(k) = bk(k) - w(k) * bk(k+1)
  end do

  return
end subroutine cbk
!> @brief subroutine cchg.
!> @return None.
!>
!> @param a [inout] Argument a.
!> @param b [in] Argument b.
!> @param z [inout] Argument z.
!> @param chg [inout] Argument chg.
subroutine cchg ( a, b, z, chg )

!*****************************************************************************80
!
!! CCHG computes the confluent hypergeometric function.
!
!  Discussion:
!
!    This function computes the confluent hypergeometric function
!    M(a,b,z) with real parameters a, b and complex argument z.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    26 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) A, B, parameter values.
!
!    Input, complex ( kind = real64 ) Z, the argument.
!
!    Output, complex ( kind = real64 ) CHG, the value of M(a,b,z).
!
  implicit none

  real(real64), intent(inout) :: a
  real(real64) a0
  real(real64) a1
  real(real64), intent(in) :: b
  real(real64) ba
  complex(real64) cfac
  complex(real64), intent(inout) :: chg
  complex(real64) chg1
  complex(real64) chg2
  complex(real64) chw
  complex(real64) ci
  complex(real64) cr
  complex(real64) cr1
  complex(real64) cr2
  complex(real64) crg
  complex(real64) cs1
  complex(real64) cs2
  complex(real64) cy0
  complex(real64) cy1
  real(real64) g1
  real(real64) g2
  real(real64) g3
  integer(int32) i
  integer(int32) j
  integer(int32) k
  integer(int32) la
  integer(int32) m
  integer(int32) n
  integer(int32) nl
  integer(int32) ns
  real(real64) phi
  real(real64) pi
  real(real64) x
  real(real64) x0
  real(real64) y
  complex(real64), intent(inout) :: z
  complex(real64) z0

  pi = 3.141592653589793e+00_real64
  ci = cmplx ( 0.0e+00_real64, 1.0e+00_real64, kind = real64 )
  a0 = a
  a1 = a
  z0 = z

  if ( b == 0.0e+00_real64 .or. b == - int ( abs ( b ) ) ) then
    chg = cmplx ( 1.0e+30_real64, 0.0e+00_real64, kind = real64 )
  else if ( a == 0.0e+00_real64 .or. z == 0.0e+00_real64 ) then
    chg = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
  else if ( a == -1.0e+00_real64 ) then
    chg = 1.0e+00_real64 - z / b
  else if ( a == b ) then
    chg = exp ( z )
  else if ( a - b == 1.0e+00_real64 ) then
    chg = ( 1.0e+00_real64 + z / b ) * exp ( z )
  else if ( a == 1.0e+00_real64 .and. b == 2.0e+00_real64 ) then
    chg = ( exp ( z ) - 1.0e+00_real64 ) / z
  else if ( a == int ( a ) .and. a < 0.0e+00_real64 ) then
    m = int ( - a )
    cr = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    chg = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    do k = 1, m
      cr = cr * ( a + k - 1.0e+00_real64 ) / k / ( b + k - 1.0e+00_real64 ) * z
      chg = chg + cr
    end do
  else

    x0 = real ( z, kind = real64 )
    if ( x0 < 0.0e+00_real64 ) then
      a = b - a
      a0 = a
      z = - z
    end if

    if ( a < 2.0e+00_real64 ) then
      nl = 0
    else
      nl = 1
      la = int ( a )
      a = a - la - 1.0e+00_real64
    end if

    do n = 0, nl

      if ( 2.0e+00_real64 <= a0 ) then
        a = a + 1.0e+00_real64
      end if

      if ( cdabs ( z ) < 20.0e+00_real64 + abs ( b ) .or. a < 0.0e+00_real64 ) then

        chg = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
        crg = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
        chw = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
        do j = 1, 500
          crg = crg * ( a + j - 1.0e+00_real64 ) / ( j * ( b + j - 1.0e+00_real64 ) ) * z
          chg = chg + crg
          if ( abs ( ( chg - chw ) / chg ) < 1.0e-15_real64 ) then
            exit
          end if
          chw = chg
        end do

      else

        call gamma ( a, g1 )
        call gamma ( b, g2 )
        ba = b - a
        call gamma ( ba, g3 )
        cs1 = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
        cs2 = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
        cr1 = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
        cr2 = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )

        do i = 1, 8
          cr1 = - cr1 * (     a + i - 1.0e+00_real64 ) * ( a - b + i ) / ( z * i )
          cr2 =   cr2 * ( b - a + i - 1.0e+00_real64 ) * ( i - a ) / ( z * i )
          cs1 = cs1 + cr1
          cs2 = cs2 + cr2
        end do

        x = real ( z, kind = real64 )
        y = imag ( z )

        if ( x == 0.0e+00_real64 .and. 0.0e+00_real64 <= y ) then
          phi = 0.5e+00_real64 * pi
        else if ( x == 0.0e+00_real64 .and. y <= 0.0e+00_real64 ) then
          phi = -0.5e+00_real64 * pi
        else
          phi = atan ( y / x )
        end if

        if ( -1.5e+00_real64 * pi < phi .and. phi <= -0.5_real64 * pi ) then
          ns = -1
        else if ( -0.5e+00_real64 * pi < phi .and. phi < 1.5e+00_real64 * pi ) then
          ns = 1
        end if

        if ( y == 0.0e+00_real64 ) then
          cfac = cos ( pi * a )
        else
          cfac = exp ( ns * ci * pi * a )
        end if

        chg1 = g2 / g3 * z ** ( - a ) * cfac * cs1
        chg2 = g2 / g1 * exp ( z ) * z ** ( a - b ) * cs2
        chg = chg1 + chg2

      end if

      if ( n == 0 ) then
        cy0 = chg
      else if ( n == 1 ) then
        cy1 = chg
      end if

    end do

    if ( 2.0e+00_real64 <= a0 ) then
      do i = 1, la - 1
        chg = ( ( 2.0e+00_real64 * a - b + z ) * cy1 + ( b - a ) * cy0 ) / a
        cy0 = cy1
        cy1 = chg
        a = a + 1.0e+00_real64
      end do
    end if

    if ( x0 < 0.0e+00_real64 ) then
      chg = chg * exp ( - z )
    end if

  end if

  a = a1
  z = z0

  return
end subroutine cchg
!> @brief subroutine cerf.
!> @return None.
!>
!> @param z [in] Argument z.
!> @param cer [inout] Argument cer.
!> @param cder [inout] Argument cder.
pure subroutine cerf ( z, cer, cder )

!*****************************************************************************80
!
!! CERF computes the error function and derivative for a complex argument.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    25 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
! 
!  Parameters:
!
!    Input, complex ( kind = real64 ), the argument.
!
!    Output, complex ( kind = real64 ) CER, CDER, the values of erf(z) and erf'(z).
!
  implicit none

  real(real64) c0
  complex(real64), intent(inout) :: cder
  complex(real64), intent(inout) :: cer
  real(real64) cs
  real(real64) ei1
  real(real64) ei2
  real(real64) eps
  real(real64) er
  real(real64) er0
  real(real64) er1
  real(real64) er2 
  real(real64) eri
  real(real64) err
  integer(int32) k
  integer(int32) n
  real(real64) pi
  real(real64) r
  real(real64) ss
  real(real64) w
  real(real64) w1
  real(real64) w2
  real(real64) x
  real(real64) x2
  real(real64) y
  complex(real64), intent(in) :: z

  eps = 1.0e-12_real64
  pi = 3.141592653589793e+00_real64
  x = real ( z, kind = real64 )
  y = imag ( z )
  x2 = x * x

  if ( x <= 3.5e+00_real64 ) then

    er = 1.0e+00_real64
    r = 1.0e+00_real64
    w = 0.0e+00_real64
    do k = 1, 100
      r = r * x2 / ( k + 0.5e+00_real64 )
      er = er + r
      if ( abs ( er - w ) <= eps * abs ( er ) ) then
        exit
      end if
      w = er
    end do

    c0 = 2.0e+00_real64 / sqrt ( pi ) * x * exp ( - x2 )
    er0 = c0 * er

  else

    er = 1.0e+00_real64
    r = 1.0e+00_real64
    do k = 1, 12
      r = - r * ( k - 0.5e+00_real64 ) / x2
      er = er + r
    end do
    c0 = exp ( - x2 ) / ( x * sqrt ( pi ) )
    er0 = 1.0e+00_real64 - c0 * er

  end if

  if ( y == 0.0e+00_real64 ) then

    err = er0
    eri = 0.0e+00_real64

  else

    cs = cos ( 2.0e+00_real64 * x * y )
    ss = sin ( 2.0e+00_real64 * x * y )
    er1 = exp ( - x2 ) * ( 1.0e+00_real64 - cs ) / ( 2.0e+00_real64 * pi * x )
    ei1 = exp ( - x2 ) * ss / ( 2.0e+00_real64 * pi * x )
    er2 = 0.0e+00_real64
    w1 = 0.0e+00_real64
    do n = 1, 100
      er2 = er2 + exp ( - 0.25e+00_real64 * n * n ) &
        / ( n * n + 4.0e+00_real64 * x2 ) * ( 2.0e+00_real64 * x &
        - 2.0e+00_real64 * x * cosh ( n * y ) * cs &
        + n * sinh ( n * y ) * ss )
      if ( abs ( ( er2 - w1 ) / er2 ) < eps ) then
        exit
      end if
      w1 = er2
    end do

    c0 = 2.0e+00_real64 * exp ( - x2 ) / pi
    err = er0 + er1 + c0 * er2
    ei2 = 0.0e+00_real64
    w2 = 0.0e+00_real64
    do n = 1, 100
      ei2 = ei2 + exp ( - 0.25e+00_real64 * n * n ) &
        / ( n * n + 4.0e+00_real64 * x2 ) * ( 2.0e+00_real64 * x &
        * cosh ( n * y ) * ss + n * sinh ( n * y ) * cs )
      if ( abs ( ( ei2 - w2 ) / ei2 ) < eps ) then
        exit
      end if
      w2 = ei2
    end do

    eri = ei1 + c0 * ei2

  end if

  cer = cmplx ( err, eri, kind = real64 )
  cder = 2.0e+00_real64 / sqrt ( pi ) * exp ( - z * z )

  return
end subroutine cerf
!> @brief subroutine cerror.
!> @return None.
!>
!> @param z [in] Argument z.
!> @param cer [inout] Argument cer.
pure subroutine cerror ( z, cer )

!*****************************************************************************80
!
!! CERROR computes the error function for a complex argument.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    15 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
! 
!  Parameters:
!
!    Input, complex ( kind = real64 ) Z, the argument.
!
!    Output, complex ( kind = real64 ) CER, the function value.
!
  implicit none

  real(real64) a0
  complex(real64) c0
  complex(real64), intent(inout) :: cer
  complex(real64) cl
  complex(real64) cr
  complex(real64) cs
  integer(int32) k
  real(real64) pi
  complex(real64), intent(in) :: z
  complex(real64) z1

  a0 = abs ( z )
  c0 = exp ( - z * z )
  pi = 3.141592653589793e+00_real64
  z1 = z

  if ( real ( z, kind = real64 ) < 0.0e+00_real64 ) then
    z1 = - z
  end if

  if ( a0 <= 5.8e+00_real64 ) then    

    cs = z1
    cr = z1
    do k = 1, 120
      cr = cr * z1 * z1 / ( k + 0.5e+00_real64 )
      cs = cs + cr
      if ( abs ( cr / cs ) < 1.0e-15_real64 ) then
        exit
      end if
    end do

    cer = 2.0e+00_real64 * c0 * cs / sqrt ( pi )

  else

    cl = 1.0e+00_real64 / z1              
    cr = cl
    do k = 1, 13
      cr = -cr * ( k - 0.5e+00_real64 ) / ( z1 * z1 )
      cl = cl + cr
      if ( abs ( cr / cl ) < 1.0e-15_real64 ) then
        exit
      end if
    end do

    cer = 1.0e+00_real64 - c0 * cl / sqrt ( pi )

  end if

  if ( real ( z, kind = real64 ) < 0.0e+00_real64 ) then
    cer = -cer
  end if

  return
end subroutine cerror
!> @brief subroutine cerzo.
!> @return None.
!>
!> @param nt [in] Argument nt.
!> @param zo [inout] Argument zo.
subroutine cerzo ( nt, zo )

!*****************************************************************************80
!
!! CERZO evaluates the complex zeros of the error function.
!
!  Discussion:
!
!    The modified Newton method is used.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    15 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) NT, the number of zeros.
!
!    Output, complex ( kind = real64 ) ZO(NT), the zeros.
!
  implicit none

  integer(int32), intent(in) :: nt

  integer(int32) i
  integer(int32) it
  integer(int32) j
  integer(int32) nr
  real(real64) pi
  real(real64) pu
  real(real64) pv
  real(real64) px
  real(real64) py
  real(real64) w
  real(real64) w0
  complex(real64) z
  complex(real64) zd
  complex(real64) zf
  complex(real64) zfd
  complex(real64) zgd
  complex(real64), intent(inout) :: zo(nt)
  complex(real64) zp
  complex(real64) zq
  complex(real64) zw

  pi = 3.141592653589793e+00_real64

  do nr = 1, nt

    pu = sqrt ( pi * ( 4.0e+00_real64 * nr - 0.5e+00_real64 ) )
    pv = pi * sqrt ( 2.0e+00_real64 * nr - 0.25e+00_real64 )
    px = 0.5e+00_real64 * pu - 0.5e+00_real64 * log ( pv ) / pu
    py = 0.5e+00_real64 * pu + 0.5e+00_real64 * log ( pv ) / pu
    z = cmplx ( px, py, kind = real64 )
    it = 0

    do

      it = it + 1
      call cerf ( z, zf, zd )
      zp = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
      do i = 1, nr - 1
        zp = zp * ( z - zo(i) )
      end do
      zfd = zf / zp

      zq = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
      do i = 1, nr - 1
        zw = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
        do j = 1, nr - 1
          if ( j /= i ) then
            zw = zw * ( z - zo(j) )
          end if
        end do
        zq = zq + zw
      end do

      zgd = ( zd - zq * zfd ) / zp
      z = z - zfd / zgd
      w0 = w
      w = abs ( z )

      if ( 50 < it .or. abs ( ( w - w0 ) / w ) <= 1.0e-11_real64 ) then
        exit
      end if

    end do

    zo(nr) = z

  end do

  return
end subroutine cerzo
!> @brief subroutine cfc.
!> @return None.
!>
!> @param z [in] Argument z.
!> @param zf [inout] Argument zf.
!> @param zd [inout] Argument zd.
pure subroutine cfc ( z, zf, zd )

!*****************************************************************************80
!
!! CFC computes the complex Fresnel integral C(z) and C'(z).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    26 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, complex ( kind = real64 ) Z, the argument.
!
!    Output, complex ( kind = real64 ) ZF, ZD, the values of C(z) and C'(z).
!
  implicit none

  complex(real64) c
  complex(real64) cf
  complex(real64) cf0
  complex(real64) cf1
  complex(real64) cg
  complex(real64) cr
  real(real64) eps
  integer(int32) k
  integer(int32) m
  real(real64) pi
  real(real64) w0
  real(real64) wa
  real(real64) wa0
  complex(real64), intent(in) :: z
  complex(real64) z0
  complex(real64), intent(inout) :: zd
  complex(real64), intent(inout) :: zf
  complex(real64) zp
  complex(real64) zp2

  eps = 1.0e-14_real64
  pi = 3.141592653589793e+00_real64
  w0 = abs ( z )
  zp = 0.5e+00_real64 * pi * z * z
  zp2 = zp * zp
  z0 = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )

  if ( z == z0 ) then

    c = z0

  else if ( w0 <= 2.5e+00_real64 ) then

    cr = z
    c = cr
    wa0 = 0.0_real64
    do k = 1, 80
      cr = -0.5e+00_real64 * cr * ( 4.0e+00_real64 * k - 3.0e+00_real64 ) &
        / k / ( 2.0e+00_real64 * k - 1.0e+00_real64 ) &
        / ( 4.0e+00_real64 * k + 1.0e+00_real64 ) * zp2
      c = c + cr
      wa = abs ( c )
      if ( abs ( ( wa - wa0 ) / wa ) < eps .and. 10 < k ) then
        exit
      end if
      wa0 = wa
    end do

  else if ( 2.5e+00_real64 < w0 .and. w0 < 4.5e+00_real64 ) then

    m = 85
    c = z0
    cf1 = z0
    cf0 = cmplx ( 1.0e-30_real64, 0.0e+00_real64, kind = real64 )
    do k = m, 0, -1
      cf = ( 2.0e+00_real64 * k + 3.0e+00_real64 ) * cf0 / zp - cf1
      if ( k == int ( k / 2 ) * 2 ) then
        c = c + cf
      end if
      cf1 = cf0
      cf0 = cf
    end do
    c = sqrt ( 2.0e+00_real64 / ( pi * zp ) ) * sin ( zp ) / cf * c

  else

    cr = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    cf = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    do k = 1, 20
      cr = - 0.25e+00_real64 * cr * ( 4.0e+00_real64 * k - 1.0e+00_real64 ) &
        * ( 4.0e+00_real64 * k - 3.0e+00_real64 ) / zp2
      cf = cf + cr
    end do
    cr = 1.0e+00_real64 / ( pi * z * z )
    cg = cr
    do k = 1, 12
      cr = - 0.25e+00_real64 * cr * ( 4.0e+00_real64 * k + 1.0e+00_real64 ) &
        * ( 4.0e+00_real64 * k - 1.0e+00_real64 ) / zp2
      cg = cg + cr
    end do
    c = 0.5e+00_real64 + ( cf * sin ( zp ) - cg * cos ( zp ) ) / ( pi * z )

  end if

  zf = c
  zd = cos ( 0.5e+00_real64 * pi * z * z )

  return
end subroutine cfc
!> @brief subroutine cfs.
!> @return None.
!>
!> @param z [in] Argument z.
!> @param zf [inout] Argument zf.
!> @param zd [inout] Argument zd.
pure subroutine cfs ( z, zf, zd )

!*****************************************************************************80
!
!! CFS computes the complex Fresnel integral S(z) and S'(z).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    24 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, complex ( kind = real64 ) Z, the argument.
!
!    Output, complex ( kind = real64 ) ZF, ZD, the values of S(z) and S'(z).
!
  implicit none

  complex(real64) cf
  complex(real64) cf0
  complex(real64) cf1
  complex(real64) cg
  complex(real64) cr
  real(real64) eps
  integer(int32) k
  integer(int32) m
  real(real64) pi
  complex(real64) s
  real(real64) w0
  real(real64) wb
  real(real64) wb0
  complex(real64), intent(in) :: z
  complex(real64) z0
  complex(real64), intent(inout) :: zd
  complex(real64), intent(inout) :: zf
  complex(real64) zp
  complex(real64) zp2

  eps = 1.0e-14_real64
  pi = 3.141592653589793e+00_real64
  w0 = abs ( z )
  zp = 0.5e+00_real64 * pi * z * z
  zp2 = zp * zp
  z0 = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )

  if ( z == z0 ) then

    s = z0

  else if ( w0 <= 2.5e+00_real64 ) then

    s = z * zp / 3.0e+00_real64
    cr = s
    wb0 = 0.0_real64
    do k = 1, 80
      cr = -0.5e+00_real64 * cr * ( 4.0e+00_real64 * k - 1.0e+00_real64 ) / k &
        / ( 2.0e+00_real64 * k + 1.0e+00_real64 ) &
        / ( 4.0e+00_real64 * k + 3.0e+00_real64 ) * zp2
      s = s + cr
      wb = abs ( s )
      if ( abs ( wb - wb0 ) < eps .and. 10 < k ) then
        exit
      end if
      wb0 = wb
    end do

  else if ( 2.5e+00_real64 < w0 .and. w0 < 4.5e+00_real64 ) then

    m = 85
    s = z0
    cf1 = z0
    cf0 = cmplx ( 1.0e-30_real64, 0.0e+00_real64, kind = real64  )
    do k = m, 0, -1
      cf = ( 2.0e+00_real64 * k + 3.0e+00_real64 ) * cf0 / zp - cf1
      if ( k /= int ( k / 2 ) * 2 ) then
        s = s + cf
      end if
      cf1 = cf0
      cf0 = cf
    end do
    s = sqrt ( 2.0e+00_real64 / ( pi * zp ) ) * sin ( zp ) / cf * s

  else

    cr = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64  )
    cf = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64  )
    do k = 1, 20
      cr = -0.25e+00_real64 * cr * ( 4.0e+00_real64 * k - 1.0e+00_real64 ) &
        * ( 4.0e+00_real64 * k - 3.0e+00_real64 ) / zp2
      cf = cf + cr
    end do
    cr = 1.0e+00_real64 / ( pi * z * z )
    cg = cr
    do k = 1, 12
      cr = -0.25e+00_real64 * cr * ( 4.0e+00_real64 * k + 1.0e+00_real64 ) &
        * ( 4.0e+00_real64 * k - 1.0e+00_real64 ) / zp2
      cg = cg + cr
    end do
    s = 0.5e+00_real64 - ( cf * cos ( zp ) + cg * sin ( zp ) ) / ( pi * z )

  end if

  zf = s
  zd = sin ( 0.5e+00_real64 * pi * z * z )

  return
end subroutine cfs
!> @brief subroutine cgama.
!> @return None.
!>
!> @param x [inout] Argument x.
!> @param y [inout] Argument y.
!> @param kf [in] Argument kf.
!> @param gr [inout] Argument gr.
!> @param gi [inout] Argument gi.
pure subroutine cgama ( x, y, kf, gr, gi )

!*****************************************************************************80
!
!! CGAMA computes the Gamma function for complex argument.
!
!  Discussion:
!
!    This procedcure computes the gamma function �(z) or ln[�(z)]
!    for a complex argument
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    26 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) X, Y, the real and imaginary parts of 
!    the argument Z.
!
!    Input, integer(int32) KF, the function code.
!    0 for ln[�(z)]
!    1 for �(z)
!
!    Output, real(real64) GR, GI, the real and imaginary parts of
!    the selected function.
!
  implicit none

  real(real64), parameter, dimension ( 10 ) :: a = [&
    8.333333333333333e-02_real64, -2.777777777777778e-03_real64, &
    7.936507936507937e-04_real64, -5.952380952380952e-04_real64, &
    8.417508417508418e-04_real64, -1.917526917526918e-03_real64, &
    6.410256410256410e-03_real64, -2.955065359477124e-02_real64, &
    1.796443723688307e-01_real64, -1.39243221690590e+00_real64]
  real(real64) g0
  real(real64), intent(inout) :: gi
  real(real64) gi1
  real(real64), intent(inout) :: gr
  real(real64) gr1
  integer(int32) j
  integer(int32) k
  integer(int32), intent(in) :: kf
  integer(int32) na
  real(real64) pi
  real(real64) si
  real(real64) sr
  real(real64) t
  real(real64) th
  real(real64) th1
  real(real64) th2
  real(real64), intent(inout) :: x
  real(real64) x0
  real(real64) x1
  real(real64), intent(inout) :: y
  real(real64) y1
  real(real64) z1
  real(real64) z2

  pi = 3.141592653589793e+00_real64

  if ( y == 0.0e+00_real64 .and. x == int ( x ) .and. x <= 0.0e+00_real64 ) then
    gr = 1.0e+300_real64
    gi = 0.0e+00_real64
    return
  else if ( x < 0.0e+00_real64 ) then
    x1 = x
    y1 = y
    x = -x
    y = -y
  end if

  x0 = x

  if ( x <= 7.0e+00_real64 ) then
    na = int ( 7 - x )
    x0 = x + na
  end if

  z1 = sqrt ( x0 * x0 + y * y )
  th = atan ( y / x0 )
  gr = ( x0 - 0.5e+00_real64 ) * log ( z1 ) - th * y - x0 &
    + 0.5e+00_real64 * log ( 2.0e+00_real64 * pi )
  gi = th * ( x0 - 0.5e+00_real64 ) + y * log ( z1 ) - y

  do k = 1, 10
    t = z1 ** ( 1 - 2 * k )
    gr = gr + a(k) * t * cos ( ( 2.0e+00_real64 * k - 1.0e+00_real64 ) * th )
    gi = gi - a(k) * t * sin ( ( 2.0e+00_real64 * k - 1.0e+00_real64 ) * th )
  end do

  if ( x <= 7.0e+00_real64 ) then
    gr1 = 0.0e+00_real64
    gi1 = 0.0e+00_real64
    do j = 0, na - 1
      gr1 = gr1 + 0.5e+00_real64 * log ( ( x + j ) ** 2 + y * y )
      gi1 = gi1 + atan ( y / ( x + j ) )
    end do
    gr = gr - gr1
    gi = gi - gi1
  end if

  if ( x1 < 0.0e+00_real64 ) then
    z1 = sqrt ( x * x + y * y )
    th1 = atan ( y / x )
    sr = - sin ( pi * x ) * cosh ( pi * y )
    si = - cos ( pi * x ) * sinh ( pi * y )
    z2 = sqrt ( sr * sr + si * si )
    th2 = atan ( si / sr )
    if ( sr < 0.0e+00_real64 ) then
      th2 = pi + th2
    end if
    gr = log ( pi / ( z1 * z2 ) ) - gr
    gi = - th1 - th2 - gi
    x = x1
    y = y1
  end if

  if ( kf == 1 ) then
    g0 = exp ( gr )
    gr = g0 * cos ( gi )
    gi = g0 * sin ( gi )
  end if

  return
end subroutine cgama
!> @brief subroutine ch12n.
!> @return None.
!>
!> @param n [in] Argument n.
!> @param z [in] Argument z.
!> @param nm [inout] Argument nm.
!> @param chf1 [inout] Argument chf1.
!> @param chd1 [inout] Argument chd1.
!> @param chf2 [inout] Argument chf2.
!> @param chd2 [inout] Argument chd2.
subroutine ch12n ( n, z, nm, chf1, chd1, chf2, chd2 )

!*****************************************************************************80
!
!! CH12N computes Hankel functions of first and second kinds, complex argument.
!
!  Discussion:
!
!    Both the Hankel functions and their derivatives are computed.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    26 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
! 
!  Parameters:
!
!    Input, integer(int32) N, the order of the functions.
!
!    Input, complex ( kind = real64 ) Z, the argument.
!
!    Output, integer(int32) NM, the highest order computed.
!
!    Output, complex ( kind = real64 ) CHF1(0:n), CHD1(0:n), CHF2(0:n), CHD2(0:n), 
!    the values of Hn(1)(z), Hn(1)'(z), Hn(2)(z), Hn(2)'(z).
!
  implicit none

  integer(int32), intent(in) :: n

  complex(real64) cbi(0:250)
  complex(real64) cbj(0:250)
  complex(real64) cbk(0:250)
  complex(real64) cby(0:250)
  complex(real64) cdi(0:250)
  complex(real64) cdj(0:250)
  complex(real64) cdk(0:250)
  complex(real64) cdy(0:250)
  complex(real64), intent(inout) :: chd1(0:n)
  complex(real64), intent(inout) :: chd2(0:n)
  complex(real64) cf1
  complex(real64) cfac
  complex(real64), intent(inout) :: chf1(0:n)
  complex(real64), intent(inout) :: chf2(0:n)
  complex(real64) ci
  integer(int32) k
  integer(int32), intent(inout) :: nm
  real(real64) pi
  complex(real64), intent(in) :: z
  complex(real64) zi

  ci = cmplx ( 0.0e+00_real64, 1.0e+00_real64, kind = real64 )
  pi = 3.141592653589793e+00_real64

  if ( imag ( z ) < 0.0e+00_real64 ) then

    call cjynb ( n, z, nm, cbj, cdj, cby, cdy )

    do k = 0, nm
      chf1(k) = cbj(k) + ci * cby(k)
      chd1(k) = cdj(k) + ci * cdy(k)
    end do

    zi = ci * z
    call ciknb ( n, zi, nm, cbi, cdi, cbk, cdk )
    cfac = -2.0e+00_real64 / ( pi * ci )

    do k = 0, nm
      chf2(k) = cfac * cbk(k)
      chd2(k) = cfac * ci * cdk(k)
      cfac = cfac * ci
    end do

  else if ( 0.0e+00_real64 < imag ( z ) ) then

    zi = - ci * z
    call ciknb ( n, zi, nm, cbi, cdi, cbk, cdk )
    cf1 = -ci
    cfac = 2.0e+00_real64 / ( pi * ci )

    do k = 0, nm
      chf1(k) = cfac * cbk(k)
      chd1(k) = -cfac * ci * cdk(k)
      cfac = cfac * cf1
    end do

    call cjynb ( n, z, nm, cbj, cdj, cby, cdy )

    do k = 0, nm
      chf2(k) = cbj(k) - ci * cby(k)
      chd2(k) = cdj(k) - ci * cdy(k)
    end do

  else

    call cjynb ( n, z, nm, cbj, cdj, cby, cdy )

    do k = 0, nm
      chf1(k) = cbj(k) + ci * cby(k)
      chd1(k) = cdj(k) + ci * cdy(k)
      chf2(k) = cbj(k) - ci * cby(k)
      chd2(k) = cdj(k) - ci * cdy(k)
    end do

  end if

  return
end subroutine ch12n
!> @brief subroutine chgm.
!> @return None.
!>
!> @param a [inout] Argument a.
!> @param b [in] Argument b.
!> @param x [inout] Argument x.
!> @param hg [inout] Argument hg.
subroutine chgm ( a, b, x, hg )

!*****************************************************************************80
!
!! CHGM computes the confluent hypergeometric function M(a,b,x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    27 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) A, B, parameters.
!
!    Input, real(real64) X, the argument.
!
!    Output, real(real64) HG, the value of M(a,b,x).
!
  implicit none

  real(real64), intent(inout) :: a
  real(real64) a0
  real(real64) a1
  real(real64), intent(in) :: b
  real(real64), intent(inout) :: hg
  real(real64) hg1
  real(real64) hg2
  integer(int32) i
  integer(int32) j
  integer(int32) k
  integer(int32) la
  integer(int32) m
  integer(int32) n
  integer(int32) nl
  real(real64) pi
  real(real64) r
  real(real64) r1
  real(real64) r2
  real(real64) rg
  real(real64) sum1
  real(real64) sum2
  real(real64) ta
  real(real64) tb
  real(real64) tba
  real(real64), intent(inout) :: x
  real(real64) x0
  real(real64) xg
  real(real64) y0
  real(real64) y1

  pi = 3.141592653589793e+00_real64
  a0 = a
  a1 = a
  x0 = x
  hg = 0.0e+00_real64

  if ( b == 0.0e+00_real64 .or. b == - abs ( int ( b ) ) ) then
    hg = 1.0e+300_real64
  else if ( a == 0.0e+00_real64 .or. x == 0.0e+00_real64 ) then
    hg = 1.0e+00_real64
  else if ( a == -1.0e+00_real64 ) then
    hg = 1.0e+00_real64 - x / b
  else if ( a == b ) then
    hg = exp ( x )
  else if ( a - b == 1.0e+00_real64 ) then
    hg = ( 1.0e+00_real64 + x / b ) * exp ( x )
  else if ( a == 1.0e+00_real64 .and. b == 2.0e+00_real64 ) then
    hg = ( exp ( x ) - 1.0e+00_real64 ) / x
  else if ( a == int ( a ) .and. a < 0.0e+00_real64 ) then
    m = int ( - a )
    r = 1.0e+00_real64
    hg = 1.0e+00_real64
    do k = 1, m
      r = r * ( a + k - 1.0e+00_real64 ) / k / ( b + k - 1.0e+00_real64 ) * x
      hg = hg + r
    end do
  end if

  if ( hg /= 0.0e+00_real64 ) then
    return
  end if

  if ( x < 0.0e+00_real64 ) then
    a = b - a
    a0 = a
    x = abs ( x )
  end if

  if ( a < 2.0e+00_real64 ) then
    nl = 0
  end if

  if ( 2.0e+00_real64 <= a ) then
    nl = 1
    la = int ( a )
    a = a - la - 1.0e+00_real64
  end if

  do n = 0, nl

    if ( 2.0e+00_real64 <= a0 ) then
      a = a + 1.0e+00_real64
    end if

    if ( x <= 30.0e+00_real64 + abs ( b ) .or. a < 0.0e+00_real64 ) then

      hg = 1.0e+00_real64
      rg = 1.0e+00_real64
      do j = 1, 500
        rg = rg * ( a + j - 1.0e+00_real64 ) &
          / ( j * ( b + j - 1.0e+00_real64 ) ) * x
        hg = hg + rg
        if ( abs ( rg / hg ) < 1.0e-15_real64 ) then
          exit
        end if
      end do

    else

      call gamma ( a, ta )
      call gamma ( b, tb )
      xg = b - a
      call gamma ( xg, tba )
      sum1 = 1.0e+00_real64
      sum2 = 1.0e+00_real64
      r1 = 1.0e+00_real64
      r2 = 1.0e+00_real64
      do i = 1, 8
        r1 = - r1 * ( a + i - 1.0e+00_real64 ) * ( a - b + i ) / ( x * i )
        r2 = - r2 * ( b - a + i - 1.0e+00_real64 ) * ( a - i ) / ( x * i )
        sum1 = sum1 + r1
        sum2 = sum2 + r2
      end do
      hg1 = tb / tba * x ** ( - a ) * cos ( pi * a ) * sum1
      hg2 = tb / ta * exp ( x ) * x ** ( a - b ) * sum2
      hg = hg1 + hg2

    end if

    if ( n == 0 ) then
      y0 = hg
    else if ( n == 1 ) then
      y1 = hg
    end if

  end do

  if ( 2.0e+00_real64 <= a0 ) then
    do i = 1, la - 1
      hg = ( ( 2.0e+00_real64 * a - b + x ) * y1 + ( b - a ) * y0 ) / a
      y0 = y1
      y1 = hg
      a = a + 1.0e+00_real64
    end do
  end if

  if ( x0 < 0.0e+00_real64 ) then
    hg = hg * exp ( x0 )
  end if

  a = a1
  x = x0

  return
end subroutine chgm
!> @brief subroutine chgu.
!> @return None.
!>
!> @param a [inout] Argument a.
!> @param b [inout] Argument b.
!> @param x [in] Argument x.
!> @param hu [inout] Argument hu.
!> @param md [inout] Argument md.
subroutine chgu ( a, b, x, hu, md )

!*****************************************************************************80
!
!! CHGU computes the confluent hypergeometric function U(a,b,x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    27 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
! 
!  Parameters:
!
!    Input, real(real64) A, B, parameters.
!
!    Input, real(real64) X, the argument.
!
!    Output, real(real64) HU, U(a,b,x).
!
!    Output, integer(int32) MD, the method code.
!
  implicit none

  real(real64), intent(inout) :: a
  real(real64) a00
  real(real64) aa
  real(real64), intent(inout) :: b
  real(real64) b00
  logical bl1
  logical bl2
  logical bl3
  logical bn
  real(real64), intent(inout) :: hu
  real(real64) hu1
  integer(int32) id
  integer(int32) id1
  logical il1
  logical il2
  logical il3
  integer(int32), intent(inout) :: md
  real(real64), intent(in) :: x

  aa = a - b + 1.0e+00_real64
  il1 = a == int ( a ) .and. a <= 0.0e+00_real64
  il2 = aa == int ( aa ) .and. aa <= 0.0e+00_real64
  il3 = abs ( a * ( a - b + 1.0e+00_real64 ) ) / x <= 2.0e+00_real64
  bl1 = x <= 5.0e+00_real64 .or. ( x <= 10.0e+00_real64 .and. a <= 2.0e+00_real64 )
  bl2 = ( 5.0e+00_real64 < x .and. x <= 12.5e+00_real64 ) .and. &
    ( 1.0e+00_real64 <= a .and. a + 4.0e+00_real64 <= b )
  bl3 = 12.5e+00_real64 < x .and. 5.0e+00_real64 <= a .and. a + 5.0e+00_real64 <= b
  bn = b == int ( b ) .and. b .ne. 0.0e+00_real64
  id1 = -100

  if ( b .ne. int ( b ) ) then
    call chgus ( a, b, x, hu, id1 )
    md = 1
    if ( 6 <= id1 ) then
      return
    end if
    hu1 = hu
  end if

  if ( il1 .or. il2 .or. il3 ) then
    call chgul ( a, b, x, hu, id )
    md = 2
    if ( 6 <= id ) then
      return
    end if
    if ( id < id1 ) then
      md = 1
      id = id1
      hu = hu1
    end if
  end if

  if ( 0.0e+00_real64 <= a ) then
    if ( bn .and. ( bl1 .or. bl2 .or. bl3 ) ) then
      call chgubi ( a, b, x, hu, id )
      md = 3
    else
      call chguit ( a, b, x, hu, id )
      md = 4
    end if
  else
    if ( b <= a ) then
      a00 = a
      b00 = b
      a = a - b + 1.0e+00_real64
      b = 2.0e+00_real64 - b
      call chguit ( a, b, x, hu, id )
      hu = x ** ( 1.0e+00_real64 - b00 ) * hu
      a = a00
      b = b00
      md = 4
    else if ( bn .and. ( .not. il1 ) ) then
      call chgubi ( a, b, x, hu, id )
      md = 3
    end if
  end if

  if ( id < 6 ) then
    write ( error_unit, '(a)' ) ' '
    write ( error_unit, '(a)' ) 'CHGU - Warning!'
    write ( error_unit, '(a)' ) '  Accurate results were not obtained.'
  end if

  return
end subroutine chgu
!> @brief subroutine chgubi.
!> @return None.
!>
!> @param a [in] Argument a.
!> @param b [in] Argument b.
!> @param x [in] Argument x.
!> @param hu [inout] Argument hu.
!> @param id [inout] Argument id.
subroutine chgubi ( a, b, x, hu, id )

!*****************************************************************************80
!
!! CHGUBI: confluent hypergeometric function with integer argument B.
!
!  Discussion:
!
!    This procedure computes the confluent hypergeometric function
!    U(a,b,x) with integer(int32) b ( b = �1,�2,... )
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    31 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) A, B, parameters.
!
!    Input, real(real64) X, the argument.
!
!    Output, real(real64) HU, the value of U(a,b,x).
!
!    Output, integer(int32) ID, the estimated number of significant
!    digits.
!
  implicit none

  real(real64), intent(in) :: a
  real(real64) a0
  real(real64) a1
  real(real64) a2
  real(real64), intent(in) :: b
  real(real64) da1
  real(real64) da2
  real(real64) db1
  real(real64) db2
  real(real64) el
  real(real64) ga
  real(real64) ga1
  real(real64) h0
  real(real64) hm1
  real(real64) hm2
  real(real64) hm3
  real(real64) hmax
  real(real64) hmin
  real(real64), intent(inout) :: hu
  real(real64) hu1
  real(real64) hu2
  real(real64) hw
  integer(int32), intent(inout) :: id
  integer(int32) id1
  integer(int32) id2
  integer(int32) j 
  integer(int32) k
  integer(int32) m
  integer(int32) n
  real(real64) ps
  real(real64) r
  real(real64) rn
  real(real64) rn1
  real(real64) s0
  real(real64) s1
  real(real64) s2
  real(real64) sa
  real(real64) sb
  real(real64) ua
  real(real64) ub
  real(real64), intent(in) :: x

  id = - 100
  id1 = 0
  id2 = 0
  el = 0.5772156649015329e+00_real64
  n = int ( abs ( b - 1 ) )
  rn1 = 1.0e+00_real64
  rn = 1.0e+00_real64
  do j = 1, n
    rn = rn * j
    if ( j == n - 1 ) then
      rn1 = rn
    end if
  end do

  call psi ( a, ps )
  call gamma ( a, ga )

  if ( 0.0e+00_real64 < b ) then
    a0 = a
    a1 = a - n
    a2 = a1
    call gamma ( a1, ga1 )
    ua = ( - 1 ) ** ( n - 1 ) / ( rn * ga1 )
    ub = rn1 / ga * x ** ( - n )
  else
    a0 = a + n
    a1 = a0
    a2 = a
    call gamma ( a1, ga1 )
    ua = ( - 1 ) ** ( n - 1 ) / ( rn * ga ) * x ** n
    ub = rn1 / ga1
  end if

  hm1 = 1.0e+00_real64
  r = 1.0e+00_real64
  hmax = 0.0e+00_real64
  hmin = 1.0e+300_real64
  h0 = 0.0e+00_real64

  do k = 1, 150
    r = r * ( a0 + k - 1.0e+00_real64 ) * x / ( ( n + k ) * k )
    hm1 = hm1 + r
    hu1 = abs ( hm1 )
    hmax = max ( hmax, hu1 )
    hmin = min ( hmin, hu1 )
    if ( abs ( hm1 - h0 ) < abs ( hm1 ) * 1.0e-15_real64 ) then
      exit
    end if
    h0 = hm1
  end do

  da1 = log10 ( hmax )
  if ( hmin /= 0.0e+00_real64 ) then
    da2 = log10 ( hmin )
  end if
  id = 15 - int ( abs ( da1 - da2 ) )
  hm1 = hm1 * log ( x )
  s0 = 0.0e+00_real64
  do m = 1, n
    if ( 0.0e+00_real64 <= b ) then
      s0 = s0 - 1.0e+00_real64 / m
    else
      s0 = s0 + ( 1.0e+00_real64 - a ) / ( m * ( a + m - 1.0e+00_real64 ) )
    end if
  end do
  hm2 = ps + 2.0e+00_real64 * el + s0
  r = 1.0e+00_real64
  hmax = 0.0e+00_real64
  hmin = 1.0e+300_real64
  h0 = 0.0e+00_real64
  do k = 1, 150
    s1 = 0.0e+00_real64
    s2 = 0.0e+00_real64
    if ( 0.0e+00_real64 < b ) then
      do m = 1, k
        s1 = s1 - ( m + 2.0e+00_real64 * a - 2.0e+00_real64 ) / ( m * ( m + a - 1.0e+00_real64 ) )
      end do
      do m = 1, n
        s2 = s2 + 1.0e+00_real64 / ( k + m )
      end do
    else
      do m = 1, k + n
        s1 = s1 + ( 1.0e+00_real64 - a ) / ( m * ( m + a - 1.0e+00_real64 ) )
      end do
      do m = 1, k
        s2 = s2 + 1.0e+00_real64 / m
      end do
    end if
    hw = 2.0e+00_real64 * el + ps + s1 - s2
    r = r * ( a0 + k - 1.0e+00_real64 ) * x / ( ( n + k ) * k )
    hm2 = hm2 + r * hw
    hu2 = abs ( hm2 )
    hmax = max ( hmax, hu2 )
    hmin = min ( hmin, hu2 )

    if ( abs ( ( hm2 - h0 ) / hm2 ) < 1.0e-15_real64 ) then
      exit
    end if

    h0 = hm2

  end do

  db1 = log10 ( hmax )
  if ( hmin /= 0.0e+00_real64 ) then
    db2 = log10 ( hmin )
  end if
  id1 = 15 - int ( abs ( db1 - db2 ) )
  id = min ( id, id1 )

  if ( n == 0 ) then
    hm3 = 0.0e+00_real64
  else
    hm3 = 1.0e+00_real64
  end if

  r = 1.0e+00_real64
  do k = 1, n - 1
    r = r * ( a2 + k - 1.0e+00_real64 ) / ( ( k - n ) * k ) * x
    hm3 = hm3 + r
  end do

  sa = ua * ( hm1 + hm2 )
  sb = ub * hm3
  hu = sa + sb

  if ( sa /= 0.0e+00_real64 ) then
    id1 = int ( log10 ( abs ( sa ) ) )
  end if

  if ( hu /= 0.0e+00_real64 ) then
    id2 = int ( log10 ( abs ( hu ) ) )
  end if

  if ( sa * sb < 0.0e+00_real64 ) then
    id = id - abs ( id1 - id2 )
  end if

  return
end subroutine chgubi
!> @brief subroutine chguit.
!> @return None.
!>
!> @param a [in] Argument a.
!> @param b [in] Argument b.
!> @param x [in] Argument x.
!> @param hu [inout] Argument hu.
!> @param id [inout] Argument id.
subroutine chguit ( a, b, x, hu, id )

!*****************************************************************************80
!
!! CHGUIT computes the hypergeometric function using Gauss-Legendre integration.
!
!  Discussion:
!
!    This procedure computes the hypergeometric function U(a,b,x) by
!    using Gaussian-Legendre integration (n = 60)
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    22 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) A, B, parameters.
!
!    Input, real(real64) X, the argument.
!
!    Output, real(real64) HU, U(a,b,z).
!
!    Output, integer ID, the estimated number of significant digits.
!
  implicit none

  real(real64), intent(in) :: a
  real(real64) a1
  real(real64), intent(in) :: b
  real(real64) b1
  real(real64) c
  real(real64) d
  real(real64) f1
  real(real64) f2
  real(real64) g
  real(real64) ga
  real(real64), intent(inout) :: hu
  real(real64) hu0
  real(real64) hu1
  real(real64) hu2
  integer, intent(inout) :: id
  integer j
  integer k
  integer m
  real(real64) s
  real(real64), parameter, dimension ( 30 ) :: t = [&
    0.259597723012478e-01_real64, 0.778093339495366e-01_real64, &
    0.129449135396945e+00_real64, 0.180739964873425e+00_real64, &
    0.231543551376029e+00_real64, 0.281722937423262e+00_real64, &
    0.331142848268448e+00_real64, 0.379670056576798e+00_real64, &
    0.427173741583078e+00_real64, 0.473525841761707e+00_real64, &
    0.518601400058570e+00_real64, 0.562278900753945e+00_real64, &
    0.604440597048510e+00_real64, 0.644972828489477e+00_real64, &
    0.683766327381356e+00_real64, 0.720716513355730e+00_real64, &
    0.755723775306586e+00_real64, 0.788693739932264e+00_real64, &
    0.819537526162146e+00_real64, 0.848171984785930e+00_real64, &
    0.874519922646898e+00_real64, 0.898510310810046e+00_real64, &
    0.920078476177628e+00_real64, 0.939166276116423e+00_real64, &
    0.955722255839996e+00_real64, 0.969701788765053e+00_real64, &
    0.981067201752598e+00_real64, 0.989787895222222e+00_real64, &
    0.995840525118838e+00_real64, 0.999210123227436e+00_real64]
  real(real64) t1
  real(real64) t2
  real(real64) t3
  real(real64) t4
  real(real64), parameter, dimension ( 30 ) :: w = [&
    0.519078776312206e-01_real64, 0.517679431749102e-01_real64, &
    0.514884515009810e-01_real64, 0.510701560698557e-01_real64, &
    0.505141845325094e-01_real64, 0.498220356905502e-01_real64, &
    0.489955754557568e-01_real64, 0.480370318199712e-01_real64, &
    0.469489888489122e-01_real64, 0.457343797161145e-01_real64, &
    0.443964787957872e-01_real64, 0.429388928359356e-01_real64, &
    0.413655512355848e-01_real64, 0.396806954523808e-01_real64, &
    0.378888675692434e-01_real64, 0.359948980510845e-01_real64, &
    0.340038927249464e-01_real64, 0.319212190192963e-01_real64, &
    0.297524915007890e-01_real64, 0.275035567499248e-01_real64, &
    0.251804776215213e-01_real64, 0.227895169439978e-01_real64, &
    0.203371207294572e-01_real64, 0.178299010142074e-01_real64, &
    0.152746185967848e-01_real64, 0.126781664768159e-01_real64, &
    0.100475571822880e-01_real64, 0.738993116334531e-02_real64, &
    0.471272992695363e-02_real64, 0.202681196887362e-02_real64]
  real(real64), intent(in) :: x

  id = 7
  a1 = a - 1.0e+00_real64
  b1 = b - a - 1.0e+00_real64
  c = 12.0e+00_real64 / x
  hu0 = 0.0e+00_real64

  do m = 10, 100, 5

    hu1 = 0.0e+00_real64
    g = 0.5e+00_real64 * c / m
    d = g
    do j = 1, m
      s = 0.0e+00_real64
      do k = 1, 30
        t1 = d + g * t(k)
        t2 = d - g * t(k)
        f1 = exp ( - x * t1 ) * t1 ** a1 * ( 1.0e+00_real64 + t1 ) ** b1
        f2 = exp ( - x * t2 ) * t2 ** a1 * ( 1.0e+00_real64 + t2 ) ** b1
        s = s + w(k) * ( f1 + f2 )
      end do
      hu1 = hu1 + s * g
      d = d + 2.0e+00_real64 * g
    end do

    if ( abs ( 1.0e+00_real64 - hu0 / hu1 ) < 1.0e-07_real64 ) then
      exit
    end if

    hu0 = hu1

  end do

  call gamma ( a, ga )
  hu1 = hu1 / ga
  hu0 = 0.0e+00_real64

  do m = 2, 10, 2
    hu2 = 0.0e+00_real64
    g = 0.5e+00_real64 / m
    d = g
    do j = 1, m
      s = 0.0e+00_real64
      do k = 1, 30
        t1 = d + g * t(k)
        t2 = d - g * t(k)
        t3 = c / ( 1.0e+00_real64 - t1 )
        t4 = c / ( 1.0e+00_real64 - t2 ) 
        f1 = t3 * t3 / c * exp ( - x * t3 ) * t3 ** a1 * ( 1.0e+00_real64 + t3 ) ** b1
        f2 = t4 * t4 / c * exp ( - x * t4 ) * t4 ** a1 * ( 1.0e+00_real64 + t4 ) ** b1
        s = s + w(k) * ( f1 + f2 )
      end do
      hu2 = hu2 + s * g
      d = d + 2.0e+00_real64 * g
    end do

    if ( abs ( 1.0e+00_real64 - hu0 / hu2 ) < 1.0e-07_real64 ) then
      exit
    end if

    hu0 = hu2

  end do

  call gamma ( a, ga )
  hu2 = hu2 / ga
  hu = hu1 + hu2

  return
end subroutine chguit
!> @brief subroutine chgul.
!> @return None.
!>
!> @param a [in] Argument a.
!> @param b [in] Argument b.
!> @param x [in] Argument x.
!> @param hu [inout] Argument hu.
!> @param id [inout] Argument id.
pure subroutine chgul ( a, b, x, hu, id )

!*****************************************************************************80
!
!! CHGUL: confluent hypergeometric function U(a,b,x) for large argument X.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    22 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
! 
!  Parameters:
!
!    Input, real(real64) A, B, parameters.
!
!    Input, real(real64) X, the argument.
!
!    Output, real(real64) HU, the value of U(a,b,x).
!
!    Output, integer(int32) ID, the estimated number of 
!    significant digits.
!
  implicit none

  real(real64), intent(in) :: a
  real(real64) aa
  real(real64), intent(in) :: b
  real(real64), intent(inout) :: hu
  integer(int32), intent(inout) :: id
  logical il1
  logical il2
  integer(int32) k
  integer(int32) nm
  real(real64) r
  real(real64) ra
  real(real64) r0
  real(real64), intent(in) :: x

  id = -100
  aa = a - b + 1.0e+00_real64
  il1 = ( a == int ( a ) ) .and. ( a <= 0.0e+00_real64 )
  il2 = ( aa == int ( aa ) ) .and. ( aa <= 0.0e+00_real64 )

  if ( il1 .or. il2 ) then

    if ( il1 ) then
      nm = int ( abs ( a ) )
    end if

    if ( il2 ) then
      nm = int ( abs ( aa ) )
    end if

    hu = 1.0e+00_real64
    r = 1.0e+00_real64
    do k = 1, nm
      r = - r * ( a + k - 1.0e+00_real64 ) * ( a - b + k ) / ( k * x )
      hu = hu + r
    end do
    hu = x ** ( - a ) * hu
    id = 10

  else

    hu = 1.0e+00_real64
    r = 1.0e+00_real64
    do k = 1, 25
      r = - r * ( a + k - 1.0e+00_real64 ) * ( a - b + k ) / ( k * x )
      ra = abs ( r )
      if ( ( 5 < k .and. r0 <= ra ) .or. ra < 1.0e-15_real64 ) then
        exit
      end if
      r0 = ra
      hu = hu + r
    end do

    id = int ( abs ( log10 ( ra ) ) )
    hu = x ** ( - a ) * hu

  end if

  return
end subroutine chgul
!> @brief subroutine chgus.
!> @return None.
!>
!> @param a [in] Argument a.
!> @param b [in] Argument b.
!> @param x [in] Argument x.
!> @param hu [inout] Argument hu.
!> @param id [inout] Argument id.
subroutine chgus ( a, b, x, hu, id )

!*****************************************************************************80
!
!! CHGUS: confluent hypergeometric function U(a,b,x) for small argument X.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    27 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) A, B, parameters.
!
!    Input, real(real64) X, the argument.
!
!    Output, real(real64) HU, U(a,b,x).
!
!    Output, integer(int32) ID, the estimated number of 
!    significant digits.
!
  implicit none

  real(real64), intent(in) :: a
  real(real64), intent(in) :: b
  real(real64) d1
  real(real64) d2
  real(real64) ga
  real(real64) gab
  real(real64) gb
  real(real64) gb2
  real(real64) h0
  real(real64) hmax
  real(real64) hmin
  real(real64), intent(inout) :: hu
  real(real64) hu0
  real(real64) hua
  integer(int32), intent(inout) :: id
  integer(int32) j
  real(real64) pi
  real(real64) r1
  real(real64) r2
  real(real64), intent(in) :: x
  real(real64) xg1
  real(real64) xg2

  id = -100
  pi = 3.141592653589793e+00_real64
  call gamma ( a, ga )
  call gamma ( b, gb )
  xg1 = 1.0e+00_real64 + a - b
  call gamma ( xg1, gab )
  xg2 = 2.0e+00_real64 - b
  call gamma ( xg2, gb2 )
  hu0 = pi / sin ( pi * b )
  r1 = hu0 / ( gab * gb )
  r2 = hu0 * x ** ( 1.0e+00_real64 - b ) / ( ga * gb2 )
  hu = r1 - r2
  hmax = 0.0e+00_real64
  hmin = 1.0e+300_real64
  h0 = 0.0e+00_real64
  do j = 1, 150
    r1 = r1 * ( a + j - 1.0e+00_real64 ) / ( j * ( b + j - 1.0e+00_real64 ) ) * x
    r2 = r2 * ( a - b + j ) / ( j * ( 1.0e+00_real64 - b + j ) ) * x
    hu = hu + r1 - r2
    hua = abs ( hu )
    hmax = max ( hmax, hua )
    hmin = min ( hmin, hua )
    if ( abs ( hu - h0 ) < abs ( hu ) * 1.0e-15_real64 ) then
      exit
    end if
    h0 = hu
  end do

  d1 = log10 ( hmax )
  if ( hmin /= 0.0e+00_real64 ) then
    d2 = log10 ( hmin )
  end if
  id = 15 - int ( abs ( d1 - d2 ) )

  return
end subroutine chgus
!> @brief subroutine cik01.
!> @return None.
!>
!> @param z [in] Argument z.
!> @param cbi0 [inout] Argument cbi0.
!> @param cdi0 [inout] Argument cdi0.
!> @param cbi1 [inout] Argument cbi1.
!> @param cdi1 [inout] Argument cdi1.
!> @param cbk0 [inout] Argument cbk0.
!> @param cdk0 [inout] Argument cdk0.
!> @param cbk1 [inout] Argument cbk1.
!> @param cdk1 [inout] Argument cdk1.
pure subroutine cik01 ( z, cbi0, cdi0, cbi1, cdi1, cbk0, cdk0, cbk1, cdk1 )

!*****************************************************************************80
!
!! CIK01: modified Bessel I0(z), I1(z), K0(z) and K1(z) for complex argument.
!
!  Discussion:
!
!    This procedure computes the modified Bessel functions I0(z), I1(z), 
!    K0(z), K1(z), and their derivatives for a complex argument.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    31 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
! 
!  Parameters:
!
!    Input, complex ( kind = real64 ) Z, the argument.
!
!    Output, complex ( kind = real64 ) CBI0, CDI0, CBI1, CDI1, CBK0, CDK0, CBK1, 
!    CDK1, the values of I0(z), I0'(z), I1(z), I1'(z), K0(z), K0'(z), K1(z), 
!    and K1'(z).
!
  implicit none

  real(real64), parameter, dimension ( 12 ) :: a = [&
    0.125e+00_real64,           7.03125e-02_real64,&
    7.32421875e-02_real64,      1.1215209960938e-01_real64,&
    2.2710800170898e-01_real64, 5.7250142097473e-01_real64,&
    1.7277275025845e+00_real64, 6.0740420012735e+00_real64,&
    2.4380529699556e+01_real64, 1.1001714026925e+02_real64,&
    5.5133589612202e+02_real64, 3.0380905109224e+03_real64]
  real(real64) a0
  real(real64), parameter, dimension ( 10 ) :: a1 = [&
    0.125e+00_real64,            0.2109375e+00_real64, &
    1.0986328125e+00_real64,     1.1775970458984e+01_real64, &
    2.1461706161499e+002_real64, 5.9511522710323e+03_real64, &
    2.3347645606175e+05_real64,  1.2312234987631e+07_real64, &
    8.401390346421e+08_real64,   7.2031420482627e+10_real64]
  real(real64), parameter, dimension ( 12 ) :: b = [&
   -0.375e+00_real64,           -1.171875e-01_real64, &
   -1.025390625e-01_real64,     -1.4419555664063e-01_real64, &
   -2.7757644653320e-01_real64, -6.7659258842468e-01_real64, &
   -1.9935317337513e+00_real64, -6.8839142681099e+00_real64, &
   -2.7248827311269e+01_real64, -1.2159789187654e+02_real64, &
   -6.0384407670507e+02_real64, -3.3022722944809e+03_real64]
  complex(real64) ca
  complex(real64) cb
  complex(real64), intent(inout) :: cbi0
  complex(real64), intent(inout) :: cbi1
  complex(real64), intent(inout) :: cbk0
  complex(real64), intent(inout) :: cbk1
  complex(real64), intent(inout) :: cdi0
  complex(real64), intent(inout) :: cdi1
  complex(real64), intent(inout) :: cdk0
  complex(real64), intent(inout) :: cdk1
  complex(real64) ci
  complex(real64) cr
  complex(real64) cs
  complex(real64) ct
  complex(real64) cw
  integer(int32) k
  integer(int32) k0
  real(real64) pi
  real(real64) w0
  complex(real64), intent(in) :: z
  complex(real64) z1
  complex(real64) z2
  complex(real64) zr
  complex(real64) zr2

  pi = 3.141592653589793e+00_real64
  ci = cmplx ( 0.0e+00_real64, 1.0e+00_real64, kind = real64 )
  a0 = abs ( z )
  z2 = z * z
  z1 = z

  if ( a0 == 0.0e+00_real64 ) then
    cbi0 = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    cbi1 = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
    cdi0 = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
    cdi1 = cmplx ( 0.5e+00_real64, 0.0e+00_real64, kind = real64 )
    cbk0 = cmplx ( 1.0e+30_real64, 0.0e+00_real64, kind = real64 )
    cbk1 = cmplx ( 1.0e+30_real64, 0.0e+00_real64, kind = real64 )
    cdk0 = - cmplx ( 1.0e+30_real64, 0.0e+00_real64, kind = real64 )
    cdk1 = - cmplx ( 1.0e+30_real64, 0.0e+00_real64, kind = real64 )
    return
  end if

  if ( real ( z, kind = real64 ) < 0.0e+00_real64 ) then
    z1 = -z
  end if

  if ( a0 <= 18.0e+00_real64 ) then

    cbi0 = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    cr = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    do k = 1, 50
      cr = 0.25e+00_real64 * cr * z2 / ( k * k )
      cbi0 = cbi0 + cr
      if ( abs ( cr / cbi0 ) < 1.0e-15_real64 ) then
        exit
      end if
    end do

    cbi1 = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    cr = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    do k = 1, 50
      cr = 0.25e+00_real64 * cr * z2 / ( k * ( k + 1 ) )
      cbi1 = cbi1 + cr
      if ( abs ( cr / cbi1 ) < 1.0e-15_real64 ) then
        exit
      end if
    end do

    cbi1 = 0.5e+00_real64 * z1 * cbi1

  else

    if ( a0 < 35.0e+00_real64 ) then
      k0 = 12
    else if ( a0 < 50.0e+00_real64 ) then
      k0 = 9
    else
      k0 = 7
    end if

    ca = exp ( z1 ) / sqrt ( 2.0e+00_real64 * pi * z1 )
    cbi0 = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    zr = 1.0e+00_real64 / z1
    do k = 1, k0
      cbi0 = cbi0 + a(k) * zr ** k
    end do
    cbi0 = ca * cbi0
    cbi1 = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    do k = 1, k0
      cbi1 = cbi1 + b(k) * zr ** k
    end do
    cbi1 = ca * cbi1

  end if

  if ( a0 <= 9.0e+00_real64 ) then

    cs = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
    cw = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
    ct = - log ( 0.5e+00_real64 * z1 ) - 0.5772156649015329e+00_real64
    w0 = 0.0e+00_real64
    cr = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    do k = 1, 50
      w0 = w0 + 1.0e+00_real64 / k
      cr = 0.25e+00_real64 * cr / ( k * k ) * z2
      cs = cs + cr * ( w0 + ct )
      if ( abs ( ( cs - cw ) / cs ) < 1.0e-15_real64 ) then
        exit
      end if
      cw = cs
    end do

    cbk0 = ct + cs

  else

    cb = 0.5e+00_real64 / z1
    zr2 = 1.0e+00_real64 / z2
    cbk0 = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    do k = 1, 10
      cbk0 = cbk0 + a1(k) * zr2 ** k
    end do
    cbk0 = cb * cbk0 / cbi0

  end if

  cbk1 = ( 1.0e+00_real64 / z1 - cbi1 * cbk0 ) / cbi0

  if ( real ( z, kind = real64 ) < 0.0e+00_real64 ) then

    if ( imag ( z ) < 0.0e+00_real64 ) then
      cbk0 = cbk0 + ci * pi * cbi0
      cbk1 = - cbk1 + ci * pi * cbi1
    else
      cbk0 = cbk0 - ci * pi * cbi0
      cbk1 = - cbk1 - ci * pi * cbi1
    end if

    cbi1 = - cbi1

  end if

  cdi0 = cbi1
  cdi1 = cbi0 - 1.0e+00_real64 / z * cbi1
  cdk0 = - cbk1
  cdk1 = - cbk0 - 1.0e+00_real64 / z * cbk1

  return
end subroutine cik01
!> @brief subroutine ciklv.
!> @return None.
!>
!> @param v [in] Argument v.
!> @param z [in] Argument z.
!> @param cbiv [inout] Argument cbiv.
!> @param cdiv [inout] Argument cdiv.
!> @param cbkv [inout] Argument cbkv.
!> @param cdkv [inout] Argument cdkv.
subroutine ciklv ( v, z, cbiv, cdiv, cbkv, cdkv )

!*****************************************************************************80
!
!! CIKLV: modified Bessel functions Iv(z), Kv(z), complex argument, large order.
!
!  Discussion:
!
!    This procedure computes modified Bessel functions Iv(z) and
!    Kv(z) and their derivatives with a complex argument and a large order.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    31 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) V, the order of Iv(z) and Kv(z).
!
!    Input, complex ( kind = real64 ) Z, the argument.
!
!    Output, real(real64) CBIV, CDIV, CBKV, CDKV, the values of
!    Iv(z), Iv'(z), Kv(z), Kv'(z).
!
  implicit none

  real(real64) a(91)
  complex(real64), intent(inout) :: cbiv
  complex(real64), intent(inout) :: cbkv
  complex(real64), intent(inout) :: cdiv
  complex(real64), intent(inout) :: cdkv
  complex(real64) ceta
  complex(real64) cf(12)
  complex(real64) cfi
  complex(real64) cfk
  complex(real64) csi
  complex(real64) csk
  complex(real64) ct
  complex(real64) ct2
  complex(real64) cws
  integer(int32) i
  integer(int32) k
  integer(int32) km
  integer(int32) l
  integer(int32) l0
  integer(int32) lf
  real(real64) pi
  real(real64), intent(in) :: v
  real(real64) v0
  real(real64) vr
  complex(real64), intent(in) :: z

  pi = 3.141592653589793e+00_real64
  km = 12
  call cjk ( km, a )

  do l = 1, 0, -1

    v0 = v - l
    cws = sqrt ( 1.0e+00_real64 + ( z / v0 ) * ( z / v0 ) )
    ceta = cws + log ( z / v0 / ( 1.0e+00_real64 + cws ) )
    ct = 1.0e+00_real64 / cws
    ct2 = ct * ct
    do k = 1, km
      l0 = k * ( k + 1 ) / 2 + 1
      lf = l0 + k
      cf(k) = a(lf)
      do i = lf - 1, l0, -1
        cf(k) = cf(k) * ct2 + a(i)
      end do
      cf(k) = cf(k) * ct ** k
    end do
    vr = 1.0e+00_real64 / v0
    csi = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    do k = 1, km
      csi = csi + cf(k) * vr ** k
    end do
    cbiv = sqrt ( ct / ( 2.0e+00_real64 * pi * v0 ) ) * exp ( v0 * ceta ) * csi
    if ( l == 1 ) then
      cfi = cbiv
    end if
    csk = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    do k = 1, km
      csk = csk + ( - 1 ) ** k * cf(k) * vr ** k
    end do
    cbkv = sqrt ( pi * ct / ( 2.0e+00_real64 * v0 ) ) * exp ( - v0 * ceta ) * csk

    if ( l == 1 ) then
      cfk = cbkv
    end if

  end do

  cdiv =   cfi - v / z * cbiv
  cdkv = - cfk - v / z * cbkv

  return
end subroutine ciklv
!> @brief subroutine cikna.
!> @return None.
!>
!> @param n [in] Argument n.
!> @param z [in] Argument z.
!> @param nm [inout] Argument nm.
!> @param cbi [inout] Argument cbi.
!> @param cdi [inout] Argument cdi.
!> @param cbk [inout] Argument cbk.
!> @param cdk [inout] Argument cdk.
subroutine cikna ( n, z, nm, cbi, cdi, cbk, cdk )

!*****************************************************************************80
!
!! CIKNA: modified Bessel functions In(z), Kn(z), derivatives, complex argument.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    30 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) N, the order of In(z) and Kn(z).
!
!    Input, complex ( kind = real64 ) Z, the argument.
!
!    Output, integer(int32) NM, the highest order computed.
!
!    Output, complex ( kind = real64 ) CBI((0:N), CDI(0:N), CBK(0:N), CDK(0:N), 
!    the values of In(z), In'(z), Kn(z), Kn'(z).
!
  implicit none

  integer(int32), intent(in) :: n

  real(real64) a0
  complex(real64), intent(inout) :: cbi(0:n)
  complex(real64) cbi0
  complex(real64) cbi1
  complex(real64), intent(inout) :: cbk(0:n)
  complex(real64) cbk0
  complex(real64) cbk1
  complex(real64), intent(inout) :: cdi(0:n)
  complex(real64) cdi0
  complex(real64) cdi1
  complex(real64), intent(inout) :: cdk(0:n)
  complex(real64) cdk0
  complex(real64) cdk1
  complex(real64) cf
  complex(real64) cf1
  complex(real64) cf2
  complex(real64) ckk
  complex(real64) cs
  integer(int32) k
  integer(int32) m
  ! integer(int32) msta1
  ! integer(int32) msta2
  integer(int32), intent(inout) :: nm
  complex(real64), intent(in) :: z

  a0 = abs ( z )
  nm = n

  if ( a0 < 1.0e-100_real64 ) then
    do k = 0, n
      cbi(k) = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
      cdi(k) = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
      cbk(k) = - cmplx ( 1.0e+30_real64, 0.0e+00_real64, kind = real64 )
      cdk(k) = cmplx ( 1.0e+30_real64, 0.0e+00_real64, kind = real64 )
    end do
    cbi(0) = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    cdi(1) = cmplx ( 0.5e+00_real64, 0.0e+00_real64, kind = real64 )
    return
  end if

  call cik01 ( z, cbi0, cdi0, cbi1, cdi1, cbk0, cdk0, cbk1, cdk1 )

  cbi(0) = cbi0
  cbi(1) = cbi1
  cbk(0) = cbk0
  cbk(1) = cbk1
  cdi(0) = cdi0
  cdi(1) = cdi1
  cdk(0) = cdk0
  cdk(1) = cdk1

  if ( n <= 1 ) then
    return
  end if

  m = msta1 ( a0, 200 )

  if ( m < n ) then
    nm = m
  else
    m = msta2 ( a0, n, 15 )
  end if

  cf2 = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
  cf1 = cmplx ( 1.0e-30_real64, 0.0e+00_real64, kind = real64 )
  do k = m, 0, -1
    cf = 2.0e+00_real64 * ( k + 1.0e+00_real64 ) / z * cf1 + cf2
    if ( k <= nm ) then
      cbi(k) = cf
    end if
    cf2 = cf1
    cf1 = cf
  end do

  cs = cbi0 / cf
  do k = 0, nm
    cbi(k) = cs * cbi(k)
  end do

  do k = 2, nm
    if ( abs ( cbi(k-2) ) < abs ( cbi(k-1) ) ) then
      ckk = ( 1.0e+00_real64 / z - cbi(k) * cbk(k-1) ) / cbi(k-1)
    else
      ckk = ( cbi(k) * cbk(k-2) + 2.0e+00_real64 * ( k - 1.0e+00_real64 ) &
        / ( z * z ) ) / cbi(k-2)
    end if
    cbk(k) = ckk
  end do

  do k = 2, nm
    cdi(k) =   cbi(k-1) - k / z * cbi(k)
    cdk(k) = - cbk(k-1) - k / z * cbk(k)
  end do

  return
end subroutine cikna
!> @brief subroutine ciknb.
!> @return None.
!>
!> @param n [in] Argument n.
!> @param z [in] Argument z.
!> @param nm [inout] Argument nm.
!> @param cbi [inout] Argument cbi.
!> @param cdi [inout] Argument cdi.
!> @param cbk [inout] Argument cbk.
!> @param cdk [inout] Argument cdk.
subroutine ciknb ( n, z, nm, cbi, cdi, cbk, cdk )

!*****************************************************************************80
!
!! CIKNB computes complex modified Bessel functions In(z) and Kn(z).
!
!  Discussion:
!
!    This procedure also evaluates the derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    30 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) N, the order of In(z) and Kn(z).
!
!    Input, complex ( kind = real64 ) Z, the argument.
!
!    Output, integer(int32) NM, the highest order computed.
!
!    Output, complex ( kind = real64 ) CB((0:N), CDI(0:N), CBK(0:N), CDK(0:N), 
!    the values of In(z), In'(z), Kn(z), Kn'(z).
!
  implicit none

  integer(int32), intent(in) :: n

  real(real64) a0
  complex(real64) ca0
  complex(real64), intent(inout) :: cbi(0:n)
  complex(real64) cbkl
  complex(real64) cbs
  complex(real64), intent(inout) :: cdi(0:n)
  complex(real64), intent(inout) :: cbk(0:n)
  complex(real64), intent(inout) :: cdk(0:n)
  complex(real64) cf
  complex(real64) cf0
  complex(real64) cf1
  complex(real64) cg
  complex(real64) cg0
  complex(real64) cg1
  complex(real64) ci
  complex(real64) cr
  complex(real64) cs0
  complex(real64) csk0
  real(real64) el
  real(real64) fac
  integer(int32) k
  integer(int32) k0
  integer(int32) l
  integer(int32) m
  ! integer(int32) msta1
  ! integer(int32) msta2
  integer(int32), intent(inout) :: nm
  real(real64) pi
  real(real64) vt
  complex(real64), intent(in) :: z
  complex(real64) z1

  pi = 3.141592653589793e+00_real64
  el = 0.57721566490153e+00_real64
  a0 = abs ( z )
  nm = n

  if ( a0 < 1.0e-100_real64 ) then
    do k = 0, n
      cbi(k) = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
      cbk(k) = cmplx ( 1.0e+30_real64, 0.0e+00_real64, kind = real64 )
      cdi(k) = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
      cdk(k) = - cmplx ( 1.0e+30_real64, 0.0e+00_real64, kind = real64 )
    end do
    cbi(0) = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    cdi(1) = cmplx ( 0.5e+00_real64, 0.0e+00_real64, kind = real64 ) 
    return
  end if

  ci = cmplx ( 0.0e+00_real64, 1.0e+00_real64, kind = real64 )

  if ( real ( z, kind = real64 ) < 0.0e+00_real64 ) then
    z1 = -z
  else
    z1 = z
  end if

  if ( n == 0 ) then
    nm = 1
  end if

  m = msta1 ( a0, 200 )

  if ( m < nm ) then
    nm = m
  else
    m = msta2 ( a0, nm, 15 )
  end if

  cbs = 0.0e+00_real64
  csk0 = 0.0e+00_real64
  cf0 = 0.0e+00_real64
  cf1 = 1.0e-100_real64

  do k = m, 0, -1
    cf = 2.0e+00_real64 * ( k + 1.0e+00_real64 ) * cf1 / z1 + cf0
    if ( k <= nm ) then
      cbi(k) = cf
    end if
    if ( k /= 0 .and. k == 2 * int ( k / 2 ) ) then
      csk0 = csk0 + 4.0e+00_real64 * cf / k
    end if
    cbs = cbs + 2.0e+00_real64 * cf
    cf0 = cf1
    cf1 = cf
  end do

  cs0 = exp ( z1 ) / ( cbs - cf )

  do k = 0, nm
    cbi(k) = cs0 * cbi(k)
  end do

  if ( a0 <= 9.0e+00_real64 ) then

    cbk(0) = - ( log ( 0.5e+00_real64 * z1 ) + el ) * cbi(0) + cs0 * csk0
    cbk(1) = ( 1.0e+00_real64 / z1 - cbi(1) * cbk(0) ) / cbi(0)

  else

    ca0 = sqrt ( pi / ( 2.0e+00_real64 * z1 ) ) * exp ( -z1 )

    if ( a0 < 25.0e+00_real64 ) then
      k0 = 16
    else if ( a0 < 80.0e+00_real64 ) then
      k0 = 10
    else if ( a0 < 200.0e+00_real64 ) then
      k0 = 8
    else
      k0 = 6
    end if

    do l = 0, 1
      cbkl = 1.0e+00_real64
      vt = 4.0e+00_real64 * l
      cr = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
      do k = 1, k0
        cr = 0.125e+00_real64 * cr &
          * ( vt - ( 2.0e+00_real64 * k - 1.0e+00_real64 ) ** 2 ) / ( k * z1 )
        cbkl = cbkl + cr
      end do
      cbk(l) = ca0 * cbkl
    end do
  end if

  cg0 = cbk(0)
  cg1 = cbk(1)
  do k = 2, nm
    cg = 2.0e+00_real64 * ( k - 1.0e+00_real64 ) / z1 * cg1 + cg0
    cbk(k) = cg
    cg0 = cg1
    cg1 = cg
  end do

  if ( real ( z, kind = real64 ) < 0.0e+00_real64 ) then
    fac = 1.0e+00_real64
    do k = 0, nm
      if ( imag ( z ) < 0.0e+00_real64 ) then
        cbk(k) = fac * cbk(k) + ci * pi * cbi(k)
      else
        cbk(k) = fac * cbk(k) - ci * pi * cbi(k)
      end if
      cbi(k) = fac * cbi(k)
      fac = - fac
    end do
  end if

  cdi(0) = cbi(1)
  cdk(0) = -cbk(1)
  do k = 1, nm
    cdi(k) = cbi(k-1) - k / z * cbi(k)
    cdk(k) = - cbk(k-1) - k / z * cbk(k)
  end do

  return
end subroutine ciknb
!> @brief subroutine cikva.
!> @return None.
!>
!> @param v [in] Argument v.
!> @param z [in] Argument z.
!> @param vm [inout] Argument vm.
!> @param cbi [inout] Argument cbi.
!> @param cdi [inout] Argument cdi.
!> @param cbk [inout] Argument cbk.
!> @param cdk [inout] Argument cdk.
subroutine cikva ( v, z, vm, cbi, cdi, cbk, cdk )

!*****************************************************************************80
!
!! CIKVA: modified Bessel functions Iv(z), Kv(z), arbitrary order, complex.
!
!  Discussion:
!
!    Compute the modified Bessel functions Iv(z), Kv(z)
!    and their derivatives for an arbitrary order and
!    complex argument
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    31 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:       
!
!    Input, real(real64) V, the order of the functions.
!
!    Input, complex ( kind = real64 ) Z, the argument.
!
!    Output, real(real64) VM, the highest order computed.
!
!    Output, real(real64) CBI(0:N), CDI(0:N), CBK(0:N), CDK(0:N),
!    the values of In+v0(z), In+v0'(z), Kn+v0(z), Kn+v0'(z).
!
  implicit none

  real(real64) a0
  complex(real64) ca
  complex(real64) ca1
  complex(real64) ca2
  complex(real64) cb
  complex(real64), intent(inout) :: cbi(0:)
  complex(real64) cbi0
  complex(real64), intent(inout) :: cdi(0:)
  complex(real64), intent(inout) :: cbk(0:)
  complex(real64) cbk0
  complex(real64) cbk1
  complex(real64), intent(inout) :: cdk(0:)
  complex(real64) cf
  complex(real64) cf1
  complex(real64) cf2
  complex(real64) cg0
  complex(real64) cg1
  complex(real64) cgk
  complex(real64) ci
  complex(real64) ci0
  complex(real64) cp
  complex(real64) cr
  complex(real64) cr1
  complex(real64) cr2
  complex(real64) cs
  complex(real64) csu
  complex(real64) ct
  complex(real64) cvk
  real(real64) gan
  real(real64) gap
  integer(int32) k
  integer(int32) k0
  integer(int32) m
  ! integer(int32) msta1
  ! integer(int32) msta2
  integer(int32) n
  real(real64) pi
  real(real64) piv
  real(real64), intent(in) :: v
  real(real64) v0
  real(real64) v0n
  real(real64) v0p
  real(real64), intent(inout) :: vm
  real(real64) vt
  real(real64) w0
  real(real64) ws
  real(real64) ws0
  complex(real64), intent(in) :: z
  complex(real64) z1
  complex(real64) z2

  pi = 3.141592653589793e+00_real64
  ci = cmplx ( 0.0e+00_real64, 1.0e+00_real64, kind = real64 )
  a0 = abs ( z )
  z1 = z
  z2 = z * z
  n = int ( v )
  v0 = v - n
  piv = pi * v0
  vt = 4.0e+00_real64 * v0 * v0

  if ( n == 0 ) then
    n = 1
  end if

  if ( a0 < 1.0e-100_real64 ) then

    do k = 0, n
      cbi(k) = 0.0e+00_real64
      cdi(k) = 0.0e+00_real64
      cbk(k) = -1.0e+300_real64
      cdk(k) = 1.0e+300_real64
    end do

    if ( v0 == 0.0e+00_real64 ) then
      cbi(0) = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
      cdi(1) = cmplx ( 0.5e+00_real64, 0.0e+00_real64, kind = real64 )
    end if

    vm = v
    return

  end if

  if ( a0 < 35.0e+00_real64 ) then
    k0 = 14
  else if ( a0 < 50.0e+00_real64 ) then
    k0 = 10
  else
    k0 = 8
  end if

  if ( real ( z, kind = real64 ) < 0.0e+00_real64 ) then
    z1 = -z
  end if

  if ( a0 < 18.0e+00_real64 ) then

    if ( v0 == 0.0e+00_real64 ) then
      ca1 = cmplx (1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    else
      v0p = 1.0e+00_real64 + v0
      call gamma ( v0p, gap )
      ca1 = ( 0.5e+00_real64 * z1 ) ** v0 / gap
    end if

    ci0 = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    cr = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    do k = 1, 50
      cr = 0.25e+00_real64 * cr * z2 / ( k * ( k + v0 ) )
      ci0 = ci0 + cr
      if ( abs ( cr ) < abs ( ci0 ) * 1.0e-15_real64 ) then
        exit
      end if
    end do

    cbi0 = ci0 * ca1

  else

    ca = exp ( z1 ) / sqrt ( 2.0e+00_real64 * pi * z1 )
    cs = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    cr = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    do k = 1, k0
      cr = - 0.125e+00_real64 * cr &
        * ( vt - ( 2.0e+00_real64 * k - 1.0e+00_real64 ) ** 2 ) / ( k * z1 )
      cs = cs + cr
    end do
    cbi0 = ca * cs

  end if

  m = msta1 ( a0, 200 )

  if ( m < n ) then
     n = m
  else
     m = msta2 ( a0, n, 15 )
  end if

  cf2 = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
  cf1 = cmplx ( 1.0e-30_real64, 0.0e+00_real64, kind = real64 )
  do k = m, 0, -1
    cf = 2.0e+00_real64 * ( v0 + k + 1.0e+00_real64 ) / z1 * cf1 + cf2
    if ( k <= n ) then
      cbi(k) = cf
    end if
    cf2 = cf1
    cf1 = cf
  end do

  cs = cbi0 / cf
  do k = 0, n
    cbi(k) = cs * cbi(k)
  end do

  if ( a0 <= 9.0e+00_real64 ) then

    if ( v0 == 0.0e+00_real64 ) then
      ct = - log ( 0.5e+00_real64 * z1 ) - 0.5772156649015329e+00_real64
      cs = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
      w0 = 0.0e+00_real64
      cr = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
      do k = 1, 50
        w0 = w0 + 1.0e+00_real64 / k
        cr = 0.25e+00_real64 * cr / ( k * k ) * z2
        cp = cr * ( w0 + ct )
        cs = cs + cp
        if ( 10 <= k .and. abs ( cp / cs ) < 1.0e-15_real64 ) then
          exit
        end if
      end do

      cbk0 = ct + cs

    else

      v0n = 1.0e+00_real64 - v0
      call gamma ( v0n, gan )
      ca2 = 1.0e+00_real64 / ( gan * ( 0.5e+00_real64 * z1 ) ** v0 )
      ca1 = ( 0.5e+00_real64 * z1 ) ** v0 / gap
      csu = ca2 - ca1
      cr1 = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
      cr2 = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
      do k = 1, 50
        cr1 = 0.25e+00_real64 * cr1 * z2 / ( k * ( k - v0 ) )
        cr2 = 0.25e+00_real64 * cr2 * z2 / ( k * ( k + v0 ) )
        csu = csu + ca2 * cr1 - ca1 * cr2
        ws = abs ( csu )
        if ( 10 <= k .and. abs ( ws - ws0 ) / ws < 1.0e-15_real64 ) then
          exit
        end if
        ws0 = ws
      end do

      cbk0 = 0.5e+00_real64 * pi * csu / sin ( piv )

    end if

  else

    cb = exp ( - z1 ) * sqrt ( 0.5e+00_real64 * pi / z1 )
    cs = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    cr = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    do k = 1, k0
      cr = 0.125e+00_real64 * cr &
        * ( vt - ( 2.0e+00_real64 * k - 1.0e+00_real64 ) ** 2 ) / ( k * z1 )
      cs = cs + cr
    end do
    cbk0 = cb * cs

  end if

  cbk1 = ( 1.0e+00_real64 / z1 - cbi(1) * cbk0 ) / cbi(0)
  cbk(0) = cbk0
  cbk(1) = cbk1
  cg0 = cbk0
  cg1 = cbk1

  do k = 2, n
    cgk = 2.0e+00_real64 * ( v0 + k - 1.0e+00_real64 ) / z1 * cg1 + cg0
    cbk(k) = cgk
    cg0 = cg1
    cg1 = cgk
  end do

  if ( real ( z, kind = real64 ) < 0.0e+00_real64 ) then
    do k = 0, n
      cvk = exp ( ( k + v0 ) * pi * ci )
      if ( imag ( z ) < 0.0e+00_real64 ) then
        cbk(k) = cvk * cbk(k) + pi * ci * cbi(k)
        cbi(k) = cbi(k) / cvk
      else if ( 0.0e+00_real64 < imag ( z ) ) then
        cbk(k) = cbk(k) / cvk - pi * ci * cbi(k)
        cbi(k) = cvk * cbi(k)
      end if
    end do
  end if

  cdi(0) = v0 / z * cbi(0) + cbi(1)
  cdk(0) = v0 / z * cbk(0) - cbk(1)
  do k = 1, n
    cdi(k) = - ( k + v0 ) / z * cbi(k) + cbi(k-1)
    cdk(k) = - ( k + v0 ) / z * cbk(k) - cbk(k-1)
  end do

  vm = n + v0

  return
end subroutine cikva
!> @brief subroutine cikvb.
!> @return None.
!>
!> @param v [in] Argument v.
!> @param z [in] Argument z.
!> @param vm [inout] Argument vm.
!> @param cbi [inout] Argument cbi.
!> @param cdi [inout] Argument cdi.
!> @param cbk [inout] Argument cbk.
!> @param cdk [inout] Argument cdk.
subroutine cikvb ( v, z, vm, cbi, cdi, cbk, cdk )

!*****************************************************************************80
!
!! CIKVB: modified Bessel functions,Iv(z), Kv(z), arbitrary order, complex.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    02 August 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) V, the order of the functions.
!
!    Input, complex ( kind = real64 ) Z, the argument.
!
!    Output, real(real64) VM, the highest order computed.
!
!    Output, real(real64) CBI(0:N), CDI(0:N), CBK(0:N), CDK(0:N),
!    the values of In+v0(z), In+v0'(z), Kn+v0(z), Kn+v0'(z).
!
  implicit none

  real(real64) a0
  complex(real64) ca
  complex(real64) ca1
  complex(real64) ca2
  complex(real64) cb
  complex(real64), intent(inout) :: cbi(0:)
  complex(real64) cbi0
  complex(real64), intent(inout) :: cdi(0:)
  complex(real64), intent(inout) :: cbk(0:)
  complex(real64) cbk0
  complex(real64), intent(inout) :: cdk(0:)
  complex(real64) cf
  complex(real64) cf1
  complex(real64) cf2
  complex(real64) ci
  complex(real64) ci0
  complex(real64) ckk
  complex(real64) cp
  complex(real64) cr
  complex(real64) cr1
  complex(real64) cr2
  complex(real64) cs
  complex(real64) csu
  complex(real64) ct
  complex(real64) cvk
  real(real64) gan
  real(real64) gap
  integer(int32) k
  integer(int32) k0
  integer(int32) m
  ! integer(int32) msta1
  ! integer(int32) msta2
  integer(int32) n
  real(real64) pi
  real(real64) piv
  real(real64), intent(in) :: v
  real(real64) v0
  real(real64) v0n
  real(real64) v0p
  real(real64), intent(inout) :: vm
  real(real64) vt
  real(real64) w0
  complex(real64), intent(in) :: z
  complex(real64) z1
  complex(real64) z2

  z1 = z
  z2 = z * z
  a0 = abs ( z )
  pi = 3.141592653589793e+00_real64
  ci = cmplx ( 0.0e+00_real64, 1.0e+00_real64, kind = real64 )
  n = int ( v )
  v0 = v - n
  piv = pi * v0
  vt = 4.0e+00_real64 * v0 * v0

  if ( n == 0 ) then
    n = 1
  end if

  if ( a0 < 1.0e-100_real64 ) then
    do k = 0, n
      cbi(k) = 0.0e+00_real64
      cdi(k) = 0.0e+00_real64
      cbk(k) = -1.0e+300_real64
      cdk(k) = 1.0e+300_real64
    end do
    if ( v0 == 0.0e+00_real64 ) then
      cbi(0) = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
      cdi(1) = cmplx ( 0.5e+00_real64, 0.0e+00_real64, kind = real64 )
    end if
    vm = v
    return
  end if

  if ( a0 < 35.0e+00_real64 ) then
    k0 = 14
  else if ( a0 < 50.0e+00_real64 ) then
    k0 = 10
  else
    k0 = 8
  end if

  if ( real ( z, kind = real64 ) < 0.0e+00_real64 ) then
    z1 = -z
  end if

  if ( a0 < 18.0e+00_real64 ) then

    if ( v0 == 0.0e+00_real64 ) then
      ca1 = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    else
      v0p = 1.0e+00_real64 + v0
      call gamma ( v0p, gap )
      ca1 = ( 0.5e+00_real64 * z1 ) ** v0 / gap
    end if

    ci0 = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    cr = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    do k = 1, 50
      cr = 0.25e+00_real64 * cr * z2 / ( k * ( k + v0 ) )
      ci0 = ci0 + cr
      if ( abs ( cr / ci0 ) < 1.0e-15_real64 ) then
        exit
      end if
    end do

    cbi0 = ci0 * ca1

  else

    ca = exp ( z1 ) / sqrt ( 2.0e+00_real64 * pi * z1 )
    cs = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    cr = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    do k = 1, k0
      cr = -0.125e+00_real64 * cr &
        * ( vt - ( 2.0e+00_real64 * k - 1.0e+00_real64 ) ** 2 ) / ( k * z1 )
      cs = cs + cr
    end do
    cbi0 = ca * cs

  end if

  m = msta1 ( a0, 200 )
  if ( m < n ) then
    n = m
  else
    m = msta2 ( a0, n, 15 )
  end if

  cf2 = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
  cf1 = cmplx ( 1.0e-30_real64, 0.0e+00_real64, kind = real64 )
  do k = m, 0, -1
    cf = 2.0e+00_real64 * ( v0 + k + 1.0e+00_real64 ) / z1 * cf1 + cf2
    if ( k <= n ) then
      cbi(k) = cf
    end if
    cf2 = cf1
    cf1 = cf
  end do
  cs = cbi0 / cf

  do k = 0, n
    cbi(k) = cs * cbi(k)
  end do

  if ( a0 <= 9.0e+00_real64 ) then

    if ( v0 == 0.0e+00_real64 ) then

      ct = - log ( 0.5e+00_real64 * z1 ) - 0.5772156649015329e+00_real64
      cs = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
      w0 = 0.0e+00_real64
      cr = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
      do k = 1, 50
        w0 = w0 + 1.0e+00_real64 / k
        cr = 0.25e+00_real64 * cr / ( k * k ) * z2
        cp = cr * ( w0 + ct )
        cs = cs + cp
        if ( 10 <= k .and. abs ( cp / cs ) < 1.0e-15_real64 ) then
          exit
        end if
      end do

      cbk0 = ct + cs

    else

      v0n = 1.0e+00_real64 - v0
      call gamma ( v0n, gan )
      ca2 = 1.0e+00_real64 / ( gan * ( 0.5e+00_real64 * z1 ) ** v0 )
      ca1 = ( 0.5e+00_real64 * z1 ) ** v0 / gap
      csu = ca2 - ca1
      cr1 = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
      cr2 = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
      do k = 1, 50
        cr1 = 0.25e+00_real64 * cr1 * z2 / ( k * ( k - v0 ) )
        cr2 = 0.25e+00_real64 * cr2 * z2 / ( k * ( k + v0 ) )
        cp = ca2 * cr1 - ca1 * cr2
        csu = csu + cp
        if ( 10 <= k .and. abs ( cp / csu ) < 1.0e-15_real64 ) then
          exit
        end if
      end do

      cbk0 = 0.5e+00_real64 * pi * csu / sin ( piv )

    end if

  else

    cb = exp ( -z1 ) * sqrt ( 0.5e+00_real64 * pi / z1 )
    cs = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    cr = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    do k = 1, k0
      cr = 0.125e+00_real64 * cr * ( vt - ( 2.0e+00_real64 * k - 1.0e+00_real64 ) ** 2 ) &
        / ( k * z1 )
      cs = cs + cr
    end do

    cbk0 = cb * cs

  end if

  cbk(0) = cbk0

  if ( real ( z, kind = real64 ) < 0.0e+00_real64 ) then
    do k = 0, n
      cvk = exp ( ( k + v0 ) * pi * ci )
      if ( imag ( z ) < 0.0e+00_real64 ) then
        cbk(k) = cvk * cbk(k) + pi * ci * cbi(k)
        cbi(k) = cbi(k) / cvk
      else if ( 0.0e+00_real64 < imag ( z ) ) then
        cbk(k) = cbk(k) / cvk - pi * ci * cbi(k)
        cbi(k) = cvk * cbi(k)
      end if
    end do
  end if

  do k = 1, n
    ckk = ( 1.0e+00_real64 / z - cbi(k) * cbk(k-1) ) / cbi(k-1)
    cbk(k) = ckk
  end do

  cdi(0) = v0 / z * cbi(0) + cbi(1)
  cdk(0) = v0 / z * cbk(0) - cbk(1)
  do k = 1, n
    cdi(k) = - ( k + v0 ) / z * cbi(k) + cbi(k-1)
    cdk(k) = - ( k + v0 ) / z * cbk(k) - cbk(k-1)
  end do 
 
  vm = n + v0

  return
end subroutine cikvb
!> @brief subroutine cisia.
!> @return None.
!>
!> @param x [in] Argument x.
!> @param ci [inout] Argument ci.
!> @param si [inout] Argument si.
pure subroutine cisia ( x, ci, si )

!*****************************************************************************80
!
!! CISIA computes cosine Ci(x) and sine integrals Si(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    03 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) X, the argument of Ci(x) and Si(x).
!
!    Output, real(real64) CI, SI, the values of Ci(x) and Si(x).
!
  implicit none

  real(real64) bj(101)
  real(real64), intent(inout) :: ci
  real(real64) el
  real(real64) eps
  integer(int32) k
  integer(int32) m
  real(real64) p2
  real(real64), intent(inout) :: si
  real(real64), intent(in) :: x
  real(real64) x2
  real(real64) xa
  real(real64) xa0
  real(real64) xa1
  real(real64) xcs
  real(real64) xf
  real(real64) xg
  real(real64) xg1
  real(real64) xg2
  real(real64) xr
  real(real64) xs
  real(real64) xss

  p2 = 1.570796326794897e+00_real64
  el = 0.5772156649015329e+00_real64
  eps = 1.0e-15_real64
  x2 = x * x

  if ( x == 0.0e+00_real64 ) then

    ci = -1.0e+300_real64
    si = 0.0e+00_real64

  else if ( x <= 16.0e+00_real64 ) then

    xr = -0.25e+00_real64 * x2
    ci = el + log ( x ) + xr
    do k = 2, 40
      xr = -0.5e+00_real64 * xr * ( k - 1 ) / ( k * k * ( 2 * k - 1 ) ) * x2
      ci = ci + xr
      if ( abs ( xr ) < abs ( ci ) * eps ) then
        exit
      end if
    end do

    xr = x
    si = x
    do k = 1, 40
      xr = -0.5e+00_real64 * xr * ( 2 * k - 1 ) / k / ( 4 * k * k + 4 * k + 1 ) * x2
      si = si + xr
      if ( abs ( xr ) < abs ( si ) * eps ) then
        return
      end if
    end do

  else if ( x <= 32.0e+00_real64 ) then

    m = int ( 47.2e+00_real64 + 0.82e+00_real64 * x )
    xa1 = 0.0e+00_real64
    xa0 = 1.0e-100_real64
    do k = m, 1, -1
      xa = 4.0e+00_real64 * k * xa0 / x - xa1
      bj(k) = xa
      xa1 = xa0
      xa0 = xa
    end do
    xs = bj(1)
    do k = 3, m, 2
      xs = xs + 2.0e+00_real64 * bj(k)
    end do
    bj(1) = bj(1) / xs
    do k = 2, m
      bj(k) = bj(k) / xs
    end do
    xr = 1.0e+00_real64
    xg1 = bj(1)
    do k = 2, m
      xr = 0.25e+00_real64 * xr * ( 2.0e+00_real64 * k - 3.0e+00_real64 ) **2 &
        / ( ( k - 1.0e+00_real64 ) * ( 2.0e+00_real64 * k - 1.0e+00_real64 ) ** 2 ) * x
      xg1 = xg1 + bj(k) * xr
    end do

    xr = 1.0e+00_real64
    xg2 = bj(1)
    do k = 2, m
      xr = 0.25e+00_real64 * xr * ( 2.0e+00_real64 * k - 5.0e+00_real64 )**2 &
        / ( ( k-1.0e+00_real64 ) * ( 2.0e+00_real64 * k - 3.0e+00_real64 ) ** 2 ) * x
      xg2 = xg2 + bj(k) * xr
    end do

    xcs = cos ( x / 2.0e+00_real64 )
    xss = sin ( x / 2.0e+00_real64 )
    ci = el + log ( x ) - x * xss * xg1 + 2.0_real64 * xcs * xg2 - 2.0_real64 * xcs * xcs
    si = x * xcs * xg1 + 2.0_real64 * xss * xg2 - sin ( x )

  else

    xr = 1.0e+00_real64
    xf = 1.0e+00_real64
    do k = 1, 9
      xr = -2.0e+00_real64 * xr * k * ( 2 * k - 1 ) / x2
      xf = xf + xr
    end do
    xr = 1.0e+00_real64 / x
    xg = xr
    do k = 1, 8
      xr = -2.0e+00_real64 * xr * ( 2 * k + 1 ) * k / x2
      xg = xg + xr
    end do
    ci = xf * sin ( x ) / x - xg * cos ( x ) / x
    si = p2 - xf * cos ( x ) / x - xg * sin ( x ) / x

  end if

  return
end subroutine cisia
!> @brief subroutine cisib.
!> @return None.
!>
!> @param x [in] Argument x.
!> @param ci [inout] Argument ci.
!> @param si [inout] Argument si.
pure subroutine cisib ( x, ci, si )

!*****************************************************************************80
!
!! CISIB computes cosine and sine integrals.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    20 March 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) X, the argument of Ci(x) and Si(x).
!
!    Output, real(real64) CI, SI, the values of Ci(x) and Si(x).
!
  implicit none

  real(real64), intent(inout) :: ci
  real(real64) fx
  real(real64) gx
  real(real64), intent(inout) :: si
  real(real64), intent(in) :: x
  real(real64) x2

  x2 = x * x

  if ( x == 0.0e+00_real64 ) then

    ci = -1.0e+300_real64
    si = 0.0e+00_real64

  else if ( x <= 1.0e+00_real64 ) then

    ci = (((( -3.0e-08_real64        * x2 &
             + 3.10e-06_real64     ) * x2 &
             - 2.3148e-04_real64   ) * x2 &
             + 1.041667e-02_real64 ) * x2 &
             - 0.25e+00_real64     ) * x2 + 0.577215665e+00_real64 + log ( x )

     si = (((( 3.1e-07_real64        * x2 &
             - 2.834e-05_real64    ) * x2 &
             + 1.66667e-03_real64  ) * x2 &
             - 5.555556e-02_real64 ) * x2 + 1.0e+00_real64 ) * x

  else

    fx = (((( x2              &
      + 38.027264e+00_real64  ) * x2 &
      + 265.187033e+00_real64 ) * x2 &
      + 335.67732e+00_real64  ) * x2 &
      + 38.102495e+00_real64  ) /    &
      (((( x2                 &
      + 40.021433e+00_real64  ) * x2 &
      + 322.624911e+00_real64 ) * x2 &
      + 570.23628e+00_real64  ) * x2 &
      + 157.105423e+00_real64 )

    gx = (((( x2               &
      + 42.242855e+00_real64  ) * x2  &
      + 302.757865e+00_real64 ) * x2  &
      + 352.018498e+00_real64 ) * x2  &
      + 21.821899e+00_real64 ) /      &
      (((( x2                  &
      + 48.196927e+00_real64   ) * x2 &
      + 482.485984e+00_real64  ) * x2 &
      + 1114.978885e+00_real64 ) * x2 &
      + 449.690326e+00_real64  ) / x

    ci = fx * sin ( x ) / x - gx * cos ( x ) / x

    si = 1.570796327e+00_real64 - fx * cos ( x ) / x - gx * sin ( x ) / x

  end if

  return
end subroutine cisib
!> @brief subroutine cjk.
!> @return None.
!>
!> @param km [in] Argument km.
!> @param a [inout] Argument a.
pure subroutine cjk ( km, a )

!*****************************************************************************80
!
!! CJK: asymptotic expansion coefficients for Bessel functions of large order.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    01 August 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) KM, the maximum value of K.
!
!    Output, real(real64) A(L), the value of Cj(k) where j and k are 
!    related to L by L = j+1+[k*(k+1)]/2; j,k = 0,1,...,Km.
!
  implicit none

  real(real64), intent(inout) :: a(:)
  real(real64) f
  real(real64) f0
  real(real64) g
  real(real64) g0
  integer(int32) j
  integer(int32) k
  integer(int32), intent(in) :: km
  integer(int32) l1
  integer(int32) l2
  integer(int32) l3
  integer(int32) l4

  a(1) = 1.0e+00_real64
  f0 = 1.0e+00_real64
  g0 = 1.0e+00_real64
  do k = 0, km - 1
    l1 = ( k + 1 ) * ( k + 2 ) / 2 + 1
    l2 = ( k + 1 ) * ( k + 2 ) / 2 + k + 2
    f = ( 0.5e+00_real64 * k + 0.125e+00_real64 / ( k + 1 ) ) * f0
    g = - ( 1.5e+00_real64 * k + 0.625e+00_real64 &
      / ( 3.0e+00_real64 * ( k + 1.0e+00_real64 ) ) ) * g0
    a(l1) = f
    a(l2) = g
    f0 = f
    g0 = g
  end do

  do k = 1, km - 1
    do j = 1, k
      l3 = k * ( k + 1 ) / 2 + j + 1
      l4 = ( k + 1 ) * ( k + 2 ) / 2 + j + 1
      a(l4) = ( j + 0.5e+00_real64 * k + 0.125e+00_real64 &
        / ( 2.0e+00_real64 * j + k + 1.0e+00_real64 ) ) * a(l3) &
        - ( j + 0.5e+00_real64 * k - 1.0e+00_real64 + 0.625e+00_real64 &
        / ( 2.0e+00_real64 * j + k + 1.0e+00_real64 ) ) * a(l3-1)
    end do
  end do

  return
end subroutine cjk
!> @brief subroutine cjy01.
!> @return None.
!>
!> @param z [in] Argument z.
!> @param cbj0 [inout] Argument cbj0.
!> @param cdj0 [inout] Argument cdj0.
!> @param cbj1 [inout] Argument cbj1.
!> @param cdj1 [inout] Argument cdj1.
!> @param cby0 [inout] Argument cby0.
!> @param cdy0 [inout] Argument cdy0.
!> @param cby1 [inout] Argument cby1.
!> @param cdy1 [inout] Argument cdy1.
pure subroutine cjy01 ( z, cbj0, cdj0, cbj1, cdj1, cby0, cdy0, cby1, cdy1 )

!*****************************************************************************80
!
!! CJY01: complexBessel functions, derivatives, J0(z), J1(z), Y0(z), Y1(z).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    02 August 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
! 
!  Parameters:
!
!    Input, complex ( kind = real64 ) Z, the argument.
!
!    Output, complex ( kind = real64 ) CBJ0, CDJ0, CBJ1, CDJ1, CBY0, CDY0, CBY1, 
!    CDY1, the values of J0(z), J0'(z), J1(z), J1'(z), Y0(z), Y0'(z), 
!    Y1(z), Y1'(z).
!
  implicit none

  real(real64), parameter, dimension ( 12 ) :: a = [&
    -0.703125e-01_real64,0.112152099609375e+00_real64, &
    -0.5725014209747314e+00_real64,0.6074042001273483e+01_real64, &
    -0.1100171402692467e+03_real64,0.3038090510922384e+04_real64, &
    -0.1188384262567832e+06_real64,0.6252951493434797e+07_real64, &
    -0.4259392165047669e+09_real64,0.3646840080706556e+11_real64, &
    -0.3833534661393944e+13_real64,0.4854014686852901e+15_real64]
  real(real64) a0
  real(real64), parameter, dimension ( 12 ) :: a1 = [&
    0.1171875e+00_real64,-0.144195556640625e+00_real64, &
    0.6765925884246826e+00_real64,-0.6883914268109947e+01_real64, &
    0.1215978918765359e+03_real64,-0.3302272294480852e+04_real64, &
    0.1276412726461746e+06_real64,-0.6656367718817688e+07_real64, &
    0.4502786003050393e+09_real64,-0.3833857520742790e+11_real64, &
    0.4011838599133198e+13_real64,-0.5060568503314727e+15_real64]
  real(real64), parameter, dimension ( 12 ) :: b = [&
    0.732421875e-01_real64,-0.2271080017089844e+00_real64, &
    0.1727727502584457e+01_real64,-0.2438052969955606e+02_real64, &
    0.5513358961220206e+03_real64,-0.1825775547429318e+05_real64, &
    0.8328593040162893e+06_real64,-0.5006958953198893e+08_real64, &
    0.3836255180230433e+10_real64,-0.3649010818849833e+12_real64, &
    0.4218971570284096e+14_real64,-0.5827244631566907e+16_real64]
  real(real64), parameter, dimension ( 12 ) :: b1 = [&
    -0.1025390625e+00_real64,0.2775764465332031e+00_real64, &
    -0.1993531733751297e+01_real64,0.2724882731126854e+02_real64, &
    -0.6038440767050702e+03_real64,0.1971837591223663e+05_real64, &
    -0.8902978767070678e+06_real64,0.5310411010968522e+08_real64, &
    -0.4043620325107754e+10_real64,0.3827011346598605e+12_real64, &
    -0.4406481417852278e+14_real64,0.6065091351222699e+16_real64]
  complex(real64), intent(inout) :: cbj0
  complex(real64), intent(inout) :: cbj1
  complex(real64), intent(inout) :: cby0
  complex(real64), intent(inout) :: cby1
  complex(real64), intent(inout) :: cdj0
  complex(real64), intent(inout) :: cdj1
  complex(real64), intent(inout) :: cdy0
  complex(real64), intent(inout) :: cdy1
  complex(real64) ci
  complex(real64) cp
  complex(real64) cp0
  complex(real64) cp1
  complex(real64) cq0
  complex(real64) cq1
  complex(real64) cr
  complex(real64) cs
  complex(real64) ct1
  complex(real64) ct2
  complex(real64) cu
  real(real64) el
  integer(int32) k
  integer(int32) k0
  real(real64) pi
  real(real64) rp2
  real(real64) w0
  real(real64) w1
  complex(real64), intent(in) :: z
  complex(real64) z1
  complex(real64) z2

  pi = 3.141592653589793e+00_real64
  el = 0.5772156649015329e+00_real64
  rp2 = 2.0e+00_real64 / pi
  ci = cmplx ( 0.0e+00_real64, 1.0e+00_real64, kind = real64 )
  a0 = abs ( z )
  z2 = z * z
  z1 = z

  if ( a0 == 0.0e+00_real64 ) then
    cbj0 = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    cbj1 = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
    cdj0 = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
    cdj1 = cmplx ( 0.5e+00_real64, 0.0e+00_real64, kind = real64 )
    cby0 = - cmplx ( 1.0e+30_real64, 0.0e+00_real64, kind = real64 )
    cby1 = - cmplx ( 1.0e+30_real64, 0.0e+00_real64, kind = real64 )
    cdy0 = cmplx ( 1.0e+30_real64, 0.0e+00_real64, kind = real64 )
    cdy1 = cmplx ( 1.0e+30_real64, 0.0e+00_real64, kind = real64 )
    return
  end if

  if ( real ( z, kind = real64 ) < 0.0e+00_real64 ) then
    z1 = -z
  end if

  if ( a0 <= 12.0e+00_real64 ) then

    cbj0 = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    cr = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    do k = 1, 40
      cr = -0.25e+00_real64 * cr * z2 / ( k * k )
      cbj0 = cbj0 + cr
      if ( abs ( cr ) < abs ( cbj0 ) * 1.0e-15_real64 ) then
        exit
      end if
    end do

    cbj1 = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    cr = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    do k = 1, 40
      cr = -0.25e+00_real64 * cr * z2 / ( k * ( k + 1.0e+00_real64 ) )
      cbj1 = cbj1 + cr
      if ( abs ( cr ) < abs ( cbj1 ) * 1.0e-15_real64 ) then
        exit
      end if
    end do

    cbj1 = 0.5e+00_real64 * z1 * cbj1
    w0 = 0.0e+00_real64
    cr = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    cs = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
    do k = 1, 40
      w0 = w0 + 1.0e+00_real64 / k
      cr = -0.25e+00_real64 * cr / ( k * k ) * z2
      cp = cr * w0
      cs = cs + cp
      if ( abs ( cp ) < abs ( cs ) * 1.0e-15_real64 ) then
        exit
      end if
    end do

    cby0 = rp2 * ( log ( z1 / 2.0e+00_real64 ) + el ) * cbj0 - rp2 * cs
    w1 = 0.0e+00_real64
    cr = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    cs = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    do k = 1, 40
      w1 = w1 + 1.0e+00_real64 / k
      cr = -0.25e+00_real64 * cr / ( k * ( k + 1 ) ) * z2
      cp = cr * ( 2.0e+00_real64 * w1 + 1.0e+00_real64 / ( k + 1.0e+00_real64 ) )
      cs = cs + cp
      if ( abs ( cp ) < abs ( cs ) * 1.0e-15_real64 ) then
        exit
      end if
    end do

    cby1 = rp2 * ( ( log ( z1 / 2.0e+00_real64 ) + el ) * cbj1 &
      - 1.0e+00_real64 / z1 - 0.25e+00_real64 * z1 * cs )

  else

    if ( a0 < 35.0e+00_real64 ) then
      k0 = 12
    else if ( a0 < 50.0e+00_real64 ) then
      k0 = 10
    else
      k0 = 8
    end if

    ct1 = z1 - 0.25e+00_real64 * pi

    cp0 = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    do k = 1, k0
      cp0 = cp0 + a(k) * z1 ** ( - 2 * k )
    end do

    cq0 = -0.125e+00_real64 / z1
    do k = 1, k0
      cq0 = cq0 + b(k) * z1 ** ( - 2 * k - 1 )
    end do

    cu = sqrt ( rp2 / z1 )
    cbj0 = cu * ( cp0 * cos ( ct1 ) - cq0 * sin ( ct1 ) )
    cby0 = cu * ( cp0 * sin ( ct1 ) + cq0 * cos ( ct1 ) )
    ct2 = z1 - 0.75e+00_real64 * pi

    cp1 = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    do k = 1, k0
      cp1 = cp1 + a1(k) * z1 ** ( - 2 * k )
    end do

    cq1 = 0.375e+00_real64 / z1
    do k = 1, k0
      cq1 = cq1 + b1(k) * z1 ** ( - 2 * k - 1 )
    end do

    cbj1 = cu * ( cp1 * cos ( ct2 ) - cq1 * sin ( ct2 ) )
    cby1 = cu * ( cp1 * sin ( ct2 ) + cq1 * cos ( ct2 ) )

  end if

  if ( real ( z, kind = real64 ) < 0.0e+00_real64 ) then
    if ( imag ( z ) < 0.0e+00_real64 ) then
      cby0 =     cby0 - 2.0e+00_real64 * ci * cbj0
      cby1 = - ( cby1 - 2.0e+00_real64 * ci * cbj1 )
    else
      cby0 =     cby0 + 2.0e+00_real64 * ci * cbj0
      cby1 = - ( cby1 + 2.0e+00_real64 * ci * cbj1 )
    end if
    cbj1 = -cbj1
  end if

  cdj0 = -cbj1
  cdj1 = cbj0 - 1.0e+00_real64 / z * cbj1
  cdy0 = -cby1
  cdy1 = cby0 - 1.0e+00_real64 / z * cby1

  return
end subroutine cjy01
!> @brief subroutine cjylv.
!> @return None.
!>
!> @param v [in] Argument v.
!> @param z [in] Argument z.
!> @param cbjv [inout] Argument cbjv.
!> @param cdjv [inout] Argument cdjv.
!> @param cbyv [inout] Argument cbyv.
!> @param cdyv [inout] Argument cdyv.
subroutine cjylv ( v, z, cbjv, cdjv, cbyv, cdyv )

!*****************************************************************************80
!
!! CJYLV: Bessel functions Jv(z), Yv(z) of complex argument and large order v.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    25 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) V, the order of Jv(z) and Yv(z).
!
!    Input, complex ( kind = real64 ) Z, the argument.
!
!    Output, complex ( kind = real64 ) CBJV, CDJV, CBYV, CDYV, the values of Jv(z), 
!    Jv'(z), Yv(z), Yv'(z).
!
  implicit none

  real(real64) a(91)
  complex(real64), intent(inout) :: cbjv
  complex(real64), intent(inout) :: cbyv
  complex(real64), intent(inout) :: cdjv
  complex(real64), intent(inout) :: cdyv
  complex(real64) ceta
  complex(real64) cf(12)
  complex(real64) cfj
  complex(real64) cfy
  complex(real64) csj
  complex(real64) csy
  complex(real64) ct
  complex(real64) ct2
  complex(real64) cws
  integer(int32) i
  integer(int32) k
  integer(int32) km
  integer(int32) l
  integer(int32) l0
  integer(int32) lf
  real(real64) pi
  real(real64), intent(in) :: v
  real(real64) v0
  real(real64) vr
  complex(real64), intent(in) :: z

  km = 12
  call cjk ( km, a )
  pi = 3.141592653589793e+00_real64

  do l = 1, 0, -1

    v0 = v - l
    cws = sqrt ( 1.0e+00_real64 - ( z / v0 ) * ( z / v0 ) )
    ceta = cws + log ( z / v0 / ( 1.0e+00_real64 + cws ) )
    ct = 1.0e+00_real64 / cws
    ct2 = ct * ct

    do k = 1, km
      l0 = k * ( k + 1 ) / 2 + 1
      lf = l0 + k
      cf(k) = a(lf)
      do i = lf - 1, l0, -1
        cf(k) = cf(k) * ct2 + a(i)
      end do
      cf(k) = cf(k) * ct ** k
    end do

    vr = 1.0e+00_real64 / v0
    csj = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    do k = 1, km
      csj = csj + cf(k) * vr ** k
    end do
    cbjv = sqrt ( ct / ( 2.0e+00_real64 * pi * v0 ) ) * exp ( v0 * ceta ) * csj
    if ( l == 1 ) then
      cfj = cbjv
    end if
    csy = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    do k = 1, km
      csy = csy + ( -1.0e+00_real64 ) ** k * cf(k) * vr ** k
    end do
    cbyv = - sqrt ( 2.0e+00_real64 * ct / ( pi * v0 ) ) * exp ( - v0 * ceta ) * csy
    if ( l == 1 ) then
      cfy = cbyv
    end if

  end do

  cdjv = - v / z * cbjv + cfj
  cdyv = - v / z * cbyv + cfy

  return
end subroutine cjylv
!> @brief subroutine cjyna.
!> @return None.
!>
!> @param n [in] Argument n.
!> @param z [in] Argument z.
!> @param nm [inout] Argument nm.
!> @param cbj [inout] Argument cbj.
!> @param cdj [inout] Argument cdj.
!> @param cby [inout] Argument cby.
!> @param cdy [inout] Argument cdy.
subroutine cjyna ( n, z, nm, cbj, cdj, cby, cdy )

!*****************************************************************************80
!
!! CJYNA: Bessel functions and derivatives, Jn(z) and Yn(z) of complex argument.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    02 August 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) N, the order of Jn(z) and Yn(z).
!
!    Input, complex ( kind = real64 ) Z, the argument of Jn(z) and Yn(z).
!
!    Output, integer(int32) NM, the highest order computed.
!
!    Output, complex ( kind = real64 ), CBJ(0:N), CDJ(0:N), CBY(0:N), CDY(0:N),
!    the values of Jn(z), Jn'(z), Yn(z), Yn'(z).
!
  implicit none

  integer(int32), intent(in) :: n

  real(real64) a0
  complex(real64), intent(inout) :: cbj(0:n)
  complex(real64) cbj0
  complex(real64) cbj1
  complex(real64), intent(inout) :: cby(0:n)
  complex(real64) cby0
  complex(real64) cby1
  complex(real64), intent(inout) :: cdj(0:n)
  complex(real64) cdj0
  complex(real64) cdj1
  complex(real64), intent(inout) :: cdy(0:n)
  complex(real64) cdy0
  complex(real64) cdy1
  complex(real64) cf
  complex(real64) cf1
  complex(real64) cf2
  complex(real64) cg0
  complex(real64) cg1
  complex(real64) ch0
  complex(real64) ch1
  complex(real64) ch2
  complex(real64) cj0
  complex(real64) cj1
  complex(real64) cjk
  complex(real64) cp11
  complex(real64) cp12
  complex(real64) cp21
  complex(real64) cp22
  complex(real64) cs
  complex(real64) cyk
  complex(real64) cyl1
  complex(real64) cyl2
  complex(real64) cylk
  integer(int32) k
  integer(int32) lb
  integer(int32) lb0
  integer(int32) m
  ! integer(int32) msta1
  ! integer(int32) msta2
  integer(int32), intent(inout) :: nm
  real(real64) pi
  real(real64) wa
  real(real64) ya0
  real(real64) ya1
  real(real64) yak
  complex(real64), intent(in) :: z

  pi = 3.141592653589793e+00_real64
  a0 = abs ( z )
  nm = n

  if ( a0 < 1.0e-100_real64 ) then
    do k = 0, n
      cbj(k) = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
      cdj(k) = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
      cby(k) = - cmplx ( 1.0e+30_real64, 0.0e+00_real64, kind = real64 )
      cdy(k) = cmplx ( 1.0e+30_real64, 0.0e+00_real64, kind = real64 )
    end do
    cbj(0) = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    cdj(1) = cmplx ( 0.5e+00_real64, 0.0e+00_real64, kind = real64 )
    return
  end if

  call cjy01 ( z, cbj0, cdj0, cbj1, cdj1, cby0, cdy0, cby1, cdy1 )
  cbj(0) = cbj0
  cbj(1) = cbj1
  cby(0) = cby0
  cby(1) = cby1
  cdj(0) = cdj0
  cdj(1) = cdj1
  cdy(0) = cdy0
  cdy(1) = cdy1

  if ( n <= 1 ) then
    return
  end if

  if ( n < int ( 0.25e+00_real64 * a0 ) ) then

    cj0 = cbj0
    cj1 = cbj1
    do k = 2, n
      cjk = 2.0e+00_real64 * ( k - 1.0e+00_real64 ) / z * cj1 - cj0
      cbj(k) = cjk
      cj0 = cj1
      cj1 = cjk
    end do

  else

    m = msta1 ( a0, 200 )

    if ( m < n ) then
      nm = m
    else
      m = msta2 ( a0, n, 15 )
    end if

    cf2 = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 ) 
    cf1 = cmplx ( 1.0e-30_real64, 0.0e+00_real64, kind = real64 )
    do k = m, 0, -1
      cf = 2.0e+00_real64 * ( k + 1.0e+00_real64 ) / z * cf1 - cf2
      if ( k <= nm ) then
        cbj(k) = cf
      end if
      cf2 = cf1
      cf1 = cf
    end do

    if ( abs ( cbj1 ) < abs ( cbj0 ) ) then
      cs = cbj0 / cf
    else
      cs = cbj1 / cf2
    end if

    do k = 0, nm
      cbj(k) = cs * cbj(k)
    end do

  end if

  do k = 2, nm
    cdj(k) = cbj(k-1) - k / z * cbj(k)
  end do
  ya0 = abs ( cby0 )
  lb = 0
  cg0 = cby0
  cg1 = cby1
  do k = 2, nm
    cyk = 2.0e+00_real64 * ( k - 1.0e+00_real64 ) / z * cg1 - cg0
    if ( abs ( cyk ) <= 1.0e+290_real64 ) then         
      yak = abs ( cyk )
      ya1 = abs ( cg0 )
      if ( yak < ya0 .and. yak < ya1 ) then
        lb = k
      end if
      cby(k) = cyk
      cg0 = cg1
      cg1 = cyk
    end if
  end do

  if ( 4 < lb  .and. imag ( z ) /= 0.0e+00_real64 ) then

    do

      if ( lb == lb0 ) then
        exit
      end if

      ch2 = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
      ch1 = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
      lb0 = lb
      do k = lb, 1, -1
        ch0 = 2.0e+00_real64 * k / z * ch1 - ch2
        ch2 = ch1
        ch1 = ch0
      end do
      cp12 = ch0
      cp22 = ch2
      ch2 = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
      ch1 = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
      do k = lb, 1, -1
        ch0 = 2.0e+00_real64 * k / z * ch1 - ch2
        ch2 = ch1
        ch1 = ch0
      end do
      cp11 = ch0
      cp21 = ch2

      if ( lb == nm ) then
        cbj(lb+1) = 2.0e+00_real64 * lb / z * cbj(lb) - cbj(lb-1)
      end if

      if ( abs ( cbj(1) ) < abs ( cbj(0) ) ) then
        cby(lb+1) = ( cbj(lb+1) * cby0 - 2.0e+00_real64 * cp11 / ( pi * z ) ) / cbj(0)
        cby(lb) = ( cbj(lb) * cby0 + 2.0e+00_real64 * cp12 / ( pi * z ) ) / cbj(0)
      else
        cby(lb+1) = ( cbj(lb+1) * cby1 - 2.0e+00_real64 * cp21 / ( pi * z ) ) / cbj(1)
        cby(lb) = ( cbj(lb) * cby1 + 2.0e+00_real64 * cp22 / ( pi * z ) )  / cbj(1)
      end if

      cyl2 = cby(lb+1)
      cyl1 = cby(lb)
      do k = lb - 1, 0, -1
        cylk = 2.0e+00_real64 * ( k + 1.0e+00_real64 ) / z * cyl1 - cyl2
        cby(k) = cylk
        cyl2 = cyl1
        cyl1 = cylk
      end do

      cyl1 = cby(lb)
      cyl2 = cby(lb+1)
      do k = lb + 1, nm - 1
        cylk = 2.0e+00_real64 * k / z * cyl2 - cyl1
        cby(k+1) = cylk
        cyl1 = cyl2
        cyl2 = cylk
      end do

      do k = 2, nm
        wa = abs ( cby(k) )
        if ( wa < abs ( cby(k-1) ) ) then
          lb = k
        end if
      end do

    end do

  end if

  do k = 2, nm
    cdy(k) = cby(k-1) - k / z * cby(k)
  end do

  return
end subroutine cjyna
!> @brief subroutine cjynb.
!> @return None.
!>
!> @param n [in] Argument n.
!> @param z [in] Argument z.
!> @param nm [inout] Argument nm.
!> @param cbj [inout] Argument cbj.
!> @param cdj [inout] Argument cdj.
!> @param cby [inout] Argument cby.
!> @param cdy [inout] Argument cdy.
subroutine cjynb ( n, z, nm, cbj, cdj, cby, cdy )

!*****************************************************************************80
!
!! CJYNB: Bessel functions, derivatives, Jn(z) and Yn(z) of complex argument.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    03 August 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
! 
!  Parameters:
!
!    Input, integer(int32) N, the order of Jn(z) and Yn(z).
!
!    Input, complex ( kind = real64 ) Z, the argument of Jn(z) and Yn(z).
!
!    Output, integer(int32) NM, the highest order computed.
!
!    Output, complex ( kind = real64 ) CBJ(0:N), CDJ(0:N), CBY(0:N), CDY(0:N), 
!    the values of Jn(z), Jn'(z), Yn(z), Yn'(z).
!
  implicit none

  integer(int32), intent(in) :: n

  real(real64), parameter, dimension ( 4 ) :: a = [&
    -0.7031250000000000e-01_real64, 0.1121520996093750e+00_real64, &
    -0.5725014209747314e+00_real64, 0.6074042001273483e+01_real64]
  real(real64) a0
  real(real64), parameter, dimension ( 4 ) :: a1 = [&
    0.1171875000000000e+00_real64,-0.1441955566406250e+00_real64, &
    0.6765925884246826e+00_real64,-0.6883914268109947e+01_real64]
  real(real64), parameter, dimension ( 4 ) :: b = [&
    0.7324218750000000e-01_real64,-0.2271080017089844e+00_real64, &
    0.1727727502584457e+01_real64,-0.2438052969955606e+02_real64]
  real(real64), parameter, dimension ( 4 ) :: b1 = [&
   -0.1025390625000000e+00_real64,0.2775764465332031e+00_real64, &
   -0.1993531733751297e+01_real64,0.2724882731126854e+02_real64]
  complex(real64), intent(inout) :: cbj(0:n)
  complex(real64) cbj0
  complex(real64) cbj1
  complex(real64) cbjk
  complex(real64) cbs
  complex(real64), intent(inout) :: cby(0:n)
  complex(real64) cby0
  complex(real64) cby1
  complex(real64), intent(inout) :: cdj(0:n)
  complex(real64), intent(inout) :: cdy(0:n)
  complex(real64) ce
  complex(real64) cf
  complex(real64) cf1
  complex(real64) cf2
  complex(real64) cp0
  complex(real64) cp1
  complex(real64) cq0
  complex(real64) cq1
  complex(real64) cs0
  complex(real64) csu
  complex(real64) csv
  complex(real64) ct1
  complex(real64) ct2
  complex(real64) cu
  complex(real64) cyy
  real(real64) el
  integer(int32) k
  integer(int32) m
  ! integer(int32) msta1
  ! integer(int32) msta2
  integer(int32), intent(inout) :: nm
  real(real64) pi
  real(real64) r2p
  real(real64) y0
  complex(real64), intent(in) :: z

  el = 0.5772156649015329e+00_real64
  pi = 3.141592653589793e+00_real64
  r2p = 0.63661977236758e+00_real64
  y0 = abs ( imag ( z ) )
  a0 = abs ( z )
  nm = n

  if ( a0 < 1.0e-100_real64 ) then
    do k = 0, n
      cbj(k) = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
      cdj(k) = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
      cby(k) = - cmplx ( 1.0e+30_real64, 0.0e+00_real64, kind = real64 )
      cdy(k) = cmplx ( 1.0e+30_real64, 0.0e+00_real64, kind = real64 )
    end do
    cbj(0) = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    cdj(1) = cmplx ( 0.5e+00_real64, 0.0e+00_real64, kind = real64 )
    return
  end if

  if ( a0 <= 300.0e+00_real64 .or. 80 < n ) then

    if ( n == 0 ) then
      nm = 1
    end if
    m = msta1 ( a0, 200 )
    if ( m < nm ) then
      nm = m
    else
      m = msta2 ( a0, nm, 15 )
    end if

    cbs = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
    csu = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
    csv = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
    cf2 = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
    cf1 = cmplx ( 1.0e-30_real64, 0.0e+00_real64, kind = real64 )

    do k = m, 0, -1
      cf = 2.0e+00_real64 * ( k + 1.0e+00_real64 ) / z * cf1 - cf2
      if ( k <= nm ) then
        cbj(k) = cf
      end if
      if ( k == 2 * int ( k / 2 ) .and. k .ne. 0 ) then
        if ( y0 <= 1.0e+00_real64 ) then
          cbs = cbs + 2.0e+00_real64 * cf
        else
          cbs = cbs + ( -1.0e+00_real64 ) ** ( k / 2 ) * 2.0e+00_real64 * cf
        end if
        csu = csu + ( -1.0e+00_real64 ) ** ( k / 2 ) * cf / k
      else if ( 1 < k ) then
        csv = csv + ( -1.0e+00_real64 ) ** ( k / 2 ) * k / ( k * k - 1.0e+00_real64 ) * cf
      end if
      cf2 = cf1
      cf1 = cf
    end do

    if ( y0 <= 1.0e+00_real64 ) then
      cs0 = cbs + cf
    else
      cs0 = ( cbs + cf ) / cos ( z )
    end if

    do k = 0, nm
      cbj(k) = cbj(k) / cs0
    end do

    ce = log ( z / 2.0e+00_real64 ) + el
    cby(0) = r2p * ( ce * cbj(0) - 4.0e+00_real64 * csu / cs0 )
    cby(1) = r2p * ( - cbj(0) / z + ( ce - 1.0e+00_real64 ) * cbj(1) &
      - 4.0e+00_real64 * csv / cs0 )

  else

    ct1 = z - 0.25e+00_real64 * pi
    cp0 = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    do k = 1, 4
      cp0 = cp0 + a(k) * z ** ( - 2 * k )
    end do
    cq0 = -0.125e+00_real64 / z
    do k = 1, 4
      cq0 = cq0 + b(k) * z ** ( - 2 * k - 1 )
    end do
    cu = sqrt ( r2p / z )
    cbj0 = cu * ( cp0 * cos ( ct1 ) - cq0 * sin ( ct1 ) )
    cby0 = cu * ( cp0 * sin ( ct1 ) + cq0 * cos ( ct1 ) )
    cbj(0) = cbj0
    cby(0) = cby0
    ct2 = z - 0.75e+00_real64 * pi
    cp1 = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    do k = 1, 4
      cp1 = cp1 + a1(k) * z ** ( - 2 * k )
    end do
    cq1 = 0.375e+00_real64 / z
    do k = 1, 4
      cq1 = cq1 + b1(k) * z ** ( - 2 * k - 1 )
    end do
    cbj1 = cu * ( cp1 * cos ( ct2 ) - cq1 * sin ( ct2 ) )
    cby1 = cu * ( cp1 * sin ( ct2 ) + cq1 * cos ( ct2 ) )
    cbj(1) = cbj1
    cby(1) = cby1
    do k = 2, nm
      cbjk = 2.0e+00_real64 * ( k - 1.0e+00_real64 ) / z * cbj1 - cbj0
      cbj(k) = cbjk
      cbj0 = cbj1
      cbj1 = cbjk
    end do
  end if

  cdj(0) = -cbj(1)
  do k = 1, nm
    cdj(k) = cbj(k-1) - k / z * cbj(k)
  end do

  if ( 1.0e+00_real64 < abs ( cbj(0) ) ) then
    cby(1) = ( cbj(1) * cby(0) - 2.0e+00_real64 / ( pi * z ) ) / cbj(0)
  end if

  do k = 2, nm
    if ( abs ( cbj(k-2) ) <= abs ( cbj(k-1) ) ) then
      cyy = ( cbj(k) * cby(k-1) - 2.0e+00_real64 / ( pi * z ) ) / cbj(k-1)
    else
      cyy = ( cbj(k) * cby(k-2) - 4.0e+00_real64 * ( k - 1.0e+00_real64 ) &
        / ( pi * z * z ) ) / cbj(k-2)
    end if
    cby(k) = cyy
  end do

  cdy(0) = -cby(1)
  do k = 1, nm
    cdy(k) = cby(k-1) - k / z * cby(k)
  end do

  return
end subroutine cjynb
!> @brief subroutine cjyva.
!> @return None.
!>
!> @param v [in] Argument v.
!> @param z [in] Argument z.
!> @param vm [inout] Argument vm.
!> @param cbj [inout] Argument cbj.
!> @param cdj [inout] Argument cdj.
!> @param cby [inout] Argument cby.
!> @param cdy [inout] Argument cdy.
subroutine cjyva ( v, z, vm, cbj, cdj, cby, cdy )

!*****************************************************************************80
!
!! CJYVA: Bessel functions and derivatives, Jv(z) and Yv(z) of complex argument.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    03 August 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) V, the order of Jv(z) and Yv(z).
!
!    Input, complex ( kind = real64 ), the argument.
!
!    Output, real(real64) VM, the highest order computed.
!
!    Output, real(real64) CBJ(0:*), CDJ(0:*), CBY(0:*), CDY(0:*), 
!    the values of Jn+v0(z), Jn+v0'(z), Yn+v0(z), Yn+v0'(z).
!
  implicit none

  real(real64) a0
  complex(real64) ca
  complex(real64) ca0
  complex(real64) cb
  complex(real64), intent(inout) :: cbj(0:)
  complex(real64), intent(inout) :: cby(0:)
  complex(real64) cck
  complex(real64), intent(inout) :: cdj(0:)
  complex(real64), intent(inout) :: cdy(0:)
  complex(real64) cec
  complex(real64) cf
  complex(real64) cf0
  complex(real64) cf1
  complex(real64) cf2
  complex(real64) cfac0
  complex(real64) cfac1
  complex(real64) cg0
  complex(real64) cg1
  complex(real64) ch0
  complex(real64) ch1
  complex(real64) ch2
  complex(real64) ci
  complex(real64) cju0
  complex(real64) cju1
  complex(real64) cjv0
  complex(real64) cjv1
  complex(real64) cjvl
  complex(real64) cp11
  complex(real64) cp12
  complex(real64) cp21
  complex(real64) cp22
  complex(real64) cpz
  complex(real64) cqz
  complex(real64) cr
  complex(real64) cr0
  complex(real64) cr1
  complex(real64) crp
  complex(real64) crq
  complex(real64) cs
  complex(real64) cs0
  complex(real64) cs1
  complex(real64) csk
  complex(real64) cyk
  complex(real64) cyl1
  complex(real64) cyl2
  complex(real64) cylk
  complex(real64) cyv0
  complex(real64) cyv1
  real(real64) ga
  real(real64) gb
  integer(int32) j
  integer(int32) k
  integer(int32) k0
  integer(int32) l
  integer(int32) lb
  integer(int32) lb0
  integer(int32) m
  ! integer(int32) msta1
  ! integer(int32) msta2
  integer(int32) n
  real(real64) pi
  real(real64) pv0
  real(real64) pv1
  real(real64) rp2
  real(real64), intent(in) :: v
  real(real64) v0
  real(real64) vg
  real(real64) vl
  real(real64), intent(inout) :: vm
  real(real64) vv
  real(real64) w0
  real(real64) w1
  real(real64) wa
  real(real64) ya0
  real(real64) ya1
  real(real64) yak
  complex(real64), intent(in) :: z
  complex(real64) z1
  complex(real64) z2
  complex(real64) zk

  pi = 3.141592653589793e+00_real64
  rp2 = 0.63661977236758e+00_real64
  ci = cmplx ( 0.0e+00_real64, 1.0e+00_real64, kind = real64 )
  a0 = abs ( z )
  z1 = z
  z2 = z * z
  n = int ( v )
  v0 = v - n
  pv0 = pi * v0
  pv1 = pi * ( 1.0e+00_real64 + v0 )

  if ( a0 < 1.0e-100_real64 ) then

    do k = 0, n
      cbj(k) = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
      cdj(k) = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
      cby(k) = - cmplx ( 1.0e+30_real64, 0.0e+00_real64, kind = real64 )
      cdy(k) = cmplx ( 1.0e+30_real64, 0.0e+00_real64, kind = real64 )
    end do

    if ( v0 == 0.0e+00_real64 ) then
      cbj(0) = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
      cdj(1) = cmplx ( 0.5e+00_real64, 0.0e+00_real64, kind = real64 )
    else
      cdj(0) = cmplx ( 1.0e+30_real64, 0.0e+00_real64, kind = real64 )
    end if

    vm = v                     
    return

  end if

  if ( real ( z, kind = real64 ) < 0.0e+00_real64 ) then
    z1 = -z
  end if

  if ( a0 <= 12.0e+00_real64 ) then

    do l = 0, 1
      vl = v0 + l
      cjvl = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
      cr = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
      do k = 1, 40
        cr = -0.25e+00_real64 * cr * z2 / ( k * ( k + vl ) )
        cjvl = cjvl + cr
        if ( abs ( cr ) < abs ( cjvl ) * 1.0e-15_real64 ) then
          exit
        end if
      end do

      vg = 1.0e+00_real64 + vl
      call gamma ( vg, ga )
      ca = ( 0.5e+00_real64 * z1 ) ** vl / ga

      if ( l == 0 ) then
        cjv0 = cjvl * ca
      else
        cjv1 = cjvl * ca
      end if

    end do

  else

    if ( a0 < 35.0e+00_real64 ) then
      k0 = 11
    else if ( a0 <50.0e+00_real64 ) then
      k0 = 10
    else
      k0 = 8
    end if

    do j = 0, 1
      vv = 4.0e+00_real64 * ( j + v0 ) * ( j + v0 )
      cpz = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
      crp = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
      do k = 1, k0
        crp = - 0.78125e-02_real64 * crp &
          * ( vv - ( 4.0e+00_real64 * k - 3.0e+00_real64 ) ** 2 ) &
          * ( vv - ( 4.0e+00_real64 * k - 1.0e+00_real64 ) ** 2 )  &
          / ( k * ( 2.0e+00_real64 * k - 1.0e+00_real64 ) * z2 )
        cpz = cpz + crp
      end do
      cqz = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
      crq = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
      do k = 1, k0
        crq = -0.78125e-02_real64 * crq &
          * ( vv - ( 4.0e+00_real64 * k - 1.0e+00_real64 ) ** 2 ) &
          * ( vv - ( 4.0e+00_real64 * k + 1.0e+00_real64 ) ** 2 ) &
          / ( k * ( 2.0e+00_real64 * k + 1.0e+00_real64 ) * z2 )
        cqz = cqz + crq
      end do
      cqz = 0.125e+00_real64 * ( vv - 1.0e+00_real64 ) * cqz / z1
      zk = z1 - ( 0.5e+00_real64 * ( j + v0 ) + 0.25e+00_real64 ) * pi
      ca0 = sqrt ( rp2 / z1 )
      cck = cos ( zk )
      csk = sin ( zk )
      if ( j == 0 ) then
        cjv0 = ca0 * ( cpz * cck - cqz * csk )
        cyv0 = ca0 * ( cpz * csk + cqz * cck )
      else if ( j == 1 ) then
        cjv1 = ca0 * ( cpz * cck - cqz * csk )
        cyv1 = ca0 * ( cpz * csk + cqz * cck )
      end if
    end do

  end if

  if ( a0 <= 12.0e+00_real64 ) then

    if ( v0 .ne. 0.0e+00_real64 ) then

      do l = 0, 1
        vl = v0 + l
        cjvl = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
        cr = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
        do k = 1, 40
          cr = -0.25e+00_real64 * cr * z2 / ( k * ( k - vl ) )
          cjvl = cjvl + cr
          if ( abs ( cr ) < abs ( cjvl ) * 1.0e-15_real64 ) then
            exit
          end if
        end do

        vg = 1.0e+00_real64 - vl
        call gamma ( vg, gb )
        cb = ( 2.0e+00_real64 / z1 ) ** vl / gb
        if ( l == 0 ) then
          cju0 = cjvl * cb
        else
          cju1 = cjvl * cb
        end if
      end do
      cyv0 = ( cjv0 * cos ( pv0 ) - cju0 ) / sin ( pv0 )
      cyv1 = ( cjv1 * cos ( pv1 ) - cju1 ) / sin ( pv1 )

    else

      cec = log ( z1 / 2.0e+00_real64 ) + 0.5772156649015329e+00_real64
      cs0 = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
      w0 = 0.0e+00_real64
      cr0 = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
      do k = 1, 30
        w0 = w0 + 1.0e+00_real64 / k
        cr0 = -0.25e+00_real64 * cr0 / ( k * k ) * z2
        cs0 = cs0 + cr0 * w0
      end do
      cyv0 = rp2 * ( cec * cjv0 - cs0 )
      cs1 = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
      w1 = 0.0e+00_real64
      cr1 = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
      do k = 1, 30
        w1 = w1 + 1.0e+00_real64 / k
        cr1 = -0.25e+00_real64 * cr1 / ( k * ( k + 1 ) ) * z2
        cs1 = cs1 + cr1 * ( 2.0e+00_real64 * w1 + 1.0e+00_real64 / ( k + 1.0e+00_real64 ) )
      end do
      cyv1 = rp2 * ( cec * cjv1 - 1.0e+00_real64 / z1 - 0.25e+00_real64 * z1 * cs1 )

    end if

  end if

  if ( real ( z, kind = real64 ) < 0.0e+00_real64 ) then

    cfac0 = exp ( pv0 * ci )
    cfac1 = exp ( pv1 * ci )

    if ( imag ( z ) < 0.0e+00_real64 ) then
      cyv0 = cfac0 * cyv0 - 2.0e+00_real64 * ci * cos ( pv0 ) * cjv0
      cyv1 = cfac1 * cyv1 - 2.0e+00_real64 * ci * cos ( pv1 ) * cjv1
      cjv0 = cjv0 / cfac0
      cjv1 = cjv1 / cfac1
    else if ( 0.0e+00_real64 < imag ( z ) ) then
      cyv0 = cyv0 / cfac0 + 2.0e+00_real64 * ci * cos ( pv0 ) * cjv0
      cyv1 = cyv1 / cfac1 + 2.0e+00_real64 * ci * cos ( pv1 ) * cjv1
      cjv0 = cfac0 * cjv0
      cjv1 = cfac1 * cjv1
    end if

  end if

  cbj(0) = cjv0
  cbj(1) = cjv1

  if ( 2 <= n .and. n <= int ( 0.25e+00_real64 * a0 ) ) then

    cf0 = cjv0
    cf1 = cjv1
    do k = 2, n
      cf = 2.0e+00_real64 * ( k + v0 - 1.0e+00_real64 ) / z * cf1 - cf0
      cbj(k) = cf
      cf0 = cf1
      cf1 = cf
    end do

  else if ( 2 <= n ) then

    m = msta1 ( a0, 200 )
    if ( m < n ) then
      n = m
    else
      m = msta2 ( a0, n, 15 )
    end if
    cf2 = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
    cf1 = cmplx ( 1.0e-30_real64, 0.0e+00_real64, kind = real64 )
    do k = m, 0, -1
      cf = 2.0e+00_real64 * ( v0 + k + 1.0e+00_real64 ) / z * cf1 - cf2
      if ( k <= n ) then
        cbj(k) = cf
      end if
      cf2 = cf1
      cf1 = cf
    end do
    if ( abs ( cjv1 ) < abs ( cjv0 ) ) then
      cs = cjv0 / cf
    else
      cs = cjv1 / cf2
    end if

    do k = 0, n
      cbj(k) = cs * cbj(k)
    end do

  end if

    cdj(0) = v0 / z * cbj(0) - cbj(1)
    do k = 1, n
      cdj(k) = - ( k + v0 ) / z * cbj(k) + cbj(k-1)
    end do

    cby(0) = cyv0
    cby(1) = cyv1
    ya0 = abs ( cyv0 )
    lb = 0
    cg0 = cyv0
    cg1 = cyv1
    do k = 2, n
      cyk = 2.0e+00_real64 * ( v0 + k - 1.0e+00_real64 ) / z * cg1 - cg0
      if ( abs ( cyk ) <= 1.0e+290_real64 ) then
        yak = abs ( cyk )
        ya1 = abs ( cg0 )
        if ( yak < ya0 .and. yak < ya1 ) then
          lb = k
        end if
        cby(k) = cyk
        cg0 = cg1
        cg1 = cyk
      end if
    end do

    if ( 4 < lb .and. imag ( z ) /= 0.0e+00_real64 ) then

      do

        if ( lb == lb0 ) then
          exit
        end if

        ch2 = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
        ch1 = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
        lb0 = lb
        do k = lb, 1, -1
          ch0 = 2.0e+00_real64 * ( k + v0 ) / z * ch1 - ch2
          ch2 = ch1
          ch1 = ch0
        end do
        cp12 = ch0
        cp22 = ch2
        ch2 = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
        ch1 = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
        do k = lb, 1, -1
          ch0 = 2.0e+00_real64 * ( k + v0 ) / z * ch1 - ch2
          ch2 = ch1
          ch1 = ch0
        end do
        cp11 = ch0
        cp21 = ch2

        if ( lb == n ) then
          cbj(lb+1) = 2.0e+00_real64 * ( lb + v0 ) / z * cbj(lb) - cbj(lb-1)
        end if

        if ( abs ( cbj(1) ) < abs ( cbj(0) ) ) then
          cby(lb+1) = ( cbj(lb+1) * cyv0 - 2.0e+00_real64 * cp11 / ( pi * z ) ) &
            / cbj(0)
          cby(lb) = ( cbj(lb) * cyv0 + 2.0e+00_real64 * cp12 / ( pi * z ) ) / cbj(0)
        else
          cby(lb+1) = ( cbj(lb+1) * cyv1 - 2.0e+00_real64 * cp21 / ( pi * z ) ) &
            / cbj(1)
          cby(lb) = ( cbj(lb) * cyv1 + 2.0e+00_real64 * cp22 / ( pi * z ) ) / cbj(1)
        end if

        cyl2 = cby(lb+1)
        cyl1 = cby(lb)
        do k = lb - 1, 0, -1
          cylk = 2.0e+00_real64 * ( k + v0 + 1.0e+00_real64 ) / z * cyl1 - cyl2
          cby(k) = cylk
          cyl2 = cyl1
          cyl1 = cylk
        end do

      cyl1 = cby(lb)
      cyl2 = cby(lb+1)
      do k = lb + 1, n - 1
        cylk = 2.0e+00_real64 * ( k + v0 ) / z * cyl2 - cyl1
        cby(k+1) = cylk
        cyl1 = cyl2
        cyl2 = cylk
      end do

      do k = 2, n
        wa = abs ( cby(k) )
        if ( wa < abs ( cby(k-1) ) ) then
          lb = k
        end if
      end do

    end do

  end if

  cdy(0) = v0 / z * cby(0) - cby(1)
  do k = 1, n
    cdy(k) = cby(k-1) - ( k + v0 ) / z * cby(k)
  end do
  vm = n + v0

  return
end subroutine cjyva
!> @brief subroutine cjyvb.
!> @return None.
!>
!> @param v [in] Argument v.
!> @param z [in] Argument z.
!> @param vm [inout] Argument vm.
!> @param cbj [inout] Argument cbj.
!> @param cdj [inout] Argument cdj.
!> @param cby [inout] Argument cby.
!> @param cdy [inout] Argument cdy.
subroutine cjyvb ( v, z, vm, cbj, cdj, cby, cdy )

!*****************************************************************************80
!
!! CJYVB: Bessel functions and derivatives, Jv(z) and Yv(z) of complex argument.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    03 August 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) V, the order of Jv(z) and Yv(z).
!
!    Input, complex ( kind = real64 ) Z, the argument.
!
!    Output, real(real64) VM, the highest order computed.
!
!    Output, real(real64) CBJ(0:*), CDJ(0:*), CBY(0:*), CDY(0:*), 
!    the values of Jn+v0(z), Jn+v0'(z), Yn+v0(z), Yn+v0'(z).
!
  implicit none

  real(real64) a0
  complex(real64) ca
  complex(real64) ca0
  complex(real64) cb
  complex(real64), intent(inout) :: cbj(0:)
  complex(real64), intent(inout) :: cby(0:)
  complex(real64) cck
  complex(real64), intent(inout) :: cdj(0:)
  complex(real64), intent(inout) :: cdy(0:)
  complex(real64) cec
  complex(real64) cf
  complex(real64) cf1
  complex(real64) cf2
  complex(real64) cfac0
  complex(real64) ci
  complex(real64) cju0
  complex(real64) cjv0
  complex(real64) cjvn
  complex(real64) cpz
  complex(real64) cqz
  complex(real64) cr
  complex(real64) cr0
  complex(real64) crp
  complex(real64) crq
  complex(real64) cs
  complex(real64) cs0
  complex(real64) csk
  complex(real64) cyv0
  complex(real64) cyy
  real(real64) ga
  real(real64) gb
  integer(int32) k
  integer(int32) k0
  integer(int32) m
  ! integer(int32) msta1
  ! integer(int32) msta2
  integer(int32) n
  real(real64) pi
  real(real64) pv0
  real(real64) rp2
  real(real64), intent(in) :: v
  real(real64) v0
  real(real64) vg
  real(real64), intent(inout) :: vm
  real(real64) vv
  real(real64) w0
  complex(real64), intent(in) :: z
  complex(real64) z1
  complex(real64) z2
  complex(real64) zk

  pi = 3.141592653589793e+00_real64
  rp2 = 0.63661977236758e+00_real64
  ci = cmplx ( 0.0e+00_real64, 1.0e+00_real64, kind = real64 )
  a0 = abs ( z )
  z1 = z
  z2 = z * z
  n = int ( v )
  v0 = v - n
  pv0 = pi * v0
  
  if ( a0 < 1.0e-100_real64 ) then

    do k = 0, n
      cbj(k) = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
      cdj(k) = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
      cby(k) = - cmplx ( 1.0e+30_real64, 0.0e+00_real64, kind = real64 )
      cdy(k) = cmplx ( 1.0e+30_real64, 0.0e+00_real64, kind = real64 )
    end do

    if ( v0 == 0.0e+00_real64 ) then
      cbj(0) = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
      cdj(1) = cmplx ( 0.5e+00_real64, 0.0e+00_real64, kind = real64 )
    else
      cdj(0) = cmplx ( 1.0e+30_real64, 0.0e+00_real64, kind = real64 )
    end if

    vm = v
    return

  end if

  if ( real ( z, kind = real64 ) < 0.0e+00_real64 ) then
    z1 = -z
  end if

  if ( a0 <= 12.0e+00_real64 ) then

    cjv0 = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    cr = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    do k = 1, 40
      cr = -0.25e+00_real64 * cr * z2 / ( k * ( k + v0 ) )
      cjv0 = cjv0 + cr
      if ( abs ( cr ) < abs ( cjv0 ) * 1.0e-15_real64 ) then
        exit
      end if
    end do

    vg = 1.0e+00_real64 + v0
    call gamma ( vg, ga )
    ca = ( 0.5e+00_real64 * z1 ) ** v0 / ga
    cjv0 = cjv0 * ca

  else

    if ( a0 < 35.0e+00_real64 ) then
      k0 = 11
    else if ( a0 < 50.0e+00_real64 ) then
      k0 = 10
    else
      k0 = 8
    end if

    vv = 4.0e+00_real64 * v0 * v0
    cpz = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    crp = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    do k = 1, k0
      crp = -0.78125e-02_real64 * crp &
        * ( vv - ( 4.0e+00_real64 * k - 3.0e+00_real64 ) ** 2 ) &
        * ( vv - ( 4.0e+00_real64 * k - 1.0e+00_real64 ) **2 ) &
        / ( k * ( 2.0e+00_real64 * k - 1.0e+00_real64 ) * z2 )
      cpz = cpz + crp
    end do
    cqz = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    crq = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    do k = 1, k0
      crq = -0.78125e-02_real64 * crq &
        * ( vv - ( 4.0e+00_real64 * k - 1.0e+00_real64 ) ** 2 ) &
        * ( vv - ( 4.0e+00_real64 * k + 1.0e+00_real64 ) ** 2 ) &
        / ( k * ( 2.0e+00_real64 * k + 1.0e+00_real64 ) * z2 )
      cqz = cqz + crq
    end do
    cqz = 0.125e+00_real64 * ( vv - 1.0e+00_real64 ) * cqz / z1
    zk = z1 - ( 0.5e+00_real64 * v0 + 0.25e+00_real64 ) * pi
    ca0 = sqrt ( rp2 / z1 )
    cck = cos ( zk )
    csk = sin ( zk )
    cjv0 = ca0 * ( cpz * cck - cqz * csk )
    cyv0 = ca0 * ( cpz * csk + cqz * cck )

  end if

  if ( a0 <= 12.0e+00_real64 ) then

    if ( v0 .ne. 0.0e+00_real64 ) then

      cjvn = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
      cr = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
      do k = 1, 40
        cr = -0.25e+00_real64 * cr * z2 / ( k * ( k - v0 ) )
        cjvn = cjvn + cr
        if ( abs ( cr ) < abs ( cjvn ) * 1.0e-15_real64 ) then
          exit
        end if
      end do

      vg = 1.0e+00_real64 - v0
      call gamma ( vg, gb )
      cb = ( 2.0e+00_real64 / z1 ) ** v0 / gb
      cju0 = cjvn * cb
      cyv0 = ( cjv0 * cos ( pv0 ) - cju0 ) / sin ( pv0 )

    else

      cec = log ( z1 / 2.0e+00_real64 ) + 0.5772156649015329e+00_real64
      cs0 = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
      w0 = 0.0e+00_real64
      cr0 = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
      do k = 1, 30
        w0 = w0 + 1.0e+00_real64 / k
        cr0 = -0.25e+00_real64 * cr0 / ( k * k ) * z2
        cs0 = cs0 + cr0 * w0
      end do
      cyv0 = rp2 * ( cec * cjv0 - cs0 )

    end if

  end if

  if ( n == 0 ) then
    n = 1
  end if

  m = msta1 ( a0, 200 )
  if ( m < n ) then
    n = m
  else
    m = msta2 ( a0, n, 15 )
  end if

  cf2 = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
  cf1 = cmplx ( 1.0e-30_real64, 0.0e+00_real64, kind = real64 )
  do k = m, 0, -1
    cf = 2.0e+00_real64 * ( v0 + k + 1.0e+00_real64 ) / z1 * cf1 - cf2
    if ( k <= n ) then
      cbj(k) = cf
    end if
    cf2 = cf1
    cf1 = cf
  end do

  cs = cjv0 / cf
  do k = 0, n
    cbj(k) = cs * cbj(k)
  end do

  if ( real ( z, kind = real64 ) < 0.0e+00_real64) then

    cfac0 = exp ( pv0 * ci )
    if ( imag ( z ) < 0.0e+00_real64 ) then
      cyv0 = cfac0 * cyv0 - 2.0e+00_real64 * ci * cos ( pv0 ) * cjv0
    else if ( 0.0e+00_real64 < imag ( z ) ) then
      cyv0 = cyv0 / cfac0 + 2.0e+00_real64 * ci * cos ( pv0 ) * cjv0
    end if

    do k = 0, n
      if ( imag ( z ) < 0.0e+00_real64) then
        cbj(k) = exp ( - pi * ( k + v0 ) * ci ) * cbj(k)
      else if ( 0.0e+00_real64 < imag ( z ) ) then
        cbj(k) = exp ( pi * ( k + v0 ) * ci ) * cbj(k)
      end if
    end do

    z1 = z1

  end if

  cby(0) = cyv0
  do k = 1, n
    cyy = ( cbj(k) * cby(k-1) - 2.0e+00_real64 / ( pi * z ) ) / cbj(k-1)
    cby(k) = cyy
  end do

  cdj(0) = v0 / z * cbj(0) - cbj(1)
  do k = 1, n
    cdj(k) = - ( k + v0 ) / z * cbj(k) + cbj(k-1)
  end do

  cdy(0) = v0 / z * cby(0) - cby(1)
  do k = 1, n
    cdy(k) = cby(k-1) - ( k + v0 ) / z * cby(k)
  end do

  vm = n + v0

  return
end subroutine cjyvb
!> @brief subroutine clpmn.
!> @return None.
!>
!> @param mm [in] Argument mm.
!> @param m [in] Argument m.
!> @param n [inout] Argument n.
!> @param x [in] Argument x.
!> @param y [in] Argument y.
!> @param cpm [inout] Argument cpm.
!> @param cpd [inout] Argument cpd.
pure subroutine clpmn ( mm, m, n, x, y, cpm, cpd )

!*****************************************************************************80
!
!! CLPMN: associated Legendre functions and derivatives for complex argument.
!
!  Discussion:
!
!    Compute the associated Legendre functions Pmn(z)   
!    and their derivatives Pmn'(z) for a complex argument
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    01 August 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) MM, the physical dimension of CPM and CPD.
!
!    Input, integer(int32) M, N, the order and degree of Pmn(z).
!
!    Input, real(real64) X, Y, the real and imaginary parts of 
!    the argument Z.
!
!    Output, complex ( kind = real64 ) CPM(0:MM,0:N), CPD(0:MM,0:N), the values of
!    Pmn(z) and Pmn'(z).
!
  implicit none

  integer(int32), intent(in) :: mm

  complex(real64), intent(inout) :: cpd(0:mm,0:n)
  complex(real64), intent(inout) :: cpm(0:mm,0:n)
  integer(int32) i
  integer(int32) j
  integer(int32) ls
  integer(int32), intent(in) :: m
  integer(int32), intent(in) :: n
  real(real64), intent(in) :: x
  real(real64), intent(in) :: y
  complex(real64) z
  complex(real64) zq
  complex(real64) zs

  z = cmplx ( x, y, kind = real64 )

  do i = 0, n
    do j = 0, m
      cpm(j,i) = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
      cpd(j,i) = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
    end do
  end do

  cpm(0,0) = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )

  if ( abs ( x ) == 1.0e+00_real64 .and. y == 0.0e+00_real64 ) then

    do i = 1, n
      cpm(0,i) = x ** i
      cpd(0,i) = 0.5e+00_real64 * i * ( i + 1 ) * x ** ( i + 1 )
    end do

    do j = 1, n
      do i = 1, m
        if ( i == 1 ) then
          cpd(i,j) = cmplx ( 1.0e+30_real64, 0.0e+00_real64, kind = real64 )
        else if ( i == 2 ) then
          cpd(i,j) = -0.25e+00_real64 &
            * ( j + 2 ) * ( j + 1 ) * j * ( j - 1 ) * x ** ( j + 1 )
        end if
      end do
    end do

    return

  end if

  if ( 1.0e+00_real64 < abs ( z ) ) then
    ls = -1
  else
    ls = 1
  end if

  zq = sqrt ( ls * ( 1.0e+00_real64 - z * z ) )
  zs = ls * ( 1.0e+00_real64 - z * z )
  do i = 1, m
    cpm(i,i) = -ls * ( 2.0e+00_real64 * i - 1.0e+00_real64 ) * zq * cpm(i-1,i-1)
  end do
  do i = 0, m
    cpm(i,i+1) = ( 2.0e+00_real64 * i + 1.0e+00_real64 ) * z * cpm(i,i)
  end do

  do i = 0, m
    do j = i + 2, n
      cpm(i,j) = ( ( 2.0e+00_real64 * j - 1.0e+00_real64 ) * z * cpm(i,j-1) &
        - ( i + j - 1.0e+00_real64 ) * cpm(i,j-2) ) / ( j - i )
    end do
  end do

  cpd(0,0) = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
  do j = 1, n
    cpd(0,j) = ls * j * ( cpm(0,j-1) - z * cpm(0,j) ) / zs
  end do 

  do i = 1, m
    do j = i, n
      cpd(i,j) = ls * i * z * cpm(i,j) / zs &
        + ( j + i ) * ( j - i + 1.0e+00_real64 ) / zq * cpm(i-1,j)
    end do
  end do

  return
end subroutine clpmn
!> @brief subroutine clpn.
!> @return None.
!>
!> @param n [in] Argument n.
!> @param x [in] Argument x.
!> @param y [in] Argument y.
!> @param cpn [inout] Argument cpn.
!> @param cpd [inout] Argument cpd.
pure subroutine clpn ( n, x, y, cpn, cpd )

!*****************************************************************************80
!
!! CLPN computes Legendre functions and derivatives for complex argument.
!
!  Discussion:
!
!    Compute Legendre polynomials Pn(z) and their derivatives Pn'(z) for 
!    a complex argument.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    15 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) N, the degree.
!
!    Input, real(real64) X, Y, the real and imaginary parts 
!    of the argument.
!
!    Output, complex ( kind = real64 ) CPN(0:N), CPD(0:N), the values of Pn(z)
!    and Pn'(z).
!
  implicit none

  integer(int32), intent(in) :: n

  complex(real64) cp0
  complex(real64) cp1
  complex(real64), intent(inout) :: cpd(0:n)
  complex(real64) cpf
  complex(real64), intent(inout) :: cpn(0:n)
  integer(int32) k
  real(real64), intent(in) :: x
  real(real64), intent(in) :: y
  complex(real64) z

  z = cmplx ( x, y, kind = real64 )

  cpn(0) = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
  cpn(1) = z
  cpd(0) = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
  cpd(1) = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )

  cp0 = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
  cp1 = z
  do k = 2, n
    cpf = ( 2.0e+00_real64 * k - 1.0e+00_real64 ) / k * z * cp1 - ( k - 1.0e+00_real64 ) / k * cp0
    cpn(k) = cpf
    if ( abs ( x ) == 1.0e+00_real64 .and. y == 0.0e+00_real64 ) then
      cpd(k) = 0.5e+00_real64 * x ** ( k + 1 ) * k * ( k + 1.0e+00_real64 )
    else
      cpd(k) = k * ( cp1 - z * cpf ) / ( 1.0e+00_real64 - z * z )
    end if
    cp0 = cp1
    cp1 = cpf
  end do

  return
end subroutine clpn
!> @brief subroutine clqmn.
!> @return None.
!>
!> @param mm [in] Argument mm.
!> @param m [in] Argument m.
!> @param n [in] Argument n.
!> @param x [in] Argument x.
!> @param y [in] Argument y.
!> @param cqm [inout] Argument cqm.
!> @param cqd [inout] Argument cqd.
pure subroutine clqmn ( mm, m, n, x, y, cqm, cqd )

!*****************************************************************************80
!
!! CLQMN: associated Legendre functions and derivatives for complex argument.
!
!  Discussion:
!
!    This procedure computes the associated Legendre functions of the second 
!    kind, Qmn(z) and Qmn'(z), for a complex argument.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    02 August 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
! 
!  Parameters:
!
!    Input, integer(int32) MM, the physical dimension of CQM and CQD.
!
!    Input, integer(int32) M, N, the order and degree of Qmn(z).
!
!    Input, real(real64) X, Y, the real and imaginary parts of the 
!    argument Z.
!
!    Output, complex ( kind = real64 ) CQM(0:MM,0:N), CQD(0:MM,0:N), the values of
!    Qmn(z) and Qmn'(z).
!
  implicit none

  integer(int32), intent(in) :: mm
  integer(int32), intent(in) :: n 

  complex(real64) cq0
  complex(real64) cq1
  complex(real64) cq10
  complex(real64) cqf
  complex(real64) cqf0
  complex(real64) cqf1
  complex(real64) cqf2
  complex(real64), intent(inout) :: cqm(0:mm,0:n)
  complex(real64), intent(inout) :: cqd(0:mm,0:n)
  integer(int32) i
  integer(int32) j
  integer(int32) k
  integer(int32) km
  integer(int32) ls
  integer(int32), intent(in) :: m
  real(real64), intent(in) :: x
  real(real64) xc
  real(real64), intent(in) :: y
  complex(real64) z
  complex(real64) zq
  complex(real64) zs

  z = cmplx ( x, y, kind = real64 )

  if ( abs ( x ) == 1.0e+00_real64 .and. y == 0.0e+00_real64 ) then
    do i = 0, m
      do j = 0, n
        cqm(i,j) = cmplx ( 1.0e+30_real64, 0.0e+00_real64, kind = real64 )
        cqd(i,j) = cmplx ( 1.0e+30_real64, 0.0e+00_real64, kind = real64 )
      end do
    end do
    return
  end if

  xc = abs ( z )

  if ( imag ( z ) == 0.0e+00_real64 .or. xc < 1.0e+00_real64 ) then
    ls = 1
  end if

  if ( 1.0e+00_real64 < xc ) then
    ls = -1
  end if

  zq = sqrt ( ls * ( 1.0e+00_real64 - z * z ) )
  zs = ls * ( 1.0e+00_real64 - z * z )
  cq0 = 0.5e+00_real64 * log ( ls * ( 1.0e+00_real64 + z ) / ( 1.0e+00_real64 - z ) )

  if ( xc < 1.0001e+00_real64 ) then

    cqm(0,0) = cq0
    cqm(0,1) = z * cq0 - 1.0e+00_real64
    cqm(1,0) = -1.0e+00_real64 / zq
    cqm(1,1) = - zq * ( cq0 + z / ( 1.0e+00_real64 - z * z ) )
    do i = 0, 1
      do j = 2, n
        cqm(i,j) = ( ( 2.0e+00_real64 * j - 1.0e+00_real64 ) * z * cqm(i,j-1) &
          - ( j + i - 1.0e+00_real64 ) * cqm(i,j-2) ) / ( j - i )
      end do
    end do

    do j = 0, n
      do i = 2, m
        cqm(i,j) = -2.0e+00_real64 * ( i - 1.0e+00_real64 ) * z / zq * cqm(i-1,j) &
          - ls * ( j + i - 1.0e+00_real64 ) * ( j - i + 2.0e+00_real64 ) * cqm(i-2,j)
      end do
    end do

  else

    if ( 1.1e+00_real64 < xc ) then
      km = 40 + m + n
    else
      km = ( 40 + m + n ) * int ( - 1.0e+00_real64 - 1.8e+00_real64 * log ( xc - 1.0e+00_real64 ) )
    end if

    cqf2 = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
    cqf1 = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    do k = km, 0, -1
      cqf0 = ( ( 2 * k + 3.0e+00_real64 ) * z * cqf1 &
        - ( k + 2.0e+00_real64 ) * cqf2 ) / ( k + 1.0e+00_real64 )
      if ( k <= n ) then
        cqm(0,k) = cqf0
      end if
      cqf2 = cqf1
      cqf1 = cqf0
    end do

    do k = 0, n
      cqm(0,k) = cq0 * cqm(0,k) / cqf0
    end do
 
    cqf2 = 0.0e+00_real64
    cqf1 = 1.0e+00_real64
    do k = km, 0, -1
      cqf0 = ( ( 2 * k + 3.0e+00_real64 ) * z * cqf1 &
        - ( k + 1.0e+00_real64 ) * cqf2 ) / ( k + 2.0e+00_real64 )
      if ( k <= n ) then
        cqm(1,k) = cqf0
      end if
      cqf2 = cqf1
      cqf1 = cqf0
    end do

    cq10 = -1.0e+00_real64 / zq
    do k = 0, n 
      cqm(1,k) = cq10 * cqm(1,k) / cqf0
    end do

    do j = 0, n
      cq0 = cqm(0,j)
      cq1 = cqm(1,j)
      do i = 0, m - 2
        cqf = -2.0e+00_real64 * ( i + 1 ) * z / zq * cq1 &
          + ( j - i ) * ( j + i + 1.0e+00_real64 ) * cq0
        cqm(i+2,j) = cqf
        cq0 = cq1
        cq1 = cqf
      end do
    end do

  end if

  cqd(0,0) = ls / zs
  do j = 1, n
    cqd(0,j) = ls * j * ( cqm(0,j-1) - z * cqm(0,j) ) / zs
  end do

  do j = 0, n
    do i = 1, m
      cqd(i,j) = ls * i * z / zs * cqm(i,j) &
        + ( i + j ) * ( j - i + 1.0e+00_real64 ) / zq * cqm(i-1,j)
    end do
  end do

  return
end subroutine clqmn
!> @brief subroutine clqn.
!> @return None.
!>
!> @param n [in] Argument n.
!> @param x [in] Argument x.
!> @param y [in] Argument y.
!> @param cqn [inout] Argument cqn.
!> @param cqd [inout] Argument cqd.
pure subroutine clqn ( n, x, y, cqn, cqd )

!*****************************************************************************80
!
!! CLQN: Legendre function Qn(z) and derivative Wn'(z) for complex argument.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    01 August 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) N, the degree of Qn(z).
!
!    Input, real(real64) X, Y, the real and imaginary parts of the 
!    argument Z.
!
!    Output, complex ( kind = real64 ) CQN(0:N), CQD(0:N), the values of Qn(z) 
!    and Qn'(z.
!
  implicit none

  integer(int32), intent(in) :: n

  complex(real64) cq0
  complex(real64) cq1
  complex(real64) cqf0
  complex(real64) cqf1
  complex(real64) cqf2
  complex(real64), intent(inout) :: cqn(0:n)
  complex(real64), intent(inout) :: cqd(0:n)
  integer(int32) k
  integer(int32) km
  integer(int32) ls
  real(real64), intent(in) :: x
  real(real64), intent(in) :: y
  complex(real64) z

  z = cmplx ( x, y, kind = real64 )

  if ( z == 1.0e+00_real64 ) then
    do k = 0, n
      cqn(k) = cmplx ( 1.0e+30_real64, 0.0e+00_real64, kind = real64 )
      cqd(k) = cmplx ( 1.0e+30_real64, 0.0e+00_real64, kind = real64 )
    end do
    return
  end if

  if ( 1.0e+00_real64 < abs ( z ) ) then
    ls = -1
  else
    ls = +1
  end if

  cq0 = 0.5e+00_real64 * log ( ls * ( 1.0e+00_real64 + z ) / ( 1.0e+00_real64 - z ) )
  cq1 = z * cq0 - 1.0e+00_real64
  cqn(0) = cq0
  cqn(1) = cq1

  if ( abs ( z ) < 1.0001e+00_real64 ) then

    cqf0 = cq0
    cqf1 = cq1
    do k = 2, n
      cqf2 = ( ( 2.0e+00_real64 * k - 1.0e+00_real64 ) * z * cqf1 &
        - ( k - 1.0e+00_real64 ) * cqf0 ) / k
      cqn(k) = cqf2
      cqf0 = cqf1
      cqf1 = cqf2
    end do

  else

    if ( 1.1e+00_real64 < abs ( z ) ) then
      km = 40 + n
    else
      km = ( 40 + n ) * int ( - 1.0e+00_real64 &
        - 1.8e+00_real64 * log ( abs ( z - 1.0e+00_real64 ) ) )
    end if

    cqf2 = 0.0e+00_real64
    cqf1 = 1.0e+00_real64
    do k = km, 0, -1
      cqf0 = ( ( 2 * k + 3.0e+00_real64 ) * z * cqf1 &
        - ( k + 2.0e+00_real64 ) * cqf2 ) / ( k + 1.0e+00_real64 )
      if ( k <= n ) then
        cqn(k) = cqf0
      end if
      cqf2 = cqf1
      cqf1 = cqf0
    end do
    do k = 0, n
      cqn(k) = cqn(k) * cq0 / cqf0
    end do
  end if

  cqd(0) = ( cqn(1) - z * cqn(0) ) / ( z * z - 1.0e+00_real64 )
  do k = 1, n
    cqd(k) = ( k * z * cqn(k) - k * cqn(k-1) ) / ( z * z - 1.0e+00_real64 )
  end do

  return
end subroutine clqn
!> @brief subroutine comelp.
!> @return None.
!>
!> @param hk [in] Argument hk.
!> @param ck [inout] Argument ck.
!> @param ce [inout] Argument ce.
pure subroutine comelp ( hk, ck, ce )

!*****************************************************************************80
!
!! COMELP computes complete elliptic integrals K(k) and E(k).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    07 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) HK, the modulus.  0 <= HK <= 1._real64
!
!    Output, real(real64) CK, CE, the values of K(HK) and E(HK).
!
  implicit none

  real(real64) ae
  real(real64) ak
  real(real64) be
  real(real64) bk
  real(real64), intent(inout) :: ce
  real(real64), intent(inout) :: ck
  real(real64), intent(in) :: hk
  real(real64) pk

  pk = 1.0e+00_real64 - hk * hk

  if ( hk == 1.0e+00_real64 ) then

    ck = 1.0e+300_real64
    ce = 1.0e+00_real64

  else

    ak = ((( &
        0.01451196212e+00_real64   * pk &
      + 0.03742563713e+00_real64 ) * pk &
      + 0.03590092383e+00_real64 ) * pk &
      + 0.09666344259e+00_real64 ) * pk &
      + 1.38629436112e+00_real64

    bk = ((( &
        0.00441787012e+00_real64   * pk &
      + 0.03328355346e+00_real64 ) * pk &
      + 0.06880248576e+00_real64 ) * pk &
      + 0.12498593597e+00_real64 ) * pk &
      + 0.5e+00_real64

    ck = ak - bk * log ( pk )

    ae = ((( &
        0.01736506451e+00_real64   * pk &
      + 0.04757383546e+00_real64 ) * pk &
      + 0.0626060122e+00_real64  ) * pk &
      + 0.44325141463e+00_real64 ) * pk &
      + 1.0e+00_real64

    be = ((( &
        0.00526449639e+00_real64   * pk &
      + 0.04069697526e+00_real64 ) * pk &
      + 0.09200180037e+00_real64 ) * pk &
      + 0.2499836831e+00_real64  ) * pk

    ce = ae - be * log ( pk )

  end if

  return
end subroutine comelp
!> @brief subroutine cpbdn.
!> @return None.
!>
!> @param n [in] Argument n.
!> @param z [in] Argument z.
!> @param cpb [inout] Argument cpb.
!> @param cpd [inout] Argument cpd.
subroutine cpbdn ( n, z, cpb, cpd )

!*****************************************************************************80
!
!! CPBDN: parabolic cylinder function Dn(z) and Dn'(z) for complex argument.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    29 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) N, the order.
!
!    Input, complex ( kind = real64 ) Z, the argument.
!
!    Output, complex ( kind = real64 ) CPB(0:N), CPD(0:N), the values of Dn(z) 
!    and Dn'(z).
!
  implicit none

  integer(int32), intent(in) :: n

  real(real64) a0
  complex(real64) c0
  complex(real64) ca0
  complex(real64) cf
  complex(real64) cf0
  complex(real64) cf1
  complex(real64) cfa
  complex(real64) cfb
  complex(real64), intent(inout) :: cpb(0:n)
  complex(real64), intent(inout) :: cpd(0:n)
  complex(real64) cs0
  integer(int32) k
  integer(int32) m
  integer(int32) n0
  integer(int32) n1
  integer(int32) nm1
  real(real64) pi
  real(real64) x
  complex(real64), intent(in) :: z
  complex(real64) z1

  pi = 3.141592653589793e+00_real64
  x = real ( z, kind = real64 )
  a0 = abs ( z )
  c0 = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
  ca0 = exp ( -0.25e+00_real64 * z * z )

  if ( 0 <= n ) then

    cf0 = ca0
    cf1 = z * ca0
    cpb(0) = cf0
    cpb(1) = cf1
    do k = 2, n
      cf = z * cf1 - ( k - 1.0e+00_real64 ) * cf0
      cpb(k) = cf
      cf0 = cf1
      cf1 = cf
    end do

  else

    n0 = -n

    if ( x <= 0.0e+00_real64 .or. abs ( z ) == 0.0e+00_real64 ) then

      cf0 = ca0
      cpb(0) = cf0
      z1 = - z
      if ( a0 <= 7.0e+00_real64 ) then
        call cpdsa ( -1, z1, cf1 )
      else
        call cpdla ( -1, z1, cf1 )
      end if
      cf1 = sqrt ( 2.0e+00_real64 * pi ) / ca0 - cf1
      cpb(1) = cf1
      do k = 2, n0
        cf = ( - z * cf1 + cf0 ) / ( k - 1.0e+00_real64 )
        cpb(k) = cf
        cf0 = cf1
        cf1 = cf
      end do

    else

      if ( a0 <= 3.0e+00_real64 ) then

        call cpdsa ( -n0, z, cfa )
        cpb(n0) = cfa
        n1 = n0 + 1
        call cpdsa ( -n1, z, cfb )
        cpb(n1) = cfb
        nm1 = n0 - 1
        do k = nm1, 0, -1
          cf = z * cfa + ( k + 1.0e+00_real64 ) * cfb
          cpb(k) = cf
          cfb = cfa
          cfa = cf
        end do

      else

        m = 100 + abs ( n )
        cfa = c0
        cfb = cmplx ( 1.0e-30_real64, 0.0e+00_real64, kind = real64 )
        do k = m, 0, -1
          cf = z * cfb + ( k + 1.0e+00_real64 ) * cfa
          if ( k <= n0 ) then
            cpb(k) = cf
          end if
          cfa = cfb
          cfb = cf
        end do
        cs0 = ca0 / cf
        do k = 0, n0
          cpb(k) = cs0 * cpb(k)
        end do

      end if

    end if

  end if

  cpd(0) = -0.5e+00_real64 * z * cpb(0)

  if ( 0 <= n ) then
    do k = 1, n
      cpd(k) = -0.5e+00_real64 * z * cpb(k) + k * cpb(k-1)
    end do
  else
    do k = 1, n0
      cpd(k) = 0.5e+00_real64 * z * cpb(k) - cpb(k-1)
    end do
  end if

  return
end subroutine cpbdn
!> @brief subroutine cpdla.
!> @return None.
!>
!> @param n [in] Argument n.
!> @param z [in] Argument z.
!> @param cdn [inout] Argument cdn.
pure subroutine cpdla ( n, z, cdn )

!****************************************************************************80
!
!! CPDLA computes complex parabolic cylinder function Dn(z) for large argument.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    07 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer N, the order.
!
!    Input, complex ( kind = real64 ) Z, the argument.
!
!    Output, complex ( kind = real64 ) CDN, the function value.
!
  implicit none

  complex(real64) cb0
  complex(real64), intent(inout) :: cdn
  complex(real64) cr
  integer(int32) k
  integer(int32), intent(in) :: n
  complex(real64), intent(in) :: z

  cb0 = z ** n * exp ( -0.25e+00_real64 * z * z )
  cr = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
  cdn = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )

  do k = 1, 16

    cr = -0.5e+00_real64 * cr * ( 2.0e+00_real64 * k - n - 1.0e+00_real64 ) &
      * ( 2.0e+00_real64 * k - n - 2.0e+00_real64 ) / ( k * z * z )

    cdn = cdn + cr

    if ( abs ( cr ) < abs ( cdn ) * 1.0e-12_real64 ) then
      exit
    end if

  end do

  cdn = cb0 * cdn

  return
end subroutine cpdla
!> @brief subroutine cpdsa.
!> @return None.
!>
!> @param n [in] Argument n.
!> @param z [in] Argument z.
!> @param cdn [inout] Argument cdn.
subroutine cpdsa ( n, z, cdn )

!*****************************************************************************80
!
!! CPDSA computes complex parabolic cylinder function Dn(z) for small argument.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    29 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) N, the order.
!
!    Input, complex ( kind = real64 ) Z, the argument.
!
!    Output, complex ( kind = real64 ) CDN, the value of DN(z).
!
  implicit none

  complex(real64) ca0
  complex(real64) cb0
  complex(real64), intent(inout) :: cdn
  complex(real64) cdw
  complex(real64) cr
  real(real64) eps
  real(real64) g0
  real(real64) g1
  real(real64) ga0
  real(real64) gm
  integer(int32) m
  integer(int32), intent(in) :: n
  real(real64) pd
  real(real64) pi
  real(real64) sq2
  real(real64) va0
  real(real64) vm
  real(real64) vt
  real(real64) xn
  complex(real64), intent(in) :: z

  eps = 1.0e-15_real64
  pi = 3.141592653589793e+00_real64
  sq2 = sqrt ( 2.0e+00_real64 )
  ca0 = exp ( - 0.25e+00_real64 * z * z )
  va0 = 0.5e+00_real64 * ( 1.0e+00_real64 - n )

  if ( n == 0 ) then

    cdn = ca0

  else

    if ( abs ( z ) == 0.0e+00_real64 ) then

      if ( va0 <= 0.0e+00_real64 .and. va0 == int ( va0 ) ) then
        cdn = 0.0e+00_real64
      else
        call gaih ( va0, ga0 )
        pd = sqrt ( pi ) / ( 2.0e+00_real64 ** ( -0.5e+00_real64 * n ) * ga0 )
        cdn = cmplx ( pd, 0.0e+00_real64, kind = real64 )
      end if

    else

      xn = - n
      call gaih ( xn, g1 )
      cb0 = 2.0e+00_real64 ** ( -0.5e+00_real64 * n - 1.0e+00_real64 ) * ca0 / g1
      vt = -0.5e+00_real64 * n
      call gaih ( vt, g0 )
      cdn = cmplx ( g0, 0.0e+00_real64, kind = real64 )
      cr = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )

      do m = 1, 250
        vm = 0.5e+00_real64 * ( m - n )
        call gaih ( vm, gm )
        cr = - cr * sq2 * z / m
        cdw = gm * cr
        cdn = cdn + cdw
        if ( abs ( cdw ) < abs ( cdn ) * eps ) then
          exit
        end if
      end do

      cdn = cb0 * cdn

    end if

  end if

  return
end subroutine cpdsa
!> @brief subroutine cpsi.
!> @return None.
!>
!> @param x [inout] Argument x.
!> @param y [inout] Argument y.
!> @param psr [inout] Argument psr.
!> @param psi [inout] Argument psi.
pure subroutine cpsi ( x, y, psr, psi )

!*****************************************************************************80
!
!! CPSI computes the psi function for a complex argument.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    16 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) X, Y, the real and imaginary parts 
!    of the argument.
!
!    Output, real(real64) PSR, PSI, the real and imaginary parts
!    of the function value.
!
  implicit none

  real(real64), parameter, dimension ( 8 ) :: a = [&
    -0.8333333333333e-01_real64, 0.83333333333333333e-02_real64, &
    -0.39682539682539683e-02_real64, 0.41666666666666667e-02_real64, &
    -0.75757575757575758e-02_real64, 0.21092796092796093e-01_real64, &
    -0.83333333333333333e-01_real64, 0.4432598039215686e+00_real64]
  real(real64) ct2
  integer(int32) k
  integer(int32) n
  real(real64) pi
  real(real64), intent(inout) :: psi
  real(real64), intent(inout) :: psr
  real(real64) ri
  real(real64) rr
  real(real64) th
  real(real64) tm
  real(real64) tn
  real(real64), intent(inout) :: x
  real(real64) x0
  real(real64) x1
  real(real64), intent(inout) :: y
  real(real64) y1
  real(real64) z0
  real(real64) z2

  pi = 3.141592653589793e+00_real64

  if ( y == 0.0e+00_real64 .and. x == int ( x ) .and. x <= 0.0e+00_real64 ) then

    psr = 1.0e+300_real64
    psi = 0.0e+00_real64

  else

    if ( x < 0.0e+00_real64 ) then
      x1 = x
      y1 = y
      x = -x
      y = -y
    end if

    x0 = x

    if ( x < 8.0e+00_real64 ) then
      n = 8 - int ( x )
      x0 = x + n
    end if

    if ( x0 == 0.0e+00_real64 ) then
      if ( y /= 0.0e+00_real64 ) then
        th = 0.5e+00_real64 * pi
      else
        th = 0.0e+00_real64
      end if
    else
      th = atan ( y / x0 )
    end if

    z2 = x0 * x0 + y * y
    z0 = sqrt ( z2 )
    psr = log ( z0 ) - 0.5e+00_real64 * x0 / z2
    psi = th + 0.5e+00_real64 * y / z2
    do k = 1, 8
      psr = psr + a(k) * z2 ** ( - k ) * cos ( 2.0e+00_real64 * k * th )
      psi = psi - a(k) * z2 ** ( - k ) * sin ( 2.0e+00_real64 * k * th )
    end do

    if ( x < 8.0e+00_real64 ) then
      rr = 0.0e+00_real64
      ri = 0.0e+00_real64
      do k = 1, n
        rr = rr + ( x0 - k ) / ( ( x0 - k ) ** 2.0e+00_real64 + y * y )
        ri = ri + y / ( ( x0 - k ) ** 2.0e+00_real64 + y * y )
      end do
      psr = psr - rr
      psi = psi + ri
    end if

    if ( x1 < 0.0e+00_real64 ) then
      tn = tan ( pi * x )
      tm = tanh ( pi * y )
      ct2 = tn * tn + tm * tm
      psr = psr + x / ( x * x + y * y ) + pi * ( tn - tn * tm * tm ) / ct2
      psi = psi - y / ( x * x + y * y ) - pi * tm * ( 1.0e+00_real64 + tn * tn ) / ct2
      x = x1
      y = y1
    end if

  end if

  return
end subroutine cpsi
!> @brief subroutine csphik.
!> @return None.
!>
!> @param n [in] Argument n.
!> @param z [in] Argument z.
!> @param nm [inout] Argument nm.
!> @param csi [inout] Argument csi.
!> @param cdi [inout] Argument cdi.
!> @param csk [inout] Argument csk.
!> @param cdk [inout] Argument cdk.
subroutine csphik ( n, z, nm, csi, cdi, csk, cdk )

!*****************************************************************************80
!
!! CSPHIK: complex modified spherical Bessel functions and derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    29 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
! 
!  Parameters:
!
!    Input, integer(int32) N, the order of in(z) and kn(z).
!
!    Input, complex ( kind = real64 ) Z, the argument.
!
!    Output, integer(int32) NM, the highest order computed.
!
!    Output, complex ( kind = real64 ) CSI(0:N), CDI(0:N), CSK(0:N), CDK(0:N),
!    the values of in(z), in'(z), kn(z), kn'(z).
!
  implicit none

  integer(int32), intent(in) :: n

  real(real64) a0
  complex(real64) ccosh1
  complex(real64), intent(inout) :: cdi(0:n)
  complex(real64), intent(inout) :: cdk(0:n)
  complex(real64) cf
  complex(real64) cf0
  complex(real64) cf1
  complex(real64) ci
  complex(real64) cs
  complex(real64), intent(inout) :: csi(0:n)
  complex(real64) csi0
  complex(real64) csi1
  complex(real64) csinh1
  complex(real64), intent(inout) :: csk(0:n)
  integer(int32) k
  integer(int32) m
  ! integer(int32) msta1
  ! integer(int32) msta2
  integer(int32), intent(inout) :: nm
  real(real64) pi
  complex(real64), intent(in) :: z

  pi = 3.141592653589793e+00_real64
  a0 = abs ( z )    
  nm = n

  if ( a0 < 1.0e-60_real64 ) then
    do k = 0, n
      csi(k) = 0.0e+00_real64
      cdi(k) = 0.0e+00_real64
      csk(k) = 1.0e+300_real64
      cdk(k) = -1.0e+300_real64
    end do
    csi(0) = 1.0e+00_real64
    cdi(1) = 0.3333333333333333e+00_real64
    return
  end if

  ci = cmplx ( 0.0e+00_real64, 1.0e+00_real64, kind = real64 )
  csinh1 = sin ( ci * z ) / ci
  ccosh1 = cos ( ci * z )
  csi0 = csinh1 / z
  csi1 = ( - csinh1 / z + ccosh1 ) / z
  csi(0) = csi0
  csi(1) = csi1

  if ( 2 <= n ) then

    m = msta1 ( a0, 200 )
    if ( m < n ) then
      nm = m
    else
      m = msta2 ( a0, n, 15 )
    end if

    cf0 = 0.0e+00_real64
    cf1 = 1.0e+00_real64-100
    do k = m, 0, -1
      cf = ( 2.0e+00_real64 * k + 3.0e+00_real64 ) * cf1 / z + cf0
      if ( k <= nm ) then
        csi(k) = cf
      end if
      cf0 = cf1
      cf1 = cf
    end do

    if ( abs ( csi0 ) <= abs ( csi1 ) ) then
      cs = csi1 / cf0
    else
      cs = csi0 / cf
    end if

    do k = 0, nm
      csi(k) = cs * csi(k)
    end do

  end if

  cdi(0) = csi(1)
  do k = 1, nm
    cdi(k) = csi(k-1) - ( k + 1.0e+00_real64 ) * csi(k) / z
  end do

  csk(0) = 0.5e+00_real64 * pi / z * exp ( - z )
  csk(1) = csk(0) * ( 1.0e+00_real64 + 1.0e+00_real64 / z )
  do k = 2, nm
    if ( abs ( csi(k-2) ) < abs ( csi(k-1) ) ) then
      csk(k) = ( 0.5e+00_real64 * pi / ( z * z ) - csi(k) * csk(k-1) ) / csi(k-1)
    else
      csk(k) = ( csi(k) * csk(k-2) + ( k - 0.5e+00_real64 ) * pi / z ** 3 ) / csi(k-2)
    end if
  end do

  cdk(0) = -csk(1)
  do k = 1, nm
    cdk(k) = - csk(k-1) - ( k + 1.0e+00_real64 ) * csk(k) / z
  end do

  return
end subroutine csphik
!> @brief subroutine csphjy.
!> @return None.
!>
!> @param n [in] Argument n.
!> @param z [in] Argument z.
!> @param nm [inout] Argument nm.
!> @param csj [inout] Argument csj.
!> @param cdj [inout] Argument cdj.
!> @param csy [inout] Argument csy.
!> @param cdy [inout] Argument cdy.
subroutine csphjy ( n, z, nm, csj, cdj, csy, cdy )

!*****************************************************************************80
!
!! CSPHJY: spherical Bessel functions jn(z) and yn(z) for complex argument.
!
!  Discussion:
!
!    This procedure computes spherical Bessel functions jn(z) and yn(z)
!    and their derivatives for a complex argument.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    01 August 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
! 
!  Parameters:
!
!    Input, integer(int32) N, the order of jn(z) and yn(z).
!
!    Input, complex ( kind = real64 ) Z, the argument.
!
!    Output, integer(int32) NM, the highest order computed.
!
!    Output, complex ( kind = real64 ) CSJ(0:N0, CDJ(0:N), CSY(0:N), CDY(0:N),
!    the values of jn(z), jn'(z), yn(z), yn'(z).
!
  implicit none

  integer(int32), intent(in) :: n

  real(real64) a0
  complex(real64), intent(inout) :: csj(0:n)
  complex(real64), intent(inout) :: cdj(0:n)
  complex(real64), intent(inout) :: csy(0:n)
  complex(real64), intent(inout) :: cdy(0:n)
  complex(real64) cf
  complex(real64) cf0
  complex(real64) cf1
  complex(real64) cs
  complex(real64) csa
  complex(real64) csb
  integer(int32) k
  integer(int32) m
  ! integer(int32) msta1
  ! integer(int32) msta2
  integer(int32), intent(inout) :: nm
  complex(real64), intent(in) :: z

  a0 = abs ( z )
  nm = n

  if ( a0 < 1.0e-60_real64 ) then
    do k = 0, n
      csj(k) = 0.0e+00_real64
      cdj(k) = 0.0e+00_real64
      csy(k) = -1.0e+300_real64
      cdy(k) = 1.0e+300_real64
    end do
    csj(0) = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    cdj(1) = cmplx ( 0.333333333333333e+00_real64, 0.0e+00_real64, kind = real64 )
    return
  end if

  csj(0) = sin ( z ) / z
  csj(1) = ( csj(0) - cos ( z ) ) / z

  if ( 2 <= n ) then
    csa = csj(0)
    csb = csj(1)
    m = msta1 ( a0, 200 )
    if ( m < n ) then
      nm = m
    else
      m = msta2 ( a0, n, 15 )
    end if
    cf0 = 0.0e+00_real64
    cf1 = 1.0e+00_real64-100
    do k = m, 0, -1
      cf = ( 2.0e+00_real64 * k + 3.0e+00_real64 ) * cf1 / z - cf0
      if ( k <= nm ) then
        csj(k) = cf
      end if
      cf0 = cf1
      cf1 = cf
    end do

    if ( abs ( csa ) <= abs ( csb ) ) then
      cs = csb / cf0
    else
      cs = csa / cf
    end if

    do k = 0, nm
      csj(k) = cs * csj(k)
    end do

  end if

  cdj(0) = ( cos ( z ) - sin ( z ) / z ) / z
  do k = 1, nm
    cdj(k) = csj(k-1) - ( k + 1.0e+00_real64 ) * csj(k) / z
  end do
  csy(0) = - cos ( z ) / z
  csy(1) = ( csy(0) - sin ( z ) ) / z
  cdy(0) = ( sin ( z ) + cos ( z ) / z ) / z
  cdy(1) = ( 2.0e+00_real64 * cdy(0) - cos ( z ) )  / z

  do k = 2, nm
    if ( abs ( csj(k-2) ) < abs ( csj(k-1) ) ) then 
      csy(k) = ( csj(k) * csy(k-1) - 1.0e+00_real64 / ( z * z ) ) / csj(k-1)
    else
      csy(k) = ( csj(k) * csy(k-2) &
        - ( 2.0e+00_real64 * k - 1.0e+00_real64 ) / z ** 3 ) / csj(k-2)
    end if
  end do

  do k = 2, nm
    cdy(k) = csy(k-1) - ( k + 1.0e+00_real64 ) * csy(k) / z
  end do

  return
end subroutine csphjy
!> @brief subroutine cv0.
!> @return None.
!>
!> @param kd [in] Argument kd.
!> @param m [in] Argument m.
!> @param q [in] Argument q.
!> @param a0 [inout] Argument a0.
subroutine cv0 ( kd, m, q, a0 )

!*****************************************************************************80
!
!! CV0 computes the initial characteristic value of Mathieu functions.
!
!  Discussion:
!
!    This procedure computes the initial characteristic value of Mathieu 
!    functions for m <= 12 or q <= 300 or q <= m*m.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!   03 August 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) KD, the case code:
!    1, for cem(x,q)  ( m = 0,2,4,...)
!    2, for cem(x,q)  ( m = 1,3,5,...)
!    3, for sem(x,q)  ( m = 1,3,5,...)
!    4, for sem(x,q)  ( m = 2,4,6,...)
!
!    Input, integer(int32) M, the order of the functions.
!
!    Input, real(real64) Q, the parameter of the functions.
!
!    Output, real(real64) A0, the characteristic value.
!
  implicit none

  real(real64), intent(inout) :: a0
  integer(int32), intent(in) :: kd
  integer(int32), intent(in) :: m
  real(real64), intent(in) :: q
  real(real64) q2

  q2 = q * q

  if ( m == 0 ) then

    if ( q <= 1.0e+00_real64 ) then

      a0 = ((( &
          0.0036392e+00_real64   * q2 &
        - 0.0125868e+00_real64 ) * q2 &
        + 0.0546875e+00_real64 ) * q2 &
        - 0.5e+00_real64 )       * q2

    else if ( q <= 10.0e+00_real64 ) then

      a0 = (( &
          3.999267e-03_real64   * q &
        - 9.638957e-02_real64 ) * q &
        - 0.88297e+00_real64 )  * q &
        + 0.5542818e+00_real64 

    else

      call cvql ( kd, m, q, a0 )

    end if

  else if ( m == 1 ) then

    if ( q <= 1.0e+00_real64 .and. kd == 2 ) then

      a0 = ((( &
        - 6.51e-04_real64 * q &
        - 0.015625e+00_real64 ) *  q &
        - 0.125e+00_real64 ) * q &
        + 1.0e+00_real64 ) * q &
        + 1.0e+00_real64 

    else if ( q <= 1.0e+00_real64 .and. kd == 3 ) then

      a0 = ((( &
        - 6.51e-04_real64 * q &
        + 0.015625e+00_real64 ) * q &
        - 0.125e+00_real64 ) * q &
        - 1.0e+00_real64 ) * q &
        + 1.0e+00_real64
 
    else if ( q <= 10.0e+00_real64 .and. kd == 2 ) then

      a0 = ((( &
        - 4.94603e-04_real64 * q &
        + 1.92917e-02_real64 ) * q &
        - 0.3089229e+00_real64 ) * q &
        + 1.33372e+00_real64 ) * q &
        + 0.811752e+00_real64 

    else if ( q <= 10.0e+00_real64 .and. kd == 3 ) then

      a0 = (( &
          1.971096e-03_real64 * q &
        - 5.482465e-02_real64 ) * q &
        - 1.152218e+00_real64 ) * q &
        + 1.10427e+00_real64 

    else

      call cvql ( kd, m, q, a0 )

    end if

  else if ( m == 2 ) then

    if ( q <= 1.0e+00_real64 .and. kd == 1 ) then

      a0 = ((( &
        - 0.0036391e+00_real64   * q2 &
        + 0.0125888e+00_real64 ) * q2 &
        - 0.0551939e+00_real64 ) * q2 &
        + 0.416667e+00_real64 )  * q2 + 4.0e+00_real64 

    else if ( q <= 1.0e+00_real64 .and. kd == 4 ) then

      a0 = (  &
          0.0003617e+00_real64 * q2  &
        - 0.0833333e+00_real64 ) * q2 + 4.0e+00_real64 

    else if ( q <= 15.0e+00_real64 .and. kd == 1 ) then

      a0 = ((( &
          3.200972e-04_real64    * q &
        - 8.667445e-03_real64 )  * q &
        - 1.829032e-04_real64 )  * q &
        + 0.9919999e+00_real64 ) * q &
        + 3.3290504e+00_real64 

    else if ( q <= 10.0e+00_real64 .and. kd == 4 ) then

      a0 = (( &
          2.38446e-03_real64 * q &
        - 0.08725329e+00_real64 ) * q &
        - 4.732542e-03_real64 ) * q &
        + 4.00909e+00_real64 

    else

      call cvql ( kd, m, q, a0 )

    end if

  else if ( m == 3 ) then

    if ( q <= 1.0e+00_real64 .and. kd == 2 ) then
      a0 = (( &
          6.348e-04_real64 * q &
        + 0.015625e+00_real64 ) * q &
        + 0.0625_real64 ) * q2  &
        + 9.0e+00_real64 
    else if ( q <= 1.0e+00_real64 .and. kd == 3 ) then
      a0 = (( &
          6.348e-04_real64 * q &
        - 0.015625e+00_real64 ) * q &
        + 0.0625e+00_real64 ) * q2 &
        + 9.0e+00_real64 
    else if ( q <= 20.0e+00_real64 .and. kd == 2 ) then
      a0 = ((( &
          3.035731e-04_real64 * q &
        - 1.453021e-02_real64 ) * q &
        + 0.19069602e+00_real64 ) * q &
        - 0.1039356e+00_real64 ) * q &
        + 8.9449274e+00_real64 
    else if ( q <= 15.0e+00_real64 .and. kd == 3 ) then
      a0 = (( &
          9.369364e-05_real64 * q &
        - 0.03569325e+00_real64 ) * q &
        + 0.2689874e+00_real64 ) * q &
        + 8.771735e+00_real64 
    else
      call cvql ( kd, m, q, a0 )
    end if

  else if ( m == 4 ) then

    if ( q <= 1.0e+00_real64 .and. kd == 1 ) then
      a0 = (( &
        - 2.1e-06_real64 * q2 &
        + 5.012e-04_real64 ) * q2 &
        + 0.0333333_real64 ) * q2 &
        + 16.0e+00_real64
    else if ( q <= 1.0e+00_real64 .and. kd == 4 ) then
      a0 = (( &
          3.7e-06_real64 * q2 &
        - 3.669e-04_real64 ) * q2 &
        + 0.0333333e+00_real64 ) * q2 &
        + 16.0e+00_real64
    else if ( q <= 25.0e+00_real64 .and. kd == 1 ) then
      a0 = ((( &
          1.076676e-04_real64 * q &
        - 7.9684875e-03_real64 ) * q &
        + 0.17344854e+00_real64 ) * q &
        - 0.5924058e+00_real64 ) * q &
        + 16.620847e+00_real64
    else if ( q <= 20.0e+00_real64 .and. kd == 4 ) then
      a0 = (( &
        - 7.08719e-04_real64 * q &
        + 3.8216144e-03_real64 ) * q &
        + 0.1907493e+00_real64 ) * q &
        + 15.744e+00_real64
    else
      call cvql ( kd, m, q, a0 )
    end if

  else if ( m == 5 ) then

    if ( q <= 1.0e+00_real64 .and. kd == 2 ) then
      a0 = (( &
          6.8e-6_real64 * q &
        + 1.42e-05_real64 ) * q2 &
        + 0.0208333e+00_real64 ) * q2 &
        + 25.0e+00_real64
    else if ( q <= 1.0e+00_real64 .and. kd == 3 ) then
      a0 = (( &
        - 6.8e-06_real64 * q &
        + 1.42e-05_real64 ) * q2 &
        + 0.0208333e+00_real64 ) * q2 &
        + 25.0e+00_real64
    else if ( q <= 35.0e+00_real64 .and. kd == 2 ) then
      a0 = ((( &
          2.238231e-05_real64 * q &
        - 2.983416e-03_real64 ) * q &
        + 0.10706975e+00_real64 ) * q &
        - 0.600205e+00_real64 ) * q &
        + 25.93515e+00_real64
    else if ( q <= 25.0e+00_real64 .and. kd == 3 ) then
      a0 = (( &
        - 7.425364e-04_real64 * q &
        + 2.18225e-02_real64 ) * q &
        + 4.16399e-02_real64 ) * q &
        + 24.897e+00_real64
    else
      call cvql ( kd, m, q, a0 )
    end if

  else if ( m == 6 ) then

    if ( q <= 1.0e+00_real64 ) then
      a0 = ( 0.4e-06_real64 * q2 + 0.0142857_real64 ) * q2 + 36.0e+00_real64
    else if ( q <= 40.0e+00_real64 .and. kd == 1 ) then
      a0 = ((( &
        - 1.66846e-05_real64 * q &
        + 4.80263e-04_real64 ) * q &
        + 2.53998e-02_real64 ) * q &
        - 0.181233e+00_real64 ) * q  &
        + 36.423e+00_real64
    else if ( q <= 35.0e+00_real64 .and. kd == 4 ) then
      a0 = (( &
        - 4.57146e-04_real64 * q &
        + 2.16609e-02_real64 ) * q &
        - 2.349616e-02_real64 ) * q &
        + 35.99251e+00_real64
    else
      call cvql ( kd, m, q, a0 )
    end if

  else if ( m == 7 ) then

    if ( q <= 10.0e+00_real64 ) then
      call cvqm ( m, q, a0 )
    else if ( q <= 50.0e+00_real64 .and. kd == 2 ) then
      a0 = ((( &
        - 1.411114e-05_real64 * q &
        + 9.730514e-04_real64 ) * q &
        - 3.097887e-03_real64 ) * q &
        + 3.533597e-02_real64 ) * q &
        + 49.0547e+00_real64
    else if ( q <= 40.0e+00_real64 .and. kd == 3 ) then
      a0 = (( &
        - 3.043872e-04_real64 * q &
        + 2.05511e-02_real64 ) * q &
        - 9.16292e-02_real64 ) * q &
        + 49.19035e+00_real64
    else
      call cvql ( kd, m, q, a0 )
    end if

  else if ( 8 <= m ) then

    if ( q <= 3.0e+00_real64 * m ) then
      call cvqm ( m, q, a0 )
    else if ( m * m .lt. q ) then
      call cvql ( kd, m, q, a0 )
    else if ( m == 8 .and. kd == 1 ) then
      a0 = ((( &
          8.634308e-06_real64 * q &
        - 2.100289e-03_real64 ) * q &
        + 0.169072e+00_real64 ) * q &
        - 4.64336e+00_real64 ) * q &
        + 109.4211e+00_real64
    else if ( m == 8 .and. kd == 4 ) then
      a0 = (( &
        - 6.7842e-05_real64 * q &
        + 2.2057e-03_real64 ) * q &
        + 0.48296e+00_real64 ) * q &
        + 56.59e+00_real64
    else if ( m == 9 .and. kd == 2 ) then
      a0 = ((( &
          2.906435e-06_real64 * q &
        - 1.019893e-03_real64 ) * q &
        + 0.1101965e+00_real64 ) * q &
        - 3.821851e+00_real64 ) * q &
        + 127.6098e+00_real64
    else if ( m == 9 .and. kd == 3 ) then
      a0 = (( &
        - 9.577289e-05_real64 * q &
        + 0.01043839e+00_real64 ) * q &
        + 0.06588934e+00_real64 ) * q &
        + 78.0198e+00_real64
    else if ( m == 10 .and. kd == 1 ) then
      a0 = ((( &
          5.44927e-07_real64 * q &
        - 3.926119e-04_real64 ) * q &
        + 0.0612099e+00_real64 ) * q &
        - 2.600805e+00_real64 ) * q &
        + 138.1923e+00_real64
    else if ( m == 10 .and. kd == 4 ) then
      a0 = (( &
        - 7.660143e-05_real64 * q &
        + 0.01132506e+00_real64 ) * q &
        - 0.09746023e+00_real64 ) * q &
        + 99.29494e+00_real64
    else if ( m == 11 .and. kd == 2 ) then
      a0 = ((( &
        - 5.67615e-07_real64 * q &
        + 7.152722e-06_real64 ) * q &
        + 0.01920291e+00_real64 ) * q &
        - 1.081583e+00_real64 ) * q &
        + 140.88e+00_real64
    else if ( m == 11 .and. kd == 3 ) then
      a0 = (( &
        - 6.310551e-05_real64 * q &
        + 0.0119247e+00_real64 ) * q &
        - 0.2681195e+00_real64 ) * q &
        + 123.667e+00_real64 
    else if ( m == 12 .and. kd == 1 ) then
      a0 = ((( &
        - 2.38351e-07_real64 * q &
        - 2.90139e-05_real64 ) * q &
        + 0.02023088e+00_real64 ) * q &
        - 1.289e+00_real64 ) * q &
        + 171.2723e+00_real64
    else if ( m == 12 .and. kd == 4 ) then
      a0 = ((( &
          3.08902e-07_real64 * q &
        - 1.577869e-04_real64 ) * q &
        + 0.0247911e+00_real64 ) * q &
        - 1.05454e+00_real64 ) * q  &
        + 161.471e+00_real64

    end if

  end if

  return
end subroutine cv0
!> @brief subroutine cva1.
!> @return None.
!>
!> @param kd [in] Argument kd.
!> @param m [in] Argument m.
!> @param q [in] Argument q.
!> @param cv [inout] Argument cv.
pure subroutine cva1 ( kd, m, q, cv )

!*****************************************************************************80
!
!! CVA1 computes a sequence of characteristic values of Mathieu functions.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    25 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) KD, the case code.
!    1, for cem(x,q)  ( m = 0,2,4,��� )
!    2, for cem(x,q)  ( m = 1,3,5,��� )
!    3, for sem(x,q)  ( m = 1,3,5,��� )
!    4, for sem(x,q)  ( m = 2,4,6,��� )
!
!    Input, integer(int32) M, the maximum order of the Mathieu functions.
!
!    Input, real(real64) Q, the parameter of the Mathieu functions.
!
!    Output, real(real64) CV(*), characteristic values.
!    For KD = 1, CV(1), CV(2), CV(3),..., correspond to
!    the characteristic values of cem for m = 0,2,4,...
!    For KD = 2, CV(1), CV(2), CV(3),..., correspond to
!    the characteristic values of cem for m = 1,3,5,...
!    For KD = 3, CV(1), CV(2), CV(3),..., correspond to
!    the characteristic values of sem for m = 1,3,5,...
!    For KD = 4, CV(1), CV(2), CV(3),..., correspond to
!    the characteristic values of sem for m = 0,2,4,...
!       
  implicit none

  real(real64), intent(inout) :: cv(200)
  real(real64) d(500)
  real(real64) e(500)
  real(real64) eps
  real(real64) f(500)
  real(real64) g(200)
  real(real64) h(200)
  integer(int32) i
  integer(int32) ic
  integer(int32) icm
  integer(int32) j
  integer(int32) k
  integer(int32) k1
  integer(int32), intent(in) :: kd
  integer(int32), intent(in) :: m
  integer(int32) nm
  integer(int32) nm1
  real(real64), intent(in) :: q
  real(real64) s
  real(real64) t
  real(real64) t1
  real(real64) x1
  real(real64) xa
  real(real64) xb

  eps = 1.0e-14_real64

  if ( kd == 4 ) then
    icm = m / 2
  else
    icm = int ( m / 2 ) + 1
  end if

  if ( q == 0.0e+00_real64 ) then

    if ( kd == 1 ) then
      do ic = 1, icm
        cv(ic) = 4.0e+00_real64 * ( ic - 1.0e+00_real64 ) ** 2
      end do
    else if ( kd /= 4 ) then
      do ic = 1, icm
        cv(ic) = ( 2.0e+00_real64 * ic - 1.0e+00_real64 ) ** 2
      end do
    else
      do ic = 1, icm
        cv(ic) = 4.0e+00_real64 * ic * ic
      end do
    end if

  else

    nm = int ( 10e+00_real64 + 1.5e+00_real64 * m + 0.5e+00_real64 * q )
    e(1) = 0.0e+00_real64
    f(1) = 0.0e+00_real64

    if ( kd == 1 ) then

      d(1) = 0.0e+00_real64
      do i = 2, nm
        d(i) = 4.0e+00_real64 * ( i - 1.0e+00_real64 ) ** 2
        e(i) = q
        f(i) = q * q
      end do
      e(2) = sqrt ( 2.0e+00_real64 ) * q
      f(2) = 2.0e+00_real64 * q * q

    else if ( kd /= 4 ) then

      d(1) = 1.0e+00_real64 + ( -1.0e+00_real64 ) ** kd * q
      do i = 2, nm
        d(i) = ( 2.0e+00_real64 * i - 1.0e+00_real64 ) ** 2
        e(i) = q
        f(i) = q * q
      end do

    else

      d(1) = 4.0e+00_real64
      do i = 2, nm
        d(i) = 4.0e+00_real64 * i * i
        e(i) = q
        f(i) = q * q
      end do

    end if

    xa = d(nm) + abs ( e(nm) )
    xb = d(nm) - abs ( e(nm) )

    nm1 = nm - 1
    do i = 1, nm1
      t = abs ( e(i) ) + abs ( e(i+1) )
      t1 = d(i) + t
      xa = max ( xa, t1 )
      t1 = d(i) - t
      xb = min ( xb, t1 )
    end do

    do i = 1, icm
      g(i) = xa
      h(i) = xb
    end do

    do k = 1, icm

      do k1 = k, icm
        if ( g(k1) < g(k) ) then
          g(k) = g(k1)
          exit
        end if
      end do

      if ( k /= 1 .and. h(k) < h(k-1) ) then
        h(k) = h(k-1)
      end if

      do

        x1 = ( g(k) + h(k) ) /2.0e+00_real64
        cv(k) = x1

        if ( abs ( ( g(k) - h(k) ) / x1 ) < eps ) then
          exit
        end if

        j = 0
        s = 1.0e+00_real64
        do i = 1, nm
          if ( s == 0.0e+00_real64 ) then
            s = s + 1.0e-30_real64
          end if
          t = f(i) / s
          s = d(i) - t - x1
          if ( s < 0.0e+00_real64 ) then
            j = j + 1
          end if
        end do

        if ( j < k ) then
          h(k) = x1
        else
          g(k) = x1
          if ( icm <= j ) then
            g(icm) = x1
          else
            h(j+1) = max ( h(j+1), x1 )
            g(j) = min ( g(j), x1 )
          end if
        end if

      end do

      cv(k) = x1

    end do

  end if

  return
end subroutine cva1
!> @brief subroutine cva2.
!> @return None.
!>
!> @param kd [in] Argument kd.
!> @param m [in] Argument m.
!> @param q [in] Argument q.
!> @param a [inout] Argument a.
subroutine cva2 ( kd, m, q, a )

!*****************************************************************************80
!
!! CVA2 computes a specific characteristic value of Mathieu functions.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    22 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) KD, the case code:
!    1, for cem(x,q)  ( m = 0,2,4,...)
!    2, for cem(x,q)  ( m = 1,3,5,...)
!    3, for sem(x,q)  ( m = 1,3,5,...)
!    4, for sem(x,q)  ( m = 2,4,6,...)
!
!    Input, integer(int32) M, the order of the Mathieu functions.
!
!    Input, real(real64) Q, the parameter of the Mathieu functions.
!
!    Output, real(real64) A, the characteristic value.
!
  implicit none

  real(real64), intent(inout) :: a
  real(real64) a1
  real(real64) a2
  real(real64) delta
  integer(int32) i
  integer(int32) iflag
  integer(int32), intent(in) :: kd
  integer(int32), intent(in) :: m
  integer(int32) ndiv
  integer(int32) nn
  real(real64), intent(in) :: q
  real(real64) q1
  real(real64) q2
  real(real64) qq

  if ( m <= 12 .or. q <= 3.0e+00_real64 * m .or. m * m < q ) then

    call cv0 ( kd, m, q, a )

    if ( q /= 0.0e+00_real64 ) then
      iflag = 1
      call refine ( kd, m, q, a, iflag )
    end if

  else

    ndiv = 10
    delta = ( m - 3.0e+00_real64 ) * m / real ( ndiv, kind = real64 )

    if ( ( q - 3.0e+00_real64 * m ) <= ( m * m - q ) ) then

      do

        nn = int ( ( q - 3.0e+00_real64 * m ) / delta ) + 1
        delta = ( q - 3.0e+00_real64 * m ) / nn
        q1 = 2.0e+00_real64 * m
        call cvqm ( m, q1, a1 )
        q2 = 3.0e+00_real64 * m
        call cvqm ( m, q2, a2 )
        qq = 3.0e+00_real64 * m
 
        do i = 1, nn

          qq = qq + delta
          a = ( a1 * q2 - a2 * q1 + ( a2 - a1 ) * qq ) / ( q2 - q1 )
 
          if ( i == nn ) then
            iflag = -1
          else
            iflag = 1
          end if

          call refine ( kd, m, qq, a, iflag )
          q1 = q2
          q2 = qq
          a1 = a2
          a2 = a

        end do

        if ( iflag /= -10 ) then
          exit
        end if

        ndiv = ndiv * 2
        delta = ( m - 3.0e+00_real64 ) * m / real ( ndiv, kind = real64 )

      end do

    else

      do

        nn = int ( ( m * m - q ) / delta ) + 1
        delta = ( m * m - q ) / nn
        q1 = m * ( m - 1.0e+00_real64 )
        call cvql ( kd, m, q1, a1 )
        q2 = m * m
        call cvql ( kd, m, q2, a2 )
        qq = m * m

        do i = 1, nn

          qq = qq - delta
          a = ( a1 * q2 - a2 * q1 + ( a2 - a1 ) * qq ) / ( q2 - q1 )

          if ( i == nn ) then
            iflag = -1
          else
            iflag = 1
          end if

          call refine ( kd, m, qq, a, iflag )
          q1 = q2
          q2 = qq
          a1 = a2
          a2 = a

        end do

        if ( iflag /= -10 ) then
          exit
        end if

        ndiv = ndiv * 2
        delta = ( m - 3.0e+00_real64 ) * m / real ( ndiv, kind = real64 )

      end do

    end if

  end if

  return
end subroutine cva2
!> @brief subroutine cvf.
!> @return None.
!>
!> @param kd [in] Argument kd.
!> @param m [in] Argument m.
!> @param q [in] Argument q.
!> @param a [in] Argument a.
!> @param mj [in] Argument mj.
!> @param f [inout] Argument f.
pure subroutine cvf ( kd, m, q, a, mj, f )

!*****************************************************************************80
!
!! CVF computes F for the characteristic equation of Mathieu functions.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    16 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) KD, the case code:
!    1, for cem(x,q)  ( m = 0,2,4,...)
!    2, for cem(x,q)  ( m = 1,3,5,...)
!    3, for sem(x,q)  ( m = 1,3,5,...)
!    4, for sem(x,q)  ( m = 2,4,6,...)
!
!    Input, integer(int32) M, the order of the Mathieu functions.
!
!    Input, real(real64) Q, the parameter of the Mathieu functions.
!
!    Input, real(real64) A, the characteristic value.
!
!    Input, integer(int32) MJ, ?
!
!    Output, real(real64) F, the value of the function for the
!    characteristic equation.
!
  implicit none

  real(real64), intent(in) :: a
  real(real64) b
  real(real64), intent(inout) :: f
  integer(int32) ic
  integer(int32) j
  integer(int32) j0
  integer(int32) jf
  integer(int32), intent(in) :: kd
  integer(int32) l
  integer(int32) l0
  integer(int32), intent(in) :: m
  integer(int32), intent(in) :: mj
  real(real64), intent(in) :: q
  real(real64) t0
  real(real64) t1
  real(real64) t2

  b = a
  ic = int ( m / 2 )
  l = 0
  l0 = 0
  j0 = 2
  jf = ic

  if ( kd == 1 ) then
    l0 = 2
    j0 = 3
  else if ( kd == 2 .or. kd == 3 ) then
    l = 1
  else if ( kd == 4 ) then
    jf = ic - 1
  end if

  t1 = 0.0e+00_real64
  do j = mj, ic + 1, -1
    t1 = - q * q / ( ( 2.0e+00_real64 * j + l ) ** 2 - b + t1 )
  end do

  if ( m <= 2 ) then

    t2 = 0.0e+00_real64

    if ( kd == 1 ) then
      if ( m == 0 ) then
        t1 = t1 + t1
      else if ( m == 2 ) then
        t1 = - 2.0e+00_real64 * q * q / ( 4.0e+00_real64 - b + t1 ) - 4.0e+00_real64
      end if
    else if ( kd == 2 ) then
      if ( m == 1 ) then
        t1 = t1 + q
      end if
    else if ( kd == 3 ) then
      if ( m == 1 ) then
        t1 = t1 - q
      end if
    end if

  else

    if ( kd == 1 ) then
      t0 = 4.0e+00_real64 - b + 2.0e+00_real64 * q * q / b
    else if ( kd == 2 ) then
      t0 = 1.0e+00_real64 - b + q
    else if ( kd == 3 ) then
      t0 = 1.0e+00_real64 - b - q
    else if ( kd == 4 ) then
      t0 = 4.0e+00_real64 - b
    end if

    t2 = - q * q / t0
    do j = j0, jf
      t2 = - q * q / ( ( 2.0e+00_real64 * j - l - l0 ) ** 2 - b + t2 )
    end do

  end if

  f = ( 2.0e+00_real64 * ic + l ) ** 2 + t1 + t2 - b

  return
end subroutine cvf
!> @brief subroutine cvql.
!> @return None.
!>
!> @param kd [in] Argument kd.
!> @param m [in] Argument m.
!> @param q [in] Argument q.
!> @param a0 [inout] Argument a0.
pure subroutine cvql ( kd, m, q, a0 )

!*****************************************************************************80
!
!! CVQL computes the characteristic value of Mathieu functions for q <= 3*m.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    10 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) KD, the case code:
!    1, for cem(x,q)  ( m = 0,2,4,...)
!    2, for cem(x,q)  ( m = 1,3,5,...)
!    3, for sem(x,q)  ( m = 1,3,5,...)
!    4, for sem(x,q)  ( m = 2,4,6,...)
!
!    Input, integer(int32) M, the order of the Mathieu functions.
!
!    Input, real(real64) Q, the parameter value.
!
!    Output, real(real64) A0, the initial characteristic value.
!
  implicit none

  real(real64), intent(inout) :: a0
  real(real64) c1
  real(real64) cv1
  real(real64) cv2
  real(real64) d1
  real(real64) d2
  real(real64) d3
  real(real64) d4
  integer(int32), intent(in) :: kd
  integer(int32), intent(in) :: m
  real(real64) p1
  real(real64) p2
  real(real64), intent(in) :: q
  real(real64) w
  real(real64) w2
  real(real64) w3
  real(real64) w4
  real(real64) w6

  if ( kd == 1 .or. kd == 2 ) then
    w = 2.0e+00_real64 * m + 1.0e+00_real64
  else
    w = 2.0e+00_real64 * m - 1.0e+00_real64
  end if

  w2 = w * w
  w3 = w * w2
  w4 = w2 * w2
  w6 = w2 * w4
  d1 = 5.0e+00_real64 + 34.0e+00_real64 / w2 + 9.0e+00_real64 / w4
  d2 = ( 33.0e+00_real64 + 410.0e+00_real64 / w2 + 405.0e+00_real64 / w4 ) / w
  d3 = ( 63.0e+00_real64 + 1260.0e+00_real64 / w2 + 2943.0e+00_real64 / w4 + 486.0e+00_real64 / w6 ) / w2
  d4 = ( 527.0e+00_real64 + 15617.0e+00_real64 / w2 + 69001.0e+00_real64 / w4 &
    + 41607.0e+00_real64 / w6 ) / w3
  c1 = 128.0e+00_real64
  p2 = q / w4
  p1 = sqrt ( p2 )
  cv1 = - 2.0e+00_real64 * q + 2.0e+00_real64 * w * sqrt ( q ) &
    - ( w2 + 1.0e+00_real64 ) / 8.0e+00_real64
  cv2 = ( w + 3.0e+00_real64 / w ) + d1 / ( 32.0e+00_real64 * p1 ) + d2 &
    / ( 8.0e+00_real64 * c1 * p2 )
  cv2 = cv2 + d3 / ( 64.0e+00_real64 * c1 * p1 * p2 ) + d4 &
    / ( 16.0e+00_real64 * c1 * c1 * p2 * p2 )
  a0 = cv1 - cv2 / ( c1 * p1 )

  return
end subroutine cvql
!> @brief subroutine cvqm.
!> @return None.
!>
!> @param m [in] Argument m.
!> @param q [in] Argument q.
!> @param a0 [inout] Argument a0.
pure subroutine cvqm ( m, q, a0 )

!*****************************************************************************80
!
!! CVQM computes the characteristic value of Mathieu functions for q <= m*m.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    07 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) M, the order of the Mathieu functions.
!
!    Input, real(real64) Q, the parameter value.
!
!    Output, real(real64) A0, the initial characteristic value.
!
  implicit none

  real(real64), intent(inout) :: a0
  real(real64) hm1
  real(real64) hm3
  real(real64) hm5
  integer(int32), intent(in) :: m
  real(real64), intent(in) :: q

  hm1 = 0.5e+00_real64 * q / ( m * m - 1.0e+00_real64 )
  hm3 = 0.25e+00_real64 * hm1 ** 3 / ( m * m - 4.0e+00_real64 )
  hm5 = hm1 * hm3 * q / ( ( m * m - 1.0e+00_real64 ) * ( m * m - 9.0e+00_real64 ) )
  a0 = m * m + q * ( hm1 + ( 5.0e+00_real64 * m * m + 7.0e+00_real64 ) * hm3 &
    + ( 9.0e+00_real64 * m ** 4 + 58.0e+00_real64 * m * m + 29.0e+00_real64 ) * hm5 )

  return
end subroutine cvqm
!> @brief subroutine cy01.
!> @return None.
!>
!> @param kf [in] Argument kf.
!> @param z [in] Argument z.
!> @param zf [inout] Argument zf.
!> @param zd [inout] Argument zd.
pure subroutine cy01 ( kf, z, zf, zd )

!*****************************************************************************80
!
!! CY01 computes complex Bessel functions Y0(z) and Y1(z) and derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    01 August 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer KF, the function choice.
!    0 for ZF = Y0(z) and ZD = Y0'(z);
!    1 for ZF = Y1(z) and ZD = Y1'(z);
!    2 for ZF = Y1'(z) and ZD = Y1''(z).
!
!    Input, complex ( kind = real64 ) Z, the argument.
!
!    Output, complex ( kind = real64 ) ZF, ZD, the values of the requested function 
!    and derivative.
!
  implicit none

  real(real64), parameter, dimension ( 12 ) :: a = [&
    -0.703125e-01_real64, 0.112152099609375e+00_real64, &
    -0.5725014209747314e+00_real64, 0.6074042001273483e+01_real64, &
    -0.1100171402692467e+03_real64, 0.3038090510922384e+04_real64, &
    -0.1188384262567832e+06_real64, 0.6252951493434797e+07_real64, &
    -0.4259392165047669e+09_real64, 0.3646840080706556e+11_real64, &
    -0.3833534661393944e+13_real64, 0.4854014686852901e+15_real64]
  real(real64) a0
  real(real64), parameter, dimension ( 12 ) :: a1 = [&
    0.1171875e+00_real64, -0.144195556640625e+00_real64, &
    0.6765925884246826e+00_real64, -0.6883914268109947e+01_real64, &
    0.1215978918765359e+03_real64, -0.3302272294480852e+04_real64, &
    0.1276412726461746e+06_real64, -0.6656367718817688e+07_real64, &
    0.4502786003050393e+09_real64, -0.3833857520742790e+11_real64, &
    0.4011838599133198e+13_real64, -0.5060568503314727e+15_real64]
  real(real64), parameter, dimension ( 12 ) :: b = [&
    0.732421875e-01_real64, -0.2271080017089844e+00_real64, &
    0.1727727502584457e+01_real64, -0.2438052969955606e+02_real64, &
    0.5513358961220206e+03_real64, -0.1825775547429318e+05_real64, &
    0.8328593040162893e+06_real64, -0.5006958953198893e+08_real64, &
    0.3836255180230433e+10_real64, -0.3649010818849833e+12_real64, &
    0.4218971570284096e+14_real64, -0.5827244631566907e+16_real64]
  real(real64), parameter, dimension ( 12 ) :: b1 = [&
    -0.1025390625e+00_real64, 0.2775764465332031e+00_real64, &
    -0.1993531733751297e+01_real64, 0.2724882731126854e+02_real64, &
    -0.6038440767050702e+03_real64, 0.1971837591223663e+05_real64, &
    -0.8902978767070678e+06_real64, 0.5310411010968522e+08_real64, &
    -0.4043620325107754e+10_real64, 0.3827011346598605e+12_real64, &
    -0.4406481417852278e+14_real64, 0.6065091351222699e+16_real64]
  complex(real64) cbj0
  complex(real64) cbj1
  complex(real64) cby0
  complex(real64) cby1
  complex(real64) cdy0
  complex(real64) cdy1
  complex(real64) ci
  complex(real64) cp
  complex(real64) cp0
  complex(real64) cp1
  complex(real64) cq0
  complex(real64) cq1
  complex(real64) cr
  complex(real64) cs
  complex(real64) ct1
  complex(real64) ct2
  complex(real64) cu
  real(real64) el
  integer(int32) k
  integer(int32) k0
  integer(int32), intent(in) :: kf
  real(real64) pi
  real(real64) rp2
  real(real64) w0
  real(real64) w1
  complex(real64), intent(in) :: z
  complex(real64) z1
  complex(real64) z2
  complex(real64), intent(inout) :: zd
  complex(real64), intent(inout) :: zf

  pi = 3.141592653589793e+00_real64
  el = 0.5772156649015329e+00_real64
  rp2 = 2.0e+00_real64 / pi
  ci = cmplx ( 0.0e+00_real64, 1.0e+00_real64, kind = real64 )
  a0 = abs ( z )
  z2 = z * z
  z1 = z

  if ( a0 == 0.0e+00_real64 ) then

    cbj0 = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    cbj1 = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
    cby0 = cmplx ( -1.0e+30_real64, 0.0e+00_real64, kind = real64 )
    cby1 = cmplx ( -1.0e+30_real64, 0.0e+00_real64, kind = real64 )
    cdy0 = cmplx ( 1.0e+30_real64, 0.0e+00_real64, kind = real64 )
    cdy1 = cmplx ( 1.0e+30_real64, 0.0e+00_real64, kind = real64 )

  else

    if ( real ( z, kind = real64 ) < 0.0e+00_real64) then
      z1 = -z
    end if

    if ( a0 <= 12.0e+00_real64 ) then

      cbj0 = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
      cr = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
      do k = 1, 40
        cr = - 0.25e+00_real64 * cr * z2 / ( k * k )
        cbj0 = cbj0 + cr
        if ( abs ( cr ) < abs ( cbj0 ) * 1.0e-15_real64 ) then
          exit
        end if
      end do

      cbj1 = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
      cr = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
      do k = 1, 40
        cr = -0.25e+00_real64 * cr * z2 / ( k * ( k + 1.0e+00_real64 ) )
        cbj1 = cbj1 + cr
        if ( abs ( cr ) < abs ( cbj1 ) * 1.0e-15_real64 ) then
          exit
        end if
      end do

      cbj1 = 0.5e+00_real64 * z1 * cbj1
      w0 = 0.0e+00_real64
      cr = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
      cs = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
      do k = 1, 40
        w0 = w0 + 1.0e+00_real64 / k
        cr = -0.25e+00_real64 * cr / ( k * k ) * z2
        cp = cr * w0
        cs = cs + cp
        if ( abs ( cp ) < abs ( cs ) * 1.0e-15_real64 ) then
          exit
        end if
      end do

      cby0 = rp2 * ( log ( z1 / 2.0e+00_real64 ) + el ) * cbj0 - rp2 * cs
      w1 = 0.0e+00_real64
      cr = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
      cs = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
      do k = 1, 40
        w1 = w1 + 1.0e+00_real64 / k
        cr = - 0.25e+00_real64 * cr / ( k * ( k + 1 ) ) * z2
        cp = cr * ( 2.0e+00_real64 * w1 + 1.0e+00_real64 / ( k + 1.0e+00_real64 ) )
        cs = cs + cp
        if ( abs ( cp ) < abs ( cs ) * 1.0e-15_real64 ) then
          exit
        end if
      end do

      cby1 = rp2 * ( ( log ( z1 / 2.0e+00_real64 ) + el ) * cbj1 &
        - 1.0e+00_real64 / z1 - 0.25e+00_real64 * z1 * cs )

    else

      if ( a0 < 35.0e+00_real64 ) then
        k0 = 12
      else if ( a0 < 50.0e+00_real64 ) then
        k0 = 10
      else
        k0 = 8
      end if

      ct1 = z1 - 0.25e+00_real64 * pi
      cp0 = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
      do k = 1, k0
        cp0 = cp0 + a(k) * z1 ** ( - 2 * k )
      end do
      cq0 = -0.125e+00_real64 / z1
      do k = 1, k0
        cq0 = cq0 + b(k) * z1 ** ( - 2 * k - 1 )
      end do
      cu = sqrt ( rp2 / z1 )
      cbj0 = cu * ( cp0 * cos ( ct1 ) - cq0 * sin ( ct1 ) )
      cby0 = cu * ( cp0 * sin ( ct1 ) + cq0 * cos ( ct1 ) )
      ct2 = z1 - 0.75e+00_real64 * pi
      cp1 = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
      do k = 1, k0
        cp1 = cp1 + a1(k) * z1 ** ( - 2 * k )
      end do
      cq1 = 0.375e+00_real64 / z1
      do k = 1, k0
        cq1 = cq1 + b1(k) * z1 ** ( - 2 * k - 1 )
      end do
      cbj1 = cu * ( cp1 * cos ( ct2 ) - cq1 * sin ( ct2 ) )
      cby1 = cu * ( cp1 * sin ( ct2 ) + cq1 * cos ( ct2 ) )
 
    end if

    if ( real ( z, kind = real64 ) < 0.0e+00_real64 ) then

      if ( imag ( z ) < 0.0e+00_real64 ) then
        cby0 = cby0 - 2.0e+00_real64 * ci * cbj0
      else
        cby0 = cby0 + 2.0e+00_real64 * ci * cbj0
      end if

      if ( imag ( z ) < 0.0e+00_real64 ) then
        cby1 = - ( cby1 - 2.0e+00_real64 * ci * cbj1 )
      else
        cby1 = - ( cby1 + 2.0e+00_real64 * ci * cbj1 )
      end if
      cbj1 = - cbj1
 
    end if

    cdy0 = - cby1
    cdy1 = cby0 - 1.0e+00_real64 / z * cby1

  end if

  if ( kf == 0 ) then
    zf = cby0
    zd = cdy0
  else if ( kf == 1 ) then
    zf = cby1
    zd = cdy1
  else if ( kf == 2 ) then
    zf = cdy1
    zd = - cdy1 / z - ( 1.0e+00_real64 - 1.0e+00_real64 / ( z * z ) ) * cby1
  end if

  return
end subroutine cy01
!> @brief subroutine cyzo.
!> @return None.
!>
!> @param nt [in] Argument nt.
!> @param kf [in] Argument kf.
!> @param kc [in] Argument kc.
!> @param zo [inout] Argument zo.
!> @param zv [inout] Argument zv.
subroutine cyzo ( nt, kf, kc, zo, zv )

!*****************************************************************************80
!
!! CYZO computes zeros of complex Bessel functions Y0(z) and Y1(z) and Y1'(z).
!
!  Parameters:
!
!    Ths procedure computes the complex zeros of Y0(z), Y1(z) and Y1'(z), 
!    and their associated values at the zeros using the modified Newton's 
!    iteration method.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    22 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) NT, the number of zeros.
!
!    Input, integer(int32) KF, the function choice.
!    0 for Y0(z) and Y1(z0);
!    1 for Y1(z) and Y0(z1);
!    2 for Y1'(z) and Y1(z1').
!
!    Input, integer(int32) KC, complex/real choice.
!    0, for complex roots;
!    1, for real roots.
!
!    Output, real(real64) ZO(NT), ZV(NT), the zeros of Y0(z) or Y1(z) 
!    or Y1'(z), and the value of Y0'(z) or Y1'(z) or Y1(z) at the L-th zero.
!
  implicit none

  integer(int32), intent(in) :: nt

  real(real64) h
  integer(int32) i
  integer(int32) it
  integer(int32) j
  integer(int32), intent(in) :: kc
  integer(int32), intent(in) :: kf
  integer(int32) nr
  real(real64) w
  real(real64) w0
  real(real64) x
  real(real64) y
  complex(real64) z
  complex(real64) zd
  complex(real64) zero
  complex(real64) zf
  complex(real64) zfd
  complex(real64) zgd
  complex(real64), intent(inout) :: zo(nt)
  complex(real64) zp
  complex(real64) zq
  complex(real64), intent(inout) :: zv(nt)
  complex(real64) zw

  if ( kc == 0 ) then
    x = -2.4e+00_real64
    y = 0.54e+00_real64
    h = 3.14e+00_real64
  else if ( kc == 1 ) then
    x = 0.89e+00_real64
    y = 0.0e+00_real64
    h = -3.14e+00_real64
  end if

  if ( kf == 1 ) then
    x = -0.503e+00_real64
  else if ( kf == 2 ) then
    x = 0.577e+00_real64
  end if

  zero = cmplx ( x, y, kind = real64 )

  do nr = 1, nt

    if ( nr == 1 ) then
      z = zero
    else
      z = zo(nr-1) - h
    end if

    it = 0

    do

      it = it + 1
      call cy01 ( kf, z, zf, zd )

      zp = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
      do i = 1, nr - 1
        zp = zp * ( z - zo(i) )
      end do

      zfd = zf / zp

      zq = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
      do i = 1, nr - 1
        zw = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
        do j = 1, nr - 1
          if ( j /= i ) then
            zw = zw * ( z - zo(j) )
          end if
        end do
        zq = zq + zw
      end do

      zgd = ( zd - zq * zfd ) / zp
      z = z - zfd / zgd
      w0 = w
      w = abs ( z )

      if ( 50 < it .or. abs ( ( w - w0 ) / w ) <= 1.0e-12_real64 ) then
        exit
      end if

    end do

    zo(nr) = z

  end do

  do i = 1, nt
    z = zo(i)
    if ( kf == 0 .or. kf == 2 ) then
      call cy01 ( 1, z, zf, zd )
      zv(i) = zf
    else if ( kf == 1 ) then
      call cy01 ( 0, z, zf, zd )
      zv(i) = zf
    end if
  end do

  return
end subroutine cyzo
!> @brief subroutine dvla.
!> @return None.
!>
!> @param va [in] Argument va.
!> @param x [in] Argument x.
!> @param pd [inout] Argument pd.
subroutine dvla ( va, x, pd )

!*****************************************************************************80
!
!! DVLA computes parabolic cylinder functions Dv(x) for large argument.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    06 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
! 
!  Parameters:
!
!    Input, real(real64) X, the argument.
!
!    Input, real(real64) VA, the order.
!
!    Output, real(real64) PD, the function value.
!
  implicit none

  real(real64) a0
  real(real64) ep
  real(real64) eps
  real(real64) gl
  integer(int32) k
  real(real64), intent(inout) :: pd
  real(real64) pi
  real(real64) r
  real(real64), intent(in) :: va
  real(real64) vl
  real(real64), intent(in) :: x
  real(real64) x1

  pi = 3.141592653589793e+00_real64
  eps = 1.0e-12_real64
  ep = exp ( -0.25e+00_real64 * x * x )
  a0 = abs ( x ) ** va * ep
  r = 1.0e+00_real64
  pd = 1.0e+00_real64
  do k = 1, 16
    r = -0.5e+00_real64 * r * ( 2.0e+00_real64 * k - va - 1.0e+00_real64 ) &
      * ( 2.0e+00_real64 * k - va - 2.0e+00_real64 ) / ( k * x * x )
    pd = pd + r
    if ( abs ( r / pd ) < eps ) then
      exit
    end if
  end do

  pd = a0 * pd

  if ( x < 0.0e+00_real64 ) then
    x1 = - x
    call vvla ( va, x1, vl )
    call gamma ( -va, gl )
    pd = pi * vl / gl + cos ( pi * va ) * pd
  end if

  return
end subroutine dvla
!> @brief subroutine dvsa.
!> @return None.
!>
!> @param va [in] Argument va.
!> @param x [in] Argument x.
!> @param pd [inout] Argument pd.
subroutine dvsa ( va, x, pd )

!*****************************************************************************80
!
!! DVSA computes parabolic cylinder functions Dv(x) for small argument.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    07 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) VA, the order.
!
!    Input, real(real64) X, the argument.
!
!    Output, real(real64) PD, the function value.
!
  implicit none

  real(real64) a0
  real(real64) ep
  real(real64) eps
  real(real64) g0
  real(real64) g1
  real(real64) ga0
  real(real64) gm
  integer(int32) m
  real(real64), intent(inout) :: pd
  real(real64) pi
  real(real64) r
  real(real64) r1
  real(real64) sq2
  real(real64), intent(in) :: va
  real(real64) va0
  real(real64) vm
  real(real64) vt
  real(real64), intent(in) :: x

  eps = 1.0e-15_real64
  pi = 3.141592653589793e+00_real64
  sq2 = sqrt ( 2.0e+00_real64 )
  ep = exp ( -0.25e+00_real64 * x * x )
  va0 = 0.5e+00_real64 * ( 1.0e+00_real64 - va )

  if ( va == 0.0e+00_real64 ) then

    pd = ep

  else

    if ( x == 0.0e+00_real64 ) then
      if ( va0 <= 0.0e+00_real64 .and. va0 == int ( va0 ) ) then
        pd = 0.0e+00_real64
      else
        call gamma ( va0, ga0 )
        pd = sqrt ( pi ) / ( 2.0e+00_real64 ** ( -0.5e+00_real64 * va ) * ga0 )
      end if

    else

      call gamma ( -va, g1 )
      a0 = 2.0e+00_real64 ** ( -0.5e+00_real64 * va - 1.0e+00_real64 ) * ep / g1
      vt = -0.5e+00_real64 * va
      call gamma ( vt, g0 )
      pd = g0
      r = 1.0e+00_real64
      do m = 1, 250
        vm = 0.5e+00_real64 * ( m - va )
        call gamma ( vm, gm )
        r = -r * sq2 * x / m
        r1 = gm * r
        pd = pd + r1
        if ( abs ( r1 ) < abs ( pd ) * eps ) then
          exit
        end if
      end do

      pd = a0 * pd

    end if

  end if

  return
end subroutine dvsa
!> @brief subroutine e1xa.
!> @return None.
!>
!> @param x [in] Argument x.
!> @param e1 [inout] Argument e1.
pure subroutine e1xa ( x, e1 )

!*****************************************************************************80
!
!! E1XA computes the exponential integral E1(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    06 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) X, the argument.
!
!    Output, real(real64) E1, the function value.
!
  implicit none

  real(real64), intent(inout) :: e1
  real(real64) es1
  real(real64) es2
  real(real64), intent(in) :: x

  if ( x == 0.0e+00_real64 ) then

    e1 = 1.0e+300_real64

  else if ( x <= 1.0e+00_real64 ) then

    e1 = - log ( x ) + (((( &
        1.07857e-03_real64 * x &
      - 9.76004e-03_real64 ) * x &
      + 5.519968e-02_real64 ) * x &
      - 0.24991055e+00_real64 ) * x &
      + 0.99999193e+00_real64 ) * x &
      - 0.57721566e+00_real64

  else

    es1 = ((( x &
      + 8.5733287401e+00_real64 ) * x &
      +18.059016973e+00_real64  ) * x &
      + 8.6347608925e+00_real64 ) * x &
      + 0.2677737343e+00_real64

    es2 = ((( x &
      +  9.5733223454e+00_real64 ) * x &
      + 25.6329561486e+00_real64 ) * x &
      + 21.0996530827e+00_real64 ) * x &
      +  3.9584969228e+00_real64

    e1 = exp ( - x ) / x * es1 / es2

  end if

  return
end subroutine e1xa
!> @brief subroutine e1xb.
!> @return None.
!>
!> @param x [in] Argument x.
!> @param e1 [inout] Argument e1.
pure subroutine e1xb ( x, e1 )

!*****************************************************************************80
!
!! E1XB computes the exponential integral E1(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    06 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) X, the argument.
!
!    Output, real(real64) E1, the function value.
!
  implicit none

  real(real64), intent(inout) :: e1
  real(real64) ga
  integer(int32) k
  integer(int32) m
  real(real64) r
  real(real64) t
  real(real64) t0
  real(real64), intent(in) :: x

  if ( x == 0.0e+00_real64 ) then

    e1 = 1.0e+300_real64

  else if ( x <= 1.0e+00_real64 ) then

    e1 = 1.0e+00_real64
    r = 1.0e+00_real64

    do k = 1, 25
      r = -r * k * x / ( k + 1.0e+00_real64 )**2
      e1 = e1 + r
      if ( abs ( r ) <= abs ( e1 ) * 1.0e-15_real64 ) then
        exit
      end if
    end do

    ga = 0.5772156649015328e+00_real64
    e1 = - ga - log ( x ) + x * e1

  else

    m = 20 + int ( 80.0e+00_real64 / x )
    t0 = 0.0e+00_real64
    do k = m, 1, -1
      t0 = k / ( 1.0e+00_real64 + k / ( x + t0 ) )
    end do
    t = 1.0e+00_real64 / ( x + t0 )
    e1 = exp ( -x ) * t

  end if

  return
end subroutine e1xb
!> @brief subroutine e1z.
!> @return None.
!>
!> @param z [in] Argument z.
!> @param ce1 [inout] Argument ce1.
pure subroutine e1z ( z, ce1 )

!*****************************************************************************80
!
!! E1Z computes the complex exponential integral E1(z).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    16 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
! 
!  Parameters:
!
!    Input, complex ( kind = real64 ) Z, the argument.
!
!    Output, complex ( kind = real64 ) CE1, the function value.
!
  implicit none

  real(real64) a0
  complex(real64), intent(inout) :: ce1
  complex(real64) cr
  complex(real64) ct
  complex(real64) ct0
  real(real64) el
  integer(int32) k
  real(real64) pi
  real(real64) x
  complex(real64), intent(in) :: z

  pi = 3.141592653589793e+00_real64
  el = 0.5772156649015328e+00_real64
  x = real ( z, kind = real64 )
  a0 = abs ( z )

  if ( a0 == 0.0e+00_real64 ) then
    ce1 = cmplx ( 1.0e+300_real64, 0.0e+00_real64, kind = real64 )
  else if ( a0 <= 10.0e+00_real64 .or. &
    ( x < 0.0e+00_real64 .and. a0 < 20.0e+00_real64 ) ) then
    ce1 = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    cr = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    do k = 1, 150
      cr = - cr * k * z / ( k + 1.0e+00_real64 )**2
      ce1 = ce1 + cr
      if ( abs ( cr ) <= abs ( ce1 ) * 1.0e-15_real64 ) then
        exit
      end if
    end do

    ce1 = - el - log ( z ) + z * ce1

  else

    ct0 = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
    do k = 120, 1, -1
      ct0 = k / ( 1.0e+00_real64 + k / ( z + ct0 ) )
    end do
    ct = 1.0e+00_real64 / ( z + ct0 )

    ce1 = exp ( - z ) * ct
    if ( x <= 0.0e+00_real64 .and. imag ( z ) == 0.0e+00_real64 ) then
      ce1 = ce1 - pi * cmplx ( 0.0e+00_real64, 1.0e+00_real64, kind = real64 )
    end if

  end if

  return
end subroutine e1z
!> @brief subroutine eix.
!> @return None.
!>
!> @param x [in] Argument x.
!> @param ei [inout] Argument ei.
pure subroutine eix ( x, ei )

!*****************************************************************************80
!
!! EIX computes the exponential integral Ei(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    10 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) X, the argument.
!
!    Output, real(real64) EI, the function value.
!
  implicit none

  real(real64), intent(inout) :: ei
  real(real64) ga
  integer(int32) k
  real(real64) r
  real(real64), intent(in) :: x

  if ( x == 0.0e+00_real64 ) then

    ei = -1.0e+300_real64

  else if ( x <= 40.0e+00_real64 ) then

    ei = 1.0e+00_real64
    r = 1.0e+00_real64
    do k = 1, 100
      r = r * k * x / ( k + 1.0e+00_real64 )**2
      ei = ei + r
      if ( abs ( r / ei ) <= 1.0e-15_real64 ) then
        exit
      end if
    end do

    ga = 0.5772156649015328e+00_real64
    ei = ga + log ( x ) + x * ei

  else

    ei = 1.0e+00_real64
    r = 1.0e+00_real64
    do k = 1, 20
      r = r * k / x
      ei = ei + r
    end do
    ei = exp ( x ) / x * ei

  end if

  return
end subroutine eix
!> @brief subroutine elit.
!> @return None.
!>
!> @param hk [in] Argument hk.
!> @param phi [in] Argument phi.
!> @param fe [inout] Argument fe.
!> @param ee [inout] Argument ee.
pure subroutine elit ( hk, phi, fe, ee )

!*****************************************************************************80
!
!! ELIT: complete and incomplete elliptic integrals F(k,phi) and E(k,phi).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    12 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) HK, the modulus, between 0 and 1._real64
!
!    Input, real(real64) PHI, the argument in degrees.
!
!    Output, real(real64) FE, EE, the values of F(k,phi) and E(k,phi).
!
  implicit none

  real(real64) a
  real(real64) a0
  real(real64) b
  real(real64) b0
  real(real64) c
  real(real64) ce
  real(real64) ck
  real(real64) d
  real(real64) d0
  real(real64), intent(inout) :: ee
  real(real64) fac
  real(real64), intent(inout) :: fe
  real(real64) g
  real(real64), intent(in) :: hk
  integer(int32) n
  real(real64), intent(in) :: phi
  real(real64) pi
  real(real64) r

  g = 0.0e+00_real64
  pi = 3.14159265358979e+00_real64
  a0 = 1.0e+00_real64
  b0 = sqrt ( 1.0e+00_real64 - hk * hk )
  d0 = ( pi / 180.0e+00_real64 ) * phi
  r = hk * hk

  if ( hk == 1.0e+00_real64 .and. phi == 90.0e+00_real64 ) then

    fe = 1.0e+300_real64
    ee = 1.0e+00_real64

  else if ( hk == 1.0e+00_real64 ) then

    fe = log ( ( 1.0e+00_real64 + sin ( d0 ) ) / cos ( d0 ) )
    ee = sin ( d0 )

  else

    fac = 1.0e+00_real64
    do n = 1, 40
      a = ( a0 + b0 ) /2.0e+00_real64
      b = sqrt ( a0 * b0 )
      c = ( a0 - b0 ) / 2.0e+00_real64
      fac = 2.0e+00_real64 * fac
      r = r + fac * c * c
      if ( phi /= 90.0e+00_real64 ) then
        d = d0 + atan ( ( b0 / a0 ) * tan ( d0 ) )
        g = g + c * sin( d )
        d0 = d + pi * int ( d / pi + 0.5e+00_real64 )
      end if
      a0 = a
      b0 = b
      if ( c < 1.0e-07_real64 ) then
        exit
      end if
    end do

    ck = pi / ( 2.0e+00_real64 * a )
    ce = pi * ( 2.0e+00_real64 - r ) / ( 4.0e+00_real64 * a )
    if ( phi == 90.0e+00_real64 ) then
      fe = ck
      ee = ce
    else
      fe = d / ( fac * a )
      ee = fe * ce / ck + g
    end if

  end if

  return
end subroutine elit
!> @brief subroutine elit3.
!> @return None.
!>
!> @param phi [in] Argument phi.
!> @param hk [in] Argument hk.
!> @param c [in] Argument c.
!> @param el3 [inout] Argument el3.
pure subroutine elit3 ( phi, hk, c, el3 )

!*****************************************************************************80
!
!! ELIT3 computes the elliptic integral of the third kind.
!
!  Discussion:
!
!    Gauss-Legendre quadrature is used.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    14 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) PHI, the argument in degrees.
!
!    Input, real(real64) HK, the modulus, between 0 and 1._real64
!
!    Input, real(real64) C, the parameter, between 0 and 1._real64
!
!    Output, real(real64) EL3, the value of the elliptic integral
!    of the third kind.
!
  implicit none

  real(real64), intent(in) :: c
  real(real64) c0
  real(real64) c1
  real(real64) c2
  real(real64), intent(inout) :: el3
  real(real64) f1
  real(real64) f2
  real(real64), intent(in) :: hk
  integer(int32) i 
  logical lb1
  logical lb2
  real(real64), intent(in) :: phi
  real(real64), parameter, dimension ( 10 ) :: t = [&
    0.9931285991850949e+00_real64, 0.9639719272779138e+00_real64, &
    0.9122344282513259e+00_real64, 0.8391169718222188e+00_real64, &
    0.7463319064601508e+00_real64, 0.6360536807265150e+00_real64, &
    0.5108670019508271e+00_real64, 0.3737060887154195e+00_real64, &
    0.2277858511416451e+00_real64, 0.7652652113349734e-01_real64]
  real(real64) t1
  real(real64) t2
  real(real64), parameter, dimension ( 10 ) :: w = [&
    0.1761400713915212e-01_real64, 0.4060142980038694e-01_real64, &
    0.6267204833410907e-01_real64, 0.8327674157670475e-01_real64, &
    0.1019301198172404e+00_real64, 0.1181945319615184e+00_real64, &
    0.1316886384491766e+00_real64, 0.1420961093183820e+00_real64, &
    0.1491729864726037e+00_real64, 0.1527533871307258e+00_real64]

  lb1 = ( hk == 1.0e+00_real64 ) .and. ( abs ( phi - 90.0e+00_real64 ) <= 1.0e-08_real64 )

  lb2 = c == 1.0e+00_real64 .and. abs ( phi - 90.0e+00_real64 ) <= 1.0e-08_real64

  if ( lb1 .or. lb2 ) then
    el3 = 1.0e+300_real64
    return
  end if

  c1 = 0.87266462599716e-02_real64 * phi
  c2 = c1

  el3 = 0.0e+00_real64
  do i = 1, 10
    c0 = c2 * t(i)
    t1 = c1 + c0
    t2 = c1 - c0
    f1 = 1.0e+00_real64 / ( ( 1.0e+00_real64 - c * sin(t1) * sin(t1) ) &
      * sqrt ( 1.0e+00_real64 - hk * hk * sin ( t1 ) * sin ( t1 ) ) )
    f2 = 1.0e+00_real64 / ( ( 1.0e+00_real64 - c * sin ( t2 ) * sin ( t2 ) ) &
      * sqrt( 1.0e+00_real64 - hk * hk * sin ( t2 ) * sin ( t2 ) ) )
    el3 = el3 + w(i) * ( f1 + f2 )
  end do

  el3 = c1 * el3

  return
end subroutine elit3
!> @brief function envj.
!> @return Function value.
!>
!> @param n [in] Argument n.
!> @param x [in] Argument x.
function envj ( n, x )

!*****************************************************************************80
!
!! ENVJ is a utility function used by MSTA1 and MSTA2.
!
!  Discussion:
!
!    ENVJ estimates -log(Jn(x)) from the estimate
!    Jn(x) approx 1/sqrt(2*pi*n) * ( e*x/(2*n))^n
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    14 January 2016
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!    Modifications suggested by Vincent Lafage, 11 January 2016._real64
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) N, the order of the Bessel function.
!
!    Input, real(real64) X, the absolute value of the argument.
!
!    Output, real(real64) ENVJ, the value.
!
  implicit none

  real(real64) envj
  real(real64) logten
  integer(int32), intent(in) :: n
  real(real64) n_r8
  real(real64) r8_gamma_log
  real(real64), intent(in) :: x
!
!  Original code
!
  if ( .true. ) then

    envj = 0.5e+00_real64 * log10 ( 6.28e+00_real64 * n ) &
      - n * log10 ( 1.36e+00_real64 * x / n )
!
!  Modification suggested by Vincent Lafage.
!
  else

    n_r8 = real ( n, kind = real64 )
    logten = log ( 10.0e+00_real64 )
    envj = r8_gamma_log ( n_r8 + 1.0e+00_real64 ) / logten - n_r8 * log10 ( x )

  end if

  return
end function envj
!> @brief subroutine enxa.
!> @return None.
!>
!> @param n [in] Argument n.
!> @param x [in] Argument x.
!> @param en [inout] Argument en.
subroutine enxa ( n, x, en )

!*****************************************************************************80
!
!! ENXA computes the exponential integral En(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    07 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) N, the order.
!
!    Input, real(real64) X, the argument.
!
!    Output, real(real64) EN(0:N), the function values.
!
  implicit none

  integer(int32), intent(in) :: n

  real(real64) e1
  real(real64) ek
  real(real64), intent(inout) :: en(0:n)
  integer(int32) k
  real(real64), intent(in) :: x

  en(0) = exp ( - x ) / x 
  call e1xb ( x, e1 )

  en(1) = e1
  do k = 2, n
    ek = ( exp ( - x ) - x * e1 ) / ( k - 1.0e+00_real64 )
    en(k) = ek
    e1 = ek
  end do

  return
end subroutine enxa
!> @brief subroutine enxb.
!> @return None.
!>
!> @param n [in] Argument n.
!> @param x [in] Argument x.
!> @param en [inout] Argument en.
pure subroutine enxb ( n, x, en )

!*****************************************************************************80
!
!! ENXB computes the exponential integral En(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    10 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) N, the order.
!
!    Input, real(real64) X, the argument.
!
!    Output, real(real64) EN(0:N), the function values.
!
  implicit none

  integer(int32), intent(in) :: n

  real(real64), intent(inout) :: en(0:n)
  real(real64) ens
  integer(int32) j
  integer(int32) k
  integer(int32) l
  integer(int32) m
  real(real64) ps
  real(real64) r
  real(real64) rp
  real(real64) s
  real(real64) s0
  real(real64) t
  real(real64) t0
  real(real64), intent(in) :: x

  if ( x == 0.0e+00_real64 ) then

    en(0) = 1.0e+300_real64
    en(1) = 1.0e+300_real64
    do k = 2, n
      en(k) = 1.0e+00_real64 / ( k - 1.0e+00_real64 )
    end do
    return

  else if ( x <= 1.0e+00_real64 ) then

    en(0) = exp ( - x ) / x
    do l = 1, n
      rp = 1.0e+00_real64
      do j = 1, l - 1
        rp = - rp * x / j
      end do
      ps = -0.5772156649015328e+00_real64
      do m = 1, l - 1
        ps = ps + 1.0e+00_real64 / m
      end do
      ens = rp * ( - log ( x ) + ps )
      s = 0.0e+00_real64
      do m = 0, 20
        if ( m /= l - 1 ) then
          r = 1.0e+00_real64
          do j = 1, m
            r = - r * x / j
          end do
          s = s + r / ( m - l + 1.0e+00_real64 )
          if ( abs ( s - s0 ) < abs ( s ) * 1.0e-15_real64 ) then
            exit
          end if
          s0 = s
        end if
      end do

      en(l) = ens - s

    end do

  else

    en(0) = exp ( - x ) / x
    m = 15 + int ( 100.0e+00_real64 / x )
    do l = 1, n
      t0 = 0.0e+00_real64
      do k = m, 1, -1
        t0 = ( l + k - 1.0e+00_real64 ) / ( 1.0e+00_real64 + k / ( x + t0 ) )
      end do
      t = 1.0e+00_real64 / ( x + t0 )
      en(l) = exp ( - x ) * t
    end do

  end if

  return
end subroutine enxb
!> @brief subroutine error.
!> @return None.
!>
!> @param x [in] Argument x.
!> @param err [inout] Argument err.
pure subroutine error ( x, err )

!*****************************************************************************80
!
!! ERROR evaluates the error function.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    07 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) X, the argument.
!
!    Output, real(real64) ERR, the function value.
!
  implicit none

  real(real64) c0
  real(real64) eps
  real(real64) er
  real(real64), intent(inout) :: err
  integer(int32) k
  real(real64) pi
  real(real64) r
  real(real64), intent(in) :: x
  real(real64) x2

  eps = 1.0e-15_real64
  pi = 3.141592653589793e+00_real64
  x2 = x * x

  if ( abs ( x ) < 3.5e+00_real64 ) then

    er = 1.0e+00_real64
    r = 1.0e+00_real64

    do k = 1, 50
      r = r * x2 / ( k + 0.5e+00_real64 )
      er = er + r
      if ( abs ( r ) <= abs ( er ) * eps ) then
        exit
      end if
    end do

    c0 = 2.0e+00_real64 / sqrt ( pi ) * x * exp ( - x2 )
    err = c0 * er

  else

    er = 1.0e+00_real64
    r = 1.0e+00_real64
    do k = 1, 12
      r = - r * ( k - 0.5e+00_real64 ) / x2
      er = er + r
    end do

    c0 = exp ( - x2 ) / ( abs ( x ) * sqrt ( pi ) )

    err = 1.0e+00_real64 - c0 * er
    if ( x < 0.0e+00_real64 ) then
      err = -err
    end if

  end if

  return
end subroutine error
!> @brief subroutine eulera.
!> @return None.
!>
!> @param n [in] Argument n.
!> @param en [inout] Argument en.
pure subroutine eulera ( n, en )

!*****************************************************************************80
!
!! EULERA computes the Euler number En.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    10 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) N, the index of the highest value to compute.
!
!    Output, real(real64) EN(0:N), the Euler numbers up to the N-th value.
!
  implicit none

  integer(int32), intent(in) :: n

  real(real64), intent(inout) :: en(0:n)
  integer(int32) j
  integer(int32) k
  integer(int32) m
  real(real64) r
  real(real64) s

  en(0) = 1.0e+00_real64

  do m = 1, n / 2
    s = 1.0e+00_real64
    do k = 1, m - 1
      r = 1.0e+00_real64
      do j = 1, 2 * k
        r = r * ( 2.0e+00_real64 * m - 2.0e+00_real64 * k + j ) / j
      end do
      s = s + r * en(2*k)
    end do
    en(2*m) = -s
  end do

  return
end subroutine eulera
!> @brief subroutine eulerb.
!> @return None.
!>
!> @param n [in] Argument n.
!> @param en [inout] Argument en.
pure subroutine eulerb ( n, en )

!*****************************************************************************80
!
!! EULERB computes the Euler number En.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    09 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) N, the index of the highest value to compute.
!
!    Output, real(real64) EN(0:N), the Euler numbers up to the N-th value.
!
  implicit none

  integer(int32), intent(in) :: n

  real(real64), intent(inout) :: en(0:n)
  real(real64) hpi
  real(real64) isgn
  integer(int32) k
  integer(int32) m
  real(real64) r1
  real(real64) r2
  real(real64) s

  hpi = 2.0e+00_real64 / 3.141592653589793e+00_real64
  en(0) = 1.0e+00_real64
  en(2) = -1.0e+00_real64
  r1 = -4.0e+00_real64 * hpi ** 3

  do m = 4, n, 2
    r1 = - r1 * ( m - 1 ) * m * hpi * hpi
    r2 = 1.0e+00_real64
    isgn = 1.0e+00_real64
    do k = 3, 1000, 2
      isgn = - isgn
      s = ( 1.0e+00_real64 / k ) ** ( m + 1 )
      r2 = r2 + isgn * s
      if ( s < 1.0e-15_real64 ) then
        exit
      end if
    end do

    en(m) = r1 * r2

    end do

  return
end subroutine eulerb
!> @brief subroutine fcoef.
!> @return None.
!>
!> @param kd [in] Argument kd.
!> @param m [in] Argument m.
!> @param q [in] Argument q.
!> @param a [in] Argument a.
!> @param fc [inout] Argument fc.
pure subroutine fcoef ( kd, m, q, a, fc )

!*****************************************************************************80
!
!! FCOEF: expansion coefficients for Mathieu and modified Mathieu functions.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    01 August 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) KD, the case code.
!    1, for cem(x,q)  ( m = 0,2,4,...)
!    2, for cem(x,q)  ( m = 1,3,5,...)
!    3, for sem(x,q)  ( m = 1,3,5,...)
!    4, for sem(x,q)  ( m = 2,4,6,...)
!
!    Input, integer(int32) M, the order of the Mathieu function.
!
!    Input, real(real64) Q, the parameter of the Mathieu functions.
!
!    Input, real(real64) A, the characteristic value of the Mathieu
!    functions for given m and q.
!
!    Output, real(real64) FC(*), the expansion coefficients of Mathieu
!    functions ( k =  1,2,...,KM ).  FC(1),FC(2),FC(3),... correspond to
!    A0,A2,A4,... for KD = 1 case, 
!    A1,A3,A5,... for KD = 2 case,
!    B1,B3,B5,... for KD = 3 case,
!    B2,B4,B6,... for KD = 4 case.
!
  implicit none

  real(real64), intent(in) :: a
  real(real64) f
  real(real64) f1
  real(real64) f2
  real(real64) f3
  real(real64), intent(inout) :: fc(251)
  integer(int32) i
  integer(int32) j
  integer(int32) k
  integer(int32) kb
  integer(int32), intent(in) :: kd
  integer(int32) km
  integer(int32) l
  integer(int32), intent(in) :: m
  real(real64), intent(in) :: q
  real(real64) qm
  real(real64) s
  real(real64) s0
  real(real64) sp
  real(real64) ss
  real(real64) u
  real(real64) v

  if ( q <= 1.0e+00_real64 ) then
    qm = 7.5e+00_real64 + 56.1e+00_real64 * sqrt ( q ) - 134.7e+00_real64 * q &
      + 90.7e+00_real64 * sqrt ( q ) * q
  else
    qm = 17.0e+00_real64 + 3.1e+00_real64 * sqrt ( q ) - 0.126e+00_real64 * q &
      + 0.0037e+00_real64 * sqrt ( q ) * q
  end if

  km = int ( qm + 0.5e+00_real64 * m )

  if ( q == 0.0e+00_real64 ) then

    do k = 1, km
      fc(k) = 0.0e+00_real64
    end do

    if ( kd == 1 ) then
      fc((m+2)/2) = 1.0e+00_real64
      if (m == 0 ) then
        fc(1) = 1.0e+00_real64 / sqrt ( 2.0e+00_real64 )
      end if
    else if ( kd == 4 ) then
      fc(m/2) = 1.0e+00_real64
    else
      fc((m+1)/2) = 1.0e+00_real64
    end if

    return

  end if

  kb = 0
  s = 0.0e+00_real64
  f = 1.0e-100_real64
  u = 0.0e+00_real64
  fc(km) = 0.0e+00_real64

  if ( kd == 1 ) then

    l = 0

    do k = km, 3, -1

      v = u
      u = f
      f = ( a - 4.0e+00_real64 * k * k ) * u / q - v

      if ( abs ( f ) < abs ( fc(k+1) ) ) then

        kb = k
        fc(1) = 1.0e-100_real64
        sp = 0.0e+00_real64
        f3 = fc(k+1)
        fc(2) = a / q * fc(1)
        fc(3) = ( a - 4.0e+00_real64 ) * fc(2) / q - 2.0e+00_real64 * fc(1)
        u = fc(2)
        f1 = fc(3)

        do i = 3, kb
          v = u
          u = f1
          f1 = ( a - 4.0e+00_real64 * ( i - 1.0e+00_real64 ) ** 2 ) * u / q - v
          fc(i+1) = f1
          if ( i == kb ) then
            f2 = f1
          else
            sp = sp + f1 * f1
          end if
        end do

        sp = sp + 2.0e+00_real64 * fc(1) ** 2 + fc(2) ** 2 + fc(3) ** 2
        ss = s + sp * ( f3 / f2 ) ** 2
        s0 = sqrt ( 1.0e+00_real64 / ss )
        do j = 1, km
          if ( j <= kb + 1 ) then
            fc(j) = s0 * fc(j) * f3 / f2
          else
            fc(j) = s0 * fc(j)
          end if
        end do
        l = 1
        exit
      else
        fc(k) = f
        s = s + f * f
      end if

    end do

    if ( l == 0 ) then
      fc(2) = q * fc(3) / ( a - 4.0e+00_real64 - 2.0e+00_real64 * q * q / a )
      fc(1) = q / a * fc(2)
      s = s + 2.0e+00_real64 * fc(1) ** 2 + fc(2) ** 2
      s0 = sqrt ( 1.0e+00_real64 / s )
      do k = 1, km
        fc(k) = s0 * fc(k)
      end do
    end if

  else if ( kd == 2 .or. kd == 3 ) then

    l = 0

    do k = km, 3, -1

      v = u
      u = f
      f = ( a - ( 2.0e+00_real64 * k - 1 ) ** 2 ) * u / q - v

      if ( abs ( fc(k) ) <= abs ( f ) ) then
        fc(k-1) = f
        s = s + f * f
      else
        kb = k
        f3 = fc(k)
        l = 1
        exit
      end if

    end do

    if ( l == 0 ) then

      fc(1) = q / ( a - 1.0e+00_real64 - ( - 1 ) ** kd * q ) * fc(2)
      s = s + fc(1) * fc(1)
      s0 = sqrt ( 1.0e+00_real64 / s )
      do k = 1, km
        fc(k) = s0 * fc(k)
      end do

    else

      fc(1) = 1.0e-100_real64
      fc(2) = ( a - 1.0e+00_real64 - ( - 1 ) ** kd * q ) / q * fc(1)
      sp = 0.0e+00_real64
      u = fc(1)
      f1 = fc(2)
      do i = 2, kb - 1
        v = u
        u = f1
        f1 = ( a - ( 2.0e+00_real64 * i - 1.0e+00_real64 ) ** 2 ) * u / q - v
        if ( i /= kb - 1 ) then
          fc(i+1) = f1
          sp = sp + f1 * f1
        else
          f2 = f1
        end if
      end do

      sp = sp + fc(1) ** 2 + fc(2) ** 2
      ss = s + sp * ( f3 / f2 ) ** 2
      s0 = 1.0e+00_real64 / sqrt ( ss )
      do j = 1, km
        if ( j < kb ) then
          fc(j) = s0 * fc(j) * f3 / f2
        else
          fc(j) = s0 * fc(j)
        end if
      end do

    end if

  else if ( kd == 4 ) then

    l = 0

    do k = km, 3, -1
      v = u
      u = f
      f = ( a - 4.0e+00_real64 * k * k ) * u / q - v
      if ( abs ( fc(k) ) <= abs ( f ) ) then
        fc(k-1) = f
        s = s + f * f
      else
        kb = k
        f3 = fc(k)
        l = 1
        exit
      end if
    end do

    if ( l == 0 ) then

      fc(1) = q / ( a - 4.0e+00_real64 ) * fc(2)
      s = s + fc(1) * fc(1)
      s0 = sqrt ( 1.0e+00_real64 / s )
      do k = 1, km
        fc(k) = s0 * fc(k)
      end do

    else

      fc(1) = 1.0e-100_real64
      fc(2) = ( a - 4.0e+00_real64 ) / q * fc(1)
      sp = 0.0e+00_real64
      u = fc(1)
      f1 = fc(2)

      do i = 2, kb - 1
        v = u
        u = f1
        f1 = ( a - 4.0e+00_real64 * i * i ) * u / q - v
        if ( i /= kb - 1 ) then
          fc(i+1) = f1
          sp = sp + f1 * f1
        else
          f2 = f1
        end if
      end do

      sp = sp + fc(1) ** 2 + fc(2) ** 2
      ss = s + sp * ( f3 / f2 ) ** 2
      s0 = 1.0e+00_real64 / sqrt ( ss )

      do j = 1, km
        if ( j < kb ) then
          fc(j) = s0 * fc(j) * f3 / f2
        else
          fc(j) = s0 * fc(j)
        end if
      end do

    end if

  end if

  if ( fc(1) < 0.0e+00_real64 ) then
    do j = 1, km
      fc(j) = -fc(j)
    end do
  end if

  return
end subroutine fcoef
!> @brief subroutine fcs.
!> @return None.
!>
!> @param x [in] Argument x.
!> @param c [inout] Argument c.
!> @param s [inout] Argument s.
pure subroutine fcs ( x, c, s )

!*****************************************************************************80
!
!! FCS computes Fresnel integrals C(x) and S(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    17 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
! 
!  Parameters:
!
!    Input, real(real64) X, the argument.
!
!    Output, real(real64) C, S, the function values.
!
  implicit none

  real(real64), intent(inout) :: c
  real(real64) eps
  real(real64) f
  real(real64) f0
  real(real64) f1
  real(real64) g
  integer(int32) k
  integer(int32) m
  real(real64) pi
  real(real64) px
  real(real64) q
  real(real64) r
  real(real64), intent(inout) :: s
  real(real64) su
  real(real64) t
  real(real64) t0
  real(real64) t2
  real(real64), intent(in) :: x
  real(real64) xa

  eps = 1.0e-15_real64
  pi = 3.141592653589793e+00_real64
  xa = abs ( x )
  px = pi * xa
  t = 0.5e+00_real64 * px * xa
  t2 = t * t

  if ( xa == 0.0e+00_real64 ) then

    c = 0.0e+00_real64
    s = 0.0e+00_real64

  else if ( xa < 2.5e+00_real64 ) then

    r = xa
    c = r
    do k = 1, 50
      r = -0.5e+00_real64 * r * ( 4.0e+00_real64 * k - 3.0e+00_real64 ) / k &
        / ( 2.0e+00_real64 * k - 1.0e+00_real64 ) / ( 4.0e+00_real64 * k + 1.0e+00_real64 ) * t2
      c = c + r
      if ( abs ( r ) < abs ( c ) * eps ) then
        exit
      end if
    end do

    s = xa * t / 3.0e+00_real64
    r = s
    do k = 1, 50
      r = - 0.5e+00_real64 * r * ( 4.0e+00_real64 * k - 1.0e+00_real64 ) / k &
        / ( 2.0e+00_real64 * k + 1.0e+00_real64 ) / ( 4.0e+00_real64 * k + 3.0e+00_real64 ) * t2
      s = s + r
      if ( abs ( r ) < abs ( s ) * eps ) then
        if ( x < 0.0e+00_real64 ) then
          c = -c
          s = -s
        end if
        return
      end if
    end do

  else if ( xa < 4.5e+00_real64 ) then

    m = int ( 42.0e+00_real64 + 1.75e+00_real64 * t )
    su = 0.0e+00_real64
    c = 0.0e+00_real64
    s = 0.0e+00_real64
    f1 = 0.0e+00_real64
    f0 = 1.0e-100_real64

    do k = m, 0, -1
      f = ( 2.0e+00_real64 * k + 3.0e+00_real64 ) * f0 / t - f1
      if ( k == int ( k / 2 ) * 2 ) then
        c = c + f
      else
        s = s + f
      end if
      su = su + ( 2.0e+00_real64 * k + 1.0e+00_real64 ) * f * f
      f1 = f0
      f0 = f
    end do

    q = sqrt ( su )
    c = c * xa / q
    s = s * xa / q

  else

    r = 1.0e+00_real64
    f = 1.0e+00_real64
    do k = 1, 20
      r = -0.25e+00_real64 * r * ( 4.0e+00_real64 * k - 1.0e+00_real64 ) &
        * ( 4.0e+00_real64 * k - 3.0e+00_real64 ) / t2
      f = f + r
    end do
    r = 1.0e+00_real64 / ( px * xa )
    g = r
    do k = 1, 12
      r = -0.25e+00_real64 * r * ( 4.0e+00_real64 * k + 1.0e+00_real64 ) &
        * ( 4.0e+00_real64 * k - 1.0e+00_real64 ) / t2
      g = g + r
    end do

    t0 = t - int ( t / ( 2.0e+00_real64 * pi ) ) * 2.0e+00_real64 * pi
    c = 0.5e+00_real64 + ( f * sin ( t0 ) - g * cos ( t0 ) ) / px
    s = 0.5e+00_real64 - ( f * cos ( t0 ) + g * sin ( t0 ) ) / px

  end if

  if ( x < 0.0e+00_real64 ) then
    c = -c
    s = -s
  end if

  return
end subroutine fcs
!> @brief subroutine fcszo.
!> @return None.
!>
!> @param kf [in] Argument kf.
!> @param nt [in] Argument nt.
!> @param zo [inout] Argument zo.
subroutine fcszo ( kf, nt, zo )

!*****************************************************************************80
!
!! FCSZO computes complex zeros of Fresnel integrals C(x) or S(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    17 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) KF, the function code.
!    1 for C(z);
!    2 for S(z)
!
!    Input, integer(int32) NT, the total number of zeros desired.
!
!    Output, complex ( kind = real64 ) Z0(NT), the zeros.
!
  implicit none

  integer(int32), intent(in) :: nt

  integer(int32) i
  integer(int32) it
  integer(int32) j
  integer(int32), intent(in) :: kf
  integer(int32) nr
  real(real64) pi
  real(real64) psq
  real(real64) px
  real(real64) py
  real(real64) w
  real(real64) w0
  complex(real64) z
  complex(real64) zd
  complex(real64) zf
  complex(real64) zfd
  complex(real64) zgd
  complex(real64), intent(inout) :: zo(nt)
  complex(real64) zp
  complex(real64) zq
  complex(real64) zw

  pi = 3.141592653589793e+00_real64

  do nr = 1, nt

    if ( kf == 1 ) then
      psq = sqrt ( 4.0e+00_real64 * nr - 1.0e+00_real64 )
    else
      psq = 2.0e+00_real64 * sqrt ( real ( nr, kind = real64 ) )
    end if

    px = psq - log ( pi * psq ) / ( pi * pi * psq ** 3.0e+00_real64 )
    py = log ( pi * psq ) / ( pi * psq )
    z = cmplx ( px, py, kind = real64 )

    if ( kf == 2 ) then
      if ( nr == 2 ) then
        z = cmplx ( 2.8334e+00_real64, 0.2443e+00_real64, kind = real64 )
      else if ( nr == 3 ) then
        z = cmplx ( 3.4674e+00_real64, 0.2185e+00_real64, kind = real64 )
      else if ( nr == 4 ) then
        z = cmplx ( 4.0025e+00_real64, 0.2008e+00_real64, kind = real64 )
      end if
    end if

    it = 0

    do

      it = it + 1

      if ( kf == 1 ) then
        call cfc ( z, zf, zd )
      else
        call cfs ( z, zf, zd )
      end if

      zp = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
      do i = 1, nr - 1
        zp = zp * ( z - zo(i) )
      end do
      zfd = zf / zp
      zq = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )

      do i = 1, nr - 1
        zw = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
        do j = 1, nr - 1
          if ( j /= i ) then
            zw = zw * ( z - zo(j) )
          end if
        end do
        zq = zq + zw
      end do

      zgd = ( zd - zq * zfd ) / zp
      z = z - zfd / zgd
      w0 = w
      w = cdabs ( z )

      if ( abs ( ( w - w0 ) / w ) <= 1.0e-12_real64 ) then
        exit
      end if

      if ( 50 < it ) then
        exit
      end if

    end do

    zo(nr) = z

  end do

  return
end subroutine fcszo
!> @brief subroutine ffk.
!> @return None.
!>
!> @param ks [in] Argument ks.
!> @param x [in] Argument x.
!> @param fr [inout] Argument fr.
!> @param fi [inout] Argument fi.
!> @param fm [inout] Argument fm.
!> @param fa [inout] Argument fa.
!> @param gr [inout] Argument gr.
!> @param gi [inout] Argument gi.
!> @param gm [inout] Argument gm.
!> @param ga [inout] Argument ga.
pure subroutine ffk ( ks, x, fr, fi, fm, fa, gr, gi, gm, ga )

!*****************************************************************************80
!
!! FFK computes modified Fresnel integrals F+/-(x) and K+/-(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    23 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) KS, the sign code.
!    0, to calculate F+(x) and K+(x);
!    1, to calculate F_(x) and K_(x).
!
!    Input, real(real64) X, the argument.
!
!    Output, real(real64) FR, FI, FM, FA, the values of
!    Re[F+/-(x)], Im[F+/-(x)], |F+/-(x)|, Arg[F+/-(x)]  (Degs.).
!
!    Output, real(real64) GR, GI, GM, GA, the values of
!    Re[K+/-(x)], Im[K+/-(x)], |K+/-(x)|, Arg[K+/-(x)]  (Degs.).
!       
  implicit none

  real(real64) c1
  real(real64) cs
  real(real64) eps
  real(real64), intent(inout) :: fa
  real(real64), intent(inout) :: fi
  real(real64) fi0
  real(real64), intent(inout) :: fm
  real(real64), intent(inout) :: fr
  real(real64), intent(inout) :: ga
  real(real64), intent(inout) :: gi
  real(real64), intent(inout) :: gm
  real(real64), intent(inout) :: gr
  integer(int32) k
  integer(int32), intent(in) :: ks
  integer(int32) m
  real(real64) p2p
  real(real64) pi
  real(real64) pp2
  real(real64) s1
  real(real64) srd
  real(real64) ss
  real(real64), intent(in) :: x
  real(real64) x2
  real(real64) x4
  real(real64) xa
  real(real64) xc
  real(real64) xf
  real(real64) xf0
  real(real64) xf1
  real(real64) xg
  real(real64) xp
  real(real64) xq
  real(real64) xq2
  real(real64) xr
  real(real64) xs
  real(real64) xsu
  real(real64) xw

  srd = 57.29577951308233e+00_real64
  eps = 1.0e-15_real64
  pi = 3.141592653589793e+00_real64
  pp2 = 1.2533141373155e+00_real64
  p2p = 0.7978845608028654e+00_real64
  xa = abs ( x )
  x2 = x * x
  x4 = x2 * x2

  if ( x == 0.0e+00_real64 ) then

    fr = 0.5e+00_real64 * sqrt ( 0.5e+00_real64 * pi )
    fi = ( -1.0e+00_real64 ) ** ks * fr
    fm = sqrt ( 0.25e+00_real64 * pi )
    fa = ( -1.0e+00_real64 ) ** ks * 45.0e+00_real64
    gr = 0.5e+00_real64
    gi = 0.0e+00_real64
    gm = 0.5e+00_real64
    ga = 0.0e+00_real64

  else

    if ( xa <= 2.5e+00_real64 ) then

      xr = p2p * xa
      c1 = xr
      do k = 1, 50
        xr = -0.5e+00_real64 * xr * ( 4.0e+00_real64 * k - 3.0e+00_real64 ) / k &
          / ( 2.0e+00_real64 * k - 1.0e+00_real64 ) &
          / ( 4.0e+00_real64 * k + 1.0e+00_real64 ) * x4
        c1 = c1 + xr
        if ( abs ( xr / c1 ) < eps ) then
          exit
        end if
      end do

      s1 = p2p * xa * xa * xa / 3.0e+00_real64
      xr = s1
      do k = 1, 50
        xr = -0.5e+00_real64 * xr * ( 4.0e+00_real64 * k - 1.0e+00_real64 ) &
          / k / ( 2.0e+00_real64 * k + 1.0e+00_real64 ) &
          / ( 4.0e+00_real64 * k + 3.0e+00_real64 ) * x4
        s1 = s1 + xr
        if ( abs ( xr / s1 ) < eps ) then
          exit
        end if
      end do

    else if ( xa < 5.5e+00_real64 ) then

      m = int ( 42.0e+00_real64 + 1.75e+00_real64 * x2 )
      xsu = 0.0e+00_real64
      xc = 0.0e+00_real64
      xs = 0.0e+00_real64
      xf1 = 0.0e+00_real64
      xf0 = 1.0e-100_real64
      do k = m, 0, -1
        xf = ( 2.0e+00_real64 * k + 3.0e+00_real64 ) * xf0 / x2 - xf1
        if ( k == 2 * int ( k / 2 ) )  then
          xc = xc + xf
        else
          xs = xs + xf
        end if
        xsu = xsu + ( 2.0e+00_real64 * k + 1.0e+00_real64 ) * xf * xf
        xf1 = xf0
        xf0 = xf
      end do
      xq = sqrt ( xsu )
      xw = p2p * xa / xq
      c1 = xc * xw
      s1 = xs * xw

    else

      xr = 1.0e+00_real64
      xf = 1.0e+00_real64
      do k = 1, 12
        xr = -0.25e+00_real64 * xr * ( 4.0e+00_real64 * k - 1.0e+00_real64 ) &
          * ( 4.0e+00_real64 * k - 3.0e+00_real64 ) / x4
        xf = xf + xr
      end do
      xr = 1.0e+00_real64 / ( 2.0e+00_real64 * xa * xa )
      xg = xr
      do k = 1, 12
        xr = -0.25e+00_real64 * xr * ( 4.0e+00_real64 * k + 1.0e+00_real64 ) &
          * ( 4.0e+00_real64 * k - 1.0e+00_real64 ) / x4
        xg = xg + xr
      end do
      c1 = 0.5e+00_real64 + ( xf * sin ( x2 ) - xg * cos ( x2 ) ) &
        / sqrt ( 2.0e+00_real64 * pi ) / xa
      s1 = 0.5e+00_real64 - ( xf * cos ( x2 ) + xg * sin ( x2 ) ) &
        / sqrt ( 2.0e+00_real64 * pi ) / xa

    end if
 
    fr = pp2 * ( 0.5e+00_real64 - c1 )
    fi0 = pp2 * ( 0.5e+00_real64 - s1 )
    fi = ( -1.0e+00_real64 ) ** ks * fi0
    fm = sqrt ( fr * fr + fi * fi )

    if ( 0.0e+00_real64 <= fr ) then
      fa = srd * atan ( fi / fr )
    else if ( 0.0e+00_real64 < fi ) then
      fa = srd * ( atan ( fi / fr ) + pi )
    else if ( fi < 0.0e+00_real64 ) then
      fa = srd * ( atan ( fi / fr ) - pi )
    end if

    xp = x * x + pi / 4.0e+00_real64
    cs = cos ( xp )
    ss = sin ( xp )
    xq2 = 1.0e+00_real64 / sqrt ( pi )
    gr = xq2 * ( fr * cs + fi0 * ss )
    gi = ( -1.0e+00_real64 ) ** ks * xq2 * ( fi0 * cs - fr * ss )
    gm = sqrt ( gr * gr + gi * gi )

    if ( 0.0e+00_real64 <= gr ) then
      ga = srd * atan ( gi / gr )
    else if ( 0.0e+00_real64 < gi ) then
      ga = srd * ( atan ( gi / gr ) + pi )
    else if ( gi < 0.0e+00_real64 ) then
      ga = srd * ( atan ( gi / gr ) - pi )
    end if

    if ( x < 0.0e+00_real64 ) then
      fr = pp2 - fr
      fi = ( -1.0e+00_real64 ) ** ks * pp2 - fi
      fm = sqrt ( fr * fr + fi * fi )
      fa = srd * atan ( fi / fr )
      gr = cos ( x * x ) - gr
      gi = - ( -1.0e+00_real64 ) ** ks * sin ( x * x ) - gi
      gm = sqrt ( gr * gr + gi * gi )
      ga = srd * atan ( gi / gr )
    end if

  end if

  return
end subroutine ffk
!> @brief subroutine gaih.
!> @return None.
!>
!> @param x [in] Argument x.
!> @param ga [inout] Argument ga.
pure subroutine gaih ( x, ga )

!*****************************************************************************80
!
!! GAIH computes the GammaH function.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    09 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) X, the argument.
!
!    Output, real(real64) GA, the function value.
!
  implicit none

  real(real64), intent(inout) :: ga
  integer(int32) k
  integer(int32) m
  integer(int32) m1
  real(real64) pi
  real(real64), intent(in) :: x

  pi = 3.141592653589793e+00_real64

  if ( x == int ( x ) .and. 0.0_real64 < x ) then
    ga = 1.0e+00_real64
    m1 = int ( x - 1.0e+00_real64 )
    do k = 2, m1
      ga = ga * k
    end do
  else if ( x + 0.5e+00_real64 == int ( x + 0.5e+00_real64) .and. 0.0e+00_real64 < x ) then
    m = int ( x )
    ga = sqrt ( pi )
    do k = 1, m
      ga = 0.5e+00_real64 * ga * ( 2.0e+00_real64 * k - 1.0e+00_real64 )
    end do
  end if

  return
end subroutine gaih
!> @brief subroutine gam0.
!> @return None.
!>
!> @param x [in] Argument x.
!> @param ga [inout] Argument ga.
subroutine gam0 ( x, ga )

!*****************************************************************************80
!
!! GAM0 computes the Gamma function for the LAMV function.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    09 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) X, the argument.
!
!    Output, real(real64) GA, the function value.
!   
  implicit none

  real(real64), dimension ( 25 ) :: g = [&
    1.0e+00_real64, &
    0.5772156649015329e+00_real64, &
   -0.6558780715202538e+00_real64, &
   -0.420026350340952e-01_real64, &
    0.1665386113822915e+00_real64, &
   -0.421977345555443e-01_real64, &
   -0.96219715278770e-02_real64, &
    0.72189432466630e-02_real64, &
   -0.11651675918591e-02_real64, &
   -0.2152416741149e-03_real64, &
    0.1280502823882e-03_real64, &
   -0.201348547807e-04_real64, &
   -0.12504934821e-05_real64, &
    0.11330272320e-05_real64, &
   -0.2056338417e-06_real64, &
    0.61160950e-08_real64, &
    0.50020075e-08_real64, &
   -0.11812746e-08_real64, &
    0.1043427e-09_real64, &
    0.77823e-11_real64, &
   -0.36968e-11_real64, &
    0.51e-12_real64, &
   -0.206e-13_real64, &
   -0.54e-14_real64, &
    0.14e-14_real64]
  real(real64), intent(inout) :: ga
  real(real64) gr
  integer(int32) k
  real(real64), intent(in) :: x

  gr = g(25)
  do k = 24, 1, -1
    gr = gr * x + g(k)
  end do

  ga = 1.0e+00_real64 / ( gr * x )

  return
end subroutine gam0
!> @brief subroutine gamma.
!> @return None.
!>
!> @param x [in] Argument x.
!> @param ga [inout] Argument ga.
subroutine gamma ( x, ga )

!*****************************************************************************80
!
!! GAMMA evaluates the Gamma function.
!
!  Licensing:
!
!    The original FORTRAN77 version of this routine is copyrighted by 
!    Shanjie Zhang and Jianming Jin.  However, they give permission to 
!    incorporate this routine into a user program that the copyright 
!    is acknowledged.
!
!  Modified:
!
!    08 September 2007
!
!  Author:
!
!    Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45
!
!  Parameters:
!
!    Input, real(real64) X, the argument.
!    X must not be 0, or any negative integer.
!
!    Output, real(real64) GA, the value of the Gamma function.
!
  implicit none

  real(real64), dimension ( 26 ) :: g = [&
    1.0e+00_real64, &
    0.5772156649015329e+00_real64, &
   -0.6558780715202538e+00_real64, &
   -0.420026350340952e-01_real64, &
    0.1665386113822915e+00_real64, &
   -0.421977345555443e-01_real64, &
   -0.96219715278770e-02_real64, &
    0.72189432466630e-02_real64, &
   -0.11651675918591e-02_real64, &
   -0.2152416741149e-03_real64, &
    0.1280502823882e-03_real64, & 
   -0.201348547807e-04_real64, &
   -0.12504934821e-05_real64, &
    0.11330272320e-05_real64, &
   -0.2056338417e-06_real64, & 
    0.61160950e-08_real64, &
    0.50020075e-08_real64, &
   -0.11812746e-08_real64, &
    0.1043427e-09_real64, & 
    0.77823e-11_real64, &
   -0.36968e-11_real64, &
    0.51e-12_real64, &
   -0.206e-13_real64, &
   -0.54e-14_real64, &
    0.14e-14_real64, &
    0.1e-15_real64]
  real(real64), intent(inout) :: ga
  real(real64) gr
  integer(int32) k
  integer(int32) m
  integer(int32) m1
  real(real64), parameter :: pi = 3.141592653589793e+00_real64
  real(real64) r
  real(real64), intent(in) :: x
  real(real64) z

  if ( x == aint ( x ) ) then

    if ( 0.0e+00_real64 < x ) then
      ga = 1.0e+00_real64
      m1 = int ( x ) - 1
      do k = 2, m1
        ga = ga * k
      end do
    else
      ga = 1.0e+300_real64
    end if

  else

    if ( 1.0e+00_real64 < abs ( x ) ) then
      z = abs ( x )
      m = int ( z )
      r = 1.0e+00_real64
      do k = 1, m
        r = r * ( z - real ( k, kind = real64 ) )
      end do
      z = z - real ( m, kind = real64 )
    else
      z = x
    end if

    gr = g(26)
    do k = 25, 1, -1
      gr = gr * z + g(k)
    end do

    ga = 1.0e+00_real64 / ( gr * z )

    if ( 1.0e+00_real64 < abs ( x ) ) then
      ga = ga * r
      if ( x < 0.0e+00_real64 ) then
        ga = - pi / ( x* ga * sin ( pi * x ) )
      end if
    end if

  end if

  return
end subroutine gamma
!> @brief subroutine gmn.
!> @return None.
!>
!> @param m [in] Argument m.
!> @param n [in] Argument n.
!> @param c [in] Argument c.
!> @param x [in] Argument x.
!> @param bk [in] Argument bk.
!> @param gf [inout] Argument gf.
!> @param gd [inout] Argument gd.
pure subroutine gmn ( m, n, c, x, bk, gf, gd )

!*****************************************************************************80
!
!! GMN computes quantities for oblate radial functions with small argument.
!
!  Discussion:
!
!    This procedure computes Gmn(-ic,ix) and its derivative for oblate
!    radial functions with a small argument.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    15 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) M, the mode parameter;  M = 0, 1, 2, ...
!
!    Input, integer(int32) N, mode parameter, N = M, M + 1, M + 2, ...
!
!    Input, real(real64) C, spheroidal parameter.
!
!    Input, real(real64) X, the argument.
!
!    Input, real(real64) BK(*), coefficients.
!
!    Output, real(real64) GF, GD, the value of Gmn(-C,X) and Gmn'(-C,X).
!
  implicit none

  real(real64), intent(in) :: bk(200)
  real(real64), intent(in) :: c
  real(real64) eps
  real(real64), intent(inout) :: gd
  real(real64) gd0
  real(real64) gd1
  real(real64), intent(inout) :: gf
  real(real64) gf0
  real(real64) gw
  integer(int32) ip
  integer(int32) k
  integer(int32), intent(in) :: m
  integer(int32), intent(in) :: n
  integer(int32) nm
  real(real64), intent(in) :: x
  real(real64) xm

  eps = 1.0e-14_real64

  if ( n - m == 2 * int ( ( n - m ) / 2 ) ) then
    ip = 0
  else
    ip = 1
  end if

  nm = 25 + int ( 0.5e+00_real64 * ( n - m ) + c )
  xm = ( 1.0e+00_real64 + x * x ) ** ( -0.5e+00_real64 * m )
  gf0 = 0.0e+00_real64
  gw = 0.0e+00_real64
  do k = 1, nm
    gf0 = gf0 + bk(k) * x ** ( 2.0e+00_real64 * k - 2.0e+00_real64 )
    if ( abs ( ( gf0 - gw ) / gf0 ) < eps .and. 10 <= k ) then
      exit
    end if
    gw = gf0
  end do

  gf = xm * gf0 * x ** ( 1 - ip )

  gd1 = - m * x / ( 1.0e+00_real64 + x * x ) * gf
  gd0 = 0.0e+00_real64
  gw = 0.0e+00_real64

  do k = 1, nm

    if ( ip == 0 ) then
      gd0 = gd0 + ( 2.0e+00_real64 * k - 1.0e+00_real64 ) * bk(k) &
        * x ** ( 2.0e+00_real64 * k - 2.0e+00_real64 )
    else
      gd0 = gd0 + 2.0e+00_real64 * k * bk(k+1) * x ** ( 2.0e+00_real64 * k - 1.0e+00_real64 )
    end if

    if ( abs ( ( gd0 - gw ) / gd0 ) < eps .and. 10 <= k ) then
      exit
    end if

    gw = gd0

  end do

  gd = gd1 + xm * gd0

  return
end subroutine gmn
!> @brief subroutine herzo.
!> @return None.
!>
!> @param n [in] Argument n.
!> @param x [inout] Argument x.
!> @param w [inout] Argument w.
pure subroutine herzo ( n, x, w )

!*****************************************************************************80
!
!! HERZO computes the zeros the Hermite polynomial Hn(x).
!
!  Discussion:
!
!    This procedure computes the zeros of Hermite polynomial Ln(x)
!    in the interval [-1,+1], and the corresponding
!    weighting coefficients for Gauss-Hermite integration.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    15 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) N, the order of the polynomial.
!
!    Output, real(real64) X(N), the zeros.
!
!    Output, real(real64) W(N), the corresponding weights.
!
  implicit none

  integer(int32), intent(in) :: n

  real(real64) f0
  real(real64) f1
  real(real64) fd
  real(real64) gd
  real(real64) hd
  real(real64) hf
  real(real64) hn
  integer(int32) i
  integer(int32) it
  integer(int32) j
  integer(int32) k
  integer(int32) nr
  real(real64) p
  real(real64) q
  real(real64) r
  real(real64) r1
  real(real64) r2
  real(real64), intent(inout) :: w(n)
  real(real64) wp
  real(real64), intent(inout) :: x(n)
  real(real64) z
  real(real64) z0
  real(real64) zl

  hn = 1.0e+00_real64 / n
  zl = -1.1611e+00_real64 + 1.46e+00_real64 * sqrt ( real ( n, kind = real64 ) )

  do nr = 1, n / 2

    if ( nr == 1 ) then
      z = zl
    else
      z = z - hn * ( n / 2 + 1 - nr )
    end if

    it = 0

    do

      it = it + 1
      z0 = z
      f0 = 1.0e+00_real64
      f1 = 2.0e+00_real64 * z
      do k = 2, n
        hf = 2.0e+00_real64 * z * f1 - 2.0e+00_real64 * ( k - 1.0e+00_real64 ) * f0
        hd = 2.0e+00_real64 * k * f1
        f0 = f1
        f1 = hf
      end do

      p = 1.0e+00_real64
      do i = 1, nr - 1
        p = p * ( z - x(i) )
      end do
      fd = hf / p

      q = 0.0e+00_real64
      do i = 1, nr - 1
        wp = 1.0e+00_real64
        do j = 1, nr - 1
          if ( j /= i ) then
            wp = wp * ( z - x(j) )
          end if
        end do
        q = q + wp
      end do

      gd = ( hd - q * fd ) / p
      z = z - fd / gd

      if ( 40 < it .or. abs ( ( z - z0 ) / z ) <= 1.0e-15_real64 ) then
        exit
      end if

    end do

    x(nr) = z
    x(n+1-nr) = -z
    r = 1.0e+00_real64
    do k = 1, n
      r = 2.0e+00_real64 * r * k
    end do
    w(nr) = 3.544907701811e+00_real64 * r / ( hd * hd )
    w(n+1-nr) = w(nr)

  end do

  if ( n /= 2 * int ( n / 2 ) ) then
    r1 = 1.0e+00_real64
    r2 = 1.0e+00_real64
    do j = 1, n
      r1 = 2.0e+00_real64 * r1 * j
      if ( ( n + 1 ) / 2 <= j ) then
        r2 = r2 * j
      end if
    end do
    w(n/2+1) = 0.88622692545276e+00_real64 * r1 / ( r2 * r2 )
    x(n/2+1) = 0.0e+00_real64
  end if

  return
end subroutine herzo
!> @brief subroutine hygfx.
!> @return None.
!>
!> @param a [inout] Argument a.
!> @param b [inout] Argument b.
!> @param c [in] Argument c.
!> @param x [inout] Argument x.
!> @param hf [inout] Argument hf.
subroutine hygfx ( a, b, c, x, hf )

!*****************************************************************************80
!
!! HYGFX evaluates the hypergeometric function F(A,B,C,X).
!
!  Licensing:
!
!    The original FORTRAN77 version of this routine is copyrighted by 
!    Shanjie Zhang and Jianming Jin.  However, they give permission to 
!    incorporate this routine into a user program that the copyright 
!    is acknowledged.
!
!  Modified:
!
!    08 September 2007
!
!  Author:
!
!    Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45
!
!  Parameters:
!
!    Input, real(real64) A, B, C, X, the arguments of the function.
!    C must not be equal to a nonpositive integer.
!    X < 1._real64
!
!    Output, real HF, the value of the function.
!
  implicit none

  real(real64), intent(inout) :: a
  real(real64) a0
  real(real64) aa
  real(real64), intent(inout) :: b
  real(real64) bb
  real(real64), intent(in) :: c
  real(real64) c0
  real(real64) c1
  real(real64), parameter :: el = 0.5772156649015329e+00_real64
  real(real64) eps
  real(real64) f0
  real(real64) f1
  real(real64) g0
  real(real64) g1
  real(real64) g2
  real(real64) g3
  real(real64) ga
  real(real64) gabc
  real(real64) gam
  real(real64) gb
  real(real64) gbm
  real(real64) gc
  real(real64) gca
  real(real64) gcab
  real(real64) gcb
  real(real64) gm
  real(real64), intent(inout) :: hf
  real(real64) hw
  integer(int32) j
  integer(int32) k
  logical l0
  logical l1
  logical l2
  logical l3
  logical l4
  logical l5
  integer(int32) m
  integer(int32) nm
  real(real64) pa
  real(real64) pb
  real(real64), parameter :: pi = 3.141592653589793e+00_real64
  real(real64) r
  real(real64) r0
  real(real64) r1
  real(real64) rm
  real(real64) rp
  real(real64) sm
  real(real64) sp
  real(real64) sp0
  real(real64), intent(inout) :: x
  real(real64) x1

  l0 = ( c == aint ( c ) ) .and. ( c < 0.0e+00_real64 )
  l1 = ( 1.0e+00_real64 - x < 1.0e-15_real64 ) .and. ( c - a - b <= 0.0e+00_real64 )
  l2 = ( a == aint ( a ) ) .and. ( a < 0.0e+00_real64 )
  l3 = ( b == aint ( b ) ) .and. ( b < 0.0e+00_real64 )
  l4 = ( c - a == aint ( c - a ) ) .and. ( c - a <= 0.0e+00_real64 )
  l5 = ( c - b == aint ( c - b ) ) .and. ( c - b <= 0.0e+00_real64 )

  if ( l0 .or. l1 ) then
    write ( error_unit, '(a)' ) ' '
    write ( error_unit, '(a)' ) 'HYGFX - Fatal error!'
    write ( error_unit, '(a)' ) '  The hypergeometric series is divergent.'
    return
  end if

  if ( 0.95e+00_real64 < x ) then 
    eps = 1.0e-08_real64
  else
    eps = 1.0e-15_real64
  end if
  hw = 0.0e+00_real64

  if ( x == 0.0e+00_real64 .or. a == 0.0e+00_real64 .or. b == 0.0e+00_real64 ) then

    hf = 1.0e+00_real64
    return

  else if ( 1.0e+00_real64 - x == eps .and. 0.0e+00_real64 < c - a - b ) then

    call gamma ( c, gc )
    call gamma ( c - a - b, gcab )
    call gamma ( c - a, gca )
    call gamma ( c - b, gcb )
    hf = gc * gcab /( gca *gcb )
    return

  else if ( 1.0e+00_real64 + x <= eps .and. abs ( c - a + b - 1.0e+00_real64 ) <= eps ) then

    g0 = sqrt ( pi ) * 2.0e+00_real64**( - a )
    call gamma ( c, g1 )
    call gamma ( 1.0e+00_real64 + a / 2.0e+00_real64 - b, g2 )
    call gamma ( 0.5e+00_real64 + 0.5e+00_real64 * a, g3 )
    hf = g0 * g1 / ( g2 * g3 )
    return

  else if ( l2 .or. l3 ) then

    if ( l2 ) then
      nm = int ( abs ( a ) )
    end if

    if ( l3 ) then
      nm = int ( abs ( b ) )
    end if

    hf = 1.0e+00_real64
    r = 1.0e+00_real64

    do k = 1, nm
      r = r * ( a + k - 1.0e+00_real64 ) * ( b + k - 1.0e+00_real64 ) &
        / ( k * ( c + k - 1.0e+00_real64 ) ) * x
      hf = hf + r
    end do

    return

  else if ( l4 .or. l5 ) then

    if ( l4 ) then
      nm = int ( abs ( c - a ) )
    end if

    if ( l5 ) then
      nm = int ( abs ( c - b ) )
    end if

    hf = 1.0e+00_real64
    r  = 1.0e+00_real64
    do k = 1, nm
      r = r * ( c - a + k - 1.0e+00_real64 ) * ( c - b + k - 1.0e+00_real64 ) &
        / ( k * ( c + k - 1.0e+00_real64 ) ) * x
      hf = hf + r
    end do
    hf = ( 1.0e+00_real64 - x )**( c - a - b ) * hf
    return

  end if

  aa = a
  bb = b
  x1 = x
!
!  WARNING: ALTERATION OF INPUT ARGUMENTS A AND B, WHICH MIGHT BE CONSTANTS.
!
  if ( x < 0.0e+00_real64 ) then
    x = x / ( x - 1.0e+00_real64 )
    if ( a < c .and. b < a .and. 0.0e+00_real64 < b ) then
      a = bb
      b = aa
    end if
    b = c - b
  end if

  if ( 0.75e+00_real64 <= x ) then

    gm = 0.0e+00_real64

    if ( abs ( c - a - b - aint ( c - a - b ) ) < 1.0e-15_real64 ) then

      m = int ( c - a - b )
      call gamma ( a, ga )
      call gamma ( b, gb )
      call gamma ( c, gc )
      call gamma ( a + m, gam )
      call gamma ( b + m, gbm )
      call psi ( a, pa )
      call psi ( b, pb )

      if ( m /= 0 ) then
        gm = 1.0e+00_real64
      end if

      do j = 1, abs ( m ) - 1
        gm = gm * j
      end do

      rm = 1.0e+00_real64
      do j = 1, abs ( m )
        rm = rm * j
      end do

      f0 = 1.0e+00_real64
      r0 = 1.0e+00_real64
      r1 = 1.0e+00_real64
      sp0 = 0.0e+00_real64
      sp = 0.0e+00_real64

      if ( 0 <= m ) then

        c0 = gm * gc / ( gam * gbm )
        c1 = - gc * ( x - 1.0e+00_real64 )**m / ( ga * gb * rm )

        do k = 1, m - 1
          r0 = r0 * ( a + k - 1.0e+00_real64 ) * ( b + k - 1.0e+00_real64 ) &
            / ( k * ( k - m ) ) * ( 1.0e+00_real64 - x )
          f0 = f0 + r0
        end do

        do k = 1, m
          sp0 = sp0 + 1.0e+00_real64 / ( a + k - 1.0e+00_real64 ) &
            + 1.0e+00_real64 / ( b + k - 1.0e+00_real64 ) - 1.0e+00_real64 / real ( k, kind = real64 )
        end do

        f1 = pa + pb + sp0 + 2.0e+00_real64 * el + log ( 1.0e+00_real64 - x )
        hw = f1

        do k = 1, 250

          sp = sp + ( 1.0e+00_real64 - a ) / ( k * ( a + k - 1.0e+00_real64 ) ) &
            + ( 1.0e+00_real64 - b ) / ( k * ( b + k - 1.0e+00_real64 ) )

          sm = 0.0e+00_real64
          do j = 1, m
            sm = sm + ( 1.0e+00_real64 - a ) &
              / ( ( j + k ) * ( a + j + k - 1.0e+00_real64 ) ) &
              + 1.0e+00_real64 / ( b + j + k - 1.0e+00_real64 )
          end do

          rp = pa + pb + 2.0e+00_real64 * el + sp + sm + log ( 1.0e+00_real64 - x )

          r1 = r1 * ( a + m + k - 1.0e+00_real64 ) * ( b + m + k - 1.0e+00_real64 ) &
            / ( k * ( m + k ) ) * ( 1.0e+00_real64 - x )

          f1 = f1 + r1 * rp

          if ( abs ( f1 - hw ) < abs ( f1 ) * eps ) then
            exit
          end if

          hw = f1

        end do

        hf = f0 * c0 + f1 * c1

      else if ( m < 0 ) then

        m = - m
        c0 = gm * gc / ( ga * gb * ( 1.0e+00_real64 - x )**m )
        c1 = - ( - 1 )**m * gc / ( gam * gbm * rm )

        do k = 1, m - 1
          r0 = r0 * ( a - m + k - 1.0e+00_real64 ) * ( b - m + k - 1.0e+00_real64 ) &
            / ( k * ( k - m ) ) * ( 1.0e+00_real64 - x )
          f0 = f0 + r0
        end do

        do k = 1, m
          sp0 = sp0 + 1.0e+00_real64 / real ( k, kind = real64 )
        end do

        f1 = pa + pb - sp0 + 2.0e+00_real64 * el + log ( 1.0e+00_real64 - x )

        do k = 1, 250

          sp = sp + ( 1.0e+00_real64 - a ) &
            / ( k * ( a + k - 1.0e+00_real64 ) ) &
            + ( 1.0e+00_real64 - b ) / ( k * ( b + k - 1.0e+00_real64 ) )

          sm = 0.0e+00_real64
          do j = 1, m
            sm = sm + 1.0e+00_real64 / real ( j + k, kind = real64 )
          end do

          rp = pa + pb + 2.0e+00_real64 * el + sp - sm + log ( 1.0e+00_real64 - x )

          r1 = r1 * ( a + k - 1.0e+00_real64 ) * ( b + k - 1.0e+00_real64 ) &
            / ( k * ( m + k ) ) * ( 1.0e+00_real64 - x )

          f1 = f1 + r1 * rp

          if ( abs ( f1 - hw ) < abs ( f1 ) * eps ) then
            exit
          end if

          hw = f1

        end do

        hf = f0 * c0 + f1 * c1

      end if

    else

      call gamma ( a, ga )
      call gamma ( b, gb )
      call gamma ( c, gc )
      call gamma ( c - a, gca )
      call gamma ( c - b, gcb )
      call gamma ( c - a - b, gcab )
      call gamma ( a + b - c, gabc )
      c0 = gc * gcab / ( gca * gcb )
      c1 = gc * gabc / ( ga * gb ) * ( 1.0e+00_real64 - x )**( c - a - b )
      hf = 0.0e+00_real64
      r0 = c0
      r1 = c1

      do k = 1, 250

        r0 = r0 * ( a + k - 1.0e+00_real64 ) * ( b + k - 1.0e+00_real64 ) &
          / ( k * ( a + b - c + k ) ) * ( 1.0e+00_real64 - x )

        r1 = r1 * ( c - a + k - 1.0e+00_real64 ) * ( c - b + k - 1.0e+00_real64 ) &
          / ( k * ( c - a - b + k ) ) * ( 1.0e+00_real64 - x )

        hf = hf + r0 + r1

        if ( abs ( hf - hw ) < abs ( hf ) * eps ) then
          exit
        end if

        hw = hf

      end do

      hf = hf + c0 + c1

    end if

  else

    a0 = 1.0e+00_real64

    if ( a < c .and. c < 2.0e+00_real64 * a .and. b < c .and. c < 2.0e+00_real64 * b ) then

      a0 = ( 1.0e+00_real64 - x )**( c - a - b )
      a = c - a
      b = c - b

    end if

    hf = 1.0e+00_real64
    r = 1.0e+00_real64

    do k = 1, 250

      r = r * ( a + k - 1.0e+00_real64 ) * ( b + k - 1.0e+00_real64 ) &
        / ( k * ( c + k - 1.0e+00_real64 ) ) * x

      hf = hf + r

      if ( abs ( hf - hw ) <= abs ( hf ) * eps ) then
        exit
      end if

      hw = hf

    end do

    hf = a0 * hf

  end if

  if ( x1 < 0.0e+00_real64 ) then
    x = x1
    c0 = 1.0e+00_real64 / ( 1.0e+00_real64 - x )**aa
    hf = c0 * hf
  end if

  a = aa
  b = bb

  if ( 120 < k ) then
    write ( error_unit, '(a)' ) ' '
    write ( error_unit, '(a)' ) 'HYGFX - Warning!'
    write ( error_unit, '(a)' ) '  A large number of iterations were needed.'
    write ( error_unit, '(a)' ) '  The accuracy of the results should be checked.'
  end if

  return
end subroutine hygfx
!> @brief subroutine hygfz.
!> @return None.
!>
!> @param a [inout] Argument a.
!> @param b [inout] Argument b.
!> @param c [inout] Argument c.
!> @param z [in] Argument z.
!> @param zhf [inout] Argument zhf.
subroutine hygfz ( a, b, c, z, zhf )

!*****************************************************************************80
!
!! HYGFZ computes the hypergeometric function F(a,b,c,x) for complex argument.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    03 August 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) A, B, C, parameters.
!
!    Input, complex ( kind = real64 ) Z, the argument.
!
!    Output, complex ( kind = real64 ) ZHF, the value of F(a,b,c,z).
!
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  implicit none

  real(real64), intent(inout) :: a
  real(real64) a0
  real(real64) aa
  real(real64), intent(inout) :: b
  real(real64) bb
  real(real64), intent(inout) :: c
  real(real64) ca
  real(real64) cb
  real(real64) el
  real(real64) eps
  real(real64) g0
  real(real64) g1
  real(real64) g2
  real(real64) g3
  real(real64) ga
  real(real64) gab
  real(real64) gabc
  real(real64) gam
  real(real64) gb
  real(real64) gba
  real(real64) gbm
  real(real64) gc
  real(real64) gca
  real(real64) gcab
  real(real64) gcb
  real(real64) gcbk
  real(real64) gm
  integer(int32) j
  integer(int32) k
  logical l0
  logical l1
  logical l2
  logical l3
  logical l4
  logical l5
  logical l6
  integer(int32) m
  integer(int32) mab
  integer(int32) mcab
  integer(int32) nca
  integer(int32) ncb
  integer(int32) nm
  real(real64) pa
  real(real64) pac
  real(real64) pb
  real(real64) pca
  real(real64) pi
  real(real64) rk1
  real(real64) rk2
  real(real64) rm
  real(real64) sj1
  real(real64) sj2
  real(real64) sm
  real(real64) sp
  real(real64) sp0
  real(real64) sq
  real(real64) t0
  real(real64) w0
  real(real64) ws
  real(real64) x
  real(real64) y
  complex(real64), intent(in) :: z
  complex(real64) z00
  complex(real64) z1
  complex(real64) zc0
  complex(real64) zc1
  complex(real64) zf0
  complex(real64) zf1
  complex(real64), intent(inout) :: zhf
  complex(real64) zp
  complex(real64) zp0
  complex(real64) zr
  complex(real64) zr0
  complex(real64) zr1
  complex(real64) zw

  x = real ( z, kind = real64 )
  y = imag ( z )
  eps = 1.0e-15_real64
  l0 = c == int ( c ) .and. c < 0.0e+00_real64
  l1 = abs ( 1.0e+00_real64 - x ) < eps .and. y == 0.0e+00_real64 .and. &
    c - a - b <= 0.0e+00_real64
  l2 = abs ( z + 1.0e+00_real64 ) < eps .and. &
    abs ( c - a + b - 1.0e+00_real64 ) < eps
  l3 = a == int ( a ) .and. a < 0.0e+00_real64
  l4 = b == int ( b ) .and. b < 0.0e+00_real64
  l5 = c - a == int ( c - a ) .and. c - a <= 0.0e+00_real64
  l6 = c - b == int ( c - b ) .and. c - b <= 0.0e+00_real64
  aa = a
  bb = b
  a0 = abs ( z )
  if ( 0.95e+00_real64 < a0 ) then
    eps = 1.0e-08_real64
  end if
  pi = 3.141592653589793e+00_real64
  el = 0.5772156649015329e+00_real64
  zw = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )

  if ( l0 .or. l1 ) then
    write ( error_unit, '(a)' ) ' '
    write ( error_unit, '(a)' ) 'HYGFZ - Fatal error!'
    write ( error_unit, '(a)' ) '  The hypergeometric series is divergent.'
    zhf = cmplx ( ieee_value ( 0.0e+00_real64, ieee_quiet_nan ), &
      ieee_value ( 0.0e+00_real64, ieee_quiet_nan ), kind = real64 )
    return
  end if

  if ( a0 == 0.0e+00_real64 .or. a == 0.0e+00_real64 .or. b == 0.0e+00_real64 ) then

    zhf = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )

  else if ( z == 1.0e+00_real64.and. 0.0e+00_real64 < c - a - b ) then

    call gamma ( c, gc )
    call gamma ( c - a - b, gcab )
    call gamma ( c - a, gca )
    call gamma ( c - b, gcb )
    zhf = gc * gcab / ( gca * gcb )

  else if ( l2 ) then

    g0 = sqrt ( pi ) * 2.0e+00_real64 ** ( - a )
    call gamma ( c, g1 )
    call gamma ( 1.0e+00_real64 + a / 2.0e+00_real64 - b, g2 )
    call gamma ( 0.5e+00_real64 + 0.5e+00_real64 * a, g3 )
    zhf = g0 * g1 / ( g2 * g3 )

  else if ( l3 .or. l4 ) then

    if ( l3 ) then
      nm = int ( abs ( a ) )
    end if

    if ( l4 ) then
      nm = int ( abs ( b ) )
    end if

    zhf = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    zr = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    do k = 1, nm
      zr = zr * ( a + k - 1.0e+00_real64 ) * ( b + k - 1.0e+00_real64 ) &
        / ( k * ( c + k - 1.0e+00_real64 ) ) * z
      zhf = zhf + zr
    end do

  else if ( l5 .or. l6 ) then

    if ( l5 ) then
      nm = int ( abs ( c - a ) )
    end if

    if ( l6 ) then
      nm = int ( abs ( c - b ) )
    end if

    zhf = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    zr = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
    do k = 1, nm
      zr = zr * ( c - a + k - 1.0e+00_real64 ) * ( c - b + k - 1.0e+00_real64 ) &
        / ( k * ( c + k - 1.0e+00_real64 ) ) * z
      zhf = zhf + zr
    end do
    zhf = ( 1.0e+00_real64 - z ) ** ( c - a - b ) * zhf

  else if ( a0 <= 1.0e+00_real64 ) then

    if ( x < 0.0e+00_real64 ) then

      z1 = z / ( z - 1.0e+00_real64 )
      if ( a < c .and. b < a .and. 0.0e+00_real64 < b ) then  
        a = bb
        b = aa
      end if
      zc0 = 1.0e+00_real64 / ( ( 1.0e+00_real64 - z ) ** a )
      zhf = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
      zr0 = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
      zw = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
      do k = 1, 500
        zr0 = zr0 * ( a + k - 1.0e+00_real64 ) * ( c - b + k - 1.0e+00_real64 ) &
          / ( k * ( c + k - 1.0e+00_real64 ) ) * z1
        zhf = zhf + zr0
        if ( abs ( zhf - zw ) < abs ( zhf ) * eps ) then
          exit
        end if
        zw = zhf
      end do

      zhf = zc0 * zhf

    else if ( 0.90e+00_real64 <= a0 ) then

      gm = 0.0e+00_real64
      mcab = int ( c - a - b + eps * sign ( 1.0e+00_real64, c - a - b ) )

      if ( abs ( c - a - b - mcab ) < eps ) then

        m = int ( c - a - b )
        call gamma ( a, ga )
        call gamma ( b, gb )
        call gamma ( c, gc )
        call gamma ( a + m, gam )
        call gamma ( b + m, gbm ) 
        call psi ( a, pa )
        call psi ( b, pb )
        if ( m /= 0 ) then
          gm = 1.0e+00_real64
        end if
        do j = 1, abs ( m ) - 1
          gm = gm * j
        end do
        rm = 1.0e+00_real64
        do j = 1, abs ( m )
          rm = rm * j
        end do
        zf0 = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
        zr0 = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
        zr1 = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
        sp0 = 0.0e+00_real64
        sp = 0.0e+00_real64

        if ( 0 <= m ) then

          zc0 = gm * gc / ( gam * gbm )
          zc1 = - gc * ( z - 1.0e+00_real64 ) ** m / ( ga * gb * rm )
          do k = 1, m - 1
            zr0 = zr0 * ( a + k - 1.0e+00_real64 ) * ( b + k - 1.0e+00_real64 ) &
              / ( k * ( k - m ) ) * ( 1.0e+00_real64 - z )
            zf0 = zf0 + zr0
          end do
          do k = 1, m
            sp0 = sp0 + 1.0e+00_real64 / ( a + k - 1.0e+00_real64 ) &
              + 1.0e+00_real64 / ( b + k - 1.0e+00_real64 ) - 1.0e+00_real64 / k
          end do
          zf1 = pa + pb + sp0 + 2.0e+00_real64 * el + log ( 1.0e+00_real64 - z )
          zw = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
          do k = 1, 500
            sp = sp + ( 1.0e+00_real64 - a ) &
              / ( k * ( a + k - 1.0e+00_real64 ) ) + ( 1.0e+00_real64 - b ) &
              / ( k * ( b + k - 1.0e+00_real64 ) )
            sm = 0.0e+00_real64
            do j = 1, m
              sm = sm + ( 1.0e+00_real64 - a ) / ( ( j + k ) &
                * ( a + j + k - 1.0e+00_real64 ) ) &
                + 1.0e+00_real64 / ( b + j + k - 1.0e+00_real64 )
            end do
            zp = pa + pb + 2.0e+00_real64 * el + sp + sm  + log ( 1.0e+00_real64 - z )
            zr1 = zr1 * ( a + m + k - 1.0e+00_real64 ) &
              * ( b + m + k - 1.0e+00_real64 ) / ( k * ( m + k ) ) &
              * ( 1.0e+00_real64 - z )
            zf1 = zf1 + zr1 * zp
            if ( abs ( zf1 - zw ) < abs ( zf1 ) * eps ) then
              exit
            end if
            zw = zf1
          end do

          zhf = zf0 * zc0 + zf1 * zc1

        else if ( m < 0 ) then

          m = - m
          zc0 = gm * gc / ( ga * gb * ( 1.0e+00_real64 - z ) ** m )
          zc1 = - ( - 1.0e+00_real64 ) ** m * gc / ( gam * gbm * rm )
          do k = 1, m - 1
            zr0 = zr0 * ( a - m + k - 1.0e+00_real64 ) &
              * ( b - m + k - 1.0e+00_real64 ) / ( k * ( k - m ) ) &
              * ( 1.0e+00_real64 - z )
            zf0 = zf0 + zr0
          end do

          do k = 1, m
            sp0 = sp0 + 1.0e+00_real64 / k
          end do

          zf1 = pa + pb - sp0 + 2.0e+00_real64 * el + log ( 1.0e+00_real64 - z )

          zw = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
          do k = 1, 500
            sp = sp + ( 1.0e+00_real64 - a ) / ( k * ( a + k - 1.0e+00_real64 ) ) &
              + ( 1.0e+00_real64 - b ) / ( k * ( b + k - 1.0e+00_real64 ) )
            sm = 0.0e+00_real64
            do j = 1, m
              sm = sm + 1.0e+00_real64 / ( j + k )
            end do
            zp = pa + pb + 2.0e+00_real64 * el + sp - sm + log ( 1.0e+00_real64 - z )
            zr1 = zr1 * ( a + k - 1.0e+00_real64 ) * ( b + k - 1.0e+00_real64 ) &
              / ( k * ( m + k ) ) * ( 1.0e+00_real64 - z )
            zf1 = zf1 + zr1 * zp
            if ( abs ( zf1 - zw ) < abs ( zf1 ) * eps ) then
              exit
            end if
            zw = zf1

          end do

          zhf = zf0 * zc0 + zf1 * zc1

        end if

      else

        call gamma ( a, ga )
        call gamma ( b, gb )
        call gamma ( c, gc )
        call gamma ( c - a, gca )
        call gamma ( c - b, gcb )
        call gamma ( c - a - b, gcab )
        call gamma ( a + b - c, gabc )
        zc0 = gc * gcab / ( gca * gcb )
        zc1 = gc * gabc / ( ga * gb ) * ( 1.0e+00_real64 - z ) ** ( c - a - b )
        zhf = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
        zr0 = zc0
        zr1 = zc1
        zw = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
        do k = 1, 500
          zr0 = zr0 * ( a + k - 1.0e+00_real64 ) * ( b + k - 1.0e+00_real64 ) &
            / ( k * ( a + b - c + k ) ) * ( 1.0e+00_real64 - z )
          zr1 = zr1 * ( c - a + k - 1.0e+00_real64 ) &
            * ( c - b + k - 1.0e+00_real64 ) / ( k * ( c - a - b + k ) ) &
            * ( 1.0e+00_real64 - z )
          zhf = zhf + zr0 + zr1
          if ( abs ( zhf - zw ) < abs ( zhf ) * eps ) then
            exit
          end if
          zw = zhf
        end do

        zhf = zhf + zc0 + zc1

      end if

    else

      z00 = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )

      if ( c - a < a .and. c - b < b ) then
        z00 = ( 1.0e+00_real64 - z ) ** ( c - a - b )
        a = c - a
        b = c - b
      end if

      zhf = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
      zr = cmplx ( 1.0e+00_real64, 0.0e+00_real64, kind = real64 )
      zw = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )

      do k = 1, 1500
        zr = zr * ( a + k - 1.0e+00_real64 ) * ( b + k - 1.0e+00_real64 ) &
          / ( k * ( c + k - 1.0e+00_real64 ) ) * z
        zhf = zhf + zr
        if ( abs ( zhf - zw ) <= abs ( zhf ) * eps ) then
          exit
        end if
        zw = zhf
      end do

      zhf = z00 * zhf

    end if

  else if ( 1.0e+00_real64 < a0 ) then

    mab = int ( a - b + eps * sign ( 1.0e+00_real64, a - b ) )

    if ( abs ( a - b - mab ) < eps .and. a0 <= 1.1e+00_real64 ) then
      b = b + eps
    end if

    if ( eps < abs ( a - b - mab ) ) then

      call gamma ( a, ga )
      call gamma ( b, gb )
      call gamma ( c, gc )
      call gamma ( a - b, gab )
      call gamma ( b - a, gba )
      call gamma ( c - a, gca )
      call gamma ( c - b, gcb )
      zc0 = gc * gba / ( gca * gb * ( - z ) ** a )
      zc1 = gc * gab / ( gcb * ga * ( - z ) ** b )
      zr0 = zc0
      zr1 = zc1
      zhf = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
      zw = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )

      do k = 1, 500
        zr0 = zr0 * ( a + k - 1.0e+00_real64 ) * ( a - c + k ) &
          / ( ( a - b + k ) * k * z )
        zr1 = zr1 * ( b + k - 1.0e+00_real64 ) * ( b - c + k ) &
          / ( ( b - a + k ) * k * z )
        zhf = zhf + zr0 + zr1
        if ( abs ( ( zhf - zw ) / zhf ) <= eps ) then
          exit
        end if
        zw = zhf
      end do

      zhf = zhf + zc0 + zc1

    else

      if ( a - b < 0.0e+00_real64 ) then
        a = bb
        b = aa
      end if

      ca = c - a
      cb = c - b
      nca = int ( ca + eps * sign ( 1.0e+00_real64, ca ) )
      ncb = int ( cb + eps * sign ( 1.0e+00_real64, cb ) )

      if ( abs ( ca - nca ) < eps .or. abs ( cb - ncb ) < eps ) then
        c = c + eps
      end if

      call gamma ( a, ga )
      call gamma ( c, gc )
      call gamma ( c - b, gcb )
      call psi ( a, pa )
      call psi ( c - a, pca )
      call psi ( a - c, pac )
      mab = int ( a - b + eps )
      zc0 = gc / ( ga * ( - z ) ** b )
      call gamma ( a - b, gm )
      zf0 = gm / gcb * zc0
      zr = zc0
      do k = 1, mab - 1
        zr = zr * ( b + k - 1.0e+00_real64 ) / ( k * z )
        t0 = a - b - k
        call gamma ( t0, g0 )
        call gamma ( c - b - k, gcbk )
        zf0 = zf0 + zr * g0 / gcbk
      end do

      if ( mab == 0 ) then
        zf0 = cmplx ( 0.0e+00_real64, 0.0e+00_real64, kind = real64 )
      end if

      zc1 = gc / ( ga * gcb * ( - z ) ** a )
      sp = -2.0e+00_real64 * el - pa - pca
      do j = 1, mab
        sp = sp + 1.0e+00_real64 / j
      end do
      zp0 = sp + log ( - z )
      sq = 1.0e+00_real64
      do j = 1, mab
        sq = sq * ( b + j - 1.0e+00_real64 ) * ( b - c + j ) / j
      end do
      zf1 = ( sq * zp0 ) * zc1
      zr = zc1
      rk1 = 1.0e+00_real64
      sj1 = 0.0e+00_real64

      do k = 1, 10000
        zr = zr / z
        rk1 = rk1 * ( b + k - 1.0e+00_real64 ) * ( b - c + k ) / ( k * k )
        rk2 = rk1
        do j = k + 1, k + mab
          rk2 = rk2 * ( b + j - 1.0e+00_real64 ) * ( b - c + j ) / j
        end do
        sj1 = sj1 + ( a - 1.0e+00_real64 ) / ( k * ( a + k - 1.0e+00_real64 ) ) &
          + ( a - c - 1.0e+00_real64 ) / ( k * ( a - c + k - 1.0e+00_real64 ) )
        sj2 = sj1
        do j = k + 1, k + mab
          sj2 = sj2 + 1.0e+00_real64 / j
        end do
        zp = -2.0e+00_real64 * el - pa - pac + sj2 - 1.0e+00_real64 / ( k + a - c ) &
          - pi / tan ( pi * ( k + a - c ) ) + log ( - z ) 
        zf1 = zf1 + rk2 * zr * zp
        ws = abs ( zf1 )
        if ( abs ( ( ws - w0 ) / ws ) < eps ) then
          exit
        end if
        w0 = ws
      end do

      zhf = zf0 + zf1

    end if

  end if

  a = aa
  b = bb
  if ( 150 < k ) then
    write ( error_unit, '(a)' ) ' '
    write ( error_unit, '(a)' ) 'HYGFZ - Warning!'
    write ( error_unit, '(a)' ) '  The solution returned may have low accuracy.'
  end if

  return
end subroutine hygfz
!> @brief subroutine ik01a.
!> @return None.
!>
!> @param x [in] Argument x.
!> @param bi0 [inout] Argument bi0.
!> @param di0 [inout] Argument di0.
!> @param bi1 [inout] Argument bi1.
!> @param di1 [inout] Argument di1.
!> @param bk0 [inout] Argument bk0.
!> @param dk0 [inout] Argument dk0.
!> @param bk1 [inout] Argument bk1.
!> @param dk1 [inout] Argument dk1.
pure subroutine ik01a ( x, bi0, di0, bi1, di1, bk0, dk0, bk1, dk1 )

!*****************************************************************************80
!
!! IK01A compute Bessel function I0(x), I1(x), K0(x), and K1(x).
!
!  Discussion:
!
!    This procedure computes modified Bessel functions I0(x), I1(x),
!    K0(x) and K1(x), and their derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    16 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) X, the argument.
!
!    Output, real(real64) BI0, DI0, BI1, DI1, BK0, DK0, BK1, DK1, the
!    values of I0(x), I0'(x), I1(x), I1'(x), K0(x), K0'(x), K1(x), K1'(x).
!
  implicit none

  real(real64), parameter, dimension ( 12 ) :: a = [&
    0.125e+00_real64, 7.03125e-02_real64, &
    7.32421875e-02_real64, 1.1215209960938e-01_real64, &
    2.2710800170898e-01_real64, 5.7250142097473e-01_real64, &
    1.7277275025845e+00_real64, 6.0740420012735e+00_real64, &
    2.4380529699556e+01_real64, 1.1001714026925e+02_real64, &
    5.5133589612202e+02_real64, 3.0380905109224e+03_real64]
  real(real64), parameter, dimension ( 8 ) :: a1 = [&
    0.125e+00_real64, 0.2109375e+00_real64, &
    1.0986328125e+00_real64, 1.1775970458984e+01_real64, &
    2.1461706161499e+02_real64, 5.9511522710323e+03_real64, &
    2.3347645606175e+05_real64, 1.2312234987631e+07_real64]
  real(real64), parameter, dimension ( 12 ) :: b = [&
    -0.375e+00_real64, -1.171875e-01_real64, &
    -1.025390625e-01_real64, -1.4419555664063e-01_real64, &
    -2.7757644653320e-01_real64, -6.7659258842468e-01_real64, &
    -1.9935317337513e+00_real64, -6.8839142681099e+00_real64, &
    -2.7248827311269e+01_real64, -1.2159789187654e+02_real64, &
    -6.0384407670507e+02_real64, -3.3022722944809e+03_real64]
  real(real64), intent(inout) :: bi0
  real(real64), intent(inout) :: bi1
  real(real64), intent(inout) :: bk0
  real(real64), intent(inout) :: bk1
  real(real64) ca
  real(real64) cb
  real(real64) ct
  real(real64), intent(inout) :: di0
  real(real64), intent(inout) :: di1
  real(real64), intent(inout) :: dk0
  real(real64), intent(inout) :: dk1
  real(real64) el
  integer(int32) k
  integer(int32) k0
  real(real64) pi
  real(real64) r
  real(real64) w0
  real(real64) ww
  real(real64), intent(in) :: x
  real(real64) x2
  real(real64) xr
  real(real64) xr2

  pi = 3.141592653589793e+00_real64
  el = 0.5772156649015329e+00_real64
  x2 = x * x

  if ( x == 0.0e+00_real64 ) then

    bi0 = 1.0e+00_real64
    bi1 = 0.0e+00_real64
    bk0 = 1.0e+300_real64
    bk1 = 1.0e+300_real64
    di0 = 0.0e+00_real64
    di1 = 0.5e+00_real64
    dk0 = -1.0e+300_real64
    dk1 = -1.0e+300_real64
    return

  else if ( x <= 18.0e+00_real64 ) then

    bi0 = 1.0e+00_real64
    r = 1.0e+00_real64
    do k = 1, 50
      r = 0.25e+00_real64 * r * x2 / ( k * k )
      bi0 = bi0 + r
      if ( abs ( r / bi0 ) < 1.0e-15_real64 ) then
        exit
      end if
    end do

    bi1 = 1.0e+00_real64
    r = 1.0e+00_real64
    do k = 1, 50
      r = 0.25e+00_real64 * r * x2 / ( k * ( k + 1 ) )
      bi1 = bi1 + r
      if ( abs ( r / bi1 ) < 1.0e-15_real64 ) then
        exit
      end if
    end do

    bi1 = 0.5e+00_real64 * x * bi1

  else

    if ( x < 35.0e+00_real64 ) then
      k0 = 12
    else if ( x < 50.0e+00_real64 ) then
      k0 = 9
    else
      k0 = 7
    end if

    ca = exp ( x ) / sqrt ( 2.0e+00_real64 * pi * x )
    bi0 = 1.0e+00_real64
    xr = 1.0e+00_real64 / x
    do k = 1, k0
      bi0 = bi0 + a(k) * xr ** k
    end do
    bi0 = ca * bi0
    bi1 = 1.0e+00_real64
    do k = 1, k0
      bi1 = bi1 + b(k) * xr ** k
    end do
    bi1 = ca * bi1

  end if

  if ( x <= 9.0e+00_real64 ) then

    ct = - ( log ( x / 2.0e+00_real64 ) + el )
    bk0 = 0.0e+00_real64
    w0 = 0.0e+00_real64
    ww = 0.0e+00_real64
    r = 1.0e+00_real64
    do k = 1, 50
      w0 = w0 + 1.0e+00_real64 / k
      r = 0.25e+00_real64 * r / ( k * k ) * x2
      bk0 = bk0 + r * ( w0 + ct )
      if ( abs ( ( bk0 - ww ) / bk0 ) < 1.0e-15_real64 ) then
        exit
      end if
      ww = bk0
    end do

    bk0 = bk0 + ct

  else

    cb = 0.5e+00_real64 / x
    xr2 = 1.0e+00_real64 / x2
    bk0 = 1.0e+00_real64
    do k = 1, 8
      bk0 = bk0 + a1(k) * xr2 ** k
    end do
    bk0 = cb * bk0 / bi0

  end if

  bk1 = ( 1.0e+00_real64 / x - bi1 * bk0 ) / bi0
  di0 = bi1
  di1 = bi0 - bi1 / x
  dk0 = - bk1
  dk1 = - bk0 - bk1 / x

  return
end subroutine ik01a
!> @brief subroutine ik01b.
!> @return None.
!>
!> @param x [in] Argument x.
!> @param bi0 [inout] Argument bi0.
!> @param di0 [inout] Argument di0.
!> @param bi1 [inout] Argument bi1.
!> @param di1 [inout] Argument di1.
!> @param bk0 [inout] Argument bk0.
!> @param dk0 [inout] Argument dk0.
!> @param bk1 [inout] Argument bk1.
!> @param dk1 [inout] Argument dk1.
pure subroutine ik01b ( x, bi0, di0, bi1, di1, bk0, dk0, bk1, dk1 )

!*****************************************************************************80
!
!! IK01B: Bessel functions I0(x), I1(x), K0(x), and K1(x) and derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    17 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) X, the argument.
!
!    Output, real(real64) BI0, DI0, BI1, DI1, BK0, DK0, BK1, DK1, the
!    values of I0(x), I0'(x), I1(x), I1'(x), K0(x), K0'(x), K1(x), K1'(x).
!
  implicit none

  real(real64), intent(inout) :: bi0
  real(real64), intent(inout) :: bi1
  real(real64), intent(inout) :: bk0
  real(real64), intent(inout) :: bk1
  real(real64), intent(inout) :: di0
  real(real64), intent(inout) :: di1
  real(real64), intent(inout) :: dk0
  real(real64), intent(inout) :: dk1
  real(real64) t
  real(real64) t2
  real(real64), intent(in) :: x

  if ( x == 0.0e+00_real64 ) then

    bi0 = 1.0e+00_real64
    bi1 = 0.0e+00_real64
    bk0 = 1.0e+300_real64
    bk1 = 1.0e+300_real64
    di0 = 0.0e+00_real64
    di1 = 0.5e+00_real64
    dk0 = -1.0e+300_real64
    dk1 = -1.0e+300_real64
    return

  else if ( x <= 3.75e+00_real64 ) then

    t = x / 3.75e+00_real64
    t2 = t * t

    bi0 = ((((( &
        0.0045813e+00_real64   * t2 &
      + 0.0360768e+00_real64 ) * t2 &
      + 0.2659732e+00_real64 ) * t2 &
      + 1.2067492e+00_real64 ) * t2 &
      + 3.0899424e+00_real64 ) * t2 &
      + 3.5156229e+00_real64 ) * t2 &
      + 1.0e+00_real64

    bi1 = x * (((((( &
        0.00032411e+00_real64   * t2 &
      + 0.00301532e+00_real64 ) * t2 &
      + 0.02658733e+00_real64 ) * t2 &
      + 0.15084934e+00_real64 ) * t2 &
      + 0.51498869e+00_real64 ) * t2 &
      + 0.87890594e+00_real64 ) * t2 &
      + 0.5e+00_real64 )

  else

    t = 3.75e+00_real64 / x

    bi0 = (((((((( &
        0.00392377e+00_real64   * t &
      - 0.01647633e+00_real64 ) * t &
      + 0.02635537e+00_real64 ) * t &
      - 0.02057706e+00_real64 ) * t &
      + 0.916281e-02_real64 ) * t &
      - 0.157565e-02_real64 ) * t &
      + 0.225319e-02_real64 ) * t &
      + 0.01328592e+00_real64 ) * t &
      + 0.39894228e+00_real64 ) * exp ( x ) / sqrt ( x )

    bi1 = (((((((( &
      - 0.420059e-02_real64     * t &
      + 0.01787654e+00_real64 ) * t &
      - 0.02895312e+00_real64 ) * t &
      + 0.02282967e+00_real64 ) * t &
      - 0.01031555e+00_real64 ) * t &
      + 0.163801e-02_real64 ) * t &
      - 0.00362018e+00_real64 ) * t &
      - 0.03988024e+00_real64 ) * t &
      + 0.39894228e+00_real64 ) * exp ( x ) / sqrt ( x )

  end if

  if ( x <= 2.0e+00_real64 ) then

    t = x / 2.0e+00_real64
    t2 = t * t

    bk0 = ((((( &
        0.0000074e+00_real64   * t2 &
      + 0.0001075e+00_real64 ) * t2 &
      + 0.00262698e+00_real64 ) * t2 &
      + 0.0348859e+00_real64 ) * t2 &
      + 0.23069756e+00_real64 ) * t2 &
      + 0.4227842e+00_real64 ) * t2 &
      - 0.57721566e+00_real64 - bi0 * log ( t )

    bk1 = (((((( &
      - 0.00004686e+00_real64   * t2 &
      - 0.00110404e+00_real64 ) * t2 &
      - 0.01919402e+00_real64 ) * t2 &
      - 0.18156897e+00_real64 ) * t2 &
      - 0.67278579e+00_real64 ) * t2 &
      + 0.15443144e+00_real64 ) * t2 &
      + 1.0e+00_real64 ) / x + bi1 * log ( t )

  else

    t = 2.0e+00_real64 / x
    t2 = t * t

    bk0 = (((((( &
        0.00053208e+00_real64   * t &
      - 0.0025154e+00_real64 )  * t &
      + 0.00587872e+00_real64 ) * t &
      - 0.01062446e+00_real64 ) * t &
      + 0.02189568e+00_real64 ) * t &
      - 0.07832358e+00_real64 ) * t &
      + 1.25331414e+00_real64 ) * exp ( - x ) / sqrt ( x )

    bk1 = (((((( &
      - 0.00068245e+00_real64   * t &
      + 0.00325614e+00_real64 ) * t &
      - 0.00780353e+00_real64 ) * t &
      + 0.01504268e+00_real64 ) * t &
      - 0.0365562e+00_real64  ) * t & 
      + 0.23498619e+00_real64 ) * t &
      + 1.25331414e+00_real64 ) * exp ( - x ) / sqrt ( x )

  end if

  di0 = bi1
  di1 = bi0 - bi1 / x
  dk0 = -bk1
  dk1 = -bk0 - bk1 / x

  return
end subroutine ik01b
!> @brief subroutine ikna.
!> @return None.
!>
!> @param n [in] Argument n.
!> @param x [in] Argument x.
!> @param nm [inout] Argument nm.
!> @param bi [inout] Argument bi.
!> @param di [inout] Argument di.
!> @param bk [inout] Argument bk.
!> @param dk [inout] Argument dk.
subroutine ikna ( n, x, nm, bi, di, bk, dk )

!*****************************************************************************80
!
!! IKNA compute Bessel function In(x) and Kn(x), and derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    16 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) N, the order of In(x) and Kn(x).
!
!    Input, real(real64) X, the argument.
!
!    Output, integer(int32) NM, the highest order computed.
!
!    Output, real(real64) BI(0:N), DI(0:N), BK(0:N), DK(0:N),
!    the values of In(x), In'(x), Kn(x), Kn'(x).
!
  implicit none

  integer(int32), intent(in) :: n

  real(real64), intent(inout) :: bi(0:n)
  real(real64) bi0
  real(real64) bi1
  real(real64), intent(inout) :: bk(0:n)
  real(real64) bk0
  real(real64) bk1
  real(real64), intent(inout) :: di(0:n)
  real(real64) di0
  real(real64) di1
  real(real64), intent(inout) :: dk(0:n)
  real(real64) dk0
  real(real64) dk1
  real(real64) f
  real(real64) f0
  real(real64) f1
  real(real64) g
  real(real64) g0
  real(real64) g1
  real(real64) h
  real(real64) h0
  real(real64) h1
  integer(int32) k
  integer(int32) m
  ! integer(int32) msta1
  ! integer(int32) msta2
  integer(int32), intent(inout) :: nm
  real(real64) s0
  real(real64), intent(in) :: x

  nm = n

  if ( x <= 1.0e-100_real64 ) then
    do k = 0, n
      bi(k) = 0.0e+00_real64
      di(k) = 0.0e+00_real64
      bk(k) = 1.0e+300_real64
      dk(k) = -1.0e+300_real64
    end do
    bi(0) = 1.0e+00_real64
    di(1) = 0.5e+00_real64
    return
  end if

  call ik01a ( x, bi0, di0, bi1, di1, bk0, dk0, bk1, dk1 )
  bi(0) = bi0
  bi(1) = bi1
  bk(0) = bk0
  bk(1) = bk1
  di(0) = di0
  di(1) = di1
  dk(0) = dk0
  dk(1) = dk1

  if ( n <= 1 ) then
    return
  end if

  if ( 40.0e+00_real64 < x .and. n < int ( 0.25e+00_real64 * x ) ) then

    h0 = bi0
    h1 = bi1
    do k = 2, n
      h = -2.0e+00_real64 * ( k - 1.0e+00_real64 ) / x * h1 + h0
      bi(k) = h
      h0 = h1
      h1 = h
    end do

  else

    m = msta1 ( x, 200 )

    if ( m < n ) then
      nm = m
    else
      m = msta2 ( x, n, 15 )
    end if

    f0 = 0.0e+00_real64
    f1 = 1.0e-100_real64
    do k = m, 0, -1
      f = 2.0e+00_real64 * ( k + 1.0e+00_real64 ) * f1 / x + f0
      if ( k <= nm ) then
        bi(k) = f
      end if
      f0 = f1
      f1 = f
    end do
    s0 = bi0 / f
    do k = 0, nm
      bi(k) = s0 * bi(k)
    end do
  end if

  g0 = bk0
  g1 = bk1
  do k = 2, nm
    g = 2.0e+00_real64 * ( k - 1.0e+00_real64 ) / x * g1 + g0
    bk(k) = g
    g0 = g1
    g1 = g
  end do

  do k = 2, nm
    di(k) = bi(k-1) - k / x * bi(k)
    dk(k) = - bk(k-1) - k / x * bk(k)
  end do

  return
end subroutine ikna
!> @brief subroutine iknb.
!> @return None.
!>
!> @param n [in] Argument n.
!> @param x [in] Argument x.
!> @param nm [inout] Argument nm.
!> @param bi [inout] Argument bi.
!> @param di [inout] Argument di.
!> @param bk [inout] Argument bk.
!> @param dk [inout] Argument dk.
subroutine iknb ( n, x, nm, bi, di, bk, dk )

!*****************************************************************************80
!
!! IKNB compute Bessel function In(x) and Kn(x).
!
!  Discussion:
!
!    Compute modified Bessel functions In(x) and Kn(x),
!    and their derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    17 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) N, the order of In(x) and Kn(x).
!
!    Input, real(real64) X, the argument.
!
!    Output, integer(int32) NM, the highest order computed.
!
!    Output, real(real64) BI(0:N), DI(0:N), BK(0:N), DK(0:N),
!    the values of In(x), In'(x), Kn(x), Kn'(x).
!
  implicit none

  integer(int32), intent(in) :: n

  real(real64) a0
  real(real64), intent(inout) :: bi(0:n)
  real(real64), intent(inout) :: bk(0:n)
  real(real64) bkl
  real(real64) bs
  real(real64), intent(inout) :: di(0:n)
  real(real64), intent(inout) :: dk(0:n)
  real(real64) el
  real(real64) f
  real(real64) f0
  real(real64) f1
  real(real64) g
  real(real64) g0
  real(real64) g1
  integer(int32) k
  integer(int32) k0
  integer(int32) l
  integer(int32) m
  ! integer(int32) msta1
  ! integer(int32) msta2
  integer(int32), intent(inout) :: nm
  real(real64) pi
  real(real64) r
  real(real64) s0
  real(real64) sk0
  real(real64) vt
  real(real64), intent(in) :: x

  pi = 3.141592653589793e+00_real64
  el = 0.5772156649015329e0_real64
  nm = n

  if ( x <= 1.0e-100_real64 ) then
    do k = 0, n
      bi(k) = 0.0e+00_real64
      di(k) = 0.0e+00_real64
      bk(k) = 1.0e+300_real64
      dk(k) = -1.0e+300_real64
    end do
    bi(0) = 1.0e+00_real64
    di(1) = 0.5e+00_real64
    return
  end if

  if ( n == 0 ) then
    nm = 1
  end if

  m = msta1 ( x, 200 )
  if ( m < nm ) then
    nm = m
  else
    m = msta2 ( x, nm, 15 )
  end if

  bs = 0.0e+00_real64
  sk0 = 0.0e+00_real64
  f0 = 0.0e+00_real64
  f1 = 1.0e-100_real64
  do k = m, 0, -1
    f = 2.0e+00_real64 * ( k + 1.0e+00_real64 ) / x * f1 + f0
    if ( k <= nm ) then
      bi(k) = f
    end if
    if ( k /= 0 .and. k == 2 * int ( k / 2 ) ) then
      sk0 = sk0 + 4.0e+00_real64 * f / k
    end if
    bs = bs + 2.0e+00_real64 * f
    f0 = f1
    f1 = f
  end do

  s0 = exp ( x ) / ( bs - f )
  do k = 0, nm
    bi(k) = s0 * bi(k)
  end do

  if ( x <= 8.0e+00_real64 ) then
    bk(0) = - ( log ( 0.5e+00_real64 * x ) + el ) * bi(0) + s0 * sk0
    bk(1) = ( 1.0e+00_real64 / x - bi(1) * bk(0) ) / bi(0)
  else
    a0 = sqrt ( pi / ( 2.0e+00_real64 * x ) ) * exp ( - x ) 

    if ( x < 25.0e+00_real64 ) then
      k0 = 16
    else if ( x < 80.0e+00_real64 ) then
      k0 = 10
    else if ( x < 200.0e+00_real64 ) then
      k0 = 8
    else
      k0 = 6
    end if

    do l = 0, 1
      bkl = 1.0e+00_real64
      vt = 4.0e+00_real64 * l
      r = 1.0e+00_real64
      do k = 1, k0
        r = 0.125e+00_real64 * r * ( vt - ( 2.0e+00_real64 * k - 1.0e+00_real64 ) ** 2 ) / ( k * x )
        bkl = bkl + r
      end do
      bk(l) = a0 * bkl
    end do
  end if

  g0 = bk(0)
  g1 = bk(1)
  do k = 2, nm
    g = 2.0e+00_real64 * ( k - 1.0e+00_real64 ) / x * g1 + g0
    bk(k) = g
    g0 = g1
    g1 = g
  end do

  di(0) = bi(1)
  dk(0) = -bk(1)
  do k = 1, nm
    di(k) = bi(k-1) - k / x * bi(k)
    dk(k) = -bk(k-1) - k / x * bk(k)
  end do

  return
end subroutine iknb
!> @brief subroutine ikv.
!> @return None.
!>
!> @param v [in] Argument v.
!> @param x [in] Argument x.
!> @param vm [inout] Argument vm.
!> @param bi [inout] Argument bi.
!> @param di [inout] Argument di.
!> @param bk [inout] Argument bk.
!> @param dk [inout] Argument dk.
subroutine ikv ( v, x, vm, bi, di, bk, dk )

!*****************************************************************************80
!
!! IKV compute modified Bessel function Iv(x) and Kv(x) and their derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    17 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) V, the order of Iv(x) and Kv(x).
!    V = N + V0.
!
!    Input, real(real64) X, the argument.
!
!    Output, real(real64) VM, the highest order computed.
!
!    Output, real(real64) BI(0:N), DI(0:N), BK(0:N), DK(0:N), the
!    values of In+v0(x), In+v0'(x), Kn+v0(x), Kn+v0'(x).
!
  implicit none

  real(real64) a1
  real(real64) a2
  real(real64), intent(inout) :: bi(0:)
  real(real64) bi0
  real(real64), intent(inout) :: bk(0:)
  real(real64) bk0
  real(real64) bk1
  real(real64) bk2
  real(real64) ca
  real(real64) cb
  real(real64) cs
  real(real64) ct
  real(real64), intent(inout) :: di(0:)
  real(real64), intent(inout) :: dk(0:)
  real(real64) f
  real(real64) f1
  real(real64) f2
  real(real64) gan
  real(real64) gap
  integer(int32) k
  integer(int32) k0
  integer(int32) m
  ! integer(int32) msta1
  ! integer(int32) msta2
  integer(int32) n
  real(real64) pi
  real(real64) piv
  real(real64) r
  real(real64) r1
  real(real64) r2
  real(real64) sum
  real(real64), intent(in) :: v
  real(real64) v0
  real(real64) v0n
  real(real64) v0p
  real(real64), intent(inout) :: vm
  real(real64) vt
  real(real64) w0
  real(real64) wa
  real(real64) ww
  real(real64), intent(in) :: x
  real(real64) x2

  pi = 3.141592653589793e+00_real64
  x2 = x * x
  n = int ( v )
  v0 = v - n
  if ( n == 0 ) then
    n = 1
  end if

  if ( x < 1.0e-100_real64 ) then

    do k = 0, n
      bi(k) = 0.0e+00_real64
      di(k) = 0.0e+00_real64
      bk(k) = -1.0e+300_real64
      dk(k) = 1.0e+300_real64
    end do

    if ( v == 0.0e+00_real64 ) then
      bi(0) = 1.0e+00_real64
      di(1) = 0.5e+00_real64
    end if

    vm = v
    return

  end if

  piv = pi * v0
  vt = 4.0e+00_real64 * v0 * v0

  if ( v0 == 0.0e+00_real64 ) then
    a1 = 1.0e+00_real64
  else
    v0p = 1.0e+00_real64 + v0
    call gamma ( v0p, gap )
    a1 = ( 0.5e+00_real64 * x ) ** v0 / gap
  end if

  if ( x < 35.0e+00_real64 ) then
    k0 = 14
  else if ( x < 50.0e+00_real64 ) then
    k0 = 10
  else
    k0 = 8
  end if
 
  if ( x <= 18.0e+00_real64 ) then

    bi0 = 1.0e+00_real64
    r = 1.0e+00_real64
    do k = 1, 30
      r = 0.25e+00_real64 * r * x2 / ( k * ( k + v0 ) )
      bi0 = bi0 + r
      if ( abs ( r / bi0 ) < 1.0e-15_real64 ) then
        exit
      end if
    end do

    bi0 = bi0 * a1

  else

    ca = exp ( x ) / sqrt ( 2.0e+00_real64 * pi * x )
    sum = 1.0e+00_real64
    r = 1.0e+00_real64
    do k = 1, k0
      r = -0.125e+00_real64 * r * ( vt - ( 2.0e+00_real64 * k - 1.0e+00_real64 ) ** 2 ) / ( k * x )
      sum = sum + r
    end do
    bi0 = ca * sum

  end if

  m = msta1 ( x, 200 )

  if ( m < n ) then
    n = m
  else
    m = msta2 ( x, n, 15 )
  end if

  f2 = 0.0e+00_real64
  f1 = 1.0e-100_real64
  do k = m, 0, -1
    f = 2.0e+00_real64 * ( v0 + k + 1.0e+00_real64 ) / x * f1 + f2
    if ( k <= n ) then
      bi(k) = f
    end if
    f2 = f1
    f1 = f
  end do

  cs = bi0 / f
  do k = 0, n
    bi(k) = cs * bi(k)
  end do

  di(0) = v0 / x * bi(0) + bi(1)
  do k = 1, n
    di(k) = - ( k + v0 ) / x * bi(k) + bi(k-1)
  end do

  if ( x <= 9.0e+00_real64 ) then

    if ( v0 == 0.0e+00_real64 ) then

      ct = - log ( 0.5e+00_real64 * x ) - 0.5772156649015329e+00_real64
      cs = 0.0e+00_real64
      w0 = 0.0e+00_real64
      ww = 0.0e+00_real64
      r = 1.0e+00_real64
      do k = 1, 50
        w0 = w0 + 1.0e+00_real64 / k
        r = 0.25e+00_real64 * r / ( k * k ) * x2
        cs = cs + r * ( w0 + ct )
        wa = abs ( cs )
        if ( abs ( ( wa - ww ) / wa ) < 1.0e-15_real64 ) then
          exit
        end if
        ww = wa
      end do

      bk0 = ct + cs

    else

      v0n = 1.0e+00_real64 - v0
      call gamma ( v0n, gan )
      a2 = 1.0e+00_real64 / ( gan * ( 0.5e+00_real64 * x ) ** v0 )
      a1 = ( 0.5e+00_real64 * x ) ** v0 / gap
      sum = a2 - a1
      r1 = 1.0e+00_real64
      r2 = 1.0e+00_real64
      ww = 0.0e+00_real64
      do k = 1, 120
        r1 = 0.25e+00_real64 * r1 * x2 / ( k * ( k - v0 ) )
        r2 = 0.25e+00_real64 * r2 * x2 / ( k * ( k + v0 ) )
        sum = sum + a2 * r1 - a1 * r2
        wa = abs ( sum )
        if ( abs ( ( wa - ww ) / wa ) < 1.0e-15_real64 ) then
          exit
        end if
        ww = wa
      end do

      bk0 = 0.5e+00_real64 * pi * sum / sin ( piv )

    end if

  else

    cb = exp ( - x ) * sqrt ( 0.5e+00_real64 * pi / x )
    sum = 1.0e+00_real64
    r = 1.0e+00_real64
    do k = 1, k0
      r = 0.125e+00_real64 * r * ( vt - ( 2.0e+00_real64 * k - 1.0e+00_real64 ) ** 2 ) / ( k * x )
      sum = sum + r
    end do
    bk0 = cb * sum

  end if

  bk1 = ( 1.0e+00_real64 / x - bi(1) * bk0 ) / bi(0)
  bk(0) = bk0
  bk(1) = bk1
  do k = 2, n
    bk2 = 2.0e+00_real64 * ( v0 + k - 1.0e+00_real64 ) / x * bk1 + bk0
    bk(k) = bk2
    bk0 = bk1
    bk1 = bk2
  end do

  dk(0) = v0 / x * bk(0) - bk(1)
  do k = 1, n
    dk(k) = - ( k + v0 ) / x * bk(k) - bk(k-1)
  end do

  vm = n + v0

  return
end subroutine ikv
!> @brief subroutine incob.
!> @return None.
!>
!> @param a [in] Argument a.
!> @param b [in] Argument b.
!> @param x [in] Argument x.
!> @param bix [inout] Argument bix.
subroutine incob ( a, b, x, bix )

!*****************************************************************************80
!
!! INCOB computes the incomplete beta function Ix(a,b).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    22 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) A, B, parameters.
!
!    Input, real(real64) X, the argument.
!
!    Output, real(real64) BIX, the function value.
!
  implicit none

  real(real64), intent(in) :: a
  real(real64), intent(in) :: b
  real(real64), intent(inout) :: bix
  real(real64) bt
  real(real64) dk(51)
  real(real64) fk(51)
  integer(int32) k
  real(real64) s0
  real(real64) t1
  real(real64) t2
  real(real64) ta
  real(real64) tb
  real(real64), intent(in) :: x

  s0 = ( a + 1.0e+00_real64 ) / ( a + b + 2.0e+00_real64 )
  call beta ( a, b, bt )

  if ( x <= s0 ) then

    do k = 1, 20
      dk(2*k) = k * ( b - k ) * x / &
        ( a + 2.0e+00_real64 * k - 1.0e+00_real64 ) / ( a + 2.0e+00_real64 * k )
    end do

    do k = 0, 20
      dk(2*k+1) = - ( a + k ) * ( a + b + k ) * x &
        / ( a + 2.0e+00_real64 * k ) / ( a + 2.0e+00_real64 * k + 1.0e+00_real64 )
    end do

    t1 = 0.0e+00_real64
    do k = 20, 1, -1
      t1 = dk(k) / ( 1.0e+00_real64 + t1 )
    end do
    ta = 1.0e+00_real64 / ( 1.0e+00_real64 + t1 )
    bix = x ** a * ( 1.0e+00_real64 - x ) ** b / ( a * bt ) * ta

  else

    do k = 1, 20
      fk(2*k) = k * ( a - k ) * ( 1.0e+00_real64 - x ) &
        / ( b + 2.0e+00_real64 * k - 1.0e+00_real64 ) / ( b + 2.0e+00_real64 * k )
    end do

    do k = 0,20
      fk(2*k+1) = - ( b + k ) * ( a + b + k ) * ( 1.0e+00_real64 - x ) &
        / ( b + 2.0e+00_real64 * k ) / ( b + 2.0e+00_real64 * k + 1.0e+00_real64 )
    end do

    t2 = 0.0e+00_real64
    do k = 20, 1, -1
      t2 = fk(k) / ( 1.0e+00_real64 + t2 )
    end do
    tb = 1.0e+00_real64 / ( 1.0e+00_real64 + t2 )
    bix = 1.0e+00_real64 - x ** a * ( 1.0e+00_real64 - x ) ** b / ( b * bt ) * tb

  end if

  return
end subroutine incob
!> @brief subroutine incog.
!> @return None.
!>
!> @param a [in] Argument a.
!> @param x [in] Argument x.
!> @param gin [inout] Argument gin.
!> @param gim [inout] Argument gim.
!> @param gip [inout] Argument gip.
subroutine incog ( a, x, gin, gim, gip )

!*****************************************************************************80
!
!! INCOG computes the incomplete gamma function r(a,x), ,(a,x), P(a,x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    22 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) A, the parameter.
!
!    Input, real(real64) X, the argument.
!
!    Output, real(real64) GIN, GIM, GIP, the values of
!    r(a,x), �(a,x), P(a,x).
!
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  implicit none

  real(real64), intent(in) :: a
  real(real64) ga
  real(real64), intent(inout) :: gim
  real(real64), intent(inout) :: gin
  real(real64), intent(inout) :: gip
  integer(int32) k
  real(real64) r
  real(real64) s
  real(real64) t0
  real(real64), intent(in) :: x
  real(real64) xam

  xam = -  x + a * log ( x )

  if ( 700.0e+00_real64 < xam .or. 170.0e+00_real64 < a ) then
    write ( error_unit, '(a)' ) ' '
    write ( error_unit, '(a)' ) 'INCOG - Fatal error!'
    write ( error_unit, '(a)' ) '  A and/or X is too large!'
    gin = ieee_value ( 0.0e+00_real64, ieee_quiet_nan )
    gim = ieee_value ( 0.0e+00_real64, ieee_quiet_nan )
    gip = ieee_value ( 0.0e+00_real64, ieee_quiet_nan )
    return
  end if

  if ( x == 0.0e+00_real64 ) then

    gin = 0.0e+00_real64
    call gamma ( a, ga )
    gim = ga
    gip = 0.0e+00_real64

  else if ( x <= 1.0e+00_real64 + a ) then

    s = 1.0e+00_real64 / a
    r = s
    do k = 1, 60
      r = r * x / ( a + k )
      s = s + r
      if ( abs ( r / s ) < 1.0e-15_real64 ) then
        exit
      end if
    end do

    gin = exp ( xam ) * s
    call gamma ( a, ga )
    gip = gin / ga
    gim = ga - gin

  else if ( 1.0e+00_real64 + a < x ) then

    t0 = 0.0e+00_real64
    do k = 60, 1, -1
      t0 = ( k - a ) / ( 1.0e+00_real64 + k / ( x + t0 ) )
    end do
    gim = exp ( xam ) / ( x + t0 )
    call gamma ( a, ga )
    gin = ga - gim
    gip = 1.0e+00_real64 - gim / ga

  end if

  return
end subroutine incog
!> @brief subroutine itairy.
!> @return None.
!>
!> @param x [inout] Argument x.
!> @param apt [inout] Argument apt.
!> @param bpt [inout] Argument bpt.
!> @param ant [inout] Argument ant.
!> @param bnt [inout] Argument bnt.
pure subroutine itairy ( x, apt, bpt, ant, bnt )

!****************************************************************************80
!
!! ITAIRY computes the integrals of Airy functions.
!
!  Discussion:
!
!    Compute the integrals of Airy functions with respect to t,
!    from 0 and x.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    19 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) X, the upper limit of the integral.
!
!    Output, real(real64) APT, BPT, ANT, BNT, the integrals, from 0 to x,
!    of Ai(t), Bi(t), Ai(-t), and Bi(-t).
!       
  implicit none

  real(real64), parameter, dimension ( 16 ) :: a = [&
    0.569444444444444e+00_real64, 0.891300154320988e+00_real64, &
    0.226624344493027e+01_real64, 0.798950124766861e+01_real64, &
    0.360688546785343e+02_real64, 0.198670292131169e+03_real64, &
    0.129223456582211e+04_real64, 0.969483869669600e+04_real64, &
    0.824184704952483e+05_real64, 0.783031092490225e+06_real64, &
    0.822210493622814e+07_real64, 0.945557399360556e+08_real64, &
    0.118195595640730e+10_real64, 0.159564653040121e+11_real64, &
    0.231369166433050e+12_real64, 0.358622522796969e+13_real64]
  real(real64), intent(inout) :: ant
  real(real64), intent(inout) :: apt
  real(real64), intent(inout) :: bnt
  real(real64), intent(inout) :: bpt
  real(real64) c1
  real(real64) c2
  real(real64) eps
  real(real64) fx
  real(real64) gx
  integer(int32) k
  integer(int32) l
  real(real64) pi
  real(real64) q0
  real(real64) q1
  real(real64) q2
  real(real64) r
  real(real64) sr3
  real(real64) su1
  real(real64) su2
  real(real64) su3
  real(real64) su4
  real(real64) su5
  real(real64) su6
  real(real64), intent(inout) :: x
  real(real64) xe
  real(real64) xp6
  real(real64) xr1
  real(real64) xr2

  eps = 1.0e-15_real64
  pi = 3.141592653589793e+00_real64
  c1 = 0.355028053887817e+00_real64
  c2 = 0.258819403792807e+00_real64
  sr3 = 1.732050807568877e+00_real64

  if ( x == 0.0e+00_real64 ) then

    apt = 0.0e+00_real64
    bpt = 0.0e+00_real64
    ant = 0.0e+00_real64
    bnt = 0.0e+00_real64

  else

    if ( abs ( x ) <= 9.25e+00_real64 ) then

      do l = 0, 1

        x = ( -1.0e+00_real64 ) ** l * x
        fx = x
        r = x

        do k = 1, 40
          r = r * ( 3.0e+00_real64 * k - 2.0e+00_real64 ) &
            / ( 3.0e+00_real64 * k + 1.0e+00_real64 ) * x / ( 3.0e+00_real64 * k ) &
            * x / ( 3.0e+00_real64 * k - 1.0e+00_real64 ) * x 
          fx = fx + r
          if ( abs ( r ) < abs ( fx ) * eps ) then
            exit
          end if
        end do

        gx = 0.5e+00_real64 * x * x
        r = gx

        do k = 1, 40
          r = r * ( 3.0e+00_real64 * k - 1.0e+00_real64 ) &
            / ( 3.0e+00_real64 * k + 2.0e+00_real64 ) * x / ( 3.0e+00_real64 * k ) * x &
            / ( 3.0e+00_real64 * k + 1.0e+00_real64 ) * x
          gx = gx + r
          if ( abs ( r ) < abs ( gx ) * eps ) then
            exit
          end if
        end do

        ant = c1 * fx - c2 * gx
        bnt = sr3 * ( c1 * fx + c2 * gx )

        if ( l == 0 ) then
          apt = ant
          bpt = bnt
        else
          ant = -ant
          bnt = -bnt
          x = -x
        end if

      end do

    else

      q2 = 1.414213562373095e+00_real64
      q0 = 0.3333333333333333e+00_real64
      q1 = 0.6666666666666667e+00_real64
      xe = x * sqrt ( x ) / 1.5e+00_real64
      xp6 = 1.0e+00_real64 / sqrt ( 6.0e+00_real64 * pi * xe )
      su1 = 1.0e+00_real64
      r = 1.0e+00_real64
      xr1 = 1.0e+00_real64 / xe
      do k = 1, 16
        r = - r * xr1
        su1 = su1 + a(k) * r
      end do
      su2 = 1.0e+00_real64
      r = 1.0e+00_real64
      do k = 1, 16
        r = r * xr1
        su2 = su2 + a(k) * r
      end do

      apt = q0 - exp ( - xe ) * xp6 * su1
      bpt = 2.0e+00_real64 * exp ( xe ) * xp6 * su2
      su3 = 1.0e+00_real64
      r = 1.0e+00_real64
      xr2 = 1.0e+00_real64 / ( xe * xe )
      do k = 1, 8
        r = - r * xr2
        su3 = su3 + a(2*k) * r
      end do
      su4 = a(1) * xr1
      r = xr1
      do k = 1, 7
        r = -r * xr2
        su4 = su4 + a(2*k+1) * r
      end do
      su5 = su3 + su4
      su6 = su3 - su4
      ant = q1 - q2 * xp6 * ( su5 * cos ( xe ) - su6 * sin ( xe ) )
      bnt = q2 * xp6 * ( su5 * sin ( xe ) + su6 * cos ( xe ) )

    end if

  end if

  return
end subroutine itairy
!> @brief subroutine itika.
!> @return None.
!>
!> @param x [in] Argument x.
!> @param ti [inout] Argument ti.
!> @param tk [inout] Argument tk.
pure subroutine itika ( x, ti, tk )

!*****************************************************************************80
!
!! ITIKA computes the integral of the modified Bessel functions I0(t) and K0(t).
!
!  Discussion:
!
!    This procedure integrates modified Bessel functions I0(t) and
!    K0(t) with respect to t from 0 to x.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    18 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) X, the upper limit of the integral.
!
!    Output, real(real64) TI, TK, the integrals of I0(t) and K0(t)
!    from 0 to X.
!
  implicit none

  real(real64), parameter, dimension ( 10 ) :: a = [&
    0.625e+00_real64,           1.0078125e+00_real64, &
    2.5927734375e+00_real64,    9.1868591308594e+00_real64, &
    4.1567974090576e+01_real64, 2.2919635891914e+02_real64, &
    1.491504060477e+03_real64,  1.1192354495579e+04_real64, &
    9.515939374212e+04_real64,  9.0412425769041e+05_real64]
  real(real64) b1
  real(real64) b2
  real(real64) e0
  real(real64) el
  integer(int32) k
  real(real64) pi
  real(real64) r
  real(real64) rc1
  real(real64) rc2
  real(real64) rs
  real(real64), intent(inout) :: ti
  real(real64), intent(inout) :: tk
  real(real64) tw
  real(real64), intent(in) :: x
  real(real64) x2

  pi = 3.141592653589793e+00_real64
  el = 0.5772156649015329e+00_real64

  if ( x == 0.0e+00_real64 ) then

    ti = 0.0e+00_real64
    tk = 0.0e+00_real64
    return

  else if ( x < 20.0e+00_real64 ) then

    x2 = x * x
    ti = 1.0e+00_real64
    r = 1.0e+00_real64
    do k = 1, 50
      r = 0.25e+00_real64 * r * ( 2 * k - 1.0e+00_real64 ) / ( 2 * k + 1.0e+00_real64 ) &
        / ( k * k ) * x2
      ti = ti + r
      if ( abs ( r / ti ) < 1.0e-12_real64 ) then
        exit
      end if
    end do

    ti = ti * x

  else

    ti = 1.0e+00_real64
    r = 1.0e+00_real64
    do k = 1, 10
      r = r / x
      ti = ti + a(k) * r
    end do
    rc1 = 1.0e+00_real64 / sqrt ( 2.0e+00_real64 * pi * x )
    ti = rc1 * exp ( x ) * ti

  end if

  if ( x < 12.0e+00_real64 ) then

    e0 = el + log ( x / 2.0e+00_real64 )
    b1 = 1.0e+00_real64 - e0
    b2 = 0.0e+00_real64
    tw = 0.0e+00_real64
    rs = 0.0e+00_real64
    r = 1.0e+00_real64
    do k = 1, 50
      r = 0.25e+00_real64 * r * ( 2 * k - 1.0e+00_real64 ) &
        / ( 2 * k + 1.0e+00_real64 ) / ( k * k ) * x2
      b1 = b1 + r * ( 1.0e+00_real64 / ( 2 * k + 1 ) - e0 )
      rs = rs + 1.0e+00_real64 / k
      b2 = b2 + r * rs
      tk = b1 + b2
      if ( abs ( ( tk - tw ) / tk ) < 1.0e-12_real64 ) then
        exit
      end if
      tw = tk
    end do

    tk = tk * x

  else

    tk = 1.0e+00_real64
    r = 1.0e+00_real64
    do k = 1, 10
      r = -r / x
      tk = tk + a(k) * r
    end do 
    rc2 = sqrt ( pi / ( 2.0e+00_real64 * x ) )
    tk = pi / 2.0e+00_real64 - rc2 * tk * exp ( - x )

  end if

  return
end subroutine itika
!> @brief subroutine itikb.
!> @return None.
!>
!> @param x [in] Argument x.
!> @param ti [inout] Argument ti.
!> @param tk [inout] Argument tk.
pure subroutine itikb ( x, ti, tk )

!*****************************************************************************80
!
!! ITIKB computes the integral of the Bessel functions I0(t) and K0(t).
!
!  Discussion:
!
!    This procedure integrates Bessel functions I0(t) and K0(t)
!    with respect to t from 0 to x.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    24 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
! 
!  Parameters:
!
!    Input, real(real64) X, the upper limit of the integral.
!
!    Output, real(real64) TI, TK, the integral of I0(t) and K0(t)
!    from 0 to X.
!
  implicit none

  real(real64) pi
  real(real64) t
  real(real64) t1
  real(real64), intent(inout) :: ti
  real(real64), intent(inout) :: tk
  real(real64), intent(in) :: x

  pi = 3.141592653589793e+00_real64

  if ( x == 0.0e+00_real64 ) then

    ti = 0.0e+00_real64

  else if ( x < 5.0e+00_real64 ) then

    t1 = x / 5.0e+00_real64
    t = t1 * t1
    ti = (((((((( &
        0.59434e-03_real64 * t &
      + 0.4500642e-02_real64 ) * t &
      + 0.044686921e+00_real64 ) * t &
      + 0.300704878e+00_real64 ) * t &
      + 1.471860153e+00_real64 ) * t &
      + 4.844024624e+00_real64 ) * t &
      + 9.765629849e+00_real64 ) * t &
      +10.416666367e+00_real64 ) * t &
      + 5.0e+00_real64 ) * t1

  else if ( 5.0e+00_real64 <= x .and. x <= 8.0e+00_real64 ) then

    t = 5.0e+00_real64 / x
    ti = ((( &
      - 0.015166e+00_real64 * t &
      - 0.0202292e+00_real64 ) * t &
      + 0.1294122e+00_real64 ) * t &
      - 0.0302912e+00_real64 ) * t &
      + 0.4161224e+00_real64
    ti = ti * exp ( x ) / sqrt ( x )

  else

    t = 8.0e+00_real64 / x
    ti = ((((( &
      - 0.0073995e+00_real64 * t &
      + 0.017744e+00_real64 ) * t &
      - 0.0114858e+00_real64 ) * t &
      + 0.55956e-02_real64 ) * t &
      + 0.59191e-02_real64 ) * t &
      + 0.0311734e+00_real64 ) * t &
      + 0.3989423e+00_real64
    ti = ti * exp ( x ) / sqrt ( x )

  end if

  if ( x == 0.0e+00_real64 ) then

    tk = 0.0e+00_real64

  else if ( x <= 2.0e+00_real64 ) then

    t1 = x / 2.0e+00_real64
    t = t1 * t1
    tk = (((((( &
        0.116e-05_real64        * t &
      + 0.2069e-04_real64 )     * t &
      + 0.62664e-03_real64 )    * t &
      + 0.01110118e+00_real64 ) * t &
      + 0.11227902e+00_real64 ) * t &
      + 0.50407836e+00_real64 ) * t &
      + 0.84556868e+00_real64 ) * t1
    tk = tk - log ( x / 2.0e+00_real64 ) * ti

  else if ( 2.0e+00_real64 < x .and. x <= 4.0e+00_real64 ) then

    t = 2.0e+00_real64 / x
    tk = ((( &
        0.0160395e+00_real64   * t &
      - 0.0781715e+00_real64 ) * t &
      + 0.185984e+00_real64 )  * t &
      - 0.3584641e+00_real64 ) * t &
      + 1.2494934e+00_real64
    tk = pi / 2.0e+00_real64 - tk * exp ( - x ) / sqrt ( x )

  else if ( 4.0e+00_real64 < x .and. x <= 7.0e+00_real64 ) then

    t = 4.0e+00_real64 / x
    tk = ((((( &
        0.37128e-02_real64 * t &
      - 0.0158449e+00_real64 ) * t &
      + 0.0320504e+00_real64 ) * t &
      - 0.0481455e+00_real64 ) * t &
      + 0.0787284e+00_real64 ) * t &
      - 0.1958273e+00_real64 ) * t &
      + 1.2533141e+00_real64
    tk = pi / 2.0e+00_real64 - tk * exp ( - x ) / sqrt ( x )

  else

    t = 7.0e+00_real64 / x
    tk = ((((( &
        0.33934e-03_real64      * t &
      - 0.163271e-02_real64 )   * t &
      + 0.417454e-02_real64 )   * t &
      - 0.933944e-02_real64 )   * t &
      + 0.02576646e+00_real64 ) * t &
      - 0.11190289e+00_real64 ) * t &
      + 1.25331414e+00_real64
    tk = pi / 2.0e+00_real64 - tk * exp ( - x ) / sqrt ( x )

  end if

  return
end subroutine itikb
!> @brief subroutine itjya.
!> @return None.
!>
!> @param x [in] Argument x.
!> @param tj [inout] Argument tj.
!> @param ty [inout] Argument ty.
pure subroutine itjya ( x, tj, ty )

!*****************************************************************************80
!
!! ITJYA computes integrals of Bessel functions J0(t) and Y0(t).
!
!  Discussion:
!
!    This procedure integrates Bessel functions J0(t) and Y0(t) with
!    respect to t from 0 to x.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    25 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) X, the upper limit of the integral.
!
!    Output, real(real64) TJ, TY, the integrals of J0(t) and Y0(t) 
!    from 0 to x.
!
  implicit none

  real(real64) a(18)
  real(real64) a0
  real(real64) a1
  real(real64) af
  real(real64) bf
  real(real64) bg
  real(real64) el
  real(real64) eps
  integer(int32) k
  real(real64) pi
  real(real64) r
  real(real64) r2
  real(real64) rc
  real(real64) rs
  real(real64), intent(inout) :: tj
  real(real64), intent(inout) :: ty
  real(real64) ty1
  real(real64) ty2
  real(real64), intent(in) :: x
  real(real64) x2
  real(real64) xp

  pi = 3.141592653589793e+00_real64
  el = 0.5772156649015329e+00_real64
  eps = 1.0e-12_real64

  if ( x == 0.0e+00_real64 ) then

    tj = 0.0e+00_real64
    ty = 0.0e+00_real64

  else if ( x <= 20.0e+00_real64 ) then

    x2 = x * x
    tj = x
    r = x
    do k = 1, 60
      r = -0.25e+00_real64 * r * ( 2 * k - 1.0e+00_real64 ) / ( 2 * k + 1.0e+00_real64 ) &
        / ( k * k ) * x2
      tj = tj + r
      if ( abs ( r ) < abs ( tj ) * eps ) then
        exit
      end if
    end do

    ty1 = ( el + log ( x / 2.0e+00_real64 ) ) * tj
    rs = 0.0e+00_real64
    ty2 = 1.0e+00_real64
    r = 1.0e+00_real64

    do k = 1, 60
      r = -0.25e+00_real64 * r * ( 2 * k - 1.0e+00_real64 ) / ( 2 * k + 1.0e+00_real64 ) &
        / ( k * k ) * x2
      rs = rs + 1.0e+00_real64 / k
      r2 = r * ( rs + 1.0e+00_real64 / ( 2.0e+00_real64 * k + 1.0e+00_real64 ) )
      ty2 = ty2 + r2
      if ( abs ( r2 ) < abs ( ty2 ) * eps ) then
        exit
      end if
    end do

    ty = ( ty1 - x * ty2 ) * 2.0e+00_real64 / pi

  else

    a0 = 1.0e+00_real64
    a1 = 5.0e+00_real64 / 8.0e+00_real64
    a(1) = a1

    do k = 1, 16
      af = ( ( 1.5e+00_real64 * ( k + 0.5e+00_real64 ) * ( k + 5.0e+00_real64 / 6.0e+00_real64 ) &
        * a1 - 0.5e+00_real64 * ( k + 0.5e+00_real64 ) * ( k + 0.5e+00_real64 )  &
        * ( k - 0.5e+00_real64 ) * a0 ) ) / ( k + 1.0e+00_real64 )
      a(k+1) = af
      a0 = a1
      a1 = af
    end do

    bf = 1.0e+00_real64
    r = 1.0e+00_real64
    do k = 1, 8
      r = -r / ( x * x )
      bf = bf + a(2*k) * r
    end do
    bg = a(1) / x
    r = 1.0e+00_real64 / x
    do k = 1, 8
      r = -r / ( x * x )
      bg = bg + a(2*k+1) * r
    end do
    xp = x + 0.25e+00_real64 * pi
    rc = sqrt ( 2.0e+00_real64 / ( pi * x ) )
    tj = 1.0e+00_real64 - rc * ( bf * cos ( xp ) + bg * sin ( xp ) )
    ty = rc * ( bg * cos ( xp ) - bf * sin ( xp ) )

  end if

  return
end subroutine itjya
!> @brief subroutine itjyb.
!> @return None.
!>
!> @param x [in] Argument x.
!> @param tj [inout] Argument tj.
!> @param ty [inout] Argument ty.
pure subroutine itjyb ( x, tj, ty )

!*****************************************************************************80
!
!! ITJYB computes integrals of Bessel functions J0(t) and Y0(t).
!
!  Discussion:
!
!    This procedure integrates Bessel functions J0(t) and Y0(t)
!    with respect to t from 0 to x.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    25 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) X, the upper limit of the integral.
!
!    Output, real(real64) TJ, TY, the integrals of J0(t) and Y0(t) 
!    from 0 to x.
!
  implicit none

  real(real64) f0
  real(real64) g0
  real(real64) pi
  real(real64) t
  real(real64), intent(inout) :: tj
  real(real64), intent(inout) :: ty
  real(real64), intent(in) :: x
  real(real64) x1
  real(real64) xt

  pi = 3.141592653589793e+00_real64

  if ( x == 0.0e+00_real64 ) then

    tj = 0.0e+00_real64
    ty = 0.0e+00_real64

  else if ( x <= 4.0e+00_real64 ) then

    x1 = x / 4.0e+00_real64
    t = x1 * x1

    tj = ((((((( &
      - 0.133718e-03_real64      * t &
      + 0.2362211e-02_real64 )   * t &
      - 0.025791036e+00_real64 ) * t &
      + 0.197492634e+00_real64 ) * t &
      - 1.015860606e+00_real64 ) * t &
      + 3.199997842e+00_real64 ) * t &
      - 5.333333161e+00_real64 ) * t &
      + 4.0e+00_real64 ) * x1

    ty = (((((((( &
        0.13351e-04_real64       * t &
      - 0.235002e-03_real64 )    * t &
      + 0.3034322e-02_real64 )   * t &
      - 0.029600855e+00_real64 ) * t &
      + 0.203380298e+00_real64 ) * t &
      - 0.904755062e+00_real64 ) * t &
      + 2.287317974e+00_real64 ) * t &
      - 2.567250468e+00_real64 ) * t &
      + 1.076611469e+00_real64 ) * x1

    ty = 2.0e+00_real64 / pi * log ( x / 2.0e+00_real64 ) * tj - ty

  else if ( x <= 8.0e+00_real64 ) then

    xt = x - 0.25e+00_real64 * pi
    t = 16.0e+00_real64 / ( x * x )

    f0 = (((((( &
        0.1496119e-02_real64     * t &
      - 0.739083e-02_real64 )    * t &
      + 0.016236617e+00_real64 ) * t &
      - 0.022007499e+00_real64 ) * t &
      + 0.023644978e+00_real64 ) * t &
      - 0.031280848e+00_real64 ) * t &
      + 0.124611058e+00_real64 ) * 4.0e+00_real64 / x

    g0 = ((((( &
        0.1076103e-02_real64     * t &
      - 0.5434851e-02_real64 )   * t &
      + 0.01242264e+00_real64 )  * t &
      - 0.018255209e+00_real64 ) * t &
      + 0.023664841e+00_real64 ) * t &
      - 0.049635633e+00_real64 ) * t &
      + 0.79784879e+00_real64

    tj = 1.0e+00_real64 - ( f0 * cos ( xt ) - g0 * sin ( xt ) ) / sqrt ( x )

    ty = - ( f0 * sin ( xt ) + g0 * cos ( xt ) ) / sqrt ( x )

  else

    t = 64.0e+00_real64 / ( x * x )
    xt = x-0.25e+00_real64 * pi

    f0 = ((((((( &
      - 0.268482e-04_real64     * t &
      + 0.1270039e-03_real64 )  * t &
      - 0.2755037e-03_real64 )  * t &
      + 0.3992825e-03_real64 )  * t &
      - 0.5366169e-03_real64 )  * t &
      + 0.10089872e-02_real64 ) * t &
      - 0.40403539e-02_real64 ) * t &
      + 0.0623347304e+00_real64 ) * 8.0e+00_real64 / x

    g0 = (((((( &
      - 0.226238e-04_real64        * t &
      + 0.1107299e-03_real64 )     * t &
      - 0.2543955e-03_real64 )     * t &
      + 0.4100676e-03_real64 )     * t &
      - 0.6740148e-03_real64 )     * t &
      + 0.17870944e-02_real64 )    * t &
      - 0.01256424405e+00_real64 ) * t &
      + 0.79788456e+00_real64

    tj = 1.0e+00_real64  - ( f0 * cos ( xt ) - g0 * sin ( xt ) ) / sqrt ( x )

    ty = - ( f0 * sin ( xt ) + g0 * cos ( xt ) ) / sqrt ( x )

  end if

  return
end subroutine itjyb
!> @brief subroutine itsh0.
!> @return None.
!>
!> @param x [in] Argument x.
!> @param th0 [inout] Argument th0.
pure subroutine itsh0 ( x, th0 )

!*****************************************************************************80
!
!! ITSH0 integrates the Struve function H0(t) from 0 to x.
!
!  Discussion:
!
!    This procedure evaluates the integral of Struve function
!    H0(t) with respect to t from 0 and x.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    25 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) X, the upper limit of the integral.
!
!    Output, real(real64) TH0, the integral of H0(t) from 0 to x.
!
  implicit none

  real(real64) a(25)
  real(real64) a0
  real(real64) a1
  real(real64) af
  real(real64) bf
  real(real64) bg
  real(real64) el
  integer(int32) k
  real(real64) pi
  real(real64) r
  real(real64) rd
  real(real64) s
  real(real64) s0
  real(real64), intent(inout) :: th0
  real(real64) ty
  real(real64), intent(in) :: x
  real(real64) xp

  pi = 3.141592653589793e+00_real64
  r = 1.0e+00_real64            

  if ( x <= 30.0e+00_real64 ) then

    s = 0.5e+00_real64

    do k = 1, 100

      if ( k == 1 ) then
        rd = 0.5e+00_real64
      else
        rd = 1.0e+00_real64
      end if

      r = - r * rd * k / ( k + 1.0e+00_real64 ) &
        * ( x / ( 2.0e+00_real64 * k + 1.0e+00_real64 ) ) ** 2
      s = s + r

      if ( abs ( r ) < abs ( s ) * 1.0e-12_real64 ) then
        exit
      end if

    end do

    th0 = 2.0e+00_real64 / pi * x * x * s

  else

    s = 1.0e+00_real64
    do k = 1, 12
      r = - r * k / ( k + 1.0e+00_real64 ) &
        * ( ( 2.0e+00_real64 * k + 1.0e+00_real64 ) / x ) ** 2
      s = s + r
      if ( abs ( r ) < abs ( s ) * 1.0e-12_real64 ) then
        exit
      end if
    end do

    el = 0.57721566490153e+00_real64
    s0 = s / ( pi * x * x ) + 2.0e+00_real64 / pi &
      * ( log ( 2.0e+00_real64 * x ) + el )
    a0 = 1.0e+00_real64
    a1 = 5.0e+00_real64 / 8.0e+00_real64
    a(1) = a1
    do k = 1, 20
      af = ( ( 1.5e+00_real64 * ( k + 0.5e+00_real64 ) &
        * ( k + 5.0e+00_real64 / 6.0e+00_real64 ) * a1 - 0.5e+00_real64 &
        * ( k + 0.5e+00_real64 ) * ( k + 0.5e+00_real64 ) &
        * ( k - 0.5e+00_real64 ) * a0 ) ) / ( k + 1.0e+00_real64 )
      a(k+1) = af
      a0 = a1
      a1 = af
    end do

    bf = 1.0e+00_real64
    r = 1.0e+00_real64
    do k = 1, 10
      r = - r / ( x * x )
      bf = bf + a(2*k) * r
    end do
    bg = a(1) / x
    r = 1.0e+00_real64 / x
    do k = 1, 10
      r = - r / ( x * x ) 
      bg = bg + a(2*k+1) * r
    end do
    xp = x + 0.25e+00_real64 * pi
    ty = sqrt ( 2.0e+00_real64 / ( pi * x ) ) &
      * ( bg * cos ( xp ) - bf * sin ( xp ) )
    th0 = ty + s0

  end if

  return
end subroutine itsh0
!> @brief subroutine itsl0.
!> @return None.
!>
!> @param x [in] Argument x.
!> @param tl0 [inout] Argument tl0.
pure subroutine itsl0 ( x, tl0 )

!*****************************************************************************80
!
!! ITSL0 integrates the Struve function L0(t) from 0 to x.
!
!  Discussion:
!
!    This procedure evaluates the integral of modified Struve function
!    L0(t) with respect to t from 0 to x.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    31 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) X, the upper limit of the integral.
!
!    Output, real(real64) TL0, the integral of L0(t) from 0 to x.
!
  implicit none

  real(real64) a(18)
  real(real64) a0
  real(real64) a1
  real(real64) af
  real(real64) el
  integer(int32) k
  real(real64) pi
  real(real64) r
  real(real64) rd
  real(real64) s
  real(real64) s0
  real(real64) ti
  real(real64), intent(inout) :: tl0
  real(real64), intent(in) :: x

  pi = 3.141592653589793e+00_real64
  r = 1.0e+00_real64

  if ( x <= 20.0e+00_real64 ) then

    s = 0.5e+00_real64
    do k = 1, 100
 
      if ( k == 1 ) then
        rd = 0.5e+00_real64
      else
        rd = 1.0e+00_real64
      end if
      r = r * rd * k / ( k + 1.0e+00_real64 ) &
        * ( x / ( 2.0e+00_real64 * k + 1.0e+00_real64 ) ) ** 2
      s = s + r
      if ( abs ( r / s ) < 1.0e-12_real64 ) then
        exit
      end if
    end do

    tl0 = 2.0e+00_real64 / pi * x * x * s

  else

    s = 1.0e+00_real64
    do k = 1, 10
      r = r * k / ( k + 1.0e+00_real64 ) &
        * ( ( 2.0e+00_real64 * k + 1.0e+00_real64 ) / x ) ** 2
      s = s + r
      if ( abs ( r / s ) < 1.0e-12_real64 ) then
        exit
      end if
    end do

    el = 0.57721566490153e+00_real64
    s0 = - s / ( pi * x * x ) + 2.0e+00_real64 / pi &
      * ( log ( 2.0e+00_real64 * x ) + el )
    a0 = 1.0e+00_real64
    a1 = 5.0e+00_real64 / 8.0e+00_real64
    a(1) = a1
    do k = 1, 10
      af = ( ( 1.5e+00_real64 * ( k + 0.50e+00_real64 ) &
        * ( k + 5.0e+00_real64 / 6.0e+00_real64 ) * a1 - 0.5e+00_real64 &
        * ( k + 0.5e+00_real64 ) ** 2 * ( k -0.5e+00_real64 ) * a0 ) ) &
        / ( k + 1.0e+00_real64 )
      a(k+1) = af
      a0 = a1
      a1 = af
    end do

    ti = 1.0e+00_real64
    r = 1.0e+00_real64
    do k = 1, 11
      r = r / x
      ti = ti + a(k) * r
    end do
    tl0 = ti / sqrt ( 2.0e+00_real64 * pi * x ) * exp ( x ) + s0

  end if

  return
end subroutine itsl0
!> @brief subroutine itth0.
!> @return None.
!>
!> @param x [in] Argument x.
!> @param tth [inout] Argument tth.
pure subroutine itth0 ( x, tth )

!*****************************************************************************80
!
!! ITTH0 integrates H0(t)/t from x to oo.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    23 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) X, the lower limit of the integral.
!
!    Output, real(real64) TTH, the integral of H0(t)/t from x to oo.
!
  implicit none

  real(real64) f0
  real(real64) g0
  integer(int32) k
  real(real64) pi
  real(real64) r
  real(real64) s
  real(real64) t
  real(real64), intent(inout) :: tth
  real(real64) tty
  real(real64), intent(in) :: x
  real(real64) xt

  pi = 3.141592653589793e+00_real64
  s = 1.0e+00_real64
  r = 1.0e+00_real64

  if ( x < 24.5e+00_real64 ) then

    do k = 1, 60
      r = - r * x * x * ( 2.0e+00_real64 * k - 1.0e+00_real64 ) &
        / ( 2.0e+00_real64 * k + 1.0e+00_real64 ) ** 3
      s = s + r
      if ( abs ( r ) < abs ( s ) * 1.0e-12_real64 ) then
        exit
      end if
    end do

    tth = pi / 2.0e+00_real64 - 2.0e+00_real64 / pi * x * s

  else

    do k = 1, 10
      r = - r * ( 2.0e+00_real64 * k - 1.0e+00_real64 ) ** 3 &
        / ( ( 2.0e+00_real64 * k + 1.0e+00_real64 ) * x * x )
      s = s + r
      if ( abs ( r ) < abs ( s ) * 1.0e-12_real64 ) then
        exit
      end if
    end do

    tth = 2.0e+00_real64 / ( pi * x ) * s
    t = 8.0e+00_real64 / x
    xt = x + 0.25e+00_real64 * pi
    f0 = ((((( &
        0.18118e-02_real64 * t &
      - 0.91909e-02_real64 ) * t &
      + 0.017033e+00_real64 ) * t &
      - 0.9394e-03_real64 ) * t &
      - 0.051445e+00_real64 ) * t &
      - 0.11e-05_real64 ) * t &
      + 0.7978846e+00_real64
    g0 = ((((( &
      - 0.23731e-02_real64 * t &
      + 0.59842e-02_real64 ) * t &
      + 0.24437e-02_real64 ) * t &
      - 0.0233178e+00_real64 ) * t &
      + 0.595e-04_real64 ) * t &
      + 0.1620695e+00_real64 ) * t
    tty = ( f0 * sin ( xt ) - g0 * cos ( xt ) ) / ( sqrt ( x ) * x )
    tth = tth + tty

    end if

  return
end subroutine itth0
!> @brief subroutine ittika.
!> @return None.
!>
!> @param x [in] Argument x.
!> @param tti [inout] Argument tti.
!> @param ttk [inout] Argument ttk.
pure subroutine ittika ( x, tti, ttk )

!*****************************************************************************80
!
!! ITTIKA integrates (I0(t)-1)/t from 0 to x, K0(t)/t from x to infinity.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    23 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) X, the integral limit.
!
!    Output, real(real64) TTI, TTK, the integrals of [I0(t)-1]/t 
!    from 0 to x, and of K0(t)/t from x to oo.
!
  implicit none

  real(real64) b1
  real(real64), parameter, dimension ( 8 ) :: c = [&
    1.625e+00_real64, 4.1328125e+00_real64, &
    1.45380859375e+01_real64, 6.553353881835e+01_real64, &
    3.6066157150269e+02_real64, 2.3448727161884e+03_real64, &
    1.7588273098916e+04_real64, 1.4950639538279e+05_real64]
  real(real64) e0
  real(real64) el
  integer(int32) k
  real(real64) pi
  real(real64) r
  real(real64) r2
  real(real64) rc
  real(real64) rs
  real(real64), intent(inout) :: tti
  real(real64), intent(inout) :: ttk
  real(real64), intent(in) :: x

  pi = 3.141592653589793e+00_real64
  el = 0.5772156649015329e+00_real64

  if ( x == 0.0e+00_real64 ) then
    tti = 0.0e+00_real64
    ttk = 1.0e+300_real64
    return
  end if

  if ( x < 40.0e+00_real64 ) then
    tti = 1.0e+00_real64
    r = 1.0e+00_real64
    do k = 2, 50
      r = 0.25e+00_real64 * r * ( k - 1.0e+00_real64 ) / ( k * k * k ) * x * x
      tti = tti + r
      if ( abs ( r / tti ) < 1.0e-12_real64 ) then
        exit
      end if
    end do

    tti = tti * 0.125e+00_real64 * x * x

  else

    tti = 1.0e+00_real64
    r = 1.0e+00_real64
    do k = 1, 8
      r = r / x
      tti = tti + c(k) * r
    end do
    rc = x * sqrt ( 2.0e+00_real64 * pi * x )
    tti = tti * exp ( x ) / rc

  end if

  if ( x <= 12.0e+00_real64 ) then

    e0 = ( 0.5e+00_real64 * log ( x / 2.0e+00_real64 ) + el ) &
      * log ( x / 2.0e+00_real64 ) + pi * pi / 24.0e+00_real64 + 0.5e+00_real64 * el * el
    b1 = 1.5e+00_real64 - ( el + log ( x / 2.0e+00_real64 ) )
    rs = 1.0e+00_real64
    r = 1.0e+00_real64
    do k = 2, 50
      r = 0.25e+00_real64 * r * ( k - 1.0e+00_real64 ) / ( k * k * k ) * x * x
      rs = rs + 1.0e+00_real64 / k
      r2 = r * ( rs + 1.0e+00_real64 / ( 2.0e+00_real64 * k ) &
        - ( el + log ( x / 2.0e+00_real64 ) ) )
      b1 = b1 + r2
      if ( abs ( r2 / b1 ) < 1.0e-12_real64 ) then
        exit
      end if
    end do

    ttk = e0 - 0.125e+00_real64 * x * x * b1

  else

    ttk = 1.0e+00_real64
    r = 1.0e+00_real64
    do k = 1, 8
      r = - r / x
      ttk = ttk + c(k) * r
    end do
    rc = x * sqrt ( 2.0e+00_real64 / pi * x )
    ttk = ttk * exp ( - x ) / rc

  end if

  return
end subroutine ittika
!> @brief subroutine ittikb.
!> @return None.
!>
!> @param x [in] Argument x.
!> @param tti [inout] Argument tti.
!> @param ttk [inout] Argument ttk.
pure subroutine ittikb ( x, tti, ttk )

!*****************************************************************************80
!
!! ITTIKB integrates (I0(t)-1)/t from 0 to x, K0(t)/t from x to infinity.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    28 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) X, the integral limit.
!
!    Output, real(real64) TTI, TTK, the integrals of
!    [I0(t)-1]/t from 0 to x, and K0(t)/t from x to oo.
!
  implicit none

  real(real64) e0
  real(real64) el
  real(real64) pi
  real(real64) t
  real(real64) t1
  real(real64), intent(inout) :: tti
  real(real64), intent(inout) :: ttk
  real(real64), intent(in) :: x
  real(real64) x1

  pi = 3.141592653589793e+00_real64
  el = 0.5772156649015329e+00_real64

  if ( x == 0.0e+00_real64 ) then

    tti = 0.0e+00_real64

  else if ( x <= 5.0e+00_real64 ) then

    x1 = x / 5.0e+00_real64
    t = x1 * x1
    tti = ((((((( &
        0.1263e-03_real64       * t &
      + 0.96442e-03_real64 )    * t &
      + 0.968217e-02_real64 )   * t &
      + 0.06615507e+00_real64 ) * t &
      + 0.33116853e+00_real64 ) * t &
      + 1.13027241e+00_real64 ) * t &
      + 2.44140746e+00_real64 ) * t &
      + 3.12499991e+00_real64 ) * t

  else

    t = 5.0e+00_real64 / x
    tti = ((((((((( &
         2.1945464e+00_real64   * t &
      -  3.5195009e+00_real64 ) * t &
      - 11.9094395e+00_real64 ) * t &
      + 40.394734e+00_real64  ) * t &
      - 48.0524115e+00_real64 ) * t &
      + 28.1221478e+00_real64 ) * t &
      -  8.6556013e+00_real64 ) * t &
      +  1.4780044e+00_real64 ) * t &
      -  0.0493843e+00_real64 ) * t &
      +  0.1332055e+00_real64 ) * t &
      +  0.3989314e+00_real64
    tti = tti * exp ( x ) / ( sqrt ( x ) * x )

  end if

  if ( x == 0.0e+00_real64 ) then

    ttk = 1.0e+300_real64

  else if ( x <= 2.0e+00_real64 ) then

    t1 = x / 2.0e+00_real64
    t = t1 * t1
    ttk = ((((( &
        0.77e-06_real64         * t &
      + 0.1544e-04_real64 )     * t &
      + 0.48077e-03_real64 )    * t &
      + 0.925821e-02_real64 )   * t &
      + 0.10937537e+00_real64 ) * t &
      + 0.74999993e+00_real64 ) * t
    e0 = el + log ( x / 2.0e+00_real64 )
    ttk = pi * pi / 24.0e+00_real64 + e0 * ( 0.5e+00_real64 * e0 + tti ) - ttk

  else if ( x <= 4.0e+00_real64 ) then

    t = 2.0e+00_real64 / x
    ttk = ((( &
        0.06084e+00_real64    * t &
      - 0.280367e+00_real64 ) * t &
      + 0.590944e+00_real64 ) * t &
      - 0.850013e+00_real64 ) * t &
      + 1.234684e+00_real64
    ttk = ttk * exp ( - x ) / ( sqrt ( x ) * x )

  else

    t = 4.0e+00_real64 / x
    ttk = ((((( &
        0.02724e+00_real64     * t &
      - 0.1110396e+00_real64 ) * t &
      + 0.2060126e+00_real64 ) * t &
      - 0.2621446e+00_real64 ) * t &
      + 0.3219184e+00_real64 ) * t &
      - 0.5091339e+00_real64 ) * t &
      + 1.2533141e+00_real64
    ttk = ttk * exp ( - x ) / ( sqrt ( x ) * x )

  end if

  return
end subroutine ittikb
!> @brief subroutine ittjya.
!> @return None.
!>
!> @param x [in] Argument x.
!> @param ttj [inout] Argument ttj.
!> @param tty [inout] Argument tty.
pure subroutine ittjya ( x, ttj, tty )

!*****************************************************************************80
!
!! ITTJYA integrates (1-J0(t))/t from 0 to x, and Y0(t)/t from x to infinity.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    28 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) X, the integral limit.
!
!    Output, real(real64) TTJ, TTY, the integrals of [1-J0(t)]/t 
!    from 0 to x and of Y0(t)/t from x to oo.
!
  implicit none

  real(real64) a0
  real(real64) b1
  real(real64) bj0
  real(real64) bj1
  real(real64) by0
  real(real64) by1
  real(real64) e0
  real(real64) el
  real(real64) g0
  real(real64) g1
  integer(int32) k
  integer(int32) l
  real(real64) pi
  real(real64) px
  real(real64) qx
  real(real64) r
  real(real64) r0
  real(real64) r1
  real(real64) r2
  real(real64) rs
  real(real64) t
  real(real64), intent(inout) :: ttj
  real(real64), intent(inout) :: tty
  real(real64) vt
  real(real64), intent(in) :: x
  real(real64) xk

  pi = 3.141592653589793e+00_real64
  el = 0.5772156649015329e+00_real64

  if ( x == 0.0e+00_real64 ) then

    ttj = 0.0e+00_real64
    tty = -1.0e+300_real64

  else if ( x <= 20.0e+00_real64 ) then

    ttj = 1.0e+00_real64
    r = 1.0e+00_real64
    do k = 2, 100
      r = - 0.25e+00_real64 * r * ( k - 1.0e+00_real64 ) / ( k * k * k ) * x * x
      ttj = ttj + r
      if ( abs ( r ) < abs ( ttj ) * 1.0e-12_real64 ) then
        exit
      end if
    end do

    ttj = ttj * 0.125e+00_real64 * x * x
    e0 = 0.5e+00_real64 * ( pi * pi / 6.0e+00_real64 - el * el ) &
      - ( 0.5e+00_real64 * log ( x / 2.0e+00_real64 ) + el ) &
      * log ( x / 2.0e+00_real64 )
    b1 = el + log ( x / 2.0e+00_real64 ) - 1.5e+00_real64
    rs = 1.0e+00_real64
    r = -1.0e+00_real64
    do k = 2, 100
      r = - 0.25e+00_real64 * r * ( k - 1.0e+00_real64 ) / ( k * k * k ) * x * x
      rs = rs + 1.0e+00_real64 / k
      r2 = r * ( rs + 1.0e+00_real64 / ( 2.0e+00_real64 * k ) &
        - ( el + log ( x / 2.0e+00_real64 ) ) ) 
      b1 = b1 + r2
      if ( abs ( r2 ) < abs ( b1 ) * 1.0e-12_real64 ) then
        exit
      end if
    end do

    tty = 2.0e+00_real64 / pi * ( e0 + 0.125e+00_real64 * x * x * b1 )

  else

    a0 = sqrt ( 2.0e+00_real64 / ( pi * x ) )

    do l = 0, 1

      vt = 4.0e+00_real64 * l * l
      px = 1.0e+00_real64
      r = 1.0e+00_real64
      do k = 1, 14
        r = - 0.0078125e+00_real64 * r &
          * ( vt - ( 4.0e+00_real64 * k - 3.0e+00_real64 ) ** 2 ) &
          / ( x * k ) * ( vt - ( 4.0e+00_real64 * k - 1.0e+00_real64 ) ** 2 ) &
          / ( ( 2.0e+00_real64 * k - 1.0e+00_real64 ) * x )
        px = px + r
        if ( abs ( r ) < abs ( px ) * 1.0e-12_real64 ) then
          exit
        end if
      end do

      qx = 1.0e+00_real64
      r = 1.0e+00_real64
      do k = 1, 14
        r = -0.0078125e+00_real64 * r &
          * ( vt - ( 4.0e+00_real64 * k - 1.0e+00_real64 ) ** 2 ) &
          / ( x * k ) * ( vt - ( 4.0e+00_real64 * k + 1.0e+00_real64 ) ** 2 ) &
          / ( 2.0e+00_real64 * k + 1.0e+00_real64 ) / x
        qx = qx + r
        if ( abs ( r ) < abs ( qx ) * 1.0e-12_real64 ) then
          exit
        end if
      end do

      qx = 0.125e+00_real64 * ( vt - 1.0e+00_real64 ) / x * qx
      xk = x - ( 0.25e+00_real64 + 0.5e+00_real64 * l ) * pi
      bj1 = a0 * ( px * cos ( xk ) - qx * sin ( xk ) )
      by1 = a0 * ( px * sin ( xk ) + qx * cos ( xk ) )
      if ( l == 0 ) then
        bj0 = bj1
        by0 = by1
      end if

    end do

    t = 2.0e+00_real64 / x
    g0 = 1.0e+00_real64
    r0 = 1.0e+00_real64
    do k = 1, 10
      r0 = - k * k * t * t *r0
      g0 = g0 + r0
    end do

    g1 = 1.0e+00_real64
    r1 = 1.0e+00_real64
    do k = 1, 10
      r1 = - k * ( k + 1.0e+00_real64 ) * t * t * r1
      g1 = g1 + r1
    end do

    ttj = 2.0e+00_real64 * g1 * bj0 / ( x * x ) - g0 * bj1 / x &
      + el + log ( x / 2.0e+00_real64 )
    tty = 2.0e+00_real64 * g1 * by0 / ( x * x ) - g0 * by1 / x

  end if

  return
end subroutine ittjya
!> @brief subroutine ittjyb.
!> @return None.
!>
!> @param x [in] Argument x.
!> @param ttj [inout] Argument ttj.
!> @param tty [inout] Argument tty.
pure subroutine ittjyb ( x, ttj, tty )

!*****************************************************************************80
!
!! ITTJYB integrates (1-J0(t))/t from 0 to x, and Y0(t)/t from x to infinity.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    01 August 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) X, the integral limit.
!
!    Output, real(real64) TTJ, TTY, the integrals of [1-J0(t)]/t 
!    from 0 to x and of Y0(t)/t from x to oo.
!
  implicit none

  real(real64) e0
  real(real64) el
  real(real64) f0
  real(real64) g0
  real(real64) pi
  real(real64) t
  real(real64) t1
  real(real64), intent(inout) :: ttj
  real(real64), intent(inout) :: tty
  real(real64), intent(in) :: x
  real(real64) x1
  real(real64) xt

  pi = 3.141592653589793e+00_real64
  el = 0.5772156649015329e+00_real64

  if ( x == 0.0e+00_real64 ) then

    ttj = 0.0e+00_real64
    tty = -1.0e+300_real64

  else if ( x <= 4.0e+00_real64 ) then

    x1 = x / 4.0e+00_real64
    t = x1 * x1

    ttj = (((((( &
        0.35817e-04_real64 * t &
      - 0.639765e-03_real64 ) * t &
      + 0.7092535e-02_real64 ) * t &
      - 0.055544803e+00_real64 ) * t &
      + 0.296292677e+00_real64 ) * t &
      - 0.999999326e+00_real64 ) * t &
      + 1.999999936e+00_real64 ) * t
 
    tty = ((((((( &
      - 0.3546e-05_real64        * t &
      + 0.76217e-04_real64 )     * t &
      - 0.1059499e-02_real64 )   * t &
      + 0.010787555e+00_real64 ) * t &
      - 0.07810271e+00_real64 )  * t &
      + 0.377255736e+00_real64 ) * t &
      - 1.114084491e+00_real64 ) * t &
      + 1.909859297e+00_real64 ) * t

    e0 = el + log ( x / 2.0e+00_real64 )
    tty = pi / 6.0e+00_real64 + e0 / pi * ( 2.0e+00_real64 * ttj - e0 ) - tty

  else if ( x <= 8.0e+00_real64 ) then

    xt = x + 0.25e+00_real64 * pi
    t1 = 4.0e+00_real64 / x
    t = t1 * t1

    f0 = ((((( &
        0.0145369e+00_real64 * t &
      - 0.0666297e+00_real64 ) * t &
      + 0.1341551e+00_real64 ) * t &
      - 0.1647797e+00_real64 ) * t &
      + 0.1608874e+00_real64 ) * t &
      - 0.2021547e+00_real64 ) * t &
      + 0.7977506e+00_real64

    g0 = (((((( &
        0.0160672e+00_real64   * t &
      - 0.0759339e+00_real64 ) * t &
      + 0.1576116e+00_real64 ) * t &
      - 0.1960154e+00_real64 ) * t &
      + 0.1797457e+00_real64 ) * t &
      - 0.1702778e+00_real64 ) * t &
      + 0.3235819e+00_real64 ) * t1

    ttj = ( f0 * cos ( xt ) + g0 * sin ( xt ) ) / ( sqrt ( x ) * x )
    ttj = ttj + el + log ( x / 2.0e+00_real64 )
    tty = ( f0 * sin ( xt ) - g0 * cos ( xt ) ) / ( sqrt ( x ) * x )

  else

    t = 8.0e+00_real64 / x
    xt = x + 0.25e+00_real64 * pi

    f0 = ((((( &
        0.18118e-02_real64    * t &
      - 0.91909e-02_real64 )  * t &
      + 0.017033e+00_real64 ) * t &
      - 0.9394e-03_real64 )   * t &
      - 0.051445e+00_real64 ) * t &
      - 0.11e-05_real64 )     * t &
      + 0.7978846e+00_real64

    g0 = ((((( &
      - 0.23731e-02_real64     * t &
      + 0.59842e-02_real64 )   * t &
      + 0.24437e-02_real64 )   * t &
      - 0.0233178e+00_real64 ) * t &
      + 0.595e-04_real64 )     * t &
      + 0.1620695e+00_real64 ) * t

    ttj = ( f0 * cos ( xt ) + g0 * sin ( xt ) )  &
      / ( sqrt ( x ) * x ) + el + log ( x / 2.0e+00_real64 )
    tty = ( f0 * sin ( xt ) - g0 * cos ( xt ) )  &
      / ( sqrt ( x ) * x )

  end if

  return
end subroutine ittjyb
!> @brief subroutine jdzo.
!> @return None.
!>
!> @param nt [in] Argument nt.
!> @param n [inout] Argument n.
!> @param m [inout] Argument m.
!> @param p [inout] Argument p.
!> @param zo [inout] Argument zo.
subroutine jdzo ( nt, n, m, p, zo )

!*****************************************************************************80
!
!! JDZO computes the zeros of Bessel functions Jn(x) and Jn'(x).
!
!  Discussion:
!
!    This procedure computes the zeros of Bessel functions Jn(x) and
!    Jn'(x), and arrange them in the order of their magnitudes.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    01 August 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) NT, the number of zeros.
!
!    Output, integer(int32) N(*), the  order of Jn(x) or Jn'(x) associated
!    with the L-th zero.
!
!    Output, integer(int32) M(*), the serial number of the zeros of Jn(x)
!    or Jn'(x) associated with the L-th zero ( L is the serial number of all the
!    zeros of Jn(x) and Jn'(x) ).
!
!    Output, character ( len = 4 ) P(L), 'TM' or 'TE', a code for designating 
!    the zeros of Jn(x)  or Jn'(x).  In the waveguide applications, the zeros
!    of Jn(x) correspond to TM modes and those of Jn'(x) correspond to TE modes.
!
!    Output, real(real64) ZO(*), the zeros of Jn(x) and Jn'(x).
!
  implicit none

  real(real64) bj(101)
  real(real64) dj(101)
  real(real64) fj(101)
  integer(int32) i
  integer(int32) j
  integer(int32) k
  integer(int32) l
  integer(int32) l0
  integer(int32) l1
  integer(int32) l2
  integer(int32), intent(inout) :: m(1400)
  integer(int32) m1(70)
  integer(int32) mm
  integer(int32), intent(inout) :: n(1400)
  integer(int32) n1(70)
  integer(int32) nm
  integer(int32), intent(in) :: nt
  character ( len = 4 ), intent(inout) :: p(1400)
  character ( len = 4 ) p1(70)
  real(real64) x
  real(real64) x0
  real(real64) x1
  real(real64) x2
  real(real64) xm
  real(real64), intent(inout) :: zo(1400)
  real(real64) zoc(70)

  if ( nt < 600 ) then
    xm = -1.0e+00_real64 + 2.248485e+00_real64 * real ( nt, kind = real64 ) ** 0.5e+00_real64 &
      - 0.0159382e+00_real64 * nt + 3.208775e-04_real64 * real ( nt, kind = real64 ) ** 1.5e+00_real64
    nm = int ( 14.5e+00_real64 + 0.05875e+00_real64 * nt )
    mm = int ( 0.02e+00_real64 * nt ) + 6
  else
    xm = 5.0e+00_real64 + 1.445389e+00_real64 * ( real ( nt, kind = real64 ) ) ** 0.5e+00_real64 &
      + 0.01889876e+00_real64 * nt &
      - 2.147763e-04_real64 * ( real ( nt, kind = real64 ) ) ** 1.5e+00_real64
    nm = int ( 27.8e+00_real64 + 0.0327e+00_real64 * nt )
    mm = int ( 0.01088e+00_real64 * nt ) + 10
  end if

  l0 = 0

  do i = 1,nm

    x1 = 0.407658e+00_real64 + 0.4795504e+00_real64 &
      * ( real ( i - 1, kind = real64 ) ) ** 0.5e+00_real64 + 0.983618e+00_real64 * ( i - 1 )
    x2 = 1.99535e+00_real64 + 0.8333883_real64 * ( real ( i - 1, kind = real64 ) ) ** 0.5e+00_real64 &
      + 0.984584e+00_real64 * ( i - 1 )
    l1 = 0

    do j = 1, mm

      if ( i == 1 .and. j == 1 ) then

        l1 = l1 + 1
        n1(l1) = i - 1
        m1(l1) = j
        if ( i == 1 ) then
          m1(l1) = j - 1
        end if
        p1(l1) = 'TE'
        zoc(l1) = x

        if ( i <= 15 ) then
          x1 = x + 3.057e+00_real64 + 0.0122e+00_real64 * ( i - 1 ) &
            + ( 1.555e+00_real64 + 0.41575e+00_real64 * ( i - 1 ) ) / ( j + 1 ) ** 2
        else
          x1 = x + 2.918e+00_real64 + 0.01924e+00_real64 * ( i - 1 ) &
            + ( 6.26e+00_real64 + 0.13205e+00_real64 * ( i - 1 ) ) / ( j + 1 ) ** 2
        end if

      else

        x = x1

        do

          call bjndd ( i, x, bj, dj, fj )
          x0 = x
          x = x - dj(i) / fj(i)

          if ( xm < x1 ) then
            exit
          end if

          if ( abs ( x - x0 ) <= 1.0e-10_real64 ) then
            l1 = l1 + 1
            n1(l1) = i - 1
            m1(l1) = j
            if ( i == 1 ) then
              m1(l1) = j - 1
            end if
            p1(l1) = 'TE'
            zoc(l1) = x

            if ( i <= 15 ) then
              x1 = x + 3.057e+00_real64 + 0.0122e+00_real64 * ( i - 1 ) &
                + ( 1.555e+00_real64 + 0.41575e+00_real64 * ( i - 1 ) ) / ( j + 1 ) ** 2
            else
              x1 = x + 2.918e+00_real64 + 0.01924e+00_real64 * ( i - 1 ) &
                + ( 6.26e+00_real64 + 0.13205e+00_real64 * ( i - 1 ) ) / ( j + 1 ) ** 2
            end if
            exit
          end if

        end do

      end if

      x = x2

      do

        call bjndd ( i, x, bj, dj, fj )
        x0 = x
        x = x - bj(i) / dj(i)

        if ( xm < x ) then
          exit
        end if

        if ( abs ( x - x0 ) <= 1.0e-10_real64 ) then
          exit
        end if

      end do

      if ( x <= xm ) then

        l1 = l1 + 1
        n1(l1) = i - 1
        m1(l1) = j
        p1(l1) = 'TM'
        zoc(l1) = x
        if ( i <= 15 ) then
          x2 = x + 3.11e+00_real64 + 0.0138e+00_real64 * ( i - 1 ) &
            + ( 0.04832e+00_real64 + 0.2804e+00_real64 * ( i - 1 ) ) / ( j + 1 ) ** 2
        else
          x2 = x + 3.001e+00_real64 + 0.0105e+00_real64 * ( i - 1 ) &
            + ( 11.52e+00_real64 + 0.48525e+00_real64 * ( i - 1 ) ) / ( j + 3 ) ** 2
        end if

      end if

    end do

    l = l0 + l1
    l2 = l

    do

      if ( l0 == 0 ) then
        do k = 1, l
          zo(k) = zoc(k)
          n(k) = n1(k)
          m(k) = m1(k)
          p(k) = p1(k)
        end do
        l1 = 0
      else if ( l0 /= 0 ) then
        if ( zoc(l1) .le. zo(l0) ) then
          zo(l0+l1) = zo(l0)
          n(l0+l1) = n(l0)
          m(l0+l1) = m(l0)
          p(l0+l1) = p(l0)
          l0 = l0 - 1
        else
          zo(l0+l1) = zoc(l1)
          n(l0+l1) = n1(l1)
          m(l0+l1) = m1(l1)
          p(l0+l1) = p1(l1)
          l1 = l1 - 1
        end if
      end if

      if ( l1 == 0 ) then
        exit 
      end if

    end do

    l0 = l2

  end do

  return
end subroutine jdzo
!> @brief subroutine jelp.
!> @return None.
!>
!> @param u [in] Argument u.
!> @param hk [in] Argument hk.
!> @param esn [inout] Argument esn.
!> @param ecn [inout] Argument ecn.
!> @param edn [inout] Argument edn.
!> @param eph [inout] Argument eph.
pure subroutine jelp ( u, hk, esn, ecn, edn, eph )

!*****************************************************************************80
!
!! JELP computes Jacobian elliptic functions SN(u), CN(u), DN(u).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    08 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) U, the argument.
!
!    Input, real(real64) HK, the modulus, between 0 and 1._real64
!
!    Output, real(real64) ESN, ECN, EDN, EPH, the values of
!    sn(u), cn(u), dn(u), and phi (in degrees).
!
  implicit none

  real(real64) a
  real(real64) a0
  real(real64) b
  real(real64) b0
  real(real64) c
  real(real64) d
  real(real64) dn
  real(real64), intent(inout) :: ecn
  real(real64), intent(inout) :: edn
  real(real64), intent(inout) :: eph
  real(real64), intent(inout) :: esn
  real(real64), intent(in) :: hk
  integer(int32) j
  integer(int32) n 
  real(real64) pi
  real(real64) r(40)
  real(real64) sa
  real(real64) t
  real(real64), intent(in) :: u

  pi = 3.14159265358979e+00_real64
  a0 = 1.0e+00_real64
  b0 = sqrt ( 1.0e+00_real64 - hk * hk )

  do n = 1, 40

    a = ( a0 + b0 ) / 2.0e+00_real64
    b = sqrt ( a0 * b0 )
    c = ( a0 - b0 ) / 2.0e+00_real64
    r(n) = c / a

    if ( c < 1.0e-07_real64 ) then
      exit
    end if

    a0 = a
    b0 = b

  end do

  dn = 2.0e+00_real64 ** n * a * u

  do j = n, 1, -1
    t = r(j) * sin ( dn )
    sa = atan ( t / sqrt ( abs ( 1.0e+00_real64 - t * t )))
    d = 0.5e+00_real64 * ( dn + sa )
    dn = d
  end do

  eph = d * 180.0e+00_real64 / pi
  esn = sin ( d )
  ecn = cos ( d )
  edn = sqrt ( 1.0e+00_real64 - hk * hk * esn * esn )

  return
end subroutine jelp
!> @brief subroutine jy01a.
!> @return None.
!>
!> @param x [in] Argument x.
!> @param bj0 [inout] Argument bj0.
!> @param dj0 [inout] Argument dj0.
!> @param bj1 [inout] Argument bj1.
!> @param dj1 [inout] Argument dj1.
!> @param by0 [inout] Argument by0.
!> @param dy0 [inout] Argument dy0.
!> @param by1 [inout] Argument by1.
!> @param dy1 [inout] Argument dy1.
pure subroutine jy01a ( x, bj0, dj0, bj1, dj1, by0, dy0, by1, dy1 )

!*****************************************************************************80
!
!! JY01A computes Bessel functions J0(x), J1(x), Y0(x), Y1(x) and derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    01 August 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) X, the argument.
!
!    Output, real(real64) BJ0, DJ0, BJ1, DJ1, BY0, DY0, BY1, DY1,
!    the values of J0(x), J0'(x), J1(x), J1'(x), Y0(x), Y0'(x), Y1(x), Y1'(x).
!
  implicit none

  real(real64), parameter, dimension ( 12 ) :: a = [&
    -0.7031250000000000e-01_real64, 0.1121520996093750e+00_real64, &
    -0.5725014209747314e+00_real64, 0.6074042001273483e+01_real64, &
    -0.1100171402692467e+03_real64, 0.3038090510922384e+04_real64, &
    -0.1188384262567832e+06_real64, 0.6252951493434797e+07_real64, &
    -0.4259392165047669e+09_real64, 0.3646840080706556e+11_real64, &
    -0.3833534661393944e+13_real64, 0.4854014686852901e+15_real64]
  real(real64), parameter, dimension ( 12 ) :: a1 = [&
    0.1171875000000000e+00_real64, -0.1441955566406250e+00_real64, &
    0.6765925884246826e+00_real64, -0.6883914268109947e+01_real64, &
    0.1215978918765359e+03_real64, -0.3302272294480852e+04_real64, &
    0.1276412726461746e+06_real64, -0.6656367718817688e+07_real64, &
    0.4502786003050393e+09_real64, -0.3833857520742790e+11_real64, &
    0.4011838599133198e+13_real64, -0.5060568503314727e+15_real64]
  real(real64), parameter, dimension ( 12 ) :: b = [&
    0.7324218750000000e-01_real64, -0.2271080017089844e+00_real64, &
    0.1727727502584457e+01_real64, -0.2438052969955606e+02_real64, &
    0.5513358961220206e+03_real64, -0.1825775547429318e+05_real64, &
    0.8328593040162893e+06_real64, -0.5006958953198893e+08_real64, &
    0.3836255180230433e+10_real64, -0.3649010818849833e+12_real64, &
    0.4218971570284096e+14_real64, -0.5827244631566907e+16_real64]
  real(real64), parameter, dimension ( 12 ) :: b1 = [&
    -0.1025390625000000e+00_real64, 0.2775764465332031e+00_real64, &
    -0.1993531733751297e+01_real64, 0.2724882731126854e+02_real64, &
    -0.6038440767050702e+03_real64, 0.1971837591223663e+05_real64, &
    -0.8902978767070678e+06_real64, 0.5310411010968522e+08_real64, &
    -0.4043620325107754e+10_real64, 0.3827011346598605e+12_real64, &
    -0.4406481417852278e+14_real64, 0.6065091351222699e+16_real64]
  real(real64), intent(inout) :: bj0
  real(real64), intent(inout) :: bj1
  real(real64), intent(inout) :: by0
  real(real64), intent(inout) :: by1
  real(real64) cs0
  real(real64) cs1
  real(real64) cu
  real(real64), intent(inout) :: dj0
  real(real64), intent(inout) :: dj1
  real(real64), intent(inout) :: dy0
  real(real64), intent(inout) :: dy1
  real(real64) ec
  integer(int32) k
  integer(int32) k0
  real(real64) p0
  real(real64) p1
  real(real64) pi
  real(real64) q0
  real(real64) q1
  real(real64) r
  real(real64) r0
  real(real64) r1
  real(real64) rp2
  real(real64) t1
  real(real64) t2
  real(real64) w0
  real(real64) w1
  real(real64), intent(in) :: x
  real(real64) x2

  pi = 3.141592653589793e+00_real64
  rp2 = 0.63661977236758e+00_real64
  x2 = x * x

  if ( x == 0.0e+00_real64 ) then
    bj0 = 1.0e+00_real64
    bj1 = 0.0e+00_real64
    dj0 = 0.0e+00_real64
    dj1 = 0.5e+00_real64
    by0 = -1.0e+300_real64
    by1 = -1.0e+300_real64
    dy0 = 1.0e+300_real64
    dy1 = 1.0e+300_real64
    return
  end if

  if ( x <= 12.0e+00_real64 ) then

    bj0 = 1.0e+00_real64
    r = 1.0e+00_real64
    do k = 1,30
      r = -0.25e+00_real64 * r * x2 / ( k * k )
      bj0 = bj0 + r
      if ( abs ( r ) < abs ( bj0 ) * 1.0e-15_real64 ) then
        exit
      end if
    end do

    bj1 = 1.0e+00_real64
    r = 1.0e+00_real64
    do k = 1, 30
      r = -0.25e+00_real64 * r * x2 / ( k * ( k + 1.0e+00_real64 ) )
      bj1 = bj1 + r
      if ( abs ( r ) < abs ( bj1 ) * 1.0e-15_real64 ) then
        exit
      end if
    end do

    bj1 = 0.5e+00_real64 * x * bj1
    ec = log ( x / 2.0e+00_real64 ) + 0.5772156649015329e+00_real64
    cs0 = 0.0e+00_real64
    w0 = 0.0e+00_real64
    r0 = 1.0e+00_real64
    do k = 1, 30
      w0 = w0 + 1.0e+00_real64 / k
      r0 = -0.25e+00_real64 * r0 / ( k * k ) * x2
      r = r0 * w0
      cs0 = cs0 + r
      if ( abs ( r ) < abs ( cs0 ) * 1.0e-15_real64 ) then
        exit
      end if
    end do

    by0 = rp2 * ( ec * bj0 - cs0 )
    cs1 = 1.0e+00_real64
    w1 = 0.0e+00_real64
    r1 = 1.0e+00_real64
    do k = 1, 30
      w1 = w1 + 1.0e+00_real64 / k
      r1 = -0.25e+00_real64 * r1 / ( k * ( k + 1 ) ) * x2
      r = r1 * ( 2.0e+00_real64 * w1 + 1.0e+00_real64 / ( k + 1.0e+00_real64 ) )
      cs1 = cs1 + r
      if ( abs ( r ) < abs ( cs1 ) * 1.0e-15_real64 ) then
        exit
      end if
    end do

    by1 = rp2 * ( ec * bj1 - 1.0e+00_real64 / x - 0.25e+00_real64 * x * cs1 )

  else

    if ( x < 35.0e+00_real64 ) then
      k0 = 12
    else if ( x < 50.0e+00_real64 ) then
      k0 = 10
    else
      k0 = 8
    end if

    t1 = x - 0.25e+00_real64 * pi
    p0 = 1.0e+00_real64
    q0 = -0.125e+00_real64 / x
    do k = 1, k0
      p0 = p0 + a(k) * x ** ( - 2 * k )
      q0 = q0 + b(k) * x ** ( - 2 * k - 1 )
    end do
    cu = sqrt ( rp2 / x )
    bj0 = cu * ( p0 * cos ( t1 ) - q0 * sin ( t1 ) )
    by0 = cu * ( p0 * sin ( t1 ) + q0 * cos ( t1 ) )
    t2 = x - 0.75e+00_real64 * pi
    p1 = 1.0e+00_real64
    q1 = 0.375e+00_real64 / x
    do k = 1, k0
      p1 = p1 + a1(k) * x ** ( - 2 * k )
      q1 = q1 + b1(k) * x ** ( - 2 * k - 1 )
    end do
    cu = sqrt ( rp2 / x )
    bj1 = cu * ( p1 * cos ( t2 ) - q1 * sin ( t2 ) )
    by1 = cu * ( p1 * sin ( t2 ) + q1 * cos ( t2 ) )

  end if

  dj0 = - bj1
  dj1 = bj0 - bj1 / x
  dy0 = - by1
  dy1 = by0 - by1 / x

  return
end subroutine jy01a
!> @brief subroutine jy01b.
!> @return None.
!>
!> @param x [in] Argument x.
!> @param bj0 [inout] Argument bj0.
!> @param dj0 [inout] Argument dj0.
!> @param bj1 [inout] Argument bj1.
!> @param dj1 [inout] Argument dj1.
!> @param by0 [inout] Argument by0.
!> @param dy0 [inout] Argument dy0.
!> @param by1 [inout] Argument by1.
!> @param dy1 [inout] Argument dy1.
pure subroutine jy01b ( x, bj0, dj0, bj1, dj1, by0, dy0, by1, dy1 )

!*****************************************************************************80
!
!! JY01B computes Bessel functions J0(x), J1(x), Y0(x), Y1(x) and derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    02 August 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) X, the argument.
!
!    Output, real(real64) BJ0, DJ0, BJ1, DJ1, BY0, DY0, BY1, DY1,
!    the values of J0(x), J0'(x), J1(x), J1'(x), Y0(x), Y0'(x), Y1(x), Y1'(x).
!
  implicit none

  real(real64) a0
  real(real64), intent(inout) :: bj0
  real(real64), intent(inout) :: bj1
  real(real64), intent(inout) :: by0
  real(real64), intent(inout) :: by1
  real(real64), intent(inout) :: dj0
  real(real64), intent(inout) :: dj1
  real(real64), intent(inout) :: dy0
  real(real64), intent(inout) :: dy1
  real(real64) p0
  real(real64) p1
  real(real64) pi
  real(real64) q0
  real(real64) q1
  real(real64) t
  real(real64) t2
  real(real64) ta0
  real(real64) ta1
  real(real64), intent(in) :: x

  pi = 3.141592653589793e+00_real64

  if ( x == 0.0e+00_real64 ) then

    bj0 = 1.0e+00_real64
    bj1 = 0.0e+00_real64
    dj0 = 0.0e+00_real64
    dj1 = 0.5e+00_real64
    by0 = -1.0e+300_real64
    by1 = -1.0e+300_real64
    dy0 = 1.0e+300_real64
    dy1 = 1.0e+300_real64
    return

  else if ( x <= 4.0e+00_real64 ) then

    t = x / 4.0e+00_real64
    t2 = t * t

    bj0 = (((((( &
      - 0.5014415e-03_real64 * t2 &
      + 0.76771853e-02_real64 ) * t2 &
      - 0.0709253492e+00_real64 ) * t2 &
      + 0.4443584263e+00_real64 ) * t2 &
      - 1.7777560599e+00_real64 ) * t2 &
      + 3.9999973021e+00_real64 ) * t2 &
      - 3.9999998721e+00_real64 ) * t2 &
      + 1.0e+00_real64

    bj1 = t * ((((((( &
      - 0.1289769e-03_real64 * t2 &
      + 0.22069155e-02_real64 ) * t2 &
      - 0.0236616773e+00_real64 ) * t2 &
      + 0.1777582922e+00_real64 ) * t2 &
      - 0.8888839649e+00_real64 ) * t2 &
      + 2.6666660544e+00_real64 ) * t2 &
      - 3.9999999710e+00_real64 ) * t2 &
      + 1.9999999998e+00_real64 )

    by0 = ((((((( &
      - 0.567433e-04_real64 * t2 &
      + 0.859977e-03_real64 ) * t2 &
      - 0.94855882e-02_real64 ) * t2 &
      + 0.0772975809e+00_real64 ) * t2 &
      - 0.4261737419e+00_real64 ) * t2 &
      + 1.4216421221e+00_real64 ) * t2 &
      - 2.3498519931e+00_real64 ) * t2 &
      + 1.0766115157e+00_real64 ) * t2 &
      + 0.3674669052e+00_real64

    by0 = 2.0e+00_real64 / pi * log ( x / 2.0e+00_real64 ) * bj0 + by0

    by1 = (((((((( &
        0.6535773e-03_real64 * t2 &
      - 0.0108175626e+00_real64 ) * t2 &
      + 0.107657606e+00_real64 ) * t2 &
      - 0.7268945577e+00_real64 ) * t2 &
      + 3.1261399273e+00_real64 ) * t2 &
      - 7.3980241381e+00_real64 ) * t2 &
      + 6.8529236342e+00_real64 ) * t2 &
      + 0.3932562018e+00_real64 ) * t2 &
      - 0.6366197726e+00_real64 ) / x

    by1 = 2.0e+00_real64 / pi * log ( x / 2.0e+00_real64 ) * bj1 + by1

  else

    t = 4.0e+00_real64 / x
    t2 = t * t
    a0 = sqrt ( 2.0e+00_real64 / ( pi * x ) )

    p0 = (((( &
      - 0.9285e-05_real64 * t2 &
      + 0.43506e-04_real64 ) * t2 &
      - 0.122226e-03_real64 ) * t2 &
      + 0.434725e-03_real64 ) * t2 &
      - 0.4394275e-02_real64 ) * t2 &
      + 0.999999997e+00_real64

    q0 = t * ((((( &
        0.8099e-05_real64 * t2 &
      - 0.35614e-04_real64 ) * t2 &
      + 0.85844e-04_real64 ) * t2 &
      - 0.218024e-03_real64 ) * t2 &
      + 0.1144106e-02_real64 ) * t2 &
      - 0.031249995e+00_real64 )

    ta0 = x - 0.25e+00_real64 * pi
    bj0 = a0 * ( p0 * cos ( ta0 ) - q0 * sin ( ta0 ) )
    by0 = a0 * ( p0 * sin ( ta0 ) + q0 * cos ( ta0 ) )

    p1 = (((( &
        0.10632e-04_real64 * t2 &
      - 0.50363e-04_real64 ) * t2 &
      + 0.145575e-03_real64 ) * t2 &
      - 0.559487e-03_real64 ) * t2 &
      + 0.7323931e-02_real64 ) * t2 &
      + 1.000000004e+00_real64

    q1 = t * ((((( &
      - 0.9173e-05_real64      * t2 &
      + 0.40658e-04_real64 )   * t2 &
      - 0.99941e-04_real64 )   * t2 &
      + 0.266891e-03_real64 )  * t2 &
      - 0.1601836e-02_real64 ) * t2 &
      + 0.093749994e+00_real64 )

    ta1 = x - 0.75e+00_real64 * pi
    bj1 = a0 * ( p1 * cos ( ta1 ) - q1 * sin ( ta1 ) )
    by1 = a0 * ( p1 * sin ( ta1 ) + q1 * cos ( ta1 ) )

  end if

  dj0 = - bj1
  dj1 = bj0 - bj1 / x
  dy0 = - by1
  dy1 = by0 - by1 / x

  return
end subroutine jy01b
!> @brief subroutine jyna.
!> @return None.
!>
!> @param n [in] Argument n.
!> @param x [in] Argument x.
!> @param nm [inout] Argument nm.
!> @param bj [inout] Argument bj.
!> @param dj [inout] Argument dj.
!> @param by [inout] Argument by.
!> @param dy [inout] Argument dy.
subroutine jyna ( n, x, nm, bj, dj, by, dy )

!*****************************************************************************80
!
!! JYNA computes Bessel functions Jn(x) and Yn(x) and derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    29 April 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) N, the order.
!
!    Input, real(real64) X, the argument.
!
!    Output, integer(int32) NM, the highest order computed.
!
!    Output, real(real64) BJ(0:N), DJ(0:N), BY(0:N), DY(0:N), the values
!    of Jn(x), Jn'(x), Yn(x), Yn'(x).
!
  implicit none

  integer(int32), intent(in) :: n

  real(real64), intent(inout) :: bj(0:n)
  real(real64) bj0
  real(real64) bj1
  real(real64) bjk
  real(real64), intent(inout) :: by(0:n)
  real(real64) by0
  real(real64) by1
  real(real64) cs
  real(real64), intent(inout) :: dj(0:n)
  real(real64) dj0
  real(real64) dj1
  real(real64), intent(inout) :: dy(0:n)
  real(real64) dy0
  real(real64) dy1
  real(real64) f
  real(real64) f0
  real(real64) f1
  real(real64) f2
  integer(int32) k
  integer(int32) m
  ! integer(int32) msta1
  ! integer(int32) msta2
  integer(int32), intent(inout) :: nm
  real(real64), intent(in) :: x

  nm = n

  if ( x < 1.0e-100_real64 ) then

    do k = 0, n
      bj(k) = 0.0e+00_real64
      dj(k) = 0.0e+00_real64
      by(k) = -1.0e+300_real64
      dy(k) = 1.0e+300_real64
    end do
    bj(0) = 1.0e+00_real64
    dj(1) = 0.5e+00_real64
    return

  end if

  call jy01b ( x, bj0, dj0, bj1, dj1, by0, dy0, by1, dy1 )
  bj(0) = bj0
  bj(1) = bj1
  by(0) = by0
  by(1) = by1
  dj(0) = dj0
  dj(1) = dj1
  dy(0) = dy0
  dy(1) = dy1

  if ( n <= 1 ) then
    return
  end if

  if ( n < int ( 0.9e+00_real64 * x) ) then

    do k = 2, n
      bjk = 2.0e+00_real64 * ( k - 1.0e+00_real64 ) / x * bj1 - bj0
      bj(k) = bjk
      bj0 = bj1
      bj1 = bjk
    end do

  else

    m = msta1 ( x, 200 )

    if ( m < n ) then
      nm = m
    else
      m = msta2 ( x, n, 15 )
    end if

    f2 = 0.0e+00_real64
    f1 = 1.0e-100_real64
    do k = m, 0, -1
      f = 2.0e+00_real64 * ( k + 1.0e+00_real64 ) / x * f1 - f2
      if ( k <= nm ) then
        bj(k) = f
      end if
      f2 = f1
      f1 = f
    end do

    if ( abs ( bj1 ) < abs ( bj0 ) ) then
      cs = bj0 / f
    else
      cs = bj1 / f2
    end if

    do k = 0, nm
      bj(k) = cs * bj(k)
    end do

  end if

  do k = 2, nm
    dj(k) = bj(k-1) - k / x * bj(k)
  end do

  f0 = by(0)
  f1 = by(1)
  do k = 2, nm
    f = 2.0e+00_real64 * ( k - 1.0e+00_real64 ) / x * f1 - f0
    by(k) = f
    f0 = f1
    f1 = f
  end do

  do k = 2, nm
    dy(k) = by(k-1) - k * by(k) / x
  end do

  return
end subroutine jyna
!> @brief subroutine jynb.
!> @return None.
!>
!> @param n [in] Argument n.
!> @param x [in] Argument x.
!> @param nm [inout] Argument nm.
!> @param bj [inout] Argument bj.
!> @param dj [inout] Argument dj.
!> @param by [inout] Argument by.
!> @param dy [inout] Argument dy.
subroutine jynb ( n, x, nm, bj, dj, by, dy )

!*****************************************************************************80
!
!! JYNB computes Bessel functions Jn(x) and Yn(x) and derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    02 August 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) N, the order.
!
!    Input, real(real64) X, the argument.
!
!    Output, integer(int32) NM, the highest order computed.
!
!    Output, real(real64) BJ(0:N), DJ(0:N), BY(0:N), DY(0:N), the values
!    of Jn(x), Jn'(x), Yn(x), Yn'(x).
!
  implicit none

  integer(int32), intent(in) :: n

  real(real64), parameter, dimension ( 4 ) :: a = [&
    -0.7031250000000000e-01_real64, 0.1121520996093750e+00_real64, &
    -0.5725014209747314e+00_real64, 0.6074042001273483e+01_real64]
  real(real64), parameter, dimension ( 4 ) :: a1 = [&
    0.1171875000000000e+00_real64, -0.1441955566406250e+00_real64, &
    0.6765925884246826e+00_real64, -0.6883914268109947e+01_real64]
  real(real64), parameter, dimension ( 4 ) :: b = [&
    0.7324218750000000e-01_real64, -0.2271080017089844e+00_real64, &
    0.1727727502584457e+01_real64, -0.2438052969955606e+02_real64]
  real(real64), parameter, dimension ( 4 ) :: b1 = [&
    -0.1025390625000000e+00_real64, 0.2775764465332031e+00_real64, &
    -0.1993531733751297e+01_real64, 0.2724882731126854e+02_real64]
  real(real64), intent(inout) :: bj(0:n)
  real(real64) bj0
  real(real64) bj1
  real(real64) bjk
  real(real64) bs
  real(real64), intent(inout) :: by(0:n)
  real(real64) by0
  real(real64) by1
  real(real64) byk
  real(real64) cu
  real(real64), intent(inout) :: dj(0:n)
  real(real64), intent(inout) :: dy(0:n)
  real(real64) ec
  real(real64) f
  real(real64) f1
  real(real64) f2
  integer(int32) k
  integer(int32) m
  ! integer(int32) msta1
  ! integer(int32) msta2
  integer(int32), intent(inout) :: nm
  real(real64) p0
  real(real64) p1
  real(real64) pi
  real(real64) q0
  real(real64) q1
  real(real64) r2p
  real(real64) s0
  real(real64) su
  real(real64) sv
  real(real64) t1
  real(real64) t2
  real(real64), intent(in) :: x

  pi = 3.141592653589793e+00_real64
  r2p = 0.63661977236758e+00_real64
  nm = n

  if ( x < 1.0e-100_real64 ) then
    do k = 0, n
      bj(k) = 0.0e+00_real64
      dj(k) = 0.0e+00_real64
      by(k) = -1.0e+300_real64
      dy(k) = 1.0e+300_real64
    end do
    bj(0) = 1.0e+00_real64
    dj(1) = 0.5e+00_real64
    return
  end if

  if ( x <= 300.0e+00_real64 .or. int ( 0.9e+00_real64 * x ) < n ) then

    if ( n == 0 ) then
      nm = 1
    end if

    m = msta1 ( x, 200 )

    if ( m < nm ) then
      nm = m
    else
      m = msta2 ( x, nm, 15 )
    end if

    bs = 0.0e+00_real64
    su = 0.0e+00_real64
    sv = 0.0e+00_real64
    f2 = 0.0e+00_real64
    f1 = 1.0e-100_real64

    do k = m, 0, -1
      f = 2.0e+00_real64 * ( k + 1.0e+00_real64 ) / x * f1 - f2
      if ( k <= nm ) then
        bj(k) = f
      end if
      if ( k == 2 * int ( k / 2 ) .and. k /= 0 ) then
        bs = bs + 2.0e+00_real64 * f
        su = su + ( -1.0e+00_real64 ) ** ( k / 2 ) * f / k
      else if ( 1 < k ) then
        sv = sv + ( -1.0e+00_real64 ) ** ( k / 2 ) * k / ( k * k - 1.0e+00_real64 ) * f
      end if
      f2 = f1
      f1 = f
    end do

    s0 = bs + f
    do k = 0, nm
      bj(k) = bj(k) / s0
    end do

    ec = log ( x / 2.0e+00_real64 ) + 0.5772156649015329e+00_real64
    by0 = r2p * ( ec * bj(0) - 4.0e+00_real64 * su / s0 )
    by(0) = by0
    by1 = r2p * ( ( ec - 1.0e+00_real64 ) * bj(1) - bj(0) / x - 4.0e+00_real64 * sv / s0 )
    by(1) = by1

  else

    t1 = x - 0.25e+00_real64 * pi
    p0 = 1.0e+00_real64
    q0 = -0.125e+00_real64 / x
    do k = 1, 4
      p0 = p0 + a(k) * x ** ( - 2 * k )
      q0 = q0 + b(k) * x ** ( - 2 * k - 1 )
    end do
    cu = sqrt ( r2p / x )
    bj0 = cu * ( p0 * cos ( t1 ) - q0 * sin ( t1 ) )
    by0 = cu * ( p0 * sin ( t1 ) + q0 * cos ( t1 ) )
    bj(0) = bj0
    by(0) = by0
    t2 = x - 0.75e+00_real64 * pi
    p1 = 1.0e+00_real64
    q1 = 0.375e+00_real64 / x
    do k = 1, 4
      p1 = p1 + a1(k) * x ** ( - 2 * k )
      q1 = q1 + b1(k) * x ** ( - 2 * k - 1 )
    end do
    bj1 = cu * ( p1 * cos ( t2 ) - q1 * sin ( t2 ) )
    by1 = cu * ( p1 * sin ( t2 ) + q1 * cos ( t2 ) )
    bj(1) = bj1
    by(1) = by1
    do k = 2, nm
      bjk = 2.0e+00_real64 * ( k - 1.0e+00_real64 ) / x * bj1 - bj0
      bj(k) = bjk
      bj0 = bj1
      bj1 = bjk
    end do
  end if

  dj(0) = -bj(1)
  do k = 1, nm
    dj(k) = bj(k-1) - k / x * bj(k)
  end do

  do k = 2, nm
    byk = 2.0e+00_real64 * ( k - 1.0e+00_real64 ) * by1 / x - by0
    by(k) = byk
    by0 = by1
    by1 = byk
  end do

  dy(0) = -by(1)
  do k = 1, nm
    dy(k) = by(k-1) - k * by(k) / x
  end do

  return
end subroutine jynb
!> @brief subroutine jyndd.
!> @return None.
!>
!> @param n [in] Argument n.
!> @param x [in] Argument x.
!> @param bjn [inout] Argument bjn.
!> @param djn [inout] Argument djn.
!> @param fjn [inout] Argument fjn.
!> @param byn [inout] Argument byn.
!> @param dyn [inout] Argument dyn.
!> @param fyn [inout] Argument fyn.
pure subroutine jyndd ( n, x, bjn, djn, fjn, byn, dyn, fyn )

!*****************************************************************************80
!
!! JYNDD: Bessel functions Jn(x) and Yn(x), first and second derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    02 August 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) N, the order.
!
!    Input, real(real64) X, the argument.
!
!    Output, real(real64) BJN, DJN, FJN, BYN, DYN, FYN, the values of
!    Jn(x), Jn'(x), Jn"(x), Yn(x), Yn'(x), Yn"(x).
!
  implicit none

  real(real64) bj(102)
  real(real64), intent(inout) :: bjn
  real(real64), intent(inout) :: byn
  real(real64) bs
  real(real64) by(102)
  real(real64), intent(inout) :: djn
  real(real64), intent(inout) :: dyn
  real(real64) e0
  real(real64) ec
  real(real64) f
  real(real64) f0
  real(real64) f1
  real(real64), intent(inout) :: fjn
  real(real64), intent(inout) :: fyn
  integer(int32) k
  integer(int32) m
  integer(int32) mt
  integer(int32), intent(in) :: n
  integer(int32) nt
  real(real64) s1
  real(real64) su
  real(real64), intent(in) :: x

  do nt = 1, 900
    mt = int ( 0.5e+00_real64 * log10 ( 6.28e+00_real64 * nt ) &
      - nt * log10 ( 1.36e+00_real64 * abs ( x ) / nt ) )
    if ( 20 < mt ) then
      exit
    end if
  end do

  m = nt
  bs = 0.0e+00_real64
  f0 = 0.0e+00_real64
  f1 = 1.0e-35_real64
  su = 0.0e+00_real64
  do k = m, 0, -1
    f = 2.0e+00_real64 * ( k + 1.0e+00_real64 ) * f1 / x - f0
    if ( k <= n + 1 ) then
      bj(k+1) = f
    end if
    if ( k == 2 * int ( k / 2 ) ) then
      bs = bs + 2.0e+00_real64 * f
      if ( k /= 0 ) then
        su = su + ( -1.0e+00_real64 ) ** ( k / 2 ) * f / k
      end if
    end if
    f0 = f1
    f1 = f
  end do

  do k = 0, n + 1
    bj(k+1) = bj(k+1) / ( bs - f )
  end do

  bjn = bj(n+1)
  ec = 0.5772156649015329e+00_real64
  e0 = 0.3183098861837907e+00_real64
  s1 = 2.0e+00_real64 * e0 * ( log ( x / 2.0e+00_real64 ) + ec ) * bj(1)
  f0 = s1 - 8.0e+00_real64 * e0 * su / ( bs - f )
  f1 = ( bj(2) * f0 - 2.0e+00_real64 * e0 / x ) / bj(1)

  by(1) = f0
  by(2) = f1
  do k = 2, n + 1 
    f = 2.0e+00_real64 * ( k - 1.0e+00_real64 ) * f1 / x - f0
    by(k+1) = f
    f0 = f1
    f1 = f
  end do

  byn = by(n+1)
  djn = - bj(n+2) + n * bj(n+1) / x
  dyn = - by(n+2) + n * by(n+1) / x
  fjn = ( n * n / ( x * x ) - 1.0e+00_real64 ) * bjn - djn / x
  fyn = ( n * n / ( x * x ) - 1.0e+00_real64 ) * byn - dyn / x

  return
end subroutine jyndd
!> @brief subroutine jyv.
!> @return None.
!>
!> @param v [in] Argument v.
!> @param x [in] Argument x.
!> @param vm [inout] Argument vm.
!> @param bj [inout] Argument bj.
!> @param dj [inout] Argument dj.
!> @param by [inout] Argument by.
!> @param dy [inout] Argument dy.
subroutine jyv ( v, x, vm, bj, dj, by, dy )

!*****************************************************************************80
!
!! JYV computes Bessel functions Jv(x) and Yv(x) and their derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    02 August 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) V, the order of Jv(x) and Yv(x).
!
!    Input, real(real64) X, the argument of Jv(x) and Yv(x).
!
!    Output, real(real64) VM, the highest order computed.
!
!    Output, real(real64) BJ(0:N), DJ(0:N), BY(0:N), DY(0:N),
!    the values of Jn+v0(x), Jn+v0'(x), Yn+v0(x), Yn+v0'(x).
!
  implicit none

  real(real64) a
  real(real64) a0
  real(real64) b
  real(real64), intent(inout) :: bj(0:)
  real(real64) bju0
  real(real64) bju1
  real(real64) bjv0
  real(real64) bjv1
  real(real64) bjvl
  real(real64), intent(inout) :: by(0:)
  real(real64) byv0
  real(real64) byv1
  real(real64) byvk
  real(real64) ck
  real(real64) cs
  real(real64) cs0
  real(real64) cs1
  real(real64), intent(inout) :: dj(0:)
  real(real64), intent(inout) :: dy(0:)
  real(real64) ec
  real(real64) el
  real(real64) f
  real(real64) f0
  real(real64) f1
  real(real64) f2
  real(real64) ga
  real(real64) gb
  integer(int32) j
  integer(int32) k
  integer(int32) k0
  integer(int32) l
  integer(int32) m
  ! integer(int32) msta1
  ! integer(int32) msta2
  integer(int32) n
  real(real64) pi
  real(real64) pv0
  real(real64) pv1
  real(real64) px
  real(real64) qx
  real(real64) r
  real(real64) r0
  real(real64) r1
  real(real64) rp
  real(real64) rp2
  real(real64) rq
  real(real64) sk
  real(real64), intent(in) :: v
  real(real64) v0
  real(real64) vg
  real(real64) vl
  real(real64), intent(inout) :: vm
  real(real64) vv
  real(real64) w0
  real(real64) w1
  real(real64), intent(in) :: x
  real(real64) x2
  real(real64) xk

  el = 0.5772156649015329e+00_real64
  pi = 3.141592653589793e+00_real64
  rp2 = 0.63661977236758e+00_real64
  x2 = x * x
  n = int ( v )
  v0 = v - n

  if ( x < 1.0e-100_real64 ) then

    do k = 0, n
      bj(k) = 0.0e+00_real64
      dj(k) = 0.0e+00_real64
      by(k) = -1.0e+300_real64
      dy(k) = 1.0e+300_real64
    end do

    if ( v0 == 0.0e+00_real64 ) then
      bj(0) = 1.0e+00_real64
      dj(1) = 0.5e+00_real64
    else
      dj(0) = 1.0e+300_real64
    end if
    vm = v  
    return

  end if

  if ( x <= 12.0e+00_real64 ) then

    do l = 0, 1
      vl = v0 + l
      bjvl = 1.0e+00_real64
      r = 1.0e+00_real64
      do k = 1, 40
        r = -0.25e+00_real64 * r * x2 / ( k * ( k + vl ) )
        bjvl = bjvl + r
        if ( abs ( r ) < abs ( bjvl ) * 1.0e-15_real64 ) then
          exit
        end if
      end do

      vg = 1.0e+00_real64 + vl
      call gamma ( vg, ga )
      a = ( 0.5e+00_real64 * x ) ** vl / ga

      if ( l == 0 ) then
        bjv0 = bjvl * a
      else
        bjv1 = bjvl * a
      end if

    end do

  else

    if ( x < 35.0e+00_real64 ) then
      k0 = 11
    else if ( x < 50.0e+00_real64 ) then
      k0 = 10
    else
      k0 = 8
    end if

    do j = 0, 1

      vv = 4.0e+00_real64 * ( j + v0 ) * ( j + v0 )
      px = 1.0e+00_real64
      rp = 1.0e+00_real64
      do k = 1, k0
        rp = -0.78125e-02_real64 * rp &
          * ( vv - ( 4.0e+00_real64 * k - 3.0e+00_real64 ) ** 2 ) &
          * ( vv - ( 4.0e+00_real64 * k - 1.0e+00_real64 ) ** 2 ) &
          / ( k * ( 2.0e+00_real64 * k - 1.0e+00_real64 ) * x2 )
        px = px + rp
      end do
      qx = 1.0e+00_real64
      rq = 1.0e+00_real64
      do k = 1, k0
        rq = -0.78125e-02_real64 * rq &
          * ( vv - ( 4.0e+00_real64 * k - 1.0e+00_real64 ) ** 2 ) &
          * ( vv - ( 4.0e+00_real64 * k + 1.0e+00_real64 ) ** 2 ) &
          / ( k * ( 2.0e+00_real64 * k + 1.0e+00_real64 ) * x2 )
        qx = qx + rq
      end do
      qx = 0.125e+00_real64 * ( vv - 1.0e+00_real64 ) * qx / x
      xk = x - ( 0.5e+00_real64 * ( j + v0 ) + 0.25e+00_real64 ) * pi
      a0 = sqrt ( rp2 / x )
      ck = cos ( xk )
      sk = sin ( xk )
      if ( j == 0 ) then
        bjv0 = a0 * ( px * ck - qx * sk )
        byv0 = a0 * ( px * sk + qx * ck )
      else if ( j == 1 ) then
        bjv1 = a0 * ( px * ck - qx * sk )
        byv1 = a0 * ( px * sk + qx * ck )
      end if

    end do

  end if

  bj(0) = bjv0
  bj(1) = bjv1
  dj(0) = v0 / x * bj(0) - bj(1)
  dj(1) = - ( 1.0e+00_real64 + v0 ) / x * bj(1) + bj(0)

  if ( 2 <= n .and. n <= int ( 0.9e+00_real64 * x ) ) then
    f0 = bjv0
    f1 = bjv1
    do k = 2, n
      f = 2.0e+00_real64 * ( k + v0 - 1.0e+00_real64 ) / x * f1 - f0
      bj(k) = f
      f0 = f1
      f1 = f
    end do
  else if ( 2 <= n ) then
    m = msta1 ( x, 200 )
    if ( m < n ) then
      n = m
    else
      m = msta2 ( x, n, 15 )
    end if
    f2 = 0.0e+00_real64
    f1 = 1.0e-100_real64
    do k = m, 0, -1
      f = 2.0e+00_real64 * ( v0 + k + 1.0e+00_real64 ) / x * f1 - f2
      if ( k <= n ) then
        bj(k) = f
      end if
      f2 = f1
      f1 = f
    end do

    if ( abs ( bjv1 ) < abs ( bjv0 ) ) then
      cs = bjv0 / f
    else
      cs = bjv1 / f2
    end if
    do k = 0, n
      bj(k) = cs * bj(k)
    end do
  end if

  do k = 2, n
    dj(k) = - ( k + v0 ) / x * bj(k) + bj(k-1)
  end do

  if ( x <= 12.0e+00_real64 ) then

    if ( v0 /= 0.0e+00_real64 ) then

      do l = 0, 1

        vl = v0 + l
        bjvl = 1.0e+00_real64
        r = 1.0e+00_real64
        do k = 1, 40
          r = -0.25e+00_real64 * r * x2 / ( k * ( k - vl ) )
          bjvl = bjvl + r
          if ( abs ( r ) < abs ( bjvl ) * 1.0e-15_real64 ) then
            exit
          end if
        end do

        vg = 1.0e+00_real64 - vl
        call gamma ( vg, gb )
        b = ( 2.0e+00_real64 / x ) ** vl / gb

        if ( l == 0 ) then
          bju0 = bjvl * b
        else
          bju1 = bjvl * b
        end if

      end do

      pv0 = pi * v0
      pv1 = pi * ( 1.0e+00_real64 + v0 )
      byv0 = ( bjv0 * cos ( pv0 ) - bju0 ) / sin ( pv0 )
      byv1 = ( bjv1 * cos ( pv1 ) - bju1 ) / sin ( pv1 )

    else

      ec = log ( x / 2.0e+00_real64 ) + el
      cs0 = 0.0e+00_real64
      w0 = 0.0e+00_real64
      r0 = 1.0e+00_real64
      do k = 1, 30
        w0 = w0 + 1.0e+00_real64 / k
        r0 = -0.25e+00_real64 * r0 / ( k * k ) * x2
        cs0 = cs0 + r0 * w0
      end do
      byv0 = rp2 * ( ec * bjv0 - cs0 )
      cs1 = 1.0e+00_real64
      w1 = 0.0e+00_real64
      r1 = 1.0e+00_real64
      do k = 1, 30
        w1 = w1 + 1.0e+00_real64 / k
        r1 = -0.25e+00_real64 * r1 / ( k * ( k + 1 ) ) * x2
        cs1 = cs1 + r1 * ( 2.0e+00_real64 * w1 + 1.0e+00_real64 / ( k + 1.0e+00_real64 ) )
      end do
      byv1 = rp2 * ( ec * bjv1 - 1.0e+00_real64 / x - 0.25e+00_real64 * x * cs1 )

    end if

  end if

  by(0) = byv0
  by(1) = byv1
  do k = 2, n
    byvk = 2.0e+00_real64 * ( v0 + k - 1.0e+00_real64 ) / x * byv1 - byv0
    by(k) = byvk
    byv0 = byv1
    byv1 = byvk
  end do

  dy(0) = v0 / x * by(0) - by(1)
  do k = 1, n
    dy(k) = - ( k + v0 ) / x * by(k) + by(k-1)
  end do

  vm = n + v0

  return
end subroutine jyv
!> @brief subroutine jyzo.
!> @return None.
!>
!> @param n [in] Argument n.
!> @param nt [in] Argument nt.
!> @param rj0 [inout] Argument rj0.
!> @param rj1 [inout] Argument rj1.
!> @param ry0 [inout] Argument ry0.
!> @param ry1 [inout] Argument ry1.
subroutine jyzo ( n, nt, rj0, rj1, ry0, ry1 )

!*****************************************************************************80
!
!! JYZO computes the zeros of Bessel functions Jn(x), Yn(x) and derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    28 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) N, the order of the Bessel functions.
!
!    Input, integer(int32) NT, the number of zeros.
!
!    Output, real(real64) RJ0(NT), RJ1(NT), RY0(NT), RY1(NT), the zeros 
!    of Jn(x), Jn'(x), Yn(x), Yn'(x).
!
  implicit none

  integer(int32), intent(in) :: nt

  real(real64) bjn
  real(real64) byn
  real(real64) djn
  real(real64) dyn
  real(real64) fjn
  real(real64) fyn
  integer(int32) l
  integer(int32), intent(in) :: n
  real(real64) n_r8
  real(real64), intent(inout) :: rj0(nt)
  real(real64), intent(inout) :: rj1(nt)
  real(real64), intent(inout) :: ry0(nt)
  real(real64), intent(inout) :: ry1(nt)
  real(real64) x
  real(real64) x0

  n_r8 = real ( n, kind = real64 )

  if ( n <= 20 ) then
    x = 2.82141e+00_real64 + 1.15859e+00_real64 * n_r8 
  else
    x = n + 1.85576e+00_real64 * n_r8 ** 0.33333e+00_real64 &
      + 1.03315e+00_real64 / n_r8 ** 0.33333e+00_real64
  end if

  l = 0

  do

    x0 = x
    call jyndd ( n, x, bjn, djn, fjn, byn, dyn, fyn )
    x = x - bjn / djn

    if ( 1.0e-09_real64 < abs ( x - x0 ) ) then
      cycle
    end if

    l = l + 1
    rj0(l) = x
    x = x + 3.1416e+00_real64 + ( 0.0972e+00_real64 + 0.0679e+00_real64 * n_r8 &
      - 0.000354e+00_real64 * n_r8 ** 2 ) / l

    if ( nt <= l ) then
      exit
    end if

  end do

  if ( n <= 20 ) then
    x = 0.961587e+00_real64 + 1.07703e+00_real64 * n_r8 
  else
    x = n_r8 + 0.80861e+00_real64 * n_r8 ** 0.33333e+00_real64 &
      + 0.07249e+00_real64 / n_r8 ** 0.33333e+00_real64
  end if

  if ( n == 0 ) then
    x = 3.8317e+00_real64
  end if

  l = 0

  do

    x0 = x
    call jyndd ( n, x, bjn, djn, fjn, byn, dyn, fyn )
    x = x - djn / fjn
    if ( 1.0e-09_real64 < abs ( x - x0 ) ) then
      cycle
    end if
    l = l + 1
    rj1(l) = x
    x = x + 3.1416e+00_real64 + ( 0.4955e+00_real64 + 0.0915e+00_real64 * n_r8 &
      - 0.000435e+00_real64 * n_r8 ** 2 ) / l

    if ( nt <= l ) then
      exit
    end if

  end do

  if ( n <= 20 ) then
    x = 1.19477e+00_real64 + 1.08933e+00_real64 * n_r8 
  else
    x = n_r8 + 0.93158e+00_real64 * n_r8 ** 0.33333e+00_real64 &
      + 0.26035e+00_real64 / n_r8 ** 0.33333e+00_real64
  end if
 
  l = 0

  do

    x0 = x
    call jyndd ( n, x, bjn, djn, fjn, byn, dyn, fyn )
    x = x - byn / dyn

    if ( 1.0e-09_real64 < abs ( x - x0 ) ) then
      cycle
    end if

    l = l + 1
    ry0(l) = x 
    x = x + 3.1416e+00_real64 + ( 0.312e+00_real64 + 0.0852e+00_real64 * n_r8 &
      - 0.000403e+00_real64 * n_r8 ** 2 ) / l

    if ( nt <= l ) then
      exit
    end if

  end do

  if ( n <= 20 ) then
    x = 2.67257e+00_real64 + 1.16099e+00_real64 * n_r8 
  else
    x = n_r8 + 1.8211e+00_real64 * n_r8 ** 0.33333e+00_real64 &
      + 0.94001e+00_real64 / n_r8 ** 0.33333e+00_real64
  end if
  
  l = 0

  do

    x0 = x
    call jyndd ( n, x, bjn, djn, fjn, byn, dyn, fyn )
    x = x - dyn / fyn

    if ( 1.0e-09_real64 < abs ( x - x0 ) ) then
      cycle
    end if

    l = l + 1
    ry1(l) = x
    x = x + 3.1416e+00_real64 + ( 0.197e+00_real64 + 0.0643e+00_real64 * n_r8 &
      -0.000286e+00_real64 * n_r8 ** 2 ) / l 

    if ( nt <= l ) then
      exit
    end if

  end do

  return
end subroutine jyzo
!> @brief subroutine klvna.
!> @return None.
!>
!> @param x [in] Argument x.
!> @param ber [inout] Argument ber.
!> @param bei [inout] Argument bei.
!> @param ger [inout] Argument ger.
!> @param gei [inout] Argument gei.
!> @param der [inout] Argument der.
!> @param dei [inout] Argument dei.
!> @param her [inout] Argument her.
!> @param hei [inout] Argument hei.
pure subroutine klvna ( x, ber, bei, ger, gei, der, dei, her, hei )

!*****************************************************************************80
!
!! KLVNA: Kelvin functions ber(x), bei(x), ker(x), and kei(x), and derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    03 August 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) X, the argument.
!
!    Output, real(real64) BER, BEI, GER, GEI, DER, DEI, HER, HEI, 
!    the values of ber x, bei x, ker x, kei x, ber'x, bei'x, ker'x, kei'x.
!
  implicit none

  real(real64), intent(inout) :: bei
  real(real64), intent(inout) :: ber
  real(real64) cn0
  real(real64) cp0
  real(real64) cs
  real(real64), intent(inout) :: dei
  real(real64), intent(inout) :: der
  real(real64) el
  real(real64) eps
  real(real64) fac
  real(real64), intent(inout) :: gei
  real(real64), intent(inout) :: ger
  real(real64) gs
  real(real64), intent(inout) :: hei
  real(real64), intent(inout) :: her
  integer(int32) k
  integer(int32) km
  integer(int32) m
  real(real64) pi
  real(real64) pn0
  real(real64) pn1
  real(real64) pp0
  real(real64) pp1
  real(real64) qn0
  real(real64) qn1
  real(real64) qp0
  real(real64) qp1
  real(real64) r
  real(real64) r0
  real(real64) r1
  real(real64) rc
  real(real64) rs
  real(real64) sn0
  real(real64) sp0
  real(real64) ss
  real(real64), intent(in) :: x
  real(real64) x2
  real(real64) x4
  real(real64) xc1
  real(real64) xc2
  real(real64) xd
  real(real64) xe1
  real(real64) xe2
  real(real64) xt

  pi = 3.141592653589793e+00_real64
  el = 0.5772156649015329e+00_real64
  eps = 1.0e-15_real64

  if ( x == 0.0e+00_real64 ) then
    ber = 1.0e+00_real64
    bei = 0.0e+00_real64
    ger = 1.0e+300_real64
    gei = -0.25e+00_real64 * pi
    der = 0.0e+00_real64
    dei = 0.0e+00_real64
    her = -1.0e+300_real64
    hei = 0.0e+00_real64
    return
  end if

  x2 = 0.25e+00_real64 * x * x
  x4 = x2 * x2

  if ( abs ( x ) < 10.0e+00_real64 ) then

    ber = 1.0e+00_real64
    r = 1.0e+00_real64
    do m = 1, 60
      r = -0.25e+00_real64 * r / ( m * m ) / ( 2.0e+00_real64 * m - 1.0e+00_real64 ) ** 2 * x4
      ber = ber + r
      if ( abs ( r ) < abs ( ber ) * eps ) then
        exit
      end if
    end do

    bei = x2
    r = x2
    do m = 1, 60
      r = -0.25e+00_real64 * r / ( m * m ) / ( 2.0e+00_real64 * m + 1.0e+00_real64 ) ** 2 * x4
      bei = bei + r
      if ( abs ( r ) < abs ( bei ) * eps ) then
        exit
      end if
    end do

    ger = - ( log ( x / 2.0e+00_real64 ) + el ) * ber + 0.25e+00_real64 * pi * bei
    r = 1.0e+00_real64
    gs = 0.0e+00_real64
    do m = 1, 60
      r = -0.25e+00_real64 * r / ( m * m ) / ( 2.0e+00_real64 * m - 1.0e+00_real64 ) ** 2 * x4
      gs = gs + 1.0e+00_real64 / ( 2.0e+00_real64 * m - 1.0e+00_real64 ) + 1.0e+00_real64 / ( 2.0e+00_real64 * m )
      ger = ger + r * gs
      if ( abs ( r * gs ) < abs ( ger ) * eps ) then
        exit
      end if
    end do

    gei = x2 - ( log ( x / 2.0e+00_real64 ) + el ) * bei - 0.25e+00_real64 * pi * ber
    r = x2
    gs = 1.0e+00_real64
    do m = 1, 60
      r = -0.25e+00_real64 * r / ( m * m ) / ( 2.0e+00_real64 * m + 1.0e+00_real64 ) ** 2 * x4
      gs = gs + 1.0e+00_real64 / ( 2.0e+00_real64 * m ) + 1.0e+00_real64 / ( 2.0e+00_real64 * m + 1.0e+00_real64 )
      gei = gei + r * gs
      if ( abs ( r * gs ) < abs ( gei ) * eps ) then
        exit
      end if
    end do

    der = -0.25e+00_real64 * x * x2
    r = der
    do m = 1, 60
      r = -0.25e+00_real64 * r / m / ( m + 1.0e+00_real64 ) &
        / ( 2.0e+00_real64 * m + 1.0e+00_real64 ) ** 2 * x4
      der = der + r
      if ( abs ( r ) < abs ( der ) * eps ) then
        exit
      end if
    end do

    dei = 0.5e+00_real64 * x
    r = dei
    do m = 1, 60
      r = -0.25e+00_real64 * r / ( m * m ) / ( 2.0e+00_real64 * m - 1.0e+00_real64 ) &
        / ( 2.0e+00_real64 * m + 1.0e+00_real64 ) * x4
      dei = dei + r
      if ( abs ( r ) < abs ( dei ) * eps ) then
        exit
      end if
    end do

    r = -0.25e+00_real64 * x * x2
    gs = 1.5e+00_real64
    her = 1.5e+00_real64 * r - ber / x &
      - ( log ( x / 2.0e+00_real64 ) + el ) * der + 0.25e+00_real64 * pi * dei
    do m = 1, 60
      r = -0.25e+00_real64 * r / m / ( m + 1.0e+00_real64 ) &
        / ( 2.0e+00_real64 * m + 1.0e+00_real64 ) ** 2 * x4
      gs = gs + 1.0e+00_real64 / ( 2 * m + 1.0e+00_real64 ) + 1.0e+00_real64 &
        / ( 2 * m + 2.0e+00_real64 )
      her = her + r * gs
      if ( abs ( r * gs ) < abs ( her ) * eps ) then
        exit
      end if
    end do

    r = 0.5e+00_real64 * x
    gs = 1.0e+00_real64
    hei = 0.5e+00_real64 * x - bei / x &
      - ( log ( x / 2.0e+00_real64 ) + el ) * dei - 0.25e+00_real64 * pi * der
    do m = 1, 60
      r = -0.25e+00_real64 * r / ( m * m ) / ( 2 * m - 1.0e+00_real64 ) &
        / ( 2 * m + 1.0e+00_real64 ) * x4
      gs = gs + 1.0e+00_real64 / ( 2.0e+00_real64 * m ) + 1.0e+00_real64 &
        / ( 2 * m + 1.0e+00_real64 )
      hei = hei + r * gs
      if ( abs ( r * gs ) < abs ( hei ) * eps ) then 
        return
      end if
    end do

  else

    pp0 = 1.0e+00_real64
    pn0 = 1.0e+00_real64
    qp0 = 0.0e+00_real64
    qn0 = 0.0e+00_real64
    r0 = 1.0e+00_real64

    if ( abs ( x ) < 40.0e+00_real64 ) then
      km = 18
    else
      km = 10
    end if

    fac = 1.0e+00_real64
    do k = 1, km
      fac = -fac
      xt = 0.25e+00_real64 * k * pi - int ( 0.125e+00_real64 * k ) * 2.0e+00_real64 * pi
      cs = cos ( xt )
      ss = sin ( xt )
      r0 = 0.125e+00_real64 * r0 * ( 2.0e+00_real64 * k - 1.0e+00_real64 ) ** 2 / k / x
      rc = r0 * cs
      rs = r0 * ss
      pp0 = pp0 + rc
      pn0 = pn0 + fac * rc
      qp0 = qp0 + rs
      qn0 = qn0 + fac * rs
    end do

    xd = x / sqrt (2.0e+00_real64 )
    xe1 = exp ( xd )
    xe2 = exp ( - xd )
    xc1 = 1.0e+00_real64 / sqrt ( 2.0e+00_real64 * pi * x )
    xc2 = sqrt ( 0.5e+00_real64 * pi / x )
    cp0 = cos ( xd + 0.125e+00_real64 * pi )
    cn0 = cos ( xd - 0.125e+00_real64 * pi )
    sp0 = sin ( xd + 0.125e+00_real64 * pi )
    sn0 = sin ( xd - 0.125e+00_real64 * pi )
    ger = xc2 * xe2 * (  pn0 * cp0 - qn0 * sp0 )
    gei = xc2 * xe2 * ( -pn0 * sp0 - qn0 * cp0 )
    ber = xc1 * xe1 * (  pp0 * cn0 + qp0 * sn0 ) - gei / pi
    bei = xc1 * xe1 * (  pp0 * sn0 - qp0 * cn0 ) + ger / pi
    pp1 = 1.0e+00_real64
    pn1 = 1.0e+00_real64
    qp1 = 0.0e+00_real64
    qn1 = 0.0e+00_real64
    r1 = 1.0e+00_real64
    fac = 1.0e+00_real64

    do k = 1, km
      fac = -fac
      xt = 0.25e+00_real64 * k * pi - int ( 0.125e+00_real64 * k ) * 2.0e+00_real64 * pi
      cs = cos ( xt )
      ss = sin ( xt )
      r1 = 0.125e+00_real64 * r1 &
        * ( 4.0e+00_real64 - ( 2.0e+00_real64 * k - 1.0e+00_real64 ) ** 2 ) / k / x
      rc = r1 * cs
      rs = r1 * ss
      pp1 = pp1 + fac * rc
      pn1 = pn1 + rc
      qp1 = qp1 + fac * rs
      qn1 = qn1 + rs
    end do

    her = xc2 * xe2 * ( - pn1 * cn0 + qn1 * sn0 )
    hei = xc2 * xe2 * (   pn1 * sn0 + qn1 * cn0 )
    der = xc1 * xe1 * (   pp1 * cp0 + qp1 * sp0 ) - hei / pi
    dei = xc1 * xe1 * (   pp1 * sp0 - qp1 * cp0 ) + her / pi

  end if

  return
end subroutine klvna
!> @brief subroutine klvnb.
!> @return None.
!>
!> @param x [in] Argument x.
!> @param ber [inout] Argument ber.
!> @param bei [inout] Argument bei.
!> @param ger [inout] Argument ger.
!> @param gei [inout] Argument gei.
!> @param der [inout] Argument der.
!> @param dei [inout] Argument dei.
!> @param her [inout] Argument her.
!> @param hei [inout] Argument hei.
pure subroutine klvnb ( x, ber, bei, ger, gei, der, dei, her, hei )

!*****************************************************************************80
!
!! KLVNB: Kelvin functions ber(x), bei(x), ker(x), and kei(x), and derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    03 August 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) X, the argument.
!
!    Output, real(real64) BER, BEI, GER, GEI, DER, DEI, HER, HEI, 
!    the values of ber x, bei x, ker x, kei x, ber'x, bei'x, ker'x, kei'x.
!
  implicit none

  real(real64), intent(inout) :: bei
  real(real64), intent(inout) :: ber
  real(real64) csn
  real(real64) csp
  real(real64), intent(inout) :: dei
  real(real64), intent(inout) :: der
  real(real64) fxi
  real(real64) fxr
  real(real64), intent(inout) :: gei
  real(real64), intent(inout) :: ger
  real(real64), intent(inout) :: hei
  real(real64), intent(inout) :: her
  integer(int32) l
  real(real64) pi
  real(real64) pni
  real(real64) pnr
  real(real64) ppi
  real(real64) ppr
  real(real64) ssn
  real(real64) ssp
  real(real64) t
  real(real64) t2
  real(real64) tni
  real(real64) tnr
  real(real64) tpi
  real(real64) tpr
  real(real64) u
  real(real64) v
  real(real64), intent(in) :: x
  real(real64) yc1
  real(real64) yc2
  real(real64) ye1
  real(real64) ye2
  real(real64) yd

  pi = 3.141592653589793e+00_real64

  if ( x == 0.0e+00_real64 ) then

    ber = 1.0e+00_real64
    bei = 0.0e+00_real64
    ger = 1.0e+300_real64
    gei = -0.25e+00_real64 * pi
    der = 0.0e+00_real64
    dei = 0.0e+00_real64
    her = -1.0e+300_real64
    hei = 0.0e+00_real64

  else if ( x < 8.0e+00_real64 ) then

    t = x / 8.0e+00_real64
    t2 = t * t
    u = t2 * t2

    ber = (((((( &
      - 0.901e-05_real64 * u &
      + 0.122552e-02_real64 ) * u &
      - 0.08349609e+00_real64 ) * u &
      + 2.64191397e+00_real64 ) * u &
      - 32.36345652e+00_real64 ) * u &
      + 113.77777774e+00_real64 ) * u &
      - 64.0e+00_real64 ) * u &
      + 1.0e+00_real64

    bei = t * t * (((((( &
        0.11346e-03_real64 * u &
      - 0.01103667e+00_real64 ) * u &
      + 0.52185615e+00_real64 ) * u &
      - 10.56765779e+00_real64 ) * u &
      + 72.81777742e+00_real64 ) * u &
      - 113.77777774e+00_real64 ) * u &
      + 16.0e+00_real64 )

    ger = (((((( &
      - 0.2458e-04_real64 * u &
      + 0.309699e-02_real64 ) * u &
      - 0.19636347e+00_real64 ) * u &
      + 5.65539121e+00_real64 ) * u &
      - 60.60977451e+00_real64 ) * u &
      + 171.36272133e+00_real64 ) * u &
      - 59.05819744e+00_real64 ) * u &
      - 0.57721566e+00_real64

    ger = ger - log ( 0.5e+00_real64 * x ) * ber + 0.25e+00_real64 * pi * bei

    gei = t2 * (((((( &
        0.29532e-03_real64 * u &
      - 0.02695875e+00_real64 ) * u &
      + 1.17509064e+00_real64 ) * u &
      - 21.30060904e+00_real64 ) * u &
      + 124.2356965e+00_real64 ) * u &
      - 142.91827687e+00_real64 ) * u &
      + 6.76454936e+00_real64 )

    gei = gei - log ( 0.5e+00_real64 * x ) * bei - 0.25e+00_real64 * pi * ber

    der = x * t2 * (((((( &
      - 0.394e-05_real64 * u &
      + 0.45957e-03_real64 ) * u &
      - 0.02609253e+00_real64 ) * u &
      + 0.66047849e+00_real64 ) * u &
      - 6.0681481e+00_real64 ) * u &
      + 14.22222222e+00_real64 ) * u &
      - 4.0e+00_real64 )

    dei = x * (((((( &
        0.4609e-04_real64 * u &
      - 0.379386e-02_real64 ) * u &
      + 0.14677204e+00_real64 ) * u &
      - 2.31167514e+00_real64 ) * u &
      + 11.37777772e+00_real64 ) * u &
      - 10.66666666e+00_real64 ) * u &
      + 0.5e+00_real64 ) 

    her = x * t2 * (((((( &
      - 0.1075e-04_real64 * u &
      + 0.116137e-02_real64 ) * u &
      - 0.06136358e+00_real64 ) * u &
      + 1.4138478e+00_real64 ) * u &
      - 11.36433272e+00_real64 ) * u &
      + 21.42034017e+00_real64 ) * u &
      - 3.69113734e+00_real64 )

    her = her - log ( 0.5e+00_real64 * x ) * der - ber / x  &
      + 0.25e+00_real64 * pi * dei

    hei = x * (((((( &
        0.11997e-03_real64 * u &
      - 0.926707e-02_real64 ) * u &
      + 0.33049424e+00_real64 ) * u &
      - 4.65950823e+00_real64 ) * u &
      + 19.41182758e+00_real64 ) * u &
      - 13.39858846e+00_real64 ) * u &
      + 0.21139217e+00_real64 )

    hei = hei - log ( 0.5e+00_real64 * x ) * dei - bei / x  &
      - 0.25e+00_real64 * pi * der

  else

    t = 8.0e+00_real64 / x

    do l = 1, 2

      v = ( -1.0e+00_real64 ) ** l * t

      tpr = (((( &
          0.6e-06_real64 * v &
        - 0.34e-05_real64 ) * v &
        - 0.252e-04_real64 ) * v &
        - 0.906e-04_real64 ) * v * v &
        + 0.0110486e+00_real64 ) * v

      tpi = (((( &
          0.19e-05_real64 * v &
        + 0.51e-05_real64 ) * v * v &
        - 0.901e-04_real64 ) * v &
        - 0.9765e-03_real64 ) * v &
        - 0.0110485e+00_real64 ) * v &
        - 0.3926991e+00_real64

      if ( l == 1 ) then
        tnr = tpr
        tni = tpi
      end if

    end do

    yd = x / sqrt ( 2.0e+00_real64 )
    ye1 = exp ( yd + tpr )
    ye2 = exp ( - yd + tnr )
    yc1 = 1.0e+00_real64 / sqrt ( 2.0e+00_real64 * pi * x )
    yc2 = sqrt ( pi / ( 2.0e+00_real64 * x ) )
    csp = cos ( yd + tpi )
    ssp = sin ( yd + tpi )
    csn = cos ( - yd + tni )
    ssn = sin ( - yd + tni )
    ger = yc2 * ye2 * csn
    gei = yc2 * ye2 * ssn
    fxr = yc1 * ye1 * csp
    fxi = yc1 * ye1 * ssp
    ber = fxr - gei / pi
    bei = fxi + ger / pi

    do l = 1, 2

      v = ( -1.0e+00_real64 ) ** l * t

      ppr = ((((( &
          0.16e-05_real64 * v &
        + 0.117e-04_real64 ) * v &
        + 0.346e-04_real64 ) * v &
        + 0.5e-06_real64 ) * v &
        - 0.13813e-02_real64 ) * v &
        - 0.0625001e+00_real64 ) * v &
        + 0.7071068e+00_real64

      ppi = ((((( &
        - 0.32e-05_real64 * v &
        - 0.24e-05_real64 ) * v &
        + 0.338e-04_real64 ) * v &
        + 0.2452e-03_real64 ) * v &
        + 0.13811e-02_real64 ) * v &
        - 0.1e-06_real64 ) * v &
        + 0.7071068e+00_real64

      if ( l == 1 ) then
        pnr = ppr
        pni = ppi
      end if

    end do

    her =     gei * pni - ger * pnr
    hei = - ( gei * pnr + ger * pni )
    der = fxr * ppr - fxi * ppi - hei / pi
    dei = fxi * ppr + fxr * ppi + her / pi

  end if

  return
end subroutine klvnb
!> @brief subroutine klvnzo.
!> @return None.
!>
!> @param nt [in] Argument nt.
!> @param kd [in] Argument kd.
!> @param zo [inout] Argument zo.
subroutine klvnzo ( nt, kd, zo )

!*****************************************************************************80
!
!! KLVNZO computes zeros of the Kelvin functions.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    15 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) NT, the number of zeros.
!
!    Input, integer(int32) KD, the function code.
!    1 for ber x, 
!    2 for bei x,
!    3 for ker x, 
!    4 for kei x,
!    5 for ber' x, 
!    6 for bei' x,
!    7 for ker' x, 
!    8 for kei' x.
!
!    Output, real(real64) ZO(NT), the zeros of the given Kelvin function.
!
  implicit none

  integer(int32), intent(in) :: nt

  real(real64) bei
  real(real64) ber
  real(real64) ddi
  real(real64) ddr
  real(real64) dei
  real(real64) der
  real(real64) gdi
  real(real64) gdr
  real(real64) gei
  real(real64) ger
  real(real64) hei
  real(real64) her
  integer(int32), intent(in) :: kd
  integer(int32) m
  real(real64) rt
  real(real64) rt0(8)
  real(real64), intent(inout) :: zo(nt)

  rt0(1) = 2.84891e+00_real64
  rt0(2) = 5.02622e+00_real64
  rt0(3) = 1.71854e+00_real64
  rt0(4) = 3.91467e+00_real64
  rt0(5) = 6.03871e+00_real64
  rt0(6) = 3.77268e+00_real64
  rt0(7) = 2.66584e+00_real64
  rt0(8) = 4.93181e+00_real64

  rt = rt0(kd)

  do m = 1, nt

    do

      call klvna ( rt, ber, bei, ger, gei, der, dei, her, hei )

      if ( kd == 1 ) then
        rt = rt - ber / der
      else if ( kd == 2 ) then
        rt = rt - bei / dei
      else if ( kd == 3 ) then
        rt = rt - ger / her
      else if ( kd == 4 ) then
        rt = rt - gei / hei
      else if ( kd == 5 ) then
        ddr = - bei - der / rt
        rt = rt - der / ddr
      else if ( kd == 6 ) then
        ddi = ber - dei / rt
        rt = rt - dei / ddi
      else if ( kd == 7 ) then
        gdr = - gei - her / rt
        rt = rt - her / gdr
      else
        gdi = ger - hei / rt
        rt = rt - hei / gdi
      end if

      if ( abs ( rt - rt0(kd) ) <= 5.0e-10_real64 ) then
        exit
      end if

      rt0(kd) = rt

    end do

    zo(m) = rt
    rt = rt + 4.44e+00_real64

  end do

  return
end subroutine klvnzo
!> @brief subroutine kmn.
!> @return None.
!>
!> @param m [in] Argument m.
!> @param n [in] Argument n.
!> @param c [in] Argument c.
!> @param cv [in] Argument cv.
!> @param kd [in] Argument kd.
!> @param df [in] Argument df.
!> @param dn [inout] Argument dn.
!> @param ck1 [inout] Argument ck1.
!> @param ck2 [inout] Argument ck2.
pure subroutine kmn ( m, n, c, cv, kd, df, dn, ck1, ck2 )

!*****************************************************************************80
!
!! KMN: expansion coefficients of prolate or oblate spheroidal functions.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    02 August 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) M, the mode parameter;  M = 0, 1, 2, ...
!
!    Input, integer(int32) N, mode parameter, N = M, M + 1, M + 2, ...
!
!    Input, real(real64) C, spheroidal parameter.
!
!    Input, real(real64) CV, the characteristic value.
!
!    Input, integer(int32) KD, the function code.
!    1, the prolate function.
!    -1, the oblate function.
!
!    Input, real(real64) DF(*), the expansion coefficients.
!
  implicit none

  real(real64), intent(in) :: c
  real(real64), intent(inout) :: ck1
  real(real64), intent(inout) :: ck2
  real(real64) cs
  real(real64), intent(in) :: cv
  real(real64), intent(in) :: df(200)
  real(real64), intent(inout) :: dn(200)
  real(real64) dnp
  real(real64) g0
  real(real64) gk0
  real(real64) gk1
  real(real64) gk2
  real(real64) gk3
  integer(int32) i
  integer(int32) ip
  integer(int32) j
  integer(int32) k
  integer(int32), intent(in) :: kd
  integer(int32) l
  integer(int32), intent(in) :: m
  integer(int32), intent(in) :: n
  integer(int32) nm
  integer(int32) nm1
  integer(int32) nn
  real(real64) r
  real(real64) r1
  real(real64) r2
  real(real64) r3
  real(real64) r4
  real(real64) r5
  real(real64) rk(200)
  real(real64) sa0
  real(real64) sb0
  real(real64) su0
  real(real64) sw
  real(real64) t
  real(real64) tp(200)
  real(real64) u(200)
  real(real64) v(200)
  real(real64) w(200)

  nm = 25 + int ( 0.5e+00_real64 * ( n - m ) + c )
  nn = nm + m
  cs = c * c * kd

  if ( n - m == 2 * int ( ( n - m ) / 2 ) ) then
    ip = 0
  else
    ip = 1
  end if

  do i = 1, nn + 3

    if ( ip == 0 ) then
      k = - 2 * ( i - 1 )
    else
      k = - ( 2 * i - 3 )
    end if

    gk0 = 2.0e+00_real64 * m + k
    gk1 = ( m + k ) * ( m + k + 1.0e+00_real64 )
    gk2 = 2.0e+00_real64 * ( m + k ) - 1.0e+00_real64
    gk3 = 2.0e+00_real64 * ( m + k ) + 3.0e+00_real64
    u(i) = gk0 * ( gk0 - 1.0e+00_real64 ) * cs / ( gk2 * ( gk2 + 2.0e+00_real64 ) )
    v(i) = gk1 - cv + ( 2.0e+00_real64 * ( gk1 - m * m ) - 1.0e+00_real64 ) * cs &
      / ( gk2 * gk3 )
    w(i) = ( k + 1.0e+00_real64 ) * ( k + 2.0e+00_real64 ) * cs / ( ( gk2 + 2.0e+00_real64 ) * gk3 )

  end do

  do k = 1, m
    t = v(m+1)
    do l = 0, m - k - 1
      t = v(m-l) - w(m-l+1) * u(m-l) / t
    end do
    rk(k) = -u(k) / t
  end do

  r = 1.0e+00_real64
  do k = 1, m
    r = r * rk(k)
    dn(k) = df(1) * r
  end do

  tp(nn) = v(nn+1)
  do k = nn - 1, m + 1,-1
    tp(k) = v(k+1) - w(k+2) * u(k+1) / tp(k+1)
    if ( m + 1 < k ) then
      rk(k) = -u(k) / tp(k)
    end if
  end do

  if ( m == 0 ) then
    dnp = df(1)
  else
    dnp = dn(m)
  end if

  dn(m+1) = ( - 1.0e+00_real64 ) ** ip * dnp * cs &
    / ( ( 2.0e+00_real64 * m - 1.0e+00_real64 ) &
    * ( 2.0e+00_real64 * m + 1.0e+00_real64 - 4.0e+00_real64 * ip ) * tp(m+1) )
  do k = m + 2, nn
    dn(k) = rk(k) * dn(k-1)
  end do

  r1 = 1.0e+00_real64
  do j = 1, ( n + m + ip ) / 2
       r1 = r1 * ( j + 0.5e+00_real64 * ( n + m + ip ) )
  end do
  nm1 = ( n - m ) / 2
  r = 1.0e+00_real64
  do j = 1, 2 * m + ip
    r = r * j
  end do
  su0 = r * df(1)
  sw = 0.0e+00_real64

  do k = 2, nm
    r = r * ( m + k - 1.0e+00_real64 ) * ( m + k + ip - 1.5e+00_real64 ) &
      / ( k - 1.0e+00_real64 ) / ( k + ip - 1.5e+00_real64 )
    su0 = su0 + r * df(k)
    if ( nm1 < k .and. &
      abs ( ( su0 - sw ) / su0 ) < 1.0e-14_real64 ) then
      exit
    end if
    sw = su0
  end do

  if ( kd /= 1 ) then

    r2 = 1.0e+00_real64
    do j = 1,m
      r2 = 2.0e+00_real64 * c * r2 * j
    end do
    r3 = 1.0e+00_real64
    do j = 1, ( n - m - ip ) / 2
      r3 = r3 * j
    end do
    sa0 = ( 2.0e+00_real64 * ( m + ip ) + 1.0e+00_real64 ) * r1 &
      / ( 2.0e+00_real64 ** n * c ** ip * r2 * r3 * df(1) )
    ck1 = sa0 * su0

    if ( kd == -1 ) then
      return
    end if

  end if

  r4 = 1.0e+00_real64
  do j = 1, ( n - m - ip ) / 2
    r4 = 4.0e+00_real64 * r4 * j
  end do
  r5 = 1.0e+00_real64
  do j = 1, m
    r5 = r5 * ( j + m ) / c
  end do

  if ( m == 0 ) then
    g0 = df(1)
  else
    g0 = dn(m)
  end if

  sb0 = ( ip + 1.0e+00_real64 ) * c ** ( ip + 1 ) &
    / ( 2.0e+00_real64 * ip * ( m - 2.0e+00_real64 ) + 1.0e+00_real64 ) &
    / ( 2.0e+00_real64 * m - 1.0e+00_real64 )

  ck2 = ( -1 ) ** ip * sb0 * r4 * r5 * g0 / r1 * su0

  return
end subroutine kmn
!> @brief subroutine lagzo.
!> @return None.
!>
!> @param n [in] Argument n.
!> @param x [inout] Argument x.
!> @param w [inout] Argument w.
pure subroutine lagzo ( n, x, w )

!*****************************************************************************80
!
!! LAGZO computes zeros of the Laguerre polynomial, and integration weights.
!
!  Discussion:
!
!    This procedure computes the zeros of Laguerre polynomial Ln(x) in the 
!    interval [0,�], and the corresponding weighting coefficients for 
!    Gauss-Laguerre integration.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    07 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) N, the order of the Laguerre polynomial.
!
!    Output, real(real64) X(N), the zeros of the Laguerre polynomial.
!
!    Output, real(real64) W(N), the weighting coefficients.
!
  implicit none

  integer(int32), intent(in) :: n

  real(real64) f0
  real(real64) f1
  real(real64) fd
  real(real64) gd
  real(real64) hn
  integer(int32) i
  integer(int32) it
  integer(int32) j
  integer(int32) k
  integer(int32) nr
  real(real64) p
  real(real64) pd
  real(real64) pf
  real(real64) q
  real(real64), intent(inout) :: w(n)
  real(real64) wp
  real(real64), intent(inout) :: x(n)
  real(real64) z
  real(real64) z0

  hn = 1.0e+00_real64 / real ( n, kind = real64 )

  do nr = 1, n

    if ( nr == 1 ) then
      z = hn
    else
      z = x(nr-1) + hn * nr ** 1.27e+00_real64
    end if

    it = 0

    do

      it = it + 1
      z0 = z
      p = 1.0e+00_real64
      do i = 1, nr - 1
        p = p * ( z - x(i) )
      end do

      f0 = 1.0e+00_real64
      f1 = 1.0e+00_real64 - z
      do k = 2, n
        pf = (( 2.0e+00_real64 * k - 1.0e+00_real64 - z ) * f1 &
          - ( k - 1.0e+00_real64 ) * f0 ) / k
        pd = k / z * ( pf - f1 )
        f0 = f1
        f1 = pf
      end do

      fd = pf / p

      q = 0.0e+00_real64
      do i = 1, nr - 1
        wp = 1.0e+00_real64
        do j = 1, nr - 1
          if ( j /= i ) then
            wp = wp * ( z - x(j) )
          end if
        end do
        q = q + wp
      end do

      gd = ( pd - q * fd ) / p
      z = z - fd / gd

      if ( 40 < it .or. abs ( ( z - z0 ) / z ) <= 1.0e-15_real64 ) then
        exit
      end if

    end do

    x(nr) = z
    w(nr) = 1.0e+00_real64 / ( z * pd * pd )

  end do

  return
end subroutine lagzo
!> @brief subroutine lamn.
!> @return None.
!>
!> @param n [in] Argument n.
!> @param x [in] Argument x.
!> @param nm [inout] Argument nm.
!> @param bl [inout] Argument bl.
!> @param dl [inout] Argument dl.
subroutine lamn ( n, x, nm, bl, dl )

!*****************************************************************************80
!
!! LAMN computes lambda functions and derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    14 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) N, the order.
!
!    Input, real(real64) X, the argument.
!
!    Output, integer(int32) NM, the highest order computed.
!
!    Output, real(real64) BL(0:N), DL(0:N), the
!    value of the lambda function and its derivative of orders 0 through N.
!
  implicit none

  integer(int32), intent(in) :: n

  real(real64) bg
  real(real64) bk
  real(real64), intent(inout) :: bl(0:n)
  real(real64) bs
  real(real64), intent(inout) :: dl(0:n)
  real(real64) f
  real(real64) f0
  real(real64) f1
  integer(int32) i
  integer(int32) k
  integer(int32) m
  ! integer(int32) msta1
  ! integer(int32) msta2
  integer(int32), intent(inout) :: nm
  real(real64) r
  real(real64) r0
  real(real64) uk
  real(real64), intent(in) :: x
  real(real64) x2

  nm = n

  if ( abs ( x ) < 1.0e-100_real64 ) then
    do k = 0, n
      bl(k) = 0.0e+00_real64
      dl(k) = 0.0e+00_real64
    end do
    bl(0) = 1.0e+00_real64
    dl(1) = 0.5e+00_real64
    return
  end if

  if ( x <= 12.0e+00_real64 ) then

    x2 = x * x

    do k = 0, n
      bk = 1.0e+00_real64
      r = 1.0e+00_real64
      do i = 1, 50
        r = -0.25e+00_real64 * r * x2 / ( i * ( i + k ) )
        bk = bk + r
        if ( abs ( r ) < abs ( bk ) * 1.0e-15_real64 ) then
          exit
        end if
      end do

      bl(k) = bk
      if ( 1 <= k ) then
        dl(k-1) = - 0.5e+00_real64 * x / k * bk
      end if

    end do

    uk = 1.0e+00_real64
    r = 1.0e+00_real64
    do i = 1, 50
      r = -0.25e+00_real64 * r * x2 / ( i * ( i + n + 1.0e+00_real64 ) )
      uk = uk + r
      if ( abs ( r ) < abs ( uk ) * 1.0e-15_real64 ) then
        exit
      end if
    end do

    dl(n) = -0.5e+00_real64 * x / ( n + 1.0e+00_real64 ) * uk
    return

  end if

  if ( n == 0 ) then
    nm = 1
  end if

  m = msta1 ( x, 200 )

  if ( m < nm ) then
    nm = m
  else
    m = msta2 ( x, nm, 15 )
  end if

  bs = 0.0e+00_real64
  f0 = 0.0e+00_real64
  f1 = 1.0e-100_real64
  do k = m, 0, -1
    f = 2.0e+00_real64 * ( k + 1.0e+00_real64 ) * f1 / x - f0
    if ( k <= nm ) then
      bl(k) = f
    end if
    if ( k == 2 * int ( k / 2 ) ) then
      bs = bs + 2.0e+00_real64 * f
    end if
    f0 = f1
    f1 = f
  end do

  bg = bs - f
  do k = 0, nm
    bl(k) = bl(k) / bg
  end do

  r0 = 1.0e+00_real64
  do k = 1, nm
    r0 = 2.0e+00_real64 * r0 * k / x
    bl(k) = r0 * bl(k)
  end do

  dl(0) = -0.5e+00_real64 * x * bl(1)
  do k = 1, nm
    dl(k) = 2.0e+00_real64 * k / x * ( bl(k-1) - bl(k) )
  end do

  return
end subroutine lamn
!> @brief subroutine lamv.
!> @return None.
!>
!> @param v [in] Argument v.
!> @param x [inout] Argument x.
!> @param vm [inout] Argument vm.
!> @param vl [inout] Argument vl.
!> @param dl [inout] Argument dl.
subroutine lamv ( v, x, vm, vl, dl )

!*****************************************************************************80
!
!! LAMV computes lambda functions and derivatives of arbitrary order.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    31 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) V, the order.
!
!    Input, real(real64) X, the argument.
!
!    Output, real(real64) VM, the highest order computed.
!
!    Output, real(real64) VL(0:*), DL(0:*), the Lambda function and 
!    derivative, of orders N+V0.
!
  implicit none

  real(real64), intent(in) :: v

  real(real64) a0
  real(real64) bjv0
  real(real64) bjv1
  real(real64) bk
  real(real64) ck
  real(real64) cs
  real(real64), intent(inout) :: dl(0:int(v))
  real(real64) f
  real(real64) f0
  real(real64) f1
  real(real64) f2
  real(real64) fac
  real(real64) ga
  integer(int32) i
  integer(int32) j
  integer(int32) k
  integer(int32) k0
  integer(int32) m
  ! integer(int32) msta1
  ! integer(int32) msta2
  integer(int32) n
  real(real64) pi
  real(real64) px
  real(real64) qx
  real(real64) r
  real(real64) r0
  real(real64) rc
  real(real64) rp
  real(real64) rp2
  real(real64) rq
  real(real64) sk
  real(real64) uk
  real(real64) v0
  real(real64) vk
  real(real64), intent(inout) :: vl(0:int(v))
  real(real64), intent(inout) :: vm
  real(real64) vv
  real(real64), intent(inout) :: x
  real(real64) x2
  real(real64) xk

  pi = 3.141592653589793e+00_real64
  rp2 = 0.63661977236758e+00_real64
  x = abs ( x )
  x2 = x * x
  n = int ( v )
  v0 = v - n
  vm = v

  if ( x <= 12.0e+00_real64 ) then

    do k = 0, n

      vk = v0 + k
      bk = 1.0e+00_real64
      r = 1.0e+00_real64

      do i = 1, 50
        r = -0.25e+00_real64 * r * x2 / ( i * ( i + vk ) )
        bk = bk + r
        if ( abs ( r ) < abs ( bk ) * 1.0e-15_real64 ) then
          exit
        end if
      end do

      vl(k) = bk
      uk = 1.0e+00_real64
      r = 1.0e+00_real64
      do i = 1, 50
        r = -0.25e+00_real64 * r * x2 / ( i * ( i + vk + 1.0e+00_real64 ))
        uk = uk + r
        if ( abs ( r ) < abs ( uk ) * 1.0e-15_real64 ) then
          exit
        end if
      end do

      dl(k) = - 0.5e+00_real64 * x / ( vk + 1.0e+00_real64 ) * uk

    end do

    return

  end if

  if ( x < 35.0e+00_real64 ) then
    k0 = 11
  else if ( x < 50.0e+00_real64 ) then
    k0 = 10
  else
    k0 = 8
  end if

  do j = 0, 1
    vv = 4.0e+00_real64 * ( j + v0 ) * ( j + v0 )
    px = 1.0e+00_real64
    rp = 1.0e+00_real64
    do k = 1, k0
      rp = - 0.78125e-02_real64 * rp * ( vv - ( 4.0e+00_real64 * k - 3.0e+00_real64 ) ** 2 ) &
        * ( vv - ( 4.0e+00_real64 * k - 1.0e+00_real64 ) ** 2 ) &
        / ( k * ( 2.0_real64 * k - 1.0e+00_real64 ) * x2 )
      px = px + rp
    end do
    qx = 1.0e+00_real64
    rq = 1.0e+00_real64
    do k = 1, k0
      rq = - 0.78125e-02_real64 * rq * ( vv - ( 4.0e+00_real64 * k - 1.0e+00_real64 ) ** 2 ) &
        * ( vv - ( 4.0e+00_real64 * k + 1.0e+00_real64 ) ** 2 ) &
        / ( k * ( 2.0e+00_real64 * k + 1.0e+00_real64 ) * x2 )
      qx = qx + rq
    end do
    qx = 0.125e+00_real64 * ( vv - 1.0e+00_real64 ) * qx / x
    xk = x - ( 0.5e+00_real64 * ( j + v0 ) + 0.25e+00_real64 ) * pi
    a0 = sqrt ( rp2 / x )
    ck = cos ( xk )
    sk = sin ( xk )
    if ( j == 0 ) then
      bjv0 = a0 * ( px * ck - qx * sk )
    else
      bjv1 = a0 * ( px * ck - qx * sk )
    end if
  end do

  if ( v0 == 0.0e+00_real64 ) then
    ga = 1.0e+00_real64
  else
    call gam0 ( v0, ga )
    ga = v0 * ga
  end if

  fac = ( 2.0e+00_real64 / x ) ** v0 * ga
  vl(0) = bjv0
  dl(0) = - bjv1 + v0 / x * bjv0
  vl(1) = bjv1
  dl(1) = bjv0 - ( 1.0e+00_real64 + v0 ) / x * bjv1
  r0 = 2.0e+00_real64 * ( 1.0e+00_real64 + v0 ) / x

  if ( n <= 1 ) then
    vl(0) = fac * vl(0)
    dl(0) = fac * dl(0) - v0 / x * vl(0)
    vl(1) = fac * r0 * vl(1)
    dl(1) = fac * r0 * dl(1) - ( 1.0e+00_real64 + v0 ) / x * vl(1)
    return
  end if

  if ( 2 <= n .and. n <= int ( 0.9e+00_real64 * x ) ) then

    f0 = bjv0
    f1 = bjv1
    do k = 2, n
      f = 2.0e+00_real64 * ( k + v0 - 1.0e+00_real64 ) / x * f1 - f0
      f0 = f1
      f1 = f
      vl(k) = f
    end do

  else if ( 2 <= n ) then

    m = msta1 ( x, 200 )
    if ( m < n ) then
      n = m
    else
      m = msta2 ( x, n, 15 )
    end if
    f2 = 0.0e+00_real64
    f1 = 1.0e-100_real64
    do k = m, 0, -1
      f = 2.0e+00_real64 * ( v0 + k + 1.0e+00_real64 ) / x * f1 - f2
      if ( k <= n ) then
        vl(k) = f
      end if
      f2 = f1
      f1 = f
    end do

    if ( abs ( bjv0 ) <= abs ( bjv1 ) ) then
      cs = bjv1 / f2
    else
      cs = bjv0 / f
    end if

    do k = 0, n
      vl(k) = cs * vl(k)
    end do

  end if

  vl(0) = fac * vl(0)
  do j = 1, n
    rc = fac * r0
    vl(j) = rc * vl(j)
    dl(j-1) = - 0.5e+00_real64 * x / ( j + v0 ) * vl(j)
    r0 = 2.0e+00_real64 * ( j + v0 + 1 ) / x * r0
  end do
  dl(n) = 2.0e+00_real64 * ( v0 + n ) * ( vl(n-1) - vl(n) ) / x
  vm = n + v0

  return
end subroutine lamv
!> @brief subroutine legzo.
!> @return None.
!>
!> @param n [in] Argument n.
!> @param x [inout] Argument x.
!> @param w [inout] Argument w.
pure subroutine legzo ( n, x, w )

!*****************************************************************************80
!
!! LEGZO computes the zeros of Legendre polynomials, and integration weights.
!
!  Discussion:
!
!    This procedure computes the zeros of Legendre polynomial Pn(x) in the 
!    interval [-1,1], and the corresponding weighting coefficients for 
!    Gauss-Legendre integration.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    13 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) N, the order of the polynomial.
!
!    Output, real(real64) X(N), W(N), the zeros of the polynomial,
!    and the corresponding weights.
!
  implicit none

  integer(int32), intent(in) :: n

  real(real64) f0
  real(real64) f1
  real(real64) fd
  real(real64) gd
  integer(int32) i
  integer(int32) j
  integer(int32) k
  integer(int32) n0
  integer(int32) nr
  real(real64) p
  real(real64) pd
  real(real64) pf
  real(real64) q
  real(real64), intent(inout) :: w(n)
  real(real64) wp
  real(real64), intent(inout) :: x(n)
  real(real64) z
  real(real64) z0

  n0 = ( n + 1 ) / 2

  do nr = 1, n0

    z = cos ( 3.1415926e+00_real64 * ( nr - 0.25e+00_real64 ) / n )

    do

      z0 = z
      p = 1.0e+00_real64
      do i = 1, nr - 1
        p = p * ( z - x(i))
      end do
      f0 = 1.0e+00_real64
      if ( nr == n0 .and. n /= 2 * int ( n / 2 ) ) then
        z = 0.0e+00_real64
      end if
      f1 = z
      do k = 2, n
        pf = ( 2.0e+00_real64 - 1.0e+00_real64 / k ) * z * f1 &
          - ( 1.0e+00_real64 - 1.0e+00_real64 / k ) * f0
        pd = k * ( f1 - z * pf ) / ( 1.0e+00_real64 - z * z )
        f0 = f1
        f1 = pf
      end do

      if ( z == 0.0e+00_real64 ) then
        exit
      end if

      fd = pf / p
      q = 0.0e+00_real64
      do i = 1, nr - 1
        wp = 1.0e+00_real64
        do j = 1, nr - 1
          if ( j /= i ) then
            wp = wp * ( z - x(j) )
          end if
        end do
        q = q + wp
      end do
      gd = ( pd - q * fd ) / p
      z = z - fd / gd

      if ( abs ( z - z0 ) < abs ( z ) * 1.0e-15_real64 ) then
        exit
      end if

    end do

    x(nr) = z
    x(n+1-nr) = - z
    w(nr) = 2.0e+00_real64 / ( ( 1.0e+00_real64 - z * z ) * pd * pd )
    w(n+1-nr) = w(nr)

  end do

  return
end subroutine legzo
!> @brief subroutine lgama.
!> @return None.
!>
!> @param kf [in] Argument kf.
!> @param x [in] Argument x.
!> @param gl [inout] Argument gl.
pure subroutine lgama ( kf, x, gl )

!*****************************************************************************80
!
!! LGAMA computes the gamma function or its logarithm.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    15 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) KF, the argument code.
!    1, for gamma(x);
!    2, for ln(gamma(x)).
!
!    Input, real(real64) X, the argument.
!
!    Output, real(real64) GL, the function value.
!
  implicit none

  real(real64), parameter, dimension ( 10 ) :: a = [&
    8.333333333333333e-02_real64, &
   -2.777777777777778e-03_real64, &
    7.936507936507937e-04_real64, &
   -5.952380952380952e-04_real64, &
    8.417508417508418e-04_real64, &
   -1.917526917526918e-03_real64, &
    6.410256410256410e-03_real64, &
   -2.955065359477124e-02_real64, &
    1.796443723688307e-01_real64, &
   -1.39243221690590e+00_real64]
  real(real64), intent(inout) :: gl
  real(real64) gl0
  integer(int32) k
  integer(int32), intent(in) :: kf
  integer(int32) n
  real(real64), intent(in) :: x
  real(real64) x0
  real(real64) x2
  real(real64) xp

  x0 = x

  if ( x == 1.0e+00_real64 .or. x == 2.0e+00_real64 ) then
    gl = 0.0e+00_real64
    if ( kf == 1 ) then
      gl = 1.0e+00_real64
    end if
    return
  else if ( x <= 7.0e+00_real64 ) then
    n = int ( 7.0e+00_real64 - x )
    x0 = x + n
  end if

  x2 = 1.0e+00_real64 / ( x0 * x0 )
  xp = 6.283185307179586477e+00_real64
  gl0 = a(10)

  do k = 9, 1, -1
    gl0 = gl0 * x2 + a(k)
  end do

  gl = gl0 / x0 + 0.5e+00_real64 * log ( xp ) + ( x0 - 0.5e+00_real64 ) * log ( x0 ) - x0

  if ( x <= 7.0e+00_real64 ) then
    do k = 1, n
      gl = gl - log ( x0 - 1.0e+00_real64 )
      x0 = x0 - 1.0e+00_real64
    end do
  end if

  if ( kf == 1 ) then
    gl = exp ( gl )
  end if

  return
end subroutine lgama
!> @brief subroutine lpmn.
!> @return None.
!>
!> @param mm [in] Argument mm.
!> @param m [in] Argument m.
!> @param n [in] Argument n.
!> @param x [in] Argument x.
!> @param pm [inout] Argument pm.
!> @param pd [inout] Argument pd.
pure subroutine lpmn ( mm, m, n, x, pm, pd )

!*****************************************************************************80
!
!! LPMN computes associated Legendre functions Pmn(X) and derivatives P'mn(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    19 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) MM, the leading dimension of PM and PD.
!
!    Input, integer(int32) M, the order of Pmn(x).
!
!    Input, integer(int32) N, the degree of Pmn(x).
!
!    Input, real(real64) X, the argument of Pmn(x).
!
!    Output, real(real64) PM(0:MM,0:N), PD(0:MM,0:N), the
!    values of Pmn(x) and Pmn'(x).
!
  implicit none

  integer(int32), intent(in) :: mm
  integer(int32), intent(in) :: n

  integer(int32) i
  integer(int32) j
  integer(int32) ls
  integer(int32), intent(in) :: m
  real(real64), intent(inout) :: pd(0:mm,0:n)
  real(real64), intent(inout) :: pm(0:mm,0:n)
  real(real64), intent(in) :: x
  real(real64) xq
  real(real64) xs

  do i = 0, n
    do j = 0, m
      pm(j,i) = 0.0e+00_real64
      pd(j,i) = 0.0e+00_real64
    end do
  end do

  pm(0,0) = 1.0e+00_real64

  if ( abs ( x ) == 1.0e+00_real64 ) then

    do i = 1, n
      pm(0,i) = x ** i
      pd(0,i) = 0.5e+00_real64 * i * ( i + 1.0e+00_real64 ) * x ** ( i + 1 )
    end do

    do j = 1, n
      do i = 1, m
        if ( i == 1 ) then
          pd(i,j) = 1.0e+300_real64
        else if ( i == 2 ) then
          pd(i,j) = -0.25e+00_real64 * ( j + 2 ) * ( j + 1 ) * j &
            * ( j - 1 ) * x ** ( j + 1 )
        end if
      end do
    end do

    return

  end if

  if ( 1.0e+00_real64 < abs ( x ) ) then
    ls = -1
  else
    ls = +1
  end if

  xq = sqrt ( ls * ( 1.0e+00_real64 - x * x ) )
  xs = ls * ( 1.0e+00_real64 - x * x )
  do i = 1, m
    pm(i,i) = - ls * ( 2.0e+00_real64 * i - 1.0e+00_real64 ) * xq * pm(i-1,i-1)
  end do

  do i = 0, m
    pm(i,i+1) = ( 2.0e+00_real64 * i + 1.0e+00_real64 ) * x * pm(i,i)
  end do

  do i = 0, m
    do j = i + 2, n
      pm(i,j) = ( ( 2.0e+00_real64 * j - 1.0e+00_real64 ) * x * pm(i,j-1) - &
        ( i + j - 1.0e+00_real64 ) * pm(i,j-2) ) / ( j - i )
    end do
  end do

  pd(0,0) = 0.0e+00_real64
  do j = 1, n
    pd(0,j) = ls * j * ( pm(0,j-1) - x * pm(0,j) ) / xs
  end do

  do i = 1, m
    do j = i, n
      pd(i,j) = ls * i * x * pm(i,j) / xs + ( j + i ) &
        * ( j - i + 1.0e+00_real64 ) / xq * pm(i-1,j)
    end do
  end do

  return
end subroutine lpmn
!> @brief subroutine lpmns.
!> @return None.
!>
!> @param m [in] Argument m.
!> @param n [in] Argument n.
!> @param x [in] Argument x.
!> @param pm [inout] Argument pm.
!> @param pd [inout] Argument pd.
pure subroutine lpmns ( m, n, x, pm, pd )

!*****************************************************************************80
!
!! LPMNS computes associated Legendre functions Pmn(X) and derivatives P'mn(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    18 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) M, the order of Pmn(x).
!
!    Input, integer(int32) N, the degree of Pmn(x).
!
!    Input, real(real64) X, the argument.
!
!    Output, real(real64) PM(0:N), PD(0:N), the values and derivatives
!    of the function from degree 0 to N.
!
  implicit none

  integer(int32), intent(in) :: n

  integer(int32) k
  integer(int32), intent(in) :: m
  real(real64), intent(inout) :: pm(0:n)
  real(real64) pm0
  real(real64) pm1
  real(real64) pm2
  real(real64) pmk
  real(real64), intent(inout) :: pd(0:n)
  real(real64), intent(in) :: x
  real(real64) x0

  do k = 0, n
    pm(k) = 0.0e+00_real64
    pd(k) = 0.0e+00_real64
  end do

  if ( abs ( x ) == 1.0e+00_real64 ) then

    do k = 0, n
      if ( m == 0 ) then
        pm(k) = 1.0e+00_real64
        pd(k) = 0.5e+00_real64 * k * ( k + 1.0e+00_real64 )
        if ( x < 0.0e+00_real64 ) then
          pm(k) = ( -1.0e+00_real64 ) ** k * pm(k)
          pd(k) = ( -1.0e+00_real64 ) ** ( k + 1 ) * pd(k)
        end if
      else if ( m == 1 ) then
        pd(k) = 1.0e+300_real64
      else if ( m == 2 ) then
        pd(k) = -0.25e+00_real64 * ( k + 2.0e+00_real64 ) * ( k + 1.0e+00_real64 ) &
          * k * ( k - 1.0e+00_real64 )
        if ( x < 0.0e+00_real64 ) then
          pd(k) = ( -1.0e+00_real64 ) ** ( k + 1 ) * pd(k)
        end if
      end if
    end do
    return
  end if

  x0 = abs ( 1.0e+00_real64 - x * x )
  pm0 = 1.0e+00_real64
  pmk = pm0
  do k = 1, m
    pmk = ( 2.0e+00_real64 * k - 1.0e+00_real64 ) * sqrt ( x0 ) * pm0
    pm0 = pmk
  end do
  pm1 = ( 2.0e+00_real64 * m + 1.0e+00_real64 ) * x * pm0
  pm(m) = pmk
  pm(m+1) = pm1
  do k = m + 2, n
    pm2 = ( ( 2.0e+00_real64 * k - 1.0e+00_real64 ) * x * pm1 &
      - ( k + m - 1.0e+00_real64 ) * pmk ) / ( k - m )
    pm(k) = pm2
    pmk = pm1
    pm1 = pm2
  end do

  pd(0) = ( ( 1.0e+00_real64 - m ) * pm(1) - x * pm(0) ) &
    / ( x * x - 1.0e+00_real64 )  
  do k = 1, n
    pd(k) = ( k * x * pm(k) - ( k + m ) * pm(k-1) ) &
      / ( x * x - 1.0e+00_real64 )
  end do

  return
end subroutine lpmns
!> @brief subroutine lpmv.
!> @return None.
!>
!> @param v [in] Argument v.
!> @param m [in] Argument m.
!> @param x [in] Argument x.
!> @param pmv [inout] Argument pmv.
subroutine lpmv ( v, m, x, pmv )

!*****************************************************************************80
!
!! LPMV computes associated Legendre functions Pmv(X) with arbitrary degree.
!
!  Discussion:
!
!    Compute the associated Legendre function Pmv(x) with an integer order 
!    and an arbitrary nonnegative degree v.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    19 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) V, the degree of Pmv(x).
!
!    Input, integer(int32) M, the order of Pmv(x).
!
!    Input, real(real64) X, the argument of Pm(x).
!
!    Output, real(real64) PMV, the value of Pm(x).
!
  implicit none

  real(real64) c0
  real(real64) el
  real(real64) eps
  integer(int32) j
  integer(int32) k
  integer(int32), intent(in) :: m
  integer(int32) nv
  real(real64) pa
  real(real64) pi
  real(real64), intent(inout) :: pmv
  real(real64) pss
  real(real64) psv
  real(real64) pv0
  real(real64) qr
  real(real64) r
  real(real64) r0
  real(real64) r1
  real(real64) r2
  real(real64) rg
  real(real64) s
  real(real64) s0
  real(real64) s1
  real(real64) s2
  real(real64), intent(in) :: v
  real(real64) v0
  real(real64) vs
  real(real64), intent(in) :: x
  real(real64) xq

  pi = 3.141592653589793e+00_real64
  el = 0.5772156649015329e+00_real64
  eps = 1.0e-14_real64
  nv = int ( v )
  v0 = v - nv

  if ( x == -1.0e+00_real64 .and. v /= nv ) then
    if ( m == 0 ) then
      pmv = -1.0e+300_real64
    else
      pmv = 1.0e+300_real64
    end if
    return
  end if

  c0 = 1.0e+00_real64

  if ( m /= 0 ) then

    rg = v * ( v + m )
    do j = 1, m - 1 
      rg = rg * ( v * v - j * j )
    end do
    xq = sqrt ( 1.0e+00_real64 - x * x )
    r0 = 1.0e+00_real64
    do j = 1, m
      r0 = 0.5e+00_real64 * r0 * xq / j
    end do
    c0 = r0 * rg

  end if

  if ( v0 == 0.0e+00_real64 ) then

    pmv = 1.0e+00_real64
    r = 1.0e+00_real64
    do k = 1, nv - m
      r = 0.5e+00_real64 * r * ( - nv + m + k - 1.0e+00_real64 ) &
        * ( nv + m + k ) / ( k * ( k + m ) ) * ( 1.0e+00_real64 + x )
      pmv = pmv + r
    end do
    pmv = ( -1.0e+00_real64 ) ** nv * c0 * pmv

  else

    if ( -0.35e+00_real64 <= x ) then

      pmv = 1.0e+00_real64
      r = 1.0e+00_real64
      do k = 1, 100
        r = 0.5e+00_real64 * r * ( - v + m + k - 1.0e+00_real64 ) &
          * ( v + m + k ) / ( k * ( m + k ) ) * ( 1.0e+00_real64 - x )
        pmv = pmv + r
        if ( 12 < k .and. abs ( r / pmv ) < eps ) then
          exit
        end if
      end do

      pmv = ( -1.0e+00_real64 ) ** m * c0 * pmv

    else

      vs = sin ( v * pi ) / pi
      pv0 = 0.0e+00_real64

      if ( m /= 0 ) then

        qr = sqrt ( ( 1.0e+00_real64 - x ) / ( 1.0e+00_real64 + x ) )
        r2 = 1.0e+00_real64
        do j = 1, m
          r2 = r2 * qr * j
        end do
        s0 = 1.0e+00_real64
        r1 = 1.0e+00_real64
        do k = 1, m - 1 
          r1 = 0.5e+00_real64 * r1 * ( - v + k - 1 ) * ( v + k ) &
            / ( k * ( k - m ) ) * ( 1.0e+00_real64 + x )
          s0 = s0 + r1
        end do
        pv0 = - vs * r2 / m * s0

      end if

      call psi ( v, psv )
      pa = 2.0e+00_real64 * ( psv + el ) + pi / tan ( pi * v ) &
        + 1.0e+00_real64 / v

      s1 = 0.0e+00_real64
      do j = 1, m
        s1 = s1 + ( j * j + v * v ) / ( j * ( j * j - v * v ) )
      end do

      pmv = pa + s1 - 1.0e+00_real64 / ( m - v ) &
        + log ( 0.5e+00_real64 * ( 1.0e+00_real64 + x ) )
      r = 1.0e+00_real64
      do k = 1, 100
        r = 0.5e+00_real64 * r * ( - v + m + k - 1.0e+00_real64 ) * ( v + m + k ) &
          / ( k * ( k + m ) ) * ( 1.0e+00_real64 + x )
        s = 0.0e+00_real64
        do j = 1, m
          s = s + ( ( k + j ) ** 2 + v * v ) &
            / ( ( k + j ) * ( ( k + j ) ** 2 - v * v ) )
        end do
        s2 = 0.0e+00_real64
        do j = 1, k
          s2 = s2 + 1.0e+00_real64 / ( j * ( j * j - v * v ) )
        end do
        pss = pa + s + 2.0e+00_real64 * v * v * s2 &
          - 1.0e+00_real64 / ( m + k - v ) &
          + log ( 0.5e+00_real64 * ( 1.0e+00_real64 + x ) )
        r2 = pss * r
        pmv = pmv + r2
        if ( abs ( r2 / pmv ) < eps ) then
          exit
        end if
      end do

      pmv = pv0 + pmv * vs * c0

    end if

  end if

  return
end subroutine lpmv
!> @brief subroutine lpn.
!> @return None.
!>
!> @param n [in] Argument n.
!> @param x [in] Argument x.
!> @param pn [inout] Argument pn.
!> @param pd [inout] Argument pd.
pure subroutine lpn ( n, x, pn, pd )

!*****************************************************************************80
!
!! LPN computes Legendre polynomials Pn(x) and derivatives Pn'(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    07 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) N, the maximum degree.
!
!    Input, real(real64) X, the argument.
!
!    Output, real(real64) PN(0:N), PD(0:N), the values and derivatives
!    of the polyomials of degrees 0 to N at X.
!
  implicit none

  integer(int32), intent(in) :: n

  integer(int32) k
  real(real64) p0
  real(real64) p1
  real(real64), intent(inout) :: pd(0:n)
  real(real64) pf
  real(real64), intent(inout) :: pn(0:n)
  real(real64), intent(in) :: x

  pn(0) = 1.0e+00_real64
  pn(1) = x
  pd(0) = 0.0e+00_real64
  pd(1) = 1.0e+00_real64
  p0 = 1.0e+00_real64
  p1 = x

  do k = 2, n

    pf = ( 2.0e+00_real64 * k - 1.0e+00_real64 ) / k * x * p1 &
      - ( k - 1.0e+00_real64 ) / k * p0
    pn(k) = pf

    if ( abs ( x ) == 1.0e+00_real64 ) then
      pd(k) = 0.5e+00_real64 * x ** ( k + 1 ) * k * ( k + 1.0e+00_real64 )
    else
      pd(k) = k * ( p1 - x * pf ) / ( 1.0e+00_real64 - x * x )
    end if

    p0 = p1
    p1 = pf

  end do

  return
end subroutine lpn
!> @brief subroutine lpni.
!> @return None.
!>
!> @param n [in] Argument n.
!> @param x [in] Argument x.
!> @param pn [inout] Argument pn.
!> @param pd [inout] Argument pd.
!> @param pl [inout] Argument pl.
pure subroutine lpni ( n, x, pn, pd, pl )

!*****************************************************************************80
!
!! LPNI computes Legendre polynomials Pn(x), derivatives, and integrals.
!
!  Discussion:
!
!    This routine computes Legendre polynomials Pn(x), Pn'(x)
!    and the integral of Pn(t) from 0 to x.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    13 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) N, the maximum degree.
!
!    Input, real(real64) X, the argument.
!
!    Output, real(real64) PN(0:N), PD(0:N), PL(0:N), the values, 
!    derivatives and integrals of the polyomials of degrees 0 to N at X.
!
  implicit none

  integer(int32), intent(in) :: n

  integer(int32) j
  integer(int32) k
  integer(int32) n1
  real(real64) p0
  real(real64) p1
  real(real64), intent(inout) :: pd(0:n)
  real(real64) pf
  real(real64), intent(inout) :: pl(0:n)
  real(real64), intent(inout) :: pn(0:n)
  real(real64) r
  real(real64), intent(in) :: x

  pn(0) = 1.0e+00_real64
  pn(1) = x
  pd(0) = 0.0e+00_real64
  pd(1) = 1.0e+00_real64
  pl(0) = x
  pl(1) = 0.5e+00_real64 * x * x
  p0 = 1.0e+00_real64
  p1 = x

  do k = 2, n

    pf = ( 2.0e+00_real64 * k - 1.0e+00_real64 ) / k * x * p1 - ( k - 1.0e+00_real64 ) / k * p0
    pn(k) = pf

    if ( abs ( x ) == 1.0e+00_real64 ) then
      pd(k) = 0.5e+00_real64 * x ** ( k + 1 ) * k * ( k + 1.0e+00_real64 )
    else
      pd(k) = k * ( p1 - x * pf ) / ( 1.0e+00_real64 - x * x )
    end if

    pl(k) = ( x * pn(k) - pn(k-1) ) / ( k + 1.0e+00_real64 )
    p0 = p1
    p1 = pf

    if ( k /= 2 * int ( k / 2 ) ) then

      r = 1.0e+00_real64 / ( k + 1.0e+00_real64 )
      n1 = ( k - 1 ) / 2
      do j = 1, n1
        r = ( 0.5e+00_real64 / j - 1.0e+00_real64 ) * r
      end do
      pl(k) = pl(k) + r

    end if

  end do

  return
end subroutine lpni
!> @brief subroutine lqmn.
!> @return None.
!>
!> @param mm [in] Argument mm.
!> @param m [in] Argument m.
!> @param n [in] Argument n.
!> @param x [in] Argument x.
!> @param qm [inout] Argument qm.
!> @param qd [inout] Argument qd.
pure subroutine lqmn ( mm, m, n, x, qm, qd )

!*****************************************************************************80
!
!! LQMN computes associated Legendre functions Qmn(x) and derivatives.
!
!  Discussion:
!
!    This routine computes the associated Legendre functions of the
!    second kind, Qmn(x) and Qmn'(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    13 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) MM, determines the leading dimension 
!    of QM and QD.
!
!    Input, integer(int32) M, the order of Qmn(x).
!
!    Input, integer(int32) N, the degree of Qmn(x).
!
!    Output, real(real64) QM(0:MM,0:N), QD(0:MM,0:N), contains the values
!    of Qmn(x) and Qmn'(x).
!
  implicit none

  integer(int32), intent(in) :: mm
  integer(int32), intent(in) :: n

  integer(int32) i
  integer(int32) j
  integer(int32) k
  integer(int32) km
  integer(int32) ls
  integer(int32), intent(in) :: m
  real(real64) q0
  real(real64) q1
  real(real64) q10
  real(real64), intent(inout) :: qd(0:mm,0:n)
  real(real64) qf
  real(real64) qf0
  real(real64) qf1
  real(real64) qf2
  real(real64), intent(inout) :: qm(0:mm,0:n)
  real(real64), intent(in) :: x
  real(real64) xq
  real(real64) xs

  if ( abs ( x ) == 1.0e+00_real64 ) then
    do i = 0, m
      do j = 0, n
        qm(i,j) = 1.0e+300_real64
        qd(i,j) = 1.0e+300_real64
      end do
    end do
    return
  end if

  if ( 1.0e+00_real64 < abs ( x ) ) then
    ls = -1
  else
    ls = 1
  end if

  xs = ls * ( 1.0e+00_real64 - x * x )
  xq = sqrt ( xs )
  q0 = 0.5e+00_real64 * log ( abs ( ( x + 1.0e+00_real64 ) / ( x - 1.0e+00_real64 ) ) )

  if ( abs ( x ) < 1.0001e+00_real64 ) then
    qm(0,0) = q0
    qm(0,1) = x * q0 - 1.0e+00_real64
    qm(1,0) = -1.0e+00_real64 / xq
    qm(1,1) = -xq * ( q0 + x / ( 1.0e+00_real64 - x * x ) )
    do i = 0, 1
      do j = 2, n
        qm(i,j) = ( ( 2.0e+00_real64 * j - 1.0e+00_real64 ) * x * qm(i,j-1) &
              - ( j + i - 1.0e+00_real64 ) * qm(i,j-2))/ ( j - i )
      end do
    end do

    do j = 0, n
      do i = 2, m
        qm(i,j) = -2.0e+00_real64 * ( i - 1.0e+00_real64 ) * x / xq * qm(i-1,j) &
          - ls * ( j + i - 1.0e+00_real64 ) * ( j - i + 2.0e+00_real64 ) * qm(i-2,j)
      end do
    end do

  else

    if ( 1.1e+00_real64 < abs ( x ) ) then
      km = 40 + m + n
    else
      km = ( 40 + m + n ) &
        * int ( -1.0e+00_real64 - 1.8e+00_real64 * log ( x - 1.0e+00_real64 ) )
    end if

    qf2 = 0.0e+00_real64
    qf1 = 1.0e+00_real64
    do k = km, 0, -1
      qf0 = ( ( 2 * k + 3.0e+00_real64 ) * x * qf1 &
        - ( k + 2.0e+00_real64 ) * qf2 ) / ( k + 1.0e+00_real64 )
      if ( k <= n ) then
        qm(0,k) = qf0
      end if
      qf2 = qf1
      qf1 = qf0
    end do

    do k = 0, n
      qm(0,k) = q0 * qm(0,k) / qf0
    end do

    qf2 = 0.0e+00_real64
    qf1 = 1.0e+00_real64
    do k = km, 0, -1
      qf0 = ( ( 2 * k + 3.0e+00_real64 ) * x * qf1 &
        - ( k + 1.0e+00_real64 ) * qf2 ) / ( k + 2.0e+00_real64 )
      if ( k <= n ) then
        qm(1,k) = qf0
      end if
      qf2 = qf1
      qf1 = qf0
    end do

    q10 = -1.0e+00_real64 / xq
    do k = 0, n
      qm(1,k) = q10 * qm(1,k) / qf0
    end do

    do j = 0, n
      q0 = qm(0,j)
      q1 = qm(1,j)
      do i = 0, m - 2
        qf = -2.0e+00_real64 * ( i + 1 ) * x / xq * q1 &
          + ( j - i ) * ( j + i + 1.0e+00_real64 ) * q0
        qm(i+2,j) = qf
        q0 = q1
        q1 = qf
      end do
    end do

  end if

  qd(0,0) = ls / xs
  do j = 1, n
    qd(0,j) = ls * j * ( qm(0,j-1) - x * qm(0,j) ) / xs
  end do

  do j = 0, n
    do i = 1, m
      qd(i,j) = ls * i * x / xs * qm(i,j) &
        + ( i + j ) * ( j - i + 1.0e+00_real64 ) / xq * qm(i-1,j)
    end do
  end do

  return
end subroutine lqmn
!> @brief subroutine lqmns.
!> @return None.
!>
!> @param m [in] Argument m.
!> @param n [in] Argument n.
!> @param x [in] Argument x.
!> @param qm [inout] Argument qm.
!> @param qd [inout] Argument qd.
pure subroutine lqmns ( m, n, x, qm, qd )

!*****************************************************************************80
!
!! LQMNS computes associated Legendre functions Qmn(x) and derivatives Qmn'(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    28 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) M, the order.
!
!    Input, integer(int32) N, the degree.
!
!    Input, real(real64) X, the argument.
!
!    Output, real(real64) QM(0:N), QD(0:N), the values of Qmn(x) 
!    and Qmn'(x).
!
  implicit none

  integer(int32), intent(in) :: n

  integer(int32) k
  integer(int32) km
  integer(int32) l
  integer(int32) ls
  integer(int32), intent(in) :: m
  real(real64) q0
  real(real64) q00
  real(real64) q01
  real(real64) q0l
  real(real64) q10
  real(real64) q11
  real(real64) q1l
  real(real64), intent(inout) :: qd(0:n)
  real(real64) qf0
  real(real64) qf1
  real(real64) qf2
  real(real64) qg0
  real(real64) qg1
  real(real64) qh0
  real(real64) qh1
  real(real64) qh2
  real(real64), intent(inout) :: qm(0:n)
  real(real64) qm0
  real(real64) qm1
  real(real64) qmk
  real(real64), intent(in) :: x
  real(real64) xq

  do k = 0, n
    qm(k) = 0.0e+00_real64
    qd(k) = 0.0e+00_real64
  end do

  if ( abs ( x ) == 1.0e+00_real64 ) then
     do k = 0, n
       qm(k) = 1.0e+300_real64
       qd(k) = 1.0e+300_real64
     end do
     return
  end if

  if ( 1.0e+00_real64 < abs ( x ) ) then
    ls = -1
  else
    ls = +1
  end if

  xq = sqrt ( ls * ( 1.0e+00_real64 - x * x ) )
  q0 = 0.5e+00_real64 * log ( abs ( ( x + 1.0e+00_real64 ) / ( x - 1.0e+00_real64 ) ) )
  q00 = q0
  q10 = -1.0e+00_real64 / xq
  q01 = x * q0 - 1.0e+00_real64
  q11 = - ls * xq * ( q0 + x / ( 1.0e+00_real64 - x * x ) )
  qf0 = q00
  qf1 = q10
  do k = 2, m
    qm0 = -2.0e+00_real64 * ( k - 1.0e+00_real64 ) / xq * x * qf1 &
      - ls * ( k - 1.0e+00_real64 ) * ( 2.0e+00_real64 - k ) * qf0
    qf0 = qf1
    qf1 = qm0
  end do

  if ( m == 0 ) then
    qm0 = q00
  else if ( m == 1 ) then
    qm0 = q10
  end if

  qm(0) = qm0

  if ( abs ( x ) < 1.0001e+00_real64 ) then

    if ( m == 0 .and. 0 < n ) then

      qf0 = q00
      qf1 = q01
      do k = 2, n
        qf2 = ( ( 2.0e+00_real64 * k - 1.0e+00_real64 ) * x * qf1 &
          - ( k - 1.0e+00_real64 ) * qf0 ) / k
        qm(k) = qf2
        qf0 = qf1
        qf1 = qf2
      end do

    end if
    qg0 = q01
    qg1 = q11
    do k = 2, m
      qm1 = - 2.0e+00_real64 * ( k - 1.0e+00_real64 ) / xq * x * qg1 &
        - ls * k * ( 3.0e+00_real64 - k ) * qg0
      qg0 = qg1
      qg1 = qm1
    end do

    if ( m == 0 ) then
      qm1 = q01
    else if ( m == 1 ) then
      qm1 = q11
    end if
    qm(1) = qm1

    if ( m == 1 .and. 1 < n ) then

      qh0 = q10
      qh1 = q11
      do k = 2, n
        qh2 = ( ( 2.0e+00_real64 * k - 1.0e+00_real64 ) * x * qh1 - k * qh0 ) &
          / ( k - 1.0e+00_real64 )
        qm(k) = qh2
        qh0 = qh1
        qh1 = qh2
      end do

    else if ( 2 <= m ) then

      qg0 = q00
      qg1 = q01
      qh0 = q10
      qh1 = q11

      do l = 2, n
        q0l = ( ( 2.0e+00_real64 * l - 1.0e+00_real64 ) * x * qg1 &
          - ( l - 1.0e+00_real64 ) * qg0 ) / l
        q1l = ( ( 2.0e+00_real64 * l - 1.0e+00_real64 ) * x * qh1 - l * qh0 ) &
          / ( l - 1.0e+00_real64 )
        qf0 = q0l
        qf1 = q1l
        do k = 2, m
          qmk = - 2.0e+00_real64 * ( k - 1.0e+00_real64 ) / xq * x * qf1 &
            - ls * ( k + l - 1.0e+00_real64 ) * ( l + 2.0e+00_real64 - k ) * qf0
          qf0 = qf1
          qf1 = qmk
        end do
        qm(l) = qmk
        qg0 = qg1
        qg1 = q0l
        qh0 = qh1
        qh1 = q1l
      end do

    end if

  else

    if ( 1.1e+00_real64 < abs ( x ) ) then
      km = 40 + m + n
    else
      km = ( 40 + m + n ) * int ( - 1.0e+00_real64 - 1.8e+00_real64 * log ( x - 1.0e+00_real64 ) )
    end if

    qf2 = 0.0e+00_real64
    qf1 = 1.0e+00_real64
    do k = km, 0, -1
      qf0 = ( ( 2.0e+00_real64 * k + 3.0e+00_real64 ) * x * qf1 &
        - ( k + 2.0e+00_real64 - m ) * qf2 ) / ( k + m + 1.0e+00_real64 )
      if ( k <= n ) then
        qm(k) = qf0
      end if
      qf2 = qf1
      qf1 = qf0
    end do

    do k = 0, n
     qm(k) = qm(k) * qm0 / qf0
    end do

  end if

  if ( abs ( x ) < 1.0e+00_real64 ) then
    do k = 0, n
      qm(k) = ( -1 ) ** m * qm(k)
    end do
  end if

  qd(0) = ( ( 1.0e+00_real64 - m ) * qm(1) - x * qm(0) )  / ( x * x - 1.0e+00_real64 )
  do k = 1, n
    qd(k) = ( k * x * qm(k) - ( k + m ) * qm(k-1) ) / ( x * x - 1.0e+00_real64 )
  end do

  return
end subroutine lqmns
!> @brief subroutine lqna.
!> @return None.
!>
!> @param n [in] Argument n.
!> @param x [in] Argument x.
!> @param qn [inout] Argument qn.
!> @param qd [inout] Argument qd.
pure subroutine lqna ( n, x, qn, qd )

!*****************************************************************************80
!
!! LQNA computes Legendre function Qn(x) and derivatives Qn'(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    19 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) N, the degree of Qn(x).
!
!    Input, real(real64) X, the argument of Qn(x).
!
!    Output, real(real64) QN(0:N), QD(0:N), the values of
!    Qn(x) and Qn'(x).
!
  implicit none

  integer(int32), intent(in) :: n

  integer(int32) k
  real(real64) q0
  real(real64) q1
  real(real64), intent(inout) :: qd(0:n)
  real(real64) qf
  real(real64), intent(inout) :: qn(0:n)
  real(real64), intent(in) :: x

  if ( abs ( x ) == 1.0e+00_real64 ) then

    do k = 0, n
      qn(k) = 1.0e+300_real64
      qd(k) = -1.0e+300_real64
    end do

  else if ( abs ( x ) < 1.0e+00_real64 ) then

    q0 = 0.5e+00_real64 * log ( ( 1.0e+00_real64 + x ) / ( 1.0e+00_real64 - x ) )
    q1 = x * q0 - 1.0e+00_real64
    qn(0) = q0
    qn(1) = q1
    qd(0) = 1.0e+00_real64 / ( 1.0e+00_real64 - x * x )
    qd(1) = qn(0) + x * qd(0)
    do k = 2, n
      qf = ( ( 2 * k - 1 ) * x * q1 - ( k - 1 ) * q0 ) / k
      qn(k) = qf
      qd(k) = ( qn(k-1) - x * qf ) * k / ( 1.0e+00_real64 - x * x )
      q0 = q1
      q1 = qf
    end do

  end if

  return
end subroutine lqna
!> @brief subroutine lqnb.
!> @return None.
!>
!> @param n [in] Argument n.
!> @param x [in] Argument x.
!> @param qn [inout] Argument qn.
!> @param qd [inout] Argument qd.
pure subroutine lqnb ( n, x, qn, qd )

!*****************************************************************************80
!
!! LQNB computes Legendre function Qn(x) and derivatives Qn'(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    19 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) N, the degree of Qn(x).
!
!    Input, real(real64) X, the argument of Qn(x).
!
!    Output, real(real64) QN(0:N), QD(0:N), the values of
!    Qn(x) and Qn'(x).
!
  implicit none

  integer(int32), intent(in) :: n

  real(real64) eps
  integer(int32) j
  integer(int32) k
  integer(int32) l
  integer(int32) nl
  real(real64) q0
  real(real64) q1
  real(real64) qc1
  real(real64) qc2
  real(real64), intent(inout) :: qd(0:n)
  real(real64) qf
  real(real64) qf0
  real(real64) qf1
  real(real64) qf2
  real(real64), intent(inout) :: qn(0:n)
  real(real64) qr
  real(real64), intent(in) :: x
  real(real64) x2

  eps = 1.0e-14_real64

  if ( abs ( x ) == 1.0e+00_real64 ) then
    do k = 0, n
      qn(k) = 1.0e+300_real64
      qd(k) = 1.0e+300_real64
    end do
    return
  end if

  if ( x <= 1.021e+00_real64 ) then

    x2 = abs ( ( 1.0e+00_real64 + x ) / ( 1.0e+00_real64 - x ) )
    q0 = 0.5e+00_real64 * log ( x2 )
    q1 = x * q0 - 1.0e+00_real64
    qn(0) = q0
    qn(1) = q1
    qd(0) = 1.0e+00_real64 / ( 1.0e+00_real64 - x * x )
    qd(1) = qn(0) + x * qd(0)
    do k = 2, n
      qf = ( ( 2.0e+00_real64 * k - 1.0e+00_real64 ) * x * q1 &
        - ( k - 1.0e+00_real64 ) * q0 ) / k
      qn(k) = qf
      qd(k) = ( qn(k-1) - x * qf ) * k / ( 1.0e+00_real64 - x * x )
      q0 = q1
      q1 = qf
    end do

  else

    qc2 = 1.0e+00_real64 / x
    do j = 1, n
      qc2 = qc2 * j / ( ( 2.0e+00_real64 * j + 1.0e+00_real64 ) * x )
      if ( j == n - 1 ) then
        qc1 = qc2
      end if
    end do

    do l = 0, 1

      nl = n + l
      qf = 1.0e+00_real64
      qr = 1.0e+00_real64
      do k = 1, 500
        qr = qr * ( 0.5e+00_real64 * nl + k - 1.0e+00_real64 ) &
          * ( 0.5e+00_real64 * ( nl - 1 ) + k ) &
          / ( ( nl + k - 0.5e+00_real64 ) * k * x * x )
        qf = qf + qr
        if ( abs ( qr / qf ) < eps ) then
          exit
        end if
      end do

      if ( l == 0 ) then
        qn(n-1) = qf * qc1
      else
        qn(n) = qf * qc2
      end if

    end do

    qf2 = qn(n)
    qf1 = qn(n-1)
    do k = n, 2, -1
      qf0 = ( ( 2.0e+00_real64 * k - 1.0e+00_real64 ) * x * qf1 - k * qf2 ) / ( k - 1.0e+00_real64 )
      qn(k-2) = qf0
      qf2 = qf1
      qf1 = qf0
    end do

    qd(0) = 1.0e+00_real64 / ( 1.0e+00_real64 - x * x )
    do k = 1, n
      qd(k) = k * ( qn(k-1) - x * qn(k) ) / ( 1.0e+00_real64 - x * x )
    end do

  end if

  return
end subroutine lqnb
!> @brief function msta1.
!> @return Function value.
!>
!> @param x [in] Argument x.
!> @param mp [in] Argument mp.
function msta1 ( x, mp )

!*****************************************************************************80
!
!! MSTA1 determines a backward recurrence starting point for Jn(x).
!
!  Discussion:
!
!    This procedure determines the starting point for backward  
!    recurrence such that the magnitude of    
!    Jn(x) at that point is about 10^(-MP).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    08 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) X, the argument.
!
!    Input, integer(int32) MP, the negative logarithm of the 
!    desired magnitude.
!
!    Output, integer(int32) MSTA1, the starting point.
!
  implicit none

  real(real64) a0
  ! real(real64) envj
  real(real64) f
  real(real64) f0
  real(real64) f1
  integer(int32) it
  integer(int32), intent(in) :: mp
  integer(int32) msta1
  integer(int32) n0
  integer(int32) n1
  integer(int32) nn
  real(real64), intent(in) :: x

  a0 = abs ( x )
  n0 = int ( 1.1e+00_real64 * a0 ) + 1
  f0 = envj ( n0, a0 ) - mp
  n1 = n0 + 5
  f1 = envj ( n1, a0 ) - mp
  do it = 1, 20       
    nn = n1 - int ( real ( n1 - n0, kind = real64 ) / ( 1.0e+00_real64 - f0 / f1 ) )               
    f = envj ( nn, a0 ) - mp
    if ( abs ( nn - n1 ) < 1 ) then
      exit
    end if
    n0 = n1
    f0 = f1
    n1 = nn
    f1 = f
  end do

  msta1 = nn

  return
end function msta1
!> @brief function msta2.
!> @return Function value.
!>
!> @param x [in] Argument x.
!> @param n [in] Argument n.
!> @param mp [in] Argument mp.
function msta2 ( x, n, mp )

!*****************************************************************************80
!
!! MSTA2 determines a backward recurrence starting point for Jn(x).
!
!  Discussion:
!
!    This procedure determines the starting point for a backward
!    recurrence such that all Jn(x) has MP significant digits.
!
!    Jianming Jin supplied a modification to this code on 12 January 2016._real64
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    14 January 2016
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) X, the argument of Jn(x).
!
!    Input, integer(int32) N, the order of Jn(x).
!
!    Input, integer(int32) MP, the number of significant digits.
!
!    Output, integer(int32) MSTA2, the starting point.
!
  implicit none

  real(real64) a0
  real(real64) ejn
  ! real(real64) envj
  real(real64) f
  real(real64) f0
  real(real64) f1
  real(real64) hmp
  integer(int32) it
  integer(int32), intent(in) :: mp
  integer(int32) msta2
  integer(int32), intent(in) :: n
  integer(int32) n0
  integer(int32) n1
  integer(int32) nn
  real(real64) obj
  real(real64), intent(in) :: x

  a0 = abs ( x )
  hmp = 0.5e+00_real64 * mp
  ejn = envj ( n, a0 )

  if ( ejn <= hmp ) then
    obj = mp
!
!  Original code:
!
!   n0 = int ( 1.1e+00_real64 * a0 )
!
!  Updated code:
!
    n0 = int ( 1.1e+00_real64 * a0 ) + 1
  else
    obj = hmp + ejn
    n0 = n
  end if

  f0 = envj ( n0, a0 ) - obj
  n1 = n0 + 5
  f1 = envj ( n1, a0 ) - obj

  do it = 1, 20
    nn = n1 - int ( real ( n1 - n0, kind = real64 ) / ( 1.0e+00_real64 - f0 / f1 ) )
    f = envj ( nn, a0 ) - obj
    if ( abs ( nn - n1 ) < 1 ) then
      exit
    end if
    n0 = n1
    f0 = f1
    n1 = nn
    f1 = f
  end do

  msta2 = nn + 10

  return
end function msta2
!> @brief subroutine mtu0.
!> @return None.
!>
!> @param kf [in] Argument kf.
!> @param m [in] Argument m.
!> @param q [in] Argument q.
!> @param x [in] Argument x.
!> @param csf [inout] Argument csf.
!> @param csd [inout] Argument csd.
subroutine mtu0 ( kf, m, q, x, csf, csd )

!*****************************************************************************80
!
!! MTU0 computes Mathieu functions CEM(x,q) and SEM(x,q) and derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    20 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) KF, the function code.
!    1 for computing cem(x,q) and cem'(x,q)
!    2 for computing sem(x,q) and sem'(x,q).
!
!    Input, integer(int32) M, the order of the Mathieu functions.
!
!    Input, real(real64) Q, the parameter of the Mathieu functions.
!
!    Input, real(real64) X, the argument of the Mathieu functions,
!    in degrees.
!
!    Output, real(real64) CSF, CSD, the values of cem(x,q) and cem'(x,q),
!    or of sem(x,q) and sem'(x,q).
!
  implicit none

  real(real64) a
  real(real64), intent(inout) :: csd
  real(real64), intent(inout) :: csf
  real(real64) eps
  real(real64) fg(251)
  integer(int32) ic
  integer(int32) k
  integer(int32) kd
  integer(int32), intent(in) :: kf
  integer(int32) km
  integer(int32), intent(in) :: m
  real(real64), intent(in) :: q
  real(real64) qm
  real(real64) rd
  real(real64), intent(in) :: x
  real(real64) xr

  eps = 1.0e-14_real64

  if ( kf == 1 ) then

    if ( m == 2 * int ( m / 2 ) ) then
      kd = 1
    else
      kd = 2
    end if

  else

    if ( m /= 2 * int ( m / 2 ) ) then
      kd = 3
    else
      kd = 4
    end if

  end if

  call cva2 ( kd, m, q, a )

  if ( q <= 1.0e+00_real64 ) then
    qm = 7.5e+00_real64 + 56.1e+00_real64 * sqrt ( q ) - 134.7e+00_real64 * q &
      + 90.7e+00_real64 * sqrt ( q ) * q
  else
    qm = 17.0e+00_real64 + 3.1e+00_real64 * sqrt ( q ) - 0.126e+00_real64 * q &
      + 0.0037e+00_real64 * sqrt ( q ) * q
  end if

  km = int ( qm + 0.5e+00_real64 * m )
  call fcoef ( kd, m, q, a, fg )
  ic = int ( m / 2 ) + 1
  rd = 1.74532925199433e-02_real64
  xr = x * rd

  csf = 0.0e+00_real64

  do k = 1, km

    if ( kd == 1 ) then
      csf = csf + fg(k) * cos ( ( 2.0e+00_real64 * k - 2.0e+00_real64 ) * xr )
    else if ( kd == 2 ) then
      csf = csf + fg(k) * cos ( ( 2.0e+00_real64 * k - 1.0e+00_real64 ) * xr )
    else if ( kd == 3 ) then
      csf = csf + fg(k) * sin ( ( 2.0e+00_real64 * k - 1.0e+00_real64 ) * xr )
    else if ( kd == 4 ) then
      csf = csf + fg(k) * sin ( 2.0e+00_real64 * k * xr )
    end if

    if ( ic <= k .and. abs ( fg(k) ) < abs ( csf ) * eps ) then
      exit
    end if

  end do

  csd = 0.0e+00_real64

  do k = 1, km

    if ( kd == 1 ) then
      csd = csd - ( 2 * k - 2 ) * fg(k) * sin ( ( 2 * k - 2 ) * xr )
    else if ( kd == 2 ) then
      csd = csd - ( 2 * k - 1 ) * fg(k) * sin ( ( 2 * k - 1 ) * xr )
    else if ( kd == 3 ) then
      csd = csd + ( 2 * k - 1 ) * fg(k) * cos ( ( 2 * k - 1 ) * xr )
    else if ( kd == 4 ) then
      csd = csd + 2.0e+00_real64 * k * fg(k) * cos ( 2 * k * xr )
    end if

    if ( ic <= k .and. abs ( fg(k) ) < abs ( csd ) * eps ) then
      exit
    end if

  end do

  return
end subroutine mtu0
!> @brief subroutine mtu12.
!> @return None.
!>
!> @param kf [in] Argument kf.
!> @param kc [in] Argument kc.
!> @param m [in] Argument m.
!> @param q [in] Argument q.
!> @param x [in] Argument x.
!> @param f1r [inout] Argument f1r.
!> @param d1r [inout] Argument d1r.
!> @param f2r [inout] Argument f2r.
!> @param d2r [inout] Argument d2r.
subroutine mtu12 ( kf, kc, m, q, x, f1r, d1r, f2r, d2r )

!*****************************************************************************80
!
!! MTU12 computes modified Mathieu functions of the first and second kind.
!
!  Discussion:
!
!    This procedure computes modified Mathieu functions of the first and
!    second kinds, Mcm(1)(2)(x,q) and Msm(1)(2)(x,q),
!    and their derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    31 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) KF, the function code.
!    1 for computing Mcm(x,q);
!    2 for computing Msm(x,q).
!
!    Input, integer(int32) KC, the function code.
!    1, for computing the first kind
!    2, for computing the second kind or Msm(2)(x,q) and Msm(2)'(x,q)
!    3, for computing both the first and second kinds.
!
!    Input, integer(int32) M, the order of the Mathieu functions.
!
!    Input, real(real64) Q, the parameter of the Mathieu functions.
!
!    Input, real(real64) X, the argument of the Mathieu functions.
!
!    Output, real(real64) F1R, D1R, F2R, D2R, the values of 
!    Mcm(1)(x,q) or Msm(1)(x,q), Derivative of Mcm(1)(x,q) or Msm(1)(x,q),
!    Mcm(2)(x,q) or Msm(2)(x,q), Derivative of Mcm(2)(x,q) or Msm(2)(x,q).
!
  implicit none

  real(real64) a
  real(real64) bj1(0:251)
  real(real64) bj2(0:251)
  real(real64) by1(0:251)
  real(real64) by2(0:251)
  real(real64) c1
  real(real64) c2
  real(real64), intent(inout) :: d1r
  real(real64), intent(inout) :: d2r
  real(real64) dj1(0:251)
  real(real64) dj2(0:251)
  real(real64) dy1(0:251)
  real(real64) dy2(0:251)
  real(real64) eps
  real(real64), intent(inout) :: f1r
  real(real64), intent(inout) :: f2r
  real(real64) fg(251)
  integer(int32) ic
  integer(int32) k
  integer(int32), intent(in) :: kc
  integer(int32) kd
  integer(int32), intent(in) :: kf
  integer(int32) km
  integer(int32), intent(in) :: m
  integer(int32) nm
  real(real64), intent(in) :: q
  real(real64) qm
  real(real64) u1
  real(real64) u2
  real(real64) w1
  real(real64) w2
  real(real64), intent(in) :: x

  eps = 1.0e-14_real64

  if ( kf == 1 ) then
    if ( m == 2 * int ( m / 2 ) ) then
      kd = 1
    else
      kd = 2
    end if
  else
    if ( m /= 2 * int ( m / 2 ) ) then
      kd = 3
    else
      kd = 4
    end if
  end if

  call cva2 ( kd, m, q, a )

  if ( q <= 1.0e+00_real64 ) then
    qm = 7.5e+00_real64 + 56.1e+00_real64 * sqrt ( q ) - 134.7e+00_real64 * q &
      + 90.7e+00_real64 * sqrt ( q ) * q
  else
    qm = 17.0e+00_real64 + 3.1e+00_real64 * sqrt ( q ) - 0.126e+00_real64 * q &
      + 0.0037e+00_real64 * sqrt ( q ) * q
  end if

  km = int ( qm + 0.5e+00_real64 * m )              
  call fcoef ( kd, m, q, a, fg )

  if ( kd == 4 ) then
    ic = m / 2
  else
    ic = int ( m / 2 ) + 1
  end if

  c1 = exp ( - x )
  c2 = exp ( x )
  u1 = sqrt ( q ) * c1
  u2 = sqrt ( q ) * c2

  call jynb ( km, u1, nm, bj1, dj1, by1, dy1 )
  call jynb ( km, u2, nm, bj2, dj2, by2, dy2 )

  if ( kc == 1 ) then

    f1r = 0.0e+00_real64

    do k = 1, km

      if ( kd == 1 ) then
        f1r = f1r + ( - 1.0e+00_real64 ) ** ( ic + k ) * fg(k) * bj1(k-1) * bj2(k-1)
      else if ( kd == 2 .or. kd == 3 ) then
        f1r = f1r + ( -1.0e+00_real64 ) ** ( ic + k ) * fg(k) * ( bj1(k-1) * bj2(k) &
          + ( - 1.0e+00_real64 ) ** kd * bj1(k) * bj2(k-1) )
      else
        f1r = f1r + ( -1.0e+00_real64 ) ** ( ic + k ) * fg(k) &
          * ( bj1(k-1) * bj2(k+1) - bj1(k+1) * bj2(k-1) )
      end if

      if ( 5 <= k .and. abs ( f1r - w1 ) < abs ( f1r ) * eps ) then
        exit
      end if

      w1 = f1r

    end do

    f1r = f1r / fg(1)
    d1r = 0.0e+00_real64
    do k = 1, km
      if ( kd == 1 ) then
        d1r = d1r + ( - 1.0e+00_real64 ) ** ( ic + k ) * fg(k) &
          * ( c2 * bj1(k-1) * dj2(k-1) - c1 * dj1(k-1) * bj2(k-1) )
      else if ( kd == 2 .or. kd == 3 ) then
        d1r = d1r + ( -1.0e+00_real64 ) ** ( ic + k ) * fg(k) &
          * ( c2 * ( bj1(k-1) * dj2(k) &
          + ( -1.0e+00_real64 ) ** kd * bj1(k) * dj2(k-1) ) &
          - c1 * ( dj1(k-1) * bj2(k) &
          + ( -1.0e+00_real64 ) ** kd * dj1(k) * bj2(k-1) ) )
      else
        d1r = d1r + ( -1.0e+00_real64 ) ** ( ic + k ) * fg(k) &
          * ( c2 * ( bj1(k-1) * dj2(k+1) - bj1(k+1) * dj2(k-1) ) &
          - c1 * ( dj1(k-1) * bj2(k+1) - dj1(k+1) * bj2(k-1) ) )
      end if
      if ( 5 <= k .and. abs ( d1r - w2 ) < abs ( d1r ) * eps ) then
        exit
      end if
      w2 = d1r
    end do

    d1r = d1r * sqrt ( q ) / fg(1)

  else

    f2r = 0.0e+00_real64

    do k = 1, km
      if ( kd == 1 ) then
        f2r = f2r + ( -1.0e+00_real64 ) ** ( ic + k ) * fg(k) &
          * bj1(k-1) * by2(k-1)
      else if ( kd == 2 .or. kd == 3 ) then
        f2r = f2r + ( -1.0e+00_real64 ) ** ( ic + k ) * fg(k) * ( bj1(k-1) * by2(k) &
          + ( -1.0e+00_real64 ) ** kd * bj1(k) * by2(k-1) )
      else
        f2r = f2r + ( -1.0e+00_real64 ) ** ( ic + k ) * fg(k) &
          * ( bj1(k-1) * by2(k+1) - bj1(k+1) * by2(k-1) )
      end if
      if ( 5 <= k .and. abs ( f2r - w1 ) < abs ( f2r ) * eps ) then
        exit
      end if
      w1 = f2r
    end do

    f2r = f2r / fg(1)
    d2r = 0.0e+00_real64

    do k = 1, km
      if ( kd == 1 ) then
        d2r = d2r + ( -1.0e+00_real64 ) ** ( ic + k ) * fg(k) &
          * ( c2 * bj1(k-1) * dy2(k-1) - c1 * dj1(k-1) * by2(k-1) )
      else if ( kd == 2 .or. kd == 3 ) then
        d2r = d2r + ( -1.0e+00_real64 ) ** ( ic + k ) * fg(k) &
          * ( c2 * ( bj1(k-1) * dy2(k) &
          + ( -1.0e+00_real64 ) ** kd * bj1(k) * dy2(k-1) ) &
          - c1 * ( dj1(k-1) * by2(k) + ( -1.0e+00_real64 ) ** kd &
          * dj1(k) * by2(k-1) ) )
      else
        d2r = d2r + ( -1.0e+00_real64 ) ** ( ic + k ) * fg(k) &
          * ( c2 * ( bj1(k-1) * dy2(k+1) - bj1(k+1) * dy2(k-1) ) &
          - c1 * ( dj1(k-1) * by2(k+1) - dj1(k+1) * by2(k-1) ) )
      end if

      if ( 5 <= k .and. abs ( d2r - w2 ) < abs ( d2r ) * eps ) then
        exit
      end if

      w2 = d2r

    end do

    d2r = d2r * sqrt ( q ) / fg(1)

  end if

  return
end subroutine mtu12
!> @brief subroutine othpl.
!> @return None.
!>
!> @param kf [in] Argument kf.
!> @param n [in] Argument n.
!> @param x [in] Argument x.
!> @param pl [inout] Argument pl.
!> @param dpl [inout] Argument dpl.
pure subroutine othpl ( kf, n, x, pl, dpl )

!*****************************************************************************80
!
!! OTHPL computes orthogonal polynomials Tn(x), Un(x), Ln(x) or Hn(x).
!
!  Discussion:
!
!    This procedure computes orthogonal polynomials: Tn(x) or Un(x),
!    or Ln(x) or Hn(x), and their derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    08 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) KT, the function code:
!    1 for Chebyshev polynomial Tn(x)
!    2 for Chebyshev polynomial Un(x)
!    3 for Laguerre polynomial Ln(x)
!    4 for Hermite polynomial Hn(x)
!
!    Input, integer(int32) N, the order.
!
!    Input, real(real64) X, the argument.
!
!    Output, real(real64) PL(0:N), DPL(0:N), the value and derivative of
!    the polynomials of order 0 through N at X.
!
  implicit none

  integer, intent(in) :: n

  real(real64) a
  real(real64) b
  real(real64) c
  real(real64), intent(inout) :: dpl(0:n)
  real(real64) dy0
  real(real64) dy1
  real(real64) dyn
  integer(int32) k
  integer(int32), intent(in) :: kf
  real(real64), intent(inout) :: pl(0:n)
  real(real64), intent(in) :: x
  real(real64) y0
  real(real64) y1
  real(real64) yn

  a = 2.0e+00_real64
  b = 0.0e+00_real64
  c = 1.0e+00_real64
  y0 = 1.0e+00_real64
  y1 = 2.0e+00_real64 * x
  dy0 = 0.0e+00_real64
  dy1 = 2.0e+00_real64
  pl(0) = 1.0e+00_real64
  pl(1) = 2.0e+00_real64 * x
  dpl(0) = 0.0e+00_real64
  dpl(1) = 2.0e+00_real64

  if ( kf == 1 ) then
    y1 = x
    dy1 = 1.0e+00_real64
    pl(1) = x
    dpl(1) = 1.0e+00_real64
  else if ( kf == 3 ) then
    y1 = 1.0e+00_real64 - x
    dy1 = -1.0e+00_real64
    pl(1) = 1.0e+00_real64 - x
    dpl(1) = -1.0e+00_real64
  end if

  do k = 2, n

    if ( kf == 3 ) then
      a = -1.0e+00_real64 / k
      b = 2.0e+00_real64 + a
      c = 1.0e+00_real64 + a
    else if ( kf == 4 ) then
      c = 2.0e+00_real64 * ( k - 1.0e+00_real64 )
    end if

    yn = ( a * x + b ) * y1 - c * y0
    dyn = a * y1 + ( a * x + b ) * dy1 - c * dy0
    pl(k) = yn
    dpl(k) = dyn
    y0 = y1
    y1 = yn
    dy0 = dy1
    dy1 = dyn

  end do

  return
end subroutine othpl
!> @brief subroutine pbdv.
!> @return None.
!>
!> @param v [inout] Argument v.
!> @param x [in] Argument x.
!> @param dv [inout] Argument dv.
!> @param dp [inout] Argument dp.
!> @param pdf [inout] Argument pdf.
!> @param pdd [inout] Argument pdd.
subroutine pbdv ( v, x, dv, dp, pdf, pdd )

!*****************************************************************************80
!
!! PBDV computes parabolic cylinder functions Dv(x) and derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    29 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) V, the order.
!
!    Input, real(real64) X, the argument.
!
!    Output, real(real64) DV(0:*), DP(0:*), the values of
!    Dn+v0(x), Dn+v0'(x).
!
!    Output, real(real64) PDF, PDD, the values of Dv(x) and Dv'(x).
!
  implicit none

  real(real64), intent(inout) :: dp(0:)
  real(real64), intent(inout) :: dv(0:)
  real(real64) ep
  real(real64) f
  real(real64) f0
  real(real64) f1
  integer(int32) ja
  integer(int32) k
  integer(int32) l
  integer(int32) m
  integer(int32) na
  integer(int32) nk
  integer(int32) nv
  real(real64) pd
  real(real64) pd0
  real(real64) pd1
  real(real64), intent(inout) :: pdd
  real(real64), intent(inout) :: pdf
  real(real64) s0
  real(real64), intent(inout) :: v
  real(real64) v0
  real(real64) v1
  real(real64) v2
  real(real64) vh
  real(real64), intent(in) :: x
  real(real64) xa

  xa = abs ( x )
  vh = v
  v = v + sign ( 1.0e+00_real64, v )
  nv = int ( v )
  v0 = v - nv
  na = abs ( nv )
  ep = exp ( -0.25e+00_real64 * x * x )

  if ( 1 <= na ) then
    ja = 1
  end if

  if ( 0.0e+00_real64 <= v ) then
    if ( v0 == 0.0e+00_real64 ) then
      pd0 = ep
      pd1 = x * ep
    else
      do l = 0, ja
        v1 = v0 + l
        if ( xa <= 5.8e+00_real64 ) then
          call dvsa ( v1, x, pd1 )
        else
          call dvla ( v1, x, pd1 )
        end if
        if ( l == 0 ) then
          pd0 = pd1
        end if
      end do
    end if

    dv(0) = pd0
    dv(1) = pd1
    do k = 2, na
      pdf = x * pd1 - ( k + v0 - 1.0e+00_real64 ) * pd0
      dv(k) = pdf
      pd0 = pd1
      pd1 = pdf
    end do

  else

    if ( x <= 0.0e+00_real64 ) then

      if ( xa <= 5.8e+00_real64 )  then
        call dvsa ( v0, x, pd0 )
        v1 = v0 - 1.0e+00_real64
        call dvsa ( v1, x, pd1 )
      else
        call dvla ( v0, x, pd0 )
        v1 = v0 - 1.0e+00_real64
        call dvla ( v1, x, pd1 )
      end if

      dv(0) = pd0
      dv(1) = pd1
      do k = 2, na
        pd = ( - x * pd1 + pd0 ) / ( k - 1.0e+00_real64 - v0 )
        dv(k) = pd
        pd0 = pd1
        pd1 = pd
      end do

    else if ( x <= 2.0e+00_real64 ) then

      v2 = nv + v0
      if ( nv == 0 ) then
        v2 = v2 - 1.0e+00_real64
      end if

      nk = int ( - v2 )
      call dvsa ( v2, x, f1 )
      v1 = v2 + 1.0e+00_real64
      call dvsa ( v1, x, f0 )
      dv(nk) = f1
      dv(nk-1) = f0
      do k = nk - 2, 0, -1
        f = x * f0 + ( k - v0 + 1.0e+00_real64 ) * f1
        dv(k) = f
        f1 = f0
        f0 = f
      end do

    else

      if ( xa <= 5.8e+00_real64 ) then
        call dvsa ( v0, x, pd0 )
      else
        call dvla ( v0, x, pd0 )
      end if

      dv(0) = pd0
      m = 100 + na
      f1 = 0.0e+00_real64
      f0 = 1.0e-30_real64
      do k = m, 0, -1
        f = x * f0 + ( k - v0 + 1.0e+00_real64 ) * f1
        if ( k <= na ) then
          dv(k) = f
        end if
        f1 = f0
        f0 = f
      end do
      s0 = pd0 / f
      do k = 0, na
        dv(k) = s0 * dv(k)
      end do

    end if

  end if

  do k = 0, na - 1
    v1 = abs ( v0 ) + k
    if ( 0.0e+00_real64 <= v ) then
      dp(k) = 0.5e+00_real64 * x * dv(k) - dv(k+1)
    else
      dp(k) = -0.5e+00_real64 * x * dv(k) - v1 * dv(k+1)
    end if
  end do

  pdf = dv(na-1)
  pdd = dp(na-1)
  v = vh

  return
end subroutine pbdv
!> @brief subroutine pbvv.
!> @return None.
!>
!> @param v [inout] Argument v.
!> @param x [in] Argument x.
!> @param vv [inout] Argument vv.
!> @param vp [inout] Argument vp.
!> @param pvf [inout] Argument pvf.
!> @param pvd [inout] Argument pvd.
subroutine pbvv ( v, x, vv, vp, pvf, pvd )

!*****************************************************************************80
!
!! PBVV computes parabolic cylinder functions Vv(x) and their derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    29 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) V, the order.
!
!    Input, real(real64) X, the argument.
!
!    Output, real(real64) VV(0:*), VP(0:*), the values of Vv(x), Vv'(x).
!
!    Output, real(real64) PVF, PVD, the values of Vv(x) and Vv'(x).
!
  implicit none

  real(real64) f
  real(real64) f0
  real(real64) f1
  integer(int32) ja
  integer(int32) k
  integer(int32) kv
  integer(int32) l
  integer(int32) m
  integer(int32) na
  integer(int32) nv
  real(real64) pi
  real(real64) pv0
  real(real64), intent(inout) :: pvd
  real(real64), intent(inout) :: pvf
  real(real64) q2p
  real(real64) qe
  real(real64) s0
  real(real64), intent(inout) :: v
  real(real64) v0
  real(real64) v1
  real(real64) v2
  real(real64) vh
  real(real64), intent(inout) :: vp(0:)
  real(real64), intent(inout) :: vv(0:)
  real(real64), intent(in) :: x
  real(real64) xa

  pi = 3.141592653589793e+00_real64
  xa = abs ( x )
  vh = v
  v = v + sign ( 1.0e+00_real64, v )
  nv = int ( v )
  v0 = v - nv
  na = abs ( nv )
  qe = exp ( 0.25e+00_real64 * x * x )
  q2p = sqrt ( 2.0e+00_real64 / pi )

  if ( 1 <= na ) then
    ja = 1
  end if

  if ( v <= 0.0e+00_real64 ) then

    if ( v0 == 0.0e+00_real64 ) then

      if ( xa <= 7.5e+00_real64 ) then 
        call vvsa ( v0, x, pv0 )
      else
        call vvla ( v0, x, pv0 )
      end if

      f0 = q2p * qe
      f1 = x * f0
      vv(0) = pv0
      vv(1) = f0
      vv(2) = f1

    else

      do l = 0, ja
        v1 = v0 - l
        if ( xa <= 7.5e+00_real64 ) then
          call vvsa ( v1, x, f1 )
        else
          call vvla ( v1, x, f1 )
        end if
        if ( l == 0 ) then
          f0 = f1
        end if
      end do

      vv(0) = f0
      vv(1) = f1

    end if

    if ( v0 == 0.0e+00_real64 ) then
      kv = 3
    else
      kv = 2
    end if

    do k = kv, na
      f = x * f1 + ( k - v0 - 2.0e+00_real64 ) * f0
      vv(k) = f
      f0 = f1
      f1 = f
    end do

  else

    if ( 0.0e+00_real64 <= x .and. x <= 7.5e+00_real64 ) then

      v2 = v
      if ( v2 < 1.0e+00_real64 ) then
        v2 = v2 + 1.0e+00_real64
      end if

      call vvsa ( v2, x, f1 )
      v1 = v2 - 1.0e+00_real64
      kv = int ( v2 )
      call vvsa ( v1, x, f0 )
      vv(kv) = f1
      vv(kv-1) = f0
      do k = kv - 2, 0, - 1
        f = x * f0 - ( k + v0 + 2.0e+00_real64 ) * f1
        if ( k <= na ) then
          vv(k) = f
        end if
        f1 = f0
        f0 = f
      end do

    else if ( 7.5e+00_real64 < x ) then

      call vvla ( v0, x, pv0 )
      m = 100 + abs ( na )
      vv(1) = pv0
      f1 = 0.0e+00_real64
      f0 = 1.0e-40_real64
      do k = m, 0, -1
        f = x * f0 - ( k + v0 + 2.0e+00_real64 ) * f1
        if ( k <= na ) then
          vv(k) = f
        end if
        f1 = f0
        f0 = f
      end do
      s0 = pv0 / f
      do k = 0, na
        vv(k) = s0 * vv(k)
      end do

    else

      if ( xa <= 7.5e+00_real64 ) then
        call vvsa ( v0, x, f0 )
        v1 = v0 + 1.0e+00_real64
        call vvsa ( v1, x, f1 )
      else
        call vvla ( v0, x, f0 )
        v1 = v0 + 1.0e+00_real64
        call vvla ( v1, x, f1 )
      end if

      vv(0) = f0
      vv(1) = f1
      do k = 2, na
        f = ( x * f1 - f0 ) / ( k + v0 )
        vv(k) = f
        f0 = f1
        f1 = f
      end do

    end if

  end if

  do k = 0, na - 1
    v1 = v0 + k
    if ( 0.0e+00_real64 <= v ) then
      vp(k) = 0.5e+00_real64 * x * vv(k) - ( v1 + 1.0e+00_real64 ) * vv(k+1)
    else
      vp(k) = - 0.5e+00_real64 * x * vv(k) + vv(k+1)
    end if
  end do

  pvf = vv(na-1)
  pvd = vp(na-1)
  v = vh

  return
end subroutine pbvv
!> @brief subroutine pbwa.
!> @return None.
!>
!> @param a [in] Argument a.
!> @param x [in] Argument x.
!> @param w1f [inout] Argument w1f.
!> @param w1d [inout] Argument w1d.
!> @param w2f [inout] Argument w2f.
!> @param w2d [inout] Argument w2d.
subroutine pbwa ( a, x, w1f, w1d, w2f, w2d )

!*****************************************************************************80
!
!! PBWA computes parabolic cylinder functions W(a,x) and derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    29 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) A, the parameter.
!
!    Input, real(real64) X, the argument.
!
!    Output, real(real64) W1F, W1D, W2F, W2D, the values of
!    W(a,x), W'(a,x), W(a,-x), W'(a,-x).
!
  implicit none

  real(real64), intent(in) :: a
  integer(int32), parameter :: iter = 100
  real(real64) d(iter)
  real(real64) d1
  real(real64) d2
  real(real64) dl
  real(real64) eps
  real(real64) f1
  real(real64) f2
  real(real64) g1
  real(real64) g2
  real(real64) h(iter)
  real(real64) h0
  real(real64) h1
  real(real64) hl
  integer(int32) k
  integer(int32) l1
  integer(int32) l2
  integer(int32) m
  real(real64) p0
  real(real64) r
  real(real64) r1
  real(real64) ugi
  real(real64) ugr
  real(real64) vgi
  real(real64) vgr
  real(real64), intent(inout) :: w1d
  real(real64), intent(inout) :: w1f
  real(real64), intent(inout) :: w2d
  real(real64), intent(inout) :: w2f
  real(real64), intent(in) :: x
  real(real64) x1
  real(real64) x2
  real(real64) y1
  real(real64) y1d
  real(real64) y1f
  real(real64) y2d
  real(real64) y2f

  eps = 1.0e-15_real64
  p0 = 0.59460355750136e+00_real64

  if ( a == 0.0e+00_real64 ) then
    g1 = 3.625609908222e+00_real64
    g2 = 1.225416702465e+00_real64
  else
    x1 = 0.25e+00_real64
    y1 = 0.5e+00_real64 * a
    call cgama ( x1, y1, 1, ugr, ugi )
    g1 = sqrt ( ugr * ugr + ugi * ugi )
    x2 = 0.75e+00_real64
    call cgama ( x2, y1, 1, vgr, vgi )
    g2 = sqrt ( vgr * vgr + vgi * vgi )
  end if

  f1 = sqrt ( g1 / g2 )
  f2 = sqrt ( 2.0e+00_real64 * g2 / g1 )
  h0 = 1.0e+00_real64
  h1 = a
  h(1) = a
  do l1 = 4, 2*iter, 2
    m = l1 / 2
    hl = a * h1 - 0.25e+00_real64 * ( l1 - 2.0e+00_real64 ) * ( l1 - 3.0e+00_real64 ) * h0
    h(m) = hl
    h0 = h1
    h1 = hl
  end do
  y1f = 1.0e+00_real64
  r = 1.0e+00_real64
  do k = 1, iter
    r = 0.5e+00_real64 * r * x * x / ( k * ( 2.0e+00_real64 * k - 1.0e+00_real64 ) )
    r1 = h(k) * r
    y1f = y1f + r1
    if ( abs ( r1 / y1f ) <= eps .and. 30 < k ) then
      exit
    end if
  end do

  y1d = a
  r = 1.0e+00_real64
  do k = 1, iter - 1
    r = 0.5e+00_real64 * r * x * x / ( k * ( 2.0e+00_real64 * k + 1.0e+00_real64 ) )
    r1 = h(k+1) * r
    y1d = y1d + r1
    if ( abs ( r1 / y1d ) <= eps .and. 30 < k ) then
      exit
    end if
  end do

  y1d = x * y1d
  d1 = 1.0e+00_real64
  d2 = a
  d(1) = 1.0e+00_real64
  d(2) = a
  do l2 = 5, 160, 2
    m = ( l2 + 1 ) / 2
    dl = a * d2 - 0.25e+00_real64 * ( l2 - 2.0e+00_real64 ) * ( l2 - 3.0e+00_real64 ) * d1
    d(m) = dl
    d1 = d2
    d2 = dl
  end do

  y2f = 1.0e+00_real64
  r = 1.0e+00_real64
  do k = 1, iter - 1
    r = 0.5e+00_real64 * r * x * x / ( k * ( 2.0e+00_real64 * k + 1.0e+00_real64 ) )
    r1 = d(k+1) * r
    y2f = y2f + r1
    if ( abs ( r1 / y2f ) <= eps .and. 30 < k ) then
      exit
    end if
  end do

  y2f = x * y2f
  y2d = 1.0e+00_real64
  r = 1.0e+00_real64
  do k = 1, iter - 1
    r = 0.5e+00_real64 * r * x * x / ( k * ( 2.0e+00_real64 * k - 1.0e+00_real64 ) )
    r1 = d(k+1) * r
    y2d = y2d + r1
    if ( abs ( r1 / y2d ) <= eps .and. 30 < k ) then
      exit
    end if
  end do

  w1f = p0 * ( f1 * y1f - f2 * y2f )
  w2f = p0 * ( f1 * y1f + f2 * y2f )
  w1d = p0 * ( f1 * y1d - f2 * y2d )
  w2d = p0 * ( f1 * y1d + f2 * y2d )

  return
end subroutine pbwa
!> @brief subroutine psi.
!> @return None.
!>
!> @param x [in] Argument x.
!> @param ps [inout] Argument ps.
pure subroutine psi ( x, ps )

!*****************************************************************************80
!
!! PSI computes the PSI function.
!
!  Licensing:
!
!    The original FORTRAN77 version of this routine is copyrighted by 
!    Shanjie Zhang and Jianming Jin.  However, they give permission to 
!    incorporate this routine into a user program that the copyright 
!    is acknowledged.
!
!  Modified:
!
!    08 September 2007
!
!  Author:
!
!    Original FORTRAN77 by Shanjie Zhang, Jianming Jin.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45
!
!  Parameters:
!
!    Input, real(real64) X, the argument.
!
!    Output, real(real64) PS, the value of the PSI function.
!
  implicit none

  real(real64), parameter :: a1 = -0.83333333333333333e-01_real64
  real(real64), parameter :: a2 =  0.83333333333333333e-02_real64
  real(real64), parameter :: a3 = -0.39682539682539683e-02_real64
  real(real64), parameter :: a4 =  0.41666666666666667e-02_real64
  real(real64), parameter :: a5 = -0.75757575757575758e-02_real64
  real(real64), parameter :: a6 =  0.21092796092796093e-01_real64
  real(real64), parameter :: a7 = -0.83333333333333333e-01_real64
  real(real64), parameter :: a8 =  0.4432598039215686e+00_real64
  real(real64), parameter :: el = 0.5772156649015329e+00_real64
  integer(int32) k
  integer(int32) n
  real(real64), parameter :: pi = 3.141592653589793e+00_real64
  real(real64), intent(inout) :: ps
  real(real64) s
  real(real64), intent(in) :: x
  real(real64) x2
  real(real64) xa

  xa = abs ( x )
  s = 0.0e+00_real64

  if ( x == aint ( x ) .and. x <= 0.0e+00_real64 ) then

    ps = 1.0e+300_real64
    return

  else if ( xa == aint ( xa ) ) then

    n = int ( xa )
    do k = 1, n - 1
      s = s + 1.0e+00_real64 / real ( k, kind = real64 )
    end do

    ps = - el + s

  else if ( xa + 0.5e+00_real64 == aint ( xa + 0.5e+00_real64 ) ) then

    n = int ( xa - 0.5e+00_real64 )

    do k = 1, n
      s = s + 1.0e+00_real64 / real ( 2 * k - 1, kind = real64 )
    end do

    ps = - el + 2.0e+00_real64 * s - 1.386294361119891e+00_real64

  else

    if ( xa < 10.0e+00_real64 ) then

      n = 10 - int ( xa )
      do k = 0, n - 1
        s = s + 1.0e+00_real64 / ( xa + real ( k, kind = real64 ) )
      end do

      xa = xa + real ( n, kind = real64 )

    end if

    x2 = 1.0e+00_real64 / ( xa * xa )

    ps = log ( xa ) - 0.5e+00_real64 / xa + x2 * ((((((( &
             a8   &
      * x2 + a7 ) &
      * x2 + a6 ) &
      * x2 + a5 ) &
      * x2 + a4 ) &
      * x2 + a3 ) &
      * x2 + a2 ) &
      * x2 + a1 )

    ps = ps - s

  end if

  if ( x < 0.0e+00_real64 ) then
    ps = ps - pi * cos ( pi * x ) / sin ( pi * x ) - 1.0e+00_real64 / x
  end if

  return
end subroutine psi
!> @brief subroutine qstar.
!> @return None.
!>
!> @param m [in] Argument m.
!> @param n [in] Argument n.
!> @param c [in] Argument c.
!> @param ck [in] Argument ck.
!> @param ck1 [in] Argument ck1.
!> @param qs [inout] Argument qs.
!> @param qt [inout] Argument qt.
pure subroutine qstar ( m, n, c, ck, ck1, qs, qt )

!*****************************************************************************80
!
!! QSTAR computes Q*mn(-ic) for oblate radial functions with a small argument.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    18 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) M, the mode parameter;  M = 0, 1, 2, ...
!
!    Input, integer(int32) N, mode parameter, N = M, M + 1, M + 2, ...
!
!    Input, real(real64) C, spheroidal parameter.
!
!    Input, real(real64) CK(*), ?
!
!    Input, real(real64) CK1, ?
!
!    Output, real(real64) QS, ?
!
!    Output, real(real64) QT, ?
!
  implicit none

  real(real64) ap(200)
  real(real64), intent(in) :: c
  real(real64), intent(in) :: ck(200)
  real(real64), intent(in) :: ck1
  integer(int32) i
  integer(int32) ip
  integer(int32) k
  integer(int32) l
  integer(int32), intent(in) :: m
  integer(int32), intent(in) :: n
  real(real64), intent(inout) :: qs
  real(real64) qs0
  real(real64), intent(inout) :: qt
  real(real64) r
  real(real64) s
  real(real64) sk

  if ( n - m == 2 * int ( ( n - m ) / 2 ) ) then
    ip = 0
  else
    ip = 1
  end if

  r = 1.0e+00_real64 / ck(1) ** 2
  ap(1) = r
  do i = 1, m
    s = 0.0e+00_real64
    do l = 1, i
      sk = 0.0e+00_real64
      do k = 0, l
        sk = sk + ck(k+1) * ck(l-k+1)
      end do
      s = s + sk * ap(i-l+1)
    end do
    ap(i+1) = -r * s
  end do 

  qs0 = ap(m+1)     
  do l = 1, m
    r = 1.0e+00_real64
    do k = 1, l
      r = r * ( 2.0e+00_real64 * k + ip ) &
        * ( 2.0e+00_real64 * k - 1.0e+00_real64 + ip ) / ( 2.0e+00_real64 * k ) ** 2
    end do
    qs0 = qs0 + ap(m-l+1) * r
  end do

  qs = ( -1.0e+00_real64 ) ** ip * ck1 * ( ck1 * qs0 ) / c
  qt = - 2.0e+00_real64 / ck1 * qs

  return
end subroutine qstar
!> @brief function r8_gamma_log.
!> @return Function value.
!>
!> @param x [in] Argument x.
function r8_gamma_log ( x )

!*****************************************************************************80
!
!! R8_GAMMA_LOG evaluates the logarithm of the gamma function.
!
!  Discussion:
!
!    This routine calculates the LOG(GAMMA) function for a positive real
!    argument X.  Computation is based on an algorithm outlined in
!    references 1 and 2._real64  The program uses rational functions that
!    theoretically approximate LOG(GAMMA) to at least 18 significant
!    decimal digits.  The approximation for X > 12 is from reference
!    3, while approximations for X < 12.0_real64 are similar to those in
!    reference 1, but are unpublished.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 2013
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody, Kenneth Hillstrom,
!    Chebyshev Approximations for the Natural Logarithm of the
!    Gamma Function,
!    Mathematics of Computation,
!    Volume 21, Number 98, April 1967, pages 198-203._real64
!
!    Kenneth Hillstrom,
!    ANL/AMD Program ANLC366S, DGAMMA/DLGAMA,
!    May 1969._real64
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
!    Charles Mesztenyi, John Rice, Henry Thatcher,
!    Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968,
!    LC: QA297.C64.
!
!  Parameters:
!
!    Input, real(real64) X, the argument of the function.
!
!    Output, real(real64) R8_GAMMA_LOG, the value of the function.
!
  implicit none

  real(real64), dimension ( 7 ) :: c = [&
    -1.910444077728e-03_real64, &
     8.4171387781295e-04_real64, &
    -5.952379913043012e-04_real64, &
     7.93650793500350248e-04_real64, &
    -2.777777777777681622553e-03_real64, &
     8.333333333333333331554247e-02_real64, &
     5.7083835261e-03_real64]
  real(real64) corr
  real(real64) :: d1 = -5.772156649015328605195174e-01_real64
  real(real64) :: d2 = 4.227843350984671393993777e-01_real64
  real(real64) :: d4 = 1.791759469228055000094023e+00_real64
  real(real64), parameter :: frtbig = 2.25e+76_real64
  integer(int32) i
  real(real64), dimension ( 8 ) :: p1 = [&
    4.945235359296727046734888e+00_real64, &
    2.018112620856775083915565e+02_real64, &
    2.290838373831346393026739e+03_real64, &
    1.131967205903380828685045e+04_real64, &
    2.855724635671635335736389e+04_real64, &
    3.848496228443793359990269e+04_real64, &
    2.637748787624195437963534e+04_real64, &
    7.225813979700288197698961e+03_real64]
  real(real64), dimension ( 8 ) :: p2 = [&
    4.974607845568932035012064e+00_real64, &
    5.424138599891070494101986e+02_real64, &
    1.550693864978364947665077e+04_real64, &
    1.847932904445632425417223e+05_real64, &
    1.088204769468828767498470e+06_real64, &
    3.338152967987029735917223e+06_real64, &
    5.106661678927352456275255e+06_real64, &
    3.074109054850539556250927e+06_real64]
  real(real64), dimension ( 8 ) :: p4 = [&
    1.474502166059939948905062e+04_real64, &
    2.426813369486704502836312e+06_real64, &
    1.214755574045093227939592e+08_real64, &
    2.663432449630976949898078e+09_real64, &
    2.940378956634553899906876e+10_real64, &
    1.702665737765398868392998e+11_real64, &
    4.926125793377430887588120e+11_real64, &
    5.606251856223951465078242e+11_real64]
  real(real64), dimension ( 8 ) :: q1 = [&
    6.748212550303777196073036e+01_real64, &
    1.113332393857199323513008e+03_real64, &
    7.738757056935398733233834e+03_real64, &
    2.763987074403340708898585e+04_real64, &
    5.499310206226157329794414e+04_real64, &
    6.161122180066002127833352e+04_real64, &
    3.635127591501940507276287e+04_real64, &
    8.785536302431013170870835e+03_real64]
  real(real64), dimension ( 8 ) :: q2 = [&
    1.830328399370592604055942e+02_real64, &
    7.765049321445005871323047e+03_real64, &
    1.331903827966074194402448e+05_real64, &
    1.136705821321969608938755e+06_real64, &
    5.267964117437946917577538e+06_real64, &
    1.346701454311101692290052e+07_real64, &
    1.782736530353274213975932e+07_real64, &
    9.533095591844353613395747e+06_real64]
  real(real64), dimension ( 8 ) :: q4 = [&
    2.690530175870899333379843e+03_real64, &
    6.393885654300092398984238e+05_real64, &
    4.135599930241388052042842e+07_real64, &
    1.120872109616147941376570e+09_real64, &
    1.488613728678813811542398e+10_real64, &
    1.016803586272438228077304e+11_real64, &
    3.417476345507377132798597e+11_real64, &
    4.463158187419713286462081e+11_real64]
  real(real64) r8_gamma_log
  real(real64) res
  real(real64), parameter :: sqrtpi = 0.9189385332046727417803297e+00_real64
  real(real64), intent(in) :: x
  real(real64), parameter :: xbig = 2.55e+305_real64
  real(real64) xden
  real(real64), parameter :: xinf = 1.79e+308_real64
  real(real64) xm1
  real(real64) xm2
  real(real64) xm4
  real(real64) xnum
  real(real64) y
  real(real64) ysq

  y = x

  if ( 0.0e+00_real64 < y .and. y <= xbig ) then

    if ( y <= epsilon ( y ) ) then

      res = - log ( y )
!
!  EPS < X <= 1.5_real64.
!
    else if ( y <= 1.5e+00_real64 ) then

      if ( y < 0.6796875e+00_real64 ) then
        corr = -log ( y )
        xm1 = y
      else
        corr = 0.0e+00_real64
        xm1 = ( y - 0.5e+00_real64 ) - 0.5e+00_real64
      end if

      if ( y <= 0.5e+00_real64 .or. 0.6796875e+00_real64 <= y ) then

        xden = 1.0e+00_real64
        xnum = 0.0e+00_real64
        do i = 1, 8
          xnum = xnum * xm1 + p1(i)
          xden = xden * xm1 + q1(i)
        end do

        res = corr + ( xm1 * ( d1 + xm1 * ( xnum / xden ) ) )

      else

        xm2 = ( y - 0.5e+00_real64 ) - 0.5e+00_real64
        xden = 1.0e+00_real64
        xnum = 0.0e+00_real64
        do i = 1, 8
          xnum = xnum * xm2 + p2(i)
          xden = xden * xm2 + q2(i)
        end do

        res = corr + xm2 * ( d2 + xm2 * ( xnum / xden ) )

      end if
!
!  1.5_real64 < X <= 4.0_real64.
!
    else if ( y <= 4.0e+00_real64 ) then

      xm2 = y - 2.0e+00_real64
      xden = 1.0e+00_real64
      xnum = 0.0e+00_real64
      do i = 1, 8
        xnum = xnum * xm2 + p2(i)
        xden = xden * xm2 + q2(i)
      end do

      res = xm2 * ( d2 + xm2 * ( xnum / xden ) )
!
!  4.0_real64 < X <= 12.0_real64.
!
    else if ( y <= 12.0e+00_real64 ) then

      xm4 = y - 4.0e+00_real64
      xden = -1.0e+00_real64
      xnum = 0.0e+00_real64
      do i = 1, 8
        xnum = xnum * xm4 + p4(i)
        xden = xden * xm4 + q4(i)
      end do

      res = d4 + xm4 * ( xnum / xden )
!
!  Evaluate for 12 <= argument.
!
    else

      res = 0.0e+00_real64

      if ( y <= frtbig ) then

        res = c(7)
        ysq = y * y

        do i = 1, 6
          res = res / ysq + c(i)
        end do

      end if

      res = res / y
      corr = log ( y )
      res = res + sqrtpi - 0.5e+00_real64 * corr
      res = res + y * ( corr - 1.0e+00_real64 )

    end if
!
!  Return for bad arguments.
!
  else

    res = xinf

  end if
!
!  Final adjustments and return.
!
  r8_gamma_log = res

  return
end function r8_gamma_log
!> @brief subroutine rctj.
!> @return None.
!>
!> @param n [in] Argument n.
!> @param x [in] Argument x.
!> @param nm [inout] Argument nm.
!> @param rj [inout] Argument rj.
!> @param dj [inout] Argument dj.
subroutine rctj ( n, x, nm, rj, dj )

!*****************************************************************************80
!
!! RCTJ computes Riccati-Bessel function of the first kind, and derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    18 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) N, the order of jn(x).
!
!    Input, real(real64) X, the argument.
!
!    Output, integer(int32) NM, the highest order computed.
!
!    Output, real(real64) RJ(0:N), the values of x jn(x).
!
!    Output, real(real64) DJ(0:N), the values of [x jn(x)]'.
!
  implicit none

  integer(int32), intent(in) :: n

  real(real64) cs
  real(real64), intent(inout) :: dj(0:n)
  real(real64) f
  real(real64) f0
  real(real64) f1
  integer(int32) k
  integer(int32) m
  ! integer(int32) msta1
  ! integer(int32) msta2
  integer(int32), intent(inout) :: nm
  real(real64), intent(inout) :: rj(0:n)
  real(real64) rj0
  real(real64) rj1
  real(real64), intent(in) :: x

  nm = n

  if ( abs ( x ) < 1.0e-100_real64 ) then
    do k = 0, n
      rj(k) = 0.0e+00_real64
      dj(k) = 0.0e+00_real64
    end do
    dj(0) = 1.0e+00_real64
    return
  end if

  rj(0) = sin ( x )
  rj(1) = rj(0) / x - cos ( x )
  rj0 = rj(0)
  rj1 = rj(1)

  if ( 2 <= n ) then

    m = msta1 ( x, 200 )

    if ( m < n ) then
      nm = m
    else
      m = msta2 ( x, n, 15 )
    end if

    f0 = 0.0e+00_real64
    f1 = 1.0e-100_real64
    do k = m, 0, -1
      f = ( 2.0e+00_real64 * k + 3.0e+00_real64 ) * f1 / x - f0
      if ( k <= nm ) then
        rj(k) = f
      end if
      f0 = f1
      f1 = f
    end do

    if ( abs ( rj1 ) < abs ( rj0 ) ) then
      cs = rj0 / f
    else
      cs = rj1 / f0
    end if

    do k = 0, nm
      rj(k) = cs * rj(k)
    end do

  end if

  dj(0) = cos ( x )
  do k = 1, nm
    dj(k) = - k * rj(k) / x + rj(k-1)
  end do

  return
end subroutine rctj
!> @brief subroutine rcty.
!> @return None.
!>
!> @param n [in] Argument n.
!> @param x [in] Argument x.
!> @param nm [inout] Argument nm.
!> @param ry [inout] Argument ry.
!> @param dy [inout] Argument dy.
pure subroutine rcty ( n, x, nm, ry, dy )

!*****************************************************************************80
!
!! RCTY computes Riccati-Bessel function of the second kind, and derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    18 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) N, the order of yn(x).
!
!    Input, real(real64) X, the argument.
!
!    Output, integer(int32) NM, the highest order computed.
!
!    Output, real(real64) RY(0:N), the values of x yn(x).
!
!    Output, real(real64) DY(0:N), the values of [x yn(x)]'.
!
  implicit none

  integer(int32), intent(in) :: n

  real(real64), intent(inout) :: dy(0:n)
  integer(int32) k
  integer(int32), intent(inout) :: nm
  real(real64) rf0
  real(real64) rf1
  real(real64) rf2
  real(real64), intent(inout) :: ry(0:n)
  real(real64), intent(in) :: x

  nm = n

  if ( x < 1.0e-60_real64 ) then
    do k = 0, n
      ry(k) = -1.0e+300_real64
      dy(k) = 1.0e+300_real64
    end do
    ry(0) = -1.0e+00_real64
    dy(0) = 0.0e+00_real64
    return
  end if

  ry(0) = - cos ( x )
  ry(1) = ry(0) / x - sin ( x )
  rf0 = ry(0)
  rf1 = ry(1)
  do k = 2, n
    rf2 = ( 2.0e+00_real64 * k - 1.0e+00_real64 ) * rf1 / x - rf0
    if ( 1.0e+300_real64 < abs ( rf2 ) ) then
      exit
    end if
    ry(k) = rf2
    rf0 = rf1
    rf1 = rf2
  end do

  nm = k - 1
  dy(0) = sin ( x )
  do k = 1, nm
    dy(k) = - k * ry(k) / x + ry(k-1)
  end do

  return
end subroutine rcty
!> @brief subroutine refine.
!> @return None.
!>
!> @param kd [in] Argument kd.
!> @param m [in] Argument m.
!> @param q [in] Argument q.
!> @param a [inout] Argument a.
!> @param iflag [inout] Argument iflag.
subroutine refine ( kd, m, q, a, iflag )

!*****************************************************************************80
!
!! REFINE refines an estimate of the characteristic value of Mathieu functions.
!
!  Discussion:
!
!    This procedure calculates the accurate characteristic value
!    by the secant method.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    20 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) KD, the case code:
!    1, for cem(x,q)  ( m = 0,2,4,...)
!    2, for cem(x,q)  ( m = 1,3,5,...)
!    3, for sem(x,q)  ( m = 1,3,5,...)
!    4, for sem(x,q)  ( m = 2,4,6,...)
!
!    Input, integer(int32) M, the order of the Mathieu functions.
!
!    Input, real(real64) Q, the parameter of the Mathieu functions.
!
!    Input/output, real(real64) A, the characteristic value, which
!    should have been refined on output.
!
  implicit none

  real(real64), intent(inout) :: a
  real(real64) ca
  real(real64) delta
  real(real64) eps
  real(real64) f
  real(real64) f0
  real(real64) f1
  integer(int32) it
  integer(int32), intent(inout) :: iflag
  integer(int32), intent(in) :: kd
  integer(int32), intent(in) :: m
  integer(int32) mj
  real(real64), intent(in) :: q
  real(real64) x
  real(real64) x0
  real(real64) x1

  eps = 1.0e-14_real64
  mj = 10 + m
  ca = a
  delta = 0.0e+00_real64
  x0 = a
  call cvf ( kd, m, q, x0, mj, f0 )
  x1 = 1.002e+00_real64 * a
  call cvf ( kd, m, q, x1, mj, f1 )

  do
  
    do it = 1, 100
      mj = mj + 1
      x = x1 - ( x1 - x0 ) / ( 1.0e+00_real64 - f0 / f1 )
      call cvf ( kd, m, q, x, mj, f )
      if ( abs ( 1.0e+00_real64 - x1 / x ) < eps .or. f == 0.0e+00_real64 ) then
        exit
      end if
      x0 = x1
      f0 = f1
      x1 = x
      f1 = f
    end do

    a = x

    if ( 0.05e+00_real64 < delta ) then
      a = ca
      if ( iflag < 0 ) then
        iflag = -10
      end if
      return
    end if

    if ( abs ( ( a - ca ) / ca )  <= 0.05e+00_real64 ) then
      exit
    end if

    x0 = ca
    delta = delta + 0.005e+00_real64
    call cvf ( kd, m, q, x0, mj, f0 )
    x1 = ( 1.0e+00_real64 + delta ) * ca
    call cvf ( kd, m, q, x1, mj, f1 )

  end do

  return
end subroutine refine
!> @brief subroutine rmn1.
!> @return None.
!>
!> @param m [in] Argument m.
!> @param n [in] Argument n.
!> @param c [in] Argument c.
!> @param x [in] Argument x.
!> @param df [in] Argument df.
!> @param kd [in] Argument kd.
!> @param r1f [inout] Argument r1f.
!> @param r1d [inout] Argument r1d.
subroutine rmn1 ( m, n, c, x, df, kd, r1f, r1d )

!*****************************************************************************80
!
!! RMN1 computes prolate and oblate spheroidal functions of the first kind.
!
!  Discussion:
!
!    This procedure computes prolate and oblate spheroidal radial
!    functions of the first kind for given m, n, c and x.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    29 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) M, the mode parameter;  M = 0, 1, 2, ...
!
!    Input, integer(int32) N, mode parameter, N = M, M + 1, M + 2, ...
!
!    Input, real(real64) C, spheroidal parameter.
!
!    Input, real(real64) X, the argument.
!
!    Input, real(real64) DF(*), the expansion coefficients.
!
!    Input, integer(int32) KD, the function code.
!    1, the prolate function.
!    -1, the oblate function.
!
!    Output, real(real64) R1F, R1D, the function and derivative.
!
  implicit none

  real(real64) a0
  real(real64) b0
  real(real64), intent(in) :: c
  real(real64) ck(200)
  real(real64) cx
  real(real64), intent(in) :: df(200)
  real(real64) dj(0:251)
  real(real64) eps
  integer(int32) ip
  integer(int32) j
  integer(int32) k
  integer(int32), intent(in) :: kd
  integer(int32) l
  integer(int32) lg
  integer(int32), intent(in) :: m
  integer(int32), intent(in) :: n
  integer(int32) nm
  integer(int32) nm1
  integer(int32) nm2
  integer(int32) np
  real(real64) r
  real(real64) r0
  real(real64) r1
  real(real64), intent(inout) :: r1d
  real(real64), intent(inout) :: r1f
  real(real64) r2
  real(real64) r3
  real(real64) reg
  real(real64) sa0
  real(real64) sj(0:251)
  real(real64) suc
  real(real64) sud
  real(real64) sum
  real(real64) sw
  real(real64) sw1
  real(real64), intent(in) :: x

  eps = 1.0e-14_real64
  nm1 = int ( ( n - m ) / 2 )
  if ( n - m == 2 * nm1 ) then
    ip = 0
  else
    ip = 1
  end if
  nm = 25 + nm1 + int ( c )
  reg = 1.0e+00_real64
  if ( 80 < m + nm ) then
    reg = 1.0e-200_real64
  end if
  r0 = reg
  do j = 1, 2 * m + ip
    r0 = r0 * j
  end do
  r = r0    
  suc = r * df(1)
  sw = 0.0e+00_real64
  do k = 2, nm
    r = r * ( m + k - 1.0e+00_real64 ) * ( m + k + ip - 1.5e+00_real64 ) &
      / ( k - 1.0e+00_real64 ) / ( k + ip - 1.5e+00_real64 )
    suc = suc + r * df(k)

    if ( nm1 < k .and. abs ( suc - sw ) < abs ( suc ) * eps ) then
      exit
    end if

    sw = suc

  end do

  if ( x == 0.0e+00_real64 ) then

    call sckb ( m, n, c, df, ck )
    sum = 0.0e+00_real64
    sw1 = 0.0e+00_real64
    do j = 1, nm
      sum = sum + ck(j)
      if ( abs ( sum - sw1 ) < abs ( sum ) * eps ) then
        exit
      end if
      sw1 = sum
    end do

    r1 = 1.0e+00_real64
    do j = 1, ( n + m + ip ) / 2
      r1 = r1 * ( j + 0.5e+00_real64 * ( n + m + ip ) )
    end do

    r2 = 1.0e+00_real64
    do j = 1, m
      r2 = 2.0e+00_real64 * c * r2 * j
    end do

    r3 = 1.0e+00_real64
    do j = 1, ( n - m - ip ) / 2
      r3 = r3 * j
    end do

    sa0 = ( 2.0e+00_real64 * ( m + ip ) + 1.0e+00_real64 ) * r1 &
      / ( 2.0e+00_real64 ** n * c ** ip * r2 * r3 )

    if ( ip == 0 ) then
      r1f = sum / ( sa0 * suc ) * df(1) * reg
      r1d = 0.0e+00_real64
    else if ( ip == 1 ) then
      r1f = 0.0e+00_real64
      r1d = sum / ( sa0 * suc ) * df(1) * reg
    end if

    return

  end if

  cx = c * x
  nm2 = 2 * nm + m
  call sphj ( nm2, cx, nm2, sj, dj )
  a0 = ( 1.0e+00_real64 - kd / ( x * x ) ) ** ( 0.5e+00_real64 * m ) / suc  
  r1f = 0.0e+00_real64
  sw = 0.0e+00_real64
  do k = 1, nm
    l = 2 * k + m - n - 2 + ip
    if ( l == 4 * int ( l / 4 ) ) then
      lg = 1
    else
      lg = -1
    end if
    if ( k == 1 ) then
      r = r0
    else
      r = r * ( m + k - 1.0e+00_real64 ) * ( m + k + ip - 1.5e+00_real64 ) &
        / ( k - 1.0e+00_real64 ) / ( k + ip - 1.5e+00_real64 )
    end if
    np = m + 2 * k - 2 + ip
    r1f = r1f + lg * r * df(k) * sj(np)
    if ( nm1 < k .and. abs ( r1f - sw ) < abs ( r1f ) * eps ) then
      exit
    end if
    sw = r1f
  end do

  r1f = r1f * a0
  b0 = kd * m / x ** 3.0e+00_real64 / ( 1.0e+00_real64 - kd / ( x * x ) ) * r1f    
  sud = 0.0e+00_real64
  sw = 0.0e+00_real64

  do k = 1, nm

    l = 2 * k + m - n - 2 + ip

    if ( l == 4 * int ( l / 4 ) ) then
      lg = 1
    else
      lg = -1
    end if

    if ( k == 1 ) then
      r = r0
    else
      r = r * ( m + k - 1.0e+00_real64 ) * ( m + k + ip - 1.5e+00_real64 ) &
        / ( k - 1.0e+00_real64 ) / ( k + ip - 1.5e+00_real64 )
    end if

    np = m + 2 * k - 2 + ip
    sud = sud + lg * r * df(k) * dj(np)
    if ( nm1 < k .and. abs ( sud - sw ) < abs ( sud ) * eps ) then
      exit
    end if
    sw = sud
  end do

  r1d = b0 + a0 * c * sud

  return
end subroutine rmn1
!> @brief subroutine rmn2l.
!> @return None.
!>
!> @param m [in] Argument m.
!> @param n [in] Argument n.
!> @param c [in] Argument c.
!> @param x [in] Argument x.
!> @param df [in] Argument df.
!> @param kd [in] Argument kd.
!> @param r2f [inout] Argument r2f.
!> @param r2d [inout] Argument r2d.
!> @param id [inout] Argument id.
subroutine rmn2l ( m, n, c, x, df, kd, r2f, r2d, id )

!*****************************************************************************80
!
!! RMN2L: prolate and oblate spheroidal functions, second kind, large CX.
!
!  Discussion:
!
!    This procedure computes prolate and oblate spheroidal radial functions 
!    of the second kind for given m, n, c and a large cx.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    30 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) M, the mode parameter;  M = 0, 1, 2, ...
!
!    Input, integer(int32) N, mode parameter, N = M, M + 1, M + 2, ...
!
!    Input, real(real64) C, spheroidal parameter.
!
!    Input, real(real64) X, the argument.
!
!    Input, real(real64) DF(*), the expansion coefficients.
!
!    Input, integer(int32) KD, the function code.
!    1, the prolate function.
!    -1, the oblate function.
!
!    Output, real(real64) R2F, R2D, the function and derivative values.
!
  implicit none

  real(real64) a0
  real(real64) b0
  real(real64), intent(in) :: c
  real(real64) cx
  real(real64), intent(in) :: df(200)
  real(real64) dy(0:251)
  real(real64) eps
  real(real64) eps1
  real(real64) eps2
  integer(int32), intent(inout) :: id
  integer(int32) id1
  integer(int32) id2
  integer(int32) ip
  integer(int32) j
  integer(int32) k
  integer(int32), intent(in) :: kd
  integer(int32) l
  integer(int32) lg
  integer(int32), intent(in) :: m
  integer(int32), intent(in) :: n
  integer(int32) nm
  integer(int32) nm1
  integer(int32) nm2
  integer(int32) np
  real(real64) r
  real(real64) r0
  real(real64), intent(inout) :: r2d
  real(real64), intent(inout) :: r2f
  real(real64) reg
  real(real64) sw
  real(real64) suc
  real(real64) sud
  real(real64) sy(0:251)
  real(real64), intent(in) :: x

  eps = 1.0e-14_real64

  nm1 = int ( ( n - m ) / 2 )

  if ( n - m == 2 * nm1 ) then
    ip = 0
  else
    ip = 1
  end if
  nm = 25 + nm1 + int ( c )

  if ( 80 < m + nm ) then
    reg = 1.0e-200_real64
  else
    reg = 1.0e+00_real64
  end if
  nm2 = 2 * nm + m
  cx = c * x
  call sphy ( nm2, cx, nm2, sy, dy )
  r0 = reg
  do j = 1, 2 * m + ip
    r0 = r0 * j
  end do
  r = r0    
  suc = r * df(1)
  sw = 0.0e+00_real64
  do k = 2, nm
    r = r * ( m + k - 1.0e+00_real64 ) * ( m + k + ip - 1.5e+00_real64 ) &
      / ( k - 1.0e+00_real64 ) / ( k + ip - 1.5e+00_real64 )
    suc = suc + r * df(k)
    if ( nm1 < k .and. abs ( suc - sw ) < abs ( suc ) * eps ) then
      exit
    end if
    sw = suc
  end do

  a0 = ( 1.0e+00_real64 - kd / ( x * x ) ) ** ( 0.5e+00_real64 * m ) / suc
  r2f = 0.0e+00_real64
  sw = 0.0e+00_real64
  do k = 1, nm
    l = 2 * k + m - n - 2 + ip
    if ( l == 4 * int ( l / 4 ) ) then
      lg = 1
    else
      lg = -1
    end if

    if ( k == 1 ) then
      r = r0
    else
      r = r * ( m + k - 1.0e+00_real64 ) * ( m + k + ip - 1.5e+00_real64 ) &
        / ( k - 1.0e+00_real64 ) / ( k + ip - 1.5e+00_real64 )
    end if

    np = m + 2 * k - 2 + ip
    r2f = r2f + lg * r * ( df(k) * sy(np) )
    eps1 = abs ( r2f - sw )
    if ( nm1 < k .and. eps1 < abs ( r2f ) * eps ) then
      exit
    end if
    sw = r2f
  end do

  id1 = int ( log10 ( eps1 / abs ( r2f ) + eps ) )
  r2f = r2f * a0

  if ( nm2 <= np ) then
    id = 10
    return
  end if

  b0 = kd * m / x ** 3.0e+00_real64 / ( 1.0e+00_real64 - kd / ( x * x ) ) * r2f
  sud = 0.0e+00_real64
  sw = 0.0e+00_real64
  do k = 1, nm
    l = 2 * k + m - n - 2 + ip
    if ( l == 4 * int ( l / 4 ) ) then
      lg = 1
    else
      lg = -1
    end if
    if (k == 1) then
      r = r0
    else
      r = r * ( m + k - 1.0e+00_real64 ) * ( m + k + ip - 1.5e+00_real64 ) &
        / ( k - 1.0e+00_real64 ) / ( k + ip - 1.5e+00_real64 )
    end if
    np = m + 2 * k - 2 + ip
    sud = sud + lg * r * ( df(k) * dy(np) )
    eps2 = abs ( sud - sw )
    if ( nm1 < k .and. eps2 < abs ( sud ) * eps ) then
      exit
    end if
    sw = sud
  end do

  r2d = b0 + a0 * c * sud
  id2 = int ( log10 ( eps2 / abs ( sud ) + eps ) )
  id = max ( id1, id2 )

  return
end subroutine rmn2l
!> @brief subroutine rmn2so.
!> @return None.
!>
!> @param m [in] Argument m.
!> @param n [in] Argument n.
!> @param c [in] Argument c.
!> @param x [in] Argument x.
!> @param cv [in] Argument cv.
!> @param df [in] Argument df.
!> @param kd [in] Argument kd.
!> @param r2f [inout] Argument r2f.
!> @param r2d [inout] Argument r2d.
subroutine rmn2so ( m, n, c, x, cv, df, kd, r2f, r2d )

!*****************************************************************************80
!
!! RMN2SO: oblate radial functions of the second kind with small argument.
!
!  Discussion:
!
!    This procedure computes oblate radial functions of the second kind
!    with a small argument, Rmn(-ic,ix) and Rmn'(-ic,ix).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    27 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
! 
!  Parameters:
!
!    Input, integer(int32) M, the mode parameter;  M = 0, 1, 2, ...
!
!    Input, integer(int32) N, mode parameter, N = M, M + 1, M + 2, ...
!
!    Input, real(real64) C, spheroidal parameter.
!
!    Input, real(real64) X, the argument.
!
!    Input, real(real64) CV, the characteristic value.
!
!    Input, real(real64) DF(*), the expansion coefficients.
!
!    Input, integer(int32) KD, the function code.
!    1, the prolate function.
!    -1, the oblate function.
!
!    Output, real(real64) R2F, R2D, the values of Rmn(-ic,ix) 
!    and Rmn'(-ic,ix).
!
  implicit none

  real(real64) bk(200)
  real(real64), intent(in) :: c
  real(real64) ck(200)
  real(real64) ck1
  real(real64) ck2
  real(real64), intent(in) :: cv
  real(real64), intent(in) :: df(200)
  real(real64) dn(200)
  real(real64) eps
  real(real64) gd
  real(real64) gf
  real(real64) h0
  integer(int32) ip
  integer(int32) j
  integer(int32), intent(in) :: kd
  integer(int32), intent(in) :: m
  integer(int32), intent(in) :: n
  integer(int32) nm
  real(real64) pi
  real(real64) qs
  real(real64) qt
  real(real64) r1d
  real(real64) r1f
  real(real64), intent(inout) :: r2d
  real(real64), intent(inout) :: r2f
  real(real64) sum
  real(real64) sw
  real(real64), intent(in) :: x

  if ( abs ( df(1) ) <= 1.0e-280_real64 ) then
    r2f = 1.0e+300_real64
    r2d = 1.0e+300_real64
    return
  end if

  eps = 1.0e-14_real64
  pi = 3.141592653589793e+00_real64
  nm = 25 + int ( ( n - m ) / 2 + c )
  if ( n - m == 2 * int ( ( n - m ) / 2 ) ) then
    ip = 0
  else
    ip = 1
  end if

  call sckb ( m, n, c, df, ck )
  call kmn ( m, n, c, cv, kd, df, dn, ck1, ck2 )
  call qstar ( m, n, c, ck, ck1, qs, qt )
  call cbk ( m, n, c, cv, qt, ck, bk )

  if ( x == 0.0e+00_real64 ) then

    sum = 0.0e+00_real64
    sw = 0.0e+00_real64
    do j = 1, nm
      sum = sum + ck(j)
      if ( abs ( sum - sw ) < abs ( sum ) * eps ) then
        exit
      end if
      sw = sum
    end do

    if ( ip == 0 ) then
      r1f = sum / ck1
      r2f = - 0.5e+00_real64 * pi * qs * r1f
      r2d = qs * r1f + bk(1)
    else if ( ip == 1 ) then
       r1d = sum / ck1
       r2f = bk(1)
       r2d = -0.5e+00_real64 * pi * qs * r1d
    end if

    return

  else

    call gmn ( m, n, c, x, bk, gf, gd )
    call rmn1 ( m, n, c, x, df, kd, r1f, r1d )
    h0 = atan ( x ) - 0.5e+00_real64 * pi
    r2f = qs * r1f * h0 + gf
    r2d = qs * ( r1d * h0 + r1f / ( 1.0e+00_real64 + x * x ) ) + gd

  end if

  return
end subroutine rmn2so
!> @brief subroutine rmn2sp.
!> @return None.
!>
!> @param m [in] Argument m.
!> @param n [in] Argument n.
!> @param c [in] Argument c.
!> @param x [in] Argument x.
!> @param cv [in] Argument cv.
!> @param df [in] Argument df.
!> @param kd [in] Argument kd.
!> @param r2f [inout] Argument r2f.
!> @param r2d [inout] Argument r2d.
subroutine rmn2sp ( m, n, c, x, cv, df, kd, r2f, r2d )

!*****************************************************************************80
!
!! RMN2SP: prolate, oblate spheroidal radial functions, kind 2, small argument.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    28 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) M, the mode parameter;  M = 0, 1, 2, ...
!
!    Input, integer(int32) N, mode parameter, N = M, M + 1, M + 2, ...
!
!    Input, real(real64) C, spheroidal parameter.
!
!    Input, real(real64) X, the argument.
!
!    Input, real(real64) CV, the characteristic value.
!
!    Input, real(real64) DF(*), the expansion coefficients.
!
!    Input, integer(int32) KD, the function code.
!    1, the prolate function.
!    -1, the oblate function.
!
!    Output, real(real64) R2F, R2D, the values of the function and 
!    its derivative.
!
  implicit none

  real(real64), intent(in) :: c
  real(real64) ck1
  real(real64) ck2
  real(real64), intent(in) :: cv
  real(real64), intent(in) :: df(200)
  real(real64) dn(200)
  real(real64) eps
  real(real64) ga
  real(real64) gb
  real(real64) gc
  integer(int32) ip
  integer(int32) j
  integer(int32) j1
  integer(int32) j2
  integer(int32) k
  integer(int32), intent(in) :: kd
  integer(int32) ki
  integer(int32) l1
  integer(int32), intent(in) :: m
  integer(int32), intent(in) :: n
  integer(int32) nm
  integer(int32) nm1
  integer(int32) nm2
  integer(int32) nm3
  real(real64) pd(0:251)
  real(real64) pm(0:251)
  real(real64) qd(0:251)
  real(real64) qm(0:251)
  real(real64) r1
  real(real64) r2
  real(real64), intent(inout) :: r2d
  real(real64), intent(inout) :: r2f
  real(real64) r3
  real(real64) r4
  real(real64) sd
  real(real64) sd0
  real(real64) sd1
  real(real64) sd2
  real(real64) sdm
  real(real64) sf
  real(real64) spd1
  real(real64) spd2
  real(real64) spl
  real(real64) su0
  real(real64) su1
  real(real64) su2
  real(real64) sum
  real(real64) sw
  real(real64), intent(in) :: x

  if ( abs ( df(1) ) < 1.0e-280_real64 ) then
    r2f = 1.0e+300_real64
    r2d = 1.0e+300_real64
    return
  end if

  eps = 1.0e-14_real64

  nm1 = int ( ( n - m ) / 2 )

  if ( n - m == 2 * nm1 ) then
    ip = 0
  else
    ip = 1
  end if

  nm = 25 + nm1 + int ( c )
  nm2 = 2 * nm + m
  call kmn ( m, n, c, cv, kd, df, dn, ck1, ck2 )
  call lpmns ( m, nm2, x, pm, pd )
  call lqmns ( m, nm2, x, qm, qd )

  su0 = 0.0e+00_real64
  sw = 0.0e+00_real64
  do k = 1, nm
    j = 2 * k - 2 + m + ip
    su0 = su0 + df(k) * qm(j)
    if ( nm1 < k .and. abs ( su0 - sw ) < abs ( su0 ) * eps ) then
      exit                                                               
    end if
    sw = su0
  end do

  sd0 = 0.0e+00_real64
  sw = 0.0e+00_real64

  do k = 1, nm
    j = 2 * k - 2 + m + ip
    sd0 = sd0 + df(k) * qd(j)
    if ( nm1 < k .and. abs ( sd0 - sw ) < abs ( sd0 ) * eps ) then
      exit
    end if
    sw = sd0
  end do

  su1 = 0.0e+00_real64
  sd1 = 0.0e+00_real64
  do k = 1, m
    j = m - 2 * k + ip
    if ( j < 0 ) then
      j = - j - 1
    end if
    su1 = su1 + dn(k) * qm(j)
    sd1 = sd1 + dn(k) * qd(j)
  end do

  ga = ( ( x - 1.0e+00_real64 ) / ( x + 1.0e+00_real64 ) ) ** ( 0.5e+00_real64 * m )

  do k = 1, m

    j = m - 2 * k + ip

    if ( 0 <= j ) then
      cycle
    end if

    if ( j < 0 ) then
      j = - j - 1
    end if 
    r1 = 1.0e+00_real64
    do j1 = 1, j
      r1 = ( m + j1 ) * r1
    end do
    r2 = 1.0e+00_real64
    do j2 = 1, m - j - 2
      r2 = j2 * r2
    end do
    r3 = 1.0e+00_real64
    sf = 1.0e+00_real64
    do l1 = 1, j
      r3 = 0.5e+00_real64 * r3 * ( - j + l1 - 1.0e+00_real64 ) * ( j + l1 ) &
        / ( ( m + l1 ) * l1 ) * ( 1.0e+00_real64 - x )
      sf = sf + r3
    end do

    if ( m - j <= 1 ) then
      gb = 1.0e+00_real64
    else
      gb = ( m - j - 1.0e+00_real64 ) * r2
    end if

    spl = r1 * ga * gb * sf
    su1 = su1 + ( -1 ) ** ( j + m ) * dn(k) * spl
    spd1 = m / ( x * x - 1.0e+00_real64 ) * spl
    gc = 0.5e+00_real64 * j * ( j + 1.0_real64 ) / ( m + 1.0e+00_real64 )
    sd = 1.0e+00_real64
    r4 = 1.0e+00_real64
    do l1 = 1, j - 1
      r4 = 0.5e+00_real64 * r4 * ( - j + l1 ) * ( j + l1 + 1.0e+00_real64 ) &
        / ( ( m + l1 + 1.0e+00_real64 ) * l1 ) * ( 1.0e+00_real64 - x )
      sd = sd + r4
    end do

    spd2 = r1 * ga * gb * gc * sd
    sd1 = sd1 + ( - 1 ) ** ( j + m ) * dn(k) * ( spd1 + spd2 )

  end do

  su2 = 0.0e+00_real64
  sw = 0.0e+00_real64
  ki = ( 2 * m + 1 + ip ) / 2
  nm3 = nm + ki
  do k = ki, nm3
    j = 2 * k - 1 - m - ip
    su2 = su2 + dn(k) * pm(j)
    if ( m < j .and. &
      abs ( su2 - sw ) < abs ( su2 ) * eps ) then
      exit
    end if
    sw = su2
  end do

  sd2 = 0.0e+00_real64
  sw = 0.0e+00_real64

  do k = ki, nm3
    j = 2 * k - 1 - m - ip
    sd2 = sd2 + dn(k) * pd(j)
    if ( m < j .and. &
      abs ( sd2 - sw ) < abs ( sd2 ) * eps ) then
      exit
    end if
    sw = sd2
  end do

  sum = su0 + su1 + su2
  sdm = sd0 + sd1 + sd2
  r2f = sum / ck2
  r2d = sdm / ck2

  return
end subroutine rmn2sp
!> @brief subroutine rswfo.
!> @return None.
!>
!> @param m [in] Argument m.
!> @param n [in] Argument n.
!> @param c [in] Argument c.
!> @param x [in] Argument x.
!> @param cv [in] Argument cv.
!> @param kf [in] Argument kf.
!> @param r1f [inout] Argument r1f.
!> @param r1d [inout] Argument r1d.
!> @param r2f [inout] Argument r2f.
!> @param r2d [inout] Argument r2d.
subroutine rswfo ( m, n, c, x, cv, kf, r1f, r1d, r2f, r2d )

!*****************************************************************************80
!
!! RSWFO computes prolate spheroidal radial function of first and second kinds.
!
!  Discussion:
!
!    This procedure computes oblate radial functions of the first
!    and second kinds, and their derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    07 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
! 
!  Parameters:
!
!    Input, integer(int32) M, the mode parameter;  M = 0, 1, 2, ...
!
!    Input, integer(int32) N, mode parameter, N = M, M + 1, M + 2, ...
!
!    Input, real(real64) C, spheroidal parameter.
!
!    Input, real(real64) X, the argument.
!
!    Input, real(real64) CV, the characteristic value.
!
!    Input, integer(int32) KF, the function code.
!    1, for the first kind
!    2, for the second kind
!    3, for both the first and second kinds.
!
!    Output, real(real64) R1F, the radial function of the first kind;
!
!    Output, real(real64) R1D, the derivative of the radial function of
!    the first kind;
!
!    Output, real(real64) R2F, the radial function of the second kind;
!
!    Output, real(real64) R2D, the derivative of the radial function of
!    the second kind;
!
  implicit none

  real(real64), intent(in) :: c
  real(real64), intent(in) :: cv
  real(real64) df(200)
  integer(int32) id
  integer(int32) kd
  integer(int32), intent(in) :: kf
  integer(int32), intent(in) :: m
  integer(int32), intent(in) :: n
  real(real64), intent(inout) :: r1d
  real(real64), intent(inout) :: r1f
  real(real64), intent(inout) :: r2d
  real(real64), intent(inout) :: r2f
  real(real64), intent(in) :: x

  kd = -1
  call sdmn ( m, n, c, cv, kd, df )

  if ( kf /= 2 ) then
    call rmn1 ( m, n, c, x, df, kd, r1f, r1d )
  end if

  if ( 1 < kf ) then
    id = 10
    if ( 1.0e-08_real64 < x ) then
      call rmn2l ( m, n, c, x, df, kd, r2f, r2d, id )
    end if
    if ( -1 < id ) then
      call rmn2so ( m, n, c, x, cv, df, kd, r2f, r2d )
    end if
  end if

  return
end subroutine rswfo

!> @brief subroutine rswfp.
!> @return None.
!>
!> @param m [in] Argument m.
!> @param n [in] Argument n.
!> @param c [in] Argument c.
!> @param x [in] Argument x.
!> @param cv [in] Argument cv.
!> @param kf [in] Argument kf.
!> @param r1f [inout] Argument r1f.
!> @param r1d [inout] Argument r1d.
!> @param r2f [inout] Argument r2f.
!> @param r2d [inout] Argument r2d.
subroutine rswfp ( m, n, c, x, cv, kf, r1f, r1d, r2f, r2d )

!*****************************************************************************80
!
!! RSWFP computes prolate spheroidal radial function of first and second kinds.
!
!  Discussion:
!
!    This procedure computes prolate spheriodal radial functions of the
!    first and second kinds, and their derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    07 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) M, the mode parameter;  M = 0, 1, 2, ...
!
!    Input, integer(int32) N, mode parameter, N = M, M + 1, M + 2, ...
!
!    Input, real(real64) C, spheroidal parameter.
!
!    Input, real(real64) X, the argument of the radial function, 1 < X.
!
!    Input, real(real64) CV, the characteristic value.
!
!    Input, integer(int32) KF, the function code.
!    1, for the first kind
!    2, for the second kind
!    3, for both the first and second kinds.
!
!    Output, real(real64) R1F, the radial function of the first kind;
!
!    Output, real(real64) R1D, the derivative of the radial function of
!    the first kind;
!
!    Output, real(real64) R2F, the radial function of the second kind;
!
!    Output, real(real64) R2D, the derivative of the radial function of
!    the second kind;
!
  implicit none

  real(real64), intent(in) :: c
  real(real64), intent(in) :: cv
  real(real64) df(200)
  integer(int32) id
  integer(int32) kd
  integer(int32), intent(in) :: kf
  integer(int32), intent(in) :: m
  integer(int32), intent(in) :: n
  real(real64), intent(inout) :: r1d
  real(real64), intent(inout) :: r1f
  real(real64), intent(inout) :: r2d
  real(real64), intent(inout) :: r2f
  real(real64), intent(in) :: x

  kd = 1
  call sdmn ( m, n, c, cv, kd, df )

  if ( kf /= 2 ) then
    call rmn1 ( m, n, c, x, df, kd, r1f, r1d )
  end if

  if ( 1 < kf ) then
    call rmn2l ( m, n, c, x, df, kd, r2f, r2d, id )
    if ( -8 < id ) then
      call rmn2sp ( m, n, c, x, cv, df, kd, r2f, r2d )
    end if
  end if

  return
end subroutine rswfp
!> @brief subroutine scka.
!> @return None.
!>
!> @param m [in] Argument m.
!> @param n [in] Argument n.
!> @param c [in] Argument c.
!> @param cv [in] Argument cv.
!> @param kd [in] Argument kd.
!> @param ck [inout] Argument ck.
pure subroutine scka ( m, n, c, cv, kd, ck )

!*****************************************************************************80
!
!! SCKA: expansion coefficients for prolate and oblate spheroidal functions.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    22 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) M, the mode parameter.
!
!    Input, integer(int32) N, the mode parameter.
!
!    Input, real(real64) C, the spheroidal parameter.
!
!    Input, real(real64) CV, the characteristic value.
!
!    Input, integer(int32) KD, the function code.
!    1, the prolate function.
!    -1, the oblate function.
!
!    Output, real(real64) CK(*), the expansion coefficients.
!    CK(1), CK(2),... correspond to c0, c2,..., and so on.
!       
  implicit none

  real(real64), intent(in) :: c
  real(real64), intent(inout) :: ck(200)
  real(real64) cc
  real(real64) cs
  real(real64), intent(in) :: cv
  real(real64) f
  real(real64) f0
  real(real64) f1
  real(real64) f2
  real(real64) fl
  real(real64) fs
  integer(int32) ip
  integer(int32) j
  integer(int32) k
  integer(int32) k1
  integer(int32) kb
  integer(int32), intent(in) :: kd
  integer(int32), intent(in) :: m
  integer(int32), intent(in) :: n
  integer(int32) nm
  real(real64) r1
  real(real64) r2
  real(real64) s0
  real(real64) su1
  real(real64) su2

  cc = max ( c, 1.0e-10_real64 )

  nm = 25 + int ( ( n - m ) / 2 + cc )
  cs = cc * cc * kd

  if ( n - m == 2 * int ( ( n - m ) / 2 ) ) then
    ip = 0
  else
    ip = 1
  end if

  fs = 1.0e+00_real64
  f1 = 0.0e+00_real64
  f0 = 1.0e-100_real64
  kb = 0
  ck(nm+1) = 0.0e+00_real64

  do k = nm, 1, -1

    f = ((( 2.0e+00_real64 * k + m + ip ) &
      * ( 2.0e+00_real64 * k + m + 1.0e+00_real64 + ip ) - cv + cs ) * f0 &
      - 4.0e+00_real64 * ( k + 1.0e+00_real64 ) * ( k + m + 1.0e+00_real64 ) * f1 ) / cs

    if ( abs ( ck(k+1) ) < abs ( f ) ) then

      ck(k) = f
      f1 = f0
      f0 = f

      if ( 1.0e+100_real64 < abs ( f ) ) then
        do k1 = nm, k, -1
          ck(k1) = ck(k1) * 1.0e-100_real64
        end do
        f1 = f1 * 1.0e-100_real64
        f0 = f0 * 1.0e-100_real64
      end if

    else

      kb = k
      fl = ck(k+1)
      f1 = 1.0e+00_real64
      f2 = 0.25e+00_real64 * ( ( m + ip ) * ( m + ip + 1.0e+00_real64 ) &
        - cv + cs ) / ( m + 1.0e+00_real64 ) * f1
      ck(1) = f1

      if ( kb == 1 ) then
        fs = f2
      else if (kb == 2 ) then
        ck(2) = f2
        fs = 0.125e+00_real64 * ( ( ( m + ip + 2.0e+00_real64 ) &
          * ( m + ip + 3.0e+00_real64 ) - cv + cs ) * f2 &
          - cs * f1 ) / ( m + 2.0e+00_real64 )
      else
        ck(2) = f2
        do j = 3, kb + 1
          f = 0.25e+00_real64 * ( ( ( 2.0e+00_real64 * j + m + ip - 4.0e+00_real64 ) &
            * ( 2.0e+00_real64 * j + m + ip - 3.0e+00_real64 ) - cv + cs ) * f2 &
            - cs * f1 ) / ( ( j - 1.0e+00_real64 ) * ( j + m - 1.0e+00_real64 ) )
          if ( j <= kb ) then
            ck(j) = f
          end if
          f1 = f2
          f2 = f
        end do
        fs = f
      end if

      exit

    end if

  end do

  su1 = 0.0e+00_real64
  do k = 1, kb
    su1 = su1 + ck(k)
  end do

  su2 = 0.0e+00_real64
  do k = kb + 1, nm
    su2 = su2 + ck(k)
  end do

  r1 = 1.0e+00_real64
  do j = 1, ( n + m + ip ) / 2
    r1 = r1 * ( j + 0.5e+00_real64 * ( n + m + ip ) )
  end do

  r2 = 1.0e+00_real64
  do j = 1, ( n - m - ip ) / 2
    r2 = - r2 * j
  end do

  if ( kb == 0 ) then
    s0 = r1 / ( 2.0e+00_real64 ** n * r2 * su2 )
  else
    s0 = r1 / ( 2.0e+00_real64 ** n * r2 * ( fl / fs * su1 + su2 ) )
  end if

  do k = 1, kb
    ck(k) = fl / fs * s0 * ck(k)
  end do

  do k = kb + 1, nm
    ck(k) = s0 * ck(k)
  end do

  return
end subroutine scka
!> @brief subroutine sckb.
!> @return None.
!>
!> @param m [in] Argument m.
!> @param n [in] Argument n.
!> @param c [in] Argument c.
!> @param df [in] Argument df.
!> @param ck [inout] Argument ck.
pure subroutine sckb ( m, n, c, df, ck )

!*****************************************************************************80
!
!! SCKB: expansion coefficients for prolate and oblate spheroidal functions.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    15 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) M, the mode parameter.
!
!    Input, integer(int32) N, the mode parameter.
!
!    Input, real(real64) C, the spheroidal parameter.
!
!    Input, real(real64) DF(*), the expansion coefficients DK.
!
!    Output, real(real64) CK(*), the expansion coefficients CK.
!
  implicit none

  real(real64), intent(in) :: c
  real(real64), intent(inout) :: ck(200)
  real(real64) cc
  real(real64) d1
  real(real64) d2
  real(real64) d3
  real(real64), intent(in) :: df(200)
  real(real64) fac
  integer(int32) i
  integer(int32) i1
  integer(int32) i2
  integer(int32) ip
  integer(int32) k
  integer(int32), intent(in) :: m
  integer(int32), intent(in) :: n
  integer(int32) nm
  real(real64) r
  real(real64) r1
  real(real64) reg
  real(real64) sum
  real(real64) sw

  cc = max ( c, 1.0e-10_real64 )

  nm = 25 + int ( 0.5e+00_real64 * ( n - m ) + cc )
 
  if ( n - m == 2 * int ( ( n - m ) / 2 ) ) then
    ip = 0
  else
    ip = 1
  end if

  if ( 80 < m + nm ) then
    reg = 1.0e-200_real64
  else
    reg = 1.0e+00_real64
  end if

  fac = - 0.5e+00_real64 ** m

  do k = 0, nm - 1

    fac = - fac
    i1 = 2 * k + ip + 1
    r = reg
    do i = i1, i1 + 2 * m - 1
      r = r * i
    end do 

    i2 = k + m + ip
    do i = i2, i2 + k - 1
      r = r * ( i + 0.5e+00_real64 )
    end do

    sum = r * df(k+1)
    do i = k + 1, nm
      d1 = 2.0e+00_real64 * i + ip
      d2 = 2.0e+00_real64 * m + d1
      d3 = i + m + ip - 0.5e+00_real64
      r = r * d2 * ( d2 - 1.0e+00_real64 ) * i * ( d3 + k ) &
        / ( d1 * ( d1 - 1.0e+00_real64 ) * ( i - k ) * d3 )
      sum = sum + r * df(i+1)
      if ( abs ( sw - sum ) < abs ( sum ) * 1.0e-14_real64 ) then
        exit
      end if
      sw = sum
    end do

    r1 = reg
    do i = 2, m + k
      r1 = r1 * i
    end do

    ck(k+1) = fac * sum / r1

  end do

  return
end subroutine sckb
!> @brief subroutine sdmn.
!> @return None.
!>
!> @param m [in] Argument m.
!> @param n [in] Argument n.
!> @param c [in] Argument c.
!> @param cv [in] Argument cv.
!> @param kd [in] Argument kd.
!> @param df [inout] Argument df.
pure subroutine sdmn ( m, n, c, cv, kd, df )

!*****************************************************************************80
!
!! SDMN: expansion coefficients for prolate and oblate spheroidal functions.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    29 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) M, the mode parameter.
!
!    Input, integer(int32) N, the mode parameter.
!
!    Input, real(real64) C, the spheroidal parameter.
!
!    Input, real(real64) CV, the characteristic value.
!
!    Input, integer(int32) KD, the function code.
!    1, the prolate function.
!    -1, the oblate function.
!
!    Output, real(real64) DF(*), expansion coefficients;
!    DF(1), DF(2), ... correspond to d0, d2, ... for even n-m and d1,
!    d3, ... for odd n-m
!
  implicit none

  real(real64) a(200)
  real(real64), intent(in) :: c
  real(real64) cs
  real(real64), intent(in) :: cv
  real(real64) d(200)
  real(real64) d2k
  real(real64), intent(inout) :: df(200)
  real(real64) dk0
  real(real64) dk1
  real(real64) dk2
  real(real64) f
  real(real64) f0
  real(real64) f1
  real(real64) f2
  real(real64) fl
  real(real64) fs
  real(real64) g(200)
  integer(int32) i
  integer(int32) ip
  integer(int32) j
  integer(int32) k
  integer(int32) k1
  integer(int32) kb
  integer(int32), intent(in) :: kd
  integer(int32), intent(in) :: m
  integer(int32), intent(in) :: n
  integer(int32) nm
  real(real64) r1
  real(real64) r3
  real(real64) r4
  real(real64) s0
  real(real64) su1
  real(real64) su2
  real(real64) sw

  nm = 25 + int ( 0.5e+00_real64 * ( n - m ) + c )

  if ( c < 1.0e-10_real64 ) then
     do i = 1, nm
       df(i) = 0e+00_real64
     end do
     df((n-m)/2+1) = 1.0e+00_real64
     return
  end if   

  cs = c * c * kd

  if ( n - m == 2 * int ( ( n - m ) / 2 ) ) then
    ip = 0
  else
    ip = 1
  end if

  do i = 1, nm + 2
    if ( ip == 0 ) then
      k = 2 * ( i - 1 )
    else
      k = 2 * i - 1
    end if
    dk0 = m + k
    dk1 = m + k + 1
    dk2 = 2 * ( m + k )
    d2k = 2 * m + k
    a(i) = ( d2k + 2.0e+00_real64 ) * ( d2k + 1.0e+00_real64 ) &
      / ( ( dk2 + 3.0e+00_real64 ) * ( dk2 + 5.0e+00_real64 ) ) * cs
    d(i) = dk0 * dk1 &
      + ( 2.0e+00_real64 * dk0 * dk1 - 2.0e+00_real64 * m * m - 1.0e+00_real64 ) &
      / ( ( dk2 - 1.0e+00_real64 ) * ( dk2 + 3.0e+00_real64 ) ) * cs
    g(i) = k * ( k - 1.0e+00_real64 ) / ( ( dk2 - 3.0e+00_real64 ) &
      * ( dk2 - 1.0e+00_real64 ) ) * cs
  end do

  fs = 1.0e+00_real64
  f1 = 0.0e+00_real64
  f0 = 1.0e-100_real64
  kb = 0
  df(nm+1) = 0.0e+00_real64

  do k = nm, 1, -1

    f = - ( ( d(k+1) - cv ) * f0 + a(k+1) * f1 ) / g(k+1)

    if ( abs ( df(k+1) ) < abs ( f ) ) then

      df(k) = f
      f1 = f0
      f0 = f
      if ( 1.0e+100_real64 < abs ( f ) ) then
        do k1 = k, nm
          df(k1) = df(k1) * 1.0e-100_real64
        end do
        f1 = f1 * 1.0e-100_real64
        f0 = f0 * 1.0e-100_real64
      end if  

    else

      kb = k
      fl = df(k+1)
      f1 = 1.0e-100_real64
      f2 = - ( d(1) - cv ) / a(1) * f1
      df(1) = f1

      if ( kb == 1 ) then

        fs = f2

      else if ( kb == 2 ) then

        df(2) = f2
        fs = - ( ( d(2) - cv ) * f2 + g(2) * f1 ) / a(2)

      else 

        df(2) = f2
        do j = 3, kb + 1
          f = - ( ( d(j-1) - cv ) * f2 + g(j-1) * f1 ) / a(j-1)
          if ( j <= kb ) then
            df(j) = f
          end if
          if ( 1.0e+100_real64 < abs ( f ) ) then
            do k1 = 1, j
              df(k1) = df(k1) * 1.0e-100_real64
            end do
            f = f * 1.0e-100_real64
            f2 = f2 * 1.0e-100_real64
          end if  
          f1 = f2
          f2 = f
        end do
        fs = f

      end if

      exit

    end if

  end do

  su1 = 0.0e+00_real64

  r1 = 1.0e+00_real64
  do j = m + ip + 1, 2 * ( m + ip )
    r1 = r1 * j
  end do

  su1 = df(1) * r1
  do k = 2, kb
    r1 = - r1 * ( k + m + ip - 1.5e+00_real64 ) / ( k - 1.0e+00_real64 )
    su1 = su1 + r1 * df(k)
  end do

  su2 = 0.0e+00_real64
  do k = kb + 1, nm
    if ( k /= 1 ) then
      r1 = - r1 * ( k + m + ip - 1.5e+00_real64 ) / ( k - 1.0e+00_real64 )
    end if
    su2 = su2 + r1 * df(k)
    if ( abs ( sw - su2 ) < abs ( su2 ) * 1.0e-14_real64 ) then
      exit
    end if
    sw = su2
  end do

  r3 = 1.0e+00_real64
  do j = 1, ( m + n + ip ) / 2
    r3 = r3 * ( j + 0.5e+00_real64 * ( n + m + ip ) )
  end do

  r4 = 1.0e+00_real64
  do j = 1, ( n - m - ip ) / 2
    r4 = -4.0e+00_real64 * r4 * j
  end do

  s0 = r3 / ( fl * ( su1 / fs ) + su2 ) / r4
  do k = 1, kb
    df(k) = fl / fs * s0 * df(k)
  end do

  do k = kb + 1, nm
    df(k) = s0 * df(k)
  end do

  return
end subroutine sdmn
!> @brief subroutine segv.
!> @return None.
!>
!> @param m [in] Argument m.
!> @param n [in] Argument n.
!> @param c [in] Argument c.
!> @param kd [in] Argument kd.
!> @param cv [inout] Argument cv.
!> @param eg [inout] Argument eg.
pure subroutine segv ( m, n, c, kd, cv, eg )

!*****************************************************************************80
!
!! SEGV computes the characteristic values of spheroidal wave functions.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    28 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) M, the mode parameter.
!
!    Input, integer(int32) N, the mode parameter.
!
!    Input, real(real64) C, the spheroidal parameter.
!
!    Input, integer(int32) KD, the function code.
!    1, the prolate function.
!    -1, the oblate function.
!
!    Output, real(real64) CV, the characteristic value.
!
!    Output, real(real64) EG(*), the characteristic value for 
!    mode parameters m and n.  ( L = n - m + 1 )
!
  implicit none

  real(real64) a(300)
  real(real64) b(100)
  real(real64), intent(in) :: c
  real(real64) cs
  real(real64), intent(inout) :: cv
  real(real64) cv0(100)
  real(real64) d(300)
  real(real64) d2k
  real(real64) dk0
  real(real64) dk1
  real(real64) dk2
  real(real64) e(300)
  real(real64), intent(inout) :: eg(200)
  real(real64) f(300)
  real(real64) g(300)
  real(real64) h(100)
  integer(int32) i
  integer(int32) icm
  integer(int32) j
  integer(int32) k
  integer(int32) k1
  integer(int32), intent(in) :: kd
  integer(int32) l
  integer(int32), intent(in) :: m
  integer(int32), intent(in) :: n
  integer(int32) nm
  integer(int32) nm1
  real(real64) s
  real(real64) t
  real(real64) t1
  real(real64) x1
  real(real64) xa
  real(real64) xb

  if ( c < 1.0e-10_real64 ) then
    do i = 1, n
      eg(i) = ( i + m ) * ( i + m - 1.0e+00_real64 )
    end do
    cv = eg(n-m+1)
    return
  end if

  icm = ( n - m + 2 ) / 2
  nm = 10 + int ( 0.5e+00_real64 * ( n - m ) + c )
  cs = c * c * kd

  do l = 0, 1

    do i = 1, nm
      if ( l == 0 ) then
        k = 2 * ( i - 1 )
      else
        k = 2 * i - 1
      end if
      dk0 = m + k
      dk1 = m + k + 1
      dk2 = 2 * ( m + k )
      d2k = 2 * m + k
      a(i) = ( d2k + 2.0e+00_real64 ) * ( d2k + 1.0e+00_real64 ) &
        / ( ( dk2 + 3.0e+00_real64 ) * ( dk2 + 5.0e+00_real64 ) ) * cs
      d(i) = dk0 * dk1 + ( 2.0e+00_real64 * dk0 * dk1 &
        - 2.0_real64 * m * m - 1.0e+00_real64 ) &
        / ( ( dk2 - 1.0e+00_real64 ) * ( dk2 + 3.0e+00_real64 ) ) * cs
      g(i) = k * ( k - 1.0e+00_real64 ) / ( ( dk2 - 3.0e+00_real64 ) &
        * ( dk2 - 1.0e+00_real64 ) ) * cs
    end do

    do k = 2, nm
      e(k) = sqrt ( a(k-1) * g(k) )
      f(k) = e(k) * e(k)
    end do

    f(1) = 0.0e+00_real64
    e(1) = 0.0e+00_real64
    xa = d(nm) + abs ( e(nm) )
    xb = d(nm) - abs ( e(nm) )
    nm1 = nm - 1
    do i = 1, nm1
      t = abs ( e(i) ) + abs ( e(i+1) )
      t1 = d(i) + t
      if ( xa < t1 ) then
        xa = t1
      end if
      t1 = d(i) - t
      if ( t1 < xb ) then
        xb = t1
      end if
    end do

    do i = 1, icm
      b(i) = xa
      h(i) = xb
    end do

    do k = 1, icm

      do k1 = k, icm
        if ( b(k1) < b(k) ) then
          b(k) = b(k1)
          exit
        end if
      end do

      if ( k /= 1 .and. h(k) < h(k-1) ) then
        h(k) = h(k-1)
      end if

      do

        x1 = ( b(k) + h(k) ) /2.0e+00_real64
        cv0(k) = x1

        if ( abs ( ( b(k) - h(k) ) / x1 ) < 1.0e-14_real64 ) then
          exit
        end if

        j = 0
        s = 1.0e+00_real64

        do i = 1, nm

          if ( s == 0.0e+00_real64 ) then
            s = s + 1.0e-30_real64
          end if
          t = f(i) / s
          s = d(i) - t - x1
          if ( s < 0.0e+00_real64 ) then
            j = j + 1
          end if
        end do
 
        if ( j < k ) then

          h(k) = x1

        else

          b(k) = x1
          if ( icm <= j ) then
            b(icm) = x1
          else
            if ( h(j+1) < x1 ) then
              h(j+1) = x1
            end if
            if ( x1 < b(j) ) then
              b(j) = x1
            end if
          end if

        end if

      end do

      cv0(k) = x1

      if ( l == 0 ) then
        eg(2*k-1) = cv0(k)
      else
        eg(2*k) = cv0(k)
      end if

    end do

  end do

  cv = eg(n-m+1)

  return
end subroutine segv
!> @brief subroutine sphi.
!> @return None.
!>
!> @param n [in] Argument n.
!> @param x [in] Argument x.
!> @param nm [inout] Argument nm.
!> @param si [inout] Argument si.
!> @param di [inout] Argument di.
subroutine sphi ( n, x, nm, si, di )

!*****************************************************************************80
!
!! SPHI computes spherical Bessel functions in(x) and their derivatives in'(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    18 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) N, the order of In(X).
!
!    Input, real(real64) X, the argument.
!
!    Output, integer(int32) NM, the highest order computed.
!
!    Output, real(real64) SI(0:N), DI(0:N), the values and derivatives
!    of the function of orders 0 through N.
!
  implicit none

  integer(int32), intent(in) :: n

  real(real64) cs
  real(real64), intent(inout) :: di(0:n)
  real(real64) f
  real(real64) f0
  real(real64) f1
  integer(int32) k
  integer(int32) m
  ! integer(int32) msta1
  ! integer(int32) msta2
  integer(int32), intent(inout) :: nm
  real(real64), intent(inout) :: si(0:n)
  real(real64) si0
  real(real64), intent(in) :: x

  nm = n

  if ( abs ( x ) < 1.0e-100_real64 ) then
    do k = 0, n
      si(k) = 0.0e+00_real64
      di(k) = 0.0e+00_real64
    end do
    si(0) = 1.0e+00_real64
    di(1) = 0.333333333333333e+00_real64
    return
  end if

  si(0) = sinh ( x ) / x
  si(1) = -( sinh ( x ) / x - cosh ( x ) ) / x
  si0 = si(0)

  if ( 2 <= n ) then

    m = msta1 ( x, 200 )
    if ( m < n ) then
      nm = m
    else
      m = msta2 ( x, n, 15 )
    end if
    f0 = 0.0e+00_real64
    f1 = 1.0e+00_real64-100
    do k = m, 0, -1
      f = ( 2.0e+00_real64 * k + 3.0e+00_real64 ) * f1 / x + f0
      if ( k <= nm ) then
        si(k) = f
      end if
      f0 = f1
      f1 = f
    end do
    cs = si0 / f
    do k = 0, nm
      si(k) = cs * si(k)
    end do

  end if

  di(0) = si(1)
  do k = 1, nm
    di(k) = si(k-1) - ( k + 1.0e+00_real64 ) / x * si(k)
  end do

  return
end subroutine sphi
!> @brief subroutine sphj.
!> @return None.
!>
!> @param n [in] Argument n.
!> @param x [in] Argument x.
!> @param nm [inout] Argument nm.
!> @param sj [inout] Argument sj.
!> @param dj [inout] Argument dj.
subroutine sphj ( n, x, nm, sj, dj )

!*****************************************************************************80
!
!! SPHJ computes spherical Bessel functions jn(x) and their derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    12 January 2016
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin.
!    Modifications suggested by Vincent Lagage, 12 January 2016._real64
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) N, the order.
!
!    Input, real(real64) X, the argument.
!
!    Output, integer(int32) NM, the highest order computed.
!
!    Output, real(real64) SJ(0:N), the values of jn(x).
!
!    Output, real(real64) DJ(0:N), the values of jn'(x).
!
  implicit none

  integer(int32), intent(in) :: n

  real(real64) cs
  real(real64), intent(inout) :: dj(0:n)
  real(real64) f
  real(real64) f0
  real(real64) f1
  integer(int32) k
  integer(int32) m
  ! integer(int32) msta1
  ! integer(int32) msta2
  integer(int32), intent(inout) :: nm
  real(real64) sa
  real(real64) sb
  real(real64), intent(inout) :: sj(0:n)
  real(real64), intent(in) :: x

  nm = n
!
!  Original code.
!
  if ( .true. ) then

    if ( abs ( x ) <= 1.0e-100_real64 ) then
      do k = 0, n
        sj(k) = 0.0e+00_real64
        dj(k) = 0.0e+00_real64
      end do
      sj(0) = 1.0e+00_real64
      dj(1) = 0.3333333333333333e+00_real64
      return
    end if
!
!  Updated code.
!
  else

    if ( abs ( x ) <= 1.0e-16_real64 ) then
      do k = 0, n
        sj(k) = 0.0e+00_real64
        dj(k) = 0.0e+00_real64
      end do
      sj(0) = 1.0e+00_real64
      if ( 0 < n ) then
        do k = 1, n
          sj(k) = sj(k-1) * x / real ( 2 * k + 1, kind = real64 )
        end do
        dj(1) = 1.0e+00_real64 / 3.0e+00_real64
      end if
      return
    end if

  end if

  sj(0) = sin ( x ) / x
  sj(1) = ( sj(0) - cos ( x ) ) / x

  if ( 2 <= n ) then

    sa = sj(0)
    sb = sj(1)
    m = msta1 ( x, 200 )
    if ( m < n ) then
      nm = m
    else
      m = msta2 ( x, n, 15 )
    end if

    f0 = 0.0e+00_real64
    f1 = 1.0e+00_real64-100
    do k = m, 0, -1
      f = ( 2.0e+00_real64 * k + 3.0e+00_real64 ) * f1 / x - f0
      if ( k <= nm ) then
        sj(k) = f
      end if
      f0 = f1
      f1 = f
    end do

    if ( abs ( sa ) <= abs ( sb ) ) then
      cs = sb / f0
    else
      cs = sa / f
    end if

    do k = 0, nm
      sj(k) = cs * sj(k)
    end do

  end if      

  dj(0) = ( cos(x) - sin(x) / x ) / x
  do k = 1, nm
    dj(k) = sj(k-1) - ( k + 1.0e+00_real64 ) * sj(k) / x
  end do

  return
end subroutine sphj
!> @brief subroutine sphk.
!> @return None.
!>
!> @param n [in] Argument n.
!> @param x [in] Argument x.
!> @param nm [inout] Argument nm.
!> @param sk [inout] Argument sk.
!> @param dk [inout] Argument dk.
pure subroutine sphk ( n, x, nm, sk, dk )

!*****************************************************************************80
!
!! SPHK computes modified spherical Bessel functions kn(x) and derivatives.
!
!  Discussion:
!
!    This procedure computes modified spherical Bessel functions
!    of the second kind, kn(x) and kn'(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    15 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) N, the order.
!
!    Input, real(real64) X, the argument.
!
!    Output, integer(int32) NM, the highest order computed.
!
!    Output, real(real64) SK(0:N), DK(0:N), the values of kn(x) and kn'(x).
!
  implicit none

  integer(int32), intent(in) :: n

  real(real64), intent(inout) :: dk(0:n)
  real(real64) f
  real(real64) f0
  real(real64) f1
  integer(int32) k
  integer(int32), intent(inout) :: nm
  real(real64), intent(inout) :: sk(0:n)
  real(real64) pi
  real(real64), intent(in) :: x

  pi = 3.141592653589793e+00_real64
  nm = n
  if ( x < 1.0e-60_real64 ) then
    do k = 0,n
      sk(k) = 1.0e+300_real64
      dk(k) = -1.0e+300_real64
    end do
    return
  end if

  sk(0) = 0.5e+00_real64 * pi / x * exp ( - x )
  sk(1) = sk(0) * ( 1.0e+00_real64 + 1.0e+00_real64 / x )
  f0 = sk(0)
  f1 = sk(1)
  do k = 2, n
    f = ( 2.0e+00_real64 * k - 1.0e+00_real64 ) * f1 / x + f0
    sk(k) = f
    if ( 1.0e+300_real64 < abs ( f ) ) then
      exit
    end if
    f0 = f1
    f1 = f
  end do

  nm = k - 1

  dk(0) = -sk(1)
  do k = 1, nm
    dk(k) = -sk(k-1) - ( k + 1.0e+00_real64 ) / x * sk(k)
  end do

  return
end subroutine sphk
!> @brief subroutine sphy.
!> @return None.
!>
!> @param n [in] Argument n.
!> @param x [in] Argument x.
!> @param nm [inout] Argument nm.
!> @param sy [inout] Argument sy.
!> @param dy [inout] Argument dy.
pure subroutine sphy ( n, x, nm, sy, dy )

!*****************************************************************************80
!
!! SPHY computes spherical Bessel functions yn(x) and their derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    15 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer(int32) N, the order.
!
!    Input, real(real64) X, the argument.
!
!    Output, integer(int32) NM, the highest order computed.
!
!    Output, real(real64) SY(0:N), DY(0:N), the values of yn(x) and yn'(x).
! 
  implicit none

  integer(int32), intent(in) :: n

  real(real64), intent(inout) :: dy(0:n)
  real(real64) f
  real(real64) f0
  real(real64) f1
  integer(int32) k
  integer(int32), intent(inout) :: nm
  real(real64), intent(inout) :: sy(0:n)
  real(real64), intent(in) :: x

  nm = n

  if ( x < 1.0e-60_real64 ) then
    do k = 0, n
      sy(k) = -1.0e+300_real64
      dy(k) = 1.0e+300_real64
    end do
    return
  end if

  sy(0) = - cos ( x ) / x
  sy(1) = ( sy(0) - sin ( x ) ) / x
  f0 = sy(0)
  f1 = sy(1)
  do k = 2, n
    f = ( 2.0e+00_real64 * k - 1.0e+00_real64 ) * f1 / x - f0
    sy(k) = f
    if ( 1.0e+300_real64 <= abs ( f ) ) then
      exit
    end if              
    f0 = f1
    f1 = f
  end do

  nm = k - 1
  dy(0) = ( sin ( x ) + cos ( x ) / x ) / x
  do k = 1, nm
    dy(k) = sy(k-1) - ( k + 1.0e+00_real64 ) * sy(k) / x
  end do

  return
end subroutine sphy
!> @brief subroutine stvh0.
!> @return None.
!>
!> @param x [in] Argument x.
!> @param sh0 [inout] Argument sh0.
pure subroutine stvh0 ( x, sh0 )

!*****************************************************************************80
!
!! STVH0 computes the Struve function H0(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    22 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) X, the argument.
!
!    Output, real(real64) SH0, the value of H0(x).
!
  implicit none

  real(real64) a0
  real(real64) by0
  integer(int32) k
  integer(int32) km
  real(real64) p0
  real(real64) pi
  real(real64) q0
  real(real64) r
  real(real64) s
  real(real64), intent(inout) :: sh0
  real(real64) t
  real(real64) t2
  real(real64) ta0
  real(real64), intent(in) :: x

  pi = 3.141592653589793e+00_real64
  s = 1.0e+00_real64
  r = 1.0e+00_real64

  if ( x <= 20.0e+00_real64 ) then
    a0 = 2.0e+00_real64 * x / pi
    do k = 1, 60
      r = - r * x / ( 2.0e+00_real64 * k + 1.0e+00_real64 ) * x &
        / ( 2.0e+00_real64 * k + 1.0e+00_real64 )
      s = s + r
      if ( abs ( r ) < abs ( s ) * 1.0e-12_real64 ) then
        exit
      end if
    end do

    sh0 = a0 * s

  else

    if ( x < 50.0e+00_real64 ) then
      km = int ( 0.5e+00_real64 * ( x + 1.0e+00_real64 ) )
    else
      km = 25
    end if

    do k = 1, km
      r = - r * ( ( 2.0e+00_real64 * k - 1.0e+00_real64 ) / x ) ** 2
      s = s + r
      if ( abs ( r ) < abs ( s ) * 1.0e-12_real64 ) then
        exit
      end if
    end do

    t = 4.0e+00_real64 / x
    t2 = t * t

    p0 = (((( &
      - 0.37043e-05_real64     * t2 &
      + 0.173565e-04_real64 )  * t2 &
      - 0.487613e-04_real64 )  * t2 &
      + 0.17343e-03_real64 )   * t2 &
      - 0.1753062e-02_real64 ) * t2 &
      + 0.3989422793e+00_real64

    q0 = t * ((((( &
        0.32312e-05_real64     * t2 &
      - 0.142078e-04_real64 )  * t2 &
      + 0.342468e-04_real64 )  * t2 &
      - 0.869791e-04_real64 )  * t2 &
      + 0.4564324e-03_real64 ) * t2 &
      - 0.0124669441e+00_real64 )

    ta0 = x - 0.25e+00_real64 * pi
    by0 = 2.0e+00_real64 / sqrt ( x ) &
      * ( p0 * sin ( ta0 ) + q0 * cos ( ta0 ) )
    sh0 = 2.0e+00_real64 / ( pi * x ) * s + by0

  end if

  return
end subroutine stvh0
!> @brief subroutine stvh1.
!> @return None.
!>
!> @param x [in] Argument x.
!> @param sh1 [inout] Argument sh1.
pure subroutine stvh1 ( x, sh1 )

!*****************************************************************************80
!
!! STVH1 computes the Struve function H1(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    22 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) X, the argument.
!
!    Output, real(real64) SH1, the value of H1(x).
!
  implicit none

  real(real64) a0
  real(real64) by1
  integer(int32) k
  integer(int32) km
  real(real64) p1
  real(real64) pi
  real(real64) q1
  real(real64) r
  real(real64) s
  real(real64), intent(inout) :: sh1
  real(real64) t
  real(real64) t2
  real(real64) ta1
  real(real64), intent(in) :: x

  pi = 3.141592653589793e+00_real64
  r = 1.0e+00_real64

  if ( x <= 20.0e+00_real64 ) then

    s = 0.0e+00_real64
    a0 = - 2.0e+00_real64 / pi
    do k = 1, 60
      r = - r * x * x / ( 4.0e+00_real64 * k * k - 1.0e+00_real64 )
      s = s + r
      if ( abs ( r ) < abs ( s ) * 1.0e-12_real64 ) then
        exit
      end if
    end do

    sh1 = a0 * s

  else

    s = 1.0e+00_real64

    if ( x <= 50.0e+00_real64 ) then
      km = int ( 0.5e+00_real64 * x )
    else
      km = 25
    end if

    do k = 1, km
      r = - r * ( 4.0e+00_real64 * k * k - 1.0e+00_real64 ) / ( x * x )
      s = s + r
      if ( abs ( r ) < abs ( s ) * 1.0e-12_real64 ) then
        exit
      end if
    end do

    t = 4.0e+00_real64 / x
    t2 = t * t

    p1 = (((( &
        0.42414e-05_real64      * t2 &
      - 0.20092e-04_real64 )    * t2 &
      + 0.580759e-04_real64 )   * t2 &
      - 0.223203e-03_real64 )   * t2 &
      + 0.29218256e-02_real64 ) * t2 &
      + 0.3989422819e+00_real64

    q1 = t * ((((( &
      - 0.36594e-05_real64     * t2 &
      + 0.1622e-04_real64 )    * t2 &
      - 0.398708e-04_real64 )  * t2 &
      + 0.1064741e-03_real64 ) * t2 &
      - 0.63904e-03_real64 )   * t2 &
      + 0.0374008364e+00_real64 )

    ta1 = x - 0.75e+00_real64 * pi
    by1 = 2.0e+00_real64 / sqrt ( x ) * ( p1 * sin ( ta1 ) + q1 * cos ( ta1 ) )
    sh1 = 2.0e+00_real64 / pi * ( 1.0e+00_real64 + s / ( x * x ) ) + by1

  end if

  return
end subroutine stvh1
!> @brief subroutine stvhv.
!> @return None.
!>
!> @param v [in] Argument v.
!> @param x [in] Argument x.
!> @param hv [inout] Argument hv.
 subroutine stvhv ( v, x, hv )

!*****************************************************************************80
!
!! STVHV computes the Struve function Hv(x) with arbitrary order v.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    24 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) V, the order of the function.
!
!    Input, real(real64) X, the argument.
!
!    Output, real(real64) HV, the value of Hv(x).
!
  implicit none

  real(real64) bf
  real(real64) bf0
  real(real64) bf1
  real(real64) by0
  real(real64) by1
  real(real64) byv
  real(real64) ga
  real(real64) gb
  real(real64), intent(inout) :: hv
  integer(int32) k
  integer(int32) l
  integer(int32) n
  real(real64) pi
  real(real64) pu0
  real(real64) pu1
  real(real64) qu0
  real(real64) qu1
  real(real64) r1
  real(real64) r2
  real(real64) s
  real(real64) s0
  real(real64) sa
  real(real64) sr
  real(real64) t0
  real(real64) t1
  real(real64) u
  real(real64) u0
  real(real64), intent(in) :: v
  real(real64) v0
  real(real64) va
  real(real64) vb
  real(real64) vt
  real(real64), intent(in) :: x

  pi = 3.141592653589793e+00_real64

  if ( x == 0.0e+00_real64 ) then
    if ( -1.0e+00_real64 < v .or. int ( v ) - v == 0.5e+00_real64 ) then
      hv = 0.0e+00_real64
    else if ( v < -1.0e+00_real64 ) then
      hv = ( -1.0e+00_real64 ) ** ( int ( 0.5e+00_real64 - v ) - 1 ) * 1.0e+300_real64
    else if ( v == -1.0e+00_real64 ) then
      hv = 2.0e+00_real64 / pi
    end if
    return
  end if

  if ( x <= 20.0e+00_real64 ) then

    v0 = v + 1.5e+00_real64
    call gamma ( v0, ga )
    s = 2.0e+00_real64 / ( sqrt ( pi ) * ga )
    r1 = 1.0e+00_real64

    do k = 1, 100
      va = k + 1.5e+00_real64
      call gamma ( va, ga )
      vb = v + k + 1.5e+00_real64
      call gamma ( vb, gb )
      r1 = -r1 * ( 0.5e+00_real64 * x ) ** 2
      r2 = r1 / ( ga * gb )
      s = s + r2
      if ( abs ( r2 ) < abs ( s ) * 1.0e-12_real64 ) then
        exit
      end if
    end do

    hv = ( 0.5e+00_real64 * x ) ** ( v + 1.0e+00_real64 ) * s

  else

    sa = ( 0.5e+00_real64 * x ) ** ( v - 1.0e+00_real64 ) / pi
    v0 = v + 0.5e+00_real64
    call gamma ( v0, ga )
    s = sqrt ( pi ) / ga
    r1 = 1.0e+00_real64

    do k = 1, 12
      va = k + 0.5e+00_real64
      call gamma ( va, ga )
      vb = - k + v + 0.5e+00_real64
      call gamma ( vb, gb )
      r1 = r1 / ( 0.5e+00_real64 * x ) ** 2
      s = s + r1 * ga / gb
    end do

    s0 = sa * s
    u = abs ( v )
    n = int ( u )
    u0 = u - n

    do l = 0, 1

      vt = 4.0e+00_real64 * ( u0 + l ) ** 2
      r1 = 1.0e+00_real64
      pu1 = 1.0e+00_real64
      do k = 1, 12
        r1 = -0.0078125e+00_real64 * r1 &
          * ( vt - ( 4.0e+00_real64 * k - 3.0e+00_real64 ) ** 2 ) &
          * ( vt - ( 4.0e+00_real64 * k - 1.0e+00_real64 ) ** 2 ) &
          / ( ( 2.0e+00_real64 * k - 1.0e+00_real64 ) * k * x * x )
        pu1 = pu1 + r1
      end do

      qu1 = 1.0e+00_real64
      r2 = 1.0e+00_real64
      do k = 1, 12
        r2 = -0.0078125e+00_real64 * r2 &
          * ( vt - ( 4.0e+00_real64 * k - 1.0e+00_real64 ) ** 2 ) &
          * ( vt - ( 4.0e+00_real64 * k + 1.0e+00_real64 ) ** 2 ) &
          / ( ( 2.0e+00_real64 * k + 1.0e+00_real64 ) * k * x * x )
        qu1 = qu1 + r2
      end do
      qu1 = 0.125e+00_real64 * ( vt - 1.0e+00_real64 ) / x * qu1

      if ( l == 0 ) then
        pu0 = pu1
        qu0 = qu1
      end if

    end do

    t0 = x - ( 0.5e+00_real64 * u0 + 0.25e+00_real64 ) * pi
    t1 = x - ( 0.5e+00_real64 * u0 + 0.75e+00_real64 ) * pi
    sr = sqrt ( 2.0e+00_real64 / ( pi * x ) )
    by0 = sr * ( pu0 * sin ( t0 ) + qu0 * cos ( t0 ) )
    by1 = sr * ( pu1 * sin ( t1 ) + qu1 * cos ( t1 ) )
    bf0 = by0
    bf1 = by1
    do k = 2, n
      bf = 2.0e+00_real64 * ( k - 1.0e+00_real64 + u0 ) / x * bf1 - bf0
      bf0 = bf1
      bf1 = bf
    end do

    if ( n == 0 ) then
      byv = by0
    else if ( n == 1 ) then
      byv = by1
    else
      byv = bf
    end if
    hv = byv + s0
  end if

  return
end subroutine stvhv
!> @brief subroutine stvl0.
!> @return None.
!>
!> @param x [in] Argument x.
!> @param sl0 [inout] Argument sl0.
pure subroutine stvl0 ( x, sl0 )

!*****************************************************************************80
!
!! STVL0 computes the modified Struve function L0(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    22 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) X, the argument.
!
!    Output, real(real64) SL0, the function value.
!
  implicit none

  real(real64) a0
  real(real64) a1
  real(real64) bi0
  integer(int32) k
  integer(int32) km
  real(real64) pi
  real(real64) r
  real(real64) s
  real(real64), intent(inout) :: sl0
  real(real64), intent(in) :: x

  pi = 3.141592653589793e+00_real64
  s = 1.0e+00_real64
  r = 1.0e+00_real64

  if ( x <= 20.0e+00_real64 ) then

    a0 = 2.0e+00_real64 * x / pi

    do k = 1, 60
      r = r * ( x / ( 2.0e+00_real64 * k + 1.0e+00_real64 ) ) ** 2
      s = s + r
      if ( abs ( r / s ) < 1.0e-12_real64 ) then
        exit
      end if
    end do

    sl0 = a0 * s

  else

    if ( x < 50.0e+00_real64 ) then
      km = int ( 0.5e+00_real64 * ( x + 1.0e+00_real64 ) )
    else
      km = 25
    end if

    do k = 1, km
      r = r * ( ( 2.0e+00_real64 * k - 1.0e+00_real64 ) / x ) ** 2
      s = s + r
      if ( abs ( r / s ) < 1.0e-12_real64 ) then
        exit
      end if
    end do

    a1 = exp ( x ) / sqrt ( 2.0e+00_real64 * pi * x )
    r = 1.0e+00_real64
    bi0 = 1.0e+00_real64
    do k = 1, 16
      r = 0.125e+00_real64 * r * ( 2.0e+00_real64 * k - 1.0e+00_real64 ) ** 2 / ( k * x )
      bi0 = bi0 + r
      if ( abs ( r / bi0 ) < 1.0e-12_real64 ) then
        exit
      end if
    end do

    bi0 = a1 * bi0
    sl0 = - 2.0e+00_real64 / ( pi * x ) * s + bi0

  end if

  return
end subroutine stvl0
!> @brief subroutine stvl1.
!> @return None.
!>
!> @param x [in] Argument x.
!> @param sl1 [inout] Argument sl1.
pure subroutine stvl1 ( x, sl1 )

!*****************************************************************************80
!
!! STVL1 computes the modified Struve function L1(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    05 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) X, the argument.
!
!    Output, real(real64) SL1, the function value.
!
  implicit none

  real(real64) a1
  real(real64) bi1
  integer(int32) k
  integer(int32) km
  real(real64) pi
  real(real64) r
  real(real64) s
  real(real64), intent(inout) :: sl1
  real(real64), intent(in) :: x

  pi = 3.141592653589793e+00_real64
  r = 1.0e+00_real64
  if ( x <= 20.0e+00_real64 ) then
    s = 0.0e+00_real64
    do k = 1, 60
      r = r * x * x / ( 4.0e+00_real64 * k * k - 1.0e+00_real64 )
      s = s + r
      if ( abs ( r / s ) < 1.0e-12_real64 ) then
        exit
      end if
    end do

    sl1 = 2.0e+00_real64 / pi * s

  else

    s = 1.0e+00_real64
    km = int ( 0.50e+00_real64 * x )
    km = min ( km, 25 )

    do k = 1, km
      r = r * ( 2.0e+00_real64 * k + 3.0e+00_real64 ) &
        * ( 2.0e+00_real64 * k + 1.0e+00_real64 ) / ( x * x )
      s = s + r
      if ( abs ( r / s ) < 1.0e-12_real64 ) then
        exit
      end if
    end do

    sl1 = 2.0e+00_real64 / pi * ( -1.0e+00_real64 + 1.0e+00_real64 &
      / ( x * x ) + 3.0e+00_real64 * s / x**4 )
    a1 = exp ( x ) / sqrt ( 2.0e+00_real64 * pi * x )
    r = 1.0e+00_real64
    bi1 = 1.0e+00_real64
    do k = 1, 16
      r = -0.125e+00_real64 * r &
        * ( 4.0e+00_real64 - ( 2.0e+00_real64 * k - 1.0e+00_real64 ) ** 2 ) / ( k * x )
      bi1 = bi1 + r
      if ( abs ( r / bi1 ) < 1.0e-12_real64 ) then
        exit
      end if
    end do

    sl1 = sl1 + a1 * bi1

  end if

  return
end subroutine stvl1
!> @brief subroutine stvlv.
!> @return None.
!>
!> @param v [in] Argument v.
!> @param x [in] Argument x.
!> @param slv [inout] Argument slv.
subroutine stvlv ( v, x, slv )

!*****************************************************************************80
!
!! STVLV computes the modified Struve function Lv(x) with arbitary order.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    04 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) V, the order of Lv(x).
!
!    Input, real(real64) X, the argument of Lv(x).
!
!    Output, real(real64) SLV, the value of Lv(x).
!
  implicit none

  real(real64) bf
  real(real64) bf0
  real(real64) bf1
  real(real64) biv
  real(real64) biv0
  real(real64) ga
  real(real64) gb
  integer(int32) k
  integer(int32) l
  integer(int32) n
  real(real64) pi
  real(real64) r
  real(real64) r1
  real(real64) r2
  real(real64) s
  real(real64) s0
  real(real64) sa
  real(real64), intent(inout) :: slv
  real(real64) u
  real(real64) u0
  real(real64), intent(in) :: v
  real(real64) v0
  real(real64) va
  real(real64) vb
  real(real64) vt
  real(real64), intent(in) :: x

  pi = 3.141592653589793e+00_real64

  if ( x == 0.0e+00_real64 ) then

    if ( -1.0e+00_real64 < v .or. int ( v ) - v == 0.5e+00_real64 ) then
      slv = 0.0e+00_real64
    else if ( v < -1.0e+00_real64 ) then
      slv = ( -1 ) ** ( int ( 0.5e+00_real64 - v ) - 1 ) * 1.0e+300_real64
    else if ( v == -1.0e+00_real64 ) then
      slv = 2.0e+00_real64 / pi
    end if

  else if ( x <= 40.0e+00_real64 ) then

    v0 = v + 1.5e+00_real64
    call gamma ( v0, ga )
    s = 2.0e+00_real64 / ( sqrt ( pi ) * ga )
    r1 = 1.0e+00_real64
    do k = 1, 100
      va = k + 1.5e+00_real64
      call gamma ( va, ga )
      vb = v + k + 1.5e+00_real64
      call gamma ( vb, gb )
      r1 = r1 * ( 0.5e+00_real64 * x ) ** 2
      r2 = r1 / ( ga * gb )
      s = s + r2
      if ( abs ( r2 / s ) < 1.0e-12_real64 ) then
        exit
      end if
    end do

    slv = ( 0.5e+00_real64 * x ) ** ( v + 1.0e+00_real64 ) * s

  else

    sa = -1.0e+00_real64 / pi * ( 0.5e+00_real64 * x ) ** ( v - 1.0e+00_real64 )
    v0 = v + 0.5e+00_real64
    call gamma ( v0, ga )
    s = - sqrt ( pi ) / ga
    r1 = -1.0e+00_real64
    do k = 1, 12
      va = k + 0.5e+00_real64
      call gamma ( va, ga )
      vb = - k + v + 0.5e+00_real64
      call gamma ( vb, gb )
      r1 = - r1 / ( 0.5e+00_real64 * x ) ** 2
      s = s + r1 * ga / gb
    end do
    s0 = sa * s
    u = abs ( v )
    n = int ( u )
    u0 = u - n
    do l = 0, 1
      vt = u0 + l
      r = 1.0e+00_real64
      biv = 1.0e+00_real64
      do k = 1, 16
        r = -0.125e+00_real64 * r * ( 4.0e+00_real64 * vt * vt - &
          ( 2.0e+00_real64 * k - 1.0e+00_real64 )**2 ) / ( k * x )
        biv = biv + r
        if ( abs ( r / biv ) < 1.0e-12_real64 ) then
          exit
        end if
      end do

      if ( l == 0 ) then
        biv0 = biv
      end if

    end do

    bf0 = biv0
    bf1 = biv
    do k = 2, n
      bf = - 2.0e+00_real64 * ( k - 1.0e+00_real64 + u0 ) / x * bf1 + bf0
      bf0 = bf1
      bf1 = bf
    end do

    if ( n == 0 ) then
      biv = biv0
    else if ( 1 < n ) then
      biv = bf
    end if

    slv = exp ( x ) / sqrt ( 2.0e+00_real64 * pi * x ) * biv + s0

  end if

  return
end subroutine stvlv
!> @brief subroutine timestamp.
!> @return None.
!>
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    May 31 2001   9:45:54.872_real64 AM
!
!  Modified:
!
!    31 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer(int32) d
  character ( len = 8 ) date
  integer(int32) h
  integer(int32) m
  integer(int32) mm
  character ( len = 9 ), parameter, dimension(12) :: month = [&
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ']
  integer(int32) n
  integer(int32) s
  character ( len = 10 ) time
  integer(int32) values(8)
  integer(int32) y
  character ( len = 5 ) zone

  call date_and_time ( date, time, zone, values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( error_unit, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end subroutine timestamp
!> @brief subroutine vvla.
!> @return None.
!>
!> @param va [in] Argument va.
!> @param x [in] Argument x.
!> @param pv [inout] Argument pv.
subroutine vvla ( va, x, pv )

!*****************************************************************************80
!
!! VVLA computes parabolic cylinder function Vv(x) for large arguments.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    04 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) X, the argument.
!
!    Input, real(real64) VA, the order nu.
!
!    Output, real(real64) PV, the value of V(nu,x).
!
  implicit none

  real(real64) a0
  real(real64) dsl
  real(real64) eps
  real(real64) gl
  integer(int32) k
  real(real64) pdl
  real(real64) pi
  real(real64), intent(inout) :: pv
  real(real64) qe
  real(real64) r
  real(real64), intent(in) :: va
  real(real64), intent(in) :: x
  real(real64) x1

  pi = 3.141592653589793e+00_real64
  eps = 1.0e-12_real64
  qe = exp ( 0.25e+00_real64 * x * x )
  a0 = abs ( x ) ** ( -va - 1.0e+00_real64 ) * sqrt ( 2.0e+00_real64 / pi ) * qe

  r = 1.0e+00_real64
  pv = 1.0e+00_real64
  do k = 1, 18
    r = 0.5e+00_real64 * r * ( 2.0e+00_real64 * k + va - 1.0e+00_real64 ) &
      * ( 2.0e+00_real64 * k + va ) / ( k * x * x )
    pv = pv + r
    if ( abs ( r / pv ) < eps ) then
      exit
    end if
  end do

  pv = a0 * pv

  if ( x < 0.0e+00_real64 ) then
    x1 = -x
    call dvla ( va, x1, pdl )
    call gamma ( -va, gl )
    dsl = sin ( pi * va ) * sin ( pi * va )
    pv = dsl * gl / pi * pdl - cos ( pi * va ) * pv
  end if

  return
end subroutine vvla
!> @brief subroutine vvsa.
!> @return None.
!>
!> @param va [inout] Argument va.
!> @param x [in] Argument x.
!> @param pv [inout] Argument pv.
subroutine vvsa ( va, x, pv )

!*****************************************************************************80
!
!! VVSA computes parabolic cylinder function V(nu,x) for small arguments.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    04 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real(real64) X, the argument.
!
!    Input, real(real64) VA, the order nu.
!
!    Output, real(real64) PV, the value of V(nu,x).
!
  implicit none

  real(real64) a0
  real(real64) ep
  real(real64) eps
  real(real64) fac
  real(real64) g1
  real(real64) ga0
  real(real64) gm
  real(real64) gw
  integer(int32) m
  real(real64) pi
  real(real64), intent(inout) :: pv
  real(real64) r
  real(real64) r1
  real(real64) sq2
  real(real64) sv
  real(real64) sv0
  real(real64) v1
  real(real64), intent(inout) :: va
  real(real64) va0
  real(real64) vb0
  real(real64) vm
  real(real64), intent(in) :: x

  eps = 1.0e-15_real64
  pi = 3.141592653589793e+00_real64
  ep = exp ( -0.25e+00_real64 * x * x )
  va0 = 1.0e+00_real64 + 0.5e+00_real64 * va

  if ( x == 0.0e+00_real64 ) then

    if ( ( va0 <= 0.0e+00_real64 .and. va0 == int ( va0 ) ) .or. &
      va == 0.0e+00_real64 ) then
      pv = 0.0e+00_real64
    else
      vb0 = -0.5e+00_real64 * va
      sv0 = sin ( va0 * pi )
      call gamma ( va0, ga0 )
      pv = 2.0e+00_real64 ** vb0 * sv0 / ga0
    end if

  else

    sq2 = sqrt ( 2.0e+00_real64 )
    a0 = 2.0e+00_real64 ** ( -0.5e+00_real64 * va ) * ep / ( 2.0e+00_real64 * pi )
    sv = sin ( - ( va + 0.5e+00_real64 ) * pi )
    v1 = -0.5e+00_real64 * va
    call gamma ( v1, g1 )
    pv = ( sv + 1.0e+00_real64 ) * g1
    r = 1.0e+00_real64
    fac = 1.0e+00_real64

    do m = 1, 250
      vm = 0.5e+00_real64 * ( m - va )
      call gamma ( vm, gm )
      r = r * sq2 * x / m
      fac = - fac
      gw = fac * sv + 1.0e+00_real64
      r1 = gw * r * gm
      pv = pv + r1
      if ( abs ( r1 / pv ) < eps .and. gw /= 0.0e+00_real64 ) then
        exit
      end if
    end do

    pv = a0 * pv

  end if

  return
end subroutine vvsa

end module mod_special_functions
