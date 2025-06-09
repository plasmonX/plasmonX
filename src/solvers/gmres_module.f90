!> GMRES module
!!
!! This module contains the subroutines for complex GMRES algorithm
!!
!! Date         : 2025
!!
module gmres_module

   use parameters_module
   use target_module 
   use output_module
   use memory_manager_module

   implicit none

contains 

   !> Solve a linear system A x = b using restarted GMRES
   !! Uses either explicit matrix or on-the-fly matrix-vector product.
   !! Performs Modified Gram-Schmidt orthogonalization and applies Givens rotations.
   !!
   !! Input  : i_freq            -- frequency index, required if onthefly = .true.
   !! Input  : n                 -- system size
   !! In/Out : x                 -- solution vector (initial guess assumed 0)
   !! Input  : matrix_iterative  -- system matrix (used only if onthefly = .false.)
   !! Input  : onthefly          -- if true, use target_%product_matrix_vector
   !! Input  : tol               -- convergence tolerance on relative residual
   !! Input  : max_iter          -- maximum number of outer iterations
   !! Input  : restart           -- restart parameter (dimension of Krylov space)
   subroutine gmres_solver(i_freq, n, x, matrix_iterative, onthefly,        &
                           tol, max_iter, restart)

      implicit none

      !input/output variables
      integer, intent(in)                        :: i_freq   ! index of frequency
      integer, intent(in)                        :: n        ! system size
      complex(dp), dimension(n), intent(inout)   :: x        ! solution array
      !matrix left hand side (optional -- if on memory)
      complex(dp), dimension(n,n), optional, intent(in) :: matrix_iterative 
      !on the fly algorithm
      logical, intent(in)                        :: onthefly     
      real(dp), intent(in)                       :: tol      ! threshold
      integer, intent(in)                        :: max_iter ! max GMRES cycles
      integer, intent(in)                        :: restart  ! restart dimension

      !internal variables
      complex(dp), dimension(:,:), allocatable :: V, H
      complex(dp), dimension(:), allocatable   :: w, r, s, y, tmp
      real(dp), dimension(:), allocatable      :: cs
      complex(dp)              :: zero_cp, one_cp
      complex(dp)              :: zdotc
      complex(dp)              :: tmp_Hjj, tmp_Hj1j
      real(dp)                 :: beta
      real(dp)                 :: normb, resid
      real(dp)                 :: dznrm2
      real(dp)                 :: temp
      integer                  :: i, j, m, iter
      logical                  :: converged

      zero_cp = dcmplx(zero, zero)
      one_cp  = dcmplx(one, zero)
      m = min(restart, n)

      call mem_man%alloc(V, n, m+1, "v_gmres")
      call mem_man%alloc(H, m+1, m+1, "h_gmres")
      call mem_man%alloc(w, n,   "w_gmres") ! tmp sol
      call mem_man%alloc(r, n,   "r_gmres") ! tmp rhs
      call mem_man%alloc(s, m+1, "s_gmres")
      call mem_man%alloc(y, m,   "y_gmres")
      call mem_man%alloc(cs, m,  "cs_gmres")

      !write r in the right hand side
      call zcopy(n, x, 1, r, 1)

      !initialization
      x = zero_cp
      iter = 0
      converged = .false.

      !let us calculate the norm of b
      normb = dznrm2(n, r, 1)

      !the solution is zero
      if (normb .eq. zero) return
   
      do while (iter < max_iter .and. .not. converged)
   
         beta = normb

         V(:,1) = zero_cp
         call zaxpy(n, dcmplx(one/beta,zero), r, 1, V(:,1), 1)

         s    = zero_cp
         cs   = zero
         H    = zero_cp
         temp = zero

         H(1,m+1) = beta
         do i = 1, m
            H(i+1,m+1) = zero
         enddo

         write(out_%iunit,'(6x,a,5x,a, 9x)') "Iteration", "Residual"
   
         do j = 1, m

            if (onthefly) then
               call target_%product_matrix_vector(i_freq, V(:,j), w)
            else
               call zgemv('N', n, n, &
                          one_cp,    &
                          matrix_iterative, &
                          n,                &
                          V(:,j),           &
                          1,                &
                          zero_cp,          &
                          w,                &
                          1)
            endif
            call target_%new_gmres_diagonal_shift(i_freq,V(:,j),w)
   
            if (j.lt.25) then
               do i = 1, j
                  H(i,j) = zdotc(n, V(1,i),1,w,1)
                  call zaxpy(n, -H(i,j), V(:,i), 1, w, 1)
               end do
            else
               !H(1:j, j) = V(:,1:j)^H * w
               call zgemv('C', n, j, one_cp, V, n, w, 1, zero_cp, H(1,j), 1)
               !w = w - V(:,1:j) * H(1:j, j)
               call zgemv('N', n, j, -one_cp, V, n, H(1,j), 1, one_cp, w, 1)
               
               !Second projection to remove numerical errors
               call mem_man%alloc(tmp, j, "tmp_gmres")
               call zgemv('C', n, j, one_cp, V, n, w, 1, zero_cp, tmp, 1)
               call zgemv('N', n, j, -one_cp, V, n, tmp, 1, one_cp, w, 1)
               do i = 1, j
                  H(i,j) = H(i,j) + tmp(i)
               end do            
               call mem_man%dealloc(tmp, "tmp_gmres")
            endif
 
            H(j+1,j) = dznrm2(n, w, 1)


            if (H(j+1,j) .ne. zero) then
               V(:,j+1) = zero_cp
               call zaxpy(n, one_cp / H(j+1,j), w, 1, V(:,j+1), 1)
            endif
   
            ! Apply Givens rotations to H(1:j+1,j)
            do i = 1, j-1
               call zrot(1, H(i,j), 1, H(i+1,j), 1, cs(i), s(i))
            end do
   
            tmp_Hjj = H(j,j)
            tmp_Hj1j= H(j+1,j)

            ! Compute i-th rotation matrix
            call zrotg(tmp_Hjj, tmp_Hj1j, cs(j), s(j))
   
            ! Apply rotation to rhs s
            call zrot(1, H(j,m+1), 1, H(j+1,m+1), 1, cs(j), s(j))
            call zrot(1, H(j,j), 1, H(j+1,j), 1, cs(j), s(j))

            H(j+1, j) = zero_cp

            resid = abs(H(j+1,m+1)) 

            write(out_%iunit,'(1x,i10,5x,e13.5)') j, resid
            flush(out_%iunit)

            if (resid < tol) then

               ! Solve upper triangular system H(1:m,1:m)*y = s(1:m)
               call zcopy(j, H(1,m+1), 1, y, 1)
               call ztrsv('U', 'N', 'N', j, H, m+1, y, 1)
               ! x = x + V(:,1:m)*y
               call zgemv('N', n, j, one_cp, V, n, y, 1, one_cp, x, 1)

               ! Check final residual
               if (onthefly) then
                  call target_%product_matrix_vector(i_freq, x, w)
               else
                  call zgemv('N', n, n,        &
                             one_cp,           &
                             matrix_iterative, &
                             n,                &
                             x,                &
                             1,                &
                             zero_cp,          &
                             w,                &
                             1)
               endif
               call target_%new_gmres_diagonal_shift(i_freq,x,w)

               !final check residual
               !$omp parallel do 
               do i=1,n
                  w(i) = r(i) - w(i)
               end do
               !$omp end parallel do
               resid = dznrm2(n, w, 1) 

               !set converged
               if(resid .lt. tol) then
                  converged = .true.
                  exit
               endif
            endif
         end do
   
         if (converged) exit
         iter = iter + 1
         call out_%error("GMRES not converged and restart GMRES NYI.")

      end do
   
      call mem_man%dealloc(V, "v_gmres")
      call mem_man%dealloc(H, "h_gmres")
      call mem_man%dealloc(w, "w_gmres")
      call mem_man%dealloc(r, "r_gmres")
      call mem_man%dealloc(s, "s_gmres")
      call mem_man%dealloc(y, "y_gmres")
      call mem_man%dealloc(cs,"cs_gmres")

   end subroutine gmres_solver

end module gmres_module

