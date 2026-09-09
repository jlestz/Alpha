--- /mnt/user-data/uploads/Alpha_transport.f90
+++ Alpha_transport.f90
@@ -84,6 +84,18 @@
 !  every real in this file is declared plain "real" under implicit none;
 !  in single precision this tolerance would never be reached.
      real(dp), parameter :: error_tol = 1.0E-12_dp
+!  JBL: print stride for the per-iteration scalar error diagnostics
+!  written to Alpha_transport.out (main-loop, alpha2, and He-loop
+!  write statements below). Writing every iteration was blowing up
+!  the log file at 1E5+ iterations; this cuts it to about
+!  n_up_loop/n_print_stride lines. ii=1 and ii=n_up_loop are always
+!  printed regardless of stride so the first and last regular-loop
+!  values are never missing. A converged early exit also always gets
+!  one more regular-format line at the exit iteration (in addition to
+!  the separate "Converged at ii=..." summary message), so a run that
+!  converges off the n_print_stride grid still ends with a normal
+!  data point, not just the summary line.
+     integer, parameter :: n_print_stride = 100
      real(dp) :: relax
      real(dp) :: relax_f
      real(dp) :: thfrac
@@ -1564,7 +1576,11 @@
   
  
 !  write(3,*) 'ii=',ii,'D_TAE=',D_TAE,'error=',error,error_rho(2),error_rho(25)
-   write(3,*) 'ii=',ii,'D_TAE=',D_TAE,'error=',error,error_f,error_f_rho(25)
+!  JBL: throttle to every n_print_stride iterations, always keeping
+!  the first and last (n_up_loop) regular-loop values.
+   if ((ii .eq. 1) .or. (mod(ii-1,n_print_stride) .eq. 0) .or. (ii .eq. n_up_loop)) then
+     write(3,*) 'ii=',ii,'D_TAE=',D_TAE,'error=',error,error_f,error_f_rho(25)
+   endif
 
 ! JBL 08/31/26 debug 
    if ((l_debug_plots .eq. 1) .and. (mod(ii,100) .eq. 0)) then
@@ -1724,7 +1740,11 @@
 
 
 
-  write(3,*) 'ii=',ii,'error2=',error2,error2_rho(1),error2_rho(25)
+!  JBL: throttle to every n_print_stride iterations, always keeping
+!  the first and last (n_up_loop) regular-loop values.
+  if ((ii .eq. 1) .or. (mod(ii-1,n_print_stride) .eq. 0) .or. (ii .eq. n_up_loop)) then
+    write(3,*) 'ii=',ii,'error2=',error2,error2_rho(1),error2_rho(25)
+  endif
 
    do i=1,n_rho_grid-1
     n_alpha2_tran_rho(i) = relax*n_alpha2_tran_rho(i)+(1.-relax)*n_alpha2_tran_p_rho(i)
@@ -1744,11 +1764,21 @@
 !  the exit.
    if (NBI_flag .eq. 2) then
      if ((error_check .lt. error_tol) .and. (error2 .lt. error_tol)) then
+!  JBL: also emit the regular per-iteration diagnostic lines on a
+!  converged exit, in the same format/columns as the throttled writes
+!  above, so the log always has a final regular-format data point
+!  even when convergence lands off the n_print_stride grid. The
+!  separate "Converged at ii=..." summary line below is kept as well.
+       write(3,*) 'ii=',ii,'D_TAE=',D_TAE,'error=',error,error_f,error_f_rho(25)
+       write(3,*) 'ii=',ii,'error2=',error2,error2_rho(1),error2_rho(25)
        write(3,*) 'Converged at ii=',ii,' error_check=',error_check,' error2=',error2
        exit
      endif
    else
      if (error_check .lt. error_tol) then
+!  JBL: same as above -- emit the regular-format diagnostic line on a
+!  converged exit, in addition to the summary "Converged at ii=..." line.
+       write(3,*) 'ii=',ii,'D_TAE=',D_TAE,'error=',error,error_f,error_f_rho(25)
        write(3,*) 'Converged at ii=',ii,' error_check=',error_check
        exit
      endif
@@ -2175,7 +2205,11 @@
 
 
 
-  write(3,*) 'ii=',ii,'error=',error,error_rho(1),error_rho(25)
+!  JBL: throttle to every n_print_stride iterations, always keeping
+!  the first and last (n_up_loop) regular-loop values.
+  if ((ii .eq. 1) .or. (mod(ii-1,n_print_stride) .eq. 0) .or. (ii .eq. n_up_loop)) then
+    write(3,*) 'ii=',ii,'error=',error,error_rho(1),error_rho(25)
+  endif
 
    do i=1,n_rho_grid-1
     n_He_tran_rho(i) = relax*n_He_tran_rho(i)+(1.-relax)*n_He_tran_p_rho(i)
@@ -2185,6 +2219,10 @@
 !  above. Placed after the relax step so the exit leaves the same
 !  fully-relaxed state a normal completed iteration would.
    if (error .lt. error_tol) then
+!  JBL: same pattern as the main loop -- emit the regular-format
+!  diagnostic line on a converged exit, in addition to the summary
+!  "He loop converged at ii=..." line.
+     write(3,*) 'ii=',ii,'error=',error,error_rho(1),error_rho(25)
      write(3,*) 'He loop converged at ii=',ii,' error=',error
      exit
    endif
