Damiens version of the SPARC code, originally written by Shravan Hanasoge. Upgraded with a lot of debugging + help from Jamon Pennicot.
Used for the study of linear perturbations in the solar atmosphere.

Damiens to-do list:
0) Quiet simulations currently read in a x=0,y=0 of a 3D array, it should read a 1D array. 
1) Re-do formatting in the initialize routine that prints simulation setup to be neater and more readable.
2) Neaten + comment. Remove as much as possible from the 'common' variables in all_modules.f90. Use INTENT(IN), INTENT(OUT) where possible to reduce chances of errors. Add public/private to the modules. Remove ifs from for loops where possible. When are nested loops quicker than array operations??
3) Put PMLs in there own module, they make the RHS update routines unreadable.
4) Move driver initialization into init_driver in initialization.f90
5) Move all IO into its own module.
6) Add a few more options to bkg + ptb saves
7) Displacement PMLs not working well, requires filtering EVERY timestep, there is probably a better way of formulating or implementing them.
8) move random-sources to use HDF5, remove fits from code.

version 1.6 03/10/2017

- Fixed the derivatives routine (I think)

Version 1.5 26/07/2016

- Flows added to Velocity, tested in both displacement and velocity
- Fixed and tested random sources.

Version 1.4 26/07/2016

- Fixed small bug with filtering
- Restart now working
- Max VA can now be set in control file, reduction is not adjusted in all_modules.f90
- Save Chunk can now save a 3D subset of the total array.

Version 1.3 07/06/2016

- Jamon added a moving source for his 'Sunquake' runs.
- Randomised pulses fixed, still to slow.
- Added read-write PMLs, fixed my HDF5 routines a bit.
- Fixed a number of random typos/bugs.
- Changed how parameters file has been read in.

Version 1.2 - 19/03/2016

- Jamon has ironed most of the bugs out of 2D. Displacement mode is still largely untested.
- Damping is in. Still not working in 2D.
- Tweaked the derivative routines to make them more readable and slightly faster. Will perform similar tweak to filtering routines one day.
- using simple 5-point differences in xy, seems to fixed some of the weird artifacts. The large stencils maybe don't work with the low ppw used?

Version 1.1 - 24/02/2016

- HDF5 Files are in, fits is out
- 2D is now 2.5D Jamon is testing + fixing the bugs.
- The wave damping is in progress, but I need to learn to use fftw.

Version 1.0

I have made a variety of changes, which can be summarised as:

1) New .trol file which allows parameter changes without recompilation.
2) New explicit derivative and digital filter scheme, although currently an 11pt, this
could easily be reduced to a 7pt for further speed increases.
3) Domain decomposition is now done in 3-dimensions. So a reasonable level of scaling is
achieved to far more processors. Communication is done in the derivative routines and a (lazy)
attempt has been made at masking it with computation.
4) The modules have been re-arranged slightly. They now behave as follows:

- Driver.f90: as before, begins the time-step and calls the correct routines.
- All_modules.f90: Contains the important arrays + variables. Also contains all I/O routines + some other stuff, like calculating norms,
interpolation for drivers etc.
- initialisation.f90 : Contains all initialisation: i) MPI, ii) Array Allocation + initialization, iii) init the RHS, iv) init the stepping variables, v) init the derivative + filters.
- step.f90 : Performs the Runge-Kutta time step integration.
- RHS.f90 : Performs the RHS updates for the linear MHD equations. Current contains a MHD + Quiet (Hydro only) Routines. Displacement + Flows versions of this will be added as I have time.
- RHS2D.f90: As above but for 2D, probably doesn't work at the moment, I don't know, I don't use it. Will fix one day.
- derivatives.f90: All derivative and filtering routines.

5) There is no longer a 'background' file. All that is read in is: 1-dimensional arrays: zaxis (in Mm), grav (cgs). 3D-arrays: pressure, density, sound speed, if magnetic: B_x, B_y, B_z. I can add an option for this back in if it needed.
6) The horizontal digital filters will kill all low-frequency oscillations, this means your nu-omega ridge diagrams will be truncated at a l-value depending on your grid resoultion. However these filters
are far more stable when dealing with large density + pressure gradients in the horizontal direction (i.e. a sunspot model with large wilson depression). I currently filter z-direction every 10 second, and xy-direction
every 30 seconds. But you will want to adjust this based on your time-step.
7) The PSR_pulse stuff is a random-ish pulse driver I am working on. It currently does not work very well.
8) Horizontal PMLS are in, the PMLs have been recoded in terms of N,RC & kappa as in the paper, so they can be tuned. 
9) HDF5 Files are in and working.
