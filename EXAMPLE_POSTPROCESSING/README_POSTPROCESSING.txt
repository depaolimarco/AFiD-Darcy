We provide here few examples of post-processing routines (Matlab file postprocessing.m).
To perform the analysis, set to 1 the flag corresponding as below:
1) Visualization of global quantities computed runtime (flag i_Nusselt).
	- Visualization of residuals of first- and second-order energy budgets (dimensionless form in convective units)
	- For a derivation of the second order energy budget, see for instance Eq. (18) in [1] and references therein.
2) Visualization of horizontally-avergaed profiles
	- time-averaged profiles (flag i_profiles_avg) 
	- instantaneous profiles (flag i_profiles_ins)
	- note the input parameters 601-607 to control output of the statistics and profiles
3) Visualization of intantaneous scalar fields
	- it produces a movie saved in the folder 'movies' (flag i_movie_1)
        - note the input parameters 601-607 to control output of the fields

To perfom the analysis, run the file from the current folder and set the desired flag(s).

References
[1] De Paoli, M. (2023). Convective mixing in porous media: a review of Darcy, pore-scale and Hele-Shaw studies. The European Physical Journal E, 46(12), 129.
