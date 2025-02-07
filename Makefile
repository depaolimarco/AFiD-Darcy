# Choose the machine being used
# Options: PC, SNELLIUS, IRENE, MARENOSTRUM, SUPERMUC, DISCOVERER
MACHINE=PC
FLAVOUR=GNU
# Modules required for each HPC system as follows:
# SNELLIUS:
#	GNU: 2022 foss/2022a HDF5/1.12.2-gompi-2022a
# 	Intel: 2022 intel/2022a FFTW/3.3.10-GCC-11.3.0 HDF5/1.12.2-iimpi-2021a
# IRENE (Intel): flavor/hdf5/parallel hdf5 fftw3/gnu
# MARENOSTRUM (Intel): fabric intel mkl impi hdf5 fftw szip
# SUPERMUC (Intel): fftw hdf5
# DISCOVERER:
#	GNU: hdf5/1/1.14/latest-gcc-openmpi fftw/3/latest-gcc-openmpi lapack
#	Intel: hdf5/1/1.14/latest-intel-openmpi fftw/3/latest-gcc-openmpi mkl

#=======================================================================
#  Compiler options
#=======================================================================

# Object and module directory:
OBJDIR=obj

ifeq ($(FLAVOUR),GNU)
	FC = h5pfc -cpp -ffree-line-length-none -fdefault-real-8 -fdefault-double-8 -fallow-argument-mismatch
	FC += -O3
else
	FC = h5pfc -fpp -r8
endif

ifeq ($(MACHINE),PC)
# GNU Debug Flags
	# FC += -O0 -g -fbacktrace -Wall -Wextra
	# FC += -Wpedantic
	# FC += -Warray-temporaries
	# FC += -fcheck=all -finit-real=snan -ffpe-trap=invalid #-std=f2018
# FC += -O0 -pg -fbacktrace -fbounds-check
# Intel Debug Flags
# FC += -O0 -g -traceback -check bounds
	ifeq ($(FLAVOUR),GNU)
		# LDFLAGS = -L/opt/homebrew/lib/ -lfftw3 -llapack -lhdf5_fortran -lhdf5 -lz -ldl -lm
		LDFLAGS = -L/usr/local/lib/ -lfftw3 -llapack -lhdf5_fortran -lhdf5 -lz -ldl -lm
	else
		LDFLAGS = -lfftw3 -qmkl=sequential
	endif
endif

ifeq ($(MACHINE),DISCOVERER)
	ifeq ($(FLAVOUR),GNU)
		LDFLAGS += -lfftw3 -llapack -ldl
	else
		LDFLAGS += -lfftw3 -qmkl=sequential
	endif
endif
ifeq ($(MACHINE),SNELLIUS)
	ifeq ($(FLAVOUR),GNU)
		FC += -O2 -march=znver1 -mtune=znver1 -mfma -mavx2 -m3dnow -fomit-frame-pointer
		LDFLAGS = -lfftw3 -lopenblas -ldl
	else
		FC += -align array64byte -fma -ftz -fomit-frame-pointer
		LDFLAGS = -lfftw3 -qmkl=sequential
	endif
endif
ifeq ($(MACHINE),IRENE)
	FC += -mtune=skylake -xCORE-AVX512 -m64 -fPIC $(FFTW3_FFLAGS)
	LDFLAGS = $(FFTW3_LDFLAGS) $(MKL_LDFLAGS) -ldl
endif
ifeq ($(MACHINE),MARENOSTRUM)
	FC += -mtune=skylake -xCORE-AVX512 -m64 -fPIC $(FFTW_FFLAGS)
	LDFLAGS = $(FFTW_LIBS) -mkl=sequential
endif
ifeq ($(MACHINE),SUPERMUC)
	FC = mpif90 -fpp -r8 -O3 $(HDF5_INC)
	LDFLAGS = $(FFTW_LIB) $(HDF5_F90_SHLIB) $(HDF5_SHLIB) -qmkl=sequential
endif

ifeq ($(FLAVOUR),GNU)
	FC += -J $(OBJDIR)
else
	FC += -module $(OBJDIR)
endif

#===============================================================================
# 2. Make Rules
#===============================================================================
#-------------------------------------------------------------------------------
# Non-module Fortran files to be compiled:
#-------------------------------------------------------------------------------
EXTRA_DIST = 	2DECOMP/transpose_z_to_x.F90 2DECOMP/transpose_x_to_z.F90 2DECOMP/transpose_x_to_y.F90\
				2DECOMP/transpose_y_to_x.F90 2DECOMP/transpose_y_to_z.F90 2DECOMP/transpose_z_to_y.F90\
				2DECOMP/factor.F90 2DECOMP/halo.F90 2DECOMP/fft_common.F90 2DECOMP/alloc.F90 2DECOMP/halo_common.F90

FFILES += 	CalcMaxCFL.F90 CalcOutput.F90\
			CalcLocalDivergence.F90 CalcFluxes.F90\
			CheckDivergence.F90 CorrectPressure.F90 CorrectVelocity.F90\
			CreateGrid.F90 CreateResultDirectories.F90 CreateInitialConditions.F90\
			DeallocateVariables.F90 ExplicitTermsTemp.F90 \
			Factorize.F90 GlobalQuantities.F90 HdfReadContinua.F90 HdfRoutines.F90\
			Interp.F90 ImplicitAndUpdateTemp.F90\
			InitPressureSolver.F90 InitTimeMarchScheme.F90 InitVariables.F90\
			LocateLargeDivergence.F90 MakeMovieCut.F90 MpiAuxRoutines.F90 QuitRoutine.F90\
			ReadFlowField.F90 ReadInputFile.F90 ResetLogs.F90 SetVelocityField.F90 SlabDumpRoutines.F90\
			SolveImpEqnUpdate_Temp.F90 SolvePressureCorrection.F90\
			StatReadReduceWrite.F90 StatRoutines.F90 TimeMarcher.F90\
			WriteFlowField.F90 WriteFlowFieldSnapshot.F90 WriteGridInfo.F90\
			Main.F90

#-------------------------------------------------------------------------------
# Files that create modules:
#-------------------------------------------------------------------------------
MFILES = Param.F90 decomp_2d.F90 AuxiliaryRoutines.F90 decomp_2d_fft.F90

# Object and module directory:
SRCDIR=src
OBJDIR=obj
OBJS  := $(FFILES:%.F90=$(OBJDIR)/%.o)
MOBJS := $(MFILES:%.F90=$(OBJDIR)/%.o)

#-------------------------------------------------------------------------------
# Make program 
#-------------------------------------------------------------------------------
PROGRAM = afid 

#Compiling 
all: directories $(PROGRAM) 
$(PROGRAM): $(MOBJS) $(OBJS) 
	$(FC) -o $@ $^ $(LDFLAGS) 

#-------------------------------------------------------------------------------
# Dependencies 
#-------------------------------------------------------------------------------
$(OBJDIR)/Param.o: $(SRCDIR)/Param.F90
	$(FC) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/decomp_2d.o: $(SRCDIR)/2DECOMP/decomp_2d.F90
	$(FC) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/AuxiliaryRoutines.o: $(SRCDIR)/AuxiliaryRoutines.F90 
	$(FC) -c -o $@ $< $(LDFLAGS) 
$(OBJDIR)/decomp_2d_fft.o: $(SRCDIR)/2DECOMP/decomp_2d_fft.F90
	$(FC) -c -o $@ $< $(LDFLAGS) 
$(OBJDIR)/ImplicitDecomp.o: $(SRCDIR)/ImplicitDecomp.F90 
	$(FC) -c -o $@ $< $(LDFLAGS)

$(OBJDIR)/%.o: $(SRCDIR)/%.F90 $(MOBJS)
	$(FC) -c -o $@ $< $(LDFLAGS)

#-------------------------------------------------------------------------------
# Clean up 
#-------------------------------------------------------------------------------
clean: 
	rm -rf $(OBJDIR)/*.o $(OBJDIR)/*.mod $(OBJDIR)/*genmod* $(OBJDIR)/*.o obj\
	rm afid

.PHONY: directories
directories: $(OBJDIR) 
$(OBJDIR): 
	mkdir -p ${OBJDIR} 
