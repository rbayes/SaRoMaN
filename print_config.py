import os

class print_config:

	def __init__(self):
		#Mind geometry
		self.MIND_type = 0   # Cylinder
		self.MIND_xdim = 6.0 # m
		self.MIND_ydim = 6.0 # m
		self.MIND_zdim = 13.0 # m
		self.MIND_vertex_xdim = 2.0 # m
		self.MIND_vertex_ydim = 2.0 # m
		self.MIND_vertex_zdim = 2.0 # m
		self.MIND_ear_xdim  = 0.4393 # m
		self.MIND_ear_ydim = 2.8994 # m
		self.MIND_bore_diameter = 0.2 # m
		#Mind internal dimensions
		self.MIND_active_mat = 'G4_POLYSTYRENE'
		self.MIND_width_active = 1.5 # cm
		self.MIND_rad_length_active = 413.1 #mm
		self.MIND_active_layers = 1
		self.MIND_passive_mat = 'G4_Fe'
		self.MIND_width_passive = 1.5 # cm
		self.MIND_rad_length_passive = 17.58 #mm
		self.MIND_width_bracing = 0.1 # cm Aluminum
		self.MIND_width_air = 0.25 # cm
		self.MIND_rad_length_air = 303.9 #mm

	def print_file(self,filename,data):
	    '''
	    Print data to file filename.
	    '''
	    outfile = open(filename,'w+')
	    outfile.write(data)
	    outfile.close()

	def print_digi_config(self,dictionary,filename):
		filedata = '''
########################################################################
#                                                                      #
#  This is a parameter file that can be read by bhep sreader class.    #
#                                                                      #
########################################################################

#############################################
#  parameters for the run type
#############################################

# number of events to be processed.
RUN nEvents I %(Nevts)d

# gausian sigma for smear (cm)
RUN Gaus_Sigma D 1.0

# energy smear (%%)
RUN Eng_Res D 0.11

# seed value for random generator
CON Gen_seed D 107311191

#############################################
# Parameters for hit_constructor
#############################################

# MIND dimensions (m)
CON MIND_x D %(MIND_xdim)d
CON MIND_y D %(MIND_ydim)d
CON MIND_z D %(MIND_zdim)d
CON VERTEX_x D %(MIND_vertex_xdim)d
CON VERTEX_y D %(MIND_vertex_ydim)d
CON VERTEX_z D %(MIND_vertex_zdim)d

# MIND internal dimensions (cm)
# ****IRON****
CON widthI D 1.5

# ****SCIN**** (cm)
CON widthS D 1.5
CON nplane I 1

# **** Aluminum **** (cm)
CON widthAl D 0.1

CON isOctagonal I 0

# ****AIR **** (cm)
CON widthA D 0.25
# radiation length (mm)
CON x0Sc D 413.1

CON rec_boxX D 2.0
CON rec_boxY D 2.0

# minimum energy at plane to be detected.(MeV)
CON min_eng D 0.000016

#############################################
#  Data to read and write
#############################################

DATA idst_files SV 1
%(out_base)s/G4_out/nd_%(part)s%(inttype)s/nd_%(part)s%(inttype)s_%(seed)d.dst.root

DATA odst_file S %(out_base)s/digi_out/nd_%(part)s%(inttype)s/nd_%(part)s%(inttype)s_%(seed)d_digi.dst.root
'''% dict(dictionary, **vars(self))

		self.print_file(filename,filedata)

	def print_mindG4_config(self,dictionary,filename):
		filedata = '''
### ---------------------------------------------------------------------------
###  $Id: example_mindG4.config 436 2010-11-03 15:06:50Z alaing $
###
###  Example configuration file for the mindG4 application.
### ---------------------------------------------------------------------------  


### JOB options ###########################################
JOB output_dst    S %(out_base)s/G4_out/nd_%(part)s%(inttype)s/nd_%(part)s%(inttype)s_%(seed)d.dst.root
JOB number_events I %(Nevts)s

JOB random_seed I 13243%(seed)d

### GEOMETRY configuration parameters #####################

### Detector dimensions (in mm)
GEOMETRY width  D 6000.
GEOMETRY height D 6000.
GEOMETRY length D 2000.
GEOMETRY vwidth D 2000.
GEOMETRY vheight D 2000.
GEOMETRY vdepth  D 2000.
GEOMETRY ear_width D 439.3
GEOMETRY ear_height D 2899.4
GEOMETRY bore_diameter D 200

### Calorimeter layers thicknesses (in mm)
GEOMETRY active_thickness  D 15.
GEOMETRY passive_thickness D 15.
GEOMETRY bracing_thickness D 1.

GEOMETRY IsOctagonal I 0

### Active and passive materials (see documentation for values)
GEOMETRY active_material  S G4_POLYSTYRENE
GEOMETRY passive_material S G4_Fe

### Uncomment if you want to simulate TASD.
## Sets all volumes as sensitive detectors.
## Remember to use G4_POLYSTYRENE as passive_material as well!!.
# GEOMETRY TASD I 0

### number of active layers per piece. (number of gaps should be 1 greater!)
GEOMETRY num_active_layers I 1

### Air gaps between layers (in mm)
GEOMETRY gap1 D 2.5
GEOMETRY gap2 D 2.5
GEOMETRY gap3 D 2.5
GEOMETRY gap4 D 2.5

### MAgnetic field. Still uniform vector.(T)
## A single bit to turn the uniform field off in favour of a toroidal field
GEOMETRY isUniform I 0

# if uniform then [bx,by,bz]. if toroidal then [br,btheta,bphi]
GEOMETRY field DV 3
0.
1.
0.

# If we wish to use a field map then the identity of the field map
# must be entered here.
# GEOMETRY FieldMap S /data/neutrino04/common_SW/MIND/mindG4/MIND_field_map_files-5_cm_grid/iron_field_halfplane_2.res
# GEOMETRY FieldMap S /data/neutrino04/common_SW/MIND/mindG4/MIND_field_map_files-5_cm_grid/iron_field_combined_z=0.res

GEOMETRY FieldScaling D %(Bfield)d

### GENERATION configuration parameters ################### 

### Three generators, SINGLE_PARTICLE, NUANCE and GENIE, with different
### configuration parameters, are available. The first is chosen
### as default. Comment out the following lines if you want to
### use the latter.

GENERATION generator S GENIE
# GENERATION generator S SINGLE_PARTICLE


# GENERATION particle_name S %(part)s

### Particle energy will be sample between these two values (in GeV)
# GENERATION energy_min    D 0.300
# GENERATION energy_max    D 3.000

### Uncomment the following lines if you want to use the NUANCE
### event generator.
#
### NUANCE/GENIE data files for passive and active materials!! Set the filenames here!!

GENERATION active_material_data S %(out_base)s/genie_samples/nd_%(part)s%(inttype)s/ev0_%(ASeed)s_%(pid)d_1000060120[0.922582],1000010010[0.077418]_%(Nevts)s.root
GENERATION passive_material_data S %(out_base)s/genie_samples/nd_%(part)s%(inttype)s/ev0_%(seed)d_%(pid)d_1000260560_%(Nevts)s.root
#
### Vertex location (RANDOM, ACTIVE, PASSIVE, FIXED).
GENERATION vertex_location S RANDOM

### Special simulation for training purpose
### Use the muon but not the hadronization from GENIE
# GENERATION FSL_Select I 0

# Vertex if FIXED requested.
#GENERATION fvert DV 3
#0.
#0.
#-18000
#

### PHYSICS configuration parameters ######################

### Production cut (in mm) for secondaries
PHYSICS production_cut D 30.
#
### Minimum Kinetic energy for a particle to be tracked (MeV).
PHYSICS minimum_kinEng D 100.
'''% dictionary

		self.print_file(filename,filedata)

	def print_rec_config(self,dictionary,filename):
		filedata = '''
########################################################################
#                                                                      #
#  This is a parameter file that can be read by bhep sreader class.    #
#                                                                      #
########################################################################

#############################################
#  parameters for the setup
#############################################

# MIND dimensions (m)
RUN MIND_x D 6.
RUN MIND_y D 6.
RUN MIND_z D 13.
RUN vertex_x D 2.
RUN vertex_y D 2.
RUN vertexDepth D 2.
RUN EAR_height D 2.8994
RUN EAR_width  D 0.4393

RUN IsOctagonal I 0

# MIND internal dimensions (cm)
# ****IRON****
RUN widthI D 1.5
# radiation length (mm)
RUN x0Fe D 17.58

# ****SCIN**** (cm)
RUN widthS D 1.5
RUN nplane I 1
# radiation length (mm)
RUN x0Sc D 413.1

# ****AIR **** (cm)
RUN widthA D 0.25
# radiation length (m)
RUN x0AIR D 303.9

# ****Aluminium **** (cm)
RUN widthAl D 0.1

# relative density, Sc/Fe, AIR/Sc.
RUN rel_denSI D 0.135
RUN rel_denAS D 0.001

# Attenuation in Wavelength shifting fibre
RUN WLSatten D 5000.

# measurement type.
RUN meas_type S xyz

# magnetic field.
RUN mag_field DV 3
0.
1.
0.

# RUN mag_field_map S /data/neutrino04/common_SW/MIND/mindG4/MIND_field_map_files-5_cm_grid/iron_field_halfplane_2.res
# RUN mag_field_map S /data/neutrino04/common_SW/MIND/mindG4/MIND_field_map_files-5_cm_grid/iron_field_combined_z=0.res

RUN fieldScale D %(seed)d
########

# energy loss (MeV/cm)
# 1.5 cm Fe, 2 cm Sc, 1.0 cm Air
# RUN de_dx_min D 0.5323 
# 2 cm Fe, 2 cm Sc
# RUN de_dx_min D 0.637  
# 1.5 cm Fe, 1.5 cm Sc, 0.5 cm Air
# RUN de_dx_min D 0.6484
# 1.5 cm Fe, 1.5 cm Sc, 0.5 cm Air, 0.1 cm Aluminum
RUN de_dx_min D 0.575
# 3 cm Fe, 2 cm Sc
# RUN de_dx_min D 0.757
RUN de_dx_scint D 0.205

# Position resolution for detector (cm).
RUN pos_res D 0.75

# Step size for track fitting (cm).
RUN StepSize D 1.

# name of detector for hit getter.
RUN detect S tracking

########
# For hit clustering.(edge in cm)
RUN do_clust I 1
RUN rec_boxX D 2.0

# min energy at plane for detection (MeV) ***must be same as in digi!!***
RUN min_eng D 0.000016

# Seed for smear on cluster position.NOT USED!! NEEDS TIDIED!
RUN Gen_seed D 373940592

# sigma (cm)
RUN pos_sig D 0.577
RUN zpos_sig D 0.433

#############################################
#  parameters for the analysis
#############################################

# monitor ntuple file
RUN out_file S %(out_base)s/rec_out/nd_%(part)s%(inttype)s/nd_%(part)s%(inttype)s_%(seed)d.root

#liklihood file.
RUN like_file S %(out_base)s/likelihoods/like_%(part)s_%(inttype)s_%(seed)d.root

# type of fit
RUN kfitter S kalman

# fitter model
RUN model S particle/helix

# verbosities for recpack services
ANA vfit I 0
ANA vmat I 0
ANA vnav I 0
ANA vmod I 0
ANA vsim I 0

# maximum chi2 for tracks
RUN chi2fit_max D 50

# maximum local chi2 for nodes (muon fit)
RUN chi2node_max D 20
RUN max_outliers I 5

# maximum local chi2 for nodes (pattern rec.)
RUN pat_rec_max_chi D 20
RUN pat_rec_max_outliers I 10
RUN max_consec_missed_planes I 3
# minimum proportion of planes used to not reject
# event failing consec planes cut
RUN min_use_prop D 0.7

# Stuff for skipper etc.
RUN maxBlobSkip D 0.2
RUN minBlobOcc D 2

# cuts on cellular automaton trajectories. separation in cm.
RUN max_sep D 7
RUN max_traj I 40
RUN accept_chi D 20
RUN max_coincidence D 0.5

# cut on the maximum number of reconstructed trajectories
RUN max_N_trajectories I 4

# lowest and highest number of hits required for fit.
RUN low_Pass_hits I 8
RUN high_Pass_hits I 500

# proportion of nodes which must be fir not to trigger backwards fit.
RUN low_fit_cut0 D 0.8

# proportion of nodes which must be fit to accept type 2 (free mu) interaction.
RUN low_fit_cut2 D 0.6

# Fiducial cuts. (reduction in z and reduction at both sides in x,y. cm)
RUN z_cut D 300
RUN x_cut D 50
RUN y_cut D 50

# fit data twice (0=false, 1=true)
ANA refit I 1

# multiplication factor for covariance in refit seed
ANA facRef D 10000.

# Relevant parameters for pattern recognition.
# Do Pattern recognition (0=false, 1=true)
ANA patRec I 1

# Min./Max. isolated hits for patRec seed.
ANA min_seed_hits I 7
ANA max_seed_hits I 20
ANA min_check_nodes I 3

# Minumum proportion of isolated hits for central p seed.
ANA min_iso_prop D 0.8

# additional variables for mutliple track analysis
ANA plane_occupancy I 10
ANA max_multOcc_plane I 10

# Bit to tell if require output of likelihood info.
RUN likeli I 1

###  data to read ###
DATA idst_files SV 1
%(out_base)s/digi_out/nd_%(part)s%(inttype)s/nd_%(part)s%(inttype)s_%(seed)d_digi.dst.root

###  data to write ###
#
# DATA odst_files SV 1
# /home/alaing/ntupElec/tnd_e.FITTED.gz
'''% dictionary

		self.print_file(filename,filedata)
