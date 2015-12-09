#######################################################################################################################
#Created by Patrik Hallsjo @ University of Glasgow
#Need automatic dating through GIT, 
#Modified on 4/11-2015
#Created on 25/9-2015
#######################################################################################################################
#General python import
#######################################################################################################################
import os
#######################################################################################################################
#Importing own python files
#######################################################################################################################

#######################################################################################################################
#Class generation
#######################################################################################################################
class print_config:

	def __init__(self,GenerationMode):
		self.GenerationMode = GenerationMode

################################################
	def print_file(self,filename,data):
	    '''
	    Print data to file filename.
	    '''
	    outfile = open(filename,'w+')
	    outfile.write(data)
	    outfile.close()

	def print_digi_config(self,dictionary,filename):

		filedata = self.print_general('CON',dictionary)
		self.print_file(filename,filedata)

	def print_mindG4_config(self,dictionary,filename):

		filedata = self.print_general('GEOMETRY',dictionary)

		if(self.GenerationMode == 'SINGLE_PARTICLE'):
			filedata += '''

### Three generators, SINGLE_PARTICLE, NUANCE and GENIE, with different
### configuration parameters, are available. The first is chosen
### as default. Comment out the following lines if you want to
### use the latter.
GENERATION generator S SINGLE_PARTICLE


GENERATION particle_name S %(part)s

### Particle energy will be sample between these two values (in GeV)
GENERATION energy_min    D 0.05
GENERATION energy_max    D 4.000

'''% dict(dictionary, **vars(self))
		elif(self.GenerationMode == 'GENIE'):
			filedata += '''
GENERATION generator S GENIE
'''% dict(dictionary, **vars(self))			

		

		filedata += '''
### JOB options ###########################################
JOB output_dst    S %(out_base)s/G4_out/nd_%(part)s%(inttype)s/nd_%(part)s%(inttype)s_%(seed)s.dst.root
JOB number_events I %(Nevts)s

JOB random_seed I 13243%(seed)s

### GEOMETRY configuration parameters #####################

### Uncomment if you want to simulate TASD.
## Sets all volumes as sensitive detectors.
## Remember to use G4_POLYSTYRENE as passive_material as well!!.
# GEOMETRY TASD I 0

### Air gaps between layers (in mm) must be half of the airgap!
GEOMETRY gap1 D %(MIND_thickness_air_mm)s   
GEOMETRY gap2 D 0
#GEOMETRY gap3 D 2.5
#GEOMETRY gap4 D 2.5

### GDML
GEOMETRY useGDML I  %(useGDML)s
GEOMETRY writeGDML I 0
GEOMETRY GDMLFileName S %(xml_file_path)s


### MAgnetic field. Still uniform vector.(T)
## A single bit to turn the uniform field off in favour of a toroidal field
GEOMETRY isUniform I 0

# if uniform then [bx,by,bz]. if toroidal then [br,btheta,bphi]
# magnetic field, not used but must be given
GEOMETRY field DV 3
0.
1.
0.

# If we wish to use a field map then the identity of the field map
# must be entered here.
GEOMETRY FieldMap S %(field_map_full_name)s

GEOMETRY FieldScaling D %(Bfield)s

### GENERATION configuration parameters ################### 

### NUANCE/GENIE data files for passive and active materials!! Set the filenames here!!

GENERATION active_material_data S %(out_base)s/genie_samples/nd_%(part)s%(inttype)s/ev0_%(ASeed)s_%(pid)d_1000060120[0.922582],1000010010[0.077418]_%(Nevts)s.root
GENERATION passive_material_data S %(out_base)s/genie_samples/nd_%(part)s%(inttype)s/ev0_%(seed)s_%(pid)d_1000260560_%(Nevts)s.root
#
### Vertex location (RANDOM, ACTIVE, PASSIVE, FIXED, GAUSS).
GENERATION vertex_location S GAUSS


### Special simulation for training purpose
### Use the muon but not the hadronization from GENIE
# GENERATION FSL_Select I 0

# Vertex if FIXED requested.
GENERATION fvert DV 3
0.
0.
-2000
 
GENERATION bspot DV 2
100.
100.

GENERATION costh_min D 1 
GENERATION costh_max D 0.9

### PHYSICS configuration parameters ######################

### Production cut (in mm) for secondaries
PHYSICS production_cut D 30.
#
### Minimum Kinetic energy for a particle to be tracked (MeV).
PHYSICS minimum_kinEng D 10. #100.
'''% dict(dictionary, **vars(self))

		self.print_file(filename,filedata)

	def print_rec_config(self,dictionary,filename):
		filedata = self.print_general('RUN',dictionary)

		filedata += '''

#############################################
#  parameters for the setup
#############################################

#Parsed xml_file
RUN xml_parsed S %(parsed_file_path)s

# relative density, Sc/Fe, AIR/Sc.
RUN rel_denSI D %(config_rec_rel_denSI)s
RUN rel_denAS D %(config_rec_rel_denAS)s

# Attenuation in Wavelength shifting fibre
RUN WLSatten D %(config_rec_WLSatten)s

# measurement type.
RUN meas_type S %(config_rec_meas_type)s

# magnetic field, not used but must be given
RUN mag_field DV 3
0.
1.
0.

RUN mag_field_map S %(field_map_full_name)s

RUN fieldScale D %(Bfield)s
########

# energy loss (MeV/mm)
RUN de_dx_scint D %(MIND_active_de_dx)s
RUN de_dx_fe D %(MIND_passive_de_dx)s
RUN de_dx_min D %(MIND_module_de_dx)s

# Position resolution for detector (cm).
RUN pos_res D %(config_rec_pos_res)s

# Step size for track fitting (cm).
RUN StepSize D %(config_rec_step_size)s

# name of detector for hit getter.
RUN detect S %(config_rec_detect)s

########
# For hit clustering.(edge in cm)
RUN do_clust I %(config_rec_do_clust)s
RUN rec_boxX D 0.966 #%(config_rec_boxX)s

# min energy at plane for detection (MeV) ***must be same as in digi!!***
RUN min_eng D %(MIND_min_eng_at_plane)s

# Seed for smear on cluster position
#RUN Gen_seed D %(config_rec_gen_seed)s

# sigma (cm)
RUN pos_sig D %(MIND_width_sigma)s 
RUN zpos_sig D %(MIND_thickness_sigma)s 

#############################################
#  parameters for the analysis
#############################################

# monitor ntuple file
RUN out_file S %(out_base)s/rec_out/nd_%(part)s%(inttype)s/nd_%(part)s%(inttype)s_%(seed)s.root

#liklihood file.
RUN like_file S %(out_base)s/likelihoods/like_%(part)s_%(inttype)s_%(seed)s.root

# type of fit
RUN kfitter S %(config_rec_kfitter)s

# fitter model
RUN model S %(config_rec_model)s

# verbosities for recpack services
ANA vfit I %(config_rec_vfit)s
ANA vmat I %(config_rec_vmat)s
ANA vnav I %(config_rec_vnav)s
ANA vmod I %(config_rec_vmod)s
ANA vsim I %(config_rec_vsim)s

# maximum chi2 for tracks
RUN chi2fit_max D %(config_rec_chi2fit_max)s

# maximum local chi2 for nodes (muon fit)
RUN chi2node_max D %(config_rec_chi2node_max)s
RUN max_outliers I %(config_rec_max_outliners)s

# maximum local chi2 for nodes (pattern rec.)
RUN pat_rec_max_chi D %(config_rec_pat_rec_max_chi)s
RUN pat_rec_max_outliers I %(config_rec_pat_rec_max_outliners)s
RUN max_consec_missed_planes I %(config_rec_max_consec_missed_planes)s
# minimum proportion of planes used to not reject
# event failing consec planes cut
RUN min_use_prop D %(config_rec_min_use_prop)s

# Stuff for skipper etc.
RUN maxBlobSkip D %(config_rec_max_blobSkip)s
RUN minBlobOcc D %(config_rec_min_blobOcc)s

# cuts on cellular automaton trajectories. separation in cm.
RUN max_sep D %(config_rec_max_sep)s
RUN max_traj I %(config_rec_max_traj)s
RUN accept_chi D %(config_rec_accept_chi)s
RUN max_coincidence D %(config_rec_max_coincidence)s

# cut on the maximum number of reconstructed trajectories
RUN max_N_trajectories I %(config_rec_max_n_trajs)s

# lowest and highest number of hits required for fit.
RUN low_Pass_hits I %(config_rec_low_pass_hits)s
RUN high_Pass_hits I %(config_rec_high_pass_hits)s

# proportion of nodes which must be fir not to trigger backwards fit.
RUN low_fit_cut0 D %(config_rec_low_fit_cut0)s

# proportion of nodes which must be fit to accept type 2 (free mu) interaction.
RUN low_fit_cut2 D %(config_rec_low_fit_cut2)s

# Fiducial cuts. (reduction in z and reduction at both sides in x,y. cm)
RUN z_cut D %(config_rec_z_cut)s
RUN x_cut D %(config_rec_x_cut)s
RUN y_cut D %(config_rec_y_cut)s

# fit data twice (0=false, 1=true)
ANA refit I %(config_rec_refit)s

# multiplication factor for covariance in refit seed
ANA facRef D %(config_rec_fac_refit)s

# Relevant parameters for pattern recognition.
# Do Pattern recognition (0=false, 1=true)
ANA patRec I %(config_rec_pat_rec)s

# Min./Max. isolated hits for patRec seed.
ANA min_seed_hits I %(config_rec_min_seed_hits)s
ANA max_seed_hits I %(config_rec_max_seed_hits)s
ANA min_check_nodes I %(config_rec_min_check_nodes)s

# Minumum proportion of isolated hits for central p seed.
ANA min_iso_prop D %(config_rec_min_iso_prop)s

# additional variables for mutliple track analysis
ANA plane_occupancy I %(config_rec_plane_occupancy)s
ANA max_multOcc_plane I %(config_rec_max_multOcc_plane)s

# Bit to tell if require output of likelihood info.
RUN likeli I %(config_rec_likelihood)s

###  data to read ###
DATA idst_files SV 1
%(out_base)s/digi_out/nd_%(part)s%(inttype)s/nd_%(part)s%(inttype)s_%(seed)s_digi.dst.root

###  data to write ###
#
# DATA odst_files SV 1
# /home/alaing/ntupElec/tnd_e.FITTED.gz
'''% dictionary

		self.print_file(filename,filedata)

	def print_general(self, inPreParam, dictionary):
		self.preParam = inPreParam


		filedata = '''
'''

		if(self.preParam == 'CON'):
			filedata+= '''
# number of events to be processed.
RUN nEvents I %(Nevts)d

# gausian sigma for smear (cm)
RUN Gaus_Sigma D %(config_digi_gaus_sigma)s

# energy smear (%%)
RUN Eng_Res D %(config_digi_eng_res)s

# seed value for random generator
CON Gen_seed D %(config_digi_seed)s

CON rec_boxX D 0.966 #%(MIND_width_active)s
CON rec_boxY D 0.666  #%(MIND_width_active)s

CON nVoxX I 47;

# Attenuation in Wavelength shifting fibre
CON WLSatten D %(config_rec_WLSatten)s

# minimum energy at plane to be detected.(MeV)
CON min_eng D %(MIND_min_eng_at_plane)s

DATA idst_files SV 1
%(out_base)s/G4_out/nd_%(part)s%(inttype)s/nd_%(part)s%(inttype)s_%(seed)s.dst.root

DATA odst_file S %(out_base)s/digi_out/nd_%(part)s%(inttype)s/nd_%(part)s%(inttype)s_%(seed)s_digi.dst.root
'''% dict(dictionary, **vars(self))

		filedata += '''

%(preParam)s IsOctagonal I %(MIND_type)s

# MIND dimensions (m)
%(preParam)s MIND_x D %(MIND_xdim)s
%(preParam)s MIND_y D %(MIND_ydim)s
%(preParam)s MIND_z D %(MIND_zdim)s
%(preParam)s vertex_x D %(MIND_vertex_xdim)s
%(preParam)s vertex_y D %(MIND_vertex_ydim)s
%(preParam)s vertex_z D %(MIND_vertex_zdim)s //vertexDepth
%(preParam)s ear_height D %(MIND_ear_ydim)s
%(preParam)s ear_width  D %(MIND_ear_xdim)s
%(preParam)s bore_diameter D %(MIND_bore_diameter)s

#Mind internal dimensions
%(preParam)s active_material S %(MIND_active_mat)s
%(preParam)s active_thickness D %(MIND_thickness_active)s //widthSc
%(preParam)s x0Sc D %(MIND_rad_length_active)s
%(preParam)s active_layers I %(MIND_active_layers)s //nlayers
%(preParam)s passive_material S %(MIND_passive_mat)s
%(preParam)s passive_thickness D %(MIND_thickness_passive)s //widthI
%(preParam)s x0Fe D %(MIND_rad_length_passive)s

%(preParam)s bracing_material S %(MIND_bracing_mat)s 
%(preParam)s bracing_thickness D %(MIND_thickness_bracing)s // widthAl

%(preParam)s air_gap D %(MIND_thickness_air)s //withA
%(preParam)s x0AIR D %(MIND_rad_length_air)s

'''% dict(dictionary, **vars(self))

		return filedata

#######################################################################################################################
#File specific functions
#######################################################################################################################




