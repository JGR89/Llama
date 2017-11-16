##############################################################
#                                                            #
#                  SET USER PARAMETERS HERE                  #                                                              #
#                                                            #
##############################################################

#Packages
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

#Source Inputs - ages
ages = [10**7, 10**7.5] #Jet ages to make mock observations of
 
#Observing inputs
source_RA = 0 #degrees
source_DEC = -45 #degrees
beam_diameter = 1.0 #Beam diameter in arcseconds
observation_date = '2017-08-04'#YYYY-MM-DD


#Define raise radioEvolution directory, and Llama output directory
import os
raise_radioEvolution_dir = '/Users/Jonathan/Desktop/Projects/PHD_code/RAiSE_files/'
output_directory = '/Users/Jonathan/Desktop/Projects/PHD_code/Llama_v3/Llama_output/'
#'/Users/Jonathan/Desktop/Projects/PHD_code/Llama_v3/Llama_output/'

#Check if the Llama_output directory exists. Raise a ValueError if it does
if os.path.isdir(output_directory) is False:
    os.mkdir(output_directory)
else:
    raise ValueError('output directory already exists')  
    
#Define location of surfbright and ldtrack folders of RAiSE
SB_dir = raise_radioEvolution_dir + 'SurfBright/'
LDtrack_dir = raise_radioEvolution_dir + 'LDtracks/'

##############################################################
#                                                            #
#                       DEFINE CLASSES                       #                                                              
#                                                            #
##############################################################
class Timer(object):    
    '''
        goo.gl/dZ4nBO
        Times imbedded code execution
        Example usage for averaging over 10 runs:

        ts = []   
        for i in range(10):
            with Timer() as t:  
                <RUN CODE>
            ts.append(t.interval)
        print 'mean was', np.mean(ts) 
    '''
    def __enter__(self):
        self.start = time.clock()
        return self

    def __exit__(self, *args):
        self.end = time.clock()
        self.interval = self.end - self.start



class raise_sim:
    """Creates a RAiSE object.
    Inputs are:
    - RAiSE folder name (i.e. where the simulation is)
    - RAiSE file name 
    - redshift of source (as above)
    - log10 of the source luminosity in [W / Hz] (as above)
    - beam size of mock observation (i.e. primary beam diameter in arcsec)
    """
    
    def LDtrack_age_luminosity_size_entry(self,
                                          LDtrack_fname,
                                          input_age,
                                          frequency):
        """
        Gets closest age match for a source in a RAiSE LD track file.

        Inputs:
        - LDtrack file
        - age
        - frequency

        When
        age(n) <= input_age < age(n+1)
        then input_age = age(n) (i.e. the smallest, closest entry)

        Returns: table_index, age, luminosity (at desired frequency)

        """
        from astropy.table import Table
        import bisect
        #LDtrack_fname = self.LDtrack_fname
        #input_age = self.input_age
        #frequency = self.frequency

        LDtrack_table = Table.read(LDtrack_fname, format='ascii')

        ages = LDtrack_table['Time (yrs)'].data
        i = bisect.bisect_left(ages,self.input_age)    
        if i == 0:
            index = i #otherwise it'll give us the age at the end of the table when index = i-1
        else:
            index = i-1

        age = ages[index]

        #Find size of the source (in kpc)
        size_col_index = LDtrack_table.colnames.index('Size (kpc)')
        size = LDtrack_table[index][size_col_index]

        #Find luminosity for this age source, at given freq

        luminosity_freq_col_index = LDtrack_table.colnames.index(('L{} (W/Hz)').format(self.frequency))
        luminosity = LDtrack_table[index][luminosity_freq_col_index]

        return index, age, luminosity, size
    
   
    def __init__(self,
                 SB_filename,
                 LD_filename,
                 redshift,
                 input_age,
                 frequency,
                 beam_diameter):
        
        import numpy as np
        from astropy.cosmology import Planck15 as cosmo
        import astropy.units as u
        from astropy.table import Table
        self.SB_filename = str(SB_filename)
        self.LD_filename = str(LD_filename)
    
        #Source info
        self.z = redshift
        self.frequency = frequency
        self.input_age = input_age
        self.DL = cosmo.luminosity_distance(self.z)
        self.LDtrack_index, self.source_age, self.luminosity, self.source_size = self.LDtrack_age_luminosity_size_entry(self.LD_filename, self.input_age, self.frequency)
        self.beam_diameter = beam_diameter
        
        #Various grids and derived source values
        
        print 'SB_filename is {}'.format(self.SB_filename)
        #self.surfBrightMap = np.loadtxt(fname=self.SB_filename,dtype=object,delimiter="  ") #mJy/beam
        self.surfBrightMap = np.loadtxt(fname=self.SB_filename)
        print 'shape is {}'.format(self.surfBrightMap.shape)
        
        self.grid_length = self.surfBrightMap.shape[0]
        self.grid_height = self.surfBrightMap.shape[1]
        self.surfBrightMap_integrated = sum(self.surfBrightMap.flatten()) #[mJy/beam * npix]
        
       
        #Spatial conversions
        self.arcsec_per_kpc = cosmo.arcsec_per_kpc_proper(redshift) #do I want proper or comoving?
        self.kpc_per_pixel = self.source_size/self.grid_length
        self.arcsec_per_pixel = self.arcsec_per_kpc*self.kpc_per_pixel
        self.cell_area = self.arcsec_per_pixel**2
        
        #Calculate source flux
        self.luminosityMap = (self.luminosity)*(self.surfBrightMap/self.surfBrightMap_integrated)
        self.flux = self.luminosityMap/(4*np.pi*(self.DL.to(u.m))**2)*(u.W/u.Hz)
        
        #Normalise by beam (i.e * by beam area/cell area)
        self.beam_area = np.pi*(self.beam_diameter/2)**2 #Beam area in arcsec^2                                               
        self.flux_beam_normalised = (self.flux/self.beam_area).to(u.Jy) #flux with each cell weighted by ratio of beam size/pixel size
                                        
        return
##############################################################
#                                                            #
#                     DEFINE FUNCTIONS                       #                                                              #
#                                                            #
##############################################################

def FluxDensity_to_fits(raise_object,source_RA,source_DEC,object_name,obs_date,output_dir,output_fname):
    #Convert a map from a raise_class into a fits file with header information.
    """
    Input formats:
    data_array is an NxM numpy array
    source RA and DEC are floats (or ints I suppose)
    object_name is a string, gets inserted into header OBJECT card
    obs_date is the date of the observation. format: YYYY-MM-DD (i think fits needs it to be this way?)
    """
    from datetime import datetime
    from astropy.io import fits
    
    hdu = fits.PrimaryHDU()
    hdu.data = raise_object.flux_beam_normalised.value
    
    time = str(datetime.now())[0:10]
    required_header_inputs = {'CTYPE1': 'RA---SIN', #need to make sure this is correct
                              'CTYPE2': 'DEC--SIN', #need to make sure this is correct
                              'CRVAL1': source_RA,
                              'CRVAL2': source_DEC,
                              'CDELT1': raise_object.arcsec_per_pixel.value/3600.0, #fits standard says it must be in degrees
                              'CDELT2': raise_object.arcsec_per_pixel.value/3600.0,
                              'DATAMAX': raise_object.flux_beam_normalised.value.max(), #broken at the moment as the grid is not calculated correctly
                              'DATAMIN': raise_object.flux_beam_normalised.value.min(), #broken at the moment as the grid is not calculated correctly                            
                              'DATE': time, #date of header generation
                              'DATE-OBS': obs_date, #might need to be user defined?
                              'EQUINOX': 2000.0, #randomly setting this at the moment
                              'OBJECT': str(object_name),
                              'TELESCOP': 'RAiSE_sim',
                              'BUNIT': 'Jy/pixel',
                              'BZERO': 0.0,
                              'BSCALE': 1.0}
    
    #Fail safe for if surface brightness array is 0 (i.e. source unable to be seen at current frequency)
    if np.isnan(raise_object.flux_beam_normalised.value.max()):
        #Might need to do something else here to stop a fits file from being written?
        required_header_inputs['DATAMAX'] = 0
    
    #Floor the min value
    if np.isnan(raise_object.flux_beam_normalised.value.min()):
        #Might need to do something else here to stop a fits file from being written?
        required_header_inputs['DATAMIN'] = 0
    
    #Source centre is always at centre pixel.
    if hdu.header['NAXIS1'] % 2 == 0:
        hdu.header['CRPIX1'] = hdu.header['NAXIS1']/2.0
        hdu.header['CRPIX2'] = hdu.header['NAXIS2']/2.0
    else:
        hdu.header['CRPIX1'] = (hdu.header['NAXIS1']+1)/2.0
        hdu.header['CRPIX2'] = (hdu.header['NAXIS2']+1)/2.0
    
    #Add the header information
    for key, value in required_header_inputs.iteritems():
        hdu.header[key] = value
    
    #Write the fits file (assues output_dir string ends in a /)
    fits.writeto('{}{}'.format(output_dir,output_fname),hdu.data,hdu.header,overwrite=True)
    return

def _gen_uvgen_string(uvgen_dict):
    """
    Generates bash code to print to a miriad shell script.
    Input is a dictionary of input:value pairs for the miriad task uvgen
    ant_numbers is a 
    """
    output_string='uvgen'
    for key,value in uvgen_dict.iteritems():
        if key == 'source' or key == 'ant':
            output_string+=' {}="{}"'.format(key,value)
        elif key == 'select':
            output_string+=' {}=\'{}\''.format(key,value)
        else:
            output_string+= ' {}=\'{}\''.format(key,value).replace('[','').replace(']','')
    return output_string

def _gen_uu_vv_string(uu_vv_dict):
    """
    It plots uu,vv for all baselines on the one plot
    """
    output_string = 'uvplt'
    for key, value in uu_vv_dict.iteritems():
        output_string+= ' {}=\'{}\''.format(key, value)
    return output_string

#Define function to plot uv tracks of source
def gen_shell_script_uvtracks(uvgen_dict,uu_vv_dict):
    fname = 'uvtracks_script.sh'
    shebang = '#!/bin/bash'
    """
    Generates a shell script to display the uv tracks to be observed.
    #Uncomment the elseif statement if you need file paths to have "" instead of ''
    Uses _genuvgen_line(), _gen_fits_line(), _gen_uvmodel_line() and _gen_uvplt_line()
    """

    uvgen_scriptEntry = _gen_uvgen_string(uvgen_dict)
    #uvgen_scriptEntry+='\n'
    uvplt_scriptEntry = _gen_uu_vv_string(uu_vv_dict)
    
    #Write each line to a file
    with open(fname,'a') as f:
        f.write(shebang + '\n')
        f.write(uvgen_scriptEntry + '\n')
        f.write(uvplt_scriptEntry)
    return fname

def run_shell_script(fname):
    """
    shell script must have correct shebang.
    For me, this is #!/bin/bash
    """
    if type(fname) != str:
        raise TypeError('Input must be a string')
    
    import subprocess
    import os
    
    #Make script executable 
    os.system('chmod +x {}'.format(fname))
    
    #Run it
    subprocess.call(['./{}'.format(fname)])
    
    return

def observe_fits(w_dir,
                 fits_image, 
                 Nant,
                 freq,
                 source_RA,
                 source_DEC):
    
    import os
    """
    Creates visibility data of a fits image
    1. Uses MIRIADs UVGEN task to generate uv tracks (i.e. sampling distribution on the sky) of a source
    these are initially done for a point source
    
    2. Uses MIRIADs FITS task to convert a flux density map with each pixel normalised by the primary beam area. 
    
    3. Uses MIRIADs UVMODEL task to replace the model source from a point source to the raise source
    
    4. Uses MIRIADs UVPLT task to get amplitude, phase, dtime (time in days), and uvdistance vectors.
    
    #Works by:
    Generating a list of telescope Baselines 'ant(1)ant(2)' for uvgen, or 'ant(1)(2)' for uvmodel
    makes a folder named with array specs
    Write and runs a MIRIAD shell scripts for each antenna pair
    
    Inputs:
    w_dir: Working directory (will be, for example, /Users/Jonathan/Desktop/Projects/PHD_code/Llama_v1/Llama_output/H=14.00_Q=37.00_z=0.10/freq=10.18)
    Nant: number of antennas in the array, e.g. 5
    freq: frequency of the observation in GHz. e.g. 1.38
    source_RA, source DEC: sin projection? of right ascention and declination of the source
    
    Outputs:
    
    Creates a folder called 'mock_observations', where all of the miriad outputs go
    
    - raise.mir file (miriad version of the fits file)
    - uvgen.vis (visibilities for a point source)
    - raise.vis (visibilities for a RAiSE source)
    - A .dat file for the telescope parameters (a print of uvgen_dict)
    
    
    Within the 'mock_observations' folder,: 
    - MIRIAD a separate folder for each baseline is created.
    A .dat of the visibilities for each baseline (phase, amp dtime (fractional days))
    
    """
    
    def _antlist(Nant):
        antlist=[]
        for i in range(Nant):
            antlist.append('ant(%s)' % str(i+1))
        return antlist  
    
    def gen_VisData(ant_numbers,
                uvgen_dict,
                fits_dict,
                uvmodel_dict,
                uvdist_phase_dict, 
                uvdist_amp_dict,
                dtime_uvdist_dict):
        fname = 'miriad_script.sh'
        shebang = '#!/bin/bash' #this line might cause issues depending on someones setup? This is for bash
        """
        Ant numbers are of the form [1,2]
        """
        
        #Define functions to generate MIRIAD strings
        def _gen_fits_string(fits_dict):
            """
            Inputs the fits image saved from running the FluxDensity_to_fits().
            Converts it to RAiSE format.
            """
            output_string = 'fits'
            for key, value in fits_dict.iteritems():
                output_string+= ' {}=\'{}\''.format(key, value)
            return output_string

        def _gen_uvmodel_string(uvmodel_dict):
            """
            Replaces an observed point source from uvgen with a model
            """
            output_string = 'uvmodel'
            for key, value in uvmodel_dict.iteritems():
                output_string+=' {}=\'{}\''.format(key, value)
            return output_string

        def _gen_uvplt_string(uvplt_dict):
            """
            Outputs a uvplt string for a miriad shell script file.
            Inputs
            """
            output_string = 'uvplt'
            for key, value in uvplt_dict.iteritems():
                output_string+= ' {}=\'{}\''.format(key, value)
            return output_string
        
        #Create shell script strings
        uvgen_scriptEntry = _gen_uvgen_string(uvgen_dict)
        fits_scriptEntry = _gen_fits_string(fits_dict)
        uvmodel_scriptEntry = _gen_uvmodel_string(uvmodel_dict)
        uvplt_uvdist_phase_scriptEntry = _gen_uvplt_string(uvdist_phase_dict)
        uvplt_uvdist_amp_scriptEntry = _gen_uvplt_string(uvdist_amp_dict)
        uvplt_dtime_uvdist_scriptEntry = _gen_uvplt_string(dtime_uvdist_dict)


        #Write the stript entries to fname (hard coded fname is miriad_script.sh)
        with open(fname,'a') as f:
            f.write(shebang + '\n')
            f.write(uvgen_scriptEntry + '\n')
            f.write(fits_scriptEntry + '\n')
            f.write(uvmodel_scriptEntry + '\n')
            f.write(uvplt_uvdist_phase_scriptEntry + '\n')
            f.write(uvplt_uvdist_amp_scriptEntry + '\n')
            f.write(uvplt_dtime_uvdist_scriptEntry)
            f.close()

        return fname
    
    #Define static MIRIAD task dictionaries
    uvgen_dict = {'source':'/Users/Jonathan/Documents/miriad/cat/point.source',
                  'ant':'/Users/Jonathan/Documents/miriad/cat/3.0c.ant',
                  'telescop':'atca',
                  'corr':[0,1,0,104],
                  'freq':freq,
                  'ellim':12,
                  'lat':-30,
                  'jyperk':12.7,
                  #'select':ant_pair can be added as a key if needed
                  'radec':[source_RA,source_DEC],
                  'harange':[-6,6,0.1],
                  'stokes':'i',
                  'cycle':[0.1,0],
                  'gnoise':0,
                  'pnoise':[0,0,0,0],
                  'systemp':0,
                  'baseunit':-51.0204,
                  'out':'point_source.vis'}       
        
    fits_dict = {'in':'{}'.format(fits_image), #location of fits file
                 'op':'xyin',
                 'out':'raise.mir'} 
       
    uvmodel_dict = {'vis':'point_source.vis',
                    'model':'raise.mir',
                    'options':'replace',
                    'out':'raise.vis'}
        
    #Special dictionary for all baselines (may be able to rewrite gen_uvplt_dict to have a special base for nobase=True so i don'thave to define this separate function)
    uvplt_uu_vv_dict = {'vis':'point_source.vis',
                        'line':'wide',
                        'axis':'uu,vv',
                        'options':'nobase, log',
                        'log':'uu_vv.log',
                        'device':'uuvv.png/PNG'}
    
    #Define dynamic uvplt dictionary generator:
    def gen_uvplt_dict(xvar,yvar,ant_numbers):
        """
        Generates a dictionary used to create MIRIAD uvplt taks inputs
        Is a function of baseline
        Possible values for xvar and yvar:
        uu, vv, uvdistance, amplitude, dtime, phase
        ant_numbers is of the form [1,2], for example.
        """
        uvplt_dict = {'vis':'raise.vis',
                      'line':'wide', #do i need this line?
                      'axis':'{},{}'.format(xvar,yvar), #This line is what gets changed to plot different things.
                      'select':'ant({})({})'.format(ant_numbers[0],ant_numbers[1]),
                      'device':'{}{}.png/PNG'.format(xvar,yvar),
                      'options':'log',
                      'log':'{}_{}.log'.format(xvar,yvar)}
        return uvplt_dict
     
        
    ############################################
    #        EXTRACT AND PLOT uu,vv data       #
    ############################################
       
    #Write and run a shell script for the uvtracks (on a point source - model replaced later)
    
    uvtracks_script = gen_shell_script_uvtracks(uvgen_dict,uvplt_uu_vv_dict)
    run_shell_script(uvtracks_script)
                  
    ############################################
    # EXTRACT AND PLOT DATA FROM EACH BASELINE #
    ############################################
        
    #Create dictionaries for each required output (visibilities: phase, amp, uvdistance)
    import itertools
    select = []
    for subset in itertools.combinations(_antlist(Nant),2): #subset is just the antennas
        select.append(subset[0]+subset[1])
    
    #Select is a list of strings of the form ['ant(1)ant(2)','ant(1)ant(3)','....']
        
    for ant_pair in select:

        #Make a directory for current set of telescope parameters
        ant_dir = '{}/{}'.format(w_dir,ant_pair.replace('(','').replace(')',''))
        os.mkdir(ant_dir)
        os.chdir(ant_dir)
        
        #Put the 2 ant numbers into a list [n, n+1], as uvplt requires select=ant(1)(2), NOT ant(1)ant(2)
        ant_numbers = [int(s) for s in list(ant_pair) if s.isdigit()]        
            
        #Define uvplt dicts    
        uvdist_phase_dict = gen_uvplt_dict('uvdistance','phase',ant_numbers) #make sure ant_numbers syntax for miriad for this function is correct
        uvdist_amp_dict = gen_uvplt_dict('uvdistance','amplitude',ant_numbers)
        dtime_uvdist_dict = gen_uvplt_dict('uvdistance','dtime',ant_numbers)
        
        #Write shell script for current antenna pair
        gen_VisData_output = gen_VisData(ant_pair,uvgen_dict,fits_dict,uvmodel_dict, uvdist_phase_dict, uvdist_amp_dict, dtime_uvdist_dict)
        
        #Execute shell script
        run_shell_script(gen_VisData_output)

        os.chdir('../')
        #print 'shell script names is {}'.format(gen_VisData_output)
     
        
        #print 'Observation complete.'
        #print 'Note to self: print telescope specs to file (i.e. uvgen)'
    return

#Step through RAiSE data and make mock observations


                            #################################################################################
                            #           ___       ___       ________  _____ ______   ________               #
                            #          |\  \     |\  \     |\   __  \|\   _ \  _   \|\   __  \              #
                            #          \ \  \    \ \  \    \ \  \|\  \ \  \\\__\ \  \ \  \|\  \             #
                            #           \ \  \    \ \  \    \ \   __  \ \  \\|__| \  \ \   __  \            #
                            #            \ \  \____\ \  \____\ \  \ \  \ \  \    \ \  \ \  \ \  \           #
                            #             \ \_______\ \_______\ \__\ \__\ \__\    \ \__\ \__\ \__\          #
                            #              \|_______|\|_______|\|__|\|__|\|__|     \|__|\|__|\|__|          #
                            #                                                                               #
                            #################################################################################

#Time the script                                                           
import time
start = time.time()

#Cycle through RAiSE output directories (SB and LD tracks)
for sb_subdir, ldtrack_subdir in zip(os.listdir(SB_dir), os.listdir(LDtrack_dir)):
    folder = '{}{}/'.format(output_directory,sb_subdir.replace('surfbright_',''))
    os.mkdir(folder)
    
    #Cycle through the LD track entries closest to each jet age as defined in list ages
    for ldtrack_array_file in os.listdir(LDtrack_dir + ldtrack_subdir):
        
        #Define LD file path
        LD_filename = ldtrack_array_file
        LDtrack_full_path = '{}{}/{}'.format(LDtrack_dir, ldtrack_subdir, LD_filename)
        
        #Make LD dictionary out of filename to get simulation parameters
        LDtrack_parameters = {}
        for entry in LD_filename.split('_')[1:-1]:
            k,v = entry.split("=")
            LDtrack_parameters.setdefault(k,[]).append(v)
        
        #Get redshift from the LD filename dictionary (not strictly necessary, just makes code neater)
        z = float(LDtrack_parameters['z'][0])

        # -- Scale the surface brightness map for each frequency -- #
        
        #Make SB dictionary out of filename to get simulation parameters
        for surfbright_map in os.listdir(SB_dir + sb_subdir):
            surfBright_map_parameters = {}
            for entry in surfbright_map.strip('.dat').strip('surfbright_').split('_'):
                k,v = entry.split("=")
                surfBright_map_parameters.setdefault(k,[]).append(v)
            
            #Get current frequency from the SB dictionary
            freq = surfBright_map_parameters['nu'][0]
            
            #Make a folder for the current frequency
            subfolder = '{}{}/'.format(folder,'freq={}'.format(freq))
            os.mkdir(subfolder)                                                                                

            #Cycle through input ages to get from the RAiSE LD track, for the current frequency
            for age in ages:
                #Define SB file path
                SB_filename = surfbright_map
                SB_full_path = '{}{}/{}'.format(SB_dir, sb_subdir, SB_filename)
                
                #Import the raise simulation
                raise_simulation = raise_sim(SB_full_path, LDtrack_full_path, z, age, freq, beam_diameter)
                
                #convert the flux density map to a fits file and add header information
                object_name = 'age={}_z={}_freq={}_beam_d={}'.format(str(np.log10(age)),z,freq,beam_diameter)
                output_fname = '{}.fits'.format(object_name)
                
                FluxDensity_to_fits(raise_simulation, source_RA, source_DEC, object_name, observation_date,subfolder,output_fname)
                
                #Create folder for the observations of each fits file (i.e, one per source age)
                observation_folder = '{}OBS_AGE{}'.format(subfolder,np.log10(age))
                os.mkdir(observation_folder)
                os.chdir(observation_folder)
                
                #Remember that RAiSE takes log10(freq), whereas MIRIAD (and fhus the observe_fits fn) takes freq in GHz                
                #Check this is correct.
                freq_GHz = (10**float(freq))/10**9
                print 'freq_GHz = ' + str(freq_GHz)
                
                #Generate mock observations using miriad
                with Timer() as t:                
                    observe_fits(os.getcwd(),'../../{}'.format(output_fname),5,freq_GHz,source_RA,source_DEC)
                print 'observe fits takes {}'.format(t.interval)
                # observe fits inputs: fits_image_path, Nant, f req, source_RA, source_DEC
                #Outputs are a folder with the uv tracks.
                #Within this folder are individual folders for each baseline, containing Vis amplitudes, phases, uvdistances and dtime
                
                #os.chdir('../') dont think i need this           
                
print 'COMPLETE'
print 'Llama took ', time.time()-start, 'seconds to run'

