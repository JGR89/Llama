#Functions and function testing
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

def ReadMIRIADLog2(location):
    """
    Returns the two columns from a MIRIAD uvplt .log file as numpy arrays.
    location can be full path, or simply the filename, eg. 'uvdistance_amplitude.log'
    If the latter is used, the os.getcwd() must be set to the folder containing the log file.
    """
    import numpy as np
    import sys  
    #reload(sys)
    #sys.setdefaultencoding('utf8')
    
    with open(location) as f:
        fileContents = f.readlines()[8:]
        f.close()
    
    colA = np.array([])
    colB = np.array([])
    
    for item in fileContents:
        colAEntry, colBEntry = item.rsplit(' ')
        colA = np.append([colA],[float(colAEntry)])
        colB = np.append([colB],[float(colBEntry)])
        
    return colA, colB

def thickAxes(ax,linewidth):
    """
    Thickens a matplotlib axes, as set by linewidth
    e.g. ax = plt.plot(x,y)
    thickAxes(ax)
    """
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(linewidth)
    return

def GenSubPlot_uvtracks(uu,vv,showplt=False,figsize=(30,10)):
    fig = plt.figure(1, figsize)
    ax = plt.subplot(111)
    thickAxes(ax,2) #make sure this function actually works
    
    line = plt.plot(uu,vv,'.',markersize=5)
    plt.title('uv tracks',fontsize=16)
    plt.xlabel('u [$k\lambda$]',fontsize=13)
    plt.ylabel('v [$k\lambda$]',fontsize=13)
    
    plt.savefig('uvtracks{}'.format('.png'),bbox_inches="tight")
    
    if showplt is True:
        plt.show()
    else:
        None
    plt.close(fig)
    
    return

def GenSubPlot(uvdist,dtime,phase,amp,showplt=False,figsize=(30,10)):
    """
    When os.getcwd() contains uvdistance_dtime.log, uvdistance_phase.log, and uvdistance_amplitude.log
    this function creates a subplot of the visibilities and saves it to file.
    Setting showplt=True will make the plot appear on the screen.
    The dimensions of these figures is set by figsize
    """
    #Init the figure
    fig = plt.figure(1, figsize)
    
    #Subplot 1
    ax = plt.subplot(131)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(2)
    
    line = plt.plot(uvdist,amp,'b',lw=2)
    plt.title('uvdist v amp', fontsize=13)
    plt.xlabel('uvdist (k$\lambda$)', fontsize=13)
    plt.ylabel('amp (Jy)', fontsize=13)
    
    #Subplot 2
    ax = plt.subplot(132)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(2)
    
    line = plt.plot(uvdist,phase,'b',lw=2)
    plt.title('uvdist v phase', fontsize=13)
    plt.xlabel('uvdist (k$\lambda$)', fontsize=13)
    plt.ylabel('phase (degrees)', fontsize=13)
    
    #Subplot 3
    ax = plt.subplot(133)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(2)
        
    line = plt.plot(dtime,uvdist,'b',lw=2)
    plt.title('dtime v uvdist', fontsize=13)
    plt.xlabel('fractional days', fontsize=13)
    plt.ylabel('uvdist (k$\lambda$)', fontsize=13)
    
    #add more subplots here if need be
    
    plt.suptitle('Visibilities', fontsize=16)
    
    plt.savefig('visibilies{}'.format('.png'),bbox_inches="tight")
    
    if showplt is True:
        plt.show()
    else:
        None
    plt.close(fig)
    return

def compute_baseline_chi_sq(x,mu,sigma,nFreePar=None):
    """
    #Details: http://astronomy.swin.edu.au/~cblake/StatsLecture3.pdf
    
    #Note: Stas told me not to write this myself incase I get it wrong.. So take it out I suppose
    
    -Inputs-
    data: x 
    error of x: sigma
    model: mu
    
    degrees of freedom: dof
        - defaults to False. If set, the probably distribution if the model is correct is also returned.
        - ask ross: how many free parameters does raise have?
    -Outputs-
    chi_sq value: goodness of fit of data (x) to the model (mu)
    
    if number of model free parameters is set, output probability dist
    P(chi**2)
    
    """
    import numpy as np
    
    dif = np.asarray(np.log10(x), dtype=float) - np.asarray(np.log10(mu), dtype=float)
    fraction = dif/np.asarray(sigma, dtype=float)
    chi_sq = np.sum(fraction**2)
    
    if nFreePar != None:
        k = 1 #because below is from P is propto. Do I need to output a dictionary where the keys are k's, and the values are P numpy vectors perhaps?
        nPoints = len(x)
        dof = nPoints - nFreePar
        P = k*(chi_sq**((dof-2.0)/2.0))*np.exp(-chi_sq/2.0)
    else:
        P = None
            
    return chi_sq, P

def GetPsudoVis(psudo_data_dir, baseline):
    """
    Pulls data from a folder that i'm pretending has real data in it.
    When I get real data, this function will be subed out for one that reads the real stuff
    
    Inputs: 
        - baseline: string of form 'ant1ant2', or 'ant1ant3' etc.
        - psudo_data_dir: path of directory containing log files
    
    Outputs:
        - Vectors of data contained in the log files (amplitude, phase, dtime and uvdistance)
    """
    fname1 = 'uvdistance_amplitude.log'
    fname2 = 'uvdistance_dtime.log'
    fname3 = 'uvdistance_phase.log'
    
    psudo_uvdist, psudo_amp = ReadMIRIADLog2('{}/{}/{}'.format(psudo_data_dir,baseline,fname1))
    _, psudo_dtime = ReadMIRIADLog2('{}/{}/{}'.format(psudo_data_dir,baseline,fname2))
    _, psudo_phase = ReadMIRIADLog2('{}/{}/{}'.format(psudo_data_dir,baseline,fname3))
    
    return psudo_amp, psudo_phase, psudo_dtime, psudo_uvdist

def checkLengths(*arg):
    """
    Raises an error if the lenth of all the input lists are not equal
    Possible issue: if the inputs are lists or numpy vectors might be an issue? I'm using lists
    
    Inputs: Lists (not tested for numpy vectors)
    """
    import numpy as np
    n = len(arg[0])
    
    if all(len(x) == n for x in arg):
        pass
    else:
        raise ValueError('Lists differ in length')
    
    return
#Next piece of Llama code needs to cycle through all the baselines, and compare them to a set of real baselines
######
######
######

import time
import os
import matplotlib.pyplot as plt
import numpy as np

#Create an empty dictionary
chi_sq_dict = dict({})

plot_suppress = True

#data_dir = Not existing at the moment as I do not have real data
psudo_data_dir = '/Users/Jonathan/Desktop/Projects/PHD_code/Llama_v3/psudo_vis/'
w_dir = '/Users/Jonathan/Desktop/Projects/PHD_code/Llama_v3/Llama_output/'

#Define main analysis function
def Llama_analysis(w_dir,print_visplots=False,pltshow=False):
    
    """
    This function cycles through the subfolders of w_dir (i.e. Llama_output) and compares the 
    generated visibilities from Llama with a set of real visibilities (data_dir, or psudo_data_dir for now)
    
    Folder structure of Llama output:
    Llama_output/
         |--H=#_Q=#_z=#/ (RAiSE output folders for halo mass [H], jet power [Q] and redshift [z])
                |--freq=#/ (RAiSE output folders for each simulated frequency)
                     |--age=#_z=#_freq=#_beam_d=#.fits (Scaled RAiSE SB map to flux density [Jy/pix? can't remember. check.])
                     |--OBS_AGE#/
                           |--point_source.vis/ (Used by MIRIAD)
                           |-uuvv.png (uv tracks image)
                           |-To Do: print A DATA FILE OF UUVV TOO SO I CAN DO MY OWN PLOTS 
                           |-uvtracks_script.sh (change name to generate_uvtracks.sh or something to be more descriptive)
                           |--ant#ant#+1/ (baseline folders of ant # and ant #+1 for the simulated array)
                                 |-miriad_script.sh (change name to sim_observation_miriad.sh or something)
                                 |-uvdistance_dtime.log
                                 |-uvdistance_phase.log
                                 |-uvdistance_amplitude.log
                                 |--assorted MIRIAD output files
    
    
    Outputs: A dictionary that gives the chi-sq values for each baseline 
    """
    
    #Define functions here
    def getBaselineAntennas(baseline):
        """
        Input: baseline = 'antNantM' where M = N+1 [e.g ant1ant2]
        Output: N and M given the string above
        """
        split_string = baseline.split('ant')
        ants = []
        for entry in split_string:
            if len(entry) is 0:
                None
            else:
                ants.append(entry)
        return ants[0], ants[1]
    
    #Cycling through the Llama_output directories:
    
    for redshift_folder in os.listdir(w_dir): #Cycle through H=#_Q=#_z=# folders
        
        #If foldername begins with H, proceed (this is a shitty way of doing this)
        if redshift_folder[0] != 'H':
            None
        else:
            print '--------- Processing {} ---------'.format(redshift_folder)
            os.chdir('{}{}'.format(w_dir,redshift_folder))

            for freq_folder in os.listdir(w_dir + redshift_folder): #Cycle through freq=# folders
                print 'Processing {}'.format(freq_folder)
                os.chdir('{}{}/{}'.format(w_dir,redshift_folder,freq_folder))

                freq_folder_contents = os.listdir(os.getcwd())
                for OBS_AGE_entry in freq_folder_contents: #Cycle through OBS_AGE=# for current freq.
                    
                    #If foldername begins with OBS_AGE, proceed.
                    if OBS_AGE_entry[0:7] == 'OBS_AGE':
                        print 'Processing {}'.format(OBS_AGE_entry)
                        os.chdir('{}{}/{}/{}'.format(w_dir,redshift_folder,freq_folder,OBS_AGE_entry))
                        
                        OBS_AGE_folder_contents = os.listdir(os.getcwd())
                        
                        #Initialise a list for the chi-sq values for the current observation 
                        bline_chisqs = list([]) #issues appending to an empty numpy array in a loop, so use a list
                        
                        for baseline in OBS_AGE_folder_contents: #cycle through each baseline, but first..                            
                             
                            #Process entire array uv tracks before actually analysing individual baselines
                            if baseline == 'uu_vv.log': 
                                uu,vv = ReadMIRIADLog2('uu_vv.log')
                                
                                if print_visplots is True:
                                    GenSubPlot_uvtracks(uu,vv,showplt=False)
                                    
                            #elif statement: add in nobase plots here for amp and phase if need be
                                                   
                            #Process each antenna file/baseline 
                            elif baseline[0:3] == 'ant': #If the folder entry is a baseline folder (e.g. ant1ant2), continue (shitty way of doing this)
                                
                                #Enter into directory #Do i need to change all of the directory cding?
                                os.chdir('{}{}/{}/{}/{}'.format(w_dir,redshift_folder,freq_folder,OBS_AGE_entry,baseline))
                                
                                ant1, ant2 = getBaselineAntennas(baseline)
                                print 'Processing baseline {}-{}'.format(ant1,ant2)
                                
                                #Pull data from logs
                                with Timer() as t:
                                    uvdist,dtime = ReadMIRIADLog2('uvdistance_dtime.log')
                                    _,phase = ReadMIRIADLog2('uvdistance_phase.log')
                                    _,amp = ReadMIRIADLog2('uvdistance_amplitude.log')
                                print 'MIRIAD log files read in {} seconds'.format(t.interval)

                                #Check lengths. If they are unequal, raise an error because there's an issue
                                checkLengths(uvdist,dtime,phase,amp)
                                
                                #======================================#
                                # -- Check against real visibilities --
                                #======================================#
                                #Read in the real data (psudo data for now)
                                with Timer() as t:
                                    psudo_amp, psudo_phase, psudo_dtime, psudo_uvdistance = GetPsudoVis(psudo_data_dir,baseline) #does this really need to happen with every loop? Think about it
                                    sigma_amp = np.sqrt(np.std(np.log10(amp))**2 + (0.2)**2)
                                print 'GetPsudoVis took {}'.format(t.interval)
                                #Calculate chi-sq value for current baseline
                                print 'computing chi sq for baseline'
                                
                                with Timer() as t:
                                    chi_sq,_ = compute_baseline_chi_sq(psudo_amp,amp,sigma_amp,nFreePar=None)
                                print 'chi sq calculated in {}'.format(t.interval)
                                #Add sim parameters to dict. key:value = simulation parameters:chi_sq value
                                bline_chisqs.append(chi_sq) #save the current baselne chisq val
                                
                                #print 'current recorded values for {}: {}'.format(baseline,bline_chisqs)                    
                                #print 'current dict sum: {}'.format(sum(bline_chisqs))
                                
                                #Make nicer plots and save them to the folder (if toggled on)
                                if print_visplots is True:
                                    #'Printing plots to {}'.format(os.getcwd())
                                    GenSubPlot(uvdist,dtime,phase,amp,showplt=False)
                                else:
                                    None
                                
                                #Return to OBS_AGE folder ready to process next baseline
                                os.chdir('../') 
                            else: #skip the folder, i.e. if it doesn't begin with 'ant'
                                #print 'Skipping {}'.format(baseline)
                                None
                        
                        #at the moment this is only for the amplitude.
                        print 'Writing to dict key: {}-{}-{}'.format(redshift_folder,freq_folder,OBS_AGE_entry)
                        
                        #Add the chi sq values together
                        chi_sq_dict['{}-{}-{}'.format(redshift_folder,freq_folder, OBS_AGE_entry)] = sum(bline_chisqs)
                        #print 'Observation chi-sq sum: {}'.format(np.sum(bline_chisqs))
                        #print '~~~~~ CURRENT DICT: {} ~~~~~'.format(chi_sq_dict)
                        
                        #Calculate probability, i.e. e^(-sum(chisq)/2), where ive already done the sum
                        
                        #Change back one directory ready to enter into the next baseline folder
                        os.chdir('../')
                        #baseline loop ends here

                    else: #Skip folder if it doesn't beigin with 'OBS_AGE'
                        #print 'Skipping {}'.format(OBS_AGE_entry)
                        None
            
            #Cd back one directory to process next frequency (i.e. go back to redshift directory)
            os.chdir('../')
            
        #Do I need an os.chdir here too? To head back to the Llama_output dir ready for the next redshift loop
        #i.e. ready for the next H=#_Q=#_z=# foldeR?
    return chi_sq_dict

start=time.time()
chi_sq = Llama_analysis(w_dir,print_visplots=False)
#print 'chi_sq keys: {}'.format(list(chi_sq.keys()))
print 'Llama_analysis took ', time.time()-start, 'seconds to run'

#We want chi-sq probabilities. Need P(chi_sq) propto e^(-chi_sq/2)
#Could probably add this into the actual function
def chi_sq_prob(chi_sq_dict):
    """
        Converts a dictionary of chi_sq values to a probability distribution
        I.e. takes values, chi_sq, and computes e^(-chisq/2) for each value

        Inputs: dictionary of chi-sq values
        Outputs: dictionary of chi_sq probabilities
    """
    import copy

    #create a new isolated dictionary
    prob_dict = copy.deepcopy(chi_sq_dict)

    #Convert the values using e**(-value/2.0)
    prob_dict.update((x, np.exp(-y/2.0)) for x, y in prob_dict.items())

    return prob_dict

prob_dict = chi_sq_prob(chi_sq)

#Order the dictionary by values
#see https://www.saltycrane.com/blog/2007/09/how-to-sort-python-dictionary-by-keys/

#Create a table of halo masses, jet powers, redshifts, freq, log10(age), chi_sq, P(chi_sq)


def interpret_analysis_key(dict_key_string):
    """
        Given a dictionary key of the form of the example (this function can be updated if more paramters are added to the name)
        Example key_list item: H=14.00_Q=37.00_z=0.75-freq=9.68-OBS_AGE7.0
        Example key_list value: 14, 37, 0.75, 9.68, 7.0
    """

    dict_key_string = dict_key_string.replace('OBS_AGE','OBSAGE=').replace('-','_')
    
    dict_key_string = dict_key_string.split('_')
    parameters = {}
    for entry in dict_key_string:
        value = entry.split('=')
        print 'value: {}, value[0]: {}'.format(value, value[0])
        parameters[value[0]]=value[1]
        
    print 'parameters: {}'.format(parameters)
    if len(parameters) < 5:
        raise ValueError('error in interpret_analysis_key usage.\nFunction return only outputs 5 variables.\ni.e. There were more than 4 RAiSE simulation parameters defined in the simuulation\nIf this is raised in error, comment this out. Otherwise, update the function.')

    H = parameters['H']
    Q = parameters['Q']
    z = parameters['z']
    freq = parameters['freq']
    age = parameters['OBSAGE']
    
    return H, Q, z, freq, age


#Put all the values into table format (easier to do Beysian statistics over)

#Create an empty table?
#CRAETE IT HERE
from astropy.table import Table
column_names = ['H', 'Q', 'z', 'freq', 'age', 'chi_sq', 'P(chi_sq']
simulation_table = Table(names=column_names)

#Need to add another loop here to cycle throgh chi_sq, and P(chi_sq)

for key,P_chi_sq in prob_dict.iteritems(): #distribute dictionary informatuon into a table
    H, Q, z, freq, age = interpret_analysis_key(key)
    
    #Get chi_sq value from dict. Key is that same as prob_dict keys
    print 'key: {}'.format(key)
    print 'chi_sq: {}'.format(chi_sq)
    chi_sq_val = chi_sq[key]
    

    row_values = [H, Q, z, freq, age, chi_sq_val, P_chi_sq]
    
    #Raise an error if the table in being incorrectly filled
    if len(column_names) != len(row_values):
        raise ValueError('more variables were read than were input into the table')
    
    #Append row to table, as well as the chi_sq value and P(chi_sq) value
    simulation_table.add_row(row_values)
    print('row_values added to table: {}'.format(row_values))
    print('column_names:             {}').format(column_names)

    #Assignmnt column_values into appropiate table columns
    #ADD VALUES TO TABLE DEFINED JUST OUTSIDE OF FOR LOOP

#Do statistics with the table now

simulation_table.pprint(max_lines=-1, max_width=-1)