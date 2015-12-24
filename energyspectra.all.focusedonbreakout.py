import numpy as np
from matplotlib import pyplot as pp
import matplotlib as mpl
import linecache
from scipy.interpolate import interp1d as intp
#import matplotlib.colors as colors
from matplotlib import cm

lwvalue = 0.8
mpl.rcParams['lines.linewidth'] = 1.8
mpl.rcParams['xtick.labelsize'] = 18
mpl.rcParams['ytick.labelsize'] = 18
mpl.rcParams['axes.linewidth'] = lwvalue
mpl.rcParams['font.family'] = 'serif'
pp.tick_params(axis='both',size=10,width=1.1)

models=[ ["wh07.12.40g",167],[ "wh07.15.40g",199], [ "wh07.15.40g.shen",209],[ "wh07.20.40g",261],[ "wh07.25.40g",267]]
models_label = [r"12 M$_{\odot}$", r"15 M$_{\odot}$", r"15 M$_{\odot}$ SEOS", r"20 M$_{\odot}$", r"25 M$_{\odot}$"]

start_time = -25
end_time   = 11
xlimright = (33,)
nbin = 40

r_line=606

mycolormap = cm.gist_rainbow #You can change the color map by changing this line
m=cm.ScalarMappable(cmap=mycolormap)
m.set_array(range(start_time,end_time))

# "../../Snproject_backup/"
s = "../../" + models[1][0] + "/dumps/dump_00164"# + str(SBOtime)
r = np.loadtxt(s, usecols=(0,), skiprows=1, unpack=True)/1e5 # converts to km

s = "../../" + "group_structure.dat"
bin_center,bin_width = np.loadtxt(s,skiprows=5,unpack=True)
skipcols = 8 + 3*nbin  #this is the 0-indexed number of 
                        #the first relevant column
                                     
for h in (0,):#range(0,3):  #Number of neutrino flavors  
  bin_center_used = bin_center[(h*nbin):((h+1)*nbin)]
  bin_center_finer_spacing = np.linspace(bin_center_used[0],bin_center_used[nbin-1],1000)
  for i in (1,):
      print  models[i][0]
      for j in range(start_time, end_time+1):
         s = "../../" + models[i][0] + "/dumps/dump_00" + str(j+models[i][1]).rjust(3, '0')
         f = open(s)
         line = linecache.getline(s,r_line+2)#+2 to compensate for comment line and 0-indexed array values but 1-indexed line numbers
         linesplit = [float(p) for p in line.split()]
         spectrum = linesplit[skipcols+h*nbin:(skipcols+(h+1)*nbin)]
         luminosity=list()
         for l in range(0,len(spectrum)):
             luminosity.append(spectrum[l] * (4 * np.pi * (r[r_line]*1e5)**2 )/bin_width[l+h*nbin] )
         if len(luminosity)!=nbin:
             print "Wrong number of bins!  It is " + str(len(luminosity))  
  
         lum2 = intp(bin_center_used,luminosity,kind='cubic')
         rgba = mycolormap( (float(j-start_time))/float(end_time-start_time) )
         pp.plot(bin_center_finer_spacing,np.divide(lum2(bin_center_finer_spacing),1.e53),color=rgba)#ls=lineslist[i])
         #pp.scatter(bin_center_used,np.divide(luminosity,1.e53),color=rgba,marker='o')
 
      cbar=pp.colorbar(m)
      #cbar.ax.tick_params(labelsize=16) 
      cbar.set_label('Time Since Peak (ms)',size=16)
      #pp.rc('text',usetex=True)
      pp.xlabel("Energy (MeV)",size=16)
      pp.xlim(0,xlimright[h])
      pp.ylim(0,.32)
      if h==0:
          pp.ylabel(r'L$_{\nu_e}$ (10$^{53}$ erg s$^{-1}$ MeV$^{-1}$)',size=16)
          nutype = "e"
          nutypeword="electron"
      else:
          print "For some reason this code thinks there are more than three neutrino types, so I don't know what label to put."


      stemp = "../figs/" + "nue_energy_spectrum_time.pdf"
      pp.savefig(stemp)
      pp.close()
