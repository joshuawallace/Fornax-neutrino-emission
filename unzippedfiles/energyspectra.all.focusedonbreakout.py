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

start_time = -40
end_time   = 80
xlimright = (33,)
nbin = 40

#r_line=606

mycolormap = cm.gist_rainbow #You can change the color map by changing this line
m=cm.ScalarMappable(cmap=mycolormap)
m.set_array(range(start_time,end_time))

s = "nue_energy_bins.txt"
e_bin_center,e_bin_width = np.loadtxt(s,unpack=True)
s = "nua_energy_bins.txt"
a_bin_center,a_bin_width = np.loadtxt(s,unpack=True)
s = "nux_energy_bins.txt"
x_bin_center,x_bin_width = np.loadtxt(s,unpack=True)


for i in range(len(models)):
  print models[i][0]

  e_spectra = []
  a_spectra = []
  x_spectra = []

  for j in range(start_time,end_time+1):
    s = models[i][0] + "/" + "nu_spectra.t_"
    if j<0:
        s = s + "n" + str(abs(j)).rjust(3, '0') + ".dat"
    else:
        s = s + str(j).rjust(3, '0') + ".dat"

    temp = np.loadtxt(s,unpack=True)
    e_spectra.append(temp[0])
    a_spectra.append(temp[1])
    x_spectra.append(temp[2])

  #e spectra
  e_x = np.linspace(e_bin_center[0],e_bin_center[-1],550)
  for j in range(len(e_spectra)):
    #print len(e_bin_center)
    #print len(e_spectra[j])
    #print e_spectra[j]
    e_cubic = intp(e_bin_center,e_spectra[j],kind='cubic')
    rgba = mycolormap( float(j)/float(end_time-start_time)  )
    pp.plot(e_x,e_cubic(e_x),color=rgba)
                                    
 
  cbar=pp.colorbar(m)
      #cbar.ax.tick_params(labelsize=16) 
  cbar.set_label('Time Since Peak (ms)',size=16)
      #pp.rc('text',usetex=True)
  pp.xlabel("Energy (MeV)",size=16)
  pp.xlim(0,30.)
  pp.ylim(bottom=0)
#      pp.ylim(0,.32)
  pp.title("Electron neutrino")
  pp.ylabel(r'Spectral enegy density (10$^{51}$ erg s$^{-1}$ MeV$^{-1}$)',size=16)

  stemp = models[i][0] + "_nue_energy_spectrum_time.pdf"
  pp.savefig(stemp)
  pp.close()

    #a spectra
  a_x = np.linspace(a_bin_center[0],a_bin_center[-1],550)
  for j in range(len(a_spectra)):
    a_cubic = intp(a_bin_center,a_spectra[j],kind='cubic')
    rgba = mycolormap( float(j+start_time)/float(end_time-start_time)  )
    pp.plot(a_x,a_cubic(a_x),color=rgba)
                                    
 
  cbar=pp.colorbar(m)
      #cbar.ax.tick_params(labelsize=16) 
  cbar.set_label('Time Since Peak (ms)',size=16)
      #pp.rc('text',usetex=True)
  pp.xlabel("Energy (MeV)",size=16)
  pp.xlim(0,40.)
  pp.ylim(bottom=0)
#      pp.ylim(0,.32)
  pp.title("Anti-electron neutrino")
  pp.ylabel(r'Spectral enegy density (10$^{51}$ erg s$^{-1}$ MeV$^{-1}$)',size=16)

  stemp = models[i][0] + "_nua_energy_spectrum_time.pdf"
  pp.savefig(stemp)
  pp.close()

      #x spectra
  x_x = np.linspace(x_bin_center[0],x_bin_center[-1],550 )
  for j in range(len(x_spectra)):
    x_cubic = intp(x_bin_center,x_spectra[j],kind='cubic')
    rgba = mycolormap( float(j+start_time)/float(end_time-start_time)  )
    pp.plot(x_x,x_cubic(x_x),color=rgba)
                                    
 
  cbar=pp.colorbar(m)
      #cbar.ax.tick_params(labelsize=16) 
  cbar.set_label('Time Since Peak (ms)',size=16)
      #pp.rc('text',usetex=True)
  pp.xlabel("Energy (MeV)",size=16)
  pp.xlim(0,65.)
  pp.ylim(bottom=0)
#      pp.ylim(0,.32)
  pp.title("Muon/tau neutrino")
  pp.ylabel(r'Spectral enegy density (10$^{51}$ erg s$^{-1}$ MeV$^{-1}$)',size=16)

  stemp = models[i][0] + "_nux_energy_spectrum_time.pdf"
  pp.savefig(stemp)
  pp.close()
