#This short and sweet script pulls out the fluxes from all the different dumps
#adds them together, 4pi r^2's them and then prints out the resulting
#luminosities as a function of radius.

#Version 2 now just prints the radius in a separate file and then big arrays for electron, anti, and mu types
#Where each row corresponds to a radius and each column corresponds to a timestep

import numpy as np
import matplotlib.pyplot as pp
from matplotlib.backends.backend_pdf import PdfPages


start_time_pre = -5
end_time_post   = 15

#r_line = 606

c = 2.998e5 #km/s

r_radii_int = [50,100,200,300,400,500,606]

#print end_time

nbin = 40
skipcols = 8 + 3*nbin

bin_center,bin_width = np.loadtxt("../group_structure.dat",skiprows=5,unpack=True)

u_bin_center = bin_center[2*nbin:3*nbin]
u_bin_width  = bin_width[2*nbin:3*nbin]


models=[ ["wh07.12.40g",167],[ "wh07.15.40g",199], [ "wh07.15.40g.shen",210],[ "wh07.20.40g",261],[ "wh07.25.40g",268]]

for iii in range(len(models)):
    print models[iii][0]
    prefix = "../" + models[iii][0] + "/"
    s = prefix + "dumps/dump_00164"# + str(SBOtime)

    r = np.loadtxt(s, usecols=(0,), skiprows=1,unpack=True)#/1e5 # converts to km
    #print r[r_line]/1.e5
    r_radii = np.divide(r[r_radii_int],1.e5)
    print "radii:"
    print r_radii
    time_difference = np.divide(np.subtract(np.divide(r[r_radii_int[-1]],1.e5),r_radii),c)
    print "time difference:"
    print time_difference


    u_lum_all=[]

    for kk in range(len(r_radii)):
        time_offset = int(round(time_difference[kk] * 1000))

        start_time = models[iii][1] + start_time_pre - time_offset #extra bit to factor in light travel time
        end_time   = models[iii][1] + end_time_post  - time_offset
        print start_time
        print end_time
    
        energyspectra_pages = PdfPages(models[iii][0] + '_energyspectra_rintarrayindex_' + str(kk) + '.pdf')

        for i in range(start_time,end_time+1):
            if i%200 == 0:
                print i
            s = prefix + "dumps/dump_00" + str(i).rjust(3, '0')
            all = np.loadtxt(s)


            u_lum=all[r_radii_int[kk],(skipcols+2*nbin):(skipcols + 3*nbin)]
    #        print e_lum


            u_lum = np.multiply(u_lum, 4 * np.pi * (np.multiply(r_radii[kk],1.e5))**2 / 1e51) #/ u_bin_width[j]


            #u_lum_all.append(np.sum(u_lum))

            u_lum = np.divide(u_lum,u_bin_width) #To put in energy densities

            pp.plot(u_bin_center,u_lum)
            pp.xlabel("Energy (MeV)")
            pp.ylabel("Spectral density (10$^{51}$ erg/s/MeV)")
            pp.title(models[iii][0])
            ymax = 2
            pp.ylim(0,ymax)
            pp.text(55,.9*ymax,"Time is " + str(i-models[iii][1]+time_offset) + " ms from maximum lum")
            pp.text(55,.8*ymax,"Radius is " + str(int(round(r_radii[kk]))) + " km")
            energyspectra_pages.savefig()
            pp.close()
        energyspectra_pages.close()



                       




