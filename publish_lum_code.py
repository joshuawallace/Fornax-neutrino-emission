#This short and sweet script pulls out the fluxes from all the different dumps
#adds them together, 4pi r^2's them and then prints out the resulting
#luminosities as a function of radius.

#Version 2 now just prints the radius in a separate file and then big arrays for electron, anti, and mu types
#Where each row corresponds to a radius and each column corresponds to a timestep

import numpy as np
import matplotlib.pyplot as pp

start_time_pre = -40
end_time_post   = 80

r_line = 606

#print end_time

nbin = 40
skipcols = 8 + 3*nbin

bin_center,bin_width = np.loadtxt("../group_structure.dat",skiprows=5,unpack=True)
e_bin_center = bin_center[0:nbin]
e_bin_width  = bin_width[0:nbin]
a_bin_center = bin_center[nbin:2*nbin]
a_bin_width  = bin_width[nbin:2*nbin]
u_bin_center = bin_center[2*nbin:3*nbin]
u_bin_width  = bin_width[2*nbin:3*nbin]


models=[ ["wh07.12.40g",167],[ "wh07.15.40g",199], [ "wh07.15.40g.shen",210],[ "wh07.20.40g",261],[ "wh07.25.40g",268]]

for iii in range(len(models)):
    print models[iii][0]
    prefix = "../" + models[iii][0] + "/"
    s = prefix + "dumps/dump_00164"# + str(SBOtime)

    r = np.loadtxt(s, usecols=(0,), skiprows=1,unpack=True)#/1e5 # converts to km
    print r[r_line]/1.e5


    e_lum_all=[]
    a_lum_all=[]
    u_lum_all=[]

    start_time = models[iii][1] + start_time_pre
    end_time   = models[iii][1] + end_time_post

    for i in range(start_time,end_time+1):
        if i%200 == 0:
            print i
        s = prefix + "dumps/dump_00" + str(i).rjust(3, '0')
        all = np.loadtxt(s)


        e_lum=all[r_line,skipcols:(skipcols + nbin)]
        a_lum=all[r_line,(skipcols+nbin):(skipcols + 2*nbin)]
        u_lum=all[r_line,(skipcols+2*nbin):(skipcols + 3*nbin)]
#        print e_lum
        
        #e_lum = np.sum(e_lum,axis=1)
        #a_lum = np.sum(a_lum,axis=1)
        #u_lum = np.sum(u_lum,axis=1)


        for j in range(0,len(e_lum)):
            e_lum[j]=e_lum[j] * 4 * np.pi * (r[r_line])**2 / 1e51 #/ e_bin_width[j]
            a_lum[j]=a_lum[j] * 4 * np.pi * (r[r_line])**2 / 1e51 #/ a_bin_width[j]
            u_lum[j]=u_lum[j] * 4 * np.pi * (r[r_line])**2 / 1e51 #/ u_bin_width[j]

        #make appropriate values in nux spectra 0
        for j in range(36,40):
            u_lum[j] = 0.

        e_lum_all.append(np.sum(e_lum))
        a_lum_all.append(np.sum(a_lum))
        u_lum_all.append(np.sum(u_lum))

        e_lum = np.divide(e_lum,e_bin_width) #To put in energy densities
        a_lum = np.divide(a_lum,a_bin_width) #To put in energy densities
        u_lum = np.divide(u_lum,u_bin_width) #To put in energy densities

        #if i%5 == 0:
        #    print np.sum(e_lum)
        #e_lum_all.append(e_lum.tolist())
        #a_lum_all.append(a_lum.tolist())
        #u_lum_all.append(u_lum.tolist())


        savefileprefix = models[iii][0] + '/'
        fileprefix = 'nu_spectra.t_'
        endfix      = '.dat'

        i_eff = i - models[iii][1]
        if i_eff >=0:
            np.savetxt(savefileprefix + fileprefix + str(i_eff).rjust(3, '0') + endfix,np.transpose([e_lum,a_lum,u_lum]),header = '\n'.join(['t = ' + str(i_eff).rjust(3, '0') + ' ms, neutrino energy spectrum for the ' + models[iii][0] + ' model','nu_e                   anti-nu_e                nux']))
        else:
            np.savetxt(savefileprefix + fileprefix + 'n' + str(abs(i_eff)).rjust(3, '0') + endfix,np.transpose([e_lum,a_lum,u_lum]),header = '\n'.join(['t = -' + str(abs(i_eff)).rjust(3, '0') + ' ms, neutrino energy spectrum for the ' + models[iii][0] + ' model','nu_e                   anti-nu_e              nux']))

                       
    xt = np.arange(start_time_pre,end_time_post+1)
    pp.plot(xt,e_lum_all,color='red',label='e')
    pp.plot(xt,a_lum_all,color='blue',label='a')
    pp.plot(xt,u_lum_all,color='green',label='u')
    pp.xlabel("Time (ms)")
    pp.ylabel("Luminosity 10$^{51}$ erg/s")
    pp.legend(loc='best')
    pp.savefig(models[iii][0] + ".lum.pdf")
    pp.close()



