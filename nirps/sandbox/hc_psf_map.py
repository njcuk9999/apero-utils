from astropy.io import fits
import numpy as np
from photutils import DAOStarFinder, morphology,data_properties
import matplotlib.pyplot as plt
import os
from astropy.table import Table

plt.ioff()


file_comp = '../simu_nirps/simu_HA/20200402085619_ramp_003_HA_HC_HC_pp.fits'

do_plot_comp = False

focus_range = np.array([-1000,-750,-500,-250,0,250,500,750,1000,np.nan])

keys = ['semimajor_axis_sigma', 'sharpness', 'orientation', 'elongation','rvcontent']
units = ['pixel', 'unitless', 'deg', 'semi-major/semi-minor','Q/Q$_0$']

# axis 0 = focus position
# axis 1 = keys
# axis 2 = region on array
# axis 3 = value/rms around value
u = np.zeros([len(focus_range), len(keys), 5,2])

for ifocus_ref in range(len(focus_range)):

    focus_ref = focus_range[ifocus_ref]

    if np.isfinite(focus_ref):
        file = 'psf_map'+str(int(focus_ref))+'.fits'
    else:
        file = file_comp



    tbl_name = file.split('.fits')[0]+'_psf.csv'


    if os.path.isfile(tbl_name) == False:
        print('reading file : '+file)
        im = fits.getdata(file)
        # constructing a mask for the edges of the array
        mask = np.ones_like(im)
        mask[0:10,:] = 0
        mask[-10:,:] = 0
        mask[:,0:10] = 0
        mask[:,-10:] = 0

        print('finding 1-sigma std')
        std = np.nanpercentile(np.abs(im),68)

        print('finding PSFs in the image, cutoff at 3 sigma and 5 std from noise')
        daofind = DAOStarFinder(fwhm=3.0, threshold=25.*std)
        sources = daofind(im*mask)
        sources['rvcontent'] = np.zeros_like(sources,dtype = float)

        print('preparing output table to accept values')
        # source properties that will be kept in the master table
        columns = ['elongation', 'orientation','semiminor_axis_sigma','semimajor_axis_sigma','gini']

        # preparing the output table for new columns
        for col in columns:
            sources[col] = np.zeros_like(sources, dtype=float)


        for ii in range(len(sources)):
            data = im[int(sources['ycentroid'][ii])-6:int(sources['ycentroid'][ii])+7,int(sources['xcentroid'][ii])-6:int(sources['xcentroid'][ii])+7]

            profil = np.sum(data/np.sum(data),axis=1)
            rv_content = np.sum(np.gradient(profil) ** 2)

            sources['rvcontent'][ii] = rv_content
            tbl = data_properties(data).to_table(columns=columns)

            for col in columns:
                sources[col][ii] =  np.array(tbl[col])[0]
            print('PSF {0} / {1}, elongation = {2}, semi-major axis = {3}'.format(ii,len(sources), np.array(tbl['elongation'])[0], np.array(tbl['semimajor_axis_sigma'])[0] ) )
        sources.write(tbl_name)

    # reading the list of sources
    sources = Table.read(tbl_name)

    x = sources['xcentroid']
    y = sources['ycentroid']

    for ikey in range(len(keys)):
        tmp = sources[keys[ikey]]

        reg1 = (x < 2048) & (y < 2048)
        reg2 = (x < 2048) & (y > 2048)
        reg3 = (x > 2048) & (y < 2048)
        reg4 = (x > 2048) & (y > 2048)

        u[ifocus_ref,ikey,0,0] = np.nanmedian(tmp[reg1])
        u[ifocus_ref,ikey,1,0] = np.nanmedian(tmp[reg2])
        u[ifocus_ref,ikey,2,0] = np.nanmedian(tmp[reg3])
        u[ifocus_ref,ikey,3,0] = np.nanmedian(tmp[reg4])

        u[ifocus_ref,ikey,0,1] = np.nanpercentile(np.abs(tmp[reg1] - np.nanmedian(tmp[reg1])),68)
        u[ifocus_ref,ikey,1,1] = np.nanpercentile(np.abs(tmp[reg2] - np.nanmedian(tmp[reg2])),68)
        u[ifocus_ref,ikey,2,1] = np.nanpercentile(np.abs(tmp[reg3] - np.nanmedian(tmp[reg3])),68)
        u[ifocus_ref,ikey,3,1] = np.nanpercentile(np.abs(tmp[reg4] - np.nanmedian(tmp[reg4])),68)


fig,ax = plt.subplots(nrows = 2, ncols = 2, figsize = [12,12])

colors = ['ro-','go-','bo-','co-']
labels = ['x<2048, y<2048','x<2048, y>2048','x>2048, y<2048', 'x>2048, y>2048']

facecolor = ['red','green','blue','cyan']
for i in range(4):
    ax[0,0].fill_between(np.array(focus_range),u[:,0,i,0]-u[:,0,i,1]/2,u[:,0,i,0]+u[:,0,i,1]/2,facecolor = facecolor[i], alpha = 0.3,label = labels[i])

    if do_plot_comp:
        g = np.isnan(focus_range)
        ax[0,0].errorbar(np.zeros_like(g[g],dtype = float)+i*50,u[g,0,i,0],yerr = u[g,0,i,1]/2,color = facecolor[i],ecolor = facecolor[i], barsabove = True,fmt = 'o')

ax[0,0].set(title = keys[0],xlabel  = 'focus ($\mu$m)',ylabel = units[0])

for i in range(4):
    ax[0,1].fill_between(np.array(focus_range),u[:,1,i,0]-u[:,1,i,1]/2,u[:,1,i,0]+u[:,0,i,1]/2,facecolor = facecolor[i], alpha = 0.3)

    if do_plot_comp:
        g = np.isnan(focus_range)
        ax[0,1].errorbar(np.zeros_like(g[g],dtype = float)+i*50,u[g,1,i,0],yerr = u[g,1,i,1]/2,color = facecolor[i],ecolor = facecolor[i], barsabove = True,fmt = 'o')


ax[0,1].set(title = keys[1],xlabel  = 'focus ($\mu$m)',ylabel = units[1])

for i in range(4):
    ax[1,0].fill_between(np.array(focus_range),u[:,2,i,0]-u[:,2,i,1]/2,u[:,2,i,0]+u[:,0,i,1]/2,facecolor = facecolor[i], alpha = 0.3)

    if do_plot_comp:
        g = np.isnan(focus_range)
        ax[1,0].errorbar(np.zeros_like(g[g],dtype = float)+i*50,u[g,2,i,0],yerr = u[g,2,i,1]/2,color = facecolor[i],ecolor = facecolor[i], barsabove = True,fmt = 'o')


ax[1,0].set(title = keys[2],xlabel  = 'focus ($\mu$m)',ylabel = units[2])

for i in range(4):
    ax[1,1].fill_between(np.array(focus_range),u[:,3,i,0]-u[:,3,i,1]/2,u[:,3,i,0]+u[:,0,i,1]/2,facecolor = facecolor[i], alpha = 0.3)

    if do_plot_comp:
        g = np.isnan(focus_range)
        ax[1,1].errorbar(np.zeros_like(g[g],dtype = float)+i*50,u[g,3,i,0],yerr = u[g,3,i,1]/2,color = facecolor[i],ecolor = facecolor[i], barsabove = True,fmt = 'o')


ax[1,1].set(title = keys[3],xlabel  = 'focus ($\mu$m)',ylabel = units[3])

#ax[1].plot(focus_ref,np.nanmed
ax[0,0].legend()
plt.tight_layout()
plt.savefig('psf_morphology2.png')
plt.show(block = True)


fig,ax = plt.subplots(nrows = 1, ncols = 1, figsize = [7,7])
for i in range(4):

    g = np.isfinite(focus_range)
    norm = np.max(u[g,4,i,0])
    ax.fill_between(np.array(focus_range),u[:,4,i,0]/norm-u[:,4,i,1]/2/norm,u[:,4,i,0]/norm+u[:,4,i,1]/2/norm,facecolor = facecolor[i], alpha = 0.5,label = labels[i])

    if do_plot_comp:
        g = np.isnan(focus_range)
        ax.errorbar(np.zeros_like(g[g],dtype = float)+i*50,u[g,4,i,0]/norm,yerr = u[g,4,i,1]/2/norm,color = facecolor[i],ecolor = facecolor[i], barsabove = True,fmt = 'o')


ax.set(title = keys[4],xlabel  = 'focus ($\mu$m)',ylabel = units[4])
ax.legend()
plt.tight_layout()
plt.savefig('rvcontent.png')
plt.show()


stop

fig_name =  file.split('.fits')[0]+'_vector_field.png'


#if os.path.isfile(fig_name) == False:
if True:
    g = (sources['elongation']-1) < np.nanmedian(sources['elongation']-1)*5
    sources = sources[g]

    x = sources['xcentroid']
    y = sources['ycentroid']


    u = (sources['elongation']-1)*np.cos(sources['orientation']/180*np.pi  )
    v = (sources['elongation']-1)*np.sin(sources['orientation']/180*np.pi  )

    plt.quiver(x,y,u,v,scale = 1)
    plt.xlabel('x field position')
    plt.ylabel('y field position')

    plt.savefig(fig_name,figsize = [12,5.5])

    plt.show(block = True)
