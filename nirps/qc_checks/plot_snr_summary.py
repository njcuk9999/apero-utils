from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time

"""
Tool to plot a few things from the SNR summary table
"""

tbl = Table.read('all_snr_apero.csv')

plt.plot(tbl['SNR_Y'],tbl['SNR_Y']/tbl['SNR_Y_APERO'],'b.',alpha = 0.3, label = 'Y band')
plt.plot(tbl['SNR_J'],tbl['SNR_J']/tbl['SNR_J_APERO'],'g.',alpha = 0.3, label = 'J band')
plt.plot(tbl['SNR_H'],tbl['SNR_H']/tbl['SNR_H_APERO'],'r.',alpha = 0.3, label = 'H band')
plt.xlabel('difference SNR')
plt.ylabel('difference SNR /APERO SNR')
plt.ylim(0,2)
plt.xlim(10,800)
plt.xscale('log')
plt.legend()
plt.grid(color = 'lightgrey', alpha = 0.8)
plt.savefig('snr_comparison.pdf',dpi=300,bbox_inches='tight')
plt.show()


tbl['plot_date'] = Time(tbl['MJDMID'],format = 'mjd').plot_date

t1 = Time('2023-04-23').plot_date
t2 = Time('2023-05-10').plot_date

tbl = tbl[(tbl['plot_date']>t1) & (tbl['plot_date']<t2)]

plt.figure(figsize = (18,8))
for uobj in np.unique(tbl['OBJECT']):
    g = np.where(tbl['OBJECT']==uobj)[0]
    plt.plot_date(tbl['plot_date'][g],tbl['OVERHEAD'][g],'o',alpha = 0.5)
plt.xlabel('date')
plt.ylabel('Frame-to-frame overhead [s]')
plt.ylim(0,200)
plt.grid(color = 'lightgrey', alpha = 0.8)

plt.xlim(t1,t2)
#plt.plot([t1,t2],[11.5,11.5],'r--',label = 'nominal overhead')
#plt.legend()
plt.tight_layout()
plt.savefig('frame_to_frame_overhead.pdf',dpi=300,bbox_inches='tight')
plt.show()
