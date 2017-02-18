"""
Use this for exploring fft and filters
"""
import numpy as np
import os
from copy import deepcopy
from matplotlib import pylab as plt
from geoNet import utils, putils

velloc="/".join([os.getcwd(), "Vol1", "data", "velBB"])
accloc="/".join([os.getcwd(), "Vol1", "data", "accBB"])

velMECS = utils.get_stat_data(velloc, 'MECS')
accMECS = utils.get_stat_data(accloc, 'MECS')

#vel_data = utils.int_stat_data(accMECS)
#acc_data = utils.diff_stat_data(velMECS)

g=981.

accMECS_fft = utils.fft_stat_data(accMECS)
fig, ax = plt.subplots(3,1)
#ax[0].plot(accMECS_fft['f'], np.abs(accMECS_fft['000']))
#ax[1].plot(accMECS_fft['f'], np.abs(accMECS_fft['090']))
#ax[2].plot(accMECS_fft['f'], np.abs(accMECS_fft['ver']))
ax[0].plot(accMECS_fft['f'], accMECS_fft['000'])
ax[1].plot(accMECS_fft['f'], accMECS_fft['090'])
ax[2].plot(accMECS_fft['f'], accMECS_fft['ver'])

accMECS_LF = deepcopy(accMECS)
accMECS_HF = deepcopy(accMECS)
tf_LF = utils.filt_stat_data(accMECS_LF,  1., 'lowpass', output='sos', order=4,
                   worN=accMECS_LF['000'].size)
tf_HF = utils.filt_stat_data(accMECS_HF, 1., 'highpass', output='sos', order=4,
                   worN=accMECS_HF['000'].size)
LF_fft = utils.fft_stat_data(accMECS_LF)
HF_fft = utils.fft_stat_data(accMECS_HF)

fig_tf, ax_tf = plt.subplots()
ax_tf.plot(tf_LF['w'], np.abs(tf_LF['h']), label='LF')
ax_tf.plot(tf_HF['w'], np.abs(tf_HF['h']), label='HF')
ifmax = LF_fft['f'].size/2
LF = np.abs(LF_fft['ver'][1:ifmax])
HF = np.abs(HF_fft['ver'][1:ifmax])
#LF = (LF_fft['000'][1:ifmax])
#HF = (HF_fft['000'][1:ifmax])
ax_tf.plot(LF_fft['f'][1:ifmax], LF/LF.max(), label='000-LF')
ax_tf.plot(HF_fft['f'][1:ifmax], HF/HF.max(), label='000-HF')
ax_tf.vlines(1., 0., 1.)
ax_tf.legend(loc='best')
ax_tf.set_xscale('log')
plt.show()
#plt.show(block=False)
