
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['OpenSans']})
rc('text', usetex=True)

shift = -4
plt.gca().add_patch(plt.Rectangle(xy=(-2+shift, -2), width=4, height=4, fc='none', ec='black', lw=2))

jx=0.5+shift
jy=-0.5
plt.gca().add_patch(plt.Circle((jx, jy),1.3, fc=(0, 1, 0, 0.5)))
plt.gca().add_patch(plt.Circle((jx, jy),0.5, fc='white'))
plt.gca().add_patch(plt.Circle((jx, jy),0.25, fc='blue'))
plt.text(0.5+shift, jy, 'A', horizontalalignment='center', verticalalignment='center', fontsize=12, color='white')

plt.text(-0.5+shift-0.9/np.sqrt(2), -0.5+0.9/np.sqrt(2), r'$v_{AV}$', horizontalalignment='center', verticalalignment='center', fontsize=18, color='black')

plt.gca().add_patch(plt.Circle((-0.5+shift+1/1.2/np.sqrt(2), -0.5-1/1.2/np.sqrt(2)),0.25, fc='black'))
plt.text(-0.5+shift+1/1.2/np.sqrt(2), -0.5-1/1.2/np.sqrt(2), 'B', horizontalalignment='center', verticalalignment='center', fontsize=12, color='white')

plt.gca().add_patch(plt.Circle((0.75+shift, 1.5),0.25, fc='black'))
plt.text(0.75+shift, 1.5, 'B', horizontalalignment='center', verticalalignment='center', fontsize=12, color='white')

plt.gca().add_patch(plt.Circle((0.5+shift+1/1.8/np.sqrt(2), -0.5+1/1.8/np.sqrt(2)),0.25, fc='black'))
plt.text(0.5+shift+1/1.8/np.sqrt(2), -0.5+1/1.8/np.sqrt(2), 'B', horizontalalignment='center', verticalalignment='center', fontsize=12, color='white')

plt.gca().annotate("", xy=(-1.8, 1), xytext=(1.8, 1), arrowprops=dict(arrowstyle="<-", lw=2, color='red'))
plt.text(0, 1.5, r'A$\rightarrow$C', horizontalalignment='center', verticalalignment='center', fontsize=16)
plt.text(0, 0.5, r'B$\rightarrow$D', horizontalalignment='center', verticalalignment='center', fontsize=16)
plt.gca().annotate("", xy=(-1.8, -1), xytext=(1.8, -1), arrowprops=dict(arrowstyle="->", lw=2, color='blue'))
plt.text(0, -0.5, r'A$\leftarrow$C', horizontalalignment='center', verticalalignment='center', fontsize=16)
plt.text(0, -1.5, r'B$\leftarrow$D', horizontalalignment='center', verticalalignment='center', fontsize=16)

shift = 4
jx=0.5+shift
jy=-0.5
plt.gca().add_patch(plt.Rectangle(xy=(-2+shift, -2), width=4, height=4, fc='none', ec='black', lw=2))
plt.gca().add_patch(plt.Circle((jx, jy),1.3, fc=(0, 1, 0, 0.5)))
plt.gca().add_patch(plt.Circle((jx, jy),0.5, fc='white'))
plt.gca().add_patch(plt.Circle((jx, jy),0.25, fc='red'))
plt.text(jx, jy, 'C', horizontalalignment='center', verticalalignment='center', fontsize=12, color='white')

plt.gca().add_patch(plt.Circle((0.75+shift, 1.5),0.25, fc='black'))
plt.text(0.75+shift, 1.5, 'B', horizontalalignment='center', verticalalignment='center', fontsize=12, color='white')

plt.gca().add_patch(plt.Circle((0.5+shift+1/1.8/np.sqrt(2), -0.5+1/1.8/np.sqrt(2)),0.25, fc='purple'))
plt.text(0.5+shift+1/1.8/np.sqrt(2), -0.5+1/1.8/np.sqrt(2), 'D', horizontalalignment='center', verticalalignment='center', fontsize=12, color='white')

plt.gca().add_patch(plt.Circle((-0.5+shift+1/1.2/np.sqrt(2), -0.5-1/1.2/np.sqrt(2)),0.25, fc='black'))
plt.text(-0.5+shift+1/1.2/np.sqrt(2), -0.5-1/1.2/np.sqrt(2), 'B', horizontalalignment='center', verticalalignment='center', fontsize=12, color='white')

plt.text(-0.5+shift-0.9/np.sqrt(2), -0.5+0.9/np.sqrt(2), r'$v_{AV}$', horizontalalignment='center', verticalalignment='center', fontsize=18, color='black')


plt.axis('scaled')
plt.axis('off')
#plt.show()
plt.savefig('rxvb.png', bbox_inches='tight', transparent=True, dpi=300)

