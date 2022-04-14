import matplotlib.pyplot as plt
#import matplotlib
#matplotlib.rcParams['mathtext.fontset'] = 'custom'
#matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
#matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
#matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
#rcParams['font.family'] = 'sans-serif'

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica"]})

circles=[
[-0.5, 1.2], [0.7, 1.2],
[-1.2, 0], [1.2, 0],
[-0.5, -1.2], [0.7, -1.2],
]

#plt.axes()
for circle in circles:
    plt.gca().add_patch(plt.Circle((circle[0], circle[1]),1, fc='lightgrey'))
for circle in circles:
    plt.gca().add_patch(plt.Circle((circle[0], circle[1]),0.5, fc='grey'))
#plt.gca().set_facecolor('black')
#plt.figure(facecolor='yellow')
plt.gca().add_patch(plt.Circle((0, 0),0.05, fc='black'))
plt.text(0.05, -0.05, r'$x_o$', fontsize=24, color='black')
plt.gca().add_patch(plt.Rectangle(xy=(-0.5, -0.5), width=1, height=1, fc='none', ec='black', lw=2))
#plt.gca().add_patch(plt.Circle((0., 0.),0.5, fc='none', ec='black', linestyle='--'))
plt.gca().add_patch(plt.Circle((0.2, 0.3),0.05, fc='blue'))
plt.text(0.25, 0.25, r'$x_n$', fontsize=24, color='blue')
plt.gca().add_patch(plt.Rectangle(xy=(-0.5+0.2, -0.5+0.3), width=1, height=1, fc='none', ec='blue', lw=2))
#plt.gca().add_patch(plt.Circle((0.2, 0.3),0.5, fc='none', ec='blue', linestyle='--'))
shift=0.02
plt.gca().annotate("", xy=(0.175+shift, 0.275), xytext=(0.025+shift, 0.025*3/2), arrowprops=dict(arrowstyle="->", lw=2))
shift=-0.02
plt.gca().annotate("", xy=(0.175+shift, 0.275), xytext=(0.025+shift, 0.025*3/2), arrowprops=dict(arrowstyle="<-", lw=2, color='blue'))
#plt.gca().annotate("", xy=(0.175, 0.275), xytext=(0.025, 0.025*3/2), arrowprops=dict(arrowstyle="<-", lw=2, color='blue'))

# delta
plt.gca().annotate("", xy=(-0.5, -0.55), xytext=(0.5, -0.55), arrowprops=dict(arrowstyle="<->", lw=2))
plt.text(-0.1, -0.75, r'$2\delta$', fontsize=24, color='black')

plt.axis('scaled')
plt.axis('off')
plt.xlim([-1., 1.])
plt.ylim([-1., 1.])
plt.savefig('dccb_draw.png', bbox_inches='tight')
#plt.show()
