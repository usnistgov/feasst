import pygraphviz as pgv
B=pgv.AGraph('monte_carlo.dot') # create a new graph from file
B.layout('dot') # layout with default (neato)
#B.layout() # layout with default (neato)
#B.draw('monte_carlo.png') # draw png
#B.draw('monte_carlo.eps') # draw eps
B.draw('monte_carlo.svg') # draw svg

B=pgv.AGraph('trial.dot') # create a new graph from file
B.layout('dot') # layout with default (neato)
B.draw('trial.svg') # draw svg
