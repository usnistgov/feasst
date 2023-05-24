import pygraphviz as pgv
A=pgv.AGraph('sim.dot') # create a new graph from file
A.layout('dot') # layout with default (neato)
A.draw('sim.svg') # draw svg

#B=pgv.AGraph('sim2.dot') # create a new graph from file
#B.layout('dot') # layout with default (neato)
#B.draw('sim2.svg') # draw svg
