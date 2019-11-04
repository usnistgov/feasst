import pygraphviz as pgv
B=pgv.AGraph('configuration.dot') # create a new graph from file
B.layout('dot') # layout with default (neato)
#B.layout() # layout with default (neato)
#B.draw('configuration.png') # draw png
#B.draw('configuration.eps') # draw eps
B.draw('configuration.svg') # draw svg
