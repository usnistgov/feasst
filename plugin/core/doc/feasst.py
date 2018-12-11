import pygraphviz as pgv
B=pgv.AGraph('feasst.dot') # create a new graph from file
B.layout('dot') # layout with default (neato)
#B.layout() # layout with default (neato)
B.draw('feasst.png') # draw png
