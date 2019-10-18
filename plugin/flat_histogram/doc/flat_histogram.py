import pygraphviz as pgv
B=pgv.AGraph('flat_histogram.dot') # create a new graph from file
B.layout('dot') # layout with default (neato)
B.draw('flat_histogram.svg') # draw svg
