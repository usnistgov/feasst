import feasst
mc = feasst.MonteCarlo()
text_input = """RandomMT19937 seed 123
Configuration cubic_side_length 8 particle_type0 /feasst/particle/lj.fstprt"""
for line in text_input.split('\n'):
    feasst.parse(mc, line)

# alternatively, line by line
feasst.parse(mc, 'RandomMT19937 seed 123')
feasst.parse(mc, 'Configuration cubic_side_length 8 particle_type0 /feasst/particle/lj.fstprt')
