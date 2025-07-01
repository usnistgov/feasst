import feasst
mc = feasst.MonteCarlo()
text_input = """RandomMT19937 seed=123
Configuration cubic_side_length=8 particle_type=lj:/feasst/particle/lj.txt"""
for line in text_input.split('\n'):
    feasst.parse(mc, line)

# alternatively, line by line
feasst.parse(mc, 'RandomMT19937 seed=123')
feasst.parse(mc, 'Configuration cubic_side_length=8 particle_type=lj:/feasst/particle/lj.txt add_num_lj_particles=1')

# Test that the first config has no particles and the second one has one at the origin
try:
    mc.configuration(0).particle(0).site(0).position(0)
    assert False
except RuntimeError:
    assert True
assert mc.configuration(1).particle(0).site(0).position(0) == 0.
