"""
Block average analysis for correlated data

In this example, a random walk is performed in one dimension with periodic boundaries.
The step size is chosen randomly uniform in [-1, 1], and the repeating cell is 25.
Use the blocking method (http://dx.doi.org/10.1063/1.457480) to find the standard deviation of the mean from uncorrelated data, and to estimate the correlation time (i.e., the number of steps until the particle could be anywhere in the repeating cell).

For 1e6 steps, the error on the average position is about 0.08 with about 68% confidence, and a correlation time of roughly 2^11 steps correspoding with roughly (25/0.5)^2.
"""

import random
import feasst

pbc_length = 25
step_size = 1
num_steps = int(1e6)

av_position = feasst.Accumulator({'max_block_operations': '20'})
position = random.uniform(-pbc_length/2, pbc_length/2)
for step in range(num_steps):
    position += step_size*random.uniform(-1, 1)
    if position < -pbc_length/2:
        position += pbc_length
    elif position > pbc_length/2:
        position -= pbc_length
    av_position.accumulate(position)
for op in range(av_position.max_block_operations()):
    print(op, av_position.block_stdev(op, 10))
