import random
import feasst

pbc_length = 25
step_size = 1
num_steps = int(1e6)

av_position = feasst.Accumulator({'max_block_operations': '20'})
#fst.args({"max_block_operations": "20"}))
#av_position2 = fst.MakeAccumulator()
#rng = fst.MakeRandomMT19937()
position = random.uniform(-pbc_length/2, pbc_length/2)
for step in range(num_steps):
    position += step_size*random.uniform(-1, 1)
    if position < -pbc_length/2:
        position += pbc_length
    elif position > pbc_length/2:
        position -= pbc_length
    av_position.accumulate(position)
#    av_position2.accumulate(position)
for op in range(av_position.max_block_operations()):
    print(op, av_position.block_stdev(op, 10))
#self.assertAlmostEqual(av_position2.block_stdev(), 0.08, delta=0.05)
