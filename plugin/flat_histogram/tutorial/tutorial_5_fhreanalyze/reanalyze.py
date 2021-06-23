import argparse
import feasst as fst
import pyfeasst

print(fst.version())
parser = argparse.ArgumentParser()
parser.add_argument("--first_processor", type=int, help="id of the first processor", default=3)
parser.add_argument("--last_processor", type=int, help="id of the last processor", default=3)
args = parser.parse_args()
print("args:", args)

for proc in range(args.first_processor, args.last_processor + 1):
    mc = fst.MonteCarlo().deserialize(pyfeasst.read_checkpoint(
        "checkpoint" + str(proc) + ".fst"))
    fh = fst.FlatHistogram(mc.criteria())
    hist = fh.macrostate().histogram()
    tm = fst.MakeTransitionMatrix(fst.args({"min_sweeps": "1",
                                            "num_blocks": "300"}))
    tm.resize(hist)
    tm.read_dump_file("dump" + str(proc) + ".txt")
    tm.infrequent_update()
    print(tm.write())
    print(tm.write_per_bin_header())
    for bn in range(hist.size()):
        print(str(int(hist.center_of_bin(bn))) + "," + tm.write_per_bin(bn))
    #tm.write(
    #for val in tm.ln_prob().values():
    #    print(val)
