import feasst

def add(monte_carlo, steps_per, proc="", log="log.txt", tolerance=1e-8):
    """Add Log, Movie, CheckEnergy and Tuner to monte_carlo

    steps_per -- perform analysis every this many steps
    proc -- append movie file name with proc
    log -- set the log file (default empty: standard output)
    """
    monte_carlo.add(feasst.MakeLog(feasst.args(
        {"steps_per" : str(steps_per),
         "file_name": str(log),
         "clear_file": "true"})))
    monte_carlo.add(feasst.MakeMovie(feasst.args(
        {"steps_per" : str(steps_per),
         "file_name" : "movie"+str(proc)+".xyz"})))
    monte_carlo.add(feasst.MakeCheckEnergy(feasst.args(
        {"steps_per" : str(steps_per),
         "tolerance" : str(tolerance)})))
    monte_carlo.add(feasst.MakeTuner(feasst.args(
        {"steps_per" : str(steps_per)})))
