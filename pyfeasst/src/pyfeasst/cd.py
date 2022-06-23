"""This module provides a convenient context manager"""

import os

class cd:
    """
    Context manager for changing the current working directory.
    Author: Brian Hunt on Stack Overflow.
    https://stackoverflow.com/questions/431684/how-do-i-change-directory-cd-in-python/13197763#13197763

    >>> import os
    >>> from pathlib import Path
    >>> import subprocess
    >>> from pyfeasst import cd
    >>> subprocess.call("mkdir -p tmp", shell=True, executable='/bin/bash')
    >>> with cd.cd("tmp"): print(Path(os.getcwd()).name)
    tmp
    """
    def __init__(self, new_path):
        self._new_path = os.path.expanduser(new_path)
        self._saved_path = None

    def __enter__(self):
        self._saved_path = os.getcwd()
        os.chdir(self._new_path)

    def __exit__(self, etype, value, traceback):
        os.chdir(self._saved_path)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
