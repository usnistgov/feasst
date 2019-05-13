# pyfeasst - collection of python utilities

import os

class cd:
    """Context manager for changing the current working directory
       Thanks to Brian Hunt on Stack Overflow
       https://stackoverflow.com/questions/431684/how-do-i-change-directory-cd-in-python/13197763#13197763 """
    def __init__(self, newPath):
        self._newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self._savedPath = os.getcwd()
        os.chdir(self._newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self._savedPath)

import subprocess

def bash_command(cmd):
    subprocess.call(cmd, shell=True, executable='/bin/bash')
    #subprocess.Popen(['/bin/bash', '-c', cmd])
