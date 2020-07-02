import os
import sys
import platform
import argparse
import subprocess
from pathlib import Path


parser = argparse.ArgumentParser()
parser.add_argument("--utility", default="",
                    help='Just use a utility and exit')
parser.add_argument("--python", default="",
                    help='Default Python interpreter')
args = parser.parse_args()  # get the arguments


def get_qhroot():
    """Gets the root path of quantum honeycomp"""
    return Path(__file__).resolve().parent.parent


qhroot = get_qhroot()  # get the root path
os.environ["QHROOT"] = str(qhroot)  # create the environmental variable


def get_command(name="python"):
    """Return the path for Anaconda Python, which has pyqt by default"""
    #dirname = os.path.dirname(os.path.realpath(__file__))
    qhroot = get_qhroot()
    path = qhroot.joinpath("pysrc", "interpreter")
    sys.path.insert(1, str(path))  # add this path
    print("Path: ", str(path))
    import pycommand
    return pycommand.get_python()


python = get_command()


if args.python != "":  # non empty string
    path = qhroot.joinpath("pysrc", "interpreter", "pythoninterpreter.py")
    with open(path, mode='wt') as config:
        config.write("mainpython = r\""+args.python+"\"\n")    


# Use a utility and exit
if args.utility != "":  # non empty string
    utility = qhroot.joinpath('utilities')
    subprocess.run(["python", str(utility), args.utility])
    #os.system(python + " "+qhroot+"/utilities/"+args.utility)
    exit()


pyqt = qhroot.joinpath('bin', "versions", "quantum-honeycomp-pyqt")
subprocess.run(["python", str(pyqt)])
