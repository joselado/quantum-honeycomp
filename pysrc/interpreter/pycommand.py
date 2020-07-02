import subprocess
import os
import sys
import platform
from pathlib import Path


def write_python_exec(name):
    """Create pythoninterpreter.py which contains the path to the target Python executable"""
    os.environ["PYTHON_EXEC"] = name
    path = Path('pythoninterpreter.py')
    with open(path, mode='wt') as config:
        config.write("mainpython = r\""+os.environ["PYTHON_EXEC"]+"\"")


def correct_python(install=False):
    """Check if a suitable Python is installed"""
    try:
        import PyQt5
        import scipy
        import numpy
        import numba
        import matplotlib
        return True
    except:
        if install:
            install_dependencies(executable=sys.executable)  # try to install
            return correct_python(install=False)  # try again
        return False


def install_python():
    """Install a correct Python distribution"""
    if correct_python(install=True):  # The current one is the right Python
        write_python_exec(sys.executable)  # write this Python distribution
        print("Found a correct python distribution")
        return  # nothing to do
    dirname = os.path.dirname(os.path.realpath(__file__))  # this directory
    os.system("rm -rf "+dirname+"/python_interpreter")  # remove the subfolder
    os.system("mkdir "+dirname+"/python_interpreter")  # create the subfolder
    os.chdir(dirname+"/python_interpreter")  # go to this directory
    pypath = dirname+"/python_interpreter/python3"  # path to Python
    if platform.system() == "Linux":
        anapath = "https://repo.anaconda.com/archive/Anaconda3-2020.02-Linux-x86_64.sh"
        os.system("wget "+anapath)  # download Anaconda Python
        anafile = "Anaconda3-2020.02-Linux-x86_64.sh"
    elif platform.system() == "Windows":
        print("Please download and install Anaconda")
        exit()
    else:  # for Mac
        anapath = "https://repo.anaconda.com/archive/Anaconda3-2020.02-MacOSX-x86_64.sh"
        os.system("curl -LO "+anapath)  # download Anaconda Python
        anafile = "Anaconda3-2020.02-MacOSX-x86_64.sh "  # file to install
    os.system("bash "+anafile+" -b -p "+pypath)  # install anaconda
    os.system("rm "+anafile)  # remove the installer
    # now get the executable
    dirname = os.path.dirname(os.path.realpath(__file__))  # this directory
    pyint = dirname + "/python_interpreter/python3/bin/python3"  # local one
    write_python_exec(pyint)  # write this Python distribution


def install_dependencies(executable=None):
    if executable is None:
        executable = get_python()
    for l in ["mayavi", "numba", "scipy", "numpy", "matplotlib"]:
        try:
            install_package(l, executable=executable)
        except:
            pass


def get_python():
    """Return the path for Anaconda Python, which has pyqt by default"""
    try:
        #from .pythoninterpreter import mainpython
        mainpython = sys.executable
        #mainpython = os.environ["PYTHON_EXEC"]
        print("Using the interpreter", mainpython)
        return mainpython
    except:
        print("No python interpreter found, exiting")
        exit()


def add_to_path():
    """Add quantum honeycomp to the PATH"""
    out = os.environ["SHELL"]
    home = os.environ["HOME"]
    if out == "/bin/bash":
        if platform.system() == "Linux":
            rcfile = home+"/.bashrc"
        else:
            rcfile = home+"/.bash_profile"
    elif out == "/bin/zsh":
        rcfile = home+"/.zshrc"
    qhpath = os.path.dirname(os.path.realpath(__file__))+"/../../bin"
    try:
        ls = open(rcfile, "r").read()  # if the file exists
    except:
        ls = ""  # otherwise
    addrc = "alias quantum-honeycomp=\"" + get_python() + " "+qhpath + \
        "/quantum-honeycomp\""
#    addrc = "\nexport PATH=\""+qhpath+"\":$PATH\n"
    open(rcfile, "w").write(ls+addrc)  # add to the bash


def install_package(package, executable=None):
    if executable is None:
        executable = sys.executable
    subprocess.check_call([executable, "-m", "pip", "install", package])


def run_qh():
    """Run Quantum Honeycomp"""
    qhpath = Path(__file__).resolve().parent.parent.parent
    qh = qhpath.joinpath('bin', "quantum-honeycomp.py")
    subprocess.run(["python", str(qh)])
