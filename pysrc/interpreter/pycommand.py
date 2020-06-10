import subprocess
import os
import sys
import platform

def correct_python():
    """CHeck if a suitable Python is installed"""
    try:
        import scipy
        import numpy
        import multiprocessing
        out,err = subprocess.Popen(['python', '--version'],
               stdout=subprocess.PIPE,
               stderr=subprocess.STDOUT).communicate()
    except: out = ""
    # check if JUlia has the correct version
    return "Python 3.7" in str(out)

def install_python():
    """Install a correct Python distribution"""
    if correct_python(): 
        print("Found a correct python distribution")
        return # nothing to do
    pwd = os.getcwd() # get the current directory
    dirname = os.path.dirname(os.path.realpath(__file__)) # this directory
    os.system("rm -rf "+dirname+"/python_interpreter") # remove the subfolder
    os.system("mkdir "+dirname+"/python_interpreter") # create the subfolder
    os.chdir(dirname+"/python_interpreter") # go to this directory
    pypath = dirname+"/python_interpreter/python3" # path to Python
    if platform.system()=="Linux":
        anapath = "https://repo.anaconda.com/archive/Anaconda3-2020.02-Linux-x86_64.sh"
        os.system("wget "+anapath) # download Anaconda Python
        anafile = "Anaconda3-2020.02-Linux-x86_64.sh"
    elif platform.system()=="Windows":
        print("Not compatible with Windows yet")
        exit()
    else: # for Mac
        anapath = "https://repo.anaconda.com/archive/Anaconda3-2020.02-MacOSX-x86_64.sh"
        os.system("curl -LO "+anapath) # download Anaconda Python
        anafile = "Anaconda3-2020.02-MacOSX-x86_64.sh " # file to install
    os.system("bash "+anafile+" -b -p "+pypath) # install anaconda
    os.system("rm "+anafile) # remove the installer
    install_dependencies() # install the dependencies

def install_dependencies():
    pip = get_pip() # pip command
    for l in ["mayavi","multiprocessing","numba"]:
        os.system(pip+" install "+l) # install this library



def get_python():
  """Return the path for Anaconda Python, which has pyqt by default"""
  if correct_python(): return "python" # default python command
  else:
    dirname = os.path.dirname(os.path.realpath(__file__)) # this directory
    return dirname +"/python_interpreter/python3/bin/python"


def get_pip():
  """Return the path for Anaconda Python, which has pyqt by default"""
  if correct_python(): return "python" # default python command
  else:
    dirname = os.path.dirname(os.path.realpath(__file__)) # this directory
    return dirname +"/python_interpreter/python3/bin/pip"


def add_to_path():
    """Add quantum honeycomp to the PATH"""
    out = os.environ["SHELL"]
    home = os.environ["HOME"]
    if out=="/bin/bash":
        if platform.system()=="Linux":  rcfile = home+"/.bashrc"
        else: rcfile = home+"/.bash_profile"
    elif out=="/bin/zsh": 
        rcfile = home+"/.zshrc"
    qhpath = os.path.dirname(os.path.realpath(__file__))+"/../../bin"
    ls = open(rcfile,"r").read()
    addrc = "\nexport PATH=\""+qhpath+"\":$PATH\n"
    open(rcfile,"w").write(ls+addrc)
    print(rcfile,qhpath)





