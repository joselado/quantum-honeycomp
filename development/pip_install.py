#!/usr/bin/python

# this script will install all the necessary python packages using pip

import os

def install(name):
  os.system("pip install "+name)

packages = ["scipy","numpy","matplotlib","multiprocess"]

for p in packages: install(p)

# just in case
os.system("python compile_fortran.py")  # compile fortran libraries

