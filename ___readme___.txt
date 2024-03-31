
 .oooooo..o  .oooooo..o oooooooooo.                          
d8P'    `Y8 d8P'    `Y8 `888'   `Y8b                         
Y88bo.      Y88bo.       888      888 oo.ooooo.  oooo    ooo 
 `"Y8888o.   `"Y8888o.   888      888  888' `88b  `88.  .8'  
     `"Y88b      `"Y88b  888      888  888   888   `88..8'   
oo     .d8P oo     .d8P  888     d88'  888   888    `888'    
8""88888P'  8""88888P'  o888bood8P'    888bod8P'     .8'     
                                       888       .o..P'      
                                      o888o      `Y8P'       
                                                             
DESCRIPTION

A collection of python packages for the structural and dynamic analysis of dynamical systems..


CONTRIBUTORS

The Structural and Stochastic Dynamics research group at the University of Liège
2023-... Vincent Denoël


STABILIZED VERSIONS
Date     ┃  Release  ┃  Comments
Mar. 23  ┃  V.1.0.0  ┃  Development of the main folder arrangement


CONTENTS
dyn    # structural dynamics : integrators, ...
misc   #
signal # some signal processing techniques : first passage time, ...
stat   # advanced statistics methods
tools  # miscellaneous tools for plotting, time-freq analysis
wind   # wind engineering : buffeting, turbulence simulator, flutter analysis, ...

INSTALLATION

In the VS code terminal, create a new environment, called SSDpy
python3.9 -m venv SSDpy


MAKE THE PACKAGE ACCESSIBLE IN ANY PROJECT : OPTION 1

Run SSDpy_update, which will copy the files from your local folder to the site-packages/SSD folder


MAKE THE PACKAGE ACCESSIBLE IN ANY PROJECT : OPTION 2 (use .pth file >> seem to not work...???)

You need to move and personalize a .pth file that contains the directory to search for, and place it in the {venv-root}/lib/{python-version}/site-packages directory.
see https://docs.python.org/3/install/index.html#modifying-python-s-search-path
"The most convenient way [to add a directory] is to add a path configuration file to a directory that’s already on Python’s path ... Path configuration files have an extension of .pth, and each line must contain a single path that will be appended to sys.path."

See file
SSDpy.pth
in /opt/anaconda3/envs/YOUR_ENV/lib/pythonX.X/site-packages

This should should contain the list of files to the permanent folder of packages that need to be accessed. For example,
/Users/vincentDenoel/Dropbox/001_ULg/001_Recherches/MyPythonToolboxes/SSDpy
/Users/vincentDenoel/Dropbox/001_ULg/001_Recherches/MyPythonToolboxes/SSDpy/dyn
/Users/vincentDenoel/Dropbox/001_ULg/001_Recherches/MyPythonToolboxes/SSDpy/stat
/Users/vincentDenoel/Dropbox/001_ULg/001_Recherches/MyPythonToolboxes/SSDpy/wind
/Users/vincentDenoel/Dropbox/001_ULg/001_Recherches/MyPythonToolboxes/SSDpy/wind/tornado_model

You can check the system path with the following lines of code (in Terminal):
import sys; print(sys.path)


