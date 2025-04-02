import pkgutil
import importlib
import os

from .constitutive.constlaws import boucWen

from .dyn.EOMsolvers import newmark, NewmarkMDDL
from .dyn.log_decrement import log_decrement
from .dyn.damping import RayleighDamping

from .gcdc import gcdc_tools

from .num import winkler

from .num.nonlin import example_structures
from .num.nonlin import nonlin_FEM

from .plot import my_colormap as cmp
from .plot.plot_tf import plot_tf
from .plot.plot_settings import plot_settings

from .signal.envelopes import envelopes
from .signal.findpeaks import findpeaks
from .signal.phase import phase
from .signal.histv import histv
from .signal.vpsd import vpsd

from .tools.clipfig import clipfig
from .tools.matlab import read_mat_file
from .tools.saveas import saveas

from .wind.bamm import bamm
from .wind.gpvm import gpvm
#from .wind.tornado_model.tornado_model import tornado_model
from .wind.wind_psd import wind_psd
#from .wind.viv import viv

'''
# Dynamically import all modules from all subdirectories in SSDpy
for loader, module_name, is_pkg in pkgutil.walk_packages(__path__, __name__ + "."):
    module = importlib.import_module(module_name)
    # Optionally, add specific attributes from these modules to the SSDpy namespace
    globals().update({name: getattr(module, name) for name in dir(module) if not name.startswith("_")})
'''