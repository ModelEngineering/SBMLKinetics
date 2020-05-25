'''Provides Information on SBML Kinetics Laws'''


from src.common import constants as cn
from src.common import util

import collections
import numpy as np
import libsbml
import zipfile


class KineticLaw(object):

  def __init__(self, libsbml_kinetics):
    """
    :param libsbml.KineticLaw libsbml_kinetics:
    """
    self.libsbml_kinetics = libsbml_kinetics
    self.formula = self.libsbml_kinetics.getFormula()
