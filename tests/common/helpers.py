from src.common import constants as cn
from src.common.simple_sbml import SimpleSBML

import libsbml
import os


DIR = os.path.dirname(os.path.abspath(__file__))
TEST_PATH = os.path.join(DIR, "test_file.xml")


def getSimple():
  return SimpleSBML(TEST_PATH)
