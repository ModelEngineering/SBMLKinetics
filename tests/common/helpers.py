from src.common import constants as cn
from src.common import simple_sbml

import libsbml
import os


DIR = os.path.dirname(os.path.abspath(__file__))
TEST_PATH = os.path.join(DIR, "test_file.xml")
TEST_PATH_1 = os.path.join(DIR, "BIOMD0000000006.xml")
TEST_PATH_56 = os.path.join(DIR, "BIOMD0000000056.xml")

def getSimple():
  return simple_sbml.SimpleSBML(TEST_PATH)

def getSimple_BIOMD6():
  return simple_sbml.SimpleSBML(TEST_PATH_1)

def getSimple_BIOMD56():
  return simple_sbml.SimpleSBML(TEST_PATH_56)
