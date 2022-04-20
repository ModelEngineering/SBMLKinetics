from SBMLKinetics.common import constants as cn
from SBMLKinetics.common import simple_sbml

import libsbml
import os


DIR = os.path.dirname(os.path.abspath(__file__))
TEST_PATH = os.path.join(DIR, "test_file.xml")
TEST_PATH_2 = os.path.join(DIR, "BIOMD0000000002.xml")
TEST_PATH_3 = os.path.join(DIR, "BIOMD0000000003.xml")
TEST_PATH_5 = os.path.join(DIR, "BIOMD0000000005.xml")
TEST_PATH_6 = os.path.join(DIR, "BIOMD0000000006.xml")
TEST_PATH_43 = os.path.join(DIR, "BIOMD0000000043.xml")
TEST_PATH_56 = os.path.join(DIR, "BIOMD0000000056.xml")
TEST_PATH_239 = os.path.join(DIR, "BIOMD0000000239.xml")

def getSimple():
  return simple_sbml.SimpleSBML(TEST_PATH)

def getSimple_BIOMD2():
  return simple_sbml.SimpleSBML(TEST_PATH_2)

def getSimple_BIOMD3():
  return simple_sbml.SimpleSBML(TEST_PATH_3)

def getSimple_BIOMD5():
  return simple_sbml.SimpleSBML(TEST_PATH_5)

def getSimple_BIOMD6():
  return simple_sbml.SimpleSBML(TEST_PATH_6)

def getSimple_BIOMD43():
  return simple_sbml.SimpleSBML(TEST_PATH_43)

def getSimple_BIOMD56():
  return simple_sbml.SimpleSBML(TEST_PATH_56)

def getSimple_BIOMD239():
  return simple_sbml.SimpleSBML(TEST_PATH_239)