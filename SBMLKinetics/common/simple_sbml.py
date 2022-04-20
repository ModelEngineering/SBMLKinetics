'''Python representation of an SBML model.'''

from SBMLKinetics.common import constants as cn
from SBMLKinetics.common.kinetic_law import KineticLaw
from SBMLKinetics.common.reaction import Reaction
from SBMLKinetics.common.function_definition import FunctionDefinition
from SBMLKinetics.common import util

import collections
import os.path
import numpy as np
import sys
import libsbml
import urllib3
import warnings
import zipfile


TYPE_MODEL = "type_model"  # libsbml model
TYPE_XML = "type_xml"  # XML string
TYPE_ANTIMONY = "type_xml"  # Antimony string
TYPE_FILE = "type_file" # File reference

# filename: name of file processed
# number: index of item
# model: libsbml.Model
IteratorItem = collections.namedtuple('IteratorItem',
    'filename number model')


class SimpleSBML(object):

  def __init__(self, model_reference):
    """
    Initializes instance variables
    :param str model_reference: string or SBML file or Roadrunner object
    """
    ##### PUBLIC #####
    self.model = None  # libsbml object
    self.reactions = []  # Python wrapper for Reaction
    self.species = []  # libsbml Species
    # Read the model
    if "RoadRunner" in str(type(model_reference)):
      model_reference = model_reference.getSBML()
    if util.isSBMLModel(model_reference):
      self.model = model_reference
    else:
      xml = util.getXML(model_reference)
      reader = libsbml.SBMLReader()
      document = reader.readSBMLFromString(xml)
      util.checkSBMLDocument(document, model_reference=model_reference)
      self.model = document.getModel()
    # Do the initializations
    self.species = [self.model.getSpecies(nn)
        for nn in range(self.model.getNumSpecies())]
    self.parameters = [self.model.getParameter(nn)
        for nn in range(self.model.getNumParameters())]
    self.function_definitions = self.getFunctionDefinitions()
    self.reactions = [Reaction(self.model.getReaction(nn), 
        function_definitions=self.function_definitions)
        for nn in range(self.model.getNumReactions())]

  def getFunctionDefinitions(self):
    """
    Returns all function definitions in the model.py

    Returns
    -------
    list-libsbml.FunctionDefinition
    """
    sbml_definitions =[self.model.getFunctionDefinition(n) for n 
        in range(self.model.getNumFunctionDefinitions())]
    function_definitions = [FunctionDefinition(s) for s in sbml_definitions]
    return function_definitions

  def getReaction(self, an_id):
    """
    Finds a reaction with the specified id.
    :param str an_id: id for the reaction
    :return Reaction/None:
    """
    return self._getInstance(self.reactions, an_id)

  def getSpecies(self, an_id):
    """
    Finds and returns the species with given name
    :param str an_id:
    Return None if there is no such species.
    """
    return self._getInstance(self.species, an_id)

  def getParameter(self, an_id):
    """
    Finds and returns the Parameter with given name
    :param str an_id:
    Return None if there is no such parameter.
    """
    return self._getInstance(self.parameters, an_id)

  def _getInstance(self, a_list, an_id):
    """
    Finds and returns the species with given name
    Return None if there is no such molecules
    :param str id:
    """
    results = [e for e in a_list if e.getId() == an_id]
    if len(results) > 1:
      raise ValueError(
          "Two instances have the same id: %s" %
          id)
    if len(results) == 0:
      return None
    return results[0]


#################### FUNCTIONS #########################
def readURL(url):
  """
  :param str url:
  :return str: file content
  """
  def do():
    http = urllib3.PoolManager()
    response = http.request('GET', url)
    return response.data.decode("utf-8") 
 # Catch bogus warnings
  with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    result = do()
  return result
  
def getZipfilePaths(data_dir=cn.BIOMODELS_DIR,
    zip_filename=cn.BIOMODELS_ZIP_FILENAME):
  """
  :param str data_dir: absolute path of the directory containing
      the xml files
  :param str zip_filename: name of the zipfile to process
  :return list-str, ZipFile: list of file paths, ZipFile object for file
  """
  path = os.path.join(data_dir, zip_filename)
  zipper = zipfile.ZipFile(path, "r")
  files = [f.filename for f in zipper.filelist]
  return files, zipper
  
def modelIterator(initial=0, final=1000,
    data_dir=cn.BIOMODELS_DIR,
    zip_filename=cn.BIOMODELS_ZIP_FILENAME):
  """
  Iterates across all models in a data directory.
  :param int initial: initial file to process
  :param int final: final file to process
  :param str data_dir: absolute path of the 
      directory containing
      the xml files
  :param str zip_filename: name of the zipfile to process. 
      If None, then looks for XML files in the directory.
  :return IteratorItem:
  """
  files, zipper = getZipfilePaths(
      data_dir=data_dir, zip_filename=zip_filename)
  # Functions for file types
  def readXML(filename):
    path = os.path.join(data_dir, filename)
    with open(path, 'r') as fd:
      lines = ''.join(fd.readlines())
    return lines
  def readZip(filename):
    with zipper.open(filename, 'r') as fid:
      lines = fid.read()
    return lines
  #
  if zip_filename is not None:
    read_func = readZip
  else:
    read_func = readXML
  begin_num = max(initial, 0)
  num = begin_num - 1
  end_num = min(len(files), final)
  for filename in files[begin_num:end_num]:
    num += 1
    lines = read_func(filename)
    if isinstance(lines, bytes):
      lines = lines.decode("utf-8") 
    try:
      #model could be invalid sbml
      model = SimpleSBML(lines)
      iterator_item = IteratorItem(filename=filename,
      model=model, number=num)
    except Exception as e:
      print(e)
      iterator_item = None
    yield iterator_item
