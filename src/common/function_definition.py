"""Representation of an SBML function defintion."""

import libsbml


class FunctionDefinition():

  def __init__(self, sbml_function_definition):
    self.sbml_function_definition = sbml_function_definition
    self.function_name = self.sbml_function_definition.getName()
    self.id = self.sbml_function_definition.getId()
    self.argument_names = [self.sbml_function_definition.getArgument(n).getName()
        for n in range(self.sbml_function_definition.getNumArguments())]
    self.body = libsbml.formulaToL3String(self.sbml_function_definition.getBody())

  def __repr__(self):
    argument_call = ",".join(self.argument_names)
    call_str = "%s(%s)" % (self.id, argument_call)
    return "%s: %s" % (call_str, self.body)
