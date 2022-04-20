"""Representation of an SBML function defintion."""

import libsbml


class FunctionDefinition():

    def __init__(self, sbml_function_definition,
          name=None, fid=None, arguments=None, body=None):
        """
        Parameters
        ----------
        sbml_function_definition: libsbml.FunctionDefinition
        name: str
        fid: str
        arguments: list-str
        body: str
        """
        self.sbml_function_definition = sbml_function_definition
        if self.sbml_function_definition is not None:
            name = self.sbml_function_definition.getName()
            fid = self.sbml_function_definition.getId()
            arguments = [
                  self.sbml_function_definition.getArgument(n).getName()
                  for n in 
                  range(self.sbml_function_definition.getNumArguments())]
            body = libsbml.formulaToL3String(
                  self.sbml_function_definition.getBody())
        self.function_name = name
        self.id = fid
        self.argument_names = arguments
        self.body = body
  
    def __repr__(self):
        argument_call = ",".join(self.argument_names)
        call_str = "%s(%s)" % (self.id, argument_call)
        return "%s: %s" % (call_str, self.body)
  
    @classmethod
    def makeBuiltinFunctions(cls):
        """
        Creates a function definition for the SBML delay function.

        Returns
        -------
        list-FunctionDefinition
        """
        functions = []
        # Delay
        functions.append(cls(None, name="delay", fid="delay",
               arguments=["a_species", "num"], body = "a_species"))
        # Exponential
        functions.append(cls(None, name="exp", fid="exp",
               arguments=["num"], body = "2.71828182**num"))
        #
        return functions
