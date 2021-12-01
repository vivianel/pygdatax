import inspect


class FunctionDescription(object):

    def __init__(self, function_handle):
        self.args_name = []
        self.kwargs = {}
        self.kwargs_type = {}
        self.fullcommand = ''
        self.module_name = ''
        self.function_name = ''
        self.decorator = ''
        self.docstring = ''
        if inspect.isfunction(function_handle):
            self.module_name = inspect.getmodule(function_handle)
            self.function_name = function_handle.__name__
            source = inspect.getsource(function_handle)
            first_line = source.splitlines()[0]
            if '@' in first_line:
                self.decorator = first_line
            self.docstring = inspect.getdoc(function_handle)
            sign = inspect.signature(function_handle)
            params = sign.parameters
            for arg_name in params:
                if params[arg_name].default is inspect._empty:
                    self.args_name.append(arg_name)
                else:  # this is a keword argument thougz
                    default_value = params[arg_name].default
                    param_type = params[arg_name].annotation
                    if param_type is inspect._empty:
                        self.kwargs_type[arg_name] = None
                    else:
                        self.kwargs_type[arg_name] = param_type
                    if default_value is None:
                        self.kwargs[arg_name] = 'None'
                    else:
                        self.kwargs[arg_name] = str(default_value)
            self._makeDefaultCommandline()

    def _makeDefaultCommandline(self):
        self.fullcommand = self.function_name+'('
        for arg in self.args_name:
            self.fullcommand += arg+', '
        for key in self.kwargs:
            self.fullcommand += key+'='+self.kwargs[key]+', '
        # remove last comma
        self.fullcommand = self.fullcommand[:-2]
        self.fullcommand += ')'

    def makeCommandLine(self, **kwargs):
        commandLine = self.function_name+'('
        for arg in self.args_name:
            commandLine += arg+', '
        for key in kwargs:
            if key in self.kwargs:
                if self.kwargs_type[key] is str:
                    commandLine += key + '=' + '"' + kwargs[key] + '"' + ', '
                else:
                    commandLine += key + '=' + str(kwargs[key]) + ', '
        commandLine = commandLine[:-2]
        commandLine += ')'
        return commandLine



def get_commandList(module, decorator='@nxlib.treatment_function'):
    commandList = []
    members = inspect.getmembers(module)
    for m in members:
        des = FunctionDescription(m[1])
        if des.fullcommand != '' and decorator in des.decorator:
            commandList.append(des.fullcommand)
    return commandList


def get_functionList(module, decorator='@nxlib.treatment_function'):
    functionList = []
    members = inspect.getmembers(module)
    for m in members:
        des = FunctionDescription(m[1])
        if des.fullcommand != '' and  decorator in des.decorator:
            functionList.append(des.function_name)
    return functionList


def get_descriptionDict(module, decorator='@nxlib.treatment_function'):
    descriptionDict = {}
    members = inspect.getmembers(module)
    for m in members:
        des = FunctionDescription(m[1])
        if des.fullcommand != '' and decorator in des.decorator:
            descriptionDict[des.function_name] = des
    return descriptionDict


if __name__ == '__main__':
    from pygdatax.instruments import sansllb

    azi = FunctionDescription(sansllb.azimutal_integration)
    print(azi.makeCommandLine(w='2'))
    l = get_commandList(sansllb)
    print(l)
