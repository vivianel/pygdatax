import inspect


class FunctionDescription(object):

    def __init__(self, function_handle):
        self.args_name = []
        self.kwargs = {}
        self.fullcommand = ''
        self.module_name = ''
        self.function_name = ''
        self.decorator = None
        if inspect.isfunction(function_handle):
            self.module_name = inspect.getmodule(function_handle)
            self.function_name = function_handle.__name__
            source = inspect.getsource(function_handle)
            first_line = source.splitlines()[0]
            if '@' in first_line:
                self.decorator = first_line
            sign = inspect.signature(function_handle)
            params = sign.parameters
            for arg_name in params:
                if params[arg_name].default is inspect._empty:
                    self.args_name.append(arg_name)
                else: # this is a keword argument thougz
                    default_value = params[arg_name].default
                    if default_value is None:
                        self.kwargs[arg_name] = 'None'
                    else:
                        self.kwargs[arg_name] = str(default_value)
            self._makeCommandline()

    def _makeCommandline(self):
        self.fullcommand = self.function_name+'('
        for arg in self.args_name:
            self.fullcommand += arg+', '
        for key in self.kwargs:
            self.fullcommand += key+'='+self.kwargs[key]+', '
        # remove last comma
        self.fullcommand = self.fullcommand[:-2]
        self.fullcommand += ')'


def get_commandList(module, decorator='@nxlib.treatment_function'):
    commandList = []
    members = inspect.getmembers(module)
    for m in members:
        des = FunctionDescription(m[1])
        if des.fullcommand != '' and des.decorator == decorator:
            commandList.append(des.fullcommand)
    return commandList

if __name__ =='__main__':
    from pygdatax import xeuss
    azi = FunctionDescription(xeuss.azimutal_integration)
    l = get_commandList(xeuss)
    print(l)
