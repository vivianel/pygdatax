import inspect


class FunctionDescription(object):
    args_name = []
    kwargs = {}
    fullcommand = ''
    module_name = ''
    function_name=''

    def __init__(self, function_handle):
        if inspect.isfunction(function_handle):
            self.module_name = inspect.getmodule(function_handle)
            self.function_name = function_handle.__name__
            source = inspect.getsource(function_handle)
            first_line = source.splitlines()
            if '@nxlib.treatment_function' in first_line:
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




if __name__ =='__main__':
    import inspect
    from pygdatax import xeuss
    azi = FunctionDescription(xeuss.azimutal_integration)
    print(azi.fullcommand)
