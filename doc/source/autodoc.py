import inspect 

class AutodocCollector:
    entries = []
    def __init__(self, module, dunder=False):
        for key in dir(module):
            thing = getattr(module, key)
            if key.startswith('__') and not dunder:
                continue
            
            if 'nanobind.nb_func' in str(type(thing)):
                self.parse_func(thing)
            elif 'enum.EnumType' in str(type(thing)):
                self.parse_enum(thing)
            elif 'nanobind.nb_type_0' in str(type(thing)):
                self.parse_class(thing)
            else:
                print(key, type(thing), thing.__name__, inspect.getmembers(thing))
            
    def add_entry(self, e):
        self.entries.append(e)
        
    def to_file(self, path, title):
        with open(path, 'w') as fp:
            fp.write(title+'\n')
            fp.write('='*len(title)+'\n\n')
            fp.write('\n\n'.join(self.entries))
                            
    def parse_func(self, thing):
        mother = inspect.getmodule(thing)
        self.add_entry(f'.. autofunction:: {mother.__name__}.{thing.__name__}\n\n    :teqpflsh:`{thing.__name__}`')
        
    def parse_enum(self, thing):
        mother = inspect.getmodule(thing)
        self.add_entry(f'.. autoclass:: {mother.__name__}.{thing.__name__}\n    :undoc-members:\n    :members:\n    :show-inheritance:\n\n    :teqpflsh:`{thing.__name__}`')
        
    def parse_class(self, thing):
        mother = inspect.getmodule(thing)
        members = inspect.getmembers(thing)
        e = f'.. autoclass:: {mother.__name__}.{thing.__name__}\n    :show-inheritance:'
        e += f'\n\n    C++ docs: :teqpflsh:`{thing.__name__}`'
        for name, member in members:
            if name.startswith('__'): continue
            if isinstance(member, property):
                e += '\n\n' + ' '*4 + f'.. autoproperty:: {name}\n'
            else:
                e += '\n\n' + ' '*4 + f'.. automethod:: {name}\n'
        self.add_entry(e)
            
if __name__ == '__main__':
    import teqpflsh
    ac = AutodocCollector(teqpflsh._teqpflsh_impl)
    ac.to_file('api/teqpflsh.rst', title='teqpflsh Package')