
def load_yaml_with_envirnomemnt_interpolation(fName):

    CF = open(fname, 'r')
    config = yaml.safe_load(CF)
    CF.close()

    ptn = re.compile(r'(^\[E\:)(.*)(\])')

    for k,v in config.items():
        if type(v) != str:
            continue
        m = ptn.match(v)
        if m:
            s = m.span()
            name = m.groups(0)[1]
            n = os.environ[name]
            if n:
                config[k] = v.replace(v[s[0]:s[1]],n)
            else:
                print('Varianble %s is not defined' % name, file=sys.stderr)
                exit(1) 
    return config  


class Project:
    """
        This class represents a snakeobject project. 
        
        The project directory is given as the ``directory`` parameter.  
        If none, the values of the SO_PROJECT is used as a project directory.
        If non existent, the current working directory and its parent directories are examined in order if they contains a so_project.yaml file. The first directory found to contain so_project.yaml file is used as a project directory.
        If no directory is found to contain a so_project.yaml files, the current workgin directory is used as a project.

        A snakeobject project attributes are:

        * ``directory``, the project directory 
        * ``parameters``, a key value dictionalry for global project lelvel parameters.
        * ``objectGraph``, the objectGraph for the project
    """

    def __init__(self,directory=None):

        pass

