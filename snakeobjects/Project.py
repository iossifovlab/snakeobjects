
class Project:
    """This class represents a snakeobject project. 
       A snakeobject project attributes are:
       - directory, the project directory 
       - parameters, a key value dictionalry for global project lelvel parameters.
       - objectGraph, the objectGraph for the project
    """

    def __init__(self,directory=None):
        """The constructor of a snakemake project object. 
        The project directory is given as a parameter.  
        If none, the values of the SO_PROJECT is used as a project directory.
        If non existent, the current working directory and its parent directories are examined in order if they contains a so_project.yaml file. The first directory found to contain so_project.yaml file is used as a project directory.
        If no directory is found to contain a so_project.yaml files, the current workgin directory is used as a project.
        """ 

        pass


