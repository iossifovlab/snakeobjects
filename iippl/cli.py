from argparse import ArgumentParser
from iippl import __version__
import importlib.resources as importlib_resources


def cli(args=None):
    p = ArgumentParser(
        description="Our pipeline infrastructure.",
        conflict_handler='resolve'
    )
    p.add_argument(
        '-V', '--version',
        action='version',
        help='Show the conda-prefix-replacement version number and exit.',
        version="iippl %s" % __version__,
    )
    p.add_argument('command')

    args = p.parse_args(args)

    if args.command == "main.snakefile":
        print(importlib_resources.read_text(__package__,'main.snakefile'),end='')
    else:
        print("Don't know the command:", args.command)
        return 1

    # No return value means no error.
    # Return a value of 1 or higher to signify an error.
    # See https://docs.python.org/3/library/sys.html#sys.exit



if __name__ == '__main__':
    import sys
    cli(sys.argv[1:])
