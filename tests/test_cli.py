from snakeobjects import cli

def test_cli_jobscript():
    assert cli.cli(["jobscript.sh"]) is None

def test_cli_header():
    assert cli.cli(["header.snakefile"]) is None
 
