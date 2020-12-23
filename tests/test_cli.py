from snakeobjects import cli

def test_cli_jobscript():
    assert cli.cli(["jobscript.sh"]) is None

