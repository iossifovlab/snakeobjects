from snakeobjects import cli

def test_cli_jobscript():
    assert cli.cli(["version"]) is None
    assert cli.cli(["help"]) is None

