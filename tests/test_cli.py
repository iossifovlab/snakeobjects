from iippl import cli

def test_cli_template():
    assert cli.cli(["main.snakefile"]) is None
