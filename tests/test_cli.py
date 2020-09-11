from iippl import cli

def test_cli_template():
    assert cli.cli(["header.snakefile"]) is None
