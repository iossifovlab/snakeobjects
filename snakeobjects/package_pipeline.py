import pkg_resources
from importlib_metadata import entry_points

print("Boris' approach")
pkg_name, data_name = "example_package", "snakefiles"
snakefile_directory = pkg_resources.resource_filename(pkg_name, data_name)
print(f"\tpkg_name: {pkg_name}")
print(f"\tdata_name: {data_name}")
print(f"\tsnakefile_directory: {snakefile_directory}")
print()

print("Ivan's approach")

discovered_plugins = entry_points(group='snakeobjects.pipeline')
for dp in discovered_plugins:
    print(dp.name, dp, dp.load(), sep="\n")
    print()
