from setuptools import setup
import versioneer

requirements = [
    # package requirements go here
]

setup(
    name='snakeobjects',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="Snakeobjects, an object-oriented workflow management system based on snakemake.",
    license="MIT",
    author="Ivan Iossifov",
    author_email='iossifov@cshl.edu',
    url='https://github.com/iossifovlab/snakeobjects',
    packages=['snakeobjects'],
    package_data={
        "snakeobjects": ["jobscript.sh","header.snakefile"],
    },
    include_package_data=True,
    scripts = [],
    entry_points={
        'console_scripts': [
            'sobjects=snakeobjects.cli:cli'
        ]
    },
    install_requires=requirements,
    keywords=['snakeobjects','snakemake','workflow managment'],
    classifiers=[
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7'
    ]
)
