from setuptools import setup
import versioneer

requirements = [
    # package requirements go here
]

setup(
    name='pipes',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="Our pipeline infrastructure.",
    license="MIT",
    author="Ivan Iossifov",
    author_email='iossifov@cshl.edu',
    url='https://github.com/iossifovlab/pipes',
    packages=['iippl'],
    entry_points={
        'console_scripts': [
            'iippl=iippl.cli:cli'
        ]
    },
    install_requires=requirements,
    keywords='pipes','iippl'
    classifiers=[
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ]
)
