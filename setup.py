from setuptools import setup, find_packages

setup(name='phylotyper',
    version='0.1.2',
    description='Phylogenetic-based prediction of subtypes',
    url='https://github.com/superphy/insilico-subtyping',
    author='Matt Whiteside',
    author_email='matthew.whiteside@phac-aspc.gc.ca',
    license='APL 2.0',
    entry_points={
        'console_scripts': [
            'phylotyper = phylotyper.__main__:main'
        ]
    },
    install_requires=[
        'biopython',
        'pyyaml',
        'rpy2<=2.8',
        'setuptools',
    ],
    install_package_data = True,
    zip_safe = False,
    packages = find_packages(),
    package_data = { 'phylotyper': ['data/*/*', 'R/*', 'subtypes_index.yaml']},
    # test_suite = 'nose.collector',
    # tests_require = ['nose']
    # setup_requires=['pytest-runner'],
    # tests_require = ['pytest']
)
