from setuptools import setup, find_packages

setup(name='phylotyper',
	version='0.1.0',
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
		'pyaml',
		'rpy2',
	],
	include_package_data = True,
	zip_safe = False,
	packages = find_packages(exclude=["*.tests", "*.tests.*", "tests.*", "tests", "*.data", "*.data.*", "data.*", "data"]),
	test_suite = 'nose.collector',
	tests_require = ['nose']
)
