from setuptools import setup

setup(name='phylotyper',
	version='0.1.0',
	packages=['phylotyper'],
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
		'Bio',
		'pkg_resources',
		'yaml',
		'rpy2',
	],
	include_package_data = True,
	zip_safe = False
)
