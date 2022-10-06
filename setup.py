from setuptools import setup, find_packages

setup(
	name='Mulguisin',
	version='0.4',
	license='BSD',
	description='Mulguisin: Cluster/group finding algorithm',
	author='Young Ju, Inkyu Park, Sungwook E. Hong, Cristiano G. Sabiu',
	author_email='youngju@uos.ac.kr',
	packages=find_packages(),
	#package_dir={'': 'Mulguisin'},
	url='https://github.com/youngju20/Mulguisin',
	install_requires=[
		'scipy','numpy','h5py','numba',
	],
		
)
