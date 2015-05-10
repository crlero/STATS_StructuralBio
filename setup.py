from distutils.core import setup


setup(
      name='statools',
      
      version='1.0',
      
      description='A sample Structural Biology-Python project installable for the Calculation\
      of Protein structure Statistical Potentials.',
      author='Jorge Roel & Cristina Leal',
      author_email='cristina.leal01@estudiant.upf.edu',
      
      classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Scientific-comunity-Students',
        'Programming Language :: Python :: 3.4',
    ],
      
			packages=['statools'],
      scripts =['statools/stats'],)
