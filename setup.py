import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='ldsc2',
    version='0.0.32',
    author='Anthony Aylward',
    author_email='aaylward@eng.ucsd.edu',
    description='Simplify stratified LD score regression',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/anthony-aylward/ldsc2.git',
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent"
    ],
    install_requires=[
        'funcgenom',
        'gitpython',
        'pandas',
        'pybedtools',
        'numpy'
    ],
    entry_points={
        'console_scripts': [
            'ldsc2-download=ldsc2.download:main',
            'ldsc2-build-on-baseline=ldsc2.build_on_baseline:main',
            'ldsc2-strat=ldsc2.strat:main'
        ]
    }
)
