from setuptools import setup, find_namespace_packages

setup(
    name='bifrost_sp_ecoli',
    version='0.0.1',
    description='Ecoli serotyping and stx typing for Bifrost',
    url='https://github.com/ssi-dk/bifrost_sp_ecoli',
    author='Rasmus',
    author_email='raah@ssi.dk',
    license='MIT',
    packages=find_namespace_packages(),
    install_requires=[
        'bifrostlib >= 2.1.9',
        'biopython>=1.77'
    ],
    python_requires='>=3.11',
    package_data={'bifrost_sp_ecoli': ['config.yaml', 'pipeline.smk']},
    include_package_data=True
)
