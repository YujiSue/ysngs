from setuptools import setup, find_packages
setup(
    name='ysngs',
    version='0.2.1',
    author='Yuji Suehiro',
    packages=find_packages(),
    package_data={'': ['ngsapp.json']},
    url='https://github.com/YujiSue/ysngs',
    description='Scripts for NGS data analysis.'
)
