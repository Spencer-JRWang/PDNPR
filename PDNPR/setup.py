from setuptools import setup, find_packages
setup(
    name='PDNPR',
    version='1.0',
    description='protein dynamic network pathway runner',
    author='Wang Jingran',
    author_email='jrwangspencer@stu.suda.edu.cn',
    packages=find_packages(),
    include_package_data=True,
    url='https://github.com/Spencer-JRWang/PDNPR',
    install_requires=[
        "pymol-open-source",
        "networkx",
        "matplotlib",
        "numpy",
        "mdtraj",
        "tkinker"
    ]
)