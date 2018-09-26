from setuptools import setup
import os

def readme():
    with open(os.path.dirname(__file__) + '/README.rst') as rst:
        return rst.read()

setup(
    name='spectre',
    version='0.1',
    description='Python package for the spectral coherence analysis of ribosome profiling sequence data.',
    long_description=readme(),
    url='https://github.com/mills-lab/spectre/',
    author='Sang Y. Chun',
    author_email='stonyc@umich.edu',
    license='BSD 3-Clause',
    install_requires=['pytest','numpy','pandas','docopts','HTSeq'],
    packages=['spectre', 'tests'],
    package_dir={'spectre': './spectre', 'tests': './tests'},
    package_data={'spectre': ['data/*', 'LICENSE', 'CONTRIBUTING.md', 'CODE_OF_CONDUCT.md']},
    classifiers=[
        'Development Status :: 1 - Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: End Users/Desktop',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Information Technology',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        'Operating System :: Unix',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: MacOS',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.0',
        'Programming Language :: Python :: 3.1',
        'Programming Language :: Python :: 3.2',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Information Analysis'
    ],
    keywords='ribosome profiling coherence translation',
    include_package_data=True,
    zip_safe=False
)