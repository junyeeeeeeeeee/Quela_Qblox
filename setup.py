from setuptools import setup, find_packages

setup(
    name='Quela_Qblox',          # Replace with your package name
    version='2024.11.13',
    author='Dai-Jia, Wu',
    author_email='porkface0301@gmail.com',
    description='AS measurement software for Qblox',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/junyeeeeeeeeee/Quela_Qblox',
    packages=find_packages(),          # Automatically find package directories
    install_requires=[                 # Optional, if you have dependencies
        'qblox-instruments==0.12.0',
        'quantify-core==0.7.4',
        'quantify-scheduler==0.20.0',
        'colorama',
        'scikit-learn',
        'scipy',
        'sympy',
        'numpy>=1.26.2',
        'matplotlib>=3.9',
        'xarray>=2024.6.0'
    ],
    classifiers=[
        'Programming Language :: Python :: 3.10',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.10',
)
