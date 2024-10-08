from setuptools import find_packages, setup

with open("README.md", "r") as f:
    long_description = f.read()

setup(
    name="ASQuelaQblox",
    version="2024.08.12",
    description="Academia Sinica Quela-Qblox toolbox",
    package_dir={"": "src",},
    # packages=find_packages(where="src"),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/junyeeeeeeeeee/Quela_Qblox",
    author= ["Dai-Jia, Wu", "Jun-Yee, Lee", "Wei-En, Lin"],
    author_email="porkface0301@gmail.com",
    license="MIT",
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.10",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        "qblox-instruments==0.12.0",
        "quantify-core==0.7.4",
        "quantify-scheduler==0.20.0",
        "colorama",
        "scikit-learn",
        "numpy==1.26.2",
        "scipy",
        "matplotlib>=3.9",
        "xarray==2023.12.0"
    ],
    # extras_require={
    #     "dev": ["pytest>=7.0", "twine>=4.0.2"],
    # },
    python_requires=">=3.10",
)