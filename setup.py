from setuptools import setup, find_packages

setup(
    name="OMKar",
    version="1.0.0",
    author="Siavash Raeisi Dehkordi, Zhaoyang Jia",
    author_email="siavash.raisi@gmail.com",
    description="A computational tool for automated karyotyping using OGM data.",
    url="https://github.com/siavashre/OMKar/",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "pandas",
        "pulp",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
)
