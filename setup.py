import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="tango",
    version="0.4.1",
    author="John Sundh",
    author_email="john.sundh@scilifelab.se",
    description="A package to assign taxonomy to metagenomic contigs",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/johnne/tango",
    packages=["tango"],
    entry_points={'console_scripts': ['tango = tango.__main__:main']},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU GENERAL PUBLIC LICENSE",
        "Operating System :: OS Independent",
    ],
)
