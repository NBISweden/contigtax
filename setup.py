import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="tango",
    version="0.0.1",
    author="John Sundh",
    author_email="john.sundh@scilifelab.se",
    description="A package to add taxonomy to metagenomic contigs",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/johnne/tango",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU GENERAL PUBLIC LICENSE",
        "Operating System :: OS Independent",
    ],
)
