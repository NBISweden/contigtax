import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="contigtax",
    version="0.5.9",
    author="John Sundh",
    author_email="john.sundh@scilifelab.se",
    description="A package to assign taxonomy to metagenomic contigs",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/NBISweden/contigtax",
    packages=["contigtax"],
    entry_points={'console_scripts': ['contigtax = contigtax.__main__:main']},
    scripts=["contigtax/evaluate_contigtax.py"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
