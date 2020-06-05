import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="bioCanon-phac", # Replace with your own username
    version="0.0.8",
    author="James Robertson, Amanda Saunders, Justin Schonfeld",
    author_email="justin.schonfeld@canada.ca",
    description="A package for generating potential biohansel schemes",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/phac-nml/bioCanon",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
