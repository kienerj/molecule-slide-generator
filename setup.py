import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="molecule-slide_generator",
    version="0.0.1",    
    author="Joos Kiener",
    author_email="joos.kiener@gmail.com",
    description="Generate images of molecules and their properties for use in presentations and reports ",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/kienerj/molecule-slide-generator",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Development Status :: 4 - Beta"
    ],
    python_requires='>=3.6',
    install_requires=[
       'rdkit>=2020.09.1',
       'numpy',
       'Pillow'
    ]
)