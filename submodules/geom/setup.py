from setuptools import setup, find_packages

# Read README.md for PyPI long description
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="geom",
    version="1.0.0",
    author="Pablo Grobas Illobre",
    author_email="pgrobasillobre@gmail.com", 
    description="GEOM: A CLI for geometry manipulation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pgrobasillobre/geom", 
    packages=find_packages(include=["geom", "geom.*"]),
    install_requires=[
        "gmsh>=4.11.1",
        "ase>=3.22.1",
        "numpy>=1.24.3",
        "pytest>=8.3.4",
        "launchpadlib>=2.1.0"
    ],
    entry_points={
        "console_scripts": [
            "geom = geom.__main__:main",
        ]
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Development Status :: 4 - Beta",  # Change to '5 - Production/Stable' when ready
        "Intended Audience :: Developers",
        "Topic :: Scientific/Engineering",
    ],
    python_requires=">=3.6",
    include_package_data=True,
    project_urls={
        "Bug Tracker": "https://github.com/pgrobasillobre/issues",  # Update GitHub URL
        "Documentation": "https://github.com/pgrobasillobre/geom/wiki",
    },
)

