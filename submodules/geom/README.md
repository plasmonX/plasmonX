# GEOM – Geometry File Management

**GEOM** is a command-line tool for creating, modifying, and analyzing XYZ geometry files. It provides **geometry transformations, nanoparticle generation, and graphene structure creation** for computational research.

## Features

- **Geometry Transformations**: Translation, rotation, merging, and specular (mirror) transformations.
- **Nanoparticle Generation**: Sphere, rod, core-shell, tip, pyramid, cone, icosahedron, and more.
- **Graphene Structures**: Ribbons, disks, rings, and triangles.
- **Advanced Options**: Alloying, dimer formation, and bowtie configurations.
- **Minimum Distance Calculation** between XYZ geometries.
- **Geometrical Center Computation**.

## Installation

GEOM requires **Python 3** and the following dependencies:

- `gmsh>=4.11.1`
- `ase>=3.22.1`
- `numpy>=1.24.3`
- `pytest>=8.3.4`
- `launchpadlib>=2.1.0`

### Install with:
```bash
./install.sh
```

## Usage

After installation, load the GEOM environment by running the `geom_load` function, which will configure the necessary environment variables and aliases:

```
geom_load
```

Once the environment is set up, run the following command to see the available options:

```
geom -h
```

This will display the help menu with all the available commands and their descriptions.

Example commands:

- **Rotate geometry 90 degrees** around the Y-axis:

```
geom -r1 90 geom.xyz origin_CM_yes +y
```

- Generate a nanoparticle sphere:

```
geom -create -sphere Ag 30
```

- Generate a graphene ribbon: 

```
geom -create -graphene rib 50 20
```

## Running Tests

After running ./install.sh, the tests are executed automatically.

To manually run the tests again:

```
./geom/tests/run_all_tests.sh
```

## License

GEOM is licensed under the **GNU General Public License v3.0**.

## Contact

For issues or contributions:

- Email: **pgrobasillobre@gmail.com**
- Github issues: https://github.com/pgrobasillobre/geom/issues
