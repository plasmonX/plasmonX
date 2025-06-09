
class parameters:
   """
   Parameters Class for Metal and Graphene Structures.

   This class provides essential lattice parameters, minimum interatomic distances, 
   and atomic arrangements for various elements, including FCC, BCC, HCP metals, 
   and graphene.

   Attributes:
       convergence (float): Convergence threshold for numerical operations.
       convergence_step (float): Step size for convergence calculations.
       min_dist_translate (float): Minimum allowable distance for translation operations.

       lattice_constant (dict): Lattice constants (in Ångströms) for various elements.
           - Keys: Atomic symbols (e.g., "Ag", "Au", "C").
           - Values: Single float (for cubic structures) or tuple (for HCP structures).

       min_dist (dict): Minimum interatomic distances (in Ångströms).
           - Keys: Atomic symbols.
           - Values: Distance thresholds based on lattice parameters.

       atomic_arrangement (dict): Crystallographic structure of elements.
           - Keys: Atomic symbols.
           - Values: String representation of atomic packing ("FCC", "BCC", "HCP").

       metal_atomtypes (list): List of valid atomic types for metal nanoparticles.
   
   Sources:
       - Chemistry LibreTexts - Crystal Structures of Metals
         https://chem.libretexts.org/Bookshelves/Inorganic_Chemistry/Introduction_to_Inorganic_Chemistry_(Wikibook)/06%3A_Metals_and_Alloys-_Structure_Bonding_Electronic_and_Magnetic_Properties/6.04%3A_Crystal_Structures_of_Metals
       - Wikipedia - Graphene
         https://en.wikipedia.org/wiki/Graphene

         (Both accessed on 21st February 2025.)
   
   Notes:
       - All values are in Ångströms.
       - The `lattice_constant` dictionary includes FCC, BCC, HCP metals and graphene.
       - The `min_dist` dictionary ensures reasonable atomic separations for simulations.
       - The class is mainly used for molecular structure generation and validation.
   """

   def __init__(self):
      """
      Initializes default parameters for molecular structures.

      Attributes Initialized:
          - `convergence` and `convergence_step`: Thresholds for numerical calculations.
          - `lattice_constant`: Dictionary of atomic lattice parameters (FCC, BCC, HCP, Graphene).
          - `min_dist`: Minimum interatomic distances to avoid structural inconsistencies.
          - `atomic_arrangement`: Defines the atomic packing type for different elements.
          - `metal_atomtypes`: List of metal elements supported in the simulations.
      """
       
      # Convergence options
      self.convergence      = 0.010
      self.convergence_step = 0.001

      # General parameters (units = Angstroms)
      self.min_dist_translate = 1.0

      # Lattice parameters of metals / graphene (units = Angstroms)
      self.lattice_constant = {
          "ag": 4.09000, "au": 4.09000, "al": 4.05000, "ni": 3.52000, "cu": 3.61000,
          "pd": 3.89000, "pt": 3.92000, "pb": 4.95000, "th": 5.08000, "ce": 3.68000, "yb": 5.49000,
          "na": 4.22500, "fe": 2.87000, "mo": 3.15000, "w": 3.16000, "v": 3.02000,
          "nb": 3.30000, "ta": 3.31000, "eu": 4.61000, "ba": 5.02000, "ra": 5.15000,
          "cr": 2.88000, "li": 3.49000, "k": 5.23000, "rb": 5.59000, "cs": 6.05000, "fr": 5.70000,
          "mg": (3.21000, 5.21000), "ti": (2.95000, 4.68000), "zn": (2.66000, 4.95000),
          "os": (2.73000, 4.32000), "sc": (3.31000, 5.27000), "zr": (3.23000, 5.15000),
          "hf": (3.20000, 5.08000), "gd": (3.63000, 5.78000), "tb": (3.62000, 5.69000),
          "dy": (3.61000, 5.65000), "ho": (3.60000, 5.61000), "er": (3.59000, 5.60000),
          "tm": (3.58000, 5.57000), "lu": 3.51000, "be": (2.29000, 3.58000),
          "re": (2.76000, 4.46000), "co": (3.54000, 4.52000), "ru": (2.71000, 4.28000),
          "c": 2.46  # Graphene
      }

      # Merge cutoff as minimum interatomic distance given the lattice constants above (units = Angstroms)
      self.min_dist = {
          "ag": 2.88, "au": 2.88, "al": 2.86378, "ni": 2.48902, "cu": 2.55266,
          "pd": 2.75065, "pt": 2.77186, "pb": 3.50018, "th": 3.59210, "ce": 2.60215, "yb": 3.88202,
          "na": 3.65, "fe": 2.48549, "mo": 2.72798, "w": 2.73664, "v": 2.61540,
          "nb": 2.85788, "ta": 2.86654, "eu": 3.99238, "ba": 4.34745, "ra": 4.46003,
          "cr": 2.49415, "li": 3.02243, "k": 4.52931, "rb": 4.84108, "cs": 5.23945, "fr": 4.93634,
          "mg": 3.21, "ti": 2.95, "zn": 2.66, "os": 2.73, "sc": 3.31, "zr": 3.23,
          "hf": 3.2, "gd": 3.63, "tb": 3.62, "dy": 3.61, "ho": 3.6, "er": 3.59,
          "tm": 3.58, "lu": 3.51, "be": 2.29, "re": 2.76, "co": 3.54, "ru": 2.71,
          "c": 1.42 
      }

      # Atomic packing of metals
      self.atomic_arrangement = {
          "ag": "FCC", "au": "FCC", "al": "FCC", "ni": "FCC", "cu": "FCC",
          "pd": "FCC", "pt": "FCC", "pb": "FCC", "th": "FCC", "ce": "FCC", "yb": "FCC",
          "na": "BCC", "fe": "BCC", "mo": "BCC", "w": "BCC", "v": "BCC",
          "nb": "BCC", "ta": "BCC", "eu": "BCC", "ba": "BCC", "ra": "BCC",
          "cr": "BCC", "li": "BCC", "k": "BCC", "rb": "BCC", "cs": "BCC", "fr": "BCC",
          "mg": "HCP", "ti": "HCP", "zn": "HCP", "os": "HCP", "sc": "HCP", "zr": "HCP",
          "hf": "HCP", "gd": "HCP", "tb": "HCP", "dy": "HCP", "ho": "HCP", "er": "HCP",
          "tm": "HCP", "lu": "HCP", "be": "HCP", "re": "HCP", "co": "HCP", "ru": "HCP"
      }

      # List of accepted atom types
      self.metal_atomtypes = list(self.lattice_constant.keys())
