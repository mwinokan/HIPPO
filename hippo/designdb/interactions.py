"""Constants for protein-ligand interaction detection (fingerprinting).

The feature families and complementary-feature / interaction-type maps are
re-exported from ``molparse`` (the single source of truth). The distance/angle
cutoffs live here, ported from the legacy ``hippo`` ``pose`` module.

Used by :class:`.InteractionService` and the :class:`.Pose` component.
"""

from molparse.rdkit.features import COMPLEMENTARY_FEATURES, FEATURE_FAMILIES, INTERACTION_TYPES

# maximum centroid-centroid distance (Angstrom) for each interaction type
INTERACTION_CUTOFF = {
    'Hydrophobic': 4.5,
    'Hydrogen Bond': 3.5,
    'Electrostatic': 4.5,
    'π-stacking': 6.0,
    'π-cation': 4.5,
    # https://pubs.acs.org/doi/full/10.1021/acs.cgd.5b01058
    'Sulfur-Sulfur': 4.0,
}

# π-stacking geometry cutoffs (Angstrom)
PI_STACK_MIN_CUTOFF = 3.8
PI_STACK_F2F_CUTOFF = 4.5
PI_STACK_E2F_CUTOFF = 6.0

# warn (rather than silently skip) when a residue mismatch (mutation) occurs
# within this distance (Angstrom) of any ligand feature
MUTATION_WARNING_DIST = 15

__all__ = [
    'FEATURE_FAMILIES',
    'COMPLEMENTARY_FEATURES',
    'INTERACTION_TYPES',
    'INTERACTION_CUTOFF',
    'PI_STACK_MIN_CUTOFF',
    'PI_STACK_F2F_CUTOFF',
    'PI_STACK_E2F_CUTOFF',
    'MUTATION_WARNING_DIST',
]
