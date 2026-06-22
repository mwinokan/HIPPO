"""Protein-ligand interaction fingerprinting.

Detects interactions for a :class:`.Pose`: extracts protein features (populating
:class:`.FeatureModel`), runs the geometric detector, resolves duplicate
interactions, and populates :class:`.InteractionModel`. Protein features come from
the pose's own ``protein_system``. Entry point: :meth:`.Pose.calculate_interactions`.
"""

import json

import mrich
import numpy as np
from designdb.interactions import (
    COMPLEMENTARY_FEATURES,
    INTERACTION_CUTOFF,
    INTERACTION_TYPES,
    PI_STACK_F2F_CUTOFF,
    PI_STACK_MIN_CUTOFF,
)
from designdb.models import FeatureModel, InteractionModel


def _norm(coords) -> np.ndarray:
    """Principal axis (first eigenvector of the covariance) of a set of points."""
    coords = np.array(coords).T
    cov = np.cov(coords)
    eig = np.linalg.eig(cov)
    return eig[1][:, 0]


def _unit_vector(vector) -> np.ndarray:
    """Unit vector in the direction of ``vector``."""
    return vector / np.linalg.norm(vector)


def _angle_between(v1, v2) -> float:
    """Angle (degrees, folded to 0-90) between two vectors."""
    v1_u = _unit_vector(v1)
    v2_u = _unit_vector(v2)
    a = 180 * np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)) / np.pi
    if a > 90:
        a = 180 - a
    return a


class InteractionService:
    """Construction and persistence of pose-protein interactions."""

    @staticmethod
    def calculate(
        pose,
        *,
        resolve: bool = True,
        distance_padding: float = 0.0,
        angle_padding: float = 0.0,
        force: bool = False,
        commit: bool = True,
        debug: bool = False,
    ) -> None:
        """Enumerate valid interactions between a pose's ligand and protein.

        :param pose: the :class:`.Pose` to fingerprint
        :param resolve: cull duplicate / less-significant interactions
        :param distance_padding: padding (Angstrom) added to all distance cutoffs
        :param angle_padding: padding (degrees) added to all angle cutoffs
        :param force: recalculate even if the pose is already fingerprinted
        :param commit: persist the results to the database
        :param debug: increase verbosity
        """

        if pose.has_fingerprint and not force:
            if debug:
                mrich.warning(f'{pose} is already fingerprinted')
            return

        protein_system = pose.protein_system
        if protein_system is None and pose.reference is not None:
            protein_system = pose.reference.protein_system

        if protein_system is None:
            raise NotImplementedError(
                f'No protein system for {pose} '
                f'(protein_link={pose.protein_link!r}, reference={pose.reference_id})'
            )

        mol = pose.mol
        if mol is None:
            mrich.error(f'Could not read molecule for {pose}')
            return

        candidates = InteractionService._detect(
            pose=pose,
            protein_system=protein_system,
            mol=mol,
            distance_padding=distance_padding,
            angle_padding=angle_padding,
            debug=debug,
        )

        if resolve:
            candidates = InteractionService._resolve(candidates, debug=debug)

        if not commit:
            return

        # replace any existing interactions for this pose
        InteractionModel.objects.filter(pose_id=pose.id).delete()

        InteractionModel.objects.bulk_create(
            [
                InteractionModel(
                    feature_id=c['feature_id'],
                    pose_id=pose.id,
                    interaction_type=c['type'],
                    interaction_family=c['family'],
                    interaction_atom_id=json.dumps(c['atom_ids']),
                    interaction_prot_coord=json.dumps(c['prot_coord']),
                    interaction_lig_coord=json.dumps(c['lig_coord']),
                    interaction_distance=c['distance'],
                    interaction_angle=c['angle'],
                    interaction_energy=None,
                )
                for c in candidates
            ]
        )

        pose.set_has_fingerprint(True, commit=commit)

        if debug:
            mrich.success(f'{pose}: {len(candidates)} interactions')

    @staticmethod
    def _protein_feature_id(target, prot_feature) -> int:
        """Get-or-create the :class:`.FeatureModel` for a molparse protein feature."""
        atom_name = ' '.join(a.name for a in prot_feature.atoms)
        feature, _ = FeatureModel.objects.get_or_create(
            feature_family=prot_feature.family,
            target=target,
            feature_chain_name=prot_feature.res_chain,
            feature_residue_name=prot_feature.res_name,
            feature_residue_number=prot_feature.res_number,
            feature_atom_name=atom_name,
        )
        return feature.id

    @staticmethod
    def _detect(
        pose, protein_system, mol, distance_padding, angle_padding, debug
    ) -> list[dict]:
        """Run the geometric detector, returning candidate interaction dicts."""

        target = pose.target

        # organise ligand features by family
        comp_features_by_family: dict[str, list] = {}
        for f in pose.features:
            comp_features_by_family.setdefault(f.family, []).append(f)

        candidates: list[dict] = []

        for prot_feature in protein_system.get_protein_features():

            prot_family = prot_feature.family

            if prot_family not in COMPLEMENTARY_FEATURES:
                continue

            prot_coords = [a.np_pos for a in prot_feature.atoms]
            if not prot_coords:
                continue
            prot_coord = np.sum(prot_coords, axis=0) / len(prot_coords)

            feature_id = InteractionService._protein_feature_id(target, prot_feature)

            for complementary_family in COMPLEMENTARY_FEATURES[prot_family]:

                interaction_type = INTERACTION_TYPES[
                    (prot_family, complementary_family)
                ]

                for lig_feature in comp_features_by_family.get(
                    complementary_family, []
                ):

                    lig_pos = np.asarray(lig_feature.position)
                    distance = float(np.linalg.norm(lig_pos - prot_coord))
                    angle = None

                    if (
                        distance
                        > INTERACTION_CUTOFF[interaction_type] + distance_padding
                    ):
                        continue

                    lig_coords = None
                    if interaction_type.startswith('π'):
                        conf = mol.GetConformer()
                        lig_coords = [
                            np.array(conf.GetAtomPosition(i - 1))
                            for i in lig_feature.atom_numbers
                        ]

                    if interaction_type == 'π-stacking':
                        # require at least one atom within the min cutoff
                        min_distance = min(
                            float(np.linalg.norm(lig_coord - p_coord))
                            for lig_coord in lig_coords
                            for p_coord in prot_coords
                        )
                        if min_distance > PI_STACK_MIN_CUTOFF + distance_padding:
                            continue

                        lig_norm = _norm([list(p) for p in lig_coords])
                        prot_norm = _norm(prot_coords)
                        angle = _angle_between(lig_norm, prot_norm)

                        # face-to-face has a stricter distance cutoff
                        if (
                            angle < 40 - angle_padding
                            and distance > PI_STACK_F2F_CUTOFF + distance_padding
                        ):
                            continue

                    elif interaction_type == 'π-cation':
                        if prot_family == 'Aromatic':
                            aromatic_norm = _norm(prot_coords)
                            cation_vec = lig_pos - prot_coord
                        else:
                            aromatic_norm = _norm([list(p) for p in lig_coords])
                            cation_vec = prot_coord - lig_pos

                        angle = _angle_between(aromatic_norm, cation_vec)
                        if angle > 30 + angle_padding:
                            continue

                    candidates.append(
                        {
                            'feature_id': feature_id,
                            'feature_family': prot_family,
                            'feature_atom_name': ' '.join(
                                a.name for a in prot_feature.atoms
                            ),
                            'type': interaction_type,
                            'family': lig_feature.family,
                            'atom_ids': [int(i) for i in lig_feature.atom_numbers],
                            'prot_coord': [float(x) for x in prot_coord],
                            'lig_coord': [float(x) for x in lig_pos],
                            'distance': distance,
                            'angle': None if angle is None else float(angle),
                        }
                    )

        if debug:
            mrich.debug(f'{pose}: {len(candidates)} candidate interactions')

        return candidates

    @staticmethod
    def _resolve(candidates: list[dict], debug: bool = False) -> list[dict]:
        """Cull duplicate / less-significant interactions.

        Keeps, per interaction type: the closest interaction per ligand-atom group
        (Hydrogen Bond, π-cation, Electrostatic), the closest per protein feature
        (π-stacking), all Sulfur-Sulfur, and -- for Hydrophobic -- de-duplicates
        lumped vs. single hydrophobes then keeps the closest per protein feature.
        """

        for i, c in enumerate(candidates):
            c['_idx'] = i

        keep: set[int] = set()

        def keep_min_per(predicate, key) -> None:
            """Keep the min-distance candidate within each ``key`` group."""
            best: dict = {}
            for c in candidates:
                if not predicate(c):
                    continue
                k = key(c)
                if k not in best or c['distance'] < best[k]['distance']:
                    best[k] = c
            keep.update(c['_idx'] for c in best.values())

        keep_min_per(
            lambda c: c['type'] == 'Hydrogen Bond', lambda c: tuple(c['atom_ids'])
        )
        keep_min_per(lambda c: c['type'] == 'π-stacking', lambda c: c['feature_id'])
        keep_min_per(
            lambda c: c['type'] == 'π-cation', lambda c: tuple(c['atom_ids'])
        )
        keep_min_per(
            lambda c: c['type'] == 'Electrostatic', lambda c: tuple(c['atom_ids'])
        )

        # Sulfur-Sulfur: keep all
        keep.update(
            c['_idx'] for c in candidates if c['type'] == 'Sulfur-Sulfur'
        )

        # Hydrophobic: de-duplicate lumped vs. single hydrophobes
        hydrophobic = [c for c in candidates if c['type'] == 'Hydrophobic']

        covered: dict = {}
        lumped_lumped: dict = {}
        for c in hydrophobic:
            families = (c['feature_family'], c['family'])
            names = c['feature_atom_name'].split()
            if families == ('LumpedHydrophobe', 'Hydrophobe'):
                for name in names:
                    covered.setdefault((name, c['atom_ids'][0]), []).append(c['_idx'])
            elif families == ('Hydrophobe', 'LumpedHydrophobe'):
                for atom_id in c['atom_ids']:
                    covered.setdefault(
                        (c['feature_atom_name'], atom_id), []
                    ).append(c['_idx'])
            elif families == ('LumpedHydrophobe', 'LumpedHydrophobe'):
                for name in names:
                    for atom_id in c['atom_ids']:
                        covered.setdefault((name, atom_id), []).append(c['_idx'])
                lumped_lumped[c['feature_atom_name']] = tuple(c['atom_ids'])

        keep_hydrophobic = {c['_idx'] for c in hydrophobic}
        for c in hydrophobic:
            families = (c['feature_family'], c['family'])
            if families == ('Hydrophobe', 'Hydrophobe'):
                if (c['feature_atom_name'], c['atom_ids'][0]) in covered:
                    keep_hydrophobic.discard(c['_idx'])
            elif families == ('LumpedHydrophobe', 'Hydrophobe'):
                value = lumped_lumped.get(c['feature_atom_name'])
                if value is not None and c['atom_ids'][0] in value:
                    keep_hydrophobic.discard(c['_idx'])

        # of the surviving hydrophobes, keep the closest per protein feature
        by_idx = {c['_idx']: c for c in candidates}
        best_per_feature: dict = {}
        for idx in keep_hydrophobic:
            c = by_idx[idx]
            k = c['feature_id']
            if k not in best_per_feature or (
                c['distance'] < best_per_feature[k]['distance']
            ):
                best_per_feature[k] = c
        keep.update(c['_idx'] for c in best_per_feature.values())

        resolved = [c for c in candidates if c['_idx'] in keep]
        for c in candidates:
            c.pop('_idx', None)

        if debug:
            mrich.debug(f'resolved {len(candidates)} -> {len(resolved)} interactions')

        return resolved
