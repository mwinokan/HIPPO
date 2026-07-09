import inspect
import json
import logging
import re
import shutil
from collections.abc import Callable
from itertools import combinations
from os.path import relpath
from pathlib import Path
from pprint import pprint
from zipfile import ZipFile

import community as louvain
import mcol
import molparse as mp
import mrich
import networkx as nx
import pandas as pd
from designdb.models import (
    CompoundModel,
    InspirationModel,
    InteractionModel,
    PoseModel,
    PoseTagJunctionModel,
    PoseTagModel,
    ScoreValueModel,
    SubsiteModel,
    SubsiteTagModel,
    TargetModel,
)
from designdb.sets.interaction import InteractionSet
from designdb.settings import DEFAULT_POSE_METHODS
from designdb.utils import ScoreSubquery, normalize_string_list
from designdb.utils_frag import generate_header
from django.conf import settings
from django.db import IntegrityError, transaction
from django.db.models import (
    Count,
    Exists,
    FloatField,
    OuterRef,
    Q,
    QuerySet,
    Subquery,
)
from django.db.models.fields.json import KeyTextTransform
from django.db.models.functions import Cast
from IPython.display import display
from ipywidgets import (
    BoundedIntText,
    Checkbox,
    GridBox,
    Layout,
    VBox,
    interactive,
    interactive_output,
)
from molparse.rdkit import draw_grid, draw_mols
from pandas import DataFrame

# from mypackage.services.compound import CompoundService
from rdkit import Chem

# from rdkit.Chem import inchi
from rdkit.Chem import PandasTools, SDWriter

if settings.MANAGE_MODELS:
    from designdb.utils import JsonGroupArray as ArrayAgg
else:
    from django.contrib.postgres.aggregates import ArrayAgg

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import pandas
    from designdb.components.pose import Pose
    from designdb.sets.compound import CompoundSet


# from .validation.compound import ValidationError, validate_compound_data

SDF_XCAv2_PATTERN = re.compile(
    r'^.*-.\d{4}_._\d*_\d_.*-.\d{4}\+.\+\d*\+\d_ligand\.sdf$'
)
SDF_XCAV3_PATTERN = re.compile(
    r'^.*-.\d{4}_._\d*_._\d_.*-.\d{4}\+.\+\d*\+.\+\d_ligand\.sdf$'
)


SDF_FRAGALYSIS_PATTERN = re.compile(r'^.*\d{4}[a-z].sdf$')
PDBID_PATTERN = re.compile(r'^[A-Za-z0-9]{4}-[a-z].sdf$')


logger = logging.getLogger(__name__)


class PoseSet:
    """Object representing a subset of the 'pose' table in the :class:`.Database`.

    .. attention::

            :class:`.PoseSet` objects should not be created directly. Instead use the
            :meth:`.HIPPO.poses` property. See :doc:`getting_started` and
            :doc:`insert_elaborations`.

    Use as an iterable
    ==================

    Iterate through :class:`.PoseModel` objects in the set:

    ::

            pset = animal.poses[:100]

            for pose in pset:
                    ...

    Check membership
    ================

    To determine if a :class:`.PoseModel` is present in the set:

    ::

            is_member = pose in cset

    Selecting compounds in the set
    ==============================

    The :class:`.PoseSet` can be indexed like standard Python lists by their indices

    ::

            pset = animal.poses[1:100]

            # indexing individual compounds
            pose = pset[0]  # get the first pose
            pose = pset[1]  # get the second pose
            pose = pset[-1] # get the last pose

            # getting a subset of compounds using a slice
            pset2 = pset[13:18] # using a slice

    """

    def __init__(
        self,
        queryset=None,
        *,
        sort: bool = True,
        name: str | None = None,
    ) -> None:
        """PoseSet initialisation"""

        # let's have a queryset no matter what.
        if queryset:
            self._queryset = queryset
        else:
            self._queryset = PoseModel.objects.none()

        self._name = name
        if sort:
            self._queryset = self._queryset.order_by('pk')

        self._interactions = None
        self._metadata_dict = None

    ### DUNDERS

    def __str__(self):
        """Unformatted string representation"""
        if self.name:
            s = f'{self._name}: '
        else:
            s = ''

        s += f'{{P × {len(self)}}}'

        return s

    def __repr__(self) -> str:
        """ANSI Formatted string representation"""
        return f'{mcol.bold}{mcol.underline}{self}{mcol.unbold}{mcol.ununderline}'

    def __rich__(self) -> str:
        """Rich Formatted string representation"""
        return f'[bold underline]{self}'

    def __len__(self) -> int:
        """The number of poses in this set"""
        return self._queryset.count()

    def __iter__(self):
        """Iterate through poses in this set as :class:`.Pose` components"""
        from designdb.components.pose import Pose

        return (Pose(p) for p in self._queryset)

    def __getitem__(
        self,
        key: int | slice,
    ) -> 'Pose | PoseSet':
        """Get poses or subsets thereof from this set

        :param key: integer index or slice of indices

        """
        from designdb.components.pose import Pose

        match key:
            case int():
                # index by position in the (ordered) set, not by pk; support
                # negative indices (e.g. pset[-1] -> last pose)
                n = len(self)
                idx = key + n if key < 0 else key
                if not 0 <= idx < n:
                    raise IndexError(f'PoseSet index out of range: {key}')
                return Pose(self._queryset[idx])

            case slice():
                # positional slice of the (ordered) members. Slice `.all()` (a
                # fresh, unevaluated clone) so this returns a queryset even when
                # self._queryset is already evaluated (which would otherwise slice
                # to a list of instances). sort=False: can't re-order a sliced qs.
                return PoseSet(self._queryset.all()[key], sort=False)

            case _:
                raise NotImplementedError

    def __add__(
        self,
        other: 'PoseSet',
    ) -> 'PoseSet':
        """Add a :class:`.PoseSet` to this set"""
        if isinstance(other, PoseSet):
            return PoseSet(
                PoseModel.objects.filter(
                    Q(pk__in=self._queryset) | Q(pk__in=other.queryset)
                ),
                sort=False,
            )
        elif isinstance(other, PoseModel):
            return PoseSet(
                PoseModel.objects.filter(Q(pk__in=self._queryset) | Q(pk=other.pk)),
                sort=False,
            )
        else:
            raise NotImplementedError

    def __sub__(
        self,
        other: 'PoseSet',
    ) -> 'PoseSet':
        """Substract a :class:`.PoseSet` from this set"""
        match other:
            case PoseSet():
                return PoseSet(
                    PoseModel.objects.filter(
                        Q(pk__in=self._queryset) & ~Q(pk__in=other.queryset)
                    ),
                    sort=False,
                )
            case int():
                return PoseSet(
                    PoseModel.objects.filter(
                        Q(pk__in=self._queryset) & ~Q(pk=other.pk)
                    ),
                    sort=False,
                )

    def __and__(self, other: 'PoseSet'):
        """AND set operation, returns only poses in both sets"""

        match other:
            case PoseSet():
                return PoseSet(
                    PoseModel.objects.filter(
                        Q(pk__in=self._queryset) & Q(pk__in=other.queryset)
                    ),
                    sort=False,
                )

            case _:
                raise NotImplementedError

    def __or__(self, other: 'PoseSet'):
        """OR set operation, returns union of both sets"""

        match other:
            case PoseSet():
                return PoseSet(
                    PoseModel.objects.filter(
                        Q(pk__in=self._queryset) | Q(pk__in=other.queryset)
                    ),
                    sort=False,
                )

            case _:
                raise NotImplementedError

    def __xor__(self, other: 'PoseSet'):
        """Exclusive OR set operation, returns all poses in either set but not both"""

        match other:
            case PoseSet():
                return PoseSet(
                    PoseModel.objects.filter(
                        Q(Q(pk__in=self._queryset) | Q(pk__in=other.queryset))
                        & ~Q(Q(pk__in=self._queryset) & Q(pk__in=other.queryset))
                    ),
                    sort=False,
                )

            case _:
                raise NotImplementedError

    def __call__(
        self,
        *,
        tag: str = None,
        pose_method: str = None,
        target: int = None,
        subsite: int = None,
    ) -> 'PoseSet':
        """Filter poses by a given tag, pose method name, SubsiteModel ID, or
        target ID. See
        :meth:`.PoseSet.get_by_tag`, :meth:`.PoseSet.get_by_method`,
        :meth:`.PoseSet.get_by_target`, and :meth:`.PoseSet.get_by_subsite`"""

        if tag:
            return self.get_by_tag(tag)
        elif pose_method:
            return self.get_by_method(pose_method)
        elif target:
            return self.get_by_target(target=TargetModel.objects.get(pk=target))
        elif subsite:
            return self.get_by_subsite(subsite=SubsiteModel.objects.get(pk=subsite))
        else:
            raise NotImplementedError

    @classmethod
    def get_by_references(cls, poseset: 'PoseSet') -> 'PoseSet':
        return PoseSet(
            PoseModel.objects.filter(pk__in=poseset._queryset.values('pose_reference'))
        )

    # there's a method get_by_inspiration
    @classmethod
    def get_by_inspirations(cls, poseset: 'PoseSet') -> 'PoseSet':
        return PoseSet(
            PoseModel.objects.filter(
                pk__in=InspirationModel.objects.filter(
                    derivative_pose__in=poseset._queryset,
                ).values(
                    'original_pose',
                ),
            ),
        )

    ### FILTERING

    def get_by_tag(
        self,
        tag: str,
        inverse: bool = False,
    ) -> 'PoseSet':
        """Get all child poses with a certain tag

        :param tag: tag to filter by
        :param inverse: return all poses *not* tagged with ``tag``
            (Default value = False)

        """
        self._queryset = self._queryset.annotate(
            has_tag=Exists(
                PoseTagJunctionModel.objects.filter(
                    pose=OuterRef('pk'),
                    pose_tag__pose_tag_name=tag,
                ),
            ),
        )
        if inverse:
            return PoseSet(self._queryset.filter(has_tag=False))
        else:
            return PoseSet(self._queryset.filter(has_tag=True))

    def get_by_method(self, method: str) -> 'PoseSet':
        """Get all poses associated with a given pose method name."""
        return PoseSet(
            self._queryset.filter(
                methods__pose_method_name=method,
            )
        )

    def get_by_metadata(
        self, key: str, value: str | None = None, debug: bool = False
    ) -> 'PoseSet':
        """Get all child poses with by their metadata. If no value is passed, then
        simply containing the key in the metadata dictionary is sufficient

        :param key: metadata key to search for
        :param value: metadata value, if ``None`` return poses with the metadata key
            regardless of value (Default value = None)

        """
        if value is None:
            # metadata stored as string
            return PoseSet(
                self._queryset.filter(pose_metadata__contains=f'"{key}"'),
            )

        else:
            if isinstance(value, str):
                value = f'"{value}"'

            return PoseSet(
                self._queryset.filter(pose_metadata__contains=f'"{key}: {value}"'),
            )

    def get_by_inspiration(self, inspiration: PoseModel, inverse: bool = False):
        """Get all child poses with with this inspiration.

        :param inspiration: inspiration :class:`.PoseModel` ID or object
        :param inverse: invert the selection (Default value = False)

        """
        # not entirely sure which way the filtering should go
        qs = (
            InspirationModel.objects.filter(
                derivative_pose=inspiration,
            ).values('original_pose'),
        )

        if inverse:
            return PoseSet(self._queryset.exclude(pk__in=qs))
        else:
            return PoseSet(self._queryset.filter(pk__in=qs))

    def get_df(
        self,
        smiles: bool = True,
        inchikey: bool = True,
        alias: bool = True,
        name: bool = True,
        compound_id: bool = False,
        target_id: bool = False,
        reference_id: bool = False,
        reference_alias: bool = False,
        path: bool = False,
        mol: bool = False,
        energy_score: bool = False,
        distance_score: bool = False,
        inspiration_score: bool = False,
        metadata: bool = False,
        expand_metadata: bool = True,
        debug: bool = True,
        inspiration_ids: bool = False,
        inspiration_aliases: bool = False,
        derivative_ids: bool = False,
        tags: bool = False,
        expand_tags: bool = False,
        subsites: bool = False,
        scoring_methods: list[tuple[str, str]] | None = None,
        pose_method: bool = False,
        # skip_no_mol=True, reference: str = "name", mol: bool = False, **kwargs
    ) -> 'pandas.DataFrame':
        """Get a DataFrame of the poses in this set.

        :param smiles: include SMILES column (Default value = True)
        :param inchikey: include InChIKey column (Default value = True)
        :param alias: include alias column (Default value = True)
        :param name: include name column (Default value = True)
        :param compound_id: include :class:`.CompoundModel` ID column
            (Default value = False)
        :param reference_id: include reference :class:`.PoseModel` ID column
            (Default value = False)
        :param target_id: include reference :class:`.TargetModel` ID column
            (Default value = False)
        :param path: include path column (Default value = False)
        :param mol: include ``rdkit.Chem.Mol`` in output (Default value = False)
        :param energy_score: include energy_score column (Default value = False)
        :param distance_score: include distance_score column (Default value = False)
        :param inspiration_score: include inspiration_score column
            (Default value = False)
        :param metadata: include metadata in output (Default value = False)
        :param expand_metadata: create separate column for each metadata key
            (Default value = True)
        :param inspiration_ids: include inspiration :class:`.PoseModel` ID column
        :param inspiration_aliases: include inspiration :class:`.PoseModel` alias column
        :param derivative_ids: include derivative :class:`.PoseModel` ID column
        :param tags: include tags column
        :param subsites: include subsites column
        """

        sig = inspect.signature(self.get_df)
        flags = {
            name: locals()[name]
            for name in sig.parameters
            if name
            not in (
                'self',
                'debug',
                'expand_tags',
                'expand_metadata',
                'scoring_methods',
            )
        }
        # need id in output
        flags['id'] = True

        # alias and name both point to same thing. prefer 'name'
        if flags.get('name', False):
            flags['alias'] = True

        # this is still not working right and I don't understand. What
        # was the original code doing here? simply adding both fields,
        # name and alias?

        # dict :: func arg: (col title, qs field lookup, queryset annotation)
        # this is going to get out of hand with multiple scoring methods
        fields = {
            'id': ('id', 'id', None),
            'smiles': ('smiles', 'pose_smiles', None),
            'inchikey': ('inchikey', 'pose_inchikey', None),
            # 'alias': ('alias', 'pose_alias', None),
            'name': ('name', 'pose_alias', None),
            'compound_id': ('compound_id', 'compound__id', None),
            'target_id': ('target_id', 'target__id', None),
            'reference_id': ('reference_id', 'pose_reference', None),
            'reference_alias': (
                'reference_alias',
                'reference_alias',
                Subquery(
                    PoseModel.objects.filter(
                        pk=OuterRef('pose_reference'),
                    ).values('pose_alias')[0:1]
                ),
            ),
            'path': ('protein_link', 'protein_link', None),
            'mol': ('mol', 'pose_mol', None),
            'energy_score': (
                'energy_score',
                'energy_score',
                ScoreSubquery('energy_score'),
            ),
            'distance_score': (
                'distance_score',
                'distance_score',
                ScoreSubquery('distance_score'),
            ),
            'inspiration_score': (
                'inspiration_score',
                'inspiration_score',
                ScoreSubquery('inspiration_score'),
            ),
            'metadata': ('metadata', 'pose_metadata', None),
            'inspiration_ids': (
                'inspiration_ids',
                'inspiration_ids',
                ArrayAgg('inspirations__id'),
                # JsonGroupArray('inspirations__id'),
            ),
            'inspiration_aliases': (
                'inspiration_aliases',
                'inspiration_aliases',
                ArrayAgg('inspirations__pose_alias'),
                # JsonGroupArray('inspirations__pose_alias'),
            ),
            'derivative_ids': (
                'derivative_ids',
                'derivative_ids',
                ArrayAgg('inspirations__id'),
                # JsonGroupArray('inspirations__id'),
            ),
            'tags': (
                'tags',
                'tag_names',
                ArrayAgg('tags__pose_tag_name'),
                # JsonGroupArray('tags__pose_tag_name'),
            ),
            'subsites': (
                'subsites',
                'subsites_names',
                ArrayAgg(
                    'subsites__subsite_name',
                    filter=Q(subsites__isnull=False),
                ),
                # JsonGroupArray('subsites__subsite_name',
                #     filter=Q(subsites__isnull=False),),
            ),
            'pose_method': (
                'pose_method',
                'pose_method_names',
                ArrayAgg('methods__pose_method_name', distinct=True),
            ),
        }

        annotations = {
            v[1]: v[2] for k, v in fields.items() if flags.get(k, False) and v[2]
        }
        values = [v[1] for k, v in fields.items() if flags.get(k, False)]
        columns = {v[1]: v[0] for k, v in fields.items() if flags.get(k, False)}

        for method_name, version in scoring_methods or []:
            score_sq = Subquery(
                ScoreValueModel.objects.filter(
                    pose=OuterRef('pk'),
                    compound=OuterRef('compound'),
                    scoring_method__method_name=method_name,
                    scoring_method__method_version=version,
                )
                .annotate(
                    score_val=Cast(
                        KeyTextTransform('score', 'score'), output_field=FloatField()
                    )
                )
                .values('score_val')[:1],
                output_field=FloatField(),
            )
            annotations[method_name] = score_sq
            values.append(method_name)
            columns[method_name] = method_name

        qs = self._queryset.annotate(**annotations).values(*values)

        df = pd.DataFrame(qs)
        df = df.rename(columns=columns)
        df = df.set_index('id')

        if alias:
            df['alias'] = df.name

        if metadata and expand_metadata:
            # TODO: code specific to my current situation. have to
            # parse string to json (does postgres handle this
            # automatically?)
            # expanded = pd.json_normalize(
            #     df["metadata"].apply(lambda x: json.loads(x) if x else {}),
            # )
            expanded = pd.json_normalize(df['metadata'])
            # dropping columns is due to confusion with scores. can't be
            # permanent solution, for now, drop the common ones
            expanded = expanded.drop(
                columns=set(expanded.columns).intersection(set(df.columns)),
            )

            df = df.drop(columns=['metadata']).join(expanded)

        if tags and expand_tags:
            # surprisingly manual compared to expand_metadata, but
            # kept running into problems
            df['tags'] = df['tags'].apply(normalize_string_list)
            # get all unique tags
            all_tags = sorted(set(tag for tags in df['tags'] for tag in tags))

            # build boolean columns
            for tag in all_tags:
                df[tag] = df['tags'].apply(lambda tags, tag=tag: tag in tags)

            df = df.drop(columns=['tags'])

        # custom aggreagte field is giving me string, parse to list
        for col in [
            'inspiration_aliases',
        ]:
            if col in df.columns:
                df[col] = df[col].apply(normalize_string_list)

        return df

    def get_by_reference(
        self,
        ref_id: int,
    ) -> 'PoseSet | None':
        """Get poses with a certain reference id

        :param ref_id: reference :class:`.PoseModel` ID

        """
        qs = self._queryset.filter(pose_reference=ref_id)
        if not qs.exists():
            # odd, but keeping now
            return None

        return PoseSet(qs)

    def get_by_compound(
        self,
        *,
        compound: 'int | CompoundModel | CompoundSet',
    ) -> 'PoseSet | None':
        """Select a subset of this :class:`.PoseSet` by the associated
        :class:`.CompoundModel`.

        :param compound: :class:`.CompoundModel` object or ID
        :returns: a :class:`.PoseSet` of the selection

        """
        if isinstance(compound, int):
            return PoseSet(self._queryset.filter(compound__id=compound))
        elif isinstance(compound, CompoundModel):
            return PoseSet(self._queryset.filter(compound=compound))
        else:
            # possible crash point: assuming CompoundSet but not
            # testing type, still trying to fiugre out circular
            # imports
            return PoseSet(self._queryset.filter(compound__in=compound.queryset))

    def get_by_target(
        self,
        *,
        target: TargetModel,
    ) -> 'PoseSet | None':
        """Select a subset of this :class:`.PoseSet` by the associated
        :class:`.TargetModel`.

        :param id: :class:`.TargetModel` ID
        :returns: a :class:`.PoseSet` of the selection

        """
        # where would you need this method?? do you ever create sets
        # of poses from different targets?
        return PoseSet(self._queryset.filter(target=target))

    def get_by_subsite(
        self,
        *,
        subsite: SubsiteModel,
    ) -> 'PoseSet | None':
        """Select a subset of this :class:`.PoseSet` by the associated
        :class:`.SubsiteModel`.

        :param id: :class:`.SubsiteModel` ID
        :returns: a :class:`.PoseSet` of the selection

        """
        qs = self._queryset.filter(
            id__in=SubsiteTagModel.objects.filter(
                subsite=subsite,
            ).values('pose'),
        )

        if self.name:
            name = f'{self.name} & subsite={subsite.pk}'
        else:
            name = None

        return PoseSet(qs, name=name)

    def set_subsites_from_metadata_field(self, field: str = 'CanonSites alias') -> None:
        """Create and assign subsite entries from a pose metadata field."""
        # local import: keeps this upward set -> service call off the module graph
        from designdb.services.subsite import SubsiteService

        SubsiteService.set_subsites_from_metadata_field(self._queryset, field)

    def get_best_scoring_poses_per_compound(
        self,
        scoring_method: str,
        version: str | None = None,
        inverse: bool = False,
    ) -> 'PoseSet':
        """Return one pose per compound with the best score for the given
        scoring method.

        :param scoring_method: ``ScoringMethodModel.method_name`` to rank by
        :param version: ``ScoringMethodModel.method_version`` — required when multiple
            versions of the same method exist
        :param inverse: if ``True``, higher score is better (default: lower is better)
        """
        score_num = Cast(KeyTextTransform('score', 'score'), output_field=FloatField())

        filters = {'scoring_method__method_name': scoring_method}
        if version is not None:
            filters['scoring_method__method_version'] = version

        rows = (
            ScoreValueModel.objects.filter(pose__in=self._queryset, **filters)
            .annotate(score_num=score_num)
            .values_list('compound_id', 'pose_id', 'score_num')
        )

        # keep exactly one pose per compound (the best score; ties broken by
        # first-seen), rather than every pose tied for the best
        best: dict[int, tuple[int, float]] = {}
        for compound_id, pose_id, score in rows:
            if score is None:
                continue
            current = best.get(compound_id)
            if current is None or (
                score > current[1] if inverse else score < current[1]
            ):
                best[compound_id] = (pose_id, score)

        best_pose_ids = [pose_id for pose_id, _ in best.values()]
        return PoseSet(PoseModel.objects.filter(pk__in=best_pose_ids))

    # def filter(
    #     self,
    #     function=None,
    #     *,
    #     key: str = None,
    #     value: str = None,
    #     operator='=',
    #     inverse: bool = False,
    # ):
    #     """Filter this :class:`.PoseSet` by selecting members where
    #     ``function(pose)`` is truthy or pass a key, value, and optional operator
    #     to search by database values

    #     :param function: callable object
    #     :param key: database field for 'pose' table ('pose_' prefix not needed)
    #     :param value: value to compare to
    #     :param operator: comparison operator (default = "=")
    #     :param inverse: invert the selection (Default value = False)

    #     """

    #     if function:
    #         ids = set()
    #         for pose in self:
    #             value = function(pose)
    #             # mrich.debug(f'{pose=} {value=}')
    #             if value and not inverse:
    #                 ids.add(pose.id)
    #             elif not value and inverse:
    #                 ids.add(pose.id)

    #         return PoseSet(self.db, ids)

    #     sql = f"""
    #     SELECT pose_id FROM {self.db.SQL_SCHEMA_PREFIX}pose
    #     WHERE pose_id IN {self.str_ids}
    #     AND pose_{key} {operator} {value}
    #     """

    #     cursor = self.db.execute(sql)

    #     ids = [i for (i,) in cursor]

    #     return PoseSet(self.db, ids)

    def add_tag(
        self,
        tag: str,
    ) -> None:
        """Add this tag to every member of the set"""

        assert isinstance(tag, str)

        pose_tag = PoseTagModel(pose_tag_name=tag)
        pose_tag.save()

        PoseTagJunctionModel.objects.bulk_create(
            [
                PoseTagJunctionModel(pose=pose, pose_tag=pose_tag)
                for pose in self._queryset
            ],
            ignore_conflicts=True,
        )

        mrich.print(f'Tagged {self} w/ "{tag}"')

        # refetch in case was evaluated
        self._queryset = PoseModel.objects.filter(pk__in=self._queryset.values('pk'))

        # NB! I'm now realizing this is potentially a huge
        # problem. with every evaluation and refretch some attributes
        # may be lost. how can this be kept clean?

    def append_to_metadata(
        self,
        key,
        value,
    ) -> None:
        """Append ``value`` to the list at ``key`` in every member's metadata dict.

        :param key: the metadata key to match (created as a new list if absent)
        :param value: the value to append to the list
        """
        for pose in self._queryset:
            metadata = pose.pose_metadata or {}
            existing = metadata.get(key)
            if existing is None:
                metadata[key] = [value]
            elif isinstance(existing, list):
                existing.append(value)
            else:
                mrich.error(f'Could not append to metadata {key=}: not a list')
                continue
            pose.pose_metadata = metadata
            pose.save()
        self._queryset = PoseModel.objects.filter(pk__in=self._queryset.values('pk'))

    # TODO: implement scores
    # def calculate_inspiration_scores(
    #     self,
    #     alpha: float = 0.95,
    #     beta: float = 0.05,
    #     score_type: str = 'combo',
    # ) -> 'pd.DataFrame':
    #     """Set inspiration_score values using MoCASSIn.calculate_mocassin_tversky

    #     :param alpha: Tversky alpha parameter
    #     :param beta: Tversky beta parameter
    #     :param score_type: Score type to add to database, choose from
    #         "combo", "shape", "colour"
    #     :returns: Pandas DataFrame with molecules and scores
    #     """

    #     from mocassin.mocassin import calculate_mocassin_tversky

    #     df = self.get_df(
    #         alias=False,
    #         smiles=False,
    #         inchikey=False,
    #         inspiration_ids=True,
    #         mol=True,
    #     )

    #     inspirations = {p.id: p for p in self.inspirations}

    #     df['inspiration_mols'] = df['inspiration_ids'].apply(
    #         lambda x: [inspirations[i].mol for i in x]
    #     )

    #     n = len(df)

    #     for j, (i, row) in mrich.track(
    #         enumerate(df.iterrows()), prefix='MoCASSIn', total=n
    #     ):
    #         mrich.set_progress_field('j', j)
    #         mrich.set_progress_field('n', n)

    #         try:
    #             combo, shape, colour = calculate_mocassin_tversky(
    #                 row['inspiration_mols'],
    #                 row['mol'],
    #                 alpha=0.95,
    #                 beta=0.05,
    #             )
    #             df.loc[i, f'mocassin_combo({alpha},{beta})'] = combo
    #             df.loc[i, f'mocassin_shape({alpha},{beta})'] = shape
    #             df.loc[i, f'mocassin_colour({alpha},{beta})'] = colour
    #         except Exception as e:
    #             mrich.error(e)

    #     tuples = df[f'mocassin_{score_type}({alpha},{beta})'].items()

    #     sql = f"""UPDATE {self.db.SQL_SCHEMA_PREFIX}pose
    #     SET pose_inspiration_score = {self.db.SQL_STRING_PLACEHOLDER}
    #     WHERE pose_id = {self.db.SQL_STRING_PLACEHOLDER}"""

    #     mrich.debug('Updating pose_inspiration_score values')
    #     self.db.executemany(sql, [(b, a) for a, b in tuples])
    #     self.db.commit()

    #     return df

    ### SPLITTING

    def split_by_reference(self) -> 'dict[int,PoseSet]':
        """Split this :class:`.PoseSet` into subsets grouped by reference ID

        :returns: a dictionary with reference :class:`.PoseModel` IDs as keys and
            :class:`.PoseSet` subsets as values

        """
        sets = {}
        for ref_id in self.reference_ids:
            sets[ref_id] = self.get_by_reference(ref_id)
        return sets

    def split_by_inspirations(
        self,
        single_set: bool = False,
    ) -> 'dict[PoseSet,PoseSet] | PoseSet':
        """Split this :class:`.PoseSet` into subsets grouped by inspirations

        :param single_set: Return a single :class:`.PoseSet` with members sorted by
            inspirations (Default value = False)
        :returns: a dictionary with tuples of inspiration :class:`.PoseSet` as keys and
            :class:`.PoseSet` derivative subsets as values

        """

        sets = {}

        for pose in self._queryset:
            insp_ids = list(pose.inspirations.distinct().values_list('pk', flat=True))
            key = tuple(insp_ids)
            sets.setdefault(key, set())
            sets[key].add(pose.pk)

        mrich.var('#unique inspiration combinations', len(sets))

        if single_set:
            return PoseSet(
                PoseModel.objects.filter(
                    pk__in=[id for s in sets.values() for id in s.ids],
                    sort=False,
                )
            )

        self._queryset = PoseModel.objects.filter(pk__in=self._queryset.values('pk'))

        return {
            PoseSet(PoseModel.objects.filter(pk__in=insp_ids)): PoseSet(
                PoseModel.objects.filter(pk__in=pose_ids)
            )
            for insp_ids, pose_ids in sets.items()
        }

    ### EXPORTING

    def write_sdf(
        self,
        out_path: str,
        name_col: str = 'alias',
        inspiration_ids: bool = False,
        inspiration_aliases: bool = False,
        **kwargs,
    ) -> None:
        """Write an SDF

        :param out_path: filepath of the output
        :param name_col: pose property to use as the name column, can be
            ``["name", "alias", "inchikey", "id"]`` (Default value = 'name')
        :param inspiration_ids: include inspiration :class:`.PoseModel` ID column
        :param inspiration_aliases: include inspiration :class:`.PoseModel` alias column
        :param fragalysis_inspirations: create inspirations column "ref_mols"
        """

        df = self.get_df(
            mol=True,
            inspiration_ids=inspiration_ids,
            inspiration_aliases=inspiration_aliases,
            name=name_col == 'name',
            **kwargs,
        )

        if name_col not in ['name', 'alias', 'inchikey', 'id']:
            # try getting name from metadata
            records = self._queryset.values('id', 'pose_metadata')

            longcode_lookup = {}
            for i, d in records:
                if d:
                    metadata = json.loads(d)
                else:
                    metadata = {}

                longcode_lookup[i] = metadata.get(name_col, None)

            values = []
            for _, row in df.iterrows():
                values.append(longcode_lookup[row['id']])

            df[name_col] = values

        df = df.rename(columns={name_col: '_Name', 'mol': 'ROMol'})

        mrich.writing(out_path)

        PandasTools.WriteSDF(df, out_path, 'ROMol', '_Name', list(df.columns))

        # keep record of export
        value = str(Path(out_path).resolve())
        # self.db.remove_metadata_list_item(table='pose', key='exports', value=value)
        self.append_to_metadata(key='exports', value=value)

    def to_fragalysis(
        self,
        out_path: str,
        *,
        method: str,
        ref_url: str = 'https://hippo.winokan.com',
        submitter_name: str,
        submitter_email: str,
        submitter_institution: str,
        metadata: bool = True,
        sort_by: str | None = None,
        sort_reverse: bool = False,
        generate_pdbs: bool = False,
        copy_reference_pdbs: bool = False,
        # ingredients: IngredientSet = None,
        skip_no_reference: bool = True,
        skip_no_inspirations: bool = True,
        skip_metadata: list[str] | None = None,
        tags: bool = True,
        subsites: bool = True,
        extra_cols: dict[str, list] = None,
        inspiration_score: bool = True,
        # name_col: str = "name",
        **kwargs,
    ):
        """Prepare an SDF for upload to the RHS of Fragalysis.

        :param out_path: the file path to write to
        :param method: method used to generate the compounds
        :param ref_url: reference URL for the method
        :param submitter_name: name of the person submitting the compounds
        :param submitter_email: email of the person submitting the compounds
        :param submitter_institution: institution name of the person submitting the
            compounds
        :param metadata: include metadata in the output? (Default value = True)
        :param skipmetadata: exclude metadata keys from output
        :param sort_by: if set will sort the SDF by this column/field
            (Default value = None)
        :param sort_reverse: reverse the sorting (Default value = False)
        :param generate_pdbs: generate accompanying protein-ligand complex PDBs
            (Default value = False)
        :param ingredients: get procurement and amount information from this
            :class:`.IngredientSet` (Default value = None)
        :param tags: include a column for tags in the output (Default value = True)
        :param subsites: include a column for subsites in the output
            (Default value = True)
        :param extra_cols: extra_cols should be a dictionary with a key for each column
            name, and list values where the first element is the field description, and
            all subsequent elements are values for each pose.

        """

        assert out_path.endswith('.sdf')

        _name_col = '_Name'
        mol_col = 'ROMol'
        mol_col = 'mol'

        # make sure references are defined:
        logger.debug('entering')

        mrich.debug(len(self), 'poses in set')
        poses = None

        if skip_no_reference:
            values = self._queryset.filter(pose_reference__isnull=False)

            if not values.exists():
                mrich.debug('no references, quitting')
                logger.warning('no references, quitting')
                return

            poses = PoseSet(values)

            mrich.debug(len(poses), 'remaining after skipping null reference')

        if skip_no_inspirations:
            if not poses:
                poses = self

            values = InspirationModel.objects.filter(
                derivative_pose__in=self._queryset,
            ).values(
                'derivative_pose',
            )

            if values.exists():
                poses = PoseSet(PoseModel.objects.filter(pk__in=values))
                mrich.debug(len(poses), 'remaining after skipping null inspirations')
            else:
                logger.warning(
                    'no inspirations found; per-pose fallback will set '
                    'inspiration to self'
                )

        if not poses:
            # huh?
            poses = PoseSet(self._queryset)

        mrich.var('#poses', len(poses))

        # TODO: this should not go through the df
        pose_df = poses.get_df(
            mol=True,
            inspiration_ids=True,
            # duplicate_name="original ID",
            name=True,
            compound_id=True,
            reference_id=True,
            metadata=metadata,
            tags=tags,
            subsites=subsites,
            energy_score=True,
            distance_score=True,
            inspiration_score=inspiration_score,
            # sanitise_null_metadata_values=True,
            expand_tags=False,
            # sanitise_tag_list_separator=";",
            # sanitise_metadata_list_separator=";",
            # skip_metadata=skip_metadata,
            # **kwargs,
        )

        pose_df = pose_df.reset_index()

        # fix inspirations and reference column (comma separated aliases).
        # reference poses may lie outside this set, so look up exactly the
        # referenced pose IDs rather than scanning the whole table.
        ref_ids = {r for r in pose_df['reference_id'].tolist() if r is not None}
        lookup = dict(
            PoseModel.objects.filter(pk__in=ref_ids).values_list('pk', 'pose_alias')
        )

        # for i, row in pose_df.iterrows():
        #     strs = []
        #     for i in normalize_string_list(row['inspiration_ids']):
        #         # this is what it did in original code
        #         alias = self._queryset.get(pk=i).pose_alias
        #         if not alias:
        #             continue
        #         strs.append(alias)
        #     inspiration_strs.append(','.join(strs))

        # comma separate subsites
        if subsites:

            def fix_subsites(subsite_list: list[str]):
                """Fix subsites"""
                if not subsite_list:
                    logger.warning('no subsite list')
                    return 'None'
                return ','.join(subsite_list)

            pose_df['subsites'] = pose_df['subsites'].apply(fix_subsites)

        if tags:
            pose_df['tags'] = pose_df['tags'].apply(
                lambda x: ','.join(v for v in x if v is not None)
            )

        # pose_df['ref_mols'] = inspiration_strs
        pose_df['ref_mols'] = 'inspiration_strs'
        pose_df['ref_pdb'] = pose_df['reference_id'].apply(lambda x: lookup.get(x))

        # add compound identifier column (inchikey?)

        drops = ['inspiration_ids', 'reference_id']

        # if ingredients:
        #     drops.pop(drops.index("compound"))

        if skip_no_reference:
            prev = len(pose_df)
            pose_df = pose_df[pose_df['reference_id'].notna()]
            if len(pose_df) < prev:
                mrich.warning(f'Skipping {prev - len(pose_df)} Poses with no reference')

        pose_df = pose_df.drop(columns=drops, errors='ignore')

        pose_df[_name_col] = pose_df['name']

        pose_df.rename(
            inplace=True,
            columns={
                'id': 'HIPPO PoseModel ID',
                'compound_id': 'HIPPO CompoundModel ID',
                'mol': mol_col,
                # "smiles": "original SMILES",
                # "compound_id": "compound inchikey",
            },
        )

        extras = {
            'HIPPO PoseModel ID': 'HIPPO PoseModel ID',
            'HIPPO CompoundModel ID': 'HIPPO CompoundModel ID',
            'smiles': 'smiles',
            'ref_pdb': 'protein reference',
            'ref_mols': 'fragment inspirations',
            'alias': 'alias',
            # "compound inchikey": "compound inchikey",
            'distance_score': 'distance_score',
            'energy_score': 'energy_score',
            'inspiration_score': 'inspiration_score',
        }

        if subsites:
            extras['subsites'] = 'subsites'

        if tags:
            extras['tags'] = 'tags'

        if extra_cols:
            for key, value in extra_cols.items():
                extras[key] = value[0]

        # if ingredients:

        #     q_entries = []
        #     q_prices = []
        #     q_lead_times = []
        #     q_amounts = []

        #     currency = None

        #     for i, row in pose_df.iterrows():

        #         compound_id = self.db.get_compound_id(
        #             inchikey=row["compound inchikey"])

        #         ingredient = ingredients(compound_id=compound_id)

        #         if isinstance(ingredient, IngredientSet):
        #             ingredient = sorted(
        #                 [i for i in ingredient], key=lambda x: x.quote.price
        #             )[0]

        #         quote = ingredient.quote
        #         if not currency:
        #             currency = quote.currency
        #         else:
        #             assert quote.currency == currency

        #         q_entries.append(quote.entry_str)
        #         q_prices.append(quote.price)
        #         q_lead_times.append(quote.lead_time)
        #         q_amounts.append(quote.amount)

        #     pose_df["Supplier Catalogue Entry"] = q_entries
        #     # pose_df['Supplier:Catalogue:Entry'] = q_entries
        #     pose_df[f"Price ({currency})"] = q_prices
        #     pose_df["Lead time (working days)"] = q_lead_times
        #     pose_df["Amount (mg)"] = q_amounts

        #     extras["Supplier Catalogue Entry"] = "Supplier Catalogue Entry string"
        #     extras[f"Price ({currency})"] = "Quoted price"
        #     extras["Lead time (working days)"] = "Quoted lead-time"
        #     extras["Amount (mg)"] = "Quoted amount"

        out_path = Path(out_path).resolve()
        mrich.var('out_path', out_path)

        if generate_pdbs:
            # output subdirectory
            out_key = Path(out_path).name.removesuffix('.sdf')
            pdb_dir = Path(out_path).parent / Path(out_key)
            pdb_dir.mkdir(exist_ok=True)
            zip_path = Path(out_path).parent / f'{out_key}_pdbs.zip'

            # create the zip archive
            with ZipFile(str(zip_path.resolve()), 'w') as z:
                # loop over poses
                for (i, row), pose in zip(pose_df.iterrows(), poses, strict=False):
                    # filenames
                    pdb_name = f'{out_key}_{row._Name}.pdb'
                    pdb_path = pdb_dir / pdb_name
                    pose_df.loc[i, 'ref_pdb'] = pdb_name

                    # generate the PL-complex
                    sys = pose.complex_system

                    # write the PDB
                    mrich.writing(pdb_path)
                    sys.write(pdb_path, verbosity=0)
                    z.write(pdb_path)

            mrich.writing(f'{out_key}_pdbs.zip')

        if copy_reference_pdbs:
            # output subdirectory
            out_key = Path(out_path).name.removesuffix('.sdf')
            pdb_dir = Path(out_path).parent / Path(out_key)
            pdb_dir.mkdir(exist_ok=True)
            zip_path = Path(out_path).parent / f'{out_key}_refs.zip'

            # ref_pdb holds reference pose aliases that may lie outside this set,
            # so look them up directly (alias -> protein PDB path)
            ref_aliases = {a for a in pose_df['ref_pdb'].tolist() if a is not None}
            lookup = dict(
                PoseModel.objects.filter(pose_alias__in=ref_aliases).values_list(
                    'pose_alias', 'protein_link'
                )
            )

            zips = set()
            for ref_alias in pose_df['ref_pdb'].values:
                source = lookup.get(ref_alias)
                if not source:
                    mrich.warning(
                        f'No protein file for reference {ref_alias!r}; skipping'
                    )
                    continue
                source_path = Path(source)

                stem = source_path.name.replace('_hippo.pdb', '.pdb')
                # current Fragalysis protein-file naming
                apo_path = source_path.parent / stem.replace(
                    '.pdb', '_delig-desolv.pdb'
                )
                # DEPRECATED(apo-naming): fall back to pre-'delig' naming if present,
                # remove once all data uses 'delig'
                legacy_path = source_path.parent / stem.replace(
                    '.pdb', '_apo-desolv.pdb'
                )
                if not apo_path.exists() and legacy_path.exists():
                    apo_path = legacy_path

                if not apo_path.exists():
                    sys = mp.parse(source_path).protein_system
                    sys.write(apo_path, verbosity=0)

                target_path = pdb_dir / f'{ref_alias}.pdb'

                if not target_path.exists():
                    mrich.writing(target_path)
                    shutil.copy(apo_path, target_path)

                zips.add(target_path)

            # create the zip archive
            with ZipFile(str(zip_path.resolve()), 'w') as z:
                for path in zips:
                    z.write(path, arcname=path.name)

            mrich.writing(f'{out_key}_refs.zip')

        # create the header molecule

        header = generate_header(
            # self[0],   # <- what does that do??
            self._queryset.first(),
            method=method,
            ref_url=ref_url,
            submitter_name=submitter_name,
            submitter_email=submitter_email,
            submitter_institution=submitter_institution,
            extras=extras,
            metadata=metadata,
        )

        # # empty properties
        # pose_df["generation_date"] = [None] * len(pose_df)
        # pose_df["submitter_name"] = [None] * len(pose_df)
        # pose_df["method"] = [None] * len(pose_df)
        # pose_df["submitter_email"] = [None] * len(pose_df)
        # pose_df["ref_url"] = [None] * len(pose_df)

        if extra_cols:
            for key, value in extra_cols.items():
                if len(value) != len(pose_df) + 1:
                    mrich.error(
                        f'extra_col "{key}" does not have the correct number of values'
                    )
                    raise ValueError(
                        f'extra_col "{key}" does not have the correct number of values'
                    )
                pose_df[key] = value[1:]

        if sort_by:
            pose_df = pose_df.sort_values(by=sort_by, ascending=not sort_reverse)

        mrich.writing(out_path)

        with open(out_path, 'w') as sdfh:
            with SDWriter(sdfh) as w:
                w.write(header)
            PandasTools.WriteSDF(
                pose_df, sdfh, mol_col, _name_col, set(pose_df.columns)
            )

        # keep record of export
        value = str(Path(out_path).resolve())

        # FIXME
        # self.db.remove_metadata_list_item(table='pose', key='exports', value=value)

        self.append_to_metadata(key='exports', value=value)

        return pose_df

    def to_pymol(self, prefix: str | None = None) -> None:
        """Group the poses by reference protein and inspirations and output relevant
        PDBs and SDFs.

        :param prefix: prefix to give all output subdirectories (Default value = None)

        """

        commands = []

        prefix = prefix or ''
        if prefix:
            prefix = f'{prefix}_'

        from pathlib import Path

        for ref_id, poses in self.split_by_reference().items():
            ref_pose = PoseModel.objects.get(id=ref_id)
            ref_name = ref_pose.pose_alias or ref_id

            # create the subdirectory
            ref_dir = Path(f'{prefix}ref_{ref_name}')
            mrich.writing(ref_dir)
            ref_dir.mkdir(parents=True, exist_ok=True)

            # write the reference protein
            ref_pdb = ref_dir / f'ref_{ref_name}.pdb'
            ref_pose.protein_system.write(ref_pdb, verbosity=0)

            # color the reference:
            commands.append(f'load {ref_pdb.resolve()}')
            commands.append('hide')
            commands.append('show lines')
            commands.append('show surface')
            commands.append('util.cbaw')
            commands.append('set surface_color, white')
            commands.append('set transparency,  0.4')

            for j, (inspirations, poses) in enumerate(  # noqa: B020
                poses.split_by_inspirations().items()
            ):
                # split_by_inspirations() keys are already inspiration PoseSets
                insp_names = '-'.join(inspirations.names)

                # create the subdirectory
                insp_dir = ref_dir / insp_names
                insp_dir.mkdir(parents=True, exist_ok=True)

                # write the inspirations
                insp_sdf = insp_dir / f'{insp_names}_frags.sdf'
                inspirations.write_sdf(insp_sdf)

                commands.append(f'load {insp_sdf.resolve()}')
                commands.append(
                    f'set all_states, on, {insp_sdf.name.removesuffix(".sdf")}'
                )
                commands.append(f'util.rainbow "{insp_sdf.name.removesuffix(".sdf")}"')

                # write the poses
                pose_sdf = insp_dir / f'{insp_names}_derivatives.sdf'
                poses.write_sdf(pose_sdf)

                commands.append(f'load {pose_sdf.resolve()}')
                commands.append(f'util.cbaw "{pose_sdf.name.removesuffix(".sdf")}"')

                if j > 0:
                    commands.append(f'disable "{insp_sdf.name.removesuffix(".sdf")}"')
                    commands.append(f'disable "{pose_sdf.name.removesuffix(".sdf")}"')

        return '; '.join(commands)

    def to_knitwork(
        self, out_path: str, path_root: str = '.', aligned_files_dir: str | None = None
    ) -> None:
        """Knitwork takes a CSV input with:

        - observation shortcode
        - smiles
        - path_to_ligand_mol
        - path_to_pdb

        :param out_path: path to output CSV
        :param path_root: paths in CSV will be relative to here

        """

        out_path = Path(out_path).resolve()
        path_root = Path(path_root).resolve()
        mrich.var('out_path', out_path)
        mrich.var('path_root', path_root)
        mrich.var('aligned_files_dir', aligned_files_dir)

        assert out_path.name.endswith('.csv')

        with open(out_path, 'w') as f:
            mrich.writing(out_path)

            for pose in self._queryset:
                assert pose.pose_alias
                assert pose.methods.filter(
                    pose_method_name__in=DEFAULT_POSE_METHODS
                ).exists()

                if aligned_files_dir:
                    mol = str(pose.mol_path)
                    pdb = str(pose.apo_path)

                    assert 'aligned_files' in mol
                    assert 'aligned_files' in pdb

                    mol = mol.split('aligned_files/')[-1]
                    pdb = pdb.split('aligned_files/')[-1]

                    aligned_files_dir = Path(aligned_files_dir)

                    mol = relpath(aligned_files_dir / mol, path_root)
                    pdb = relpath(aligned_files_dir / pdb, path_root)

                else:
                    mol = relpath(pose.mol_path, path_root)
                    pdb = relpath(pose.apo_path, path_root)

                data = [pose.pose_alias, pose.compound.compound_smiles, mol, pdb]

                f.write(','.join(data))
                f.write('\n')

    def to_syndirella(
        self, out_key: 'str | Path', separate: bool = False
    ) -> 'DataFrame':
        """Create syndirella inputs"""

        out_key = Path('.') / out_key

        out_dir = out_key.parent
        out_key = out_key.name

        mrich.var('out_key', out_key)
        mrich.var('#poses', len(self))

        out_dir.mkdir(parents=True, exist_ok=True)

        ### Prepare Syndirella CSV data

        df = self.get_df(
            inchikey=False, alias=False, reference_alias=True, inspiration_aliases=True
        )
        df = df.rename(columns={'reference_alias': 'template'})

        # compound_set

        if separate:
            df['compound_set'] = df.apply(
                lambda row: f'{out_key}_{row["name"]}', axis=1
            )

        else:
            df['compound_set'] = out_key

        # template

        null_template = df['template'].isnull()
        if null_template.any():
            mrich.warning(
                len(null_template), 'poses have no reference. Setting to self'
            )
            mrich.print(df.loc[null_template, 'name'].values)
            df['template'] = df['template'].fillna(df['name'])

        # inspirations

        null_inspirations = df['inspiration_aliases'].apply(lambda x: not x)

        if null_inspirations.any():
            mrich.warning(
                len(null_inspirations), 'poses have no inspirations. Setting to self'
            )
            mrich.print(df.loc[null_inspirations, 'name'].values)
            df.loc[null_inspirations, 'inspiration_aliases'] = df.loc[
                null_inspirations
            ].apply(lambda row: set([row['name']]), axis=1)

        for i, row in df.iterrows():
            for j, inspiration in enumerate(row['inspiration_aliases']):
                df.loc[i, f'hit{j + 1}'] = inspiration

        # this from original code. looking at the data type I have,
        # this cannot possibly work. did I get something wrong filling
        # the df?
        # all_inspirations = set.union(*list(df['inspiration_aliases'].values))
        all_inspirations = set().union(*df['inspiration_aliases'])

        df = df.drop(columns=['name', 'inspiration_aliases'])

        ### Copy Templates

        template_dir = out_dir / 'templates'
        mrich.writing(template_dir)
        template_dir.mkdir(parents=True, exist_ok=True)

        templates = df['template'].unique()

        # records = self._queryset.filter(pose_alias__in=templates)
        records = PoseModel.objects.filter(
            target__in=self.targets,
            pose_alias__in=templates,
        )

        templates = PoseSet(records)

        for ref in templates:
            template = template_dir / ref.apo_path.name
            if not template.exists():
                mrich.writing(template)
                shutil.copy(ref.apo_path, template)

        ### Inspirations
        # isn't this overwriting the one few lines above??
        records = PoseModel.objects.filter(
            target__in=self.targets, pose_alias__in=all_inspirations
        )

        all_inspirations = PoseSet(records)

        ### Write CSV

        if separate:
            for _, row in df.iterrows():
                csv_name = out_dir / f'{row["compound_set"]}_syndirella_input.csv'
                mrich.writing(csv_name)
                row.to_frame().T.to_csv(csv_name, index=False)

        else:
            csv_name = out_dir / f'{out_key}_syndirella_input.csv'
            mrich.writing(csv_name)
            df.to_csv(csv_name, index=False)

        ### Write Inspirations

        sdf_name = out_dir / f'{out_key}_syndirella_inspiration_hits.sdf'
        all_inspirations.write_sdf(
            sdf_name,
            tags=False,
            metadata=False,
            name_col='name',
        )

        return df

    ### OUTPUT

    def interactive(
        self,
        print_name: str = True,
        method: str | None = None,
        function: Callable | None = None,
        **kwargs,
    ):
        """Interactive widget to navigate compounds in the table

        :param print_name: print the :class:`.PoseModel` name  (Default value = True)
        :param method: pass the name of a :class:`.PoseModel` method to interactively
            display. Keyword arguments to interactive() will be passed through
            (Default value = None)
        :param function: pass a callable which will be called as `function(pose)`

        """

        if method:

            def widget(i):
                """Method widget"""
                pose = self[i]
                if print_name:
                    print(repr(pose))
                value = getattr(pose, method)(**kwargs)
                if value:
                    display(value)

            return interactive(
                widget,
                i=BoundedIntText(
                    value=0,
                    min=0,
                    max=len(self) - 1,
                    step=1,
                    description='PoseModel:',
                    disabled=False,
                ),
            )

        elif function:

            def widget(i):
                """Function widget"""
                pose = self[i]
                if print_name:
                    display(pose)
                function(pose)

            return interactive(
                widget,
                i=BoundedIntText(
                    value=0,
                    min=0,
                    max=len(self) - 1,
                    step=1,
                    description='PoseModel:',
                    disabled=False,
                ),
            )

        else:
            a = BoundedIntText(
                value=0,
                min=0,
                max=len(self) - 1,
                step=1,
                description=f'PoseModel (/{len(self)}):',
                disabled=False,
            )

            b = Checkbox(description='Name', value=True)
            c = Checkbox(description='Summary', value=False)
            h = Checkbox(description='Tags', value=False)
            i = Checkbox(description='Subsites', value=False)
            d = Checkbox(description='2D (Comp.)', value=False)
            e = Checkbox(description='2D (PoseModel)', value=False)
            f = Checkbox(description='3D', value=True)
            g = Checkbox(description='Metadata', value=False)

            ui1 = GridBox(
                [b, c, d, h],
                layout=Layout(grid_template_columns='repeat(4, 100px)'),
            )
            ui2 = GridBox(
                [e, f, g, i],
                layout=Layout(grid_template_columns='repeat(4, 100px)'),
            )
            ui = VBox([a, ui1, ui2])

            def widget(
                i,
                name: bool = True,
                summary: bool = True,
                grid: bool = True,
                draw2d: bool = True,
                draw: bool = True,
                tags: bool = True,
                subsites: bool = True,
                metadata: bool = True,
            ):
                """Default widget"""
                pose = self._queryset.get(pk=i)
                if name:
                    print(repr(pose))

                if summary:
                    pose.summary(metadata=False, tags=False, subsites=False)
                if tags:
                    print(pose.tags)
                if subsites:
                    print(pose.subsites)
                if grid:
                    pose.grid()
                if draw2d:
                    pose.draw2d()
                if draw:
                    pose.draw()
                if metadata:
                    mrich.title('Metadata:')
                    pprint(pose.metadata)

            out = interactive_output(
                widget,
                {
                    'i': a,
                    'name': b,
                    'summary': c,
                    'grid': d,
                    'draw2d': e,
                    'draw': f,
                    'metadata': g,
                    'tags': h,
                    'subsites': i,
                },
            )

            display(ui, out)

    def summary(self) -> None:
        """Print a summary of this pose set"""
        mrich.header('PoseSet()')
        mrich.var('#poses', len(self))
        mrich.var('#compounds', self.num_compounds)
        mrich.var('tags', self.tags)

    def draw(self) -> None:
        """Render this pose set with Py3Dmol"""

        mols = [p.mol for p in self]

        drawing = draw_mols(mols)
        display(drawing)

    def grid(self) -> None:
        """Draw a grid of all contained molecules"""

        data = [(p.name, p.compound.mol) for p in self]

        mols = [d[1] for d in data]
        labels = [d[0] for d in data]

        drawing = draw_grid(mols, labels=labels)
        display(drawing)

    def subsite_summary(self) -> 'pd.DataFrame':
        """Print a table counting poses by subsite"""

        rows = (
            SubsiteTagModel.objects.filter(pose__in=self._queryset)
            .values('subsite_id', 'subsite__subsite_name')
            .annotate(num_poses=Count('pose', distinct=True))
        )

        df = DataFrame(
            [
                dict(
                    id=r['subsite_id'],
                    subsite=r['subsite__subsite_name'],
                    num_poses=r['num_poses'],
                )
                for r in rows
            ]
        )

        if len(df):
            df = df.set_index('id').sort_values(by='num_poses', ascending=False)

        mrich.print(df)

        return df

    def get_interaction_overlaps(self, return_pairs: bool = False) -> int:
        """Count the number of member pose pairs which share at least one but not all
        interactions"""

        records = InteractionModel.objects.filter(
            pose__in=self._queryset,
        ).values(
            'pose',
            'feature',
            'interaction_type',
        )

        ISETS = {}
        for r in records:
            pose_id = r['pose']
            feature_id = r['feature']
            interaction_type = r['interaction_type']
            values = ISETS.get(pose_id, set())
            values.add((interaction_type, feature_id))
            ISETS[pose_id] = values

        ids = [i for i in self.ids if i in ISETS]

        count = 0

        pairs = set()

        for pose_j, pose_k in combinations(ids, 2):
            iset_j = ISETS[pose_j]
            iset_k = ISETS[pose_k]

            intersection = iset_j & iset_k
            diff1 = iset_j - iset_k
            diff2 = iset_k - iset_j

            if intersection and diff1 and diff2:
                count += 1
                pairs.add((pose_j, pose_k))

        if return_pairs:
            return [PoseSet(PoseModel.objects.filter(pk__in=[a, b])) for a, b in pairs]

        return count

    def get_interaction_clusters(self) -> 'dict[int, PoseSet]':
        """Cluster poses based on shared interactions."""

        # get interaction records
        records = InteractionModel.objects.filter(
            pose__in=self._queryset,
        ).values(
            'pose',
            'feature__feature_residue_name',
            'feature__feature_residue_number',
            'interaction_type',
        )

        ISETS = {}
        for r in records:
            pose_id = r['pose']
            feature_residue_name = r['feature__feature_residue_name']
            feature_residue_number = r['feature__feature_residue_number']
            interaction_type = r['interaction_type']
            values = ISETS.get(pose_id, set())
            values.add((interaction_type, feature_residue_name, feature_residue_number))
            ISETS[pose_id] = values

        pairs = combinations(ISETS.keys(), 2)

        # construct overlap dictionary

        OVERLAPS = {}
        for id1, id2 in pairs:
            iset1 = ISETS[id1]
            iset2 = ISETS[id2]
            OVERLAPS[(id1, id2)] = len(iset1 & iset2)

        # make the graph
        G = nx.Graph()

        for (id1, id2), count in OVERLAPS.items():
            G.add_edge(id1, id2, weight=count)

        # partition the graph

        partition = louvain.best_partition(G, weight='weight')

        # find the clusters

        clusters = {}
        for node, cluster_id in partition.items():
            clusters.setdefault(cluster_id, set()).add(node)

        # create the PoseSets

        psets = {
            i: PoseSet(PoseModel.objects.filter(pk__in=ids), name=f'Cluster {i}')
            for i, ids in enumerate(clusters.values())
        }

        all_ids = set(sum((pset.ids for pset in psets.values()), []))

        # calculate modal interactions

        for cluster in psets.values():
            mrich.var(cluster.name, len(cluster), unit='poses')

            df = cluster.interactions.df

            unique_counts = df.groupby(['type', 'residue_name', 'residue_number'])[
                'pose_id'
            ].nunique()

            max_count = unique_counts.max()
            max_pairs = unique_counts[unique_counts == max_count]

            for (
                interaction_type,
                residue_name,
                residue_number,
            ) in max_pairs.index.values:
                mrich.print(interaction_type, 'w/', residue_name, residue_number)

        # unclustered
        unclustered = set(i for i in self.ids if i not in all_ids)
        psets[None] = PoseSet(
            PoseModel.objects.filter(pk__in=unclustered), name='Unclustered'
        )

        return psets

    ### PROPERTIES

    @property
    def queryset(self) -> QuerySet[PoseModel]:
        """Returns the ids of poses in this set"""
        return self._queryset

    @property
    def indices(self) -> list[int]:
        """Returns the ids of poses in this set"""
        return self.queryset.values_list('id', flat=True)

    @property
    def ids(self) -> list[int]:
        """Returns the ids of poses in this set"""
        return self.indices

    @property
    def name(self) -> str | None:
        """Returns the name of set"""
        return self._name

    @property
    def names(self) -> list[str]:
        """Returns the aliases of poses in this set"""
        return self._queryset.values_list('pose_alias', flat=True)

    @property
    def aliases(self) -> list[str]:
        """Returns the aliases of child poses"""
        return self._queryset.values_list('pose_alias', flat=True)

    @property
    def inchikeys(self) -> list[str]:
        """Returns the inchikeys of child poses"""
        return self._queryset.values_list('pose_inchikey', flat=True)

    @property
    def id_name_dict(self) -> dict:
        """Return a dictionary mapping pose ID's to their name"""
        return dict(self._queryset.values_list('pk', 'pose_alias'))

    @property
    def smiles(self) -> list[str]:
        """Returns the smiles of poses in this set"""
        return self._queryset.values_list('pose_smiles', flat=True)

    @property
    def tags(self) -> set[str]:
        """Returns the set of unique tags present in this pose set"""
        return self._queryset.values_list('tags__pose_tag_name', flat=True).distinct()

    @property
    def num_fingerprinted(self) -> int:
        """Count the number of fingerprinted poses"""
        # that's one field suspect not in use
        return self._queryset.filter(pose_fingerprint=1).count()

    @property
    def compounds(self) -> 'CompoundSet':
        """Get the compounds associated to this set of poses"""
        # local import to avoid the PoseSet <-> CompoundSet cycle
        from designdb.sets.compound import CompoundSet

        ids = list(self._queryset.values_list('compound_id', flat=True).distinct())
        return CompoundSet(ids)

    @property
    def mols(self) -> list[Chem.rdchem.Mol]:
        """Get the rdkit Molecules contained in this set"""
        return self._queryset.values_list('pose_mol', flat=True)

    @property
    def num_compounds(self) -> int:
        """Count the compounds associated to this set of poses"""
        return self._queryset.values('compound').distinct().count()

    @property
    def df(self) -> pd.DataFrame:
        """Get a DataFrame of the poses in this set"""
        return self.get_df(mol=True)

    @property
    def references(self) -> 'PoseSet':
        """Return a :class:`.PoseSet` of the all the distinct references in this
        :class:`.PoseSet`"""
        # TODO: call through proper factory method
        return self.get_by_references(self)

    @property
    def reference_ids(self) -> set[int]:
        """Return a set of :class:`.PoseModel` ID's of the all the distinct references
        in this :class:`.PoseSet`"""
        return self.get_by_references(self).values_list('pk', flat=True)

    @property
    def inspiration_sets(self) -> set[tuple[int, ...]]:
        """Return the unique sets of inspiration :class:`.PoseModel` IDs"""

        # group original (inspiration) pose IDs by derivative pose ID -- use the
        # FK ids (sortable ints, no extra queries) rather than the model objects
        data: dict[int, set[int]] = {}
        for deriv_id, orig_id in InspirationModel.objects.filter(
            derivative_pose__in=self._queryset
        ).values_list('derivative_pose_id', 'original_pose_id'):
            data.setdefault(deriv_id, set()).add(orig_id)

        return {tuple(sorted(v)) for v in data.values()}

    @property
    def num_inspiration_sets(self) -> int:
        """Return the number of unique sets of inspirations"""
        return len(self.inspiration_sets)

    @property
    def num_inspirations(self) -> int:
        """Return the number of unique inspirations for poses in this set"""
        # fmt: off
        return InspirationModel.objects.filter(
            derivative_pose__in=self._queryset,
        ).values(
            'original_pose',
        ).distinct().count()
        # fmt: on

    @property
    def inspirations(self) -> int:
        """Return the number of unique inspirations for poses in this set"""
        return self.get_by_inspirations(self._queryset)

    # @property
    # def str_ids(self) -> str:
    #     """Return an SQL formatted tuple string of the :class:`.PoseModel` IDs"""
    #     return str(tuple(self.ids)).replace(',)', ')')

    @property
    def targets(self) -> QuerySet[TargetModel]:
        """Returns the :class:`.TargetModel` objects of poses in this set"""
        return TargetModel.objects.filter(pk__in=self._queryset.values('target'))

    @property
    def target_names(self) -> list[str]:
        """Returns the :class:`.TargetModel` objects of poses in this set"""
        return self.targets.values_list('target_name', flat=True)

    @property
    def target_ids(self) -> list[int]:
        """Returns the :class:`.TargetModel` objects ID's of poses in this set"""
        return self.targets.values_list('id', flat=True)

    @property
    def best_placed_pose(self) -> PoseModel:
        """Returns the pose with the best distance_score in this subset"""
        return self._queryset.get(pk=self.best_placed_pose_id)

    @property
    def best_placed_pose_id(self) -> int:
        """Get the id of the pose with the best distance_score in this subset"""

        # if len(self) == 1:
        #     return self.ids[0]

        # query = 'pose_id, MIN(pose_distance_score)'
        # query = self.db.select_where(
        #     table='pose', query=query,
        #     key=f'pose_id in {self.str_ids}', multiple=False
        # )
        # return query[0]

        # TODO: scoring not implemented yet
        return self.queryset.first().pk

    @property
    def interactions(self) -> 'InteractionSet':
        """Get a :class:`.InteractionSet` for this :class:`.PoseModel`"""
        if self._interactions is None:
            self._interactions = InteractionSet.from_pose(self)
        return self._interactions

    @property
    def pose_id_metadata_dict(self) -> dict[int, dict]:
        """Get a dictionary mapping pose_ids to metadata dicts"""
        if self._metadata_dict is None:
            metadata = {}
            for p in self._queryset:
                metadata[p.pk] = p.pose_metadata
            self._metadata_dict = metadata
        return self._metadata_dict

    @property
    def fraction_fingerprinted(self) -> float:
        """Return the fraction of fingerprinted poses in this set"""
        return self.num_fingerprinted / len(self)

    @property
    def num_subsites(self) -> int:
        """Count the number of subsites that poses in this set come into contact with"""
        return (
            SubsiteModel.objects.filter(posemodels__in=self._queryset)
            .distinct()
            .count()
        )

    @property
    def subsite_balance(self) -> float:
        """Measure of how evenly subsite counts are distributed across poses in this
        set"""
        # TODO: subsites not implemented yet
        # from numpy import std

        # sql = f"""
        # SELECT COUNT(DISTINCT subsite_tag_ref)
        # FROM {self.db.SQL_SCHEMA_PREFIX}subsite_tag
        # WHERE subsite_tag_pose IN {self.str_ids}
        # GROUP BY subsite_tag_pose
        # """

        # counts = self.db.execute(sql).fetchall()

        # counts = [c for (c,) in counts] + [0 for _ in range(len(self) - len(counts))]

        # return -std(counts)
        return 4

    @property
    def subsite_ids(self) -> set[int]:
        """Return a list of subsite id's of member poses"""
        return SubsiteModel.objects.filter(
            pk__in=SubsiteTagModel.objects.filter(
                pose__in=self._queryset,
            ).values('subsite'),
        ).values_list('pk', flat=True)

    @property
    def avg_energy_score(self) -> float:
        """Average energy score of poses in this set"""
        # TODO: scores not implemented
        # from numpy import mean

        # sql = f"""
        # SELECT pose_energy_score
        # FROM {self.db.SQL_SCHEMA_PREFIX}pose
        # WHERE pose_id IN {self.str_ids}
        # """

        # scores = self.db.execute(sql).fetchall()
        # return mean([s for (s,) in scores if s is not None])
        return 4

    @property
    def avg_distance_score(self) -> float:
        """Average distance score of poses in this set"""
        # TODO: scores not implemented yet
        # from numpy import mean

        # sql = f"""
        # SELECT pose_distance_score
        # FROM {self.db.SQL_SCHEMA_PREFIX}pose
        # WHERE pose_id IN {self.str_ids}
        # """

        # scores = self.db.execute(sql).fetchall()

        # return mean([s for (s,) in scores if s is not None])
        return 4

    @property
    def derivatives(self) -> 'PoseSet':
        """Get the :class:`.PoseSet` of derivatives"""
        return PoseSet(
            PoseModel.objects.filter(
                pk__in=InspirationModel.objects.filter(
                    original_pose__in=self._queryset,
                ).values(
                    'derivative_pose',
                ),
            ),
        )

    @property
    def reference(self):
        """Bulk set the references for poses in this set"""
        raise NotImplementedError(
            'This attribute only allows setting, ``PoseSet.reference = ...``'
        )

    @reference.setter
    def reference(self, r) -> None:
        """Bulk set the references for poses in this set"""
        self._queryset.update(pose_reference=r)

    ### PRIVATE

    def _delete(self, *, force: bool = False) -> None:
        """Delete poses in this set"""

        if not force:
            mrich.warning('Deleting Poses is risky! Set force=True to continue')
            return

        try:
            with transaction.atomic():
                InspirationModel.objects.filter(
                    original_pose__in=self._queryset
                ).delete()
                InspirationModel.objects.filter(
                    derivative_pose__in=self._queryset
                ).delete()
                SubsiteTagModel.objects.filter(pose__in=self._queryset).delete()
                InteractionModel.objects.filter(pose__in=self._queryset).delete()
                self._queryset.delete()
        except IntegrityError as exc:
            mrich.error(exc)
