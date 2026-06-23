import json
from typing import TYPE_CHECKING

import mcol
import mrich
from designdb.models import ComponentModel, RouteModel
from designdb.sets.compound import CompoundSet

if TYPE_CHECKING:
    from pathlib import Path

    from designdb.recipe import Route


class RouteSet:
    """A set of Route objects"""

    def __init__(self, routes: 'list[Route]') -> None:
        """RouteSet initialisation"""

        data = {}
        for route in routes:
            # assert isinstance(route, Route)
            data[route.id] = route

        self._data = data
        self._cluster_map = None
        self._permitted_clusters = None
        self._current_cluster = None

    ### FACTORIES

    @classmethod
    def from_ids(cls, ids: list | set, progress: bool = True):
        """Generate a routeset from a set of :class:`.RouteModel` IDs

        :param db: database to link
        :param ids: :class:`.RouteModel` database IDs
        :param progress: show progress bar
        """

        # this gets stuck
        # if progress:
        #     ids = mrich.track(ids, prefix='Getting routes')

        # avoiding circular reference
        # avoiding name conflict
        from designdb.recipe import Route

        routes = [Route.get_route(id=r) for r in ids]

        # self = cls.__new__(cls)
        return RouteSet(routes)

    @classmethod
    def from_product_ids(cls, ids: list | set, progress: bool = True):
        """Generate a routeset from a set of product :class:`.CompoundModel` IDs

        :param db: database to link
        :param ids: :class:`.CompoundModel` database IDs
        """

        # str_ids = str(tuple(ids)).replace(',)', ')')

        # records = db.select_where(
        #     table='route',
        #     query='route_id',
        #     key=f'route_product IN {str_ids}',
        #     multiple=True,
        # )
        records = RouteModel.objects.filter(
            product_compound__pk__in=ids,
        )

        # route_ids = [i for (i,) in records]

        return cls.from_ids(records.values_list('id', flat=True), progress=progress)

    @classmethod
    def from_json(cls, path: 'str | Path', data: dict = None) -> 'RouteSet':
        """Load a serialised routeset from a JSON file

        :param db: database to link
        :param path: path to JSON
        :param data: serialised data (Default value = None)

        """

        from designdb.recipe import Route

        self = cls.__new__(cls)

        if data is None:
            data = json.load(open(path))

        new_data = {}
        for d in mrich.track(data['routes'].values(), prefix='Loading Routes...'):
            route_id = d['id']
            new_data[route_id] = Route.from_json(path=None, data=d)

        self._data = new_data
        self._cluster_map = None
        self._permitted_clusters = None
        self._current_cluster = None

        return self

    ### PROPERTIES

    @property
    def data(self) -> 'dict[int, RouteModel]':
        """Get internal data dictionary"""
        return self._data

    @property
    def db(self):
        """Deprecated: the modern ORM RouteSet has no associated ``db`` handle.

        Raises to surface any remaining legacy ``self.db`` callers explicitly.
        """
        raise NotImplementedError(
            'RouteSet no longer holds a `db` handle; use the Django ORM directly'
        )

    @property
    def routes(self) -> 'list[Route]':
        """Get route objects"""
        return self.data.values()

    @property
    def product_ids(self) -> list[int]:
        """Get the :class:`.CompoundModel` ID's of the products"""
        return RouteModel.objects.values_list(
            'product_compound__id', flat=True
        ).distinct()

    @property
    def reactant_ids(self) -> list[int]:
        """Get the :class:`.CompoundModel` ID's of the reactants"""

        return ComponentModel.objects.filter(
            route__in=self.ids,
            component_type=2,
        ).values_list('component_ref', flat=True)

    @property
    def products(self) -> 'CompoundSet':
        """Return a :class:`.CompoundSet` of all the route products"""
        return CompoundSet(self.product_ids)

    @property
    def reactants(self) -> 'CompoundSet':
        """Return a :class:`.CompoundSet` of all the route reactants"""
        return CompoundSet(self.reactant_ids)

    @property
    def str_ids(self) -> str:
        """Return an SQL formatted tuple string of the :class:`.RouteModel` ID's"""
        return str(tuple(self.ids)).replace(',)', ')')

    @property
    def ids(self) -> list[int]:
        """Return the :class:`.RouteModel` IDs"""
        return self.data.keys()

    @property
    def cluster_map(self) -> dict[tuple, set]:
        """Create a dictionary grouping routes by their scaffold/base cluster.

        :returns: A dictionary mapping a tuple of scaffold :class:`.CompoundModel` IDs
            to a set of :class:`.RouteModel` ID's to their superstructures.
        """

        # NOTE: scaffold-cluster grouping depends on the not-yet-ported
        # `get_compound_cluster_dict` machinery. This property is only consumed by
        # `balanced_pop`, which in turn is only used by the (out-of-scope) rgen
        # RandomRecipeSelectionGenerator. Port alongside rgen.
        raise NotImplementedError(
            'RouteSet.cluster_map requires the unported compound-clustering helper '
            '(get_compound_cluster_dict); port it together with the rgen subsystem'
        )

    ### METHODS

    def copy(self) -> 'RouteSet':
        """Copy this RouteSet"""
        return RouteSet(list(self.data.values()))

    def get_dict(self):
        """Get serialisable dictionary"""

        data = dict(routes={})

        # populate with routes
        for route_id, route in self.data.items():
            data['routes'][route_id] = route.get_dict()

        return data

    def prune_unavailable(self, suppliers: list[str]):
        """Remove routes that don't have all reactants available from given suppliers.

        Keeps only routes where every reactant component (``component_type == 2``)
        has at least one catalogue price from one of the given ``suppliers``.
        """

        from designdb.models import CataloguePriceCompoundJunctionModel

        # compound IDs quotable from the permitted suppliers
        quotable = set(
            CataloguePriceCompoundJunctionModel.objects.filter(
                catalogue_price__supplier__in=list(suppliers),
            ).values_list('compound_id', flat=True)
        )

        # reactant components grouped by route
        reactant_rows = ComponentModel.objects.filter(
            route_id__in=list(self.ids),
            component_type=2,
        ).values_list('route_id', 'component_ref')

        route_reactants: dict[int, set[int]] = {}
        for route_id, ref in reactant_rows:
            route_reactants.setdefault(route_id, set()).add(ref)

        kept = [
            route_id
            for route_id, reactants in route_reactants.items()
            if reactants <= quotable
        ]

        mrich.var('#routes before pruning', len(self))
        mrich.var('#routes after pruning', len(kept))

        return RouteSet.from_ids(kept)

    def pop_id(self) -> int:
        """Pop the last route from the set and return it's id"""
        route_id, route = self.data.popitem()
        return route_id

    def pop(self) -> 'RouteModel':
        """Pop the last route from the set and return it's object"""
        route_id, route = self.data.popitem()
        return route

    def balanced_pop(
        self, permitted_clusters: set[tuple] | None = None, debug: bool = False
    ) -> 'RouteModel':
        """Pop a route from this set, while maintaining the balance of scaffold
        clusters populations"""

        if not self._data:
            mrich.print('RouteSet depleted')
            return None

        if not self.cluster_map:
            # mrich.warning("RouteSet.cluster_map depleted but _data isn't...")
            return self.pop()

        # store the permitted clusters (or all clusters) list as property

        if self._permitted_clusters is None:
            if permitted_clusters:
                permitted_clusters = set(
                    (cluster,) if isinstance(cluster, int) else cluster
                    for cluster in permitted_clusters
                )

                self._permitted_clusters = []
                for cluster in permitted_clusters:
                    if cluster not in self.cluster_map:
                        mrich.warning(
                            cluster, 'in permitted_clusters but not cluster_map'
                        )
                    else:
                        self._permitted_clusters.append(cluster)

            else:
                self._permitted_clusters = list(self.cluster_map.keys())

        if self._current_cluster is None:
            self._current_cluster = self._permitted_clusters[0]

        ### pop a RouteModel

        if debug:
            mrich.debug(f'Would pop RouteModel from {self._current_cluster=}')

        cluster = self._current_cluster

        # pop the last route id from the given cluster

        try:
            route_id = self.cluster_map[cluster].pop()
        except IndexError:
            mrich.print(self._permitted_clusters)
            mrich.print(self.cluster_map)
            raise
        except AttributeError:
            mrich.print(cluster)
            mrich.print(self.cluster_map)
            raise
        except KeyError:
            mrich.print('cluster', cluster)
            mrich.print('self._permitted_clusters', self._permitted_clusters)
            mrich.print('self.cluster_map.keys()', self.cluster_map.keys())
            raise

        # clean up empty clusters

        if debug:
            mrich.debug('Popped route', route_id)

        # get the RouteModel object

        if route_id in self._data:
            route = self._data[route_id]
            del self._data[route_id]
        else:
            # if debug:
            mrich.debug('RouteModel not present')
            return self.balanced_pop()

        ### increment cluster

        # def increment_cluster(cluster):
        n = len(self._permitted_clusters)
        if n > 1:
            for i, cluster in enumerate(self._permitted_clusters):
                if cluster == self._current_cluster:
                    if i == n - 1:
                        self._current_cluster = self._permitted_clusters[0]
                    else:
                        self._current_cluster = self._permitted_clusters[i + 1]
                    break
            else:
                raise IndexError('This should never be reached...')

        # increment_cluster()

        if not self.cluster_map[cluster]:
            del self.cluster_map[cluster]
            if not self.cluster_map:
                mrich.debug('RouteSet.cluster_map depleted')
            self._permitted_clusters = [
                c for c in self._permitted_clusters if c != cluster
            ]
            # if debug:
            mrich.debug('Depleted cluster', cluster)

            if not self._permitted_clusters:
                mrich.debug('Depleted all permitted clusters', cluster)
                mrich.debug('Removing cluster restriction', cluster)
                self._permitted_clusters = list(self.cluster_map.keys())
                self._current_cluster = None

        if debug:
            mrich.debug('#Routes in set', len(self._data))

        return route

    def shuffle(self):
        """Randomly shuffle the routes in this set"""
        import random

        items = list(self.data.items())
        random.shuffle(items)
        self._data = dict(items)

        ### shuffle the cluster map as well

        for cluster, routes in self.cluster_map.items():
            random.shuffle(routes)
            self.cluster_map[cluster] = routes

    ### DUNDERS

    def __len__(self) -> int:
        """Number of routes in this set"""
        return len(self.data)

    def __str__(self) -> str:
        """Unformatted string representation"""
        return f'{{RouteModel × {len(self)}}}'

    def __repr__(self) -> str:
        """ANSI Formatted string representation"""
        return f'{mcol.bold}{mcol.underline}{self}{mcol.unbold}{mcol.ununderline}'

    def __rich__(self) -> str:
        """Rich Formatted string representation"""
        return f'[bold underline]{self}'

    def __iter__(self):
        """Iterate over routes in this set"""
        return iter(self.data.values())

    def __getitem__(self, key):
        """Get a specific route in this set"""
        return list(self.data.values())[key]
