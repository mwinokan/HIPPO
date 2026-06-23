"""Classes for working with Recipes (reaction networks).

A :class:`.Recipe` is a lean *aggregate*: it holds the products, reactants,
intermediates, reactions and (no-chem) compounds that make up a synthetic recipe,
and exposes price/serialisation/presentation on top of them.

Construction and DB-traversal *orchestration* lives in :class:`.RecipeService`.
The ``from_*`` and export methods on :class:`.Recipe` are **deprecated shims** that
delegate to it (see the ``DEPRECATED`` banner below) via a local import.
"""

from typing import TYPE_CHECKING

import mcol
import mrich
from designdb.components.compound import Ingredient
from designdb.components.reaction import Reaction
from designdb.models import ComponentModel, CompoundModel, ReactionModel, RouteModel
from designdb.sets.compound import CompoundSet
from designdb.sets.ingredient import IngredientSet
from designdb.sets.reaction import ReactionSet

if TYPE_CHECKING:
    from pathlib import Path

    import pandas
    from designdb.components.price import Price
    from designdb.sets.interaction import InteractionSet
    from designdb.sets.pose import PoseSet
    from designdb.sets.route import RouteSet
    from plotly import graph_objects


class Recipe:
    """A Recipe stores data corresponding to a specific synthetic recipe involving
    several products, reactants, intermediates, and reactions."""

    def __init__(
        self,
        *,
        products: 'IngredientSet | None' = None,
        reactants: 'IngredientSet | None' = None,
        intermediates: 'IngredientSet | None' = None,
        reactions: 'ReactionSet | None' = None,
        compounds: 'IngredientSet | None' = None,
    ) -> None:
        """Recipe initialisation"""

        if products is None:
            products = IngredientSet()
        if reactants is None:
            reactants = IngredientSet()
        if intermediates is None:
            intermediates = IngredientSet()
        if compounds is None:
            compounds = IngredientSet()
        if reactions is None:
            reactions = ReactionSet()

        # check typing
        assert isinstance(products, IngredientSet)
        assert isinstance(reactants, IngredientSet)
        assert isinstance(intermediates, IngredientSet)
        assert isinstance(compounds, IngredientSet)
        assert isinstance(reactions, ReactionSet)

        self._products = products
        self._reactants = reactants
        self._intermediates = intermediates
        self._reactions = reactions
        self._compounds = compounds
        self._hash = None

        self._score = None

        # caches
        self._product_compounds = None
        self._poses = None
        self._interactions = None
        self._combined_compounds = None

    ### DEPRECATED — construction shims (relocated to RecipeService)
    # These delegate to designdb.services.recipe.RecipeService and exist only to
    # keep legacy `Recipe.from_*(...)` call sites working during the migration.

    @classmethod
    def from_reaction(cls, *args, **kwargs):
        """DEPRECATED: use :meth:`.RecipeService.from_reaction`."""
        from designdb.services.recipe import RecipeService

        return RecipeService.from_reaction(*args, **kwargs)

    @classmethod
    def from_reactions(cls, *args, **kwargs):
        """DEPRECATED: use :meth:`.RecipeService.from_reactions`."""
        from designdb.services.recipe import RecipeService

        return RecipeService.from_reactions(*args, **kwargs)

    @classmethod
    def from_compounds(cls, *args, **kwargs):
        """DEPRECATED: use :meth:`.RecipeService.from_compounds`."""
        from designdb.services.recipe import RecipeService

        return RecipeService.from_compounds(*args, **kwargs)

    @classmethod
    def from_reactants(cls, *args, **kwargs):
        """DEPRECATED: use :meth:`.RecipeService.from_reactants`."""
        from designdb.services.recipe import RecipeService

        return RecipeService.from_reactants(*args, **kwargs)

    ### FACTORIES

    @classmethod
    def from_json(
        cls,
        path: 'str | Path | None' = None,
        *,
        data: dict | None = None,
        clear_quotes: bool = False,
        debug: bool = False,
    ) -> 'Recipe':
        """Load a serialised recipe from a JSON file (see :meth:`.Recipe.get_dict`).

        :param path: path to JSON (ignored if ``data`` is provided)
        :param data: pre-loaded serialised data (Default value = None)
        :param clear_quotes: ignore stored reactant/compound quotes
        :param debug: increase verbosity
        """

        import json

        if not data:
            if debug:
                mrich.reading(path)
            data = json.load(open(path))

        if debug and 'timestamp' in data:
            mrich.print(f'Recipe was generated at: {data["timestamp"]}')

        # IngredientSets are stored column-oriented (df.to_dict(orient='list'))
        products = IngredientSet.from_json(path=None, data=data['products'])
        intermediates = IngredientSet.from_json(path=None, data=data['intermediates'])
        reactants = IngredientSet.from_json(
            path=None, data=data['reactants'], supplier=data.get('reactant_supplier')
        )

        if 'compounds' in data:
            compounds = IngredientSet.from_json(
                path=None,
                data=data['compounds'],
                supplier=data.get('compound_supplier'),
            )
        else:
            compounds = IngredientSet()

        if clear_quotes:
            for iset in (reactants, compounds):
                iset.df['quote_id'] = None
                iset.df['quoted_amount'] = None

        reactions = ReactionSet(data['reaction_ids'], sort=False)

        if debug:
            mrich.var('reactants', reactants)
            mrich.var('intermediates', intermediates)
            mrich.var('products', products)
            mrich.var('reactions', reactions)
            mrich.var('compounds', compounds)

        return cls(
            products=products,
            reactants=reactants,
            intermediates=intermediates,
            reactions=reactions,
            compounds=compounds,
        )

    ### PROPERTIES

    @property
    def products(self) -> 'IngredientSet':
        """Product :class:`.IngredientSet`"""
        return self._products

    @products.setter
    def products(self, a: 'IngredientSet'):
        """Set the products"""
        self._products = a
        self.__flag_modification()

    @property
    def compounds(self) -> 'IngredientSet':
        """No-chem (directly purchased) :class:`.IngredientSet`"""
        return self._compounds

    @compounds.setter
    def compounds(self, a: 'IngredientSet'):
        """Set the compounds"""
        self._compounds = a
        self.__flag_modification()

    @property
    def reactants(self) -> 'IngredientSet':
        """Reactant :class:`.IngredientSet`"""
        return self._reactants

    @reactants.setter
    def reactants(self, a: 'IngredientSet'):
        """Set the reactants"""
        self._reactants = a
        self.__flag_modification()

    @property
    def intermediates(self) -> 'IngredientSet':
        """Intermediate :class:`.IngredientSet`"""
        return self._intermediates

    @intermediates.setter
    def intermediates(self, a: 'IngredientSet'):
        """Set the intermediates"""
        self._intermediates = a

    @property
    def reactions(self) -> 'ReactionSet':
        """:class:`.ReactionSet` for this recipe"""
        return self._reactions

    @reactions.setter
    def reactions(self, a: 'ReactionSet'):
        """Set the reactions"""
        self._reactions = a
        self.__flag_modification()

    @property
    def product(self) -> 'Ingredient':
        """Return the single product (if there's only one)"""
        assert len(self.products) == 1
        return self.products[0]

    @property
    def product_compounds(self) -> 'CompoundSet':
        """Product compounds"""
        if self._product_compounds is None:
            self._product_compounds = self.products.compounds
            self._product_compounds._name = f'products of {self}'
        return self._product_compounds

    @property
    def combined_compound_ids(self) -> set[int]:
        """Combined :class:`.CompoundModel` IDs from :meth:`.Recipe.product_compounds`
        and :meth:`.Recipe.compounds`"""
        return set(self.product_compounds.ids) | set(self.compounds.ids)

    @property
    def combined_compounds(self) -> 'CompoundSet':
        """Combined product and no-chem compounds"""
        if self._combined_compounds is None:
            self._combined_compounds = CompoundSet(list(self.combined_compound_ids))
            self._combined_compounds._name = f'combined compounds of {self}'
        return self._combined_compounds

    @property
    def poses(self) -> 'PoseSet':
        """Poses of the combined compounds"""
        if self._poses is None:
            self._poses = self.combined_compounds.poses
            self._poses._name = f'poses of {self}'
        return self._poses

    @property
    def interactions(self) -> 'InteractionSet':
        """Product pose interactions"""
        if self._interactions is None:
            self._interactions = self.poses.interactions
        return self._interactions

    @property
    def price(self) -> 'Price':
        """Total price of the reactants and no-chem compounds"""
        return self.reactants.get_price() + self.compounds.get_price()

    @property
    def num_products(self) -> int:
        """Number of products"""
        return len(self.products)

    @property
    def num_compounds(self) -> int:
        """Number of combined compounds"""
        return len(self.combined_compound_ids)

    @property
    def num_reactions(self) -> int:
        """Number of reactions"""
        return len(self.reactions)

    @property
    def num_reaction_types(self) -> int:
        """Number of distinct reaction types"""
        return self.reactions.num_types

    @property
    def num_reactants(self) -> int:
        """Number of reactants"""
        return len(self.reactants)

    @property
    def num_intermediates(self) -> int:
        """Number of intermediates"""
        return len(self.intermediates)

    @property
    def hash(self) -> str:
        """Unique hash string (set when loaded from a RecipeSet)"""
        return self._hash

    @property
    def score(self):
        """Recipe score"""
        return self._score

    @property
    def type(self) -> str:
        """Recipe type (EMPTY/MIXED/CHEM/NOCHEM)"""

        if self.empty:
            return 'EMPTY'

        chem = bool(self.reactions)
        nochem = bool(self.compounds)

        if chem and nochem:
            return 'MIXED'
        if chem and not nochem:
            return 'CHEM'
        if nochem and not chem:
            return 'NOCHEM'

    @property
    def empty(self) -> bool:
        """Is this Recipe empty?"""
        return not any(
            (
                self.reactants,
                self.products,
                self.intermediates,
                self.reactions,
                self.compounds,
            )
        )

    ### METHODS

    def get_price(self, supplier: str | None = None) -> 'Price':
        """Get the reactants price. See :meth:`.IngredientSet.get_price`

        :param supplier: restrict quotes to this supplier
        """
        return self.reactants.get_price(supplier=supplier)

    def get_ingredient(self, id) -> 'Ingredient':
        """Get an ingredient by its compound ID

        :param id: compound ID
        """
        matches = [r for r in self.reactants if r.id == id]
        if not matches:
            matches = [r for r in self.intermediates if r.id == id]
        if not matches:
            matches = [r for r in self.products if r.id == id]

        assert len(matches) == 1
        return matches[0]

    def add_ingredient(self, ingredient: 'Ingredient', amount: float = 1):
        """Add an :class:`.Ingredient` for direct purchase (no associated reactions)"""
        self.compounds.add(ingredient)

    def add_to_all_reactants(self, amount: float = 20) -> None:
        """Increment all reactants by this amount

        :param amount: amount in ``mg`` (Default value = 20)
        """
        self.reactants.df['amount'] += amount

    def copy(self) -> 'Recipe':
        """Copy this recipe"""
        return Recipe(
            products=self.products.copy(),
            reactants=self.reactants.copy(),
            intermediates=self.intermediates.copy(),
            reactions=self.reactions.copy(),
            compounds=self.compounds.copy(),
        )

    def check_integrity(self, debug: bool = False) -> bool:
        """Verify the internal integrity of this recipe."""

        if debug:
            mrich.debug('Checking integrity:', self)
            mrich.debug('Checking for duplicate compounds')

        for label, iset in (
            ('Reactant', self.reactants),
            ('Intermediate', self.intermediates),
            ('Product', self.products),
        ):
            if len(iset.compound_ids) != len(set(iset.compound_ids)):
                mrich.error(f"{label} compound ID's are not unique")
                return False

        if debug:
            mrich.debug('Checking for missing references')

        # all references should exist in the database
        if ReactionModel.objects.filter(pk__in=self.reactions.ids).count() < len(
            self.reactions
        ):
            mrich.error('Not all Reactions in Database')
            return False

        checks = (
            ('product', self.product_compounds.ids, len(self.products)),
            ('reactant', self.reactants.compounds.ids, len(self.reactants)),
            ('intermediate', self.intermediates.compounds.ids, len(self.intermediates)),
        )
        for label, ids, expected in checks:
            if CompoundModel.objects.filter(pk__in=list(ids)).count() < expected:
                mrich.error(f'Not all {label} Compounds in Database')
                return False

        reaction_intermediates = self.reactions.intermediates
        reaction_products = self.reactions.products
        reaction_reactants = self.reactions.reactants

        if debug:
            mrich.debug('Checking for missing reactions')

        for product in self.products:
            if product not in reaction_products:
                mrich.error(f'Product: {product} does not have associated reaction')
                return False

        for intermediate in self.intermediates:
            if intermediate not in reaction_intermediates:
                mrich.error(
                    f'Intermediate: {intermediate} is not in '
                    f'self.reactions.intermediates'
                )
                return False

        for reactant in self.reactants:
            if reactant not in reaction_reactants:
                mrich.error(f'Reactant: {reactant} is not in self.reactions.reactants')
                return False

        if debug:
            mrich.debug('Checking reactant quantities')

        for reaction in (Reaction(r) for r in self.reactions):
            product_ingredient = self.products(compound_id=reaction.product.id)
            if product_ingredient is None:
                product_ingredient = self.intermediates(compound_id=reaction.product.id)

            if debug and reaction.product_yield < 1.0:
                mrich.debug(f'{reaction}.product_yield={reaction.product_yield}')

            for reactant in reaction.reactants:
                reactant_ingredient = self.intermediates(compound_id=reactant.id)
                if reactant_ingredient is None:
                    reactant_ingredient = self.reactants(compound_id=reactant.id)

                required_amount = product_ingredient.amount / reaction.product_yield

                if reactant_ingredient.amount < required_amount:
                    mrich.error(
                        f'Not enough of {reactant_ingredient.compound}: '
                        f'{reactant_ingredient.amount} < {required_amount}'
                    )
                    return False

        if debug:
            mrich.success(self, 'OK')

        return True

    ### SERIALISATION

    def get_dict(
        self,
        *,
        price: bool = True,
        reactant_supplier: bool = True,
        compound_supplier: bool = True,
        timestamp: bool = True,
        compound_ids_only: bool = False,
        products: bool = True,
        serialise_price: bool = False,
    ) -> dict:
        """Serialise this recipe to a dictionary.

        :param price: include the price
        :param reactant_supplier: include the reactant supplier
        :param compound_supplier: include the compound supplier
        :param timestamp: add a timestamp
        :param compound_ids_only: store IDs only (instead of full ingredient dataframes)
        :param products: include products
        :param serialise_price: serialise the :class:`.Price` object
        """

        from datetime import datetime

        data = {}

        if timestamp:
            data['timestamp'] = str(datetime.now())

        try:
            if price and serialise_price:
                data['price'] = self.price.get_dict()
            elif price:
                data['price'] = self.price
        except AssertionError as e:
            mrich.warning(f'Could not get price: {e}')
            data['price'] = None

        if reactant_supplier:
            data['reactant_supplier'] = self.reactants.supplier
        if compound_supplier:
            data['compound_supplier'] = self.compounds.supplier

        if compound_ids_only:
            data['reactant_ids'] = self.reactants.compound_ids
            data['intermediate_ids'] = self.intermediates.compound_ids
            if products:
                data['products_ids'] = self.products.compound_ids
            data['compound_ids'] = self.compounds.compound_ids
        else:
            data['reactants'] = self.reactants.df.to_dict(orient='list')
            data['intermediates'] = self.intermediates.df.to_dict(orient='list')
            if products:
                data['products'] = self.products.df.to_dict(orient='list')
            data['compounds'] = self.compounds.df.to_dict(orient='list')

        data['reaction_ids'] = self.reactions.ids

        return data

    def write_json(
        self,
        file: 'str | Path',
        *,
        extra: dict | None = None,
        indent: str = '\t',
        **kwargs,
    ) -> None:
        """Serialise this recipe and write it to disk.

        :param file: write to this path
        :param extra: extra data to serialise
        :param indent: indentation whitespace (Default value = '\\t')
        """
        import json
        from pathlib import Path

        file = Path(file).resolve()
        assert file.parent.exists(), f'Directory does not exist: {file.parent}'

        data = self.get_dict(serialise_price=True, **kwargs)
        if extra:
            data.update(extra)

        mrich.writing(file)
        json.dump(data, open(file, 'w'), indent=indent)

    ### PRESENTATION

    def summary(self, price: bool = True) -> None:
        """Print a summary of this recipe

        :param price: print the price (Default value = True)
        """

        mrich.h1(str(self))

        if price:
            price = self.price
            if price:
                mrich.var('\nprice', price.amount, price.currency)

        if self.products:
            mrich.h3(f'{len(self.products)} products')
            if len(self.products) < 100:
                for product in self.products:
                    mrich.var(str(product.compound), f'{product.amount:.2f}', 'mg')

        if self.intermediates:
            mrich.h3(f'{len(self.intermediates)} intermediates')
            if len(self.intermediates) < 100:
                for intermediate in self.intermediates:
                    mrich.var(
                        str(intermediate.compound), f'{intermediate.amount:.2f}', 'mg'
                    )

        if self.reactants:
            mrich.h3(f'{len(self.reactants)} reactants')
            if len(self.reactants) < 100:
                for reactant in self.reactants:
                    mrich.var(str(reactant.compound), f'{reactant.amount:.2f}', 'mg')

        if self.reactions:
            mrich.h3(f'{len(self.reactions)} reactions')
            if len(self.reactions) < 100:
                for reaction in (Reaction(r) for r in self.reactions):
                    mrich.var(str(reaction), reaction.reaction_str, reaction.type)

        if self.compounds:
            mrich.h3(f'{len(self.compounds)} compounds')
            if len(self.compounds) < 100:
                for compound in self.compounds:
                    mrich.var(str(compound.compound), f'{compound.amount:.2f}', 'mg')

    def draw(self, color_mapper=None, node_size=300, graph_only=False):
        """Draw a graph of the reaction network

        :param color_mapper:  (Default value = None)
        :param node_size:  (Default value = 300)
        :param graph_only:  (Default value = False)
        """

        import networkx as nx

        color_mapper = color_mapper or {}
        colors = {}
        sizes = {}

        graph = nx.DiGraph()

        for reaction in (Reaction(r) for r in self.reactions):
            for reactant in reaction.reactants:
                key = str(reactant)
                ingredient = self.get_ingredient(id=reactant.id)
                graph.add_node(
                    key,
                    id=reactant.id,
                    smiles=reactant.smiles,
                    amount=ingredient.amount,
                    price=str(ingredient.price),
                    lead_time=ingredient.lead_time,
                )
                if not graph_only:
                    sizes[key] = ingredient.amount
                    colors[key] = color_mapper.get(key, (0.7, 0.7, 0.7))

        for product in self.products:
            key = str(product.compound)
            ingredient = self.get_ingredient(id=product.id)
            graph.add_node(
                key,
                id=product.id,
                smiles=product.smiles,
                amount=ingredient.amount,
                price=str(ingredient.price),
                lead_time=ingredient.lead_time,
            )
            if not graph_only:
                sizes[key] = product.amount
                colors[key] = color_mapper.get(key, (0.7, 0.7, 0.7))

        for reaction in (Reaction(r) for r in self.reactions):
            for reactant in reaction.reactants:
                graph.add_edge(
                    str(reactant),
                    str(reaction.product),
                    id=reaction.id,
                    type=reaction.type,
                    product_yield=reaction.product_yield,
                )

        if graph_only:
            return graph

        s_min = min(sizes.values())
        sizes = [s / s_min * node_size for s in sizes.values()]
        pos = nx.spring_layout(graph)
        return nx.draw(
            graph,
            pos=pos,
            with_labels=True,
            font_weight='bold',
            node_color=list(colors.values()),
            node_size=sizes,
        )

    def sankey(self, title: str | None = None) -> 'graph_objects.Figure':
        """Draw a plotly Sankey diagram

        :param title:  (Default value = None)
        """

        graph = self.draw(graph_only=True)

        import plotly.graph_objects as go

        nodes = {}
        for edge in graph.edges:
            for c in edge:
                if c not in nodes:
                    nodes[c] = len(nodes)

        source = [nodes[a] for a, b in graph.edges]
        target = [nodes[b] for a, b in graph.edges]
        value = [1 for _ in graph.edges]
        labels = list(nodes.keys())

        hoverkeys = None
        customdata = []
        for key in nodes.keys():
            n = graph.nodes[key]
            if not hoverkeys:
                hoverkeys = list(n.keys())
            if not n:
                mrich.error(f'problem w/ node {key=}')
                customdata.append((int(key[1:]), None))
            else:
                customdata.append(
                    tuple(v if v is not None else 'N/A' for v in n.values())
                )

        hoverkeys_edges = None
        customdata_edges = []
        for s, t in graph.edges.keys():
            edge = graph.edges[s, t]
            if not hoverkeys_edges:
                hoverkeys_edges = list(edge.keys())
            customdata_edges.append(
                tuple(v if v is not None else 'N/A' for v in edge.values())
            )

        hoverlines = [f'{key}=%{{customdata[{i}]}}' for i, key in enumerate(hoverkeys)]
        hovertemplate = 'Compound ' + '<br>'.join(hoverlines) + '<extra></extra>'

        hoverlines_edges = [
            f'{key}=%{{customdata[{i}]}}' for i, key in enumerate(hoverkeys_edges)
        ]
        hovertemplate_edges = (
            'Reaction ' + '<br>'.join(hoverlines_edges) + '<extra></extra>'
        )

        fig = go.Figure(
            data=[
                go.Sankey(
                    node=dict(
                        label=labels,
                        customdata=customdata,
                        hovertemplate=hovertemplate,
                    ),
                    link=dict(
                        customdata=customdata_edges,
                        hovertemplate=hovertemplate_edges,
                        source=source,
                        target=target,
                        value=value,
                    ),
                )
            ]
        )

        if not title:
            try:
                title = f'Recipe<br><sup>price={self.price}</sup>'
            except AssertionError:
                title = 'Recipe'

        fig.update_layout(title=title)
        return fig

    ### DEPRECATED — traversal/export shims (relocated to RecipeService)

    def get_routes(self, return_ids: bool = False) -> 'RouteSet':
        """DEPRECATED: use :meth:`.RecipeService.get_routes`."""
        from designdb.services.recipe import RecipeService

        return RecipeService.get_routes(self, return_ids=return_ids)

    def register_missing_routes(
        self, missing_only: bool = True, supplier: str = 'Enamine'
    ) -> None:
        """DEPRECATED: use :meth:`.RecipeService.register_missing_routes`."""
        from designdb.services.recipe import RecipeService

        return RecipeService.register_missing_routes(
            self, missing_only=missing_only, supplier=supplier
        )

    def write_CAR_csv(self, file: 'str | Path', return_df: bool = False):
        """DEPRECATED: use :meth:`.RecipeService.write_CAR_csv`."""
        from designdb.services.recipe import RecipeService

        return RecipeService.write_CAR_csv(self, file, return_df=return_df)

    def write_reactant_csv(
        self, file: 'str | Path', reaction_type_counts: bool = True, return_df=False
    ):
        """DEPRECATED: use :meth:`.RecipeService.write_reactant_csv`."""
        from designdb.services.recipe import RecipeService

        return RecipeService.write_reactant_csv(
            self, file, reaction_type_counts=reaction_type_counts, return_df=return_df
        )

    def write_product_csv(self, file: 'str | Path', return_df: bool = False):
        """DEPRECATED: use :meth:`.RecipeService.write_product_csv`."""
        from designdb.services.recipe import RecipeService

        return RecipeService.write_product_csv(self, file, return_df=return_df)

    def to_syndirella(self, out_key: 'str | Path', poses: 'PoseSet', *, separate=False):
        """DEPRECATED: use :meth:`.RecipeService.to_syndirella`."""
        from designdb.services.recipe import RecipeService

        return RecipeService.to_syndirella(self, out_key, poses, separate=separate)

    ### INTERNALS

    def __flag_modification(self) -> None:
        """Invalidate cached derived data after a mutation"""
        self._interactions = None
        self._score = None
        self._product_compounds = None
        self._poses = None
        self._combined_compounds = None

    ### DUNDERS

    def __str__(self) -> str:
        """Unformatted string representation"""
        s = f'(score={self.score:.3f})' if self.score else ''
        if self.hash:
            return f'Recipe_{self.hash}{s}'
        return f'Recipe{s}'

    def __longstr(self) -> str:
        """Long unformatted string representation"""

        if self.empty:
            return 'Empty Recipe()'

        if self.reactions:
            if self.intermediates:
                s = (
                    f'{self.reactants} --> {self.intermediates} --> '
                    f'{self.products} via {self.reactions}'
                )
            else:
                s = f'{self.reactants} --> {self.products} via {self.reactions}'

            if self.score:
                s += f', score={self.score:.3f}'
            if self.hash:
                return f'Recipe_{self.hash}({s})'
            return f'Recipe({s})'

        if self.hash:
            return f'Recipe_{self.hash}({self.compounds})'
        return f'Recipe(#compounds={self.num_compounds} [no-chem])'

    def __repr__(self) -> str:
        """ANSI Formatted string representation"""
        return (
            f'{mcol.bold}{mcol.underline}{self.__longstr()}'
            f'{mcol.unbold}{mcol.ununderline}'
        )

    def __rich__(self) -> str:
        """Rich Formatted string representation"""
        return f'[bold underline]{self.__longstr()}'

    def __add__(self, other: 'Recipe') -> 'Recipe':
        """Add another :class:`.Recipe` to this one"""
        result = self.copy()
        result.reactants += other.reactants
        result.intermediates += other.intermediates
        result.reactions += other.reactions
        result.products += other.products
        result.compounds += other.compounds
        return result


# name conflict with RouteModel. Trying to get rid of this entirely
class Route(Recipe):
    """A recipe with a single product, that is stored in the database"""

    def __init__(
        self,
        *,
        route_id: int,
        product: 'IngredientSet',
        reactants: 'IngredientSet',
        intermediates: 'IngredientSet',
        reactions: 'ReactionSet',
    ) -> None:
        """Route initialisation"""

        # check typing
        assert isinstance(product, IngredientSet)
        assert isinstance(reactants, IngredientSet)
        assert isinstance(intermediates, IngredientSet)
        assert isinstance(reactions, ReactionSet)

        assert len(product) == 1
        assert isinstance(route_id, int)
        assert route_id

        self._id = route_id
        self._products = product
        self._product_id = product.ids[0]
        self._reactants = reactants
        self._intermediates = intermediates
        self._reactions = reactions
        self._compounds = IngredientSet()
        self._hash = None
        self._score = None
        self._product_compounds = None
        self._poses = None
        self._interactions = None
        self._combined_compounds = None

    ### FACTORIES

    @classmethod
    def from_json(cls, path: 'str | Path | None' = None, data: dict = None) -> 'Route':
        """Load a serialised route from a JSON file

        :param path: path to JSON
        :param data: serialised data (Default value = None)
        """

        import json

        if data is None:
            data = json.load(open(path))

        self = cls.__new__(cls)

        self._id = data['id']
        self._product_id = data['product_id']
        self._products = IngredientSet.from_compounds(ids=[self._product_id])
        self._reactants = IngredientSet.from_json(
            path=None,
            data=data['reactants']['data'],
            supplier=data['reactants']['supplier'],
        )
        self._intermediates = IngredientSet.from_json(
            path=None,
            data=data['intermediates']['data'],
            supplier=data['intermediates']['supplier'],
        )
        self._reactions = ReactionSet(data['reactions']['indices'])
        self._compounds = IngredientSet()
        self._hash = None
        self._score = None
        self._product_compounds = None
        self._poses = None
        self._interactions = None
        self._combined_compounds = None

        return self

    @classmethod
    def get_route(
        cls, *, id: int, get_quote: bool = True, debug: bool = False
    ) -> 'Route':
        """Fetch a :class:`.RouteModel` stored in the database and wrap it.

        :param id: the ID of the :class:`.RouteModel` to retrieve
        :param get_quote: fetch catalogue quotes for the reactants so the route is
            priced, defaults to ``True``
        :param debug: increase verbosity for debugging
        """

        route = RouteModel.objects.get(pk=id)

        if debug:
            mrich.var('product_id', route.product_compound)

        qs = ComponentModel.objects.filter(route=route).order_by('id')

        reaction_ids = []
        reactant_ids = []
        reactant_amounts = []
        intermediate_ids = []
        intermediate_amounts = []

        for k in qs:
            ref = k.component_ref
            c_type = k.component_type
            amount = k.component_amount
            match c_type:
                case 1:
                    reaction_ids.append(ref)
                case 2:
                    reactant_ids.append(ref)
                    reactant_amounts.append(amount)
                case 3:
                    intermediate_ids.append(ref)
                    intermediate_amounts.append(amount)
                case _:
                    raise ValueError(f'Unknown component type {c_type}')

        if debug:
            mrich.var('components', qs)

        def _ingredients(ids, amounts, quote=False):
            """Build an IngredientSet, returning an empty one for no ids.

            When ``quote`` is set, each ingredient fetches its cheapest catalogue
            quote (via the ``compound_catalogue_map`` junction) so the route can be
            priced; otherwise ingredients are left unquoted.
            """
            iset = IngredientSet()
            for cid, amount in zip(ids, amounts, strict=False):
                iset.add(
                    Ingredient.from_compound(
                        compound=CompoundModel.objects.get(pk=cid),
                        amount=amount,
                        get_quote=quote,
                    )
                )
            return iset

        # products are made, not purchased, so they are never quoted
        products = IngredientSet.from_compounds(
            ids=[route.product_compound_id], amount=1
        )
        # reactants are the building blocks that get bought -> fetch quotes
        reactants = _ingredients(reactant_ids, reactant_amounts, quote=get_quote)
        intermediates = _ingredients(intermediate_ids, intermediate_amounts)

        reactions = ReactionSet(reaction_ids)

        recipe = Route(
            route_id=id,
            product=products,
            reactants=reactants,
            intermediates=intermediates,
            reactions=reactions,
        )

        if debug:
            mrich.var('recipe', recipe)

        return recipe

    ### PROPERTIES

    @property
    def product(self) -> 'Ingredient':
        """Product ingredient"""
        return self._products[0]

    @property
    def product_compound(self) -> 'CompoundModel':
        """Product compound"""
        return self.product.compound

    @property
    def id(self) -> int:
        """Route ID"""
        return self._id

    @property
    def price(self) -> 'Price':
        """Get the price of the reactants"""
        return self.reactants.price

    ### METHODS

    def get_dict(self) -> dict:
        """Serialisable dictionary"""
        return {
            'id': self.id,
            'product_id': self.product.id,
            'reactants': self.reactants.get_dict(),
            'intermediates': self.intermediates.get_dict(),
            'reactions': self.reactions.get_dict,
        }

    ### DUNDERS

    def __str__(self) -> str:
        """Unformatted string representation"""
        return f'Route #{self.id}: {self.product_compound}'

    def __repr__(self) -> str:
        """ANSI Formatted string representation"""
        return f'{mcol.bold}{mcol.underline}{self}{mcol.unbold}{mcol.ununderline}'

    def __rich__(self) -> str:
        """Rich Formatted string representation"""
        return f'[bold underline]{self}'


class RecipeSet:
    """A set of :class:`.Recipe` objects stored on disk as JSON."""

    def __init__(
        self,
        directory: 'str | Path',
        pattern: str = '*.json',
    ):
        """Load all recipes matching ``pattern`` in ``directory``."""

        from json import JSONDecodeError
        from pathlib import Path

        self._json_directory = Path(directory)
        self._json_pattern = pattern

        self._json_paths = {}
        for path in self._json_directory.glob(self._json_pattern):
            key = path.name.removeprefix('Recipe_').removesuffix('.json')
            self._json_paths[key] = path.resolve()

        mrich.reading(f'{directory}/{pattern}')

        self._recipes = {}
        for key, path in mrich.track(
            self._json_paths.items(), prefix='Loading recipes'
        ):
            try:
                recipe = Recipe.from_json(path=path, debug=False)
            except JSONDecodeError:
                mrich.error(f'Bad JSON in {path}')
                continue
            recipe._hash = key
            self._recipes[key] = recipe

        mrich.success('Loaded', len(self), 'Recipes')

    ### METHODS

    def get_values(
        self, key: str, progress: bool = False, serialise_price: bool = False
    ):
        """Get the value of attribute ``key`` for each member recipe."""
        values = []
        recipes = self._recipes.values()
        if progress:
            recipes = mrich.track(recipes, prefix=f'Calculating {self} values...')
        for recipe in recipes:
            value = getattr(recipe, key)
            if serialise_price and key == 'price':
                value = value.amount
            values.append(value)
        return values

    def get_df(self, **kwargs) -> 'pandas.DataFrame':
        """Get a dataframe of recipe dictionaries. See :meth:`.Recipe.get_dict`."""
        from pandas import DataFrame

        data = [recipe.get_dict(timestamp=False, **kwargs) for recipe in self]
        return DataFrame(data)

    def items(self) -> 'list[tuple[str, Recipe]]':
        """Data dictionary items"""
        return self._recipes.items()

    def keys(self) -> list[str]:
        """Data dictionary keys (recipe hashes)"""
        return self._recipes.keys()

    ### DUNDERS

    def __len__(self) -> int:
        """Number of recipes in this set"""
        return len(self._recipes)

    def __getitem__(self, key: int | str) -> Recipe:
        """Get a :class:`.Recipe` by index or hash"""
        match key:
            case int():
                return list(self._recipes.values())[key]
            case str():
                return self._recipes[key]
            case _:
                mrich.error(f'Unsupported RecipeSet key: {key=} {type(key)}')
        return None

    def __iter__(self):
        """Iterate over member recipes"""
        return iter(self._recipes.values())

    def __contains__(self, key: str) -> bool:
        """Is this hash present in the set?"""
        assert isinstance(key, str)
        return key in self._recipes

    def __str__(self) -> str:
        """Unformatted string representation"""
        return f'{{Recipe × {len(self)}}}'

    def __repr__(self) -> str:
        """ANSI Formatted string representation"""
        return f'{mcol.bold}{mcol.underline}{self}{mcol.unbold}{mcol.ununderline}'

    def __rich__(self) -> str:
        """Rich Formatted string representation"""
        return f'[bold underline]{self}'
