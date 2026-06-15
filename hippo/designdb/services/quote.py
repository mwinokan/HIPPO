"""Service layer for catalogue-price quoting.

In the modern DesignDB, catalogue prices live in the same database as the design
compounds (``catalogue_compounds`` / ``catalogue_prices``) and are linked to
compounds automatically by database triggers that match the registration hash
(``compounds.compound_hash == catalogue_compounds.catalogue_hash``), populating
the ``compound_catalogue_map`` junction.

Quoting is therefore a read against that junction: a compound is "quoted" if it
has at least one linked catalogue price. There is no longer any cross-database
transfer (the legacy ``HIPPO.quote_compounds(ref_animal)`` behaviour).
"""

from designdb.models import CataloguePriceCompoundJunctionModel


class QuoteService:
    """Report catalogue-price quoting for compounds using the current database."""

    @staticmethod
    def quoted_compound_ids(compound_ids: 'list[int] | set[int]') -> set[int]:
        """Return the subset of ``compound_ids`` that have at least one linked
        catalogue price.

        :param compound_ids: compound IDs to check
        :returns: the IDs that are quoted
        """
        return set(
            CataloguePriceCompoundJunctionModel.objects.filter(
                compound_id__in=list(compound_ids)
            )
            .values_list('compound_id', flat=True)
            .distinct()
        )

    @classmethod
    def partition_quoted(
        cls, compound_ids: 'list[int] | set[int]'
    ) -> tuple[set[int], set[int]]:
        """Split ``compound_ids`` into ``(quoted, unquoted)`` ID sets.

        :param compound_ids: compound IDs to partition
        :returns: ``(quoted_ids, unquoted_ids)``
        """
        ids = set(compound_ids)
        quoted = cls.quoted_compound_ids(ids)
        return quoted, ids - quoted
