from django.core.management.base import BaseCommand
from django.conf import settings
from django.utils import timezone
from django.db import transaction

import mrich
from mrich import print
from pathlib import Path

from rdkit.Chem import PandasTools, MolToSmiles
from hippo.tools import inchikey_from_smiles, sanitise_smiles
from hippo.models import *

# import pandas as pd
import re

# from rdkit.Chem import PandasTools
from rdkit.Chem.rdmolfiles import SDMolSupplier
from molparse.rdkit import mol_to_smiles

CATALOGS = {
    "Enamine_BB_catalog_EU_EUR": {
        "name": "EU stock",
        "lead_time": "2",
        "currency": "EUR",
    },
    "Enamine-FullCatalogue_EUR": {
        "name": "Global catalog",
        "lead_time": "10",
        "currency": "EUR",
    },
}

SUPPLIER = "Enamine"

PRICE_TO_MG_AMOUNT = {
    "Price_50mg": 50,
    "Price_0.1g": 100,
    "Price_0.25g": 250,
    "Price_0.5g": 500,
    "Price_1g": 1000,
    "Price_2.5g": 2500,
    "Price_5g": 5000,
    "Price_10g": 10000,
    "Price_EUR_100mg": 100,
    "Price_EUR_250mg": 250,
    "Price_EUR_500mg": 500,
    "Price_EUR_1g": 1000,
    "Price_EUR_2.5g": 2500,
    "Price_EUR_5g": 5000,
    "Price_EUR_10g": 10000,
}


class Command(BaseCommand):
    help = "Import Enamine quotes from catalogue SDF"

    def add_arguments(self, parser):
        parser.add_argument("file_path", type=str)

    def handle(
        self,
        file_path: str,
        debug: bool = False,
        *args,
        **kwargs,
    ):

        mrich.h1("Load SDF")
        mrich.var("file_path", file_path)

        infile = Path(file_path)
        key = infile.name.removesuffix(".sdf")

        catalog = CATALOGS[key]
        catalog_name = catalog["name"]
        lead_time = catalog["lead_time"]
        currency = catalog["currency"]
        mrich.var("catalog_name", catalog_name)
        mrich.var("lead_time", lead_time)
        mrich.var("currency", currency)

        compound_count = Compound.objects.count()
        quote_count = Quote.objects.count()

        suppl = SDMolSupplier(infile)
        n = len(suppl)
        mrich.var("Expecting", n, "compounds")

        batch_size = 10_000

        data = None

        supplier, _ = Supplier.objects.get_or_create(name=SUPPLIER)

        for i, mol in mrich.track(enumerate(suppl), total=n):

            if i % batch_size == 0:

                if data:
                    register(unique_smiles, data, supplier)

                data = []
                unique_smiles = set()

            if mol is None:
                mrich.error(f"Could not load molecule, {i=}")
                continue

            smiles = mol_to_smiles(mol)

            props = mol.GetPropsAsDict()

            if i == 0:
                mrich.print(props)

            entry = props["ID"]

            unique_smiles.add(smiles)

            for price_col in (p for p in props if p.startswith("Price_")):
                amount_in_mg = PRICE_TO_MG_AMOUNT[price_col]

                price = props[price_col]

                quote_data = dict(
                    smiles=smiles,
                    supplier="Enamine",
                    catalogue=catalog_name,
                    entry=entry,
                    amount=amount_in_mg,
                    purity=None,
                    lead_time=lead_time,
                    price=price,
                    currency=currency,
                )

                data.append(quote_data)

            mrich.set_progress_field("i", i)

        register(unique_smiles, data, supplier)

        mrich.success("Done.")

        mrich.var("#new compounds", Compound.objects.count() - compound_count)
        mrich.var("#new quotes", Quote.objects.count() - quote_count)
        mrich.success("Finished loading from", f"[file]{infile}")


def register(smiles, data, supplier):

    with mrich.loading("Registering compounds..."):
        smiles_map = Compound.bulk_register(smiles=smiles)

    mrich.var("#new quotes", len(data))

    quotes = []

    compounds = {
        c.smiles: c for c in Compound.objects.filter(smiles__in=smiles_map.values())
    }

    compounds = {
        s_old: compounds[smiles]
        for s_old, smiles in smiles_map.items()
        if smiles in compounds
    }

    with mrich.loading("Processing quotes..."):
        for d in data:

            compound = compounds.get(d["smiles"])

            if not compound:
                mrich.error("Could not register smiles:", d["smiles"])
                continue

            quotes.append(
                Quote(
                    supplier=supplier,
                    compound=compound,
                    smiles=d["smiles"],
                    catalogue=d["catalogue"],
                    entry=d["entry"],
                    amount=d["amount"],
                    lead_time=d["lead_time"],
                    price=d["price"],
                    currency=d["currency"],
                    date=timezone.now(),
                )
            )

    with mrich.loading("Bulk creating quotes..."):
        Quote.objects.bulk_create(
            quotes,
            ignore_conflicts=True,
            # unique_fields=["amount", "supplier", "catalogue", "entry"],
            # update_fields=["purity", "price", "date", "lead_time", "currency"],
        )
