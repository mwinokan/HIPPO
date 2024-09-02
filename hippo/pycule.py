import pycule
from .quote import Quote
import requests

import logging

# logging.getLogger("requests").setLevel(logging.WARNING)
# logging.getLogger("mcule:core").setLevel(logging.WARNING)
logging.getLogger("urllib3").setLevel(logging.WARNING)

ENAMINE_CATALOGUES = [
    "BB",
    "SCR",
    "MADE",
    "REAL",
]

ENAMINE_V2_CATALOGUES = {
    "building-blocks": "BB",
    "screening-compounds": "SCR",
    "made-building-blocks": "MADE",
    "real": "REAL",
    "bioactive-compounds": "BIO",
}

ENAMINE_V2_BASE_URL = "https://new.enaminestore.com/api/v2"

ENAMINE_V2_LEAD_TIME = {
    "BB": {
        "IN_STOCK": 5,
        "SYNTHESIS": 7,
    },
    "SCR": {
        "IN_STOCK": 5,
        "SYNTHESIS": 20,
    },
    "MADE": {
        "SYNTHESIS": 20,
    },
    "REAL": {
        "SYNTHESIS": 20,
    },
    "BIO": {
        "IN_STOCK": 5,
    },
}

import logging

logger = logging.getLogger("HIPPO")


class Quoter:
    """Class to scrape catalogue data"""

    def __init__(
        self,
        supplier: str,
        username: str = None,
        password: str = None,
        token: str = None,
    ):
        """Currently 'Enamine' and 'MCule' are supported. Once initialised"""

        if supplier == "Enamine":

            self.supplier = "Enamine"
            self.batch_size = 2000
            self.catalogues = ENAMINE_CATALOGUES
            self._get_quote = self.get_enamine_quote
            self._get_batch_quote = self.get_enamine_batch_quote

            assert username
            assert password

            self.wrapper = pycule.core.EnamineWrapper(
                username=username, password=password
            )

        elif supplier == "MCule":

            self.supplier = "MCule"
            self.batch_size = 1000
            self._get_quote = self.get_mcule_quote
            self._get_batch_quote = self.get_mcule_batch_quote

            assert token

            self.wrapper = pycule.core.MCuleWrapper(authorisation_token=token)

        else:
            logger.error(f'Unsupported supplier: "{supplier}"')
            raise NotImplementedError

    def get_quote(self, compound, **kwargs):
        """Get quotes for a compound

        :param compound:

        """
        return self._get_quote(compound, **kwargs)

    def get_batch_quote(self, compounds, **kwargs):
        """Get quotes for a compound set

        :param compounds:

        """
        return self._get_batch_quote(compounds, **kwargs)

    def __call__(self, compound):
        """Get quotes for a compound"""
        return self.get_quote(compound)

    ### ENAMINE

    def get_enamine_batch_quote(
        self,
        compounds,
        currency="USD",
        catalogues=None,
        exact=False,
        forms=False,
        analogues=False,
        equivalents=False,
    ):
        """

        :param compounds:
        :param currency:  (Default value = "USD")
        :param catalogues:  (Default value = None)
        :param exact:  (Default value = False)
        :param forms:  (Default value = False)
        :param analogues:  (Default value = False)
        :param equivalents:  (Default value = False)

        """

        import time
        from .cset import CompoundSet

        n_comps = len(compounds)
        logger.var("batch size", n_comps)
        assert n_comps <= 2000

        unmatched = compounds.copy()
        matched = CompoundSet(compounds.db)
        errors = []

        # perform an exact search on the batch of smiles (for each catalog)

        payload = dict(compounds=compounds.smiles, currency=currency)

        include = []
        if forms:
            include.append("FORMS")
        if analogues:
            include.append("ANALOGS")
        if equivalents:
            include.append("EQUIVALENTS")
        if include:
            payload["include"] = include

        catalogues = catalogues or ENAMINE_V2_CATALOGUES

        t_start = time.perf_counter()

        # return payload

        for catalogue in catalogues:

            if exact:
                logger.debug(
                    f"Exact search in Enamine:{ENAMINE_V2_CATALOGUES[catalogue]}..."
                )
            else:
                logger.debug(
                    f"Similarity==1.0 search in Enamine:{ENAMINE_V2_CATALOGUES[catalogue]}..."
                )

            # send the request
            t_catalogue_start = time.perf_counter()
            if exact:
                url = (
                    ENAMINE_V2_BASE_URL
                    + f"/catalog/search/in/{catalogue}/by/SMILEs/EXACT"
                )
            else:
                url = (
                    ENAMINE_V2_BASE_URL
                    + f"/catalog/search/in/{catalogue}/by/SMILEs/SIMILARITY/1.00"
                )
            response = requests.post(url=url, json=payload)
            logger.var(
                f"{catalogue} request time = ", time.perf_counter() - t_catalogue_start
            )

            # process the request
            if (status := response.status_code) != 200:
                logger.error(f"Request error: status={status}")
                errors.append(
                    dict(catalogue=catalogue, status=status, response=response)
                )
                continue

            data = response.json()

            results = data["results"]
            products = results["products"]

            # register the quotes in the database
            this_matched = self.parse_enamine_bulk_products(compounds.db, products)

            if this_matched:
                logger.success(f"Added quotes for {len(this_matched)} compounds")

            matched += this_matched
            unmatched -= this_matched

            return data

            # break

        logger.var("Total time = ", time.perf_counter() - t_start)

        if unmatched:
            logger.error(f"Did not find quotes for {len(unmatched)} compounds")

        if errors:
            logger.error(f"{len(errors)} failed requests")

        return dict(matched=matched, unmatched=unmatched, errors=errors)

    def parse_enamine_bulk_products(self, db, products):
        """

        :param db:
        :param products:

        """

        from .tools import flat_inchikey
        from .cset import CompoundSet

        matched = CompoundSet(db)

        for product in products:

            prices = product["prices"]
            product = product["product"]

            # quote-level data from product
            entry = product["code"]
            smiles = product["smile"]
            currency = prices["currency"]
            catalogue = ENAMINE_V2_CATALOGUES[product["catalog"]]

            # check if the product has a name
            if not entry:
                logger.warning(f"Product in {catalogue} but no entry-name. {smiles=}")

            if p := product["purity"]:
                purity = product["purity"] / 100
            else:
                purity = None

            assert catalogue, product
            assert currency, product

            # identify the compound
            compound_id = db.get_compound_id(inchikey=flat_inchikey(smiles))
            matched.add(compound_id)

            # loop through prices
            for pack in prices["all"]:

                assert pack["weight"]["measure"] == "mg"
                amount = pack["weight"]["amount"]
                price = pack["price"]

                availability_str = pack["weight"]["available"]

                try:
                    lead_time = ENAMINE_V2_LEAD_TIME[catalogue][availability_str]
                except KeyError:
                    logger.error(
                        f"Unsupported {availability_str=} [{compound_id=}, {entry=}]"
                    )
                    print(pack)
                    continue

                # register the quote
                quote_id = db.insert_quote(
                    compound=compound_id,
                    supplier="Enamine",
                    catalogue=catalogue,
                    entry=entry,
                    amount=amount,
                    price=price,
                    currency=currency,
                    purity=purity,
                    lead_time=lead_time,
                    smiles=smiles,
                    commit=False,
                )

            # break

        logger.debug("Committing changes to the database...")
        db.commit()

        return matched

    def get_enamine_quote(self, compound, currency="USD"):
        """

        :param compound:
        :param currency:  (Default value = "USD")

        """

        assert compound._table == "compound"

        try:

            smiles = compound.smiles

            logger.header(f"Exact search: {smiles=}")

            for catalogue in self.catalogues:

                logger.header(f"Searching in {self.supplier} {catalogue}...")

                try:
                    result = self.wrapper.exactsearch(
                        smiles=compound.smiles, currency=currency, catalogue=catalogue
                    )
                except requests.ConnectionError as e:
                    logger.error(f"ConnectionError: {e}")
                    return None

                # try again but with a similarity search
                if result["response"]["result"]["code"] != 0:
                    try:
                        result = self.wrapper.similaritysearch(
                            smiles, similarity_value=1.0, catalogue=catalogue
                        )
                    except requests.ConnectionError as e:
                        logger.error(f"ConnectionError: {e}")
                        return None

                if result["response"]["result"]["code"] != 0:
                    continue

                data = result["response"]["data"]

                if len(data) < 1:
                    print(result["response"])
                    raise NotImplementedError

                data = data[0]

                name = data["Id"]

                logger.success(f"Found in {catalogue} w/ id={name}")

                break

            else:

                logger.error("No match found")
                return None

            ### get the information

            logger.header(f"Searching by ID in {self.supplier} {catalogue}...")

            for catalogue in self.catalogues:

                result = self.wrapper.compoundidsearch(name, catalogue=catalogue)

                if result["response"]["result"]["code"] == 0:
                    try:
                        logger.success(f"Found in {catalogue}")
                        return (
                            self.parse_enamine_response(result, compound, catalogue),
                            result,
                        )
                    except NoDataInReponse:
                        pass

            else:

                logger.error("No catalog entry")
                return None

        except KeyboardInterrupt:
            logger.warning("Interrupted quoting")

    def pick_enamine_exact_data(self, smiles, data):
        """

        :param smiles:
        :param data:

        """

        if len(data) == 1:
            return data[0]["Id"]

        if len(data) == 0:
            raise NoDataInResponse(smiles)

        # sort the results by increasing lead_time, and increasing stereo-complexity
        data = sorted(
            data,
            key=lambda x: (
                self.parse_enamine_delivery_string(x["deliveryDays"]),
                stereo_complexity(x["smile"]),
            ),
        )
        mout.warning(
            f'Picking quickest and least stereo complex result from: {[x["Id"] for x in data]}'
        )
        return data[0]["Id"]

    def parse_enamine_response(self, result, compound, catalogue):
        """

        :param result:
        :param compound:
        :param catalogue:

        """

        quotes = []

        for i, data in enumerate(result["response"]["data"]):

            if i == 0:

                # update compound metadata
                metadata = dict()
                metadata["enamine_id"] = data["Id"]
                metadata["alias"] = data["name"]
                metadata["cas"] = data["cas"]
                metadata["mfcd"] = data["mfcd"]
                metadata["formula"] = data["formula"]
                metadata["molecular_weight"] = data["mw"]
                compound.metadata.update(metadata)

                """
				do something with: 
				* availability
				* storageCond
				* productUrl
				* lastUpdate
				"""

            lead_time = self.parse_enamine_delivery_string(
                delivery := data["deliveryDays"]
            )
            entry = data["Id"]
            smiles = data["smile"]

            if purity := data["purity"]:
                purity /= 100

            # if not (packs := data['packs']):

            if not (packs := data["packs"]) and "synthesis" in delivery:

                logger.warning("Hardcoded REAL quote")

                compound.db.insert_quote(
                    compound=compound,
                    supplier="Enamine",
                    catalogue=catalogue,
                    entry=entry,
                    amount=20,
                    price=330,
                    currency="USD",
                    purity=None,
                    lead_time=15,
                    smiles=smiles,
                )

            else:
                for pack in packs:
                    self.parse_enamine_pack(
                        compound, entry, purity, catalogue, pack, lead_time, smiles
                    )

    def pick_enamine_data(self, comp_id, data):
        """

        :param comp_id:
        :param data:

        """

        data = [d for d in data if d["Id"] == comp_id]

        if len(data) == 0:
            raise NoDataInReponse

        elif len(data) == 1:
            return data[0]

        data = [d for d in data if not any([p["price"] == 0.0 for p in d["packs"]])]

        if len(data) > 1:
            logger.warning(
                f"Taking first data entry after stripping non-priced ({comp_id})"
            )
        return data[0]

    def parse_enamine_delivery_string(self, string):
        """

        :param string:

        """

        if string.startswith("regular delivery, ") and string.endswith("days"):
            return int(string.removeprefix("regular delivery, ").removesuffix(" days"))

        if string.startswith("backorder, ") and string.endswith("days"):
            return int(string.removeprefix("backorder, ").removesuffix(" days"))

        if string == "3 weeks synthesis time":
            return 15

        raise Exception(f"Unexpected delivery string: {string}")

    def parse_enamine_pack(
        self, compound, entry, purity, catalogue, pack, lead_time, smiles
    ):
        """

        :param compound:
        :param entry:
        :param purity:
        :param catalogue:
        :param pack:
        :param lead_time:
        :param smiles:

        """

        supplier = "Enamine"

        match pack["measure"]:
            case "g":
                amount = pack["amount"] * 1000
            case "mg":
                amount = pack["amount"]
            case _:
                raise NotImplementedError(
                    f'Not supported pack measure: {pack["measure"]}'
                )

        price = pack["price"]
        currency = pack["currencyName"]

        if pack["status"] not in ["Normal", "Ice pack", ""]:
            print(pack)
            logger.warning(f"{pack['status']=}")

        if not price:
            logger.warning(f"Skipping price-less Enamine pack ({entry})")
            return None

        compound.db.insert_quote(
            compound=compound,
            supplier=supplier,
            catalogue=catalogue,
            entry=entry,
            amount=amount,
            price=price,
            currency=currency,
            purity=purity,
            lead_time=lead_time,
            smiles=smiles,
        )

    ### MCULE

    def get_mcule_quote(self, compound, exact=False):
        """

        :param compound:
        :param exact:  (Default value = False)

        """

        try:

            # single exact query

            smiles = compound.smiles

            if exact:
                logger.header(f"Exact search: {smiles=}")
                result = self.wrapper.singlequerysearch(smiles)
            else:
                logger.header(f"Similarity==1.0 search: {smiles=}")
                result = self.wrapper.similaritysearch(smiles, threshold=1.0)

            if result["response"]["results"]:

                results = result["response"]["results"]

                logger.header(results)

                mcule_ids = [data["mcule_id"] for data in results]
                smiles = [data["smiles"] for data in results]

                logger.success(f"Found w/ ids={mcule_ids}")

                for mcule_id, smile in zip(mcule_ids, smiles):

                    # logger.header(f'{self.wrapper.compoundavailability(mcule_id)=}')
                    # logger.header(f'{self.wrapper.compounddetails(mcule_id)=}')
                    result = self.wrapper.compoundprices(mcule_id)

                    for pack in result["response"]["best_prices"]:
                        self.parse_mcule_pack(
                            compound.db, compound, mcule_id, pack, smile
                        )

                    # logger.header(f'{self.wrapper.compoundpricesamount(mcule_id)=}')

            return result

        except KeyboardInterrupt:
            logger.warning("Interrupted quoting")

    def get_mcule_batch_quote(self, compounds):
        """

        :param compounds:

        """

        import time
        from .cset import CompoundSet
        from .tools import flat_inchikey
        from tqdm import tqdm

        db = compounds.db

        n_comps = len(compounds)

        matched = CompoundSet(db)

        logger.var("batch size", n_comps)
        assert n_comps <= 1000

        # bulk query compound availability

        t_start = time.perf_counter()
        logger.title("MCule bulk availability query...")
        data = self.wrapper.multiplequerieswithavailability(compounds.smiles)

        logger.var("availability query time", time.perf_counter() - t_start)

        results = data["response"]["results"]

        logger.var("#results", len(results))

        logger.title("MCule single price queries...")
        for result in tqdm(results):

            entry = result["mcule_id"]
            smiles = result["smiles"]

            # match to compound
            compound_id = db.get_compound_id(inchikey=flat_inchikey(smiles))

            if not compound_id:
                logger.error(f"Could not find compound matching {smiles=}")
                continue

            matched.add(compound_id)

            prices = self.wrapper.compoundprices(entry)

            try:
                best_prices = prices["response"]["best_prices"]
            except KeyError as e:
                print(prices)
                logger.error(f"Unsupported prices {e}")
                continue

            for pack in prices["response"]["best_prices"]:
                try:
                    self.parse_mcule_pack(
                        db, compound_id, entry, pack, smiles, commit=False
                    )
                except KeyError as e:
                    print(pack)
                    logger.error(f"Unsupported pack {e}")
                    continue

        logger.var("Total time = ", time.perf_counter() - t_start)

        unmatched = compounds - matched

        if unmatched:
            logger.error(f"Did not find quotes for {len(unmatched)} compounds")

        logger.success("Committing changes to the database...")
        db.commit()

        return dict(matched=matched, unmatched=unmatched)

    def parse_mcule_pack(self, db, compound, entry, pack, smiles, commit=True):
        """

        :param db:
        :param compound:
        :param entry:
        :param pack:
        :param smiles:
        :param commit:  (Default value = True)

        """

        supplier = "MCule"

        match pack["unit"]:
            case "g":
                amount = pack["amount"] * 1000
            case "mg":
                amount = pack["amount"]
            case _:
                raise NotImplementedError(f'Not supported pack measure: {pack["unit"]}')

        price = pack["price"]
        currency = pack["currency"]

        if "purity" in pack:
            purity = pack["purity"] / 100
        else:
            purity = None

        lead_time = pack["delivery_time_working_days"]

        db.insert_quote(
            compound=compound,
            supplier=supplier,
            catalogue=None,
            entry=entry,
            amount=amount,
            price=price,
            currency=currency,
            purity=purity,
            lead_time=lead_time,
            smiles=smiles,
            commit=commit,
        )


class NoDataInReponse(Exception):
    """ """

    ...


class NotInCatalogues(Exception):
    """ """

    ...
