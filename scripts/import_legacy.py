import hippo_legacy
import hippo
from pathlib import Path
from typer import Typer
from hippo.legacy import LEGACY_MAP
import mrich
from django.apps import apps
from pandas import DataFrame

app = Typer()


@app.command()
def migrate(source_db: str, target_db: str, skip: str = "", debug: bool = False):

    if skip:
        skip = skip.split(",")
    else:
        skip = []

    source_animal = hippo_legacy.HIPPO("source", source_db)

    target_animal = hippo.HIPPO("target", target_db)

    mrich.print(LEGACY_MAP)

    ID_MAP = {}

    for mapping in LEGACY_MAP:

        table = mapping["table"]

        if table in skip:
            continue

        mrich.h1(f"importing {table} table")

        debug = debug or mapping.get("debug", False)

        for warning in mapping.get("warnings", []):
            mrich.warning(warning)

        ID_MAP[table] = {}

        # get the django model
        model_name = mapping["model"]
        mrich.var("Django Model name", model_name)
        # model = apps.get_model(app_label="hippo", model_name=model_name)
        # mrich.var("Django Model", model)

        # get the registration method
        method_name = mapping["method"]
        method = getattr(target_animal, method_name)

        # print info
        field_map = mapping["fields"]

        for source_field, target_field in field_map.items():
            mrich.print(source_field, "-->", target_field)

        source_fields = list(field_map.keys())
        source_field_str = ",".join([f"{table}_id"] + source_fields)
        target_fields = list(field_map.values())

        related_fields = mapping.get("related_fields", {})

        # query for legacy data
        with mrich.loading("Querying legacy database"):
            records = source_animal.db.select(
                query=source_field_str, table=table, multiple=True
            )

        mrich.var(f"#{table} records", len(records))

        # bulk flag
        bulk = mapping.get("bulk", False)

        bulk_data = {}
        for record in mrich.track(records, prefix=f"registering {model_name} records"):

            source_id = record[0]
            record = record[1:]

            if debug:
                mrich.var("source_id", source_id)
                mrich.var("record", record)

            # create instance dictionary
            kwargs = {k: v for k, v in zip(target_fields, record)}

            for field, value in mapping.get("default_fields", {}).items():
                kwargs[field] = value

            for custom_field, custom_method in mapping.get("custom_fields", {}).items():
                kwargs = custom_method(kwargs)

            # get related object ID's
            for field in kwargs:
                if field in related_fields:
                    related_table = related_fields[field]
                    kwargs[field] = ID_MAP[related_table][kwargs[field]]

            if debug:
                mrich.debug("kwargs", kwargs)

            # store data for bulk insertion
            if bulk:
                bulk_data[source_id] = kwargs
                continue

            # instance = model(**kwargs)
            instance = method(debug=debug, **kwargs)
            instance.full_clean()
            instance.save()

            if debug:
                mrich.debug("Registered", instance)

            # store primary key mapping
            ID_MAP[table][source_id] = instance.id

        if bulk:

            # create dataframe
            df = DataFrame(bulk_data.values())

            # df = df[:50]

            # prepare bulk insertion data
            kwargs = {}
            for col in df.columns:
                kwargs[col] = df[col].values

            # run the bulk registration method
            target_objs = method(**kwargs)

            # insert PK mapping
            for source_id, target_obj in zip(bulk_data, target_objs):
                ID_MAP[table][source_id] = target_obj.id


def main():
    app()


if __name__ == "__main__":
    main()
