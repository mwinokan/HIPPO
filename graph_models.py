import os
from pathlib import Path
from hippo_django_orm.database import Database

db = Database("test.sqlite", include_contenttypes=True, include_templates=True)

# from django.template import loader
from django_extensions.management.modelviz import ModelGraph, generate_dot

app_label = "hippo"


def main():

    graph_models = ModelGraph([app_label])

    graph_models.generate_graph_data()

    graph_data = graph_models.get_graph_data(as_json=True)

    import mrich

    mrich.print(graph_data)

    # theme = "django2018"

    # import django_extensions

    # template_path = Path(django_extensions.__file__).parent / "templates" / "django_extensions" / "graph_models" / theme / "digraph.dot"

    # from django.template.base import Template

    # Engine()
    # with open(template_path, "rt") as f:
    # template = Template(f.read())

    # template = loader.get_template(str(template_path.resolve()))

    # dotdata = generate_dot(graph_data, template=template)

    # with open("hippo_django_orm.dot", "wt") as f:
    #     f.write(dotdata)


if __name__ == "__main__":
    main()
