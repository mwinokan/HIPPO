from pathlib import Path

import logging

logger = logging.getLogger("HIPPO")


class ProposalPage:
    """
    Recipe proposal web page

    """

    def __init__(
        self,
        output_dir: str | Path,
        *,
        animal: "HIPPO",
        scaffolds: "CompoundSet | None" = None,
        suppliers: list[str] | None = None,
        starting_recipe: "Recipe | None" = None,
        rgen: "RandomRecipeGenerator | None" = None,
        scorer: "Scorer | None" = None,
        title: str | None = None,
    ):

        # setup directories
        self._output_dir = Path(output_dir)
        self.make_directories()

        # project objects
        self._animal = animal
        self._scaffolds = scaffolds
        self._suppliers = suppliers
        self._starting_recipe = starting_recipe
        self._rgen = rgen
        self._scorer = scorer

        self._title = title or animal.name

        # yattag setup
        from yattag import Doc

        doc, tag, text, line = Doc().ttl()

        self._doc = doc
        self._tag = tag
        self._text = text
        self._line = line

        self.setup_page()

    ### FACTORIES

    ### PROPERTIES

    @property
    def animal(self) -> "HIPPO":
        """associated :class:`.HIPPO` object"""
        return self._animal

    @property
    def db(self) -> "Database":
        """associated :class:`.Database` object"""
        return self.animal.db

    @property
    def doc(self) -> "yattag.Doc":
        """yattag.Doc"""
        return self._doc

    @property
    def tag(self) -> "yattag.tag":
        """yattag.tag"""
        return self._tag

    @property
    def text(self) -> "yattag.text":
        """yattag.text"""
        return self._text

    @property
    def line(self) -> "yattag.line":
        """yattag.line"""
        return self._line

    @property
    def title(self) -> str:
        """Page title"""
        return self._title

    @property
    def output_dir(self) -> "Path":
        """Output directory"""
        return self._output_dir

    @property
    def resource_dir(self) -> "Path":
        """Output directory"""
        return self.output_dir / "web_resources"

    @property
    def index_path(self) -> "Path":
        """index.html Path"""
        return self.output_dir / "index.html"

    ### METHODS

    def write_html(self) -> None:
        """Write the index.html file"""

        from yattag import indent

        path = self.index_path

        with open(path, "wt") as f:
            logger.writing(path)
            f.writelines(indent(self.doc.getvalue()))

    ### INTERNAL HTML STUFF

    def make_directories(self) -> None:
        """Create output directories"""

        logger.writing(self.output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        logger.writing(self.resource_dir)
        (self.resource_dir).mkdir(parents=True, exist_ok=True)

    def setup_page(self) -> None:
        """Create the yattag page content"""

        self.doc.asis("<!DOCTYPE html>")

        with self.tag("html"):

            self.header()

            with self.tag("body", klass="w3-content", style="max-width:none"):

                with self.tag("div", klass="w3-container w3-teal"):
                    with self.tag("h1"):
                        self.text(self.title)

                with self.tag("div", klass="w3-container w3-padding"):
                    self.section(self.sec_targets)
                    self.section(self.sec_hits)

    def header(self) -> None:
        """Create the page header"""

        with self.tag("head"):

            self.doc.stag("meta", charset="UTF-8")
            self.doc.stag(
                "meta", name="viewport", content="width=device-width, initial-scale=1"
            )
            self.doc.stag(
                "link",
                rel="stylesheet",
                href="https://www.w3schools.com/w3css/4/w3.css",
            )
            self.doc.stag(
                "link",
                rel="stylesheet",
                href="https://fonts.googleapis.com/css?family=Oswald",
            )
            self.doc.stag(
                "link",
                rel="stylesheet",
                href="https://fonts.googleapis.com/css?family=Open Sans",
            )
            self.doc.stag(
                "link",
                rel="stylesheet",
                href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/4.7.0/css/font-awesome.min.css",
            )

            with self.tag("script", src="https://cdn.plot.ly/plotly-latest.min.js"):
                ...

            self.style()

    def style(self) -> None:
        """Create the page style"""

        # change to a .css file and use doc.stag("link", rel="stylesheet", href="style.css")

        with self.tag("style"):
            self.doc.asis(
                """h1,h2,h3,h4,h5,h6 {font-family: "Oswald"}body {font-family: "Open Sans"}"""
            )

    def section_header(self, title: str, tag: str = "h2") -> None:
        """section header"""
        with self.tag(tag):
            self.text(str(title))

    def accordion(self) -> None:
        """sub-content accordion"""
        # https://www.w3schools.com/w3css/w3css_accordions.asp
        raise NotImplementedError

    def sidebar(self) -> None:
        # https://www.w3schools.com/w3css/w3css_sidebar.asp
        raise NotImplementedError

    def var(self, key, value, tag=None, separator=": ") -> None:
        """sub-content accordion"""
        text = f"{key}{separator}{value}"

        if not tag:
            with self.tag("b"):
                self.text(key)
            self.text(separator)
            self.text(str(value))
        else:
            with self.tag(tag):
                self.var(key, value, separator=separator)

    def section(self, function) -> None:
        """create section div"""

        with self.tag("div", klass="w3-panel w3-border"):
            # with self.tag("div", klass="w3-border w3-padding"):
            function()
        # self.doc.stag("br")

    def plotly_graph(self, figure, filename):

        # from plotly.offline import plot

        f"""<div> 
            <iframe src="iter1b_go/scoring_newfeat_coverage.html" width="100%" height="500" style="border:none;"></iframe>
        </div>"""

        path = self.resource_dir / filename
        rel_path = Path(self.resource_dir.name) / filename

        # write the html file
        from molparse import write

        logger.writing(path)
        write(path, figure, verbosity=0)

        # embed the graph
        with self.tag("div"):
            with self.tag(
                "iframe",
                src=str(rel_path),
                width="100%",
                height="500",
                style="border:none",
            ):
                ...

    def table(self, data, **kwargs):
        """Embed some data as a table"""
        from pandas import DataFrame

        df = DataFrame(data)
        html = df.to_html(**kwargs, classes="w3-table-all", index=False, escape=False)
        self.doc.asis(html)
        self.doc.asis("<br>")

    ### SECTION CONTENT

    def sec_targets(self) -> None:
        """Section on targets"""

        title = "Protein Target"
        targets = self.animal.targets

        if len(targets) > 1:
            title += "s"

        self.section_header(title)

        for target in targets:
            self.section_header(target.name, "h3")

            self.var("name", target.name)

            subsites = target.subsites

            if subsites:

                self.section_header("Subsites", "h4")
                with self.tag("ul"):
                    for subsite in subsites:
                        self.var(f"Site {subsite.id}", subsite.name, tag="li")

    def sec_hits(self) -> None:
        """Section on experimental hits"""

        title = "Experimental hits"
        hit_compounds = self.animal.compounds(tag="hits")
        hit_poses = self.animal.poses(tag="hits")

        self.section_header(title)

        with self.tag("ul"):
            self.var("#compounds", len(hit_compounds), tag="li")
            self.var("#observations", len(hit_poses), tag="li")

        from .animal import GENERATED_TAG_COLS

        # tag statistics
        fig = self.animal.plot_tag_statistics(
            show_compounds=False,
            poses=hit_poses,
            logo=None,
            title="Tags",
            skip=["Pose", "hits"] + GENERATED_TAG_COLS,
        )

        self.plotly_graph(fig, "hit_tags.html")

        # files
        self.section_header("Downloads", "h3")
        path = self.resource_dir / "hit_poses.sdf"
        rel_path = Path(self.resource_dir.name) / "hit_poses.sdf"
        hit_poses.write_sdf(path, inspirations=False)
        table_data = [
            dict(
                Name="hit_poses.sdf",
                Description="SDF of the experimental hits",
                Download=f'<a href="{rel_path}" download>SDF</a>',
            )
        ]
        self.table(table_data)

    ### DUNDERS
