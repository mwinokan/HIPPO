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
        # scaffolds: "CompoundSet",
        # suppliers: list[str],
        # starting_recipe: "Recipe",
        # rgen: "RandomRecipeGenerator",
        # scorer: "Scorer",
        title: str | None = None,
    ):

        # setup directories
        self._output_dir = Path(output_dir)
        self.make_directories()

        # project objects
        self._animal = animal
        # self._scaffolds = scaffolds
        # self._suppliers = suppliers
        # self._starting_recipe = starting_recipe
        # self._rgen = rgen
        # self._scorer = scorer

        self._title = title or animal.name

        # yattag setup
        from yattag import Doc

        doc, tag, text = Doc().tagtext()

        self._doc = doc
        self._tag = tag
        self._text = text

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
    def title(self) -> str:
        """Page title"""
        return self._title

    @property
    def output_dir(self) -> "Path":
        """Output directory"""
        return self._output_dir

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

        logger.writing(self.output_dir / "web_resources")
        (self.output_dir / "web_resources").mkdir(parents=True, exist_ok=True)

    def setup_page(self) -> None:
        """Create the yattag page content"""

        self.doc.asis("<!DOCTYPE html>")

        with self.tag("html"):

            self.header()

            with self.tag("body"):

                with self.tag("h1"):
                    self.text(self.title)

                self.text("Hello world!")

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

    ### DUNDERS
