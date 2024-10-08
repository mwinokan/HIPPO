from pathlib import Path

import logging

logger = logging.getLogger("HIPPO")
logging.getLogger("PIL").setLevel(logging.WARNING)


class ProjectPage:
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
        proposals: "list[Recipe] | None" = None,
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
        self._proposals = proposals

        self._title = title or animal.name

        # yattag setup
        from yattag import Doc

        doc, tag, text, line = Doc().ttl()

        self._doc = doc
        self._tag = tag
        self._text = text
        self._line = line

        self.setup_page()
        self.write_html()

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
    def mol_image_dir(self) -> "Path":
        """Output directory"""
        return self.resource_dir / "mol_images"

    @property
    def pose_sdf_dir(self) -> "Path":
        """Output directory"""
        return self.resource_dir / "pose_sdfs"

    @property
    def index_path(self) -> "Path":
        """index.html Path"""
        return self.output_dir / "index.html"

    @property
    def proposals(self):
        return self._proposals

    @property
    def scaffolds(self):
        return self._scaffolds

    @property
    def suppliers(self):
        return self._suppliers

    @property
    def starting_recipe(self):
        return self._starting_recipe

    @property
    def rgen(self):
        return self._rgen

    @property
    def scorer(self):
        return self._scorer

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

        logger.writing(self.mol_image_dir)
        (self.mol_image_dir).mkdir(parents=True, exist_ok=True)

        logger.writing(self.pose_sdf_dir)
        (self.pose_sdf_dir).mkdir(parents=True, exist_ok=True)

    def setup_page(self) -> None:
        """Create the yattag page content"""

        self.doc.asis("<!DOCTYPE html>")

        with self.tag("html"):

            self.header()

            with self.tag("body", klass="w3-content", style="max-width:none"):

                with self.tag("div", klass="w3-bar w3-teal"):
                    with self.tag("div", klass="w3-bar-item"):
                        src = "https://github.com/mwinokan/HIPPO/raw/main/logos/hippo_assets-02.png?raw=true"
                        self.doc.stag(
                            "img", src=src, style="max-height:75px"
                        )  # , klass="w3-image")

                    with self.tag("div", klass="w3-bar-item"):
                        with self.tag("h1"):
                            self.text(self.title)

                with self.tag("div", klass="w3-container w3-dark-gray w3-padding"):
                    self.section(self.sec_targets)
                    self.section(self.sec_hits)

                    # placeholders
                    if self.scaffolds:
                        self.section(self.sec_scaffolds)
                    if self.scaffolds:
                        self.section(self.sec_elaborations)
                    # if self.quoting: self.section(self.sec_quoting)
                    # if self.product_pool: self.section(self.sec_product_pool)
                    # if self.route_pool: self.section(self.sec_route_pool)
                    if self.rgen:
                        self.section(self.sec_rgen)
                    if self.scorer:
                        self.section(self.sec_scorer)
                    if self.proposals:
                        self.section(self.sec_proposals)

                with self.tag("div", klass="w3-container w3-teal w3-padding"):
                    with self.tag("div", klass="w3-center"):
                        src = "https://github.com/mwinokan/HIPPO/raw/main/logos/hippo_logo_tightcrop.png?raw=true"
                        self.doc.stag("img", src=src, style="max-height:150px")

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

            with self.tag("script", src="https://3Dmol.org/build/3Dmol-min.js"):
                ...

            with self.tag("script", src="https://3Dmol.org/build/3Dmol.ui-min.js"):
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
        with self.tag("div", klass="w3-panel w3-border w3-white"):
            function()

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

    def table(self, data, style: str = "w3-table-all", **kwargs):
        """Embed some data as a table"""
        from pandas import DataFrame

        df = DataFrame(data)
        html = df.to_html(**kwargs, classes=style, index=False, escape=False)
        self.doc.asis(html)
        self.doc.asis("<br>")

    # def mol_grid_svg(self, cset, **kwargs):

    #     mols = [c.mol for c in cset]

    #     from rdkit.Chem.Draw import MolsToGridImage

    #     return MolsToGridImage(mols,
    #         molsPerRow=3,
    #         subImgSize=(200,
    #         200),
    #         legends=None,
    #         # highlightAtomLists=None,
    #         # highlightBondLists=None,
    #         useSVG=True,
    #         returnPNG=False,
    #         **kwargs)

    def save_compound_image(self, compound):
        from rdkit.Chem.Draw import MolToImage

        image = MolToImage(compound.mol)
        path = self.mol_image_dir / f"C{compound.id}.png"
        logger.writing(path)
        image.save(path)

    def save_pose_sdf(self, pose):
        path = self.pose_sdf_dir / f"P{pose.id}.sdf"
        self.animal.poses([pose.id]).write_sdf(path, inspirations=False)

    def save_pset_sdf(self, name, pset):
        path = self.pose_sdf_dir / f"{name}.sdf"
        pset.write_sdf(path, inspirations=False)

    def compound_image(self, compound, max_height="250px"):
        self.save_compound_image(compound)

        src = str(
            Path(self.resource_dir.name)
            / Path(self.mol_image_dir.name)
            / f"C{compound.id}.png"
        )

        self.doc.stag("img", src=src, style=f"max-height:{max_height}")

    def pose_3d_view(self, pose):

        self.save_pose_sdf(compound)

        src = str(
            Path(self.resource_dir.name)
            / Path(self.mol_image_dir.name)
            / f"C{compound.id}.png"
        )

        self.doc.stag("img", src=src, style=f"max-height:{max_height}")

    def compound_grid(self, compounds, style="w3-center", pose_modal: bool = False):

        id_num_poses_dict = compounds.id_num_poses_dict

        with self.tag("div", klass="w3-row"):
            for compound in compounds:

                with self.tag("div", klass=f"w3-col s12 m6 l4 {style}"):
                    with self.tag("p"):
                        with self.tag("b"):
                            self.text(f"{compound}")

                    self.compound_image(compound)

                    with self.tag("p", klass="w3-small w3-monospace"):
                        self.text(f"{compound.inchikey}")
                        self.doc.asis("<br>")
                        self.text(f"{compound.smiles}")

                    num_poses = id_num_poses_dict[compound.id]
                    self.button(f"{num_poses} poses", disable=num_poses == 0)

                    if pose_modal:
                        poses = compound.poses

                        modal_name = f"modal_c{compound.id}_poses"

                        # MOLECULE MODAL
                        self.modal_button(
                            f"view {len(poses)} poses",
                            modal_name,
                            disable=len(poses) == 0,
                        )

                        if poses:

                            self.save_pset_sdf(modal_name, poses)

                            def modal_content():
                                with self.tag("p"):
                                    self.text("TEXT TEXT TEXT")

                            self.modal(modal_name, modal_content)

    def modal_button(
        self, text, modal_name, disable: bool = False, style: str = "w3-teal"
    ):
        onclick = f"document.getElementById('{modal_name}').style.display='block'"
        self.button(text, style=style, onclick=onclick, disable=disable)

    def button(
        self, text: str, onclick: str = "", style="w3-teal", disable: bool = False
    ):

        assert text

        klass = f"w3-btn {style}"

        if disable:
            klass += " w3-disabled"
            onclick = ""

        with self.tag("button", klass=klass, onclick=onclick):
            self.text(text)

    def modal(self, modal_name, content_function):
        with self.tag("div", id=modal_name, klass="w3-modal"):
            with self.tag("div", klass="w3-modal-content"):
                with self.tag("div", klass="w3-container"):
                    with self.tag(
                        "span",
                        onclick=f"document.getElementById('{modal_name}').style.display='none'",
                        klass="w3-button w3-display-topright",
                    ):
                        self.doc.asis("&times")
                    content_function()

    def recipe_subsection(self, recipe, title, sankey: bool = False):

        self.section_header(title, "h3")

        recipe_name = title.lower().replace(" ", "_")

        if sankey:
            fig = recipe.sankey()
            self.plotly_graph(fig, f"{recipe_name}.html")

        # self.section_header("Products", "h4")

        df = recipe.products.compounds.get_df()

        def modal_content():
            self.table(df, style="w3-table-all w3-small")

        modal_name = f"{recipe_name}_products"
        self.modal_button("products", modal_name)
        self.modal(modal_name, modal_content)

        if intermediates := recipe.intermediates:
            # self.section_header("Intermediates", "h4")
            df = intermediates.compounds.get_df()

            def modal_content():
                self.table(df, style="w3-table-all w3-small")

            modal_name = f"{recipe_name}_intermediates"
            self.modal_button("intermediates", modal_name)
            self.modal(modal_name, modal_content)

        # self.section_header("Reactants", "h4")

        df = recipe.reactants.df

        def modal_content():
            self.table(df, style="w3-table-all w3-small")

        modal_name = f"{recipe_name}_reactants"
        self.modal_button("reactants", modal_name)
        self.modal(modal_name, modal_content)

        # self.section_header("Reactions", "h4")

        df = recipe.reactions.get_df(mols=False)

        def modal_content():
            self.table(df, style="w3-table-all w3-small")

        modal_name = f"{recipe_name}_reactions"
        self.modal_button("reactions", modal_name)
        self.modal(modal_name, modal_content)

    def scorer_attribute(self, attribute, histogram: bool = True):

        from .scoring import DEFAULT_ATTRIBUTES

        key = attribute.key

        self.section_header(f'{attribute._type}: "{key}"')

        with self.tag("ul"):
            self.var("weight", attribute.weight, tag="li")
            self.var("inverse", attribute.inverse, tag="li")
            self.var("min", attribute.min, tag="li")
            self.var("max", attribute.max, tag="li")
            self.var("mean", attribute.mean, tag="li")
            self.var("std", attribute.std, tag="li")

            if key in DEFAULT_ATTRIBUTES:
                self.var(
                    "Description", DEFAULT_ATTRIBUTES[key]["description"], tag="li"
                )

        if histogram:
            fig = attribute.histogram(progress=True)
            self.plotly_graph(fig, f"attribute_{key}_hist.html")

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

    def sec_scaffolds(self) -> None:
        """Section on scaffolds"""

        title = "Scaffolds"
        self.section_header(title)

        self.section_header("All scaffolds", "h3")

        cset = self.animal.compounds(tag="Syndirella base")

        with self.tag("ul"):
            self.var("#compounds", len(cset), tag="li")

        # route dict?
        df = cset.get_df(mol=False, num_poses=True)  # , routes=True)

        def modal_content():
            self.table(df, style="w3-table-all w3-small")

        modal_name = "all_scaffolds_modal"
        self.modal_button(f"all scaffolds table", modal_name)
        self.modal(modal_name, modal_content)

        # quoting?

        self.section_header("Selected scaffolds", "h3")

        cset = self.scaffolds

        self.compound_grid(cset, pose_modal=False)

        # files
        self.section_header("Downloads", "h3")

        table_data = []

        path = self.resource_dir / "all_scaffold_smiles.csv"
        rel_path = Path(self.resource_dir.name) / path.name
        cset.write_smiles_csv(path)
        table_data.append(
            dict(
                Name="All scaffolds",
                Description="CSV of scaffold SMILES",
                Download=f'<a href="{rel_path}" download>CSV</a>',
            )
        )

        path = self.resource_dir / "selected_scaffold_smiles.csv"
        rel_path = Path(self.resource_dir.name) / path.name
        cset.write_smiles_csv(path)
        table_data.append(
            dict(
                Name="Selected scaffolds",
                Description="CSV of scaffold SMILES",
                Download=f'<a href="{rel_path}" download>CSV</a>',
            )
        )

        self.table(table_data)

    def sec_elaborations(self) -> None:
        """Section on elaborations"""

        elabs = self.scaffolds.elabs

        if not elabs:
            logger.warning("No elaborations")
            return None

        self._elaborations = elabs

        title = "Elaborations"
        self.section_header(title)

        fig = self.animal.plot_reaction_funnel(
            title="Syndirella elaboration space", logo=False
        )
        self.plotly_graph(fig, "reaction_funnel.html")

    def sec_quoting(self) -> None:
        """Section on quoting"""

        title = "Quoting"
        self.section_header(title)

    def sec_product_pool(self) -> None:
        """Section on product_pool"""

        title = "product_pool"
        self.section_header(title)

    def sec_route_pool(self) -> None:
        """Section on route_pool"""

        title = "route_pool"
        self.section_header(title)

    def sec_rgen(self) -> None:
        """Section on rgen"""

        title = "Random Recipe Generation"
        self.section_header(title)

        rgen = self.rgen

        with self.tag("ul"):
            self.var("suppliers", rgen.suppliers, tag="li")
            self.var("max_lead_time", rgen.max_lead_time, tag="li")
            self.var("route_pool", len(rgen.route_pool), tag="li")

        self.recipe_subsection(rgen.starting_recipe, "Starting Recipe", sankey=False)

    def sec_scorer(self) -> None:
        """Section on scorer"""

        title = "Recipe Selection"
        self.section_header(title)

        scorer = self.scorer

        with self.tag("ul"):
            self.var("#recipes", len(scorer.recipes), tag="li")
            self.var("#attributes", len(scorer.attributes), tag="li")

        fig = scorer.plot(["price", "score"])
        self.plotly_graph(fig, "scorer_scatter.html")

        for attribute in scorer.attributes:
            self.scorer_attribute(attribute)

    def sec_proposals(self) -> None:
        """Section on proposals"""

        title = "Proposal Recipes"
        self.section_header(title)

        for proposal in self.proposals:

            self.section_header(str(proposal), "h3")

            from .plotting import plot_compound_tsnee

            fig = plot_compound_tsnee(proposal.products.compounds)
            self.plotly_graph(fig, f"proposal_{proposal.hash}.html")

    ### DUNDERS
