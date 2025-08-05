import molparse as mp
import plotly.express as px
import pandas as pd
import mout
import plotly.graph_objects as go

import mrich

"""

	ALL GRAPHS DEFINED HERE SHOULD:

	* Have a HIPPO logo
	* Include the target name in the title

"""

import functools


# hippo_graph decorator
def hippo_graph(func):
    """

    :param func:

    """

    @functools.wraps(func)
    def wrapper(animal, *args, logo="top right", **kwargs):
        """

        :param animal:
        :param *args:
        :param logo:  (Default value = 'top right')
        :param **kwargs:

        """

        wrapper_kwargs = {}
        wrapper_keys = ["show", "html", "pdf", "png"]
        for key in wrapper_keys:
            wrapper_kwargs[key] = kwargs.pop(key, None)

        fig = func(animal, *args, **kwargs)

        if not isinstance(fig, go.Figure):
            return fig

        if wrapper_kwargs["show"]:
            fig.show()

        if wrapper_kwargs["html"]:
            file = wrapper_kwargs["html"]
            if not file.endswith(".html"):
                file = f"{file}.html"
            mp.write(file, fig)

        if wrapper_kwargs["pdf"]:
            file = wrapper_kwargs["pdf"]
            if not file.endswith(".pdf"):
                file = f"{file}.pdf"
            mp.write(file, fig)

        if wrapper_kwargs["png"]:
            file = wrapper_kwargs["png"]
            if not file.endswith(".png"):
                file = f"{file}.png"
            mp.write(file, fig)

        if not fig.layout.images and logo:
            add_hippo_logo(fig, position=logo)

        return fig

    return wrapper


@hippo_graph
def plot_tag_statistics(
    animal,
    color="type",
    subtitle=None,
    log_y=False,
    show_compounds=True,
    show_poses=True,
    compounds=None,
    poses=None,
    title: str | None = None,
    skip: list[str] | None = None,
):
    """

    :param animal:
    :param color:  (Default value = 'type')
    :param subtitle:  (Default value = None)
    :param log_y:  (Default value = False)
    :param compounds:  (Default value = True)
    :param poses:  (Default value = True)

    """

    compounds = compounds or animal.compounds
    poses = poses or animal.poses
    skip = skip or []

    plot_data = []

    for tag in animal.tags.unique:

        if tag in skip:
            continue

        if show_compounds:
            num_compounds = len(compounds.get_by_tag(tag=tag))
            if num_compounds:
                data = dict(tag=tag, number=num_compounds, type="compounds")
                plot_data.append(data)

        if show_poses:
            num_poses = len(poses.get_by_tag(tag=tag))
            if num_poses:
                data = dict(tag=tag, number=num_poses, type="poses")
                plot_data.append(data)

    from pandas import DataFrame

    df = DataFrame(plot_data)

    df.sort_values(by="tag", inplace=True)

    fig = px.bar(df, x="tag", y="number", color=color, log_y=log_y)

    if not title:
        title = "Tag Statistics"

        if subtitle:
            title = f"<b>{animal.name}</b>: {title}<br><sup><i>{subtitle}</i></sup>"
        else:
            title = f"<b>{animal.name}</b>: {title}"

    fig.update_layout(
        title=title, title_automargin=False, title_yref="container", barmode="group"
    )

    fig.update_layout(xaxis_title="Tag", yaxis_title="#")

    return fig


@hippo_graph
def plot_interaction_histogram(
    animal,
    poses,
    feature_metadata,
    subtitle=None,
):
    """

    :param animal:
    :param poses:
    :param feature_metadata:
    :param subtitle:  (Default value = None)

    """

    raise NotImplementedError

    df = animal._fingerprint_df(poses)

    plot_data = []

    for key in df.columns:
        count = int(df[key].sum())

        if not count:
            continue

        data = dict(str=key, count=count)

        data["family"] = feature_metadata[key]["family"]
        data["res_name"] = feature_metadata[key]["res_name"]
        data["res_number"] = feature_metadata[key]["res_number"]
        data["res_chain"] = feature_metadata[key]["res_chain"]
        data["atom_numbers"] = feature_metadata[key]["atom_numbers"]

        data["res_name_number_chain_str"] = (
            f"{feature_metadata[key]['res_name']} {feature_metadata[key]['res_number']} {feature_metadata[key]['res_chain']}"
        )

        plot_data.append(data)

    plot_df = pd.DataFrame(plot_data)
    plot_df.sort_values(["res_chain", "res_number", "family"], inplace=True)
    plot_df

    fig = px.bar(
        plot_df,
        x="res_name_number_chain_str",
        y="count",
        color="family",
        hover_data=plot_df.columns,
    )

    title = "Leveraged protein features"

    if subtitle:
        title = f"<b>{animal.name}</b>: {title}<br><sup><i>{subtitle}</i></sup>"
    else:
        title = f"<b>{animal.name}</b>: {title}"

    fig.update_layout(title=title, title_automargin=False, title_yref="container")

    fig.update_layout(xaxis_title="Residue")
    fig.update_layout(yaxis_title="#Interactions")

    return fig


@hippo_graph
def plot_interaction_punchcard(
    animal,
    poses=None,
    subtitle=None,
    opacity=1.0,
    ignore_chains=False,
):
    """

    :param animal:
    :param poses:  (Default value = None)
    :param subtitle:  (Default value = None)
    :param opacity:  (Default value = 1.0)
    :param group:  (Default value = 'pose_name')
    :param ignore_chains:  (Default value = False)

    """

    import plotly
    from .pset import PoseTable

    poses = poses or animal.poses

    if isinstance(poses, PoseTable):
        iset = animal.interactions
    else:
        iset = poses.interactions

    mrich.var("#poses", len(poses))
    mrich.var("#interactions", len(iset))

    plot_data = iset.df

    name_lookup = poses.id_name_dict

    names = []
    for pose_id in plot_data["pose_id"].values:
        names.append(name_lookup[pose_id])
    plot_data["pose_name"] = names

    if ignore_chains:
        x = "res_name_number"
        plot_data[x] = plot_data[["residue_name", "residue_number"]].agg(
            lambda x: " ".join([str(i) for i in x]), axis=1
        )
        plot_data = plot_data.sort_values([group, x, "residue_number"])
        sort_key = lambda x: x[1]
    else:
        x = "chain_res_name_number_str"
        plot_data[x] = plot_data[["chain_name", "residue_name", "residue_number"]].agg(
            lambda x: " ".join([str(i) for i in x]), axis=1
        )
        sort_key = lambda x: (x[2], x[1])

    title = "Interaction Punch-Card"

    if subtitle:
        title = f"<b>{animal.name}</b>: {title}<br><sup><i>{subtitle}</i></sup>"
    else:
        title = f"<b>{animal.name}</b>: {title}"

    fig = px.scatter(
        plot_data,
        x=x,
        y="type",
        marginal_x="histogram",
        marginal_y="histogram",
        hover_data=plot_data.columns,
        color=group,
        title=title,
    )

    fig.update_layout(title=title, title_automargin=False, title_yref="container")

    fig.update_layout(xaxis_title="Residue", yaxis_title="Feature Family")

    # x-axis sorting
    categoryarray = plot_data[[x, "residue_number", "chain_name"]].agg(tuple, axis=1)
    categoryarray = sorted([v for v in categoryarray.values], key=sort_key)
    categoryarray = [v[0] for v in categoryarray]

    # sort axes
    fig.update_xaxes(categoryorder="array", categoryarray=categoryarray)
    fig.update_yaxes(categoryorder="category descending")

    for trace in fig.data:
        if type(trace) == plotly.graph_objs._histogram.Histogram:
            trace.opacity = 1
            trace.xbins.size = 1
        else:
            trace["marker"]["size"] = 10
            trace["marker"]["opacity"] = opacity

    fig.update_layout(barmode="stack")
    fig.update_layout(scattermode="group", scattergap=0.75)

    return add_punchcard_logo(fig)


# @hippo_graph
def plot_interaction_punchcard_by_tags(
    animal,
    tags: dict[str, str] | list[str],
    permitted_residues: dict[str, list[int]] | None = None,
    yaxis_title: str = "Tag",
    subtitle=None,
    opacity=0.7,
    group="type",
    marginal_histogram_x: bool = True,
    # marginal_histogram_y: bool = False,
    sizeref=0.08,
    counts: bool = True,
    ignore_chains=True,
    backbone_only: bool = False,
    sidechain_only: bool = False,
    return_plot_data: bool = False,
):
    """

    :param animal:
    :param poses:  (Default value = None)
    :param subtitle:  (Default value = None)
    :param opacity:  (Default value = 1.0)
    :param group:  (Default value = 'pose_name')
    :param ignore_chains:  (Default value = False)
    :param permitted_residues: dictionary mapping with keys matching `tags` mapping to list of IDs for residues to include for each tag

    """

    import plotly
    import numpy as np

    if isinstance(tags, list):
        tags = {v: v for v in tags}

    permitted_residues = permitted_residues or {}

    dfs = []
    for group_name, tag in tags.items():
        poses = animal.poses(tag=tag)
        mrich.debug(group_name, poses)

        name_lookup = poses.id_name_dict

        with mrich.loading("Getting interactions dataframe"):
            df = poses.interactions.df

        df["tag"] = tag
        df["group_name"] = group_name
        df["pose_name"] = [name_lookup[pose_id] for pose_id in df["pose_id"].values]

        if group_name in permitted_residues:
            subset = df[df["residue_number"].isin(permitted_residues[group_name])]
            diff = len(df) - len(subset)
            if diff:
                mrich.warning(
                    "Skipping",
                    diff,
                    "markers due unpermitted residue numbers for: ",
                    group_name,
                )
                df = subset

        dfs.append(df)

    mrich.debug("Concatenating dataframes")
    plot_data = pd.concat(dfs, ignore_index=True)

    ### add permitted residues

    if permitted_residues and ignore_chains:

        permitted_df = []
        for group_name, ids in permitted_residues.items():

            unique_combinations = plot_data[plot_data["group_name"] == group_name][
                ["residue_name", "residue_number"]
            ].drop_duplicates()

            for resid in ids:

                try:
                    resname = unique_combinations[
                        unique_combinations["residue_number"] == resid
                    ]["residue_name"].values[0]
                except IndexError:
                    continue

                # if plot_data[['residue_name', 'residue_number']]

                permitted_df.append(
                    dict(res_name_number=f"{resname} {resid}", group_name=group_name)
                )

        permitted_df = pd.DataFrame(permitted_df)

    mrich.debug("Building plot_data")
    if backbone_only:
        plot_data = plot_data[plot_data["backbone"] == True]
    if sidechain_only:
        plot_data = plot_data[plot_data["sidechain"] == True]

    if ignore_chains:
        x = "res_name_number"
        plot_data[x] = plot_data[["residue_name", "residue_number"]].agg(
            lambda x: " ".join([str(i) for i in x]), axis=1
        )
        plot_data = plot_data.sort_values([group, x, "residue_number"])
        sort_key = lambda x: x[1]
    else:
        x = "chain_res_name_number_str"
        plot_data[x] = plot_data[["chain_name", "residue_name", "residue_number"]].agg(
            lambda x: " ".join([str(i) for i in x]), axis=1
        )
        sort_key = lambda x: (x[2], x[1])

    if counts:
        mrich.debug("Summing by residue")
        orig_data = plot_data.copy()

        if ignore_chains:
            plot_data = (
                plot_data.groupby(["group_name", "type", x, "residue_number"])
                .size()
                .reset_index(name="count")
            )
        else:
            plot_data = (
                plot_data.groupby(
                    ["group_name", "type", x, "residue_number", "chain_name"]
                )
                .size()
                .reset_index(name="count")
            )

        plot_data["size"] = np.sqrt(plot_data["count"])

    # add a size reference

    type_str = "Size"
    sizes = [1, 50, 100, 250]

    dicts = []
    for group_name, size in zip(tags.keys(), sizes):
        dicts.append(
            dict(
                group_name=group_name,
                type="type_str",
                res_name_number="",
                residue_number=999,
                count=size,
                size=np.sqrt(size),
                text=size,
            )
        )

    plot_data = pd.concat([plot_data, pd.DataFrame(dicts)])

    mrich.debug("Making scatter plot")
    fig = px.scatter(
        plot_data,
        x=x,
        y="group_name",
        hover_data=plot_data.columns,
        color=group,
        size="size" if counts else None,
        text="text",
        # color_discrete_sequence=px.colors.qualitative.Dark2
    )

    fig.update_traces(textposition="middle right")

    if return_plot_data:
        data_snapshot1 = plot_data.copy()

    # fig.update_layout(title=title, title_automargin=False, title_yref="container")

    fig.update_layout(xaxis_title="Residue", yaxis_title=yaxis_title)

    # x-axis sorting
    if ignore_chains:
        categoryarray = plot_data[[x, "residue_number"]].agg(tuple, axis=1)
    else:
        categoryarray = plot_data[[x, "residue_number", "chain_name"]].agg(
            tuple, axis=1
        )
    categoryarray = sorted([v for v in categoryarray.values], key=sort_key)
    categoryarray = [v[0] for v in categoryarray]

    # sort axes
    fig.update_xaxes(categoryorder="array", categoryarray=categoryarray)

    for trace in fig.data:
        if type(trace) == plotly.graph_objs._histogram.Histogram:
            trace.opacity = 1
            trace.xbins.size = 1
        else:
            trace["marker"]["opacity"] = opacity

    if marginal_histogram_x:

        from plotly.subplots import make_subplots

        subplot_fig = make_subplots(
            rows=2,
            cols=1,
            specs=[[{"type": "histogram"}], [{"type": "scatter"}]],
            shared_xaxes=True,
            shared_yaxes=False,
            vertical_spacing=0.02,
            horizontal_spacing=0.02,
        )

        # add in scatter traces
        for trace in fig.data:
            trace.yaxis = "y2"
            trace.showlegend = False
            trace.marker.sizeref = sizeref
            trace.marker.line.width = 0
            subplot_fig.add_trace(trace)

        # aggregate data for histogram plot
        plot_data = (
            orig_data.groupby(["type", x, "residue_number"])
            .size()
            .reset_index(name="count")
        )

        # generate histogram
        fig2 = px.histogram(plot_data, x=x, y="count", color="type")

        # add in histogram traces
        for trace in fig2.data:
            trace.yaxis = "y1"
            subplot_fig.add_trace(trace)

        # if permitted_residues and ignore_chains:
        #     trace = go.Scatter(
        #         name="Subsite Residues",
        #         x=permitted_df[x],
        #         y=permitted_df["group_name"],
        #         mode="markers",
        #         marker_symbol="square-open",
        #         marker_size=25,
        #         marker_color="grey")

        #     trace.yaxis="y2"
        #     subplot_fig.add_trace(trace)

        # clean up the axes
        subplot_fig.update_layout(
            xaxis=dict(anchor="y2", visible=True, showticklabels=True, side="bottom"),
            # xaxis2=dict(visible=True, showticklabels=True, side="bottom"),
            yaxis2=dict(
                anchor="x",
                overlaying="x",
                side="top",
                categoryorder="category descending",
            ),  # Secondary y-axis for x marginal
        )

        # x-axis sorting
        if ignore_chains:
            categoryarray = plot_data[[x, "residue_number"]].agg(tuple, axis=1)
        else:
            raise NotImplementedError
        categoryarray = sorted([v for v in categoryarray.values], key=sort_key)
        categoryarray = [v[0] for v in categoryarray]
        subplot_fig.update_xaxes(categoryorder="array", categoryarray=categoryarray)

        # y-axis sorting
        categoryarray = list(reversed(tags.keys()))
        subplot_fig.update_yaxes(categoryorder="array", categoryarray=categoryarray)

        # stack histogram bars on top of each other
        subplot_fig.update_layout(barmode="stack")

        subplot_fig.update_layout(
            margin=dict(l=20, r=20, t=20, b=20),
        )

        if return_plot_data:
            return subplot_fig, [orig_data, plot_data, data_snapshot1]

        return subplot_fig

    # fig.update_layout(scattermode="group", scattergap=0.75)

    if return_plot_data:
        return fig, plot_data

    # return add_punchcard_logo(fig)
    return fig


@hippo_graph
def plot_residue_interactions(
    animal, poses, residue_number, subtitle=None, chain=None, target=1
):
    """

    :param animal:
    :param poses:
    :param residue_number:
    :param subtitle:  (Default value = None)
    :param chain:  (Default value = None)

    """

    mrich.var("#poses", len(poses))

    from .iset import InteractionSet

    iset = InteractionSet.from_residue(
        animal.db, residue_number=residue_number, chain=chain
    )

    # return iset

    mrich.var("#interactions", len(iset))

    plot_data = iset.df

    # name_lookup = {i:n for i,n in zip(poses.ids,poses.names)}
    name_lookup = poses.id_name_dict

    # print(name_lookup)

    names = []
    for pose_id in plot_data["pose_id"].values:
        names.append(name_lookup[pose_id])
    plot_data["pose_name"] = names

    fig = px.histogram(plot_data, x="pose_name", color="type")

    # return plot_data[plot_data['pose_name'] == ' x1762b']

    fig.update_xaxes(categoryorder="total descending")

    if not subtitle:
        subtitle = f"#Poses={len(poses)}"

    residue_name = animal.db.get_feature(
        id=plot_data["feature_id"].values[0]
    ).residue_name

    title = f"Interactions w/ {residue_name} {residue_number}"

    if chain:
        title += f" {chain}"

    title = f"<b>{animal.name}</b>: {title}<br><sup><i>{subtitle}</i></sup>"

    fig.update_layout(title=title, title_automargin=False, title_yref="container")

    fig.update_xaxes(title="Pose")
    fig.update_yaxes(title="#Interactions")

    return fig


@hippo_graph
# def plot_building_blocks(animal, subtitle=None, cset='elabs', color='name_is_smiles'):
def plot_reactant_amounts(
    animal, subtitle=None, color="has_price_picker", named_only=False, most_common=None
):
    """

    :param animal:
    :param subtitle:  (Default value = None)
    :param color:  (Default value = 'has_price_picker')
    :param named_only:  (Default value = False)
    :param most_common:  (Default value = None)

    """

    # cset = animal.compound_sets[cset]

    # bbs = cset.get_building_blocks()

    bbs = animal.building_blocks

    mout.debug("making plot_data")
    plot_data = []
    for bb in bbs:
        d = bb.dict
        # if most_common and d['amount'] is not None:

        if not named_only or not d["name_is_smiles"]:
            plot_data.append(d)

    if most_common:
        mout.debug("sorting")
        plot_data = sorted(plot_data, key=lambda x: x["amount"], reverse=True)[
            :most_used_number
        ]

    fig = px.bar(
        plot_data, x="name", y="amount", color=color, hover_data=plot_data[0].keys()
    )
    # fig = px.bar(plot_data, x='smiles', y='amount', color=color)

    title = "Building Blocks"

    if not subtitle:
        # subtitle = f'"{cset.name}": #BBs={len(bbs)}'
        subtitle = f"#BBs={len(bbs)}"

    title = f"<b>{animal.name}</b>: {title}<br><sup><i>{subtitle}</i></sup>"

    fig.update_layout(title=title, title_automargin=False, title_yref="container")

    fig.update_layout(xaxis_title="Reactant", yaxis_title="#Reactions")

    return fig


@hippo_graph
# def plot_building_blocks(animal, subtitle=None, cset='elabs', color='name_is_smiles'):
def plot_reactant_price(animal, subtitle=None, amount=20):
    """

    :param animal:
    :param subtitle:  (Default value = None)
    :param amount:  (Default value = 20)

    """

    # cset = animal.compound_sets[cset]

    # bbs = cset.get_building_blocks()

    bbs = animal.building_blocks

    plot_data = []
    for bb in bbs:
        d = bb.dict
        if not d["has_price_picker"]:
            continue

        d[f"price_{amount}mg"] = bb.get_price(amount)
        d[f"min_amount"] = bb.price_picker.min_amount

        plot_data.append(d)

    # fig = px.bar(plot_data, x='name', y=f'price_{amount}mg', color='lead_time', log_y=True, hover_data=plot_data[0].keys())
    fig = px.histogram(
        plot_data,
        x=f"price_{amount}mg",
        color="lead_time",
        hover_data=plot_data[0].keys(),
    )
    # fig = px.bar(plot_data, x='smiles', y='amount', color=color)

    title = "Reactant Pricing"

    if not subtitle:
        # subtitle = f'"{cset.name}": #BBs={len(bbs)}'
        subtitle = f"#BBs={len(bbs)}"

    title = f"<b>{animal.name}</b>: {title}<br><sup><i>{subtitle}</i></sup>"

    fig.update_layout(title=title, title_automargin=False, title_yref="container")

    fig.update_layout(
        yaxis_title="Number of reactants", xaxis_title=f"Price for {amount}mg [$USD]"
    )

    return fig


@hippo_graph
# def plot_building_blocks(animal, subtitle=None, cset='elabs', color='name_is_smiles'):
def plot_reactants_2d(animal, subtitle=None, amount=20):
    """

    :param animal:
    :param subtitle:  (Default value = None)
    :param amount:  (Default value = 20)

    """

    # cset = animal.compound_sets[cset]

    # bbs = cset.get_building_blocks()

    bbs = animal.building_blocks

    plot_data = []
    for bb in bbs:
        d = bb.dict
        if not d["has_price_picker"]:
            continue

        d[f"price_{amount}mg"] = bb.get_price(amount)
        d[f"min_amount"] = bb.price_picker.min_amount

        plot_data.append(d)

    fig = px.scatter(
        plot_data,
        y="amount",
        x=f"price_{amount}mg",
        color="name",
        hover_data=plot_data[0].keys(),
    )
    # fig = px.bar(plot_data, x='smiles', y='amount', color=color)

    title = "Building Blocks"

    if not subtitle:
        # subtitle = f'"{cset.name}": #BBs={len(bbs)}'
        subtitle = f"#BBs={len(bbs)}"

    title = f"<b>{animal.name}</b>: {title}<br><sup><i>{subtitle}</i></sup>"

    fig.update_layout(title=title, title_automargin=False, title_yref="container")

    fig.update_layout(scattermode="group", scattergap=0.75)

    fig.update_layout(
        yaxis_title="Quantity [mg]", xaxis_title=f"Price for {amount}mg [$USD]"
    )

    return fig


@hippo_graph
# def plot_building_blocks(animal, subtitle=None, cset='elabs', color='name_is_smiles'):
def plot_building_blocks(animal, subtitle=None, color="name_is_smiles"):
    """

    :param animal:
    :param subtitle:  (Default value = None)
    :param color:  (Default value = 'name_is_smiles')

    """

    # cset = animal.compound_sets[cset]

    # bbs = cset.get_building_blocks()

    bbs = animal.building_blocks

    plot_data = []
    for bb in bbs:
        plot_data.append(bb.dict)

    fig = px.scatter(plot_data, x="name", y="max", color="amount")
    # fig = px.bar(plot_data, x='smiles', y='amount', color=color)

    title = "Building Blocks"

    if not subtitle:
        # subtitle = f'"{cset.name}": #BBs={len(bbs)}'
        subtitle = f"#BBs={len(bbs)}"

    title = f"<b>{animal.name}</b>: {title}<br><sup><i>{subtitle}</i></sup>"

    fig.update_layout(title=title, title_automargin=False, title_yref="container")

    fig.update_layout(xaxis_title="Reactant", yaxis_title="Quantity")

    return fig


@hippo_graph
def plot_synthetic_routes(animal, subtitle=None, cset="elabs", color="num_reactants"):
    """

    :param animal:
    :param subtitle:  (Default value = None)
    :param cset:  (Default value = 'elabs')
    :param color:  (Default value = 'num_reactants')

    """

    cset = animal.compound_sets[cset]

    plot_data = []
    for reax in cset.reactions:
        plot_data.append(reax.dict)

    # fig = px.bar(plot_data, x='name', y='amount', color=color)
    fig = px.histogram(plot_data, x="type", color=color)

    title = "Synthetic Routes"

    if not subtitle:
        subtitle = f'"{cset.name}": #compounds={len(cset)}'

    title = f"<b>{animal.name}</b>: {title}<br><sup><i>{subtitle}</i></sup>"

    fig.update_layout(title=title, title_automargin=False, title_yref="container")

    fig.update_layout(xaxis_title="Compound", yaxis_title="#Routes")

    return fig


@hippo_graph
def plot_numbers(animal, subtitle=None):
    """
            - y-axis: numbers
            - x-categories
                    * hits
                    * hit poses
                    * bases
                    * base poses
                    * elabs
                    * elab poses
                    * BBs (total)
                    * BBs (in enamine)

    :param animal:
    :param subtitle:  (Default value = None)

    """

    # cset = animal.compound_sets[cset]

    plot_data = [
        dict(category="Experimental Hits", number=len(animal.hits), type="compound"),
        dict(category="Experimental Hits", number=len(animal.hits.poses), type="poses"),
        dict(category="Base compounds", number=len(animal.bases), type="compound"),
        dict(category="Base compounds", number=len(animal.bases.poses), type="poses"),
        dict(
            category="Syndirella Elaborations",
            number=len(animal.elabs),
            type="compound",
        ),
        dict(
            category="Syndirella Elaborations",
            number=len(animal.elabs.poses),
            type="poses",
        ),
        dict(
            category="Unique Reactants",
            number=len(animal.building_blocks),
            type="compound",
        ),
    ]

    fig = px.bar(plot_data, x="category", y="number", log_y=True, color="type")

    title = "Compounds & Poses"

    title = f"<b>{animal.name}</b>: {title}<br>"

    fig.update_layout(
        title=title, title_automargin=False, title_yref="container", barmode="group"
    )

    fig.update_layout(xaxis_title=None, yaxis_title="Log(Quantity)")

    return fig


@hippo_graph
def plot_compound_property(
    animal,
    prop,
    compounds=None,
    style="bar",
    null=None,
    hover_data=None,
    custom_data=None,
):
    """Get an arbitrary property from all the compounds in animal.compounds

            If one property, plot a 1D histogram
            If 2D plot a bar/scatter

    :param animal:
    :param prop:
    :param compounds:  (Default value = None)
    :param style:  (Default value = 'bar')
    :param null:  (Default value = None)

    """

    if not isinstance(prop, list):
        prop = [prop]

    if hover_data is None:
        hover_data = []

    plot_data = []

    if not compounds:
        compounds = animal.compounds

    if len(compounds) > 1000:
        compounds = mrich.track(compounds, prefix="Generating plot data")

    for comp in compounds:

        data = comp.dict

        for p in prop:

            # has attr

            if p not in data:

                # get attr
                if hasattr(comp, p):
                    v = getattr(comp, p)

                elif p in (m := comp.metadata):
                    v = m[p]

                else:
                    v = null

                data[p] = v

        plot_data.append(data)

    if len(prop) == 1:

        title = f"Compound {prop[0]}"

        fig = px.histogram(plot_data, x=prop[0])

        fig.update_layout(xaxis_title=prop[0], yaxis_title="Quantity")

    elif len(prop) == 2:

        hover_data = prop + ["smiles"] + hover_data

        title = f"Compound {prop[0]} vs {prop[1]}"

        func = eval(f"px.{style}")
        fig = func(
            plot_data,
            x=prop[0],
            y=prop[1],
            hover_data=hover_data,
            custom_data=custom_data,
        )

        fig.update_layout(xaxis_title=prop[0], yaxis_title=prop[1])

    else:
        mout.error("Unsupported", code="plotting.plot_compound_property.1")

    title = f"<b>{animal.name}</b>: {title}<br>"

    fig.update_layout(
        title=title, title_automargin=False, title_yref="container", barmode="group"
    )

    return fig


@hippo_graph
def plot_pose_property(
    animal,
    prop,
    poses=None,
    style="scatter",
    title=None,
    null=None,
    color=None,
    log_y=False,
    subtitle=None,
    data_only=False,
    custom_data=None,
    **kwargs,
):
    """Get an arbitrary property from all the poses in animal.poses

            If one property, plot a 1D histogram
            If 2D plot a scatter plot

    :param animal:
    :param prop:
    :param poses:  (Default value = None)
    :param style:  (Default value = 'scatter')
    :param title:  (Default value = None)
    :param null:  (Default value = None)
    :param color:  (Default value = None)
    :param log_y:  (Default value = False)
    :param subtitle:  (Default value = None)
    :param data_only:  (Default value = False)
    :param **kwargs:

    """

    # fetch these directly from the database
    if prop in ["energy_score", "distance_score"]:

        p = prop

        if not poses:
            # great all poses!
            poses = animal.poses
            n_poses = len(poses)
            mrich.out(f"Querying database for {n_poses} poses...")
            field = f"pose_{p}"
            title = title or f"{p} of all poses"
            query = animal.db.select_where(
                table="pose", query=field, key=f"{field} is not NULL", multiple=True
            )

        else:

            # subset of poses
            assert poses.table == "pose", f"{poses=} is not a set of Pose objects"
            n_poses = len(poses)
            mrich.out(f"Querying database for {n_poses} poses...")
            field = f"pose_{p}"
            title = title or f"{p} of pose subset"
            query = animal.db.select_where(
                table="pose",
                query=field,
                key=f"{field} is not NULL and pose_id in {poses.str_ids}",
                multiple=True,
            )

        plot_data = [{p: v} for v, in query]

        if p == "energy_score":
            subtitle = (
                subtitle
                or f"#poses={n_poses}, energy_score < 0 = {len([None for d in plot_data if d[p] < 0])/n_poses:.1%}"
            )
        elif p == "distance_score":
            subtitle = (
                subtitle
                or f"#poses={n_poses}, distance_score < 2 = {len([None for d in plot_data if d[p] < 2])/n_poses:.1%}"
            )
        else:
            subtitle = subtitle or f"#poses={n_poses}"

        prop = [prop]

    elif prop == ["energy_score", "distance_score"] or prop == [
        "distance_score",
        "energy_score",
    ]:

        query = f"pose_id, pose_distance_score, pose_energy_score"

        # hardcoded errorbars
        distance_score_err = 0.03
        energy_score_err = 6

        if color:
            query += f", {color}"

        if not poses:
            # great all poses!
            poses = animal.poses
            n_poses = len(poses)
            mrich.out(f"Querying database for {n_poses} poses...")
            title = f"distance & energy scores of all poses"
            query = animal.db.select(table="pose", query=query, multiple=True)

        else:

            # subset of poses
            n_poses = len(poses)
            mrich.out(f"Querying database for {n_poses} poses...")
            title = f"distance & energy scores of pose subset"
            query = animal.db.select_where(
                table="pose",
                query=query,
                key=f"pose_id in {poses.str_ids}",
                multiple=True,
            )

        plot_data = []
        for q in query:
            d = {
                "id": q[0],
                "distance_score": q[1],
                "energy_score": q[2],
                "distance_score_err": distance_score_err,
                "energy_score_err": energy_score_err,
            }
            if color:
                d[color] = q[-1]
                if color == "pose_compound":
                    d[color] = f"C{d[color]}"

            plot_data.append(d)

        kwargs["error_x"] = "energy_score_err"
        kwargs["error_y"] = "distance_score_err"

        subtitle = subtitle or f"#poses={n_poses}"

    # elif prop == ['num_atoms_added', 'energy_score'] or prop == ['num_atoms_added', 'energy_score']:
    # 	mrich.error('Use animal.plot_pose_risk_vs_placement')
    # 	raise NotImplementedError

    # 	query = f'pose_id, , pose_distance_score, pose_energy_score'

    # 	if not poses:
    # 		# great all poses!
    # 		poses = animal.poses
    # 		n_poses = len(poses)
    # 		mrich.out(f'Querying database for {n_poses} poses...')
    # 		title = f'distance & energy scores of all poses'
    # 		query = animal.db.select(table='pose', query=query, multiple=True)

    # 	else:

    # 		# subset of poses
    # 		n_poses = len(poses)
    # 		mrich.out(f'Querying database for {n_poses} poses...')
    # 		title = f'distance & energy scores of pose subset'
    # 		query = animal.db.select_where(table='pose', query=query, key=f'pose_id in {poses.str_ids}', multiple=True)

    # 	plot_data = [{'id':id, 'distance_score':v1, 'energy_score':v2 } for id,v1,v2 in query]

    # 	subtitle = f'#poses={n_poses}'

    else:

        if not poses:
            poses = animal.poses

        if prop == "tags":

            plot_data = []

            for tag in poses.tags:

                num_poses = len(poses.get_by_tag(tag=tag))
                data = dict(tag=tag, number=num_poses)
                plot_data.append(data)

            fig = px.bar(plot_data, x="tag", y="number", color=color, log_y=log_y)

            title = "Tag Statistics"

            fig.update_layout(
                title=title,
                title_automargin=False,
                title_yref="container",
                barmode="group",
            )

            fig.update_layout(xaxis_title="Tag", yaxis_title="#")

            return fig

        if not isinstance(prop, list):
            prop = [prop]

        plot_data = []
        if len(poses) > 1000:
            poses = mrich.track(poses, prefix="Generating plot data")

        for pose in poses:
            if len(prop) > 1:
                data = dict(id=pose.id)
            else:
                data = {}

            for p in prop:
                if p not in data:
                    # get attr
                    if hasattr(pose, p):
                        v = getattr(pose, p)

                    elif p in (m := pose.metadata):
                        v = m[p]

                    else:
                        v = null

                    data[p] = v

            if color:
                if color not in data:
                    # get attr
                    if hasattr(pose, color):
                        v = getattr(pose, color)

                    elif color in (m := pose.metadata):
                        v = m[color]

                    else:
                        v = null

                    data[color] = v

            plot_data.append(data)

    if data_only:
        return plot_data

    hover_data = ["id"]  # , 'alias', 'inchikey'] #, 'tags', 'inspirations']

    if len(prop) == 1:

        title = title or f"Pose {prop[0]}"

        fig = px.histogram(plot_data, x=prop[0], hover_data=None, color=color, **kwargs)

        fig.update_layout(xaxis_title=prop[0], yaxis_title="Quantity")

    elif len(prop) == 2:

        if style == "histogram":

            x = [d[prop[0]] for d in plot_data]
            y = [d[prop[1]] for d in plot_data]

            fig = go.Figure(go.Histogram2d(x=x, y=y, **kwargs))

        else:

            # if style == "bar":
            # style = "scatter"

            func = eval(f"px.{style}")
            fig = func(
                plot_data,
                x=prop[0],
                y=prop[1],
                color=color,
                hover_data=hover_data,
                custom_data=custom_data,
                **kwargs,
            )

        title = title or f"Pose {prop[0]} vs {prop[1]}"
        fig.update_layout(xaxis_title=prop[0], yaxis_title=prop[1])

    else:
        mout.error("Unsupported", code="plotting.plot_pose_property.1")

    title = title or f"<b>{animal.name}</b>: {title}<br>"

    if subtitle:
        title = f"{title}<br><sup><i>{subtitle}</i></sup>"

    fig.update_layout(
        title=title, title_automargin=False, title_yref="container", barmode="group"
    )

    return fig


@hippo_graph
def plot_compound_availability(animal, compounds=None, title=None, subtitle=None):
    """

    :param animal:
    :param compounds:  (Default value = None)
    :param title:  (Default value = None)
    :param subtitle:  (Default value = None)

    """

    from .cset import CompoundTable, CompoundSet

    compounds = compounds or animal.compounds

    match compounds:
        case CompoundTable():
            pairs = animal.db.select(
                table="quote",
                query="DISTINCT quote_supplier, quote_catalogue",
                multiple=True,
            )

            plot_data = []
            for supplier, catalogue in pairs:

                if catalogue is None:
                    catalogue = "None"
                    cat_str = "NULL"
                else:
                    cat_str = f'"{catalogue}"'

                (count,) = animal.db.select_where(
                    table="quote",
                    query="COUNT(DISTINCT quote_compound)",
                    key=f'quote_supplier IS "{supplier}" AND quote_catalogue IS {cat_str}',
                )

                plot_data.append(
                    dict(supplier=supplier, catalogue=catalogue, count=count)
                )

        case CompoundSet():

            pairs = animal.db.select(
                table="quote",
                query="DISTINCT quote_supplier, quote_catalogue",
                multiple=True,
            )

            plot_data = []
            for supplier, catalogue in pairs:

                if catalogue is None:
                    catalogue = "None"
                    cat_str = "NULL"
                else:
                    cat_str = f'"{catalogue}"'

                (count,) = animal.db.select_where(
                    table="quote",
                    query="COUNT(DISTINCT quote_compound)",
                    key=f'quote_supplier IS "{supplier}" AND quote_catalogue IS {cat_str} AND quote_compound IN {compounds.str_ids}',
                )

                if not count:
                    continue

                plot_data.append(
                    dict(supplier=supplier, catalogue=catalogue, count=count)
                )

        case _:
            raise NotImplementedError

    fig = px.bar(plot_data, x="catalogue", y="count", color="supplier")

    title = "Compound availability"

    title = title or f"<b>{animal.name}</b>: Compound availability<br>"

    if subtitle:
        title = f"{title}<br><sup><i>{subtitle}</i></sup>"

    fig.update_layout(
        title=title
    )  # ,title_automargin=False, title_yref='container', barmode='group')

    return fig


# @hippo_graph
def plot_compound_availability_venn(animal, compounds):
    """

    :param animal:
    :param compounds:  (Default value = None)
    :param title:  (Default value = None)
    :param subtitle:  (Default value = None)

    """

    from venn import venn
    from .cset import CompoundTable, CompoundSet

    pairs = animal.db.select(
        table="quote",
        query="DISTINCT quote_supplier, quote_catalogue",
        multiple=True,
    )

    plot_data = {}

    for supplier, catalogue in pairs:

        if catalogue is None:
            catalogue = "None"
            cat_str = "NULL"
        else:
            cat_str = f'"{catalogue}"'

        if (supplier, catalogue) not in plot_data:
            plot_data[(supplier, catalogue)] = set()

        records = animal.db.select_where(
            table="quote",
            query="quote_compound",
            key=f'quote_supplier IS "{supplier}" AND quote_catalogue IS {cat_str} AND quote_compound IN {compounds.str_ids}',
            multiple=True,
            none="quiet",
        )

        if not records:
            continue

        for (i,) in records:
            plot_data[(supplier, catalogue)].add(i)

    plot_data = {k: v for k, v in plot_data.items() if v}

    # return plot_data

    return venn(plot_data)


@hippo_graph
def plot_compound_price(
    animal,
    compounds=None,
    min_amount=1,
    subtitle=None,
    title=None,
    style="histogram",
    **kwargs,
):
    """

    :param animal:
    :param compounds:  (Default value = None)
    :param min_amount:  (Default value = 1)
    :param subtitle:  (Default value = None)
    :param title:  (Default value = None)
    :param style:  (Default value = 'histogram')
    :param **kwargs:

    """

    from .cset import CompoundTable, CompoundSet
    import numpy as np

    compounds = compounds or animal.compounds

    match compounds:
        case CompoundTable():

            if style == "scatter":

                sql = """
				SELECT quote_compound, quote_amount, MIN(quote_price), quote_lead_time, compound_smiles, COUNT(DISTINCT reactant_reaction)
				FROM quote 
				INNER JOIN compound ON quote.quote_compound = compound.compound_id
				INNER JOIN reactant ON quote.quote_compound = reactant.reactant_compound
				WHERE quote_amount >= {min_amount}
				GROUP BY quote_compound
				""".format(
                    min_amount=min_amount
                )

                results = animal.db.execute(sql).fetchall()

                n_compounds = len(results)

                plot_data = []
                for (
                    compound_id,
                    amount,
                    price,
                    lead_time,
                    smiles,
                    num_reactions,
                ) in results:
                    plot_data.append(
                        dict(
                            compound_id=compound_id,
                            min_price=price,
                            quoted_amount=amount,
                            lead_time=lead_time,
                            smiles=smiles,
                            num_reactions=num_reactions,
                            log_price_per_reaction=np.log(price / num_reactions),
                            price_per_reaction=price / num_reactions,
                        )
                    )

            else:
                data = animal.db.select_where(
                    table="quote",
                    query="quote_amount, MIN(quote_price)",
                    key=f"quote_amount >= {min_amount} GROUP BY quote_compound",
                    multiple=True,
                )

                n_compounds = len(data)

                plot_data = []
                for amount, price in data:
                    plot_data.append(dict(min_price=price, quoted_amount=amount))

        case CompoundSet():

            if style == "scatter":

                sql = """
				SELECT quote_compound, quote_amount, MIN(quote_price), quote_lead_time, compound_smiles, COUNT(DISTINCT reactant_reaction)
				FROM quote 
				INNER JOIN compound ON quote.quote_compound = compound.compound_id
				INNER JOIN reactant ON quote.quote_compound = reactant.reactant_compound
				WHERE quote_amount >= {min_amount}
				AND quote_compound IN {str_ids}
				GROUP BY quote_compound
				""".format(
                    min_amount=min_amount, str_ids=compounds.str_ids
                )

                results = animal.db.execute(sql).fetchall()

                n_compounds = len(results)

                plot_data = []
                for (
                    compound_id,
                    amount,
                    price,
                    lead_time,
                    smiles,
                    num_reactions,
                ) in results:
                    plot_data.append(
                        dict(
                            compound_id=compound_id,
                            min_price=price,
                            quoted_amount=amount,
                            lead_time=lead_time,
                            smiles=smiles,
                            num_reactions=num_reactions,
                            log_price_per_reaction=np.log(price / num_reactions),
                            price_per_reaction=price / num_reactions,
                        )
                    )

            else:
                data = animal.db.select_where(
                    table="quote",
                    query="quote_amount, MIN(quote_price)",
                    key=f"quote_amount >= {min_amount} AND quote_compound IN {compounds.str_ids} GROUP BY quote_compound",
                    multiple=True,
                )

                n_compounds = len(data)

                plot_data = []
                for amount, price in data:
                    plot_data.append(dict(min_price=price, quoted_amount=amount))

        case _:
            raise NotImplementedError("CompoundSet not yet supported")

    plot_data = sorted(plot_data, key=lambda x: x["quoted_amount"])

    match style:

        case "histogram":
            fig = px.histogram(
                plot_data, color="quoted_amount", x="min_price", **kwargs
            )

        case "violin":
            fig = px.violin(plot_data, color="quoted_amount", x="min_price", **kwargs)

        case "scatter":
            fig = px.scatter(
                plot_data,
                color="log_price_per_reaction",
                x="min_price",
                y="lead_time",
                hover_data=plot_data[0].keys(),
                **kwargs,
            )

        case _:
            raise NotImplementedError(f"{style=}")

    subtitle = subtitle or f"#compounds={n_compounds}, {min_amount=} mg"

    title = title or f"<b>{animal.name}</b>: Compound price<br>"

    if subtitle:
        title = f"{title}<br><sup><i>{subtitle}</i></sup>"

    fig.update_layout(
        title=title
    )  # ,title_automargin=False, title_yref='container', barmode='group')

    return fig


@hippo_graph
def plot_reaction_funnel(animal, title=None, subtitle=None):
    """

    :param animal:
    :param title:  (Default value = None)
    :param subtitle:  (Default value = None)

    """

    compounds = animal.compounds

    data = dict(
        number=[
            compounds.num_reactants,
            compounds.num_intermediates,
            compounds.num_products,
        ],
        category=["Reactants", "Intermediates", "Products"],
    )

    fig = px.funnel(data, x="category", y="number")

    title = title or f"<b>{animal.name}</b>: Reaction statistics"

    if subtitle:
        title = f"{title}<br><sup><i>{subtitle}</i></sup>"

    fig.update_layout(title=title, title_automargin=False, title_yref="container")

    return fig


HIPPO_LOGO_URL = "https://raw.githubusercontent.com/mwinokan/HIPPO/main/logos/hippo_logo_tightcrop.png"
HIPPO_HEAD_URL = (
    "https://raw.githubusercontent.com/mwinokan/HIPPO/main/logos/hippo_assets-02.png"
)


def plot_pose_interactions(
    animal: "HIPPO", pose: "Pose"
) -> "plotly.graph_objects.Figure":
    """3d figure showing the interactions between a :class:`.Pose` and the protein. In a Jupyter notebook this figure may be unusable, instead write it as a HTML file and open it in your browser:

    ::

            import molparse as mp
            fig = animal.plot_pose_interactions(pose)
            mp.write(f'{pose}_interactions.html', fig)

    :param pose: the :class:`.Pose` whose interactions are to be rendered

    """

    import molparse as mp

    # get interactions
    iset = pose.interactions

    # get the protein
    protein = pose.protein_system

    # get interacting residues
    pairs = iset.residue_number_chain_pairs

    # get residues
    residues = []
    for resnum, chain in pairs:
        residues.append(protein.get_chain(chain).residues[f"n{resnum}"])

    # get ligand
    lig_group = mp.rdkit.mol_to_AtomGroup(pose.mol)

    # create the combined plotting group
    plot_group = mp.AtomGroup.from_any(str(pose), residues + [lig_group])

    # interaction labels and vectors
    labels = []
    extras = []
    for interaction in iset:
        extras.append([interaction.prot_coord, interaction.lig_coord])
        labels.append(interaction.description)

    # create the figure
    fig = plot_group.plot3d(show=False, extra=extras, extra_labels=labels)

    return fig


@hippo_graph
def plot_compound_tsnee(
    compounds,
    df: "pd.DataFrame | None" = None,
    title: str | None = None,
    subtitle: str | None = None,
    legend: bool = False,
    color="cluster",
    symbol: str = "type",
    sort_by: str = "type",
    color: str = "cluster",
    **kwargs,
) -> "plotly.graph_objects.Figure":

    from .pca import get_cfps
    from sklearn.decomposition import PCA
    import numpy as np

    mrich.var("#compounds", len(compounds))

    if df is None:
        
        with mrich.loading("Getting Compound DataFrame"):
            df = compounds.get_df(mol=True, bases=True)
    
        def get_cluster(row):
    
            bases = row["bases"]
    
            if not bases:
                return row["id"]
    
            if len(bases) == 1:
                return list(bases)[0]
    
            return tuple(bases)
    
        def get_type(row):
    
            if row["bases"] is None:
                return "base"
    
            return "elaboration"
    
        with mrich.loading("Adding columns"):
            df["cluster"] = df.apply(get_cluster, axis=1)
            df["type"] = df.apply(get_type, axis=1)

    with mrich.loading("Getting Compound fingerprints"):
        df["FP"] = df["mol"].map(get_cfps)

<<<<<<< HEAD
=======
    df["bases"] = df["bases"].map(lambda x: x if not isinstance(x, float) else None)

    def get_cluster(row):

        bases = row["bases"]

        if not bases:
            return row["id"]

        if len(bases) == 1:
            return list(bases)[0]

        return tuple(bases)

    def get_type(row):

        if row["bases"] is None:
            return "scaffold"

        return "elaboration"

    with mrich.loading("Adding columns"):
        df["cluster"] = df.apply(get_cluster, axis=1)
        df["type"] = df.apply(get_type, axis=1)

    if sort_by:
        df = df.sort_values(by=sort_by)

>>>>>>> 325bedd8ef361af8c7441d54c4ad2c97b241a8e3
    X = np.array([x.fp for x in df["FP"]])

    with mrich.loading("Computing PCA"):
        pca = PCA(n_components=2, random_state=0)
        pca_fit = pca.fit_transform(X)

    df["PC1"] = pca_fit.T[0]
    df["PC2"] = pca_fit.T[1]

    hover_data = [
        "id",
        "smiles",
        # "alias",
        "inchikey",
        "PC1",
        "PC2",
        "bases",
        "cluster",
        "type",
    ]

    df["bases"] = df["bases"].astype(str)
    df["cluster"] = df["cluster"].astype(str)

    # return df

    # mrich.debug(df.columns)

    with mrich.loading("Creating figure"):
        fig = px.scatter(
            df,
            x="PC1",
            y="PC2",
            hover_data=hover_data,
            color=color,
            symbol=symbol,
            **kwargs,
        )

    subtitle = subtitle or f"#compounds={len(compounds)}"

    title = title or f"{compounds} PCA<br>"

    if subtitle:
        title = f"{title}<br><sup><i>{subtitle}</i></sup>"

    fig.update_layout(title=title)

    if not legend:
        fig.update_layout(showlegend=False)

    return fig


def add_hippo_logo(fig, in_plot=True, position="top right"):
    """

    :param fig:
    :param in_plot:  (Default value = True)
    :param position:  (Default value = 'top right')

    """

    assert fig.layout.title.text, "Figure must have a title to add the HIPPO logo"

    if in_plot:

        sizex = 0.3
        sizey = 0.3

        if "top" in position:
            yanchor = "top"
            y = 0.95
        elif "bottom" in position:
            yanchor = "bottom"
            y = 0.05
        else:
            yanchor = "middle"
            y = 0.50

        if "left" in position:
            xanchor = "left"
            x = 0.05
        elif "right" in position:
            xanchor = "right"
            x = 0.95
        else:
            xanchor = "center"
            x = 0.50

        fig.add_layout_image(
            dict(
                source=HIPPO_LOGO_URL,
                xref="paper",
                yref="paper",
                # layer='below',
                x=x,
                y=y,
                sizex=sizex,
                sizey=sizey,
                xanchor=xanchor,
                yanchor=yanchor,
            )
        )

        return fig

    has_legend = fig.layout.legend.title.text is not None
    fig.layout.margin.t = None

    if has_legend:

        fig.add_layout_image(
            dict(
                source="",
                xref="paper",
                yref="paper",
                x=1,
                y=1.05,
                sizex=0.4,
                sizey=0.4,
                xanchor="left",
                yanchor="bottom",
            )
        )

    else:

        fig.add_layout_image(
            dict(
                source=HIPPO_LOGO_URL,
                xref="paper",
                yref="paper",
                x=1,
                y=1.05,
                sizex=0.3,
                sizey=0.3,
                xanchor="right",
                yanchor="bottom",
            )
        )

    return fig


def add_punchcard_logo(fig):
    """Add the HIPPO logo to a punchcard figure"""

    fig.add_layout_image(
        dict(
            source=HIPPO_HEAD_URL,
            xref="paper",
            yref="paper",
            x=1,
            y=1,
            sizex=0.25,
            sizey=0.25,
            xanchor="right",
            yanchor="top",
        )
    )

    return fig
