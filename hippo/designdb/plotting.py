"""Plotting helpers (plotly figures) operating on sets."""

import mrich

HIPPO_HEAD_URL = (
    'https://raw.githubusercontent.com/mwinokan/HIPPO/main/logos/hippo_assets-02.png'
)


def add_punchcard_logo(fig):
    """Add the HIPPO logo to a punch-card figure."""
    fig.add_layout_image(
        dict(
            source=HIPPO_HEAD_URL,
            xref='paper',
            yref='paper',
            x=1,
            y=1,
            sizex=0.25,
            sizey=0.25,
            xanchor='right',
            yanchor='top',
        )
    )
    return fig


def plot_interaction_punchcard(
    poses,
    *,
    title_prefix: str | None = None,
    subtitle: str | None = None,
    opacity: float = 1.0,
    group: str = 'pose_name',
    ignore_chains: bool = False,
):
    """Interaction punch-card (residue vs. feature family) for a :class:`.PoseSet`.

    :param poses: the :class:`.PoseSet` to plot
    :param title_prefix: bold prefix for the title (e.g. the target name)
    :param subtitle: optional subtitle
    :param opacity: marker opacity
    :param group: column to colour points by (default: per-pose)
    :param ignore_chains: drop the chain from the residue axis label
    """
    import plotly.express as px
    import plotly.graph_objects as go

    iset = poses.interactions
    mrich.var('#poses', len(poses))
    mrich.var('#interactions', len(iset))

    plot_data = iset.df
    if plot_data.empty:
        mrich.warning('No interactions to plot')
        return None

    name_lookup = poses.id_name_dict
    plot_data['pose_name'] = [name_lookup.get(i) for i in plot_data['pose_id'].values]

    if ignore_chains:
        x = 'res_name_number'
        plot_data[x] = plot_data[['residue_name', 'residue_number']].agg(
            lambda r: ' '.join(str(i) for i in r), axis=1
        )
        sort_key = lambda v: v[1]
    else:
        x = 'chain_res_name_number_str'
        plot_data[x] = plot_data[['chain_name', 'residue_name', 'residue_number']].agg(
            lambda r: ' '.join(str(i) for i in r), axis=1
        )
        sort_key = lambda v: (v[2], v[1])

    title = 'Interaction Punch-Card'
    if title_prefix:
        title = f'<b>{title_prefix}</b>: {title}'
    if subtitle:
        title += f'<br><sup><i>{subtitle}</i></sup>'

    fig = px.scatter(
        plot_data,
        x=x,
        y='type',
        marginal_x='histogram',
        marginal_y='histogram',
        hover_data=plot_data.columns,
        color=group,
        title=title,
    )

    fig.update_layout(title=title, title_automargin=False, title_yref='container')
    fig.update_layout(xaxis_title='Residue', yaxis_title='Feature Family')

    # order the residue axis by (residue_number, chain)
    categoryarray = plot_data[[x, 'residue_number', 'chain_name']].agg(tuple, axis=1)
    categoryarray = [v[0] for v in sorted(categoryarray.values, key=sort_key)]
    fig.update_xaxes(categoryorder='array', categoryarray=categoryarray)
    fig.update_yaxes(categoryorder='category descending')

    for trace in fig.data:
        if isinstance(trace, go.Histogram):
            trace.opacity = 1
            trace.xbins.size = 1
        else:
            trace['marker']['size'] = 10
            trace['marker']['opacity'] = opacity

    fig.update_layout(barmode='stack')
    fig.update_layout(scattermode='group', scattergap=0.75)

    return add_punchcard_logo(fig)
