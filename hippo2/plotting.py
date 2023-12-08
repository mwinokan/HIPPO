
import plotly.express as px
import pandas as pd

'''

	ALL GRAPHS DEFINED HERE SHOULD:

	* Have a HIPPO logo
	* Include the target name in the title

'''

def plot_tag_statistics(animal, color=None, subtitle=None):

	plot_data = []

	for tag in animal.all_tags:

		num_poses = len(animal.get_poses(tag))

		data = dict(tag=tag, num_poses=num_poses)

		plot_data.append(data)

	fig = px.bar(plot_data, x='tag', y='num_poses', color=None)

	title = 'Tag Statistics'

	if subtitle:
		title = f'<b>{animal.target_name}</b>: {title}<br><sup><i>{subtitle}</i></sup>'
	else:
		title = f'<b>{animal.target_name}</b>: {title}'
	
	fig.update_layout(title=title,title_automargin=False, title_yref='container')

	fig.update_layout(xaxis_title='Tag', yaxis_title='#Poses')

	return add_hippo_logo(fig)

def plot_interaction_histogram(animal, poses, feature_metadata, subtitle=None,):

	df = animal._fingerprint_df(poses)

	plot_data = []
	
	for key in df.columns:
		count = int(df[key].sum())

		if not count:
			continue
		
		data = dict(str=key,count=count)

		data['family'] = feature_metadata[key]['family']
		data['res_name'] = feature_metadata[key]['res_name']
		data['res_number'] = feature_metadata[key]['res_number']
		data['res_chain'] = feature_metadata[key]['res_chain']
		data['atom_numbers'] = feature_metadata[key]['atom_numbers']

		data['res_name_number_chain_str'] = f"{feature_metadata[key]['res_name']} {feature_metadata[key]['res_number']} {feature_metadata[key]['res_chain']}"
	
		plot_data.append(data)

	plot_df = pd.DataFrame(plot_data)
	plot_df.sort_values(['res_chain', 'res_number', 'family'], inplace=True)
	plot_df

	fig = px.bar(plot_df, x='res_name_number_chain_str', y='count', color='family', hover_data=plot_df.columns)

	title='Leveraged protein features'

	if subtitle:
		title = f'<b>{animal.target_name}</b>: {title}<br><sup><i>{subtitle}</i></sup>'
	else:
		title = f'<b>{animal.target_name}</b>: {title}'
	
	fig.update_layout(title=title,title_automargin=False, title_yref='container')

	fig.update_layout(xaxis_title='Residue')
	fig.update_layout(yaxis_title='#Interactions')

	return add_hippo_logo(fig)

def plot_interaction_punchcard(animal, poses, feature_metadata, subtitle=None, opacity=1.0):

	import plotly

	plot_data = []

	for pose in poses:

		fingerprint = pose._fingerprint

		# loop over each interaction in the pose
		for key, value in fingerprint.items():
			if not value:
				continue
			
			data = dict(str=key,count=value)

			data['family'] = feature_metadata[key]['family']
			data['res_name'] = feature_metadata[key]['res_name']
			data['res_number'] = feature_metadata[key]['res_number']
			data['res_chain'] = feature_metadata[key]['res_chain']
			data['atom_numbers'] = feature_metadata[key]['atom_numbers'] 
			data['pose'] = pose.longname 

			data['res_name_number_chain_str'] = f"{feature_metadata[key]['res_name']} {feature_metadata[key]['res_number']} {feature_metadata[key]['res_chain']}"

			plot_data.append(data)

	plot_df = pd.DataFrame(plot_data)

	title = 'Interaction Punch-Card'

	if subtitle:
		title = f'<b>{animal.target_name}</b>: {title}<br><sup><i>{subtitle}</i></sup>'
	else:
		title = f'<b>{animal.target_name}</b>: {title}'

	fig = px.scatter(plot_df, x='res_name_number_chain_str', y='family',marginal_x='histogram',marginal_y='histogram', hover_data=plot_df.columns, color='pose', title=title)
	
	fig.update_layout(title=title,title_automargin=False, title_yref='container')

	fig.update_layout(xaxis_title='Residue', yaxis_title='Feature Family')

	for trace in fig.data:
		if type(trace) == plotly.graph_objs._histogram.Histogram:
			trace.opacity = 1
			trace.xbins.size = 1
		else:
			trace['marker']['size'] = 10
			trace['marker']['opacity'] = opacity

	fig.update_layout(barmode='stack')
	fig.update_layout(scattermode='group', scattergap=0.75)

	# return add_hippo_logo(fig, in_plot=False)
	return add_punchcard_logo(fig)

HIPPO_LOGO_URL = 'https://raw.githubusercontent.com/mwinokan/HIPPO/main/logos/hippo_logo_tightcrop.png'
HIPPO_HEAD_URL = 'https://raw.githubusercontent.com/mwinokan/HIPPO/main/logos/hippo_assets-02.png'

def add_hippo_logo(fig, in_plot=True):

	assert fig.layout.title.text, 'Figure must have a title to add the HIPPO logo'

	if in_plot:

		fig.add_layout_image(dict(
			source=HIPPO_LOGO_URL,
			xref="paper", yref="paper",
			# layer='below',
			
			x=0.95, y=0.95,
			sizex=0.3, sizey=0.3,
			xanchor="right", yanchor="top",
		))

		return fig

	has_legend = fig.layout.legend.title.text is not None
	fig.layout.margin.t = None

	if has_legend:

		fig.add_layout_image(dict(
			source="",
			xref="paper", yref="paper",
			
			x=1, y=1.05,
			sizex=0.4, sizey=0.4,
			xanchor="left", yanchor="bottom",
		))

	else:

		fig.add_layout_image(dict(
			source=HIPPO_LOGO_URL,
			xref="paper", yref="paper",
			
			x=1, y=1.05,
			sizex=0.3, sizey=0.3,
			xanchor="right", yanchor="bottom",
		))
	
	return fig

def add_punchcard_logo(fig):

	fig.add_layout_image(dict(
		source=HIPPO_HEAD_URL,
		xref="paper", yref="paper",
		
		x=1, y=1,
		sizex=0.25, sizey=0.25,
		xanchor="right", yanchor="top",
	))

	return fig