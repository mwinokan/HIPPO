
import molparse as mp
import plotly.express as px
import pandas as pd
import mout
import plotly.graph_objects as go

import logging
logger = logging.getLogger('HIPPO')

'''

	ALL GRAPHS DEFINED HERE SHOULD:

	* Have a HIPPO logo
	* Include the target name in the title

'''

import functools

# hippo_graph decorator
def hippo_graph(func):

	@functools.wraps(func)
	def wrapper(animal, *args, **kwargs):

		wrapper_kwargs = {}
		wrapper_keys = ['show', 'html', 'pdf', 'png']
		for key in wrapper_keys:
			wrapper_kwargs[key] = kwargs.pop(key,None)

		fig = func(animal, *args, **kwargs)

		if wrapper_kwargs['show']:
				fig.show()

		if wrapper_kwargs['html']:
			file = wrapper_kwargs['html']
			if not file.endswith('.html'):
				file = f'{file}.html'
			mp.write(file, fig)

		if wrapper_kwargs['pdf']:
			file = wrapper_kwargs['pdf']
			if not file.endswith('.pdf'):
				file = f'{file}.pdf'
			mp.write(file, fig)

		if wrapper_kwargs['png']:
			file = wrapper_kwargs['png']
			if not file.endswith('.png'):
				file = f'{file}.png'
			mp.write(file, fig)

		if not fig.layout.images:
			add_hippo_logo(fig)

		return fig
	
	return wrapper

@hippo_graph
def plot_tag_statistics(animal, color='type', subtitle=None, log_y=False):

	plot_data = []

	for tag in animal.tags.unique:

		num_compounds = len(animal.compounds.get_by_tag(tag=tag))
		data = dict(tag=tag, number=num_compounds, type='compounds')
		plot_data.append(data)

		num_poses = len(animal.poses.get_by_tag(tag=tag))
		data = dict(tag=tag, number=num_poses, type='poses')
		plot_data.append(data)
		
	fig = px.bar(plot_data, x='tag', y='number', color=color, log_y=log_y)

	title = 'Tag Statistics'

	if subtitle:
		title = f'<b>{animal.name}</b>: {title}<br><sup><i>{subtitle}</i></sup>'
	else:
		title = f'<b>{animal.name}</b>: {title}'
	
	fig.update_layout(title=title,title_automargin=False, title_yref='container', barmode='group')

	fig.update_layout(xaxis_title='Tag', yaxis_title='#')

	return fig

@hippo_graph
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
		title = f'<b>{animal.name}</b>: {title}<br><sup><i>{subtitle}</i></sup>'
	else:
		title = f'<b>{animal.name}</b>: {title}'
	
	fig.update_layout(title=title,title_automargin=False, title_yref='container')

	fig.update_layout(xaxis_title='Residue')
	fig.update_layout(yaxis_title='#Interactions')

	return fig

@hippo_graph
def plot_interaction_punchcard(animal, poses, subtitle=None, opacity=1.0, group='pose_id'):

	import plotly

	plot_data = []

	categoryarray = {}

	logger.debug(f'{poses=}')

	for pose in poses:

		fingerprint = pose.fingerprint

		if not fingerprint:
			continue

		# loop over each interaction in the pose
		for key, value in fingerprint.items():
			if not value:
				continue
			
			data = dict(str=key,count=value)

			f = animal.db.get_feature(id=key)

			data['id'] = f.id
			data['family'] = f.family
			data['residue_name'] = f.residue_name
			data['residue_number'] = f.residue_number
			data['chain_name'] = f.chain_name
			data['atom_names'] = f.atom_names
			
			data['pose_name'] = pose.name
			data['pose_longname'] = pose.longname
			data['pose_id'] = str(pose)

			data['chain_res_name_number_str'] = f.chain_res_name_number_str

			if data['chain_res_name_number_str'] not in categoryarray:
				categoryarray[data['chain_res_name_number_str']] = (data['chain_res_name_number_str'], f.chain_name, f.residue_number)

			plot_data.append(data)

	plot_df = pd.DataFrame(plot_data)

	plot_df = plot_df.sort_values(group)

	title = 'Interaction Punch-Card'

	if subtitle:
		title = f'<b>{animal.name}</b>: {title}<br><sup><i>{subtitle}</i></sup>'
	else:
		title = f'<b>{animal.name}</b>: {title}'

	fig = px.scatter(plot_df, x='chain_res_name_number_str', y='family',marginal_x='histogram',marginal_y='histogram', hover_data=plot_df.columns, color=group, title=title)
	
	fig.update_layout(title=title,title_automargin=False, title_yref='container')

	fig.update_layout(xaxis_title='Residue', yaxis_title='Feature Family')

	categoryarray = sorted([v for v in categoryarray.values()], key=lambda x: (x[1], x[2]))
	categoryarray = [v[0] for v in categoryarray]

	fig.update_xaxes(categoryorder='array', categoryarray=categoryarray)

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

@hippo_graph
# def plot_building_blocks(animal, subtitle=None, cset='elabs', color='name_is_smiles'):
def plot_reactant_amounts(animal, subtitle=None, color='has_price_picker', named_only=False, most_common=None):

	# cset = animal.compound_sets[cset]

	# bbs = cset.get_building_blocks()

	bbs = animal.building_blocks

	mout.debug('making plot_data')
	plot_data = []
	for bb in bbs:
		d = bb.dict
		# if most_common and d['amount'] is not None:
			
		if not named_only or not d['name_is_smiles']:
			plot_data.append(d)

	if most_common:
		mout.debug('sorting')
		plot_data = sorted(plot_data, key=lambda x: x['amount'], reverse=True)[:most_used_number]

	fig = px.bar(plot_data, x='name', y='amount', color=color, hover_data=plot_data[0].keys())
	# fig = px.bar(plot_data, x='smiles', y='amount', color=color)

	title = 'Building Blocks'

	if not subtitle:
		# subtitle = f'"{cset.name}": #BBs={len(bbs)}'
		subtitle = f'#BBs={len(bbs)}'

	title = f'<b>{animal.name}</b>: {title}<br><sup><i>{subtitle}</i></sup>'

	fig.update_layout(title=title,title_automargin=False, title_yref='container')

	fig.update_layout(xaxis_title='Reactant', yaxis_title='#Reactions')

	return fig

@hippo_graph
# def plot_building_blocks(animal, subtitle=None, cset='elabs', color='name_is_smiles'):
def plot_reactant_price(animal, subtitle=None, amount=20):

	# cset = animal.compound_sets[cset]

	# bbs = cset.get_building_blocks()

	bbs = animal.building_blocks

	plot_data = []
	for bb in bbs:
		d = bb.dict
		if not d['has_price_picker']:
			continue

		d[f'price_{amount}mg'] = bb.get_price(amount)
		d[f'min_amount'] = bb.price_picker.min_amount

		plot_data.append(d)

	# fig = px.bar(plot_data, x='name', y=f'price_{amount}mg', color='lead_time', log_y=True, hover_data=plot_data[0].keys())
	fig = px.histogram(plot_data, x=f'price_{amount}mg', color='lead_time', hover_data=plot_data[0].keys())
	# fig = px.bar(plot_data, x='smiles', y='amount', color=color)

	title = 'Reactant Pricing'

	if not subtitle:
		# subtitle = f'"{cset.name}": #BBs={len(bbs)}'
		subtitle = f'#BBs={len(bbs)}'

	title = f'<b>{animal.name}</b>: {title}<br><sup><i>{subtitle}</i></sup>'

	fig.update_layout(title=title,title_automargin=False, title_yref='container')

	fig.update_layout(yaxis_title='Number of reactants', xaxis_title=f'Price for {amount}mg [$USD]')

	return fig

@hippo_graph
# def plot_building_blocks(animal, subtitle=None, cset='elabs', color='name_is_smiles'):
def plot_reactants_2d(animal, subtitle=None, amount=20):

	# cset = animal.compound_sets[cset]

	# bbs = cset.get_building_blocks()

	bbs = animal.building_blocks

	plot_data = []
	for bb in bbs:
		d = bb.dict
		if not d['has_price_picker']:
			continue

		d[f'price_{amount}mg'] = bb.get_price(amount)
		d[f'min_amount'] = bb.price_picker.min_amount

		plot_data.append(d)

	fig = px.scatter(plot_data, y='amount', x=f'price_{amount}mg', color='name', hover_data=plot_data[0].keys())
	# fig = px.bar(plot_data, x='smiles', y='amount', color=color)

	title = 'Building Blocks'

	if not subtitle:
		# subtitle = f'"{cset.name}": #BBs={len(bbs)}'
		subtitle = f'#BBs={len(bbs)}'

	title = f'<b>{animal.name}</b>: {title}<br><sup><i>{subtitle}</i></sup>'

	fig.update_layout(title=title,title_automargin=False, title_yref='container')

	fig.update_layout(scattermode='group', scattergap=0.75)

	fig.update_layout(yaxis_title='Quantity [mg]', xaxis_title=f'Price for {amount}mg [$USD]')

	return fig

@hippo_graph
# def plot_building_blocks(animal, subtitle=None, cset='elabs', color='name_is_smiles'):
def plot_building_blocks(animal, subtitle=None, color='name_is_smiles'):

	# cset = animal.compound_sets[cset]

	# bbs = cset.get_building_blocks()

	bbs = animal.building_blocks

	plot_data = []
	for bb in bbs:
		plot_data.append(bb.dict)

	fig = px.scatter(plot_data, x='name', y='max', color='amount')
	# fig = px.bar(plot_data, x='smiles', y='amount', color=color)

	title = 'Building Blocks'

	if not subtitle:
		# subtitle = f'"{cset.name}": #BBs={len(bbs)}'
		subtitle = f'#BBs={len(bbs)}'

	title = f'<b>{animal.name}</b>: {title}<br><sup><i>{subtitle}</i></sup>'

	fig.update_layout(title=title,title_automargin=False, title_yref='container')

	fig.update_layout(xaxis_title='Reactant', yaxis_title='Quantity')

	return fig

@hippo_graph
def plot_synthetic_routes(animal, subtitle=None, cset='elabs', color='num_reactants'):

	cset = animal.compound_sets[cset]

	plot_data = []
	for reax in cset.reactions:
		plot_data.append(reax.dict)

	# fig = px.bar(plot_data, x='name', y='amount', color=color)
	fig = px.histogram(plot_data, x='type', color=color)

	title = 'Synthetic Routes'

	if not subtitle:
		subtitle = f'"{cset.name}": #compounds={len(cset)}'

	title = f'<b>{animal.name}</b>: {title}<br><sup><i>{subtitle}</i></sup>'

	fig.update_layout(title=title,title_automargin=False, title_yref='container')

	fig.update_layout(xaxis_title='Compound', yaxis_title='#Routes')

	return fig

@hippo_graph
def plot_numbers(animal, subtitle=None):

	'''
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
	'''

	# cset = animal.compound_sets[cset]

	plot_data = [
		dict(category='Experimental Hits', number=len(animal.hits), type='compound'),
		dict(category='Experimental Hits', number=len(animal.hits.poses), type='poses'),
		dict(category='Base compounds', number=len(animal.bases), type='compound'),
		dict(category='Base compounds', number=len(animal.bases.poses), type='poses'),
		dict(category='Syndirella Elaborations', number=len(animal.elabs), type='compound'),
		dict(category='Syndirella Elaborations', number=len(animal.elabs.poses), type='poses'),
		dict(category='Unique Reactants', number=len(animal.building_blocks), type='compound'),
	]

	fig = px.bar(plot_data, x='category', y='number', log_y=True, color='type')

	title = 'Compounds & Poses'

	title = f'<b>{animal.name}</b>: {title}<br>'

	fig.update_layout(title=title,title_automargin=False, title_yref='container', barmode='group')

	fig.update_layout(xaxis_title=None, yaxis_title='Log(Quantity)')

	return fig

@hippo_graph
def plot_compound_property(animal, prop, style='bar', null=None):

	"""
	Get an arbitrary property from all the compounds in animal.compounds
	
	If one property, plot a 1D histogram
	If 2D plot a bar/scatter

	"""

	if not isinstance(prop, list):
		prop = [prop]

	plot_data = []

	for comp in animal.compounds:

		data = comp.dict

		for p in prop:

			# has attr

			if p not in data:

				# get attr
				if hasattr(comp,p):
					v = getattr(comp,p)
				elif p in comp.metadata:
					v = comp.metadata[p]
				else:
					v = null

				data[p] = v

		plot_data.append(data)

	if len(prop) == 1:

		title = f'Compound {prop[0]}'

		fig = px.histogram(plot_data, x=prop[0])
		
		fig.update_layout(xaxis_title=prop[0], yaxis_title='Quantity')

	elif len(prop) == 2:

		title = f'Compound {prop[0]} vs {prop[1]}'
		
		func = eval(f'px.{style}')
		fig = func(plot_data, x=prop[0], y=prop[1])

		fig.update_layout(xaxis_title=prop[1], yaxis_title=prop[1])

	else:
		mout.error('Unsupported', code='plotting.plot_compound_property.1')

	title = f'<b>{animal.name}</b>: {title}<br>'

	fig.update_layout(title=title,title_automargin=False, title_yref='container', barmode='group')

	return fig

@hippo_graph
def plot_pose_property(animal, prop, poses=None, style='bar', null=None):

	"""
	Get an arbitrary property from all the poses in animal.poses
	
	If one property, plot a 1D histogram
	If 2D plot a bar/scatter

	"""

	if not isinstance(prop, list):
		prop = [prop]

	if not poses:
		poses = animal.poses

	plot_data = []

	for pose in poses:

		data = pose.dict

		for p in prop:

			# has attr

			if p not in data:

				# get attr
				if hasattr(pose,p):
					v = getattr(pose,p)
				elif p in pose.metadata:
					v = pose.metadata[p]
				else:
					v = null

				data[p] = v

		plot_data.append(data)

	if len(prop) == 1:

		title = f'Pose {prop[0]}'

		fig = px.histogram(plot_data, x=prop[0])
		
		fig.update_layout(xaxis_title=prop[0], yaxis_title='Quantity')

	elif len(prop) == 2:

		title = f'Pose {prop[0]} vs {prop[1]}'
		
		func = eval(f'px.{style}')
		fig = func(plot_data, x=prop[0], y=prop[1])

		fig.update_layout(xaxis_title=prop[1], yaxis_title=prop[1])

	else:
		mout.error('Unsupported', code='plotting.plot_pose_property.1')

	title = f'<b>{animal.name}</b>: {title}<br>'

	fig.update_layout(title=title,title_automargin=False, title_yref='container', barmode='group')

	return fig

@hippo_graph
def plot_reactant_sankey(animal, subtitle):

	'''
		BBs (total)
		BBs (in enamine)
		BBs (within budget)
		BBs (within lead-time)
	'''

	total = animal.building_blocks
	quoted = [bb for bb in animal.building_blocks if bb.price_picker is not None]
	lead_time = [bb for bb in quoted if bb.lead_time <= animal.max_lead_time]
	budget = [bb for bb in quoted if bb.get_price(animal.min_bb_quantity) <= animal.max_bb_price]
	bad = [bb for bb in quoted if bb not in lead_time and bb not in budget]
	
	n_total = len(total)
	n_quoted = len(quoted)
	n_quote_not_attempted = len([bb for bb in total if 'quote_attempted' not in bb.tags])
	n_quote_attempted = n_total - n_quote_not_attempted
	n_ok = len([bb for bb in quoted if bb in lead_time and bb in budget])
	# n_bad = len(bad)
	n_expensive = len([bb for bb in quoted if bb not in bad and bb not in budget and bb in lead_time])
	n_slow = len([bb for bb in quoted if bb not in bad and bb in budget and bb not in lead_time])
	n_bad = len(bad)

	print(n_quote_not_attempted)

	labels = [
		f'Total = {n_total}',																# 0
		f'Quote Succeeded = {n_quoted}',														# 1
		f'Quote Failed = {n_quote_attempted - n_quoted}',											# 2
		f'OK = {n_ok}',																			# 3
		f'Too expensive ({animal.min_bb_quantity}mg > ${animal.max_bb_price}) = {n_expensive}',	# 4
		f'Too slow (lead_time > {animal.max_lead_time} days) = {n_slow}',						# 5
		f'Too expensive & slow = {n_bad}',													# 6
		f'Quote not attempted = {n_quote_not_attempted}', # 7
		f'Quote attempted = {n_quote_attempted}' #8
	]

	links = [
		[8, 1, n_quoted], # total --> quote success
		[8, 2, n_quote_attempted - n_quoted], # total --> not quoted
		[1, 3, n_ok], # quoted --> ok
		[1, 4, n_expensive], # quoted --> too expensive
		[1, 5, n_slow], # quoted --> too slow
		[1, 6, len(bad)], # quoted --> bad
		[0, 7, n_quote_not_attempted], # total --> quote_not_attempted
		[0, 8, n_quote_attempted], # total --> quote_attempted
	]

	source = [l[0] for l in links]
	target = [l[1] for l in links]
	value = [l[2] for l in links]

	fig = go.Figure(data=[go.Sankey(
			node = dict(
			# pad = 15,
			# thickness = 20,
			# line = dict(color = "black", width = 0.5),
			label = labels,
			# color = "blue"
		),
			link = dict(
			source = source,
			target = target,
			value = value,
	))])

	fig.update_layout(font_size=10)

	title = 'Reactant Quoting'

	subtitle = animal.reactant_catalog

	title = f'<b>{animal.name}</b>: {title}<br><sup><i>{subtitle}</i></sup>'

	fig.update_layout(title=title,title_automargin=False, title_yref='container')

	return fig

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