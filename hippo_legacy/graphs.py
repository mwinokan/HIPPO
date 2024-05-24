
import plotly.express as px
import plotly.graph_objects as go

def hits_by_site(animal):
	fig = px.histogram(animal.compound_df,x='_site_name')
	return fig
