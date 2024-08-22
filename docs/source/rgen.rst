
========================
Random Recipe Generation
========================

1. Start with a recipe

::

	recipe = hippo.Recipe.from_json()

2. Create the generator object

::

	gen = hippo.RandomRecipeGenerator()

hippo.RandomRecipeGenerator.get_route_pool() pseudocode:

::

	route_ids = SELECT route_id FROM route WHERE route_product in RandomRecipeGenerator.product_pool
	routes = [db.get_route(id=route_id) from route_id, in route_ids]
	return RouteSet(db, routes)

3. 