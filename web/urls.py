"""
URL configuration for hippo project.

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/5.1/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""

from django.contrib import admin
from django.urls import include, path
from . import views

urlpatterns = [
    # index
    path("", views.index, name="index"),
    # admin
    path("accounts/", include("django.contrib.auth.urls")),
    path("admin/", admin.site.urls),
    path("pill_demo/", views.pill_demo, name="pill_demo"),
    # custom views
    # path("model_pill/<str:class_name>/<int:pk>/", views.model_pill, name="model_pill"),
    path("structure_pdb/<int:pk>/", views.structure_pdb, name="structure_pdb"),
    path("pose_sdf/<int:pk>/", views.pose_sdf, name="pose_sdf"),
    path("pose_compare/<str:pks>/", views.pose_compare, name="pose_compare"),
    path("pose_compare_3d/<str:pks>/", views.pose_compare_3d, name="pose_compare_3d"),
    path("fragalysis_download/", views.fragalysis_download, name="fragalysis_download"),
    path("sdf_upload/", views.sdf_upload, name="sdf_upload"),
]

for model_name, model_views in views.GENERATED_VIEWS.items():

    list_view = model_views.get("list_view", None)

    if list_view:
        urlpatterns.append(
            path(
                f"{model_name}/",
                model_views["list_view"].as_view(),
                name=f"{model_name}_list",
            )
        )

    detail_view = model_views.get("detail_view", None)

    if detail_view:

        uri_params = detail_view._uri_param_str
        uri = f"{model_name}/{uri_params}/"

        detail_view = detail_view.as_view()

        urlpatterns.append(
            path(
                uri,
                detail_view,
                name=f"{model_name}_detail",
            )
        )
