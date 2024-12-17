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
    # path("", include("hippo.urls")),
    path("", views.index, name="index"),
    # admin
    path("admin/", admin.site.urls),
    # path("targets/", views.TargetListView.as_view(), name="target-list"),
    # path("target/<int:pk>/", views.TargetDetailView.as_view(), name="target-detail"),
]

for model_name, model_views in views.GENERATED_VIEWS.items():
    urlpatterns += [
        path(
            f"{model_name}/",
            model_views["list_view"].as_view(),
            name=f"{model_name}_list",
        ),
        path(
            f"{model_name}/<int:pk>/",
            model_views["detail_view"].as_view(),
            name=f"{model_name}_detail",
        ),
    ]
