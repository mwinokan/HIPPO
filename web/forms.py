from django import forms
from hippo.models import FragalysisDownload, SdfUpload, PoseReview, ScoreUpload


class FragalysisDownloadForm(forms.ModelForm):
    class Meta:
        model = FragalysisDownload
        fields = [
            "target_name",
            "target_access_string",
            "access_token",
            "stack",
        ]


class SdfUploadForm(forms.ModelForm):
    class Meta:
        model = SdfUpload
        fields = [
            "target",
            "input_file",
            "protein_field_name",
            "inspirations_field_name",
            "pose_origins",
            "compute_embedding",
            "binding_energy_field_name",
            "inspiration_distance_field_name",
        ]


class ScoreUploadForm(forms.ModelForm):
    class Meta:
        model = ScoreUpload
        fields = [
            "target",
            "score_type",
            "related_model",
            "identifier_column_name",
            "input_file",
            "separator",
            "fields",
        ]


class PoseReviewForm(forms.ModelForm):
    class Meta:
        model = PoseReview
        fields = [
            "review",
        ]


class SearchForm(forms.Form):

    query = forms.CharField(label="Query", max_length=200)
    target = forms.BooleanField(label="Target", initial=True, required=False)
    structure = forms.BooleanField(label="Structure", initial=True, required=False)
    compound = forms.BooleanField(label="Compound", initial=True, required=False)
    pose = forms.BooleanField(label="Pose", initial=True, required=False)
