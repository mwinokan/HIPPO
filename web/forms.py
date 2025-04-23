from django import forms
from hippo.models import FragalysisDownload, SdfUpload


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
        ]
