from django import forms
from hippo.models import FragalysisDownload


class FragalysisDownloadForm(forms.ModelForm):
    class Meta:
        model = FragalysisDownload
        fields = [
            "target_name",
            "target_access_string",
            "access_token",
            "stack",
        ]
