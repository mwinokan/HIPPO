from ..models import FileModel
from .file_set import FileSet
from pathlib import Path


class File(FileModel):

    _objects = FileSet.as_manager()
    _parent_module = "files"

    @property
    def name(self):
        return Path(self.path).name

    @property
    def as_path(self):
        return Path(self.path)


def guess_file_format(path) -> str:
    for suffix in File.FILE_FORMATS:
        if path.name.endswith(suffix):
            format_type = suffix
            break
        else:
            raise ValueError(f"Unknown file extension: {path}")
    return format_type
