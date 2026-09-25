"""A selected venv must retain its installed packages across tool adapters."""

import json
import subprocess
import sys
import venv

from fasta_lake.annotation import read_annotation_settings
from fasta_lake.resources import local_interpreter


def test_python_symlink_keeps_the_selected_environment_and_settings(tmp_path):
    environment = tmp_path / "tool environment"
    venv.EnvBuilder(with_pip=False, symlinks=sys.platform != "win32").create(environment)
    python = environment / ("Scripts/python.exe" if sys.platform == "win32" else "bin/python")
    selected = local_interpreter(python)
    prefix = subprocess.check_output(
        [str(selected), "-c", "import sys; print(sys.prefix)"], text=True
    )
    assert prefix.strip() == str(environment)
    (tmp_path / "mapper.py").write_text("# external mapper fixture\n")
    (tmp_path / "data").mkdir()
    config = tmp_path / "settings.json"
    config.write_text(
        json.dumps(
            {
                "emapper": "mapper.py",
                "data_dir": "data",
                "emapper_python": str(python.relative_to(tmp_path)),
            }
        )
    )
    settings = read_annotation_settings(config)
    assert settings["emapper_python"] == python
    prefix = subprocess.check_output(
        [str(settings["emapper_python"]), "-c", "import sys; print(sys.prefix)"], text=True
    )
    assert prefix.strip() == str(environment)
