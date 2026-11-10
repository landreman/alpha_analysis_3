"""Install into an isolated temporary environment and import outside the source tree."""

from pathlib import Path
import subprocess
import sys
import tempfile


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    with tempfile.TemporaryDirectory(prefix="alpha-analysis-3-smoke-") as directory:
        environment = Path(directory) / "venv"
        subprocess.run([sys.executable, "-m", "venv", str(environment)], check=True)
        python = environment / "bin" / "python"
        subprocess.run(
            [str(python), "-m", "pip", "install", "--quiet", str(root)], check=True
        )
        subprocess.run(
            [
                str(python),
                "-c",
                "import alpha_analysis_3 as a; "
                "from alpha_analysis_3.boozer_field import BoozerField; "
                "from alpha_analysis_3.conventions import numerical_error_target; "
                "assert numerical_error_target(0.0) == 0.01; "
                "import numba, contourpy, matplotlib; print(a.__file__)",
            ],
            cwd=directory,
            check=True,
        )


if __name__ == "__main__":
    main()
