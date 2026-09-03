#!/usr/bin/env python3
"""Generate a HELM gallery with the locally built standalone Sketcher."""

import argparse
import base64
import importlib.util
import os
from pathlib import Path
import subprocess


def load_gallery_helpers():
    helper_path = Path.home() / "bin" / "generate_helm_gallery.py"
    spec = importlib.util.spec_from_file_location("helm_gallery_helpers", helper_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class LocalRenderer:
    def __init__(self):
        build_root = Path(os.environ["SCHRODINGER"])
        executable = (
            build_root / "sketcher" / "sketcher_app" / "helm_gallery_renderer"
        )
        env = os.environ.copy()
        env.setdefault("QT_QPA_PLATFORM", "offscreen")
        self.process = subprocess.Popen(
            [executable],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            text=True,
            env=env,
        )

    def render(self, helm):
        if helm.endswith("$$") and not helm.endswith("$$$"):
            helm += "$"
        request = base64.b64encode(helm.encode()).decode()
        self.process.stdin.write(request + "\n")
        self.process.stdin.flush()
        response = self.process.stdout.readline()
        if not response:
            raise RuntimeError("Local Sketcher renderer exited unexpectedly")
        success, coords_valid, payload = response.rstrip("\n").split("\t", 2)
        decoded = base64.b64decode(payload).decode()
        return success == "1", decoded, coords_valid == "1"

    def close(self):
        self.process.stdin.close()
        self.process.wait()


def main():
    parser = argparse.ArgumentParser(
        description="Generate a HELM gallery using the local Sketcher build"
    )
    parser.add_argument("csv_files", nargs="+", help="CSV files with HELM entries")
    parser.add_argument("--output", default="helm_gallery.html")
    args = parser.parse_args()

    gallery = load_gallery_helpers()
    entries = gallery.parse_csv_files(args.csv_files)
    print(f"Processing {len(entries)} entries with the local Sketcher build...")
    successes = []
    failures = []
    renderer = LocalRenderer()
    try:
        for entry_num, entry in enumerate(entries, 1):
            success, result, coords_valid = renderer.render(entry["helm_string"])
            if success:
                successes.append((entry_num, entry, result, coords_valid))
            else:
                failures.append(
                    {
                        "number": entry_num,
                        "helmString": entry["helm_string"],
                        "description": entry["description"],
                        "origin": entry["origin"],
                        "error": result,
                    }
                )
    finally:
        renderer.close()

    gallery.create_html_gallery(entries, successes, failures, args.output)


if __name__ == "__main__":
    main()
