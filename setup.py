#!/usr/bin/env python3
"""IRDB"""

from distutils.core import setup
import os
from os import path as pth
import shutil
import zipfile


ZIPPED_DATA_FILES = (
    ("MaunaKea", "TER_maunakea_tapas.dat"),
)


def ensure_zipped_data_files():
    """Create large ignored data files from their checked-in zip archives."""
    for folder, filename in ZIPPED_DATA_FILES:
        target = pth.join(folder, filename)
        archive = f"{target}.zip"
        if not pth.exists(archive):
            if not pth.exists(target):
                raise FileNotFoundError(f"Missing required data archive: {archive}")
            continue

        with zipfile.ZipFile(archive) as zip_file:
            try:
                zip_info = zip_file.getinfo(filename)
            except KeyError as err:
                raise FileNotFoundError(f"{archive} does not contain {filename}") from err

            if pth.exists(target) and pth.getsize(target) == zip_info.file_size:
                continue

            tmp_target = f"{target}.tmp"
            try:
                with zip_file.open(zip_info) as source, open(tmp_target, "wb") as dest:
                    shutil.copyfileobj(source, dest)
                os.replace(tmp_target, target)
            except Exception:
                if pth.exists(tmp_target):
                    os.remove(tmp_target)
                raise
            print(f"Extracted {archive} -> {target}")


def create_manifest():
    pkgs_list = [pkg for pkg in os.listdir("./") if \
                 pth.isdir(pkg) and pth.exists(pth.join(pkg, f"{pkg}.yaml"))]
    with open("MANIFEST.in", "w") as f:
        for pkg_name in pkgs_list:
            f.write(f"include {pkg_name}/*\n")
            print(f"Including {pkg_name}")


def setup_package():
    ensure_zipped_data_files()
    setup(name = 'IRDB',
          description = "Instrument package database",
          author = "Kieran Leschinski",
          author_email = "kieran.leschinski@unive.ac.at",
          url = "http://homepage.univie.ac.at/kieran.leschinski/",
          package_dir={'scopesim': 'scopesim'},
          packages = ["irdb/tests"],
          install_requires=["numpy>=1.16",
                            "scipy>=1.0.0",
                            "astropy>=2.0",
                            "matplotlib>=1.5",

                            "docutils",
                            "pyyaml>5.1",
                            "paramiko",

                            "scopesim>=0.10.1",
                            ],
          )


if __name__ == '__main__':
    create_manifest()
    setup_package()
