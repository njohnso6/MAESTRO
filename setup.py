# -*- coding: utf-8 -*-
# Time-stamp: <2022-03-02 12:14:09 ta32852>

import sys,os
# from distutils.core import setup


try:
    from setuptools import setup, find_packages
except ImportError:
    print("Could not load setuptools. Please install the setuptools package.")

def main():
    setup(
        name = "MAESTRO",
        version = "1.5.4",
        package_dir = {'MAESTRO':'MAESTRO'},
        packages = ['MAESTRO'],
        package_data={'MAESTRO':['Snakemake/Multiome/*', 'Snakemake/scRNA/*', 'Snakemake/integrate/*', 'Snakemake/scATAC/*', 'Snakemake/scATAC/rules/*', 'Snakemake/scRNA/rules/*', 'R/*', 'utils/*','annotations/*', 'html/*', '']},
        #data_files = [(os.path.join(sys.prefix,'bin'), ['refpkg/giggle/bin/giggle'])],
        scripts = ['MAESTRO/MAESTRO'],
        include_package_data = True,

        author = "Chenfei Wang, Dongqing Sun, Tao Liu, Changxin Wan, Ming (Tommy) Tang, Gali Bai",
        author_email = "gali.bai@hotmail.com",
        description = "MAESTRO(Model-based AnalysEs of Single-cell Transcriptome and RegulOme) is a comprehensive "
        "single-cell RNA-seq and ATAC-seq analysis suit built using snakemake. (MACS3 fork)",
        license = "GPL-3.0",
        url = "https://github.com/macs3-project/MAESTRO",

        # entry_points = {"console_scripts": ["strap = strap:main"]},
        classifiers = [
            "Development Status :: 4 - Beta",
            #"Development Status :: 5 - Production/Stable",
            "Environment :: Console",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: GPL License",
            "Natural Language :: English",
            "Programming Language :: Python :: 3",
            "Topic :: Scientific/Engineering :: Bio-Informatics"
        ],
        #install_requires=["sinto>=0.7.1",],
        #setup_requires=["sinto>=0.7.1",],
    )


if __name__ == "__main__":
    main()
