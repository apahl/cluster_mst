from setuptools import setup

setup(
    name="cluster_mst",
    version="0.1.0",
    description="Panel web app for Minimum Spanning Tree clustering of HTS results.",
    url="https://github.com/apahl/cluster_mst",
    author="Axel Pahl",
    author_email="",
    license="MIT",
    packages=["cluster_mst"],
    install_requires=[
        "holoviews",
        "panel",
    ],
    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.11",
    ],
)
