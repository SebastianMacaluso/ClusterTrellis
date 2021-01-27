import setuptools


setuptools.setup(
    name="ClusterTrellis",
    version="0.0.3",
    description="Hierarchical Cluster Trellis",
    url="https://github.com/SebastianMacaluso/ClusterTrellis",
    author="",
    author_email="sm4511@nyu.edu",
    license="MIT",
    packages=setuptools.find_packages(where="src"),
    package_dir={"": "src"},
    zip_safe=False,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
