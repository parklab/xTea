from setuptools import setup, find_packages

setup(
    name="xtea",
    version="0.1.9",
    author="Simon (Chong) Chu",
    author_email="chong.simon.chu@gmail.com",
    maintainer="Corinne Sexton",
    maintainer_email="corinne_sexton@hms.harvard.edu",
    description="Comprehensive TE insertion identification with WGS/WES data from multiple sequencing techniques",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/your-repo-url",  # Replace with your project's URL
    packages=find_packages(),
    python_requires=">=3.6",
    install_requires=[
        "pysam",
        "sortedcontainers",
        "numpy",
        "pandas",
        "pyranges",
        "configargparse",
        "onnxruntime==1.16.3"
    ],
    license="LICENSE",
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    entry_points={
        "console_scripts": [
            "xtea = xtea.run_xtea:main",
        ],
    },
    setup_requires=["setuptools"],
    tests_require=["pytest"],
)

