from setuptools import setup

setup(name="swlocal",
    version="1.0.0",
    author="Keyi Cui",
    author_email="keyi.cui@yale.edu",
    description="Generates best local sequence alignment using Smith-Waterman Algorithm",
    url="https://github.com/keyicui0717/swlocal.git",
    license="GPL",
    packages=["swlocal"],
    install_requires=["numpy", "argparse"]
)



